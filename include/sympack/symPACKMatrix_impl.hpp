
#ifndef _SUPERNODAL_MATRIX_IMPL_HP_
#define _SUPERNODAL_MATRIX_IMPL_HPP_

#include <sympack/symPACKMatrix.hpp>

#include "sympack/blas.hpp"
#include "sympack/Ordering.hpp"
#include "sympack/LoadBalancer.hpp"
#include "sympack/mpi_interf.hpp"

#include  "sympack/DistSparseMatrixGraph.hpp"

#include <queue>
#include <stdlib.h>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

namespace symPACK{
  extern "C" {
    void FORTRAN(ordsup) (int * ordflag, int *  altflag, int *  NEQNS, int *  nofsub, int *  nsuper, 
        int * xsuper, int *  xlindx, int *  lindx , int *  snode , int *  perm  , 
        int * invp  , int *  freeforw, int *  freeback, int *  sforw, int *  sback, 
        int * setseg_forw, int *  setseg_back, int *  nodehead, 
        int * nodeforw, int *  nodeback, 
        int *  setsnode, int *  supperm, int *  mark, int *  set  , int *  compset,
        int *  invp2 , int *  heap                             );
  }
}




namespace symPACK{

  supernodalTaskGraph::supernodalTaskGraph(){
    localTaskCount_=0;
  }
  supernodalTaskGraph::supernodalTaskGraph( const supernodalTaskGraph& g ){
    (*this) = g;
  }

  supernodalTaskGraph& supernodalTaskGraph::operator=( const supernodalTaskGraph& g ){
    localTaskCount_ = g.localTaskCount_;
    taskLists_.resize(g.taskLists_.size(),NULL);
    for(int i = 0; i<taskLists_.size(); ++i){
      if(g.taskLists_[i] != NULL){
        taskLists_[i] = new std::list<FBTask>(); 
        taskLists_[i]->insert(taskLists_[i]->end(),g.taskLists_[i]->begin(),g.taskLists_[i]->end());
      }
    }

#ifdef MULTITHREADING
    superLocks_.resize( g.taskLists_.size() );
    std::for_each(superLocks_.begin(), superLocks_.end(), [] (omp_lock_t & lock) { omp_init_lock(&lock);});
#endif

  }

  supernodalTaskGraph::~supernodalTaskGraph(){
    for(int i = 0; i<taskLists_.size(); ++i){
      if(taskLists_[i] != NULL){
        delete taskLists_[i];
      }
    }

#ifdef MULTITHREADING
    std::for_each(superLocks_.begin(), superLocks_.end(), [] (omp_lock_t & lock) { omp_destroy_lock(&lock);});
#endif
  }        

  void supernodalTaskGraph::removeTask(std::list<FBTask>::iterator & taskit){
    taskLists_[taskit->tgt_snode_id-1]->erase(taskit);
    //decreaseTaskCount();
  }

  std::list<FBTask>::iterator supernodalTaskGraph::addTask(FBTask & task){
    if(taskLists_[task.tgt_snode_id-1] == NULL){
      taskLists_[task.tgt_snode_id-1]=new std::list<FBTask>();
    }
    taskLists_[task.tgt_snode_id-1]->push_back(task);
    increaseTaskCount();

    std::list<FBTask>::iterator taskit = --taskLists_[task.tgt_snode_id-1]->end();
    return taskit;
  }



  Int supernodalTaskGraph::getTaskCount()
  {
    Int val = localTaskCount_;
    return val;
  }

  Int supernodalTaskGraph::setTaskCount(Int value)
  {
    localTaskCount_ = value;
    return value;
  }

  Int supernodalTaskGraph::increaseTaskCount()
  {
    Int val;
#pragma omp atomic capture
    val = ++localTaskCount_;
    return val;
  }

  Int supernodalTaskGraph::decreaseTaskCount()
  {
    Int val;
#pragma omp atomic capture
    val = --localTaskCount_;
    return val;
  }

  std::list<FBTask>::iterator supernodalTaskGraph::find_task(Int src, Int tgt, TaskType type )
  {
    scope_timer(a,FB_FIND_TASK);
    //find task corresponding to curUpdate
    auto taskit = taskLists_[tgt-1]->begin();
    for(;
        taskit!=taskLists_[tgt-1]->end();
        taskit++){
      if(taskit->src_snode_id==src && taskit->tgt_snode_id==tgt && taskit->type==type){
        break;
      }
    }
    return taskit;
  }


  template<typename T> void symPACKMatrix<T>::generateTaskGraph(supernodalTaskGraph & taskGraph,
      std::vector<Int> & AggregatesToRecv,  std::vector<Int>& LocalAggregates)
  {
    //we will need to communicate if only partial xlindx_, lindx_
    //idea: build tasklist per processor and then exchange
    std::map<Idx, std::list<std::pair<Idx,Idx> >  > Updates;
    std::vector<int> marker(np,0);
    Int numLocSnode = XsuperDist_[iam+1]-XsuperDist_[iam];
    Int firstSnode = XsuperDist_[iam];
    for(Int locsupno = 1; locsupno<locXlindx_.size(); ++locsupno){
      Idx I = locsupno + firstSnode-1;
      Int iOwner = Mapping_->Map(I-1,I-1);

      Int first_col = Xsuper_[I-1];
      Int last_col = Xsuper_[I]-1;

      //Create the factor task on the owner
      Updates[iOwner].push_back(std::make_pair(I,I));

      Int J = -1; 
      Ptr lfi = locXlindx_[locsupno-1];
      Ptr lli = locXlindx_[locsupno]-1;
      Idx prevSnode = -1;
      for(Ptr sidx = lfi; sidx<=lli;sidx++){
        Idx row = locLindx_[sidx-1];
        J = SupMembership_[row-1];
        if(J!=prevSnode){
          //J = locSupLindx_[sidx-1];
          Int iUpdater = Mapping_->Map(J-1,I-1);

          if(J>I){
            if(marker[iUpdater]!=I){
              //create the task on iUpdater
              Updates[iUpdater].push_back(std::make_pair(I,J));
              //create only one task if I don't own the factor
              if(iUpdater!=iOwner){
                marker[iUpdater]=I;
              }
            }
          }
        }
        prevSnode = J;
      }
    }

    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(std::pair<Idx,Idx>), MPI_BYTE, &type );
    MPI_Type_commit(&type);

    //then do an alltoallv
    //compute send sizes
    vector<int> ssizes(np,0);
    for(auto itp = Updates.begin();itp!=Updates.end();itp++){
      ssizes[itp->first] = itp->second.size();//*sizeof(std::pair<Idx,Idx>);
    }

    //compute send displacements
    vector<int> sdispls(np+1,0);
    sdispls[0] = 0;
    std::partial_sum(ssizes.begin(),ssizes.end(),&sdispls[1]);

    //Build the contiguous array of pairs
    vector<std::pair<Idx,Idx> > sendbuf;
    sendbuf.reserve(sdispls.back());///sizeof(std::pair<Idx,Idx>));

    for(auto itp = Updates.begin();itp!=Updates.end();itp++){
      sendbuf.insert(sendbuf.end(),itp->second.begin(),itp->second.end());
    }

    //gather receive sizes
    vector<int> rsizes(np,0);
    MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,CommEnv_->MPI_GetComm());

    //compute receive displacements
    vector<int> rdispls(np+1,0);
    rdispls[0] = 0;
    std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);


    //Now do the alltoallv
    vector<std::pair<Idx,Idx> > recvbuf(rdispls.back());///sizeof(std::pair<Idx,Idx>));
    MPI_Alltoallv(&sendbuf[0],&ssizes[0],&sdispls[0],type,&recvbuf[0],&rsizes[0],&rdispls[0],type,CommEnv_->MPI_GetComm());

    MPI_Type_free(&type);

    //now process recvbuf and add the tasks to my tasklists
    taskGraph.taskLists_.resize(TotalSupernodeCnt(),NULL);
    Int localTaskCount = 0;
    for(auto it = recvbuf.begin();it!=recvbuf.end();it++){
      FBTask curUpdate;
      curUpdate.src_snode_id = it->first;
      curUpdate.tgt_snode_id = it->second;
      curUpdate.type=(it->first==it->second)?FACTOR:UPDATE;

      Int J = curUpdate.tgt_snode_id;

      //create the list if needed
      if(taskGraph.taskLists_[J-1] == NULL){
        taskGraph.taskLists_[J-1]=new std::list<FBTask>();
      }
      taskGraph.taskLists_[J-1]->push_back(curUpdate);
      localTaskCount++;
    }
    taskGraph.localTaskCount_=localTaskCount;


#ifdef MULTITHREADING
    taskGraph.superLocks_.resize( TotalSupernodeCnt() );
    std::for_each(taskGraph.superLocks_.begin(), taskGraph.superLocks_.end(), [] (omp_lock_t & lock) { omp_init_lock(&lock);});
#endif

    for(int i = 0; i<taskGraph.taskLists_.size(); ++i){
      if(taskGraph.taskLists_[i] != NULL){
        for(auto taskit = taskGraph.taskLists_[i]->begin(); taskit!=taskGraph.taskLists_[i]->end();taskit++){
          FBTask & curUpdate = *taskit;
          if(curUpdate.type==FACTOR){
            curUpdate.remote_deps = AggregatesToRecv[curUpdate.tgt_snode_id-1];
            curUpdate.local_deps = LocalAggregates[curUpdate.tgt_snode_id-1];
          }
          else if(curUpdate.type==UPDATE){
            Int iOwner = Mapping_->Map(curUpdate.src_snode_id-1,curUpdate.src_snode_id-1);
            //If I own the factor, it is a local dependency
            if(iam==iOwner){
              curUpdate.remote_deps = 0;
              curUpdate.local_deps = 1;
            }
            else{
              curUpdate.remote_deps = 1;
              curUpdate.local_deps = 0;
            }
          }
        }
      }
    }


    std::vector<Int> levels;
    levels.resize(Xsuper_.size());
    levels[0]=0;
    Int numLevel = 0; 
    for(Int i=Xsuper_.size()-1-1; i>=0; i-- ){     
      Int fcol = Xsuper_[i];
      bassert(fcol-1>=0);
      Int pcol =ETree_.PostParent(fcol-1);
      Int supno = pcol>0?SupMembership_[pcol-1]:0;
      levels[i] = levels[supno]+1;
    }

    for(int i = 0; i<taskGraph.taskLists_.size(); ++i){
      if(taskGraph.taskLists_[i] != NULL){
        for(auto taskit = taskGraph.taskLists_[i]->begin(); taskit!=taskGraph.taskLists_[i]->end();taskit++){
          taskit->rank = levels[taskit->src_snode_id];
        }
      }
    }


  }


  template <typename T> void symPACKMatrix<T>::Factorize(){
    SYMPACK_TIMER_START(NUMERICAL_FACT);
    if(iam<np){
      switch(options_.factorization){
        case FANBOTH:
          FanBoth();
          break;
        case FANBOTH_STATIC:
          FanBoth_Static();
          break;
        default:
          FanBoth();
          break;
      }
    }
    SYMPACK_TIMER_STOP(NUMERICAL_FACT);

#ifdef PROFILE_COMM
    if(iam<np){
      logfileptr->OFS()<<"Local volume of communication: "<<gVolComm<<std::endl;
      logfileptr->OFS()<<"Local number of messages: "<<gNumMsg<<std::endl;

      size_t totalVolComm = 0;
      MPI_Reduce(&gVolComm,&totalVolComm,1,MPI_UINT64_T,MPI_SUM,0,CommEnv_->MPI_GetComm());
      size_t totalNumMsg = 0;
      MPI_Reduce(&gNumMsg,&totalNumMsg,1,MPI_UINT64_T,MPI_SUM,0,CommEnv_->MPI_GetComm());

      if(iam==0){
        std::cout<<"Total volume of communication: "<<totalVolComm<<std::endl;
        std::cout<<"Total number of messages: "<<totalNumMsg<<std::endl;
      }

      gVolComm=0;
      gNumMsg=0;
    }
#endif
  }


  //Solve related routines
  template <typename T> void symPACKMatrix<T>::Solve(T * RHS, int nrhs,  T * Xptr) {
    SYMPACK_TIMER_START(SPARSE_SOLVE);

    Int n = iSize_;
    //Int iam = CommEnv_->MPI_Rank();
    //Int np  = CommEnv_->MPI_Size();

    if(iam<np){
      std::vector<Int> children(Xsuper_.size());
      SetValue(children,0);

      for(Int I=1;I<Xsuper_.size()-1;I++){
        Int fc = Xsuper_[I-1];
        Int lc = Xsuper_[I]-1;

        Int parent = ETree_.PostParent(lc-1);
        if(parent!=0){
          ++children[SupMembership_[parent-1]-1];
        }
      }

#ifdef _DEBUG_
      logfileptr->OFS()<<"Children std::vector is"<<children<<std::endl;
#endif


      std::vector<Int> UpdatesToDo = children;

      Contributions_.resize(LocalSupernodes_.size());
      std::vector<std::stack<Int> > LocalUpdates(LocalSupernodes_.size());

      AsyncComms outgoingSend;

      //gdb_lock(3);
      //This corresponds to the k loop in dtrsm
      for(Int I=1;I<Xsuper_.size();I++){
        Int iOwner = this->Mapping_->Map(I-1,I-1);
        //If I own the column, factor it
        if( iOwner == iam ){
          //Create the blocks of my contrib with the same nz structure as L
          //and the same width as the final solution
          //MEMORY CONSUMPTION TOO HIGH ?
          //Int iLocalI = (I-1) / np +1 ;
          Int iLocalI = snodeLocalIndex(I);
          SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];
          SuperNode<T,MallocAllocator> * contrib = CreateSuperNode<MallocAllocator>(options_.decomposition,I,1,nrhs, cur_snode->NRowsBelowBlock(0) ,iSize_);

          Contributions_[iLocalI-1] = contrib;


          for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){
            NZBlockDesc & cur_desc = cur_snode->GetNZBlockDesc(blkidx);
            contrib->AddNZBlock(cur_snode->NRows(blkidx),nrhs,cur_desc.GIndex);
          }
          Int nRows = contrib->NRowsBelowBlock(0);
          std::fill(contrib->GetNZval(0),contrib->GetNZval(0)+nRows*nrhs,ZERO<T>());

          contrib->Shrink();
        }
      }

      CommList ContribsToSend; 
      DownCommList ContribsToSendDown; 


      std::vector<char> src_blocks;


      //forward-substitution phase
      //Sending contrib up the tree
      //Start from the leaves of the tree
      SYMPACK_TIMER_START(SPARSE_FWD_SUBST);

      Int I =1;
      Int iLocalI =1;
      while(iLocalI<=LocalSupernodes_.size() || !ContribsToSend.empty() || !outgoingSend.empty()){
        //Check for completion of outgoing communication
        AdvanceOutgoing(outgoingSend);

        //process some of the delayed send

        if(iLocalI>0 && iLocalI<=LocalSupernodes_.size()){

          //If I own the column, factor it
          SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];
          I = cur_snode->Id();
          Int parent = ETree_.PostParent(cur_snode->LastCol()-1);
          SuperNode<T,MallocAllocator> * contrib = Contributions_[iLocalI-1];


          //Do all my updates (Local and remote)
          //Local updates
          while(!LocalUpdates[iLocalI-1].empty()){
            Int contrib_snode_id = LocalUpdates[iLocalI-1].top();
            LocalUpdates[iLocalI-1].pop();

            SuperNode<T,MallocAllocator> * dist_contrib = snodeLocal(contrib_snode_id,Contributions_);

#ifdef _DEBUG_
            logfileptr->OFS()<<"LOCAL Supernode "<<I<<" is updated by contrib of Supernode "<<contrib_snode_id<<std::endl;
#endif

            Int iOwner = iam;
            contrib->forward_update(dist_contrib,iOwner,iam);
            --UpdatesToDo[I-1];
          }

          //do remote updates
          size_t max_bytes;
          Int nz_cnt;
          while(UpdatesToDo[I-1]>0){
            //receive children contrib
#ifdef _DEBUG_
            logfileptr->OFS()<<UpdatesToDo[I-1]<<" contribs left"<<std::endl;
#endif


            MPI_Status recv_status;
            int bytes_received = 0;

            MPI_Probe(MPI_ANY_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
            MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
            src_blocks.resize(bytes_received);
            MPI_Recv(&src_blocks[0],bytes_received,MPI_BYTE,recv_status.MPI_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
            SuperNode<T,MallocAllocator> * dist_contrib = CreateSuperNode<MallocAllocator>(options_.decomposition,&src_blocks[0],bytes_received);
            //TODO Put this back
            //Deserialize(&src_blocks[0],dist_contrib);
#ifdef _DEBUG_
            logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib.Id()<<std::endl;
#endif

            Int iOwner = recv_status.MPI_SOURCE;
            contrib->forward_update(dist_contrib,iOwner,iam);

            --UpdatesToDo[I-1];
            delete dist_contrib;
          }

          assert(UpdatesToDo[I-1]==0);

          if(UpdatesToDo[I-1]==0){

            //now compute MY contribution
            contrib->forward_update_contrib(RHS,cur_snode,this->Order_.perm);

            //send to my parent
            if(parent!=0){
              Int parent_snode_id = SupMembership_[parent-1];

              Int iTarget = this->Mapping_->Map(parent_snode_id-1,parent_snode_id-1);

              if(iTarget!=iam){
#ifdef _DEBUG_
                logfileptr->OFS()<<"Remote Supernode "<<parent_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif

                Int tgt_first_col = Xsuper_[parent_snode_id-1];
                Int tgt_last_col = Xsuper_[parent_snode_id]-1;
                Int src_nzblk_idx = 1;
                NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);

                Int src_first_row = pivot_desc.GIndex;

                bool isSkipped= false;

                Int next_local_contrib = (iLocalI < Contributions_.size())?Contributions_[iLocalI]->Id():Xsuper_.size();
                if(next_local_contrib< parent_snode_id || true){
                  //need to push the prev src_last_row
                  ContribsToSend.push(DelayedComm(contrib->Id(),parent_snode_id,1,src_first_row));
#ifdef _DEBUG_DELAY_
                  std::cout<<"P"<<iam<<" has delayed update from Contrib "<<I<<" to "<<parent_snode_id<<" from row "<<src_first_row<<std::endl;
#endif
                  isSkipped= true;
                }

                if(!isSkipped){
                  //Create a new Icomm buffer, serialize the contribution
                  // in it and add it to the outgoing comm list

                  Icomm * send_buffer = new Icomm();
                  contrib->Serialize(*send_buffer,src_nzblk_idx,src_first_row);
                  AddOutgoingComm(outgoingSend,send_buffer);

                  //{
                  //  char * buf = outgoingSend.back()->front();
                  //  size_t bytes_received = outgoingSend.back()->size();
                  //          SuperNode<T,MallocAllocator> * dist_contrib = CreateSuperNode<MallocAllocator>(options_.decomposition,&buf[0],bytes_received);
                  //  delete dist_contrib;
                  //}



                  if( outgoingSend.size() > maxIsend_){
                    MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,parent_snode_id,CommEnv_->MPI_GetComm());
                    //that deletes the Icomm
                    outgoingSend.pop_back();
                  }
                  else{
                    MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,parent_snode_id,CommEnv_->MPI_GetComm(),&outgoingSend.back()->Request);
                  }
#ifdef _DEBUG_            
                  logfileptr->OFS()<<"     Send contribution "<<I<<" to Supernode "<<parent_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
                  //                logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
                  //                logfileptr->OFS()<<"Sending "<<nzblk_cnt*sizeof(NZBlockDesc2)+nz_cnt*sizeof(T) + 3*sizeof(Int)<<" bytes"<< std::endl;
#endif
                }
              }
              else{
#ifdef _DEBUG_
                logfileptr->OFS()<<"Local Supernode "<<parent_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif

                Int iLocalJ = snodeLocalIndex(parent_snode_id);
                LocalUpdates[iLocalJ-1].push((Int)I);
              }
            }

          }
        }


        SendDelayedMessagesUp(iLocalI,ContribsToSend,outgoingSend,Contributions_);

        ++iLocalI;
      }

      while(!outgoingSend.empty()){
        AdvanceOutgoing(outgoingSend);
      } 
      MPI_Barrier(CommEnv_->MPI_GetComm());

      SYMPACK_TIMER_STOP(SPARSE_FWD_SUBST);

      //DumpContrib();

      //Back-substitution phase
      SYMPACK_TIMER_START(SPARSE_BACK_SUBST);

      //start from the root of the tree
      iLocalI = LocalSupernodes_.size() ;
      while(iLocalI>0|| !ContribsToSend.empty() || !outgoingSend.empty()){

        //Check for completion of outgoing communication
        AdvanceOutgoing(outgoingSend);

        if(iLocalI>0 && iLocalI<=LocalSupernodes_.size()){
          SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];
          I = cur_snode->Id();

          SuperNode<T,MallocAllocator> * contrib = Contributions_[iLocalI-1];

          Int parent = ETree_.PostParent(cur_snode->LastCol()-1);

          std::vector<Int> isBlockUpdated(contrib->NZBlockCnt(),0);

          if(parent!=0){
            Int parent_snode_id = SupMembership_[parent-1];
            Int iTarget = this->Mapping_->Map(parent_snode_id-1,parent_snode_id-1);
            //Do all my updates (Local and remote)
            //Local updates
            SuperNode<T,MallocAllocator> * dist_contrib;
            if(!LocalUpdates[iLocalI-1].empty()){
              Int contrib_snode_id = LocalUpdates[iLocalI-1].top();
              LocalUpdates[iLocalI-1].pop();
              dist_contrib = snodeLocal(contrib_snode_id,Contributions_);
              contrib->back_update(dist_contrib);
              //this->back_update(dist_contrib,contrib);
            }
            else{
              //Receive parent contrib
              MPI_Status recv_status;
              Int bytes_received = 0;
              MPI_Probe(iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);
              MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
              src_blocks.resize(bytes_received);

              MPI_Recv(&src_blocks[0],bytes_received,MPI_BYTE,iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);

              dist_contrib = CreateSuperNode<MallocAllocator>(options_.decomposition,&src_blocks[0],bytes_received);
              //TODO Replace this
              //Deserialize(&src_blocks[0],*dist_contrib); 
              dist_contrib->InitIdxToBlk();


#ifdef _DEBUG_
              logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib->Id()<<std::endl;
#endif
              contrib->back_update(dist_contrib);
              //this->back_update(dist_contrib,contrib);
              delete dist_contrib;
            }
          }

          //now compute MY contribution
          contrib->back_update_contrib(cur_snode);

          //send to my children
          Int colIdx = cur_snode->FirstCol()-1;
          if(colIdx>0){
            Int children_found = 0;
            while(children_found<children[I-1]){
              Int child_snode_id = SupMembership_[colIdx-1];

              Int parent = ETree_.PostParent(colIdx-1);
              if(parent!=0){
                if(SupMembership_[parent-1]==cur_snode->Id()){
                  Int iTarget = this->Mapping_->Map(child_snode_id-1,child_snode_id-1);

                  if(iTarget!=iam){

                    bool isSkipped= false;

                    Int next_local_contrib = (iLocalI >1)?Contributions_[iLocalI-2]->Id():0;
                    if(next_local_contrib > child_snode_id || true){
                      //need to push the prev src_last_row
                      //Send
                      Int src_nzblk_idx = 0;
                      NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);
                      Int src_first_row = pivot_desc.GIndex;
                      ContribsToSendDown.push(DelayedComm(contrib->Id(),child_snode_id,0,src_first_row));
#ifdef _DEBUG_DELAY_
                      std::cout<<"P"<<iam<<" has delayed update from Contrib "<<I<<" to "<<child_snode_id<<" from row "<<src_first_row<<std::endl;
#endif
                      isSkipped= true;
                    }


                    if(!isSkipped){
#ifdef _DEBUG_
                      logfileptr->OFS()<<"Remote Supernode "<<child_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif

                      Int src_nzblk_idx = 0;
                      NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);
                      Int src_first_row = pivot_desc.GIndex;
                      //Create a new Icomm buffer, serialize the contribution
                      // in it and add it to the outgoing comm list
                      Icomm * send_buffer = new Icomm();
                      //TODO replace this
                      contrib->Serialize(*send_buffer,src_nzblk_idx,src_first_row);
                      AddOutgoingComm(outgoingSend,send_buffer);


                      if( outgoingSend.size() > maxIsend_){
                        MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm());
                        outgoingSend.pop_back();
                      }
                      else{
                        MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm(),&outgoingSend.back()->Request);
                      }

#ifdef _DEBUG_            
                      logfileptr->OFS()<<"     Send contribution "<<I<<" to Supernode "<<child_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
                      //                    logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
                      //                    logfileptr->OFS()<<"Sending "<<nzblk_cnt*sizeof(NZBlockDesc2)<<" and "<<nz_cnt*sizeof(T)<<" bytes during BS"<<std::endl;
#endif
                    }
                  }
                  else{

#ifdef _DEBUG_
                    logfileptr->OFS()<<"Local Supernode "<<child_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif
                    Int iLocalJ = snodeLocalIndex(child_snode_id);
                    LocalUpdates[iLocalJ-1].push((Int)I);
                  }
                  children_found++;
                }
              }
              //last column of the prev supernode
              colIdx = Xsuper_[child_snode_id-1]-1;
              if(colIdx==0){
                break;
              }
            }


          }
        }

        SendDelayedMessagesDown(iLocalI,ContribsToSendDown,outgoingSend,Contributions_);
        --iLocalI;
      }
      SYMPACK_TIMER_STOP(SPARSE_BACK_SUBST);

      MPI_Barrier(CommEnv_->MPI_GetComm());

    }
    //DumpContrib();
    SYMPACK_TIMER_STOP(SPARSE_SOLVE);

  }

  template<typename T> void symPACKMatrix<T>::GetSolution(T * B, int nrhs){
    Int n = iSize_;
    //Int iam = CommEnv_->MPI_Rank();
    //Int np  = CommEnv_->MPI_Size();

    //if(iam<np)
    {
      std::fill(B,B+n*nrhs,T(0.0));

      //    Int nrhs = B.n();
      //Gather B from everybody and put it in the original matrix order
      std::vector<T> tmp_nzval;
      for(Int I=1; I<Xsuper_.size();++I){
        Int iOwner = this->Mapping_->Map(I-1,I-1);
        T * data;
        Int snode_size = Xsuper_[I] - Xsuper_[I-1];
        Int nzcnt = snode_size * nrhs;

        if( iOwner == iam ){
          SuperNode<T,MallocAllocator> * contrib = snodeLocal(I,Contributions_);
          //logfileptr->OFS()<<*contrib<<std::endl;
          data = contrib->GetNZval(0);
        }
        else{
          tmp_nzval.resize(nzcnt);
          data = &tmp_nzval[0];
        }

        //MPI_Bcast(data,nzcnt*sizeof(T),MPI_BYTE,iOwner,CommEnv_->MPI_GetComm());

        if(iam==iOwner){
          for(Int i = 0; i<snode_size; ++i){ 
            for(Int j = 0; j<nrhs; ++j){
              Int destRow = Xsuper_[I-1] + i;
              destRow = Order_.perm[destRow - 1];
              //B(destRow-1, j) = data[i*nrhs + j];
              B[destRow-1 + j*n] = data[i*nrhs + j];
            }
          }
        }
      }

      mpi::Allreduce((T*)MPI_IN_PLACE,&B[0],n*nrhs,MPI_SUM,fullcomm_);
      MPI_Barrier(CommEnv_->MPI_GetComm());
    }
  }


  template <typename T>
    template< class Alloc>
    void symPACKMatrix<T>::SendDelayedMessagesUp(Int iLocalI, CommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode<T,Alloc> *> & snodeColl){
      if(snodeColl.empty() || MsgToSend.empty()) { return;}

      //Index of the last global snode to do
      Int last_snode_id = Xsuper_.size()-1;
      //Index of the last local supernode
      Int last_local_id = snodeColl.back()->Id();
      //Index of the last PROCESSED supernode
      Int prev_snode_id = iLocalI<=snodeColl.size()?snodeColl[iLocalI-1]->Id():last_local_id;
      //Index of the next local supernode
      Int next_snode_id = prev_snode_id>=last_local_id?last_snode_id+1:snodeColl[iLocalI]->Id();

      bool is_last = prev_snode_id>=last_local_id;

      while( MsgToSend.size()>0){
        //Pull the highest priority message
        const DelayedComm & comm = MsgToSend.top();
        Int src_snode_id = comm.src_snode_id;
        Int tgt_snode_id = comm.tgt_snode_id;
        Int src_nzblk_idx = comm.src_nzblk_idx;
        Int src_first_row = comm.src_first_row;

#ifdef _DEBUG_DELAY_
        logfileptr->OFS()<<"Picked { "<<src_snode_id<<" -> "<<tgt_snode_id<<" }"<<std::endl;
#endif

        if(tgt_snode_id < next_snode_id || is_last /*|| OutgoingSend.size() <= maxIsend_*/){
          SendMessage(comm, OutgoingSend, snodeColl);
          //remove from the list
          MsgToSend.pop();
        }
        else{
          break;
        }
      }
    }


  template <typename T>
    template< class Alloc>
    void symPACKMatrix<T>::SendDelayedMessagesDown(Int iLocalI, DownCommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode<T,Alloc> *> & snodeColl){
      if(snodeColl.empty() || MsgToSend.empty()) { return;}

      //Index of the first local supernode
      Int first_local_id = snodeColl.front()->Id();
      //Index of the last PROCESSED supernode
      Int prev_snode_id = iLocalI>=1?snodeColl[iLocalI-1]->Id():first_local_id;
      //Index of the next local supernode
      Int next_snode_id = iLocalI<=1?0:snodeColl[iLocalI-2]->Id();

      bool is_last = prev_snode_id<=1;

      while( MsgToSend.size()>0){
        //Pull the highest priority message
        const DelayedComm & comm = MsgToSend.top();

        Int src_snode_id = comm.src_snode_id;
        Int tgt_snode_id = comm.tgt_snode_id;
        Int src_nzblk_idx = comm.src_nzblk_idx;
        Int src_first_row = comm.src_first_row;

#ifdef _DEBUG_DELAY_
        logfileptr->OFS()<<"Picked { "<<src_snode_id<<" -> "<<tgt_snode_id<<" }"<<std::endl;
#endif

        if(tgt_snode_id>next_snode_id || is_last /*|| OutgoingSend.size() <= maxIsend_*/){

          SuperNode<T,Alloc> * prev_src_snode = snodeLocal(src_snode_id,snodeColl);
          //this can be sent now

          Int iTarget = this->Mapping_->Map(tgt_snode_id-1,tgt_snode_id-1);
          if(iTarget != iam){
#ifdef _DEBUG_DELAY_
            logfileptr->OFS()<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode->Id()<<" to Supernode "<<tgt_snode_id<<std::endl;
            std::cout<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode->Id()<<" to Supernode "<<tgt_snode_id<<std::endl;
#endif

#ifdef _DEBUG_
            logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated by Supernode "<<prev_src_snode->Id()<<" rows "<<src_first_row/*<<" to "<<src_last_row*/<<" "<<src_nzblk_idx<<std::endl;
#endif

            NZBlockDesc & pivot_desc = prev_src_snode->GetNZBlockDesc(src_nzblk_idx);
            //Create a new Icomm buffer, serialize the contribution
            // in it and add it to the outgoing comm list
            Icomm * send_buffer = new Icomm();
            //TODO replace this
            prev_src_snode->Serialize(*send_buffer,src_nzblk_idx,src_first_row);
            AddOutgoingComm(OutgoingSend,send_buffer);


            if( OutgoingSend.size() > maxIsend_){
              MPI_Send(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm());
              OutgoingSend.pop_back();
            }
            else{
              MPI_Isend(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm(),&OutgoingSend.back()->Request);
            }
          }

          //remove from the list
          MsgToSend.pop();
        }
        else{
          break;
        }
      }
    }

  template <typename T>
    template< class Alloc>
    void symPACKMatrix<T>::SendMessage(const DelayedComm & comm, AsyncComms & OutgoingSend, std::vector<SuperNode<T,Alloc> *> & snodeColl){
      Int src_snode_id = comm.src_snode_id;
      Int tgt_snode_id = comm.tgt_snode_id;
      Int src_nzblk_idx = comm.src_nzblk_idx;
      Int src_first_row = comm.src_first_row;


      SuperNode<T,Alloc> * prev_src_snode = snodeLocal(src_snode_id,snodeColl);

      //this can be sent now
      Int iTarget = this->Mapping_->Map(tgt_snode_id-1,tgt_snode_id-1);
      if(iTarget != iam){
#ifdef _DEBUG_DELAY_
        logfileptr->OFS()<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode->Id()<<" to Supernode "<<tgt_snode_id<<std::endl;
        std::cout<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode->Id()<<" to Supernode "<<tgt_snode_id<<std::endl;
#endif

#ifdef _DEBUG_
        logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated by Supernode "<<prev_src_snode->Id()<<" rows "<<src_first_row/*<<" to "<<src_last_row*/<<" "<<src_nzblk_idx<<std::endl;
#endif
        NZBlockDesc & pivot_desc = prev_src_snode->GetNZBlockDesc(src_nzblk_idx);
        //Create a new Icomm buffer, serialize the contribution
        // in it and add it to the outgoing comm list
        Icomm * send_buffer = new Icomm();
        //TODO replace this
        prev_src_snode->Serialize(*send_buffer,src_nzblk_idx,src_first_row);
        AddOutgoingComm(OutgoingSend,send_buffer);


        if( OutgoingSend.size() > maxIsend_){
          MPI_Send(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm());
          OutgoingSend.pop_back();
        }
        else{
          MPI_Isend(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm(),&OutgoingSend.back()->Request);
        }

#ifdef _STAT_COMM_
        maxFactors_ = std::max(maxFactors_,send_buffer->size());
        sizesFactors_ += send_buffer->size();
        countFactors_++;
#endif


      }
    }




  template <typename T> void symPACKMatrix<T>::AdvanceOutgoing(AsyncComms & outgoingSend){
    scope_timer(a,ADVANCE_OUTGOING_COMM);
    //Check for completion of outgoing communication
    if(!outgoingSend.empty()){
      AsyncComms::iterator it = outgoingSend.begin();
      while(it != outgoingSend.end()){
        int flag = 0;
        int error_code = MPI_Test(&(*it)->Request,&flag,MPI_STATUS_IGNORE);
        if(flag){
          it = outgoingSend.erase(it);
        }
        else{
          it++;
        }
      }
    }
  }


  template<typename T> void symPACKMatrix<T>::AddOutgoingComm(AsyncComms & outgoingSend, Icomm * send_buffer){
    outgoingSend.push_back(send_buffer);
  }


  template <typename T> void symPACKMatrix<T>::GetUpdatingSupernodeCount(std::vector<Int> & sc,std::vector<Int> & mw, std::vector<Int> & mh, std::vector<Int> & numBlk){
    sc.resize(Xsuper_.size(),I_ZERO);
    std::vector<Int> marker(Xsuper_.size(),I_ZERO);
    mw.resize(Xsuper_.size(),I_ZERO);
    mh.resize(Xsuper_.size(),I_ZERO);
    numBlk.resize(Xsuper_.size(),I_ZERO);

    //Int numLocSnode = ( (Xsuper_.size()-1) / np);
    //Int firstSnode = iam*numLocSnode + 1;

    Int numLocSnode = XsuperDist_[iam+1]-XsuperDist_[iam];
    Int firstSnode = XsuperDist_[iam];

    for(Int locsupno = 1; locsupno<locXlindx_.size(); ++locsupno){
      Idx s = locsupno + firstSnode-1;

      Int first_col = Xsuper_[s-1];
      Int last_col = Xsuper_[s]-1;

      Ptr lfi = locXlindx_[locsupno-1];
      Ptr lli = locXlindx_[locsupno]-1;
      mh[s-1] = lli-lfi+1;



      //count number of contiguous blocks (at least one diagonal block)
      Idx iPrevRow = locLindx_[lfi-1]-1;
      Idx iFirstRow = locLindx_[lfi-1];

      Idx width = mh[s-1]; 


      Idx nzBlockCnt = 1;
      Idx prevSnode = -1;
      for(Ptr sidx = lfi; sidx<=lli;sidx++){

        Idx row = locLindx_[sidx-1];


        //enforce the first block to be a square diagonal block
        if(nzBlockCnt==1 && row>last_col){
          nzBlockCnt++;
        }
        else if(row!=iPrevRow+1){
          nzBlockCnt++;
        }
        iPrevRow=row;

        Int supno = SupMembership_[row-1];

        if(prevSnode!=supno){
          //Idx supno = locSupLindx_[sidx-1];
          if(marker[supno-1]!=s && supno!=s){
            marker[supno-1] = s;
            ++sc[supno-1];
            mw[supno-1] = std::max(mw[supno-1],last_col - first_col+1);
          }
        }
        prevSnode = supno;
      }
      numBlk[s-1] = nzBlockCnt;
    }

    //do an allreduce on sc, mw and mh
    //    MPI_Allreduce(MPI_IN_PLACE,&sc[0],sc.size(),MPI_INT,MPI_SUM,CommEnv_->MPI_GetComm());
    //    MPI_Allreduce(MPI_IN_PLACE,&mw[0],mw.size(),MPI_INT,MPI_MAX,CommEnv_->MPI_GetComm());
    //    MPI_Allreduce(MPI_IN_PLACE,&mh[0],mh.size(),MPI_INT,MPI_MAX,CommEnv_->MPI_GetComm());
    //    MPI_Allreduce(MPI_IN_PLACE,&numBlk[0],numBlk.size(),MPI_INT,MPI_SUM,CommEnv_->MPI_GetComm());

    MPI_Allreduce(MPI_IN_PLACE,&sc[0],sc.size(),MPI_INT,MPI_SUM,fullcomm_);
    MPI_Allreduce(MPI_IN_PLACE,&mw[0],mw.size(),MPI_INT,MPI_MAX,fullcomm_);
    MPI_Allreduce(MPI_IN_PLACE,&mh[0],mh.size(),MPI_INT,MPI_MAX,fullcomm_);
    MPI_Allreduce(MPI_IN_PLACE,&numBlk[0],numBlk.size(),MPI_INT,MPI_SUM,fullcomm_);
  }


  //  template <typename T> SparseMatrixStructure symPACKMatrix<T>::GetLocalStructure() const {
  //    return *Local_;
  //  }
  //
  //  template <typename T> SparseMatrixStructure symPACKMatrix<T>::GetGlobalStructure(){
  //    if(isGlobStructAllocated_){
  //      Global_=new SparseMatrixStructure();
  //      Local_->ToGlobal(*Global_,CommEnv_->MPI_GetComm());
  //      isGlobStructAllocated_= true;
  //    }
  //    return *Global_;
  //  }




  template<typename T>
    void symPACKMatrix<T>::DumpMatlab(){
      logfileptr->OFS()<<"+sparse([";
      for(Int I=1;I<Xsuper_.size();I++){
        Int src_first_col = Xsuper_[I-1];
        Int src_last_col = Xsuper_[I]-1;
        Int iOwner = this->Mapping_->Map(I-1,I-1);
        if( iOwner == iam ){
          SuperNode<T> * src_snode = snodeLocal(I);
          for(int blkidx=0;blkidx<src_snode->NZBlockCnt();++blkidx){
            NZBlockDesc & desc = src_snode->GetNZBlockDesc(blkidx);
            T * val = src_snode->GetNZval(desc.Offset);
            Int nRows = src_snode->NRows(blkidx);

            Int row = desc.GIndex;
            for(Int i = 0; i< nRows; ++i){
              for(Int j = 0; j< src_snode->Size(); ++j){
                logfileptr->OFS()<<row+i<<" ";
              }
            }
          }
        }
      }
      logfileptr->OFS()<<"],[";
      for(Int I=1;I<Xsuper_.size();I++){
        Int src_first_col = Xsuper_[I-1];
        Int src_last_col = Xsuper_[I]-1;
        Int iOwner = this->Mapping_->Map(I-1,I-1);
        if( iOwner == iam ){
          SuperNode<T> * src_snode = snodeLocal(I);
          for(int blkidx=0;blkidx<src_snode->NZBlockCnt();++blkidx){
            NZBlockDesc & desc = src_snode->GetNZBlockDesc(blkidx);
            T * val = src_snode->GetNZval(desc.Offset);
            Int nRows = src_snode->NRows(blkidx);

            Int row = desc.GIndex;
            for(Int i = 0; i< nRows; ++i){
              for(Int j = 0; j< src_snode->Size(); ++j){
                logfileptr->OFS()<<src_first_col+j<<" ";
              }
            }
          }
        }
      }
      logfileptr->OFS()<<"],[";
      logfileptr->OFS().precision(std::numeric_limits< T >::max_digits10);
      for(Int I=1;I<Xsuper_.size();I++){
        Int src_first_col = Xsuper_[I-1];
        Int src_last_col = Xsuper_[I]-1;
        Int iOwner = this->Mapping_->Map(I-1,I-1);
        if( iOwner == iam ){
          SuperNode<T> * src_snode = snodeLocal(I);
          for(int blkidx=0;blkidx<src_snode->NZBlockCnt();++blkidx){
            NZBlockDesc & desc = src_snode->GetNZBlockDesc(blkidx);
            T * val = src_snode->GetNZval(desc.Offset);
            Int nRows = src_snode->NRows(blkidx);

            Int row = desc.GIndex;
            for(Int i = 0; i< nRows; ++i){
              for(Int j = 0; j< src_snode->Size(); ++j){
                //                logfileptr->OFS()<<std::scientific<<val[i*src_snode->Size()+j]<<" ";
                logfileptr->OFS()<<std::scientific<<ToMatlabScalar(val[i*src_snode->Size()+j])<<" ";
              }
            }
          }
        }
      }
      logfileptr->OFS()<<"],"<<iSize_<<","<<iSize_<<")"<<std::endl;



    }










  template<typename T>
    void symPACKMatrix<T>::Dump(){
      for(Int I=1;I<Xsuper_.size();I++){
        Int src_first_col = Xsuper_[I-1];
        Int src_last_col = Xsuper_[I]-1;
        Int iOwner = this->Mapping_->Map(I-1,I-1);
        if( iOwner == iam ){
          SuperNode<T> * src_snode = snodeLocal(I);


          logfileptr->OFS()<<"+++++++++++++"<<I<<"++++++++"<<std::endl;

          logfileptr->OFS()<<"cols: ";
          for(Int i = 0; i< src_snode->Size(); ++i){
            logfileptr->OFS()<<" "<<Order_.perm[src_first_col+i-1];
          }
          logfileptr->OFS()<<std::endl;
          logfileptr->OFS()<<std::endl;
          for(int blkidx=0;blkidx<src_snode->NZBlockCnt();++blkidx){

            NZBlockDesc & desc = src_snode->GetNZBlockDesc(blkidx);
            T * val = src_snode->GetNZval(desc.Offset);
            Int nRows = src_snode->NRows(blkidx);

            Int row = desc.GIndex;
            for(Int i = 0; i< nRows; ++i){
              logfileptr->OFS()<<row+i<<" | "<<Order_.perm[row+i-1]<<":   ";
              for(Int j = 0; j< src_snode->Size(); ++j){
                logfileptr->OFS()<<val[i*src_snode->Size()+j]<<" ";
              }
              logfileptr->OFS()<<std::endl;
            }

            logfileptr->OFS()<<"_______________________________"<<std::endl;
          }
        }
      }

    }

  template<typename T>
    void symPACKMatrix<T>::DumpContrib(){
      for(Int I=1;I<Xsuper_.size();I++){
        Int src_first_col = Xsuper_[I-1];
        Int src_last_col = Xsuper_[I]-1;
        Int iOwner = this->Mapping_->Map(I-1,I-1);
        if( iOwner == iam ){
          SuperNode<T,MallocAllocator> * src_snode = snodeLocal(I,Contributions_);


          logfileptr->OFS()<<"+++++++++++++"<<I<<"++++++++"<<std::endl;

          logfileptr->OFS()<<"cols: ";
          for(Int i = 0; i< src_snode->Size(); ++i){
            logfileptr->OFS()<<" "<<Order_.perm[src_first_col+i-1];
          }
          logfileptr->OFS()<<std::endl;
          logfileptr->OFS()<<std::endl;
          for(int blkidx=0;blkidx<src_snode->NZBlockCnt();++blkidx){

            NZBlockDesc & desc = src_snode->GetNZBlockDesc(blkidx);
            T * val = src_snode->GetNZval(desc.Offset);
            Int nRows = src_snode->NRows(blkidx);

            Int row = desc.GIndex;
            for(Int i = 0; i< nRows; ++i){
              logfileptr->OFS()<<row+i<<" | "<<Order_.perm[row+i-1]<<":   ";
              for(Int j = 0; j< src_snode->Size(); ++j){
                logfileptr->OFS()<<val[i*src_snode->Size()+j]<<" ";
              }
              logfileptr->OFS()<<std::endl;
            }

            logfileptr->OFS()<<"_______________________________"<<std::endl;
          }
        }
      }

    }

  template <typename T> void symPACKMatrix<T>::Init(DistSparseMatrix<T> & pMat, symPACKOptions & options ){
    //
    //
    //#ifdef _STAT_COMM_
    //    maxAggreg_ =0;
    //    sizesAggreg_ =0;
    //    countAggreg_=0;
    //    maxFactors_ =0;
    //    sizesFactors_ =0;
    //    countFactors_=0;
    //
    //    maxAggregRecv_ =0;
    //    sizesAggregRecv_ =0;
    //    countAggregRecv_=0;
    //    maxFactorsRecv_ =0;
    //    sizesFactorsRecv_ =0;
    //    countFactorsRecv_=0;
    //#endif
    //
    //    SYMPACK_TIMER_START(SUPERMATRIX_INIT);
    //    //Create the CommEnvironment object if necessary
    //    if(CommEnv_!=NULL){
    //      delete CommEnv_;
    //    }
    //
    //    if(options.commEnv == NULL){
    //      throw std::runtime_error("The communication environment must be initialized in the options");
    //    }
    //
    //    options_ = options;
    //
    //    CommEnv_ = new CommEnvironment(*options.commEnv);
    //    iSize_ = pMat.size;
    //
    //    SYMPACK_TIMER_START(MAPPING);
    //
    //    //create the mapping
    //    Int pmapping = options.used_procs(np);
    //    Int pr = (Int)sqrt((double)pmapping);
    //    if(options.mappingTypeStr ==  "MODWRAP2DTREE"){
    //      this->Mapping_ = new ModWrap2DTreeMapping(pmapping, pr, pr, 1);
    //    }
    //    else if(options.mappingTypeStr ==  "WRAP2D"){
    //      this->Mapping_ = new Wrap2D(pmapping, pr, pr, 1);
    //    }
    //    else if(options.mappingTypeStr ==  "WRAP2DFORCED"){
    //      this->Mapping_ = new Wrap2DForced(pmapping, pr, pr, 1);
    //    }
    //    else if(options.mappingTypeStr ==  "MODWRAP2DNS"){
    //      this->Mapping_ = new Modwrap2DNS(pmapping, pr, pr, 1);
    //    }
    //    else if(options.mappingTypeStr ==  "ROW2D"){
    //      this->Mapping_ = new Row2D(pmapping, pmapping, pmapping, 1);
    //    }
    //    else if(options.mappingTypeStr ==  "COL2D"){
    //      this->Mapping_ = new Col2D(pmapping, pmapping, pmapping, 1);
    //    }
    //    else{
    //      this->Mapping_ = new Modwrap2D(pmapping, pr, pr, 1);
    //    }
    //
    //    SYMPACK_TIMER_STOP(MAPPING);
    //
    //    //Options
    //    maxIsend_ = options.maxIsend;
    //    maxIrecv_ = options.maxIrecv;
    //
    //    gMaxIrecv = maxIrecv_+maxIsend_;
    //
    //    SYMPACK_TIMER_START(SCHEDULER);
    //    switch(options_.scheduler){
    //      case DL:
    //        scheduler_ = new DLScheduler<std::list<FBTask>::iterator>();
    //        scheduler2_ = new DLScheduler<FBTask>();
    //        break;
    //      case MCT:
    //        scheduler_ = new MCTScheduler<std::list<FBTask>::iterator>();
    //        scheduler2_ = new MCTScheduler<FBTask>();
    //        break;
    //      case PR:
    //        scheduler_ = new PRScheduler<std::list<FBTask>::iterator>();
    //        scheduler2_ = new PRScheduler<FBTask>();
    //        break;
    //      case FIFO:
    //        scheduler_ = new FIFOScheduler<std::list<FBTask>::iterator>();
    //        scheduler2_ = new FIFOScheduler<FBTask>();
    //        break;
    //      default:
    //        scheduler_ = new DLScheduler<std::list<FBTask>::iterator>();
    //        scheduler2_ = new DLScheduler<FBTask>();
    //        break;
    //    }
    //    SYMPACK_TIMER_STOP(SCHEDULER);
    //
    //#ifdef EXPLICIT_PERMUTE
    //    //permute then expand is better but not doable for now
    //    //    graph_.SetComm(CommEnv_->MPI_GetComm());
    //
    //    //    logfileptr->OFS()<<"*********************************"<<std::endl;
    //    //    pMat.DumpMatlab();
    //    //    logfileptr->OFS()<<"*********************************"<<std::endl;
    //
    //    pMat.Localg_.SetSorted(0);
    //    pMat.ExpandSymmetric();
    //
    //    //    logfileptr->OFS()<<"*********************************"<<std::endl;
    //    //    pMat.DumpMatlab();
    //    //    logfileptr->OFS()<<"*********************************"<<std::endl;
    //
    //
    //    graph_ = pMat.GetLocalGraph();
    //    graph_.SetBaseval(1);
    //    graph_.SetKeepDiag(0);
    //    graph_.SetSorted(1);
    //    //graph_.ExpandSymmetric();
    //
    //#else
    //    Local_ = new SparseMatrixStructure();
    //    *Local_ = pMat.GetLocalStructure();
    //    graph_.SetComm(CommEnv_->MPI_GetComm());
    //    graph_.SetBaseval(1);
    //    graph_.SetKeepDiag(0);
    //    graph_.SetSorted(1);
    //    graph_.FromStructure(*Local_);
    //    graph_.ExpandSymmetric();
    //#endif    
    //
    //
    //    logfileptr->OFS()<<"Matrix structure expanded"<<std::endl;
    //
    //    SparseMatrixGraph * sgraph = NULL;
    //    {
    //      {
    //        //Create an Ordering object to hold the permutation
    //        Order_.SetCommEnvironment(CommEnv_);
    //
    //        double timeSta = get_time();
    //        SYMPACK_TIMER_START(ORDERING);
    //        switch(options_.ordering){
    //          case MMD:
    //            {
    //              sgraph = new SparseMatrixGraph();
    //              graph_.GatherStructure(*sgraph,0);
    //              sgraph->SetBaseval(1);
    //              sgraph->SetKeepDiag(0);
    //              Order_.MMD(*sgraph);
    //            }
    //            break;
    //
    //          case AMD:
    //            {
    //              sgraph = new SparseMatrixGraph();
    //              graph_.GatherStructure(*sgraph,0);
    //              sgraph->SetBaseval(1);
    //              sgraph->SetKeepDiag(0);
    //              Order_.AMD(*sgraph);
    //            }
    //            break;
    //
    //          case NDBOX:
    //            Order_.NDBOX(Size());
    //            break;
    //
    //          case NDGRID:
    //            Order_.NDGRID(Size());
    //            break;
    //
    //#ifdef USE_SCOTCH
    //          case SCOTCH:
    //            {
    //
    //              sgraph = new SparseMatrixGraph();
    //              graph_.GatherStructure(*sgraph,0);
    //              sgraph->SetKeepDiag(0);
    //              Order_.SCOTCH(*sgraph);
    //            }
    //            break;
    //#endif
    //#ifdef USE_METIS
    //          case METIS:
    //            {
    //              sgraph = new SparseMatrixGraph();
    //              graph_.GatherStructure(*sgraph,0);
    //              sgraph->SetKeepDiag(0);
    //              Order_.METIS(*sgraph);
    //            }
    //            break;
    //#endif
    //#ifdef USE_PARMETIS
    //          case PARMETIS:
    //            {
    //              Order_.PARMETIS(graph_);
    //            }
    //            break;
    //#endif
    //#ifdef USE_PTSCOTCH
    //          case PTSCOTCH:
    //            {
    //              Order_.PTSCOTCH(graph_);
    //            }
    //            break;
    //#endif
    //          case NATURAL:
    //            Order_.perm.resize(iSize_);
    //            for(int i=0;i<Order_.perm.size();i++){Order_.perm[i]=i+1;} 
    //            Order_.invp = Order_.perm;
    //          default:
    //
    //            //do nothing: either natural or user provided ordering
    //            break;
    //        }
    //        SYMPACK_TIMER_STOP(ORDERING);
    //
    //        double timeStop = get_time();
    //        if(iam==0){
    //          std::cout<<"Ordering time: "<<timeStop - timeSta<<std::endl;
    //        }
    //        logfileptr->OFS()<<"Ordering done"<<std::endl;
    //      }
    //    }
    //
    //
    //    std::vector<Int> cc,rc;
    //    {
    //
    //      double timeSta = get_time();
    //      //ETree_.ConstructETree(graph_,Order_);
    //      {
    //        if(sgraph==NULL){ 
    //          sgraph = new SparseMatrixGraph();
    //          graph_.GatherStructure(*sgraph,0);
    //        }
    //        //sgraph->SetKeepDiag(1);
    //        sgraph->SetBaseval(1);
    //        //        ETree_.ConstructETree(*sgraph,Order_,team_);
    //        ETree_.ConstructETree(*sgraph,Order_,CommEnv_->MPI_GetComm());
    //      }
    //
    //      ETree_.PostOrderTree(Order_);
    //      logfileptr->OFS()<<"ETREE done"<<std::endl;
    //      {
    //        double timeStart = get_time();
    //        if(sgraph==NULL){ 
    //          sgraph = new SparseMatrixGraph();
    //          graph_.GatherStructure(*sgraph,0);
    //        }
    //        //sgraph->SetKeepDiag(1);
    //        sgraph->SetBaseval(1);
    //        this->getLColRowCount(*sgraph,cc,rc);
    //        //bcast rc, cc
    //        double timeStop = get_time();
    //        if(iam==0){
    //          std::cout<<"Column count (gather + serial + bcast) construction time: "<<timeStop - timeStart<<std::endl;
    //        }
    //      }
    //
    //
    //      if(options_.ordering != NATURAL){
    //        ETree_.SortChildren(cc,Order_);
    //      }
    //
    //      double timeStop = get_time();
    //      if(iam==0){
    //        std::cout<<"Elimination tree construction time: "<<timeStop - timeSta<<std::endl;
    //      }
    //    }
    //
    //    if(sgraph!=NULL){delete sgraph;}
    //
    //
    //
    //#ifdef _DEBUG_
    //    logfileptr->OFS()<<"colcnt "<<cc<<std::endl;
    //    logfileptr->OFS()<<"rowcnt "<<rc<<std::endl;
    //#endif
    //    double flops = 0.0;
    //    for(Int i = 0; i<cc.size();++i){
    //      flops+= (double)pow((double)cc[i],2.0);
    //    }
    //
    //    if(iam==0){
    //      std::cout<<"Flops: "<<flops<<std::endl;
    //    }
    //
    //#if 1
    //    int64_t NNZ = 0;
    //    for(Int i = 0; i<cc.size();++i){
    //      NNZ+=cc[i];
    //    }
    //    if(iam==0){
    //      std::cout<<"NNZ in L factor: "<<NNZ<<std::endl;
    //    }
    //#endif
    //
    //#ifndef _NO_COMPUTATION_
    //    {
    //      double timeSta = get_time();
    //      this->findSupernodes(ETree_,Order_,cc,SupMembership_,Xsuper_,options.relax.maxSize);
    //      //logfileptr->OFS()<<"xsuper: "<<Xsuper_<<std::endl;
    //
    //      //Compute XsuperDist_
    //      {
    //        Idx supPerProc = (Xsuper_.size()-1) / np;
    //        XsuperDist_.resize(np+1,0);
    //        for(int p =0; p<np; p++){
    //          XsuperDist_[p]= p*supPerProc+1;
    //        }
    //        XsuperDist_[np] = Xsuper_.size();
    //      }
    //
    //      logfileptr->OFS()<<"Supernodes found"<<std::endl;
    //
    //      if(options_.relax.nrelax0>0){
    //        this->relaxSupernodes(ETree_, cc,SupMembership_, Xsuper_, options_.relax );
    //
    //        logfileptr->OFS()<<"Relaxation done"<<std::endl;
    //        //Refresh XsuperDist_
    //        {
    //          Idx supPerProc = (Xsuper_.size()-1) / np;
    //          XsuperDist_.resize(np+1,0);
    //          for(int p =0; p<np; p++){
    //            XsuperDist_[p]= p*supPerProc+1;
    //          }
    //          XsuperDist_[np] = Xsuper_.size();
    //        }
    //      }
    //
    //#ifdef EXPLICIT_PERMUTE
    //      {
    //        //gdb_lock(0);
    //        //recompute vertexDist based on XsuperDist
    //        std::vector<Idx> newVertexDist(np+1);
    //        for(int p = 0; p < np; p++){
    //          newVertexDist[p] = Xsuper_[XsuperDist_[p]-1];
    //          newVertexDist[p] = Xsuper_[XsuperDist_[p]-1] + (pMat.GetLocalGraph().GetBaseval()-1);
    //        }
    //        newVertexDist[np] = this->iSize_+1;
    //
    //        pMat.Permute(&Order_.invp[0],&newVertexDist[0]);
    //        graph_ = pMat.GetLocalGraph();
    //        graph_.SetBaseval(1);
    //        graph_.SetKeepDiag(1);
    //        graph_.SetSorted(1);
    //      }
    //#endif
    //
    //      {
    //        double timeSta = get_time();
    //        this->symbolicFactorizationRelaxedDist(cc);
    //
    //        double timeStop = get_time();
    //        if(iam==0){
    //          std::cout<<"Symbolic factorization time: "<<timeStop - timeSta<<std::endl;
    //        }
    //      }
    //
    //      logfileptr->OFS()<<"Symbfact done"<<std::endl;
    //      if(options.order_refinement_str == "SET")
    //      {
    //        this->refineSupernodes(2,1);
    //        if(iam==0){
    //          std::cout<<"Supernode reodering done"<<std::endl;
    //        }
    //
    //        if(options_.relax.nrelax0>0){
    //          this->relaxSupernodes(ETree_, cc,SupMembership_, Xsuper_, options_.relax );
    //          logfileptr->OFS()<<"Relaxation done"<<std::endl;
    //          //Refresh XsuperDist_
    //          {
    //            Idx supPerProc = (Xsuper_.size()-1) / np;
    //            XsuperDist_.resize(np+1,0);
    //            for(int p =0; p<np; p++){
    //              XsuperDist_[p]= p*supPerProc+1;
    //            }
    //            XsuperDist_[np] = Xsuper_.size();
    //          }
    //        }
    //
    //
    //        {
    //          //permute matrix again ?
    //
    //          double timeSta = get_time();
    //          this->symbolicFactorizationRelaxedDist(cc);
    //          double timeStop = get_time();
    //          if(iam==0){
    //            std::cout<<"Symbolic factorization time: "<<timeStop - timeSta<<std::endl;
    //          }
    //        }
    //        logfileptr->OFS()<<"Symbfact done"<<std::endl;
    //      }
    //
    //      double timeStop = get_time();
    //      if(iam==0){
    //        std::cout<<"Total symbolic factorization time: "<<timeStop - timeSta<<std::endl;
    //      }
    //    }
    //
    //#ifdef _DEBUG_
    //    logfileptr->OFS()<<"Membership list is "<<SupMembership_<<std::endl;
    //    logfileptr->OFS()<<"xsuper "<<Xsuper_<<std::endl;
    //#endif
    //
    //
    //#ifdef _OUTPUT_ETREE_
    //    logfileptr->OFS()<<"ETree is "<<ETree_<<std::endl;
    //    logfileptr->OFS()<<"Supernodal ETree is "<<ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_)<<std::endl;
    //#endif
    //
    //
    //    {
    //      scope_timer(a,LOAD_BALANCE);
    //      double timeSta = get_time();
    //
    //      SYMPACK_TIMER_START(LOAD_BALANCE);
    //      if(options.load_balance_str=="SUBCUBE-FI"){
    //        if(iam==0){ std::cout<<"Subtree to subcube FI mapping used"<<std::endl;}
    //        ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
    //        this->Balancer_ = new SubtreeToSubcube(np,SupETree,Xsuper_,XsuperDist_,SupMembership_,locXlindx_,locLindx_,cc,CommEnv_,true);
    //#ifdef _DEBUG_MAPPING_
    //        logfileptr->OFS()<<"Proc Mapping: "<<this->Balancer_->GetMap()<<std::endl;
    //        logfileptr->OFS()<<"SupETree: "<<SupETree<<std::endl;
    //#endif
    //
    //
    //        TreeMapping * test = dynamic_cast<TreeMapping*>(this->Mapping_);
    //        if(test==NULL){
    //          this->Mapping_->Update(this->Balancer_->GetMap());
    //        }
    //        else{
    //          test->Update((TreeLoadBalancer*)this->Balancer_);
    //        }
    //
    //      }
    //      else if(options.load_balance_str=="SUBCUBE-FO"){
    //        if(iam==0){ std::cout<<"Subtree to subcube FO mapping used"<<std::endl;}
    //        ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
    //        this->Balancer_ = new SubtreeToSubcube(np,SupETree,Xsuper_,XsuperDist_,SupMembership_,locXlindx_,locLindx_,cc,CommEnv_,false);
    //#ifdef _DEBUG_MAPPING_
    //        logfileptr->OFS()<<"Proc Mapping: "<<this->Balancer_->GetMap()<<std::endl;
    //        logfileptr->OFS()<<"SupETree: "<<SupETree<<std::endl;
    //#endif
    //
    //        TreeMapping * test = dynamic_cast<TreeMapping*>(this->Mapping_);
    //        if(test==NULL){
    //          this->Mapping_->Update(this->Balancer_->GetMap());
    //        }
    //        else{
    //          test->Update((TreeLoadBalancer*)this->Balancer_);
    //        }
    //
    //
    //      }
    //      else if(options.load_balance_str=="SUBCUBE-VOLUME-FI"){
    //        if(iam==0){ std::cout<<"Subtree to subcube volume FI mapping used"<<std::endl;}
    //        ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
    //        this->Balancer_ = new SubtreeToSubcubeVolume(np,SupETree,Xsuper_,XsuperDist_,SupMembership_,locXlindx_,locLindx_,cc,CommEnv_,true);
    //#ifdef _DEBUG_MAPPING_
    //        logfileptr->OFS()<<"Proc Mapping: "<<this->Balancer_->GetMap()<<std::endl;
    //        logfileptr->OFS()<<"SupETree: "<<SupETree<<std::endl;
    //#endif
    //
    //        TreeMapping * test = dynamic_cast<TreeMapping*>(this->Mapping_);
    //        if(test==NULL){
    //          this->Mapping_->Update(this->Balancer_->GetMap());
    //        }
    //        else{
    //          test->Update((TreeLoadBalancer*)this->Balancer_);
    //        }
    //
    //
    //      }
    //      else if(options.load_balance_str=="SUBCUBE-VOLUME-FO"){
    //        if(iam==0){ std::cout<<"Subtree to subcube volume FO mapping used"<<std::endl;}
    //        ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
    //        this->Balancer_ = new SubtreeToSubcubeVolume(np,SupETree,Xsuper_,XsuperDist_,SupMembership_,locXlindx_,locLindx_,cc,CommEnv_,false);
    //#ifdef _DEBUG_MAPPING_
    //        logfileptr->OFS()<<"Proc Mapping: "<<this->Balancer_->GetMap()<<std::endl;
    //        logfileptr->OFS()<<"SupETree: "<<SupETree<<std::endl;
    //#endif
    //
    //        TreeMapping * test = dynamic_cast<TreeMapping*>(this->Mapping_);
    //        if(test==NULL){
    //          this->Mapping_->Update(this->Balancer_->GetMap());
    //        }
    //        else{
    //          test->Update((TreeLoadBalancer*)this->Balancer_);
    //        }
    //
    //
    //      }
    //      else if(options.load_balance_str=="NNZ"){
    //        if(iam==0){ std::cout<<"Load Balancing on NNZ used"<<std::endl;}
    //        this->Balancer_ = new NNZBalancer(np,Xsuper_,cc);
    //#ifdef _DEBUG_MAPPING_
    //        logfileptr->OFS()<<"Proc Mapping: "<<this->Balancer_->GetMap()<<std::endl;
    //#endif
    //        this->Mapping_->Update(this->Balancer_->GetMap());
    //      }
    //      double timeStop = get_time();
    //      if(iam==0){
    //        std::cout<<"Load balancing time: "<<timeStop - timeSta<<std::endl;
    //      }
    //    }
    //
    //    SYMPACK_TIMER_START(Get_UpdateCount);
    //    GetUpdatingSupernodeCount(UpdateCount_,UpdateWidth_,UpdateHeight_,numBlk_);
    //    SYMPACK_TIMER_STOP(Get_UpdateCount);
    //
    //
    //
    //
    //    Int iam = CommEnv_->MPI_Rank();
    //    Int np  = CommEnv_->MPI_Size();
    //
    //    {
    //      double timeSta = get_time();
    //      scope_timer(b,DISTRIBUTING_MATRIX);
    //      //Icomm buffer;
    //      Icomm recv_buffer;
    //
    //
    //      AsyncComms incomingRecv;
    //      AsyncComms outgoingSend;
    //      Int numColFirst = std::max(1,iSize_ / np);
    //
    //      logfileptr->OFS()<<"Starting Send"<<std::endl;
    //
    //#ifndef NOTRY
    //      try
    //#endif
    //      {
    //#ifdef EXPLICIT_PERMUTE
    //        pMat.ToLowerTriangular();
    //        pMat.SortGraph();
    //
    //        //    logfileptr->OFS()<<"*********************************"<<std::endl;
    //        //    pMat.DumpMatlab();
    //        //    logfileptr->OFS()<<"*********************************"<<std::endl;
    //
    //        vector<std::pair<size_t,Icomm *> > send_map(np);
    //
    //        SYMPACK_TIMER_START(DISTRIBUTE_COUNTING);
    //
    //        Idx FirstLocalCol = pMat.Localg_.vertexDist[iam];
    //        Idx LastLocalCol = pMat.Localg_.vertexDist[iam+1];
    //        Int snodeCount = 0;
    //        for(Int I=1;I<Xsuper_.size();I++){
    //          Idx fc = Xsuper_[I-1];
    //          Idx lc = Xsuper_[I]-1;
    //          Int iWidth = lc-fc+1;
    //          Int iHeight = UpdateHeight_[I-1];
    //
    //          Int iDest = this->Mapping_->Map(I-1,I-1);
    //
    //          if(iDest==iam){
    //            ++snodeCount;
    //          }
    //
    //          //post all the recv and sends
    //          for(Idx col = fc;col<=lc;col++){
    //            //corresponding column in the unsorted matrix A
    //            if(col>= FirstLocalCol && col < LastLocalCol){
    //              Int nrows = 0;
    //              Idx local_col = col - FirstLocalCol+1;//1 based
    //              send_map[iDest].first += (pMat.Localg_.colptr[local_col]-pMat.Localg_.colptr[local_col-1])*(sizeof(Idx)+sizeof(T)) + sizeof(Idx) + sizeof(Ptr);
    //            }
    //          }
    //        }
    //        SYMPACK_TIMER_STOP(DISTRIBUTE_COUNTING);
    //
    //
    //        //Resize the local supernodes array
    //        LocalSupernodes_.reserve(snodeCount);
    //
    //        std::vector<upcxx::global_ptr<SuperNodeDesc > > localFactors;
    //        remoteFactors_.resize(Xsuper_.size()-1);
    //        std::fill((char*)&remoteFactors_[0],(char*)&remoteFactors_[0]+remoteFactors_.size()*sizeof(std::tuple<upcxx::global_ptr<SuperNodeDesc>,Int> ),0);
    //        //logfileptr->OFS()<<"My usable global memory size is: "<<upcxx::my_usable_global_memory_size()<<std::endl;
    //
    //        SYMPACK_TIMER_START(DISTRIBUTE_CREATE_SNODES);
    //        for(Int I=1;I<Xsuper_.size();I++){
    //          Int iDest = this->Mapping_->Map(I-1,I-1);
    //          //parse the first column to create the supernode structure
    //          if(iam==iDest){
    //            Int fc = Xsuper_[I-1];
    //            Int lc = Xsuper_[I]-1;
    //            Int iWidth = lc-fc+1;
    //            Int iHeight = UpdateHeight_[I-1];
    //            Int nzBlockCnt = numBlk_[I-1];
    //#ifndef ITREE2
    //            globToLocSnodes_.push_back(I-1);
    //#else 
    //            ITree::Interval snode_inter = { I, I, LocalSupernodes_.size() };
    //            globToLocSnodes_.Insert(snode_inter);
    //#endif
    //            SuperNode<T> * newSnode = CreateSuperNode(options_.decomposition,I,fc,lc,iHeight,iSize_,nzBlockCnt);
    //            LocalSupernodes_.push_back(newSnode);
    //          }
    //        }
    //        SYMPACK_TIMER_STOP(DISTRIBUTE_CREATE_SNODES);
    //
    //        //first create the structure of every supernode
    //        {
    //          scope_timer(a,DISTRIBUTE_CREATE_REMOTE_NODES);
    //          std::vector< int > superStructure(np);
    //          std::vector< int > rdisplsStructure;
    //          std::vector< int > sSuperStructure(np);
    //          std::vector< int > ssizes(np,0);
    //          std::vector< int > sdispls(np+1,0);
    //          std::vector< int > rsizes(np,0);
    //
    //          Int numLocSnode = XsuperDist_[iam+1]-XsuperDist_[iam];
    //          Int firstSnode = XsuperDist_[iam];
    //
    //          for(Int locsupno = 1; locsupno<locXlindx_.size(); ++locsupno){
    //            Idx I = locsupno + firstSnode-1;
    //            Int iDest = this->Mapping_->Map(I-1,I-1);
    //            ssizes[iDest] += 1 + 2*numBlk_[I-1]; //1 for supno index + numBlk startrows + numBlk number of rows
    //          }
    //
    //          sdispls[0] = 0;
    //          std::partial_sum(ssizes.begin(),ssizes.end(),&sdispls[1]);
    //
    //          rdisplsStructure = sdispls;
    //          sSuperStructure.resize(sdispls.back());
    //          for(Int locsupno = 1; locsupno<locXlindx_.size(); ++locsupno){
    //            Idx I = locsupno + firstSnode-1;
    //            Int fc = Xsuper_[I-1];
    //            Int lc = Xsuper_[I]-1;
    //            Int iWidth = lc - fc + 1;
    //            Ptr lfi = locXlindx_[locsupno-1];
    //            Ptr lli = locXlindx_[locsupno]-1;
    //            Int iDest = this->Mapping_->Map(I-1,I-1);
    //            int & tail = rdisplsStructure[iDest];
    //
    //            sSuperStructure[tail++] = I;
    //            //count number of contiguous rows
    //            for(Ptr sidx = lfi; sidx<=lli;sidx++){
    //              Idx iStartRow = locLindx_[sidx-1];
    //              Idx iPrevRow = iStartRow;
    //              Int iContiguousRows = 1;
    //              for(Int idx2 = sidx+1; idx2<=lli;idx2++){
    //                Idx iCurRow = locLindx_[idx2-1];
    //                if(iStartRow == locLindx_[lfi-1]){
    //                  if(iCurRow>iStartRow+iWidth-1){
    //                    //enforce the first block to be a square diagonal block
    //                    break;
    //                  }
    //                }
    //
    //                if(iCurRow==iPrevRow+1){
    //                  sidx++;
    //                  ++iContiguousRows;
    //                  iPrevRow=iCurRow;
    //                }
    //                else{
    //                  break;
    //                }
    //              }
    //
    //              sSuperStructure[tail++] = iStartRow;
    //              sSuperStructure[tail++] = iContiguousRows;
    //            }
    //          }
    //
    //          MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,pMat.comm);
    //          rdisplsStructure[0] = 0;
    //          std::partial_sum(rsizes.begin(),rsizes.end(),&rdisplsStructure[1]);
    //          superStructure.resize(rdisplsStructure.back());
    //
    //          //turn everything into byte sizes
    //          //      for(int p = 0; p<ssizes.size();p++){ ssizes[p]*=sizeof(int); }
    //          //      for(int p = 0; p<rsizes.size();p++){ rsizes[p]*=sizeof(int); }
    //          //      for(int p = 0; p<sdispls.size();p++){ sdispls[p]*=sizeof(int); }
    //          //      for(int p = 0; p<rdisplsStructure.size();p++){ rdisplsStructure[p]*=sizeof(int); }
    //
    //          //Do the alltoallv to get the structures        
    //          MPI_Alltoallv(&sSuperStructure[0], &ssizes[0], &sdispls[0], MPI_INT,
    //              &superStructure[0], &rsizes[0], &rdisplsStructure[0], MPI_INT,
    //              pMat.comm);
    //
    //          //loop through received structure and create supernodes
    //          for(Int p = 0; p<np; p++){
    //            int pos = rdisplsStructure[p];
    //            int end = rdisplsStructure[p+1];
    //            //        int pos = rdisplsStructure[p]/sizeof(int);
    //            //        int end = rdisplsStructure[p+1]/sizeof(int);
    //            while(pos<end){
    //              Int I = superStructure[pos++];
    //              Int nzBlockCnt = numBlk_[I-1];
    //
    //              Int fc = Xsuper_[I-1];
    //              Int lc = Xsuper_[I]-1;
    //              Int iWidth = lc-fc+1;
    //              SuperNode<T> * snode = snodeLocal(I);
    //              if(snode->NZBlockCnt()==0){
    //                for(Int i = 0; i<nzBlockCnt;i++){
    //                  Int iStartRow = superStructure[pos++];
    //                  Int iContiguousRows = superStructure[pos++];
    //                  snode->AddNZBlock(iContiguousRows , iWidth,iStartRow);
    //                }
    //                snode->Shrink();
    //              }
    //            }
    //          } 
    //        }
    //
    //        //check if the Alltoallv needs to be performed in multiple steps
    //        size_t total_send_size = 0;
    //        std::vector<size_t> stotcounts(np,0);
    //        for(Int i = 0; i< send_map.size(); i++){
    //          total_send_size += send_map[i].first;
    //          stotcounts[i] = send_map[i].first;
    //        }
    //
    //        Icomm * IsendPtr = new Icomm(total_send_size,MPI_REQUEST_NULL);
    //
    //        std::vector<size_t> spositions(np+1,0);
    //        spositions[0] = 0;
    //        std::partial_sum(stotcounts.begin(),stotcounts.end(),&spositions[1]);
    //
    //
    //        SYMPACK_TIMER_SPECIAL_START(serializing); 
    //        auto baseval = pMat.Localg_.GetBaseval();
    //        for(Int I=1;I<Xsuper_.size();I++){
    //          Idx fc = Xsuper_[I-1];
    //          Idx lc = Xsuper_[I]-1;
    //          Int iWidth = lc-fc+1;
    //          Int iHeight = UpdateHeight_[I-1];
    //
    //          Int iDest = this->Mapping_->Map(I-1,I-1);
    //
    //          //post all the recv and sends
    //          for(Idx col = fc;col<=lc;col++){
    //            //corresponding column in the unsorted matrix A
    //            if(col>= FirstLocalCol && col < LastLocalCol){
    //              Idx local_col = col - FirstLocalCol+1;//1 based
    //              Ptr colbeg = pMat.Localg_.colptr[local_col-1]-baseval;//0-based
    //              Ptr colend = pMat.Localg_.colptr[local_col]-baseval;//0-based
    //              Ptr nrows = colend-colbeg;
    //              IsendPtr->setHead(spositions[iDest]);
    //              (*IsendPtr)<<(Idx)col<<nrows;
    //              Serialize( *IsendPtr,  &pMat.Localg_.rowind[colbeg], nrows);
    //              Serialize( *IsendPtr,  &pMat.nzvalLocal[colbeg], nrows);
    //              spositions[iDest] = IsendPtr->head;
    //            }
    //          }
    //        }
    //        SYMPACK_TIMER_SPECIAL_STOP(serializing);
    //
    //
    //        Icomm * IrecvPtr = new Icomm(0,MPI_REQUEST_NULL);
    //
    //        spositions[0] = 0;
    //        std::partial_sum(stotcounts.begin(),stotcounts.end(),&spositions[1]);
    //        std::function<void(Icomm&,size_t)> resize_lambda =
    //          [](Icomm & container, size_t sz){
    //            container.resize(sz); 
    //            container.head = 0;
    //          };
    //
    //        IsendPtr->setHead(0);
    //        mpi::Alltoallv((*IsendPtr), &stotcounts[0], &spositions[0], MPI_BYTE,
    //            (*IrecvPtr),CommEnv_->MPI_GetComm(), resize_lambda);
    //
    //
    //        //    //do the all to all
    //        //    std::vector<size_t> rtotcounts(np,0);
    //        //    MPI_Alltoall(&stotcounts[0],sizeof(size_t),MPI_BYTE,&rtotcounts[0],sizeof(size_t),MPI_BYTE,CommEnv_->MPI_GetComm());
    //        //
    //        //
    //        //
    //        //
    //        //
    //        //    int intmax = std::numeric_limits<int>::max();
    //        //    int split_atav = 0;
    //        //    int chunk_atav = 0;
    //        //    for(Int i = 0; i< rtotcounts.size(); i++){
    //        //        total_recv_size += rtotcounts[i];
    //        //        //for the recv
    //        //        int cur_split_atav = (int)std::ceil((double)rtotcounts[i]/(double)intmax);
    //        //        if(split_atav<cur_split_atav){
    //        //          split_atav=cur_split_atav;
    //        //          chunk_atav=(int)std::floor((double)rtotcounts[i]/(double)split_atav);
    //        //        }
    //        //        //for the send
    //        //        cur_split_atav = (int)std::ceil((double)stotcounts[i]/(double)intmax);
    //        //        if(split_atav<cur_split_atav){
    //        //          split_atav=cur_split_atav;
    //        //          chunk_atav=(int)std::floor((double)stotcounts[i]/(double)split_atav);
    //        //        }
    //        //    }
    //        //
    //        //    Icomm * IrecvPtr = new Icomm(total_recv_size,MPI_REQUEST_NULL);
    //        //    std::vector<int> sdispls(np+1,0);
    //        //    std::vector<int> scounts(np,0);
    //        //    std::vector<int> rdispls(np+1,0);
    //        //    std::vector<int> rcounts(np,0);
    //        //    //now we know how many alltoallv are required
    //        //    for(int step_atav = 0; step_atav<split_atav;step_atav++){
    //        //      size_t offset_buffer = step_atav*chunk_atav;
    //        //      //TODO finish this
    //        //        
    //        //
    //        //      for(Int i = 0; i< scounts.size(); i++){
    //        //        scounts[i] = stotcounts[i];
    //        //      }
    //        //      for(Int i = 0; i< rcounts.size(); i++){
    //        //        rcounts[i] = rtotcounts[i];
    //        //      }
    //        //
    //        //    //compute send displacements
    //        //    sdispls[0] = 0;
    //        //    std::partial_sum(scounts.begin(),scounts.end(),&sdispls[1]);
    //        //    total_send_size = sdispls.back();
    //        //
    //        //    //compute receive displacements
    //        //    rdispls[0] = 0;
    //        //    std::partial_sum(rcounts.begin(),rcounts.end(),&rdispls[1]);
    //        //    total_recv_size = rdispls.back();
    //        //
    //        //
    //        //    MPI_Alltoallv(IsendPtr->front(), &scounts[0], &sdispls[0], MPI_BYTE,
    //        //        IrecvPtr->front(), &rcounts[0], &rdispls[0], MPI_BYTE,
    //        //        CommEnv_->MPI_GetComm());
    //        //
    //        //
    //        //    }
    //
    //        SYMPACK_TIMER_SPECIAL_START(deserializing);      
    //        IrecvPtr->setHead(0);
    //        while(IrecvPtr->head < IrecvPtr->capacity()){ 
    //          char * buffer = IrecvPtr->back();
    //          //Deserialize
    //          Idx col = *(Int*)&buffer[0];
    //          Int I = SupMembership_[col-1];
    //
    //          Int iDest = this->Mapping_->Map(I-1,I-1);
    //          bassert(iam==iDest);
    //
    //          Int fc = Xsuper_[I-1];
    //          Int lc = Xsuper_[I]-1;
    //          Int iWidth = lc-fc+1;
    //
    //          SuperNode<T> * snode = snodeLocal(I);
    //
    //          //nrows of column col sent by processor p
    //          Ptr nrows = *((Ptr*)&buffer[sizeof(Idx)]);
    //          Idx * rowind = (Idx*)(&buffer[sizeof(Ptr)+sizeof(Idx)]);
    //          T * nzvalA = (T*)(&buffer[nrows*sizeof(Idx)+sizeof(Ptr)+sizeof(Idx)]);
    //          //advance in the buffer
    //          IrecvPtr->setHead(IrecvPtr->head + sizeof(Idx) +sizeof(Ptr) + nrows*(sizeof(Idx)+sizeof(T)));
    //
    //          //Here, do a linear search instead for the blkidx
    //          Ptr colbeg = 1;
    //          Ptr colend = nrows;
    //
    //          if(colbeg<=colend){
    //            //sort rowind and nzvals
    //            //            std::vector<size_t> lperm = sort_permutation(&rowind[colbeg-1],&rowind[colend-1]+1,std::less<Int>());
    //            //            apply_permutation(&rowind[colbeg-1],&rowind[colend-1]+1,lperm);
    //            //            apply_permutation(&nzvalA[colbeg-1],&nzvalA[colend-1]+1,lperm);
    //            Idx firstRow = rowind[colbeg-1];
    //            Int blkidx = snode->FindBlockIdx(firstRow);
    //            assert(blkidx!=-1);
    //            Int blk_nrows = snode->NRows(blkidx);    
    //            NZBlockDesc * blk_desc = &snode->GetNZBlockDesc(blkidx);
    //
    //            for(Ptr rowidx = colbeg; rowidx<=colend; ++rowidx){
    //              Idx row = rowind[rowidx-1];
    //              while(row>blk_desc->GIndex+blk_nrows-1){
    //                blkidx++;
    //                blk_nrows = snode->NRows(blkidx);    
    //                blk_desc = &snode->GetNZBlockDesc(blkidx);
    //              }
    //
    //              Int local_row = row - blk_desc->GIndex + 1;
    //              Int local_col = col - fc + 1;
    //              T * nzval = snode->GetNZval(blk_desc->Offset);
    //              nzval[(local_row-1)*iWidth+local_col-1] = nzvalA[rowidx-1];
    //            }
    //          }
    //        }
    //        SYMPACK_TIMER_SPECIAL_STOP(deserializing);
    //
    //        //after the alltoallv, cleanup
    //        delete IrecvPtr;
    //        delete IsendPtr;
    //
    //
    //#else
    //        vector<std::pair<size_t,Icomm *> > send_map(np);
    //        for(Int p=0;p<np;++p){
    //          send_map[p].first = 0;
    //        }
    //        SYMPACK_TIMER_START(DISTRIBUTE_COUNTING);
    //
    //        Int snodeCount = 0;
    //        for(Int I=1;I<Xsuper_.size();I++){
    //          Idx fc = Xsuper_[I-1];
    //          Idx lc = Xsuper_[I]-1;
    //          Int iWidth = lc-fc+1;
    //          Int iHeight = UpdateHeight_[I-1];
    //
    //          Int iDest = this->Mapping_->Map(I-1,I-1);
    //
    //          if(iDest==iam){
    //            ++snodeCount;
    //          }
    //
    //          //look at the owner of the first column of the supernode
    //          Idx numColFirst = std::max(1,iSize_ / np);
    //
    //          //post all the recv and sends
    //          for(Idx col = fc;col<=lc;col++){
    //            //corresponding column in the unsorted matrix A
    //            Idx orig_col = Order_.perm[col-1];
    //            Idx iOwnerCol = std::min((orig_col-1)/numColFirst,(Idx)np-1);
    //            size_t & send_bytes = send_map[iDest].first;
    //            if(iam == iOwnerCol){
    //              Int nrows = 0;
    //              Idx local_col = (orig_col-(numColFirst)*iOwnerCol);
    //              for(Ptr rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
    //                Idx orig_row = pMat.Local_.rowind[rowidx-1];
    //                Idx row = Order_.invp[orig_row-1];
    //
    //                if(row<col){
    //                  //add the pair (col,row) to processor owning column row
    //                  Int J = SupMembership_[row-1];
    //                  Int iDestJ = this->Mapping_->Map(J-1,J-1);
    //                  send_map[iDestJ].first += sizeof(Idx)+sizeof(Int)+1*(sizeof(Idx)+sizeof(T));
    //                }
    //                else{nrows++;}
    //              }
    //              send_bytes += sizeof(Idx)+sizeof(Int)+nrows*(sizeof(Idx)+sizeof(T));
    //            }
    //          }
    //        }
    //
    //        SYMPACK_TIMER_STOP(DISTRIBUTE_COUNTING);
    //        //Resize the local supernodes array
    //        LocalSupernodes_.reserve(snodeCount);
    //
    //        std::vector<upcxx::global_ptr<SuperNodeDesc > > localFactors;
    //        remoteFactors_.resize(Xsuper_.size()-1);
    //        std::fill((char*)&remoteFactors_[0],(char*)&remoteFactors_[0]+remoteFactors_.size()*sizeof(std::tuple<upcxx::global_ptr<SuperNodeDesc>,Int> ),0);
    //        //logfileptr->OFS()<<"My usable global memory size is: "<<upcxx::my_usable_global_memory_size()<<std::endl;
    //
    //        SYMPACK_TIMER_START(DISTRIBUTE_CREATE_SNODES);
    //
    //        for(Int I=1;I<Xsuper_.size();I++){
    //
    //          Int iDest = this->Mapping_->Map(I-1,I-1);
    //
    //          //parse the first column to create the supernode structure
    //          if(iam==iDest){
    //            Int fc = Xsuper_[I-1];
    //            Int lc = Xsuper_[I]-1;
    //            Int iWidth = lc-fc+1;
    //            Int iHeight = UpdateHeight_[I-1];
    //            Int nzBlockCnt = numBlk_[I-1];
    //#ifndef ITREE2
    //            globToLocSnodes_.push_back(I-1);
    //#else 
    //            ITree::Interval snode_inter = { I, I, LocalSupernodes_.size() };
    //            globToLocSnodes_.Insert(snode_inter);
    //#endif
    //            SuperNode<T> * newSnode = CreateSuperNode(options_.decomposition,I,fc,lc,iHeight,iSize_,nzBlockCnt);
    //            LocalSupernodes_.push_back(newSnode);
    //          }
    //        }
    //        SYMPACK_TIMER_STOP(DISTRIBUTE_CREATE_SNODES);
    //
    //        //first create the structure of every supernode
    //        {
    //          std::vector< int > superStructure(np);
    //          std::vector< int > rdisplsStructure;
    //          std::vector< int > sSuperStructure(np);
    //          std::vector< int > ssizes(np,0);
    //          std::vector< int > sdispls(np+1,0);
    //          std::vector< int > rsizes(np,0);
    //
    //          Int numLocSnode = XsuperDist_[iam+1]-XsuperDist_[iam];
    //          Int firstSnode = XsuperDist_[iam];
    //
    //          for(Int locsupno = 1; locsupno<locXlindx_.size(); ++locsupno){
    //            Idx I = locsupno + firstSnode-1;
    //            Int iDest = this->Mapping_->Map(I-1,I-1);
    //            ssizes[iDest] += 1 + 2*numBlk_[I-1]; //1 for supno index + numBlk startrows + numBlk number of rows
    //          }
    //
    //          sdispls[0] = 0;
    //          std::partial_sum(ssizes.begin(),ssizes.end(),&sdispls[1]);
    //
    //          rdisplsStructure = sdispls;
    //          sSuperStructure.resize(sdispls.back());
    //          for(Int locsupno = 1; locsupno<locXlindx_.size(); ++locsupno){
    //            Idx I = locsupno + firstSnode-1;
    //            Int iDest = this->Mapping_->Map(I-1,I-1);
    //
    //            int & tail = rdisplsStructure[iDest];
    //
    //
    //            Int fc = Xsuper_[I-1];
    //            Int lc = Xsuper_[I]-1;
    //            Int iWidth = lc - fc + 1;
    //            Ptr lfi = locXlindx_[locsupno-1];
    //            Ptr lli = locXlindx_[locsupno]-1;
    //
    //            sSuperStructure[tail++] = I;
    //            //count number of contiguous rows
    //            for(Ptr sidx = lfi; sidx<=lli;sidx++){
    //              Idx iStartRow = locLindx_[sidx-1];
    //              Idx iPrevRow = iStartRow;
    //              Int iContiguousRows = 1;
    //              for(Int idx2 = sidx+1; idx2<=lli;idx2++){
    //                Idx iCurRow = locLindx_[idx2-1];
    //                if(iStartRow == locLindx_[lfi-1]){
    //                  if(iCurRow>iStartRow+iWidth-1){
    //                    //enforce the first block to be a square diagonal block
    //                    break;
    //                  }
    //                }
    //
    //                if(iCurRow==iPrevRow+1){
    //                  sidx++;
    //                  ++iContiguousRows;
    //                  iPrevRow=iCurRow;
    //                }
    //                else{
    //                  break;
    //                }
    //              }
    //
    //              sSuperStructure[tail++] = iStartRow;
    //              sSuperStructure[tail++] = iContiguousRows;
    //            }
    //          }
    //
    //          MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,pMat.comm);
    //          rdisplsStructure[0] = 0;
    //          std::partial_sum(rsizes.begin(),rsizes.end(),&rdisplsStructure[1]);
    //          superStructure.resize(rdisplsStructure.back());
    //
    //          //turn everything into byte sizes
    //          for(int p = 0; p<ssizes.size();p++){ ssizes[p]*=sizeof(int); }
    //          for(int p = 0; p<rsizes.size();p++){ rsizes[p]*=sizeof(int); }
    //          for(int p = 0; p<sdispls.size();p++){ sdispls[p]*=sizeof(int); }
    //          for(int p = 0; p<rdisplsStructure.size();p++){ rdisplsStructure[p]*=sizeof(int); }
    //
    //          //Do the alltoallv to get the structures        
    //          MPI_Alltoallv(&sSuperStructure[0], &ssizes[0], &sdispls[0], MPI_BYTE,
    //              &superStructure[0], &rsizes[0], &rdisplsStructure[0], MPI_BYTE,
    //              pMat.comm);
    //
    //          //loop through received structure and create supernodes
    //          for(Int p = 0; p<np; p++){
    //            int pos = rdisplsStructure[p]/sizeof(int);
    //            int end = rdisplsStructure[p+1]/sizeof(int);
    //            while(pos<end){
    //              Int I = superStructure[pos++];
    //              Int nzBlockCnt = numBlk_[I-1];
    //
    //              Int fc = Xsuper_[I-1];
    //              Int lc = Xsuper_[I]-1;
    //              Int iWidth = lc-fc+1;
    //              SuperNode<T> * snode = snodeLocal(I);
    //              if(snode->NZBlockCnt()==0){
    //                for(Int i = 0; i<nzBlockCnt;i++){
    //                  Int iStartRow = superStructure[pos++];
    //                  Int iContiguousRows = superStructure[pos++];
    //                  snode->AddNZBlock(iContiguousRows , iWidth,iStartRow);
    //                }
    //                snode->Shrink();
    //              }
    //            }
    //          } 
    //        }
    //
    //        {
    //          //allocate one buffer for every remote processor
    //          //compute the send structures
    //          size_t total_send_size = 0;
    //          std::vector<int> sdispls(np+1,0);
    //          std::vector<int> scounts(np,0);
    //          for(auto it = send_map.begin(); it!=send_map.end();it++){
    //            Int iCurDest = it->first;
    //            size_t & send_bytes = it->second.first;
    //            scounts[iCurDest] = send_bytes;
    //#ifdef _DEBUG_
    //            logfileptr->OFS()<<"P"<<iam<<" ---- "<<send_bytes<<" ---> P"<<iCurDest<<std::endl;
    //#endif
    //          }
    //          //compute send displacements
    //          sdispls[0] = 0;
    //          std::partial_sum(scounts.begin(),scounts.end(),&sdispls[1]);
    //          total_send_size = sdispls.back();
    //          //Fill up the send buffer
    //          SYMPACK_TIMER_SPECIAL_START(serializing);      
    //          for(Int I=1;I<Xsuper_.size();I++){
    //            Idx fc = Xsuper_[I-1];
    //            Idx lc = Xsuper_[I]-1;
    //            Int iWidth = lc-fc+1;
    //            Int iHeight = cc[fc-1];
    //
    //            Icomm & Isend = *IsendPtr;
    //            Int iDest = this->Mapping_->Map(I-1,I-1);
    //
    //#ifdef _DEBUG_
    //            logfileptr->OFS()<<"Supernode "<<I<<" is owned by P"<<iDest<<std::endl;
    //#endif
    //
    //            //look at the owner of the first column of the supernode
    //            Int numColFirst = std::max(1,iSize_ / np);
    //
    //            //Serialize
    //            for(Idx col = fc;col<=lc;col++){
    //              Idx orig_col = Order_.perm[col-1];
    //              Int iOwnerCol = std::min(Int((orig_col-1)/numColFirst),np-1);
    //              if(iam == iOwnerCol){
    //                Int nrows = 0;
    //                Int local_col = (orig_col-(numColFirst)*iOwnerCol);
    //                for(Ptr rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
    //                  Idx orig_row = pMat.Local_.rowind[rowidx-1];
    //                  Idx row = Order_.invp[orig_row-1];
    //
    //
    //                  if(row<col){
    //                    //add the pair (col,row) to processor owning column row
    //                    Int J = SupMembership_[row-1];
    //                    Int iDestJ = this->Mapping_->Map(J-1,J-1);
    //
    //                    T val = pMat.nzvalLocal[rowidx-1];
    //
    //                    //we need to set head to the proper sdispls
    //                    Isend.setHead(sposition[iDestJ]);
    //                    Isend<<row;
    //                    Isend<<1;
    //                    Isend<<col;
    //                    Isend<<val;
    //                    //backup new position for processor iDestJ
    //                    sposition[iDestJ] = Isend.head;
    //                  }
    //                  else{
    //                    nrows++;
    //                  }
    //                }
    //
    //
    //                //restore head
    //                Isend.setHead(sposition[iDest]);
    //                Isend<<col;
    //                Isend<<nrows;
    //
    //                for(Int rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
    //                  Int orig_row = pMat.Local_.rowind[rowidx-1];
    //                  Int row = Order_.invp[orig_row-1];
    //                  if(row>=col){
    //                    Isend<<row;
    //                  }
    //                }
    //
    //                for(Int rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
    //                  Int orig_row = pMat.Local_.rowind[rowidx-1];
    //                  Int row = Order_.invp[orig_row-1];
    //                  if(row>=col){
    //                    Isend<<pMat.nzvalLocal[rowidx-1];
    //                  }
    //                }
    //
    //                //backup new position for processor iDest
    //                sposition[iDest] = Isend.head;
    //
    //
    //              }
    //            }
    //          }
    //          SYMPACK_TIMER_SPECIAL_STOP(serializing);      
    //
    //          //compute the receive structures
    //          size_t total_recv_size = 0;
    //          std::vector<int> rdispls(np+1,0);
    //          std::vector<int> rcounts(np,0);
    //
    //          MPI_Alltoall(&scounts[0],sizeof(int),MPI_BYTE,&rcounts[0],sizeof(int),MPI_BYTE,pMat.comm);
    //
    //          //compute receive displacements
    //          rdispls[0] = 0;
    //          std::partial_sum(rcounts.begin(),rcounts.end(),&rdispls[1]);
    //          total_recv_size = rdispls.back();
    //
    //          Icomm * IrecvPtr = new Icomm(total_recv_size,MPI_REQUEST_NULL);
    //
    //          MPI_Alltoallv(IsendPtr->front(), &scounts[0], &sdispls[0], MPI_BYTE,
    //              IrecvPtr->front(), &rcounts[0], &rdispls[0], MPI_BYTE,
    //              pMat.comm);
    //
    //          //Need to parse the structure sent from the processor owning the first column of the supernode
    //
    //          SYMPACK_TIMER_SPECIAL_START(deserializing);      
    //          for(Int p=0;p<np;++p){
    //            IrecvPtr->setHead(rdispls[p]);
    //            while(IrecvPtr->head < rdispls[p]+rcounts[p]){ 
    //
    //              char * buffer = IrecvPtr->back();
    //              //Deserialize
    //              Idx col = *(Int*)&buffer[0];
    //              Int I = SupMembership_[col-1];
    //              //do a global to local mapping
    //
    //
    //              Int iDest = this->Mapping_->Map(I-1,I-1);
    //              bassert(iam==iDest);
    //
    //              Int fc = Xsuper_[I-1];
    //              Int lc = Xsuper_[I]-1;
    //              Int iWidth = lc-fc+1;
    //
    //              SuperNode<T> * snode = snodeLocal(I);
    //
    //              //nrows of column col sent by processor p
    //              Int nrows = *((Int*)&buffer[sizeof(Int)]);
    //              Idx * rowind = (Idx*)(&buffer[sizeof(Int)+sizeof(Idx)]);
    //              T * nzvalA = (T*)(&buffer[(1+nrows)*sizeof(Idx)+sizeof(Int)]);
    //              //advance in the buffer
    //              IrecvPtr->setHead(IrecvPtr->head + sizeof(Idx) +sizeof(Int) + nrows*(sizeof(Idx)+sizeof(T)));
    //
    //              //Here, do a linear search instead for the blkidx
    //
    //              Ptr colbeg = 1;
    //              Ptr colend = nrows;
    //              //            for(Ptr rowidx = colbeg; rowidx<=colend; ++rowidx){
    //              //              Idx row = rowind[rowidx-1];
    //              //              logfileptr->OFS()<<row<<" ";
    //              //            }
    //              //              logfileptr->OFS()<<std::endl;
    //
    //              if(colbeg<=colend){
    //                //sort rowind and nzvals
    //
    //#if 1 
    //                std::vector<size_t> lperm = sort_permutation(&rowind[colbeg-1],&rowind[colend-1]+1,std::less<Int>());
    //                apply_permutation(&rowind[colbeg-1],&rowind[colend-1]+1,lperm);
    //                apply_permutation(&nzvalA[colbeg-1],&nzvalA[colend-1]+1,lperm);
    //                Idx firstRow = rowind[colbeg-1];
    //                Int blkidx = snode->FindBlockIdx(firstRow);
    //                assert(blkidx!=-1);
    //                Int blk_nrows = snode->NRows(blkidx);    
    //                NZBlockDesc * blk_desc = &snode->GetNZBlockDesc(blkidx);
    //#else
    //                Idx prevRow = 0;
    //                Int blkidx = 0;
    //                NZBlockDesc * blk_desc = &snode->GetNZBlockDesc(blkidx);
    //                Int blk_nrows = snode->NRows(blkidx);
    //#endif
    //
    //                for(Ptr rowidx = colbeg; rowidx<=colend; ++rowidx){
    //                  Idx row = rowind[rowidx-1];
    //
    //#if 1
    //                  while(row>blk_desc->GIndex+blk_nrows-1){
    //                    blkidx++;
    //                    blk_nrows = snode->NRows(blkidx);    
    //                    blk_desc = &snode->GetNZBlockDesc(blkidx);
    //                  }
    //#else
    //                  if(row!=prevRow){
    //                    if(row>blk_desc->GIndex+blk_nrows || row<blk_desc->GIndex){
    //                      blkidx = snode->FindBlockIdx(row);
    //                      blk_desc = &snode->GetNZBlockDesc(blkidx);
    //                      blk_nrows = snode->NRows(blkidx);
    //                    }
    //                    prevRow=row;
    //                  }
    //#endif
    //
    //                  Int local_row = row - blk_desc->GIndex + 1;
    //                  Int local_col = col - fc + 1;
    //                  T * nzval = snode->GetNZval(blk_desc->Offset);
    //                  nzval[(local_row-1)*iWidth+local_col-1] = nzvalA[rowidx-1];
    //                }
    //                //for(Int rowidx = colbeg; rowidx<=colend; ++rowidx){
    //                //  Int row = rowind[rowidx-1];
    //                //  Int blkidx = snode->FindBlockIdx(row);
    //                //  assert(blkidx!=-1);
    //                //  NZBlockDesc & blk_desc = snode->GetNZBlockDesc(blkidx);
    //                //  Int local_row = row - blk_desc.GIndex + 1;
    //                //  Int local_col = col - fc + 1;
    //                //  T * nzval = snode->GetNZval(blk_desc.Offset);
    //                //  nzval[(local_row-1)*iWidth+local_col-1] = nzvalA[rowidx-1];
    //                //}
    //              }
    //            }
    //          }
    //          SYMPACK_TIMER_SPECIAL_STOP(deserializing);      
    //
    //          //after the alltoallv, cleanup
    //          delete IrecvPtr;
    //          delete IsendPtr;
    //        }
    //
    //#endif
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //        for(Int I=1;I<Xsuper_.size();I++){
    //          Int iDest = this->Mapping_->Map(I-1,I-1);
    //          //parse the first column to create the supernode structure
    //          if(iam==iDest){
    //            SuperNode<T> * newSnode = snodeLocal(I);
    //            SuperNodeDesc * meta = newSnode->GetMeta();
    //            remoteFactors_[I-1] = std::make_tuple( upcxx::global_ptr<SuperNodeDesc>( meta ), meta->blocks_cnt_) ;
    //          }
    //        }
    //        MPI_Allreduce( MPI_IN_PLACE, &remoteFactors_[0], remoteFactors_.size()*sizeof( std::tuple<upcxx::global_ptr<SuperNodeDesc>,Int> ), MPI_BYTE, MPI_BOR, CommEnv_->MPI_GetComm());
    //
    //
    //      }
    //#ifndef NOTRY
    //      catch(const std::bad_alloc& e){
    //        std::cout << "Allocation failed: " << e.what() << '\n';
    //        abort();
    //      }
    //#endif
    //      logfileptr->OFS()<<"Send Done"<<std::endl;
    //
    //      double timeStop = get_time();
    //      if(iam==0){
    //        std::cout<<"Distribution time: "<<timeStop - timeSta<<std::endl;
    //      }
    //    }
    //
    //#if 0
    //#ifdef _DEBUG_
    //    for(Int I=1;I<Xsuper_.size();I++){
    //      Int src_first_col = Xsuper_[I-1];
    //      Int src_last_col = Xsuper_[I]-1;
    //      Int iOwner = this->Mapping_->Map(I-1,I-1);
    //      //If I own the column, factor it
    //      if( iOwner == iam ){
    //        SuperNode<T> * src_snode = snodeLocal(I);
    //
    //
    //        logfileptr->OFS()<<"+++++++++++++"<<I<<"++++++++"<<std::endl;
    //
    //        logfileptr->OFS()<<"cols: ";
    //        for(Int i = 0; i< src_snode->Size(); ++i){
    //          logfileptr->OFS()<<" "<<Order_.perm[src_first_col+i-1];
    //        }
    //        logfileptr->OFS()<<std::endl;
    //        logfileptr->OFS()<<std::endl;
    //        for(int blkidx=0;blkidx<src_snode->NZBlockCnt();++blkidx){
    //
    //          NZBlockDesc & desc = src_snode->GetNZBlockDesc(blkidx);
    //          T * val = src_snode->GetNZval(desc.Offset);
    //          Int nRows = src_snode->NRows(blkidx);
    //
    //          Int row = desc.GIndex;
    //          for(Int i = 0; i< nRows; ++i){
    //            logfileptr->OFS()<<row+i<<" | "<<Order_.perm[row+i-1]<<":   ";
    //            for(Int j = 0; j< src_snode->Size(); ++j){
    //              logfileptr->OFS()<<val[i*src_snode->Size()+j]<<" ";
    //            }
    //            logfileptr->OFS()<<std::endl;
    //          }
    //
    //          logfileptr->OFS()<<"_______________________________"<<std::endl;
    //        }
    //      }
    //    }
    //#endif
    //#endif
    //    //Dump();
    //
    //#if 0
    //
    //    //allocate backup buffer
    //    {
    //      auto maxwit = max_element(&UpdateWidth_[0],&UpdateWidth_[0]+UpdateWidth_.size());
    //      auto maxhit = max_element(&UpdateHeight_[0],&UpdateHeight_[0]+UpdateHeight_.size());
    //
    //      int maxh = *maxhit;
    //      int maxw = *maxwit;
    //      for(Int I=1;I<Xsuper_.size();++I){
    //        int width = Xsuper_[I] - Xsuper_[I-1];
    //        maxw = max(maxw,width);
    //      }
    //
    //
    //      Int max_bytes = 6*sizeof(Int); 
    //      Int nrows = maxh;
    //      Int ncols = maxw;
    //      Int nz_cnt = nrows * ncols;
    //      Int nblocks = nrows;
    //      max_bytes += (nblocks)*sizeof(NZBlockDesc);
    //      max_bytes += nz_cnt*sizeof(T);
    //
    //      bool success = backupBuffer_.AllocLocal(max_bytes);
    //      assert(success);
    //
    //    }
    //
    //#endif
    //
    //
    //#ifdef PREALLOC_IRECV
    //#endif
    //
    //
    //
    //#endif
    //
    //
    //    {
    //      double timeSta = get_time();
    //      origTaskLists_.resize(Xsuper_.size(),NULL);
    //      std::vector<Int> AggregatesToRecv;
    //      std::vector<Int> LocalAggregates;
    //      FBGetUpdateCount(UpdatesToDo_,AggregatesToRecv,LocalAggregates);
    //
    //      generateTaskGraph(taskGraph_, AggregatesToRecv, LocalAggregates);
    //
    //      double timeStop = get_time();
    //      if(iam==0){
    //        std::cout<<"Task graph generation time: "<<timeStop - timeSta<<std::endl;
    //      }
    //    }
    //    MPI_Barrier(CommEnv_->MPI_GetComm());
    //
    //
    //    SYMPACK_TIMER_STOP(SUPERMATRIX_INIT);
    //
  }













  //#define BETTER_DISTRIB
  //#define BETTER_DISTRIB3
  //        //#define BETTER_DISTRIB2
  //
  //#ifndef BETTER_DISTRIB
  //        //first, count 
  //        map<Int,std::pair<size_t,Icomm *> > send_map;
  //#else
  //
  //#ifdef BETTER_DISTRIB3
  //
  //
  //        //map<Int, unordered_map<Int, std::set<triplet<T>, sortTriplet<T> > > > bufSup;
  //
  //        //   auto sortTriplet = [](const triplet<T> & a, const triplet<T> & b)->bool{
  //        //               bool retval = a.row<b.row;
  //        //               if(a.row==b.row){
  //        //               retval = a.col<b.col;
  //        //               }
  //        //               return retval;
  //        //               };
  //
  //
  //        map<Int, unordered_map<Int, std::priority_queue<triplet<T>, std::vector<triplet<T> >, sortTripletInv<T> > > > bufSup;
  //#endif
  //
  //#ifndef BETTER_DISTRIB2
  //        //map< proc , map< Sup, map<col,std::tuple<count,posRow,posNzval>  > > >
  //        map<Int, unordered_map<Int, unordered_map<Idx, std::tuple<int,int,int> > > > bufSupDesc;
  //#else
  //        Int NSuper = Xsuper_.size()-1;
  //        vector<Int> bufSupCnt(NSuper);
  //        vector<bool> bufSupMarker(NSuper,false);
  //        vector<size_t> bufSupPos(NSuper,0);
  //        //vector<Idx> lastLocalSup(np,0);
  //#endif
  //#endif
  //
  //
  //#ifndef BETTER_DISTRIB
  //        for(Int p=0;p<np;++p){
  //          send_map[p].first = 0;
  //        }
  //#endif
  //
  //        SYMPACK_TIMER_START(DISTRIBUTE_COUNTING);
  //
  //        Int snodeCount = 0;
  //        for(Int I=1;I<Xsuper_.size();I++){
  //          Idx fc = Xsuper_[I-1];
  //          Idx lc = Xsuper_[I]-1;
  //          Int iWidth = lc-fc+1;
  //          Int iHeight = UpdateHeight_[I-1];
  //
  //          Int iDest = this->Mapping_->Map(I-1,I-1);
  //
  //          if(iDest==iam){
  //            ++snodeCount;
  //          }
  //
  //          //look at the owner of the first column of the supernode
  //          Idx numColFirst = std::max(1,iSize_ / np);
  //
  //          //post all the recv and sends
  //          for(Idx col = fc;col<=lc;col++){
  //            //corresponding column in the unsorted matrix A
  //            Idx orig_col = Order_.perm[col-1];
  //            Idx iOwnerCol = std::min((orig_col-1)/numColFirst,(Idx)np-1);
  //#ifndef BETTER_DISTRIB
  //            size_t & send_bytes = send_map[iDest].first;
  //#endif
  //            if(iam == iOwnerCol){
  //              Int nrows = 0;
  //              Idx local_col = (orig_col-(numColFirst)*iOwnerCol);
  //              for(Ptr rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
  //                Idx orig_row = pMat.Local_.rowind[rowidx-1];
  //                Idx row = Order_.invp[orig_row-1];
  //
  //                if(row<col){
  //                  //add the pair (col,row) to processor owning column row
  //                  Int J = SupMembership_[row-1];
  //                  Int iDestJ = this->Mapping_->Map(J-1,J-1);
  //#ifndef BETTER_DISTRIB
  //                  send_map[iDestJ].first += sizeof(Idx)+sizeof(Int)+1*(sizeof(Idx)+sizeof(T));
  //#else
  //
  //#ifndef BETTER_DISTRIB2
  //                  auto & desc = bufSupDesc[iDestJ][J][row];
  //                  std::get<0>(desc)++;
  //#else
  //                  if(!bufSupMarker[J-1]){
  //                    bufSupMarker[J-1] = true;//++lastLocalSup[iDestJ];
  //                  }
  //                  bufSupCnt[J-1]++;
  //#endif
  //
  //#endif
  //                }
  //                else{
  //#ifndef BETTER_DISTRIB2
  //                  //add the pair (row,col) to iDest
  //                  nrows++;
  //#else
  //                  if(!bufSupMarker[I-1]){
  //                    bufSupMarker[I-1] = true;//++lastLocalSup[iDest];
  //                  }
  //                  bufSupCnt[I-1]++;
  //#endif
  //                }
  //              }
  //#ifndef BETTER_DISTRIB
  //              send_bytes += sizeof(Idx)+sizeof(Int)+nrows*(sizeof(Idx)+sizeof(T));
  //#else
  //
  //#ifndef BETTER_DISTRIB2
  //              auto & supMap = bufSupDesc[iDest][I];
  //              std::get<0>(supMap[col])+=nrows;
  //#endif
  //
  //#endif
  //            }
  //          }
  //        }
  //
  //        SYMPACK_TIMER_STOP(DISTRIBUTE_COUNTING);
  //
  //        //Resize the local supernodes array
  //        LocalSupernodes_.reserve(snodeCount);
  //
  //        //logfileptr->OFS()<<"INITIALIZING THE SHARED ARRAY"<<std::endl;
  //
  //        std::vector<upcxx::global_ptr<SuperNodeDesc > > localFactors;
  //        //localFactors.reserve(snodeCount);
  //        remoteFactors_.resize(Xsuper_.size()-1);
  //        std::fill((char*)&remoteFactors_[0],(char*)&remoteFactors_[0]+remoteFactors_.size()*sizeof(std::tuple<upcxx::global_ptr<SuperNodeDesc>,Int> ),0);
  //
  //        //logfileptr->OFS()<<"My usable global memory size is: "<<upcxx::my_usable_global_memory_size()<<std::endl;
  //
  //        SYMPACK_TIMER_START(DISTRIBUTE_CREATE_SNODES);
  //
  //        for(Int I=1;I<Xsuper_.size();I++){
  //
  //          Int iDest = this->Mapping_->Map(I-1,I-1);
  //
  //          //parse the first column to create the supernode structure
  //          if(iam==iDest){
  //            Int fc = Xsuper_[I-1];
  //            Int lc = Xsuper_[I]-1;
  //            Int iWidth = lc-fc+1;
  //            Int iHeight = UpdateHeight_[I-1];
  //            Int nzBlockCnt = numBlk_[I-1];
  //#ifndef ITREE2
  //            globToLocSnodes_.push_back(I-1);
  //#else 
  //            ITree::Interval snode_inter = { I, I, LocalSupernodes_.size() };
  //            globToLocSnodes_.Insert(snode_inter);
  //#endif
  //            SuperNode<T> * newSnode = CreateSuperNode(options_.decomposition,I,fc,lc,iHeight,iSize_,nzBlockCnt);
  //            LocalSupernodes_.push_back(newSnode);
  //          }
  //        }
  //
  //
  //        SYMPACK_TIMER_STOP(DISTRIBUTE_CREATE_SNODES);
  //
  //        //first create the structure of every supernode
  //        {
  //          std::vector< int > superStructure(np);
  //          std::vector< int > rdisplsStructure;
  //          std::vector< int > sSuperStructure(np);
  //          std::vector< int > ssizes(np,0);
  //          std::vector< int > sdispls(np+1,0);
  //          std::vector< int > rsizes(np,0);
  //
  //          Int numLocSnode = XsuperDist_[iam+1]-XsuperDist_[iam];
  //          Int firstSnode = XsuperDist_[iam];
  //
  //          for(Int locsupno = 1; locsupno<locXlindx_.size(); ++locsupno){
  //            Idx I = locsupno + firstSnode-1;
  //            Int iDest = this->Mapping_->Map(I-1,I-1);
  //            ssizes[iDest] += 1 + 2*numBlk_[I-1]; //1 for supno index + numBlk startrows + numBlk number of rows
  //          }
  //
  //          sdispls[0] = 0;
  //          std::partial_sum(ssizes.begin(),ssizes.end(),&sdispls[1]);
  //
  //          rdisplsStructure = sdispls;
  //          sSuperStructure.resize(sdispls.back());
  //          for(Int locsupno = 1; locsupno<locXlindx_.size(); ++locsupno){
  //            Idx I = locsupno + firstSnode-1;
  //            Int iDest = this->Mapping_->Map(I-1,I-1);
  //
  //            int & tail = rdisplsStructure[iDest];
  //
  //
  //            Int fc = Xsuper_[I-1];
  //            Int lc = Xsuper_[I]-1;
  //            Int iWidth = lc - fc + 1;
  //            Ptr lfi = locXlindx_[locsupno-1];
  //            Ptr lli = locXlindx_[locsupno]-1;
  //
  //            sSuperStructure[tail++] = I;
  //            //count number of contiguous rows
  //            for(Ptr sidx = lfi; sidx<=lli;sidx++){
  //              Idx iStartRow = locLindx_[sidx-1];
  //              Idx iPrevRow = iStartRow;
  //              Int iContiguousRows = 1;
  //              for(Int idx2 = sidx+1; idx2<=lli;idx2++){
  //                Idx iCurRow = locLindx_[idx2-1];
  //                if(iStartRow == locLindx_[lfi-1]){
  //                  if(iCurRow>iStartRow+iWidth-1){
  //                    //enforce the first block to be a square diagonal block
  //                    break;
  //                  }
  //                }
  //
  //                if(iCurRow==iPrevRow+1){
  //                  sidx++;
  //                  ++iContiguousRows;
  //                  iPrevRow=iCurRow;
  //                }
  //                else{
  //                  break;
  //                }
  //              }
  //
  //              sSuperStructure[tail++] = iStartRow;
  //              sSuperStructure[tail++] = iContiguousRows;
  //            }
  //          }
  //
  //          MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,pMat.comm);
  //          rdisplsStructure[0] = 0;
  //          std::partial_sum(rsizes.begin(),rsizes.end(),&rdisplsStructure[1]);
  //          superStructure.resize(rdisplsStructure.back());
  //
  //          //turn everything into byte sizes
  //          for(int p = 0; p<ssizes.size();p++){ ssizes[p]*=sizeof(int); }
  //          for(int p = 0; p<rsizes.size();p++){ rsizes[p]*=sizeof(int); }
  //          for(int p = 0; p<sdispls.size();p++){ sdispls[p]*=sizeof(int); }
  //          for(int p = 0; p<rdisplsStructure.size();p++){ rdisplsStructure[p]*=sizeof(int); }
  //
  //          //Do the alltoallv to get the structures        
  //          MPI_Alltoallv(&sSuperStructure[0], &ssizes[0], &sdispls[0], MPI_BYTE,
  //              &superStructure[0], &rsizes[0], &rdisplsStructure[0], MPI_BYTE,
  //              pMat.comm);
  //
  //          //loop through received structure and create supernodes
  //          for(Int p = 0; p<np; p++){
  //            int pos = rdisplsStructure[p]/sizeof(int);
  //            int end = rdisplsStructure[p+1]/sizeof(int);
  //            while(pos<end){
  //              Int I = superStructure[pos++];
  //              Int nzBlockCnt = numBlk_[I-1];
  //
  //              Int fc = Xsuper_[I-1];
  //              Int lc = Xsuper_[I]-1;
  //              Int iWidth = lc-fc+1;
  //              SuperNode<T> * snode = snodeLocal(I);
  //              if(snode->NZBlockCnt()==0){
  //                for(Int i = 0; i<nzBlockCnt;i++){
  //                  Int iStartRow = superStructure[pos++];
  //                  Int iContiguousRows = superStructure[pos++];
  //                  snode->AddNZBlock(iContiguousRows , iWidth,iStartRow);
  //                }
  //                snode->Shrink();
  //              }
  //            }
  //          } 
  //        }
  //
  //        {
  //          //allocate one buffer for every remote processor
  //          //compute the send structures
  //          size_t total_send_size = 0;
  //          std::vector<int> sdispls(np+1,0);
  //          std::vector<int> scounts(np,0);
  //#ifndef BETTER_DISTRIB
  //          for(auto it = send_map.begin(); it!=send_map.end();it++){
  //            Int iCurDest = it->first;
  //            size_t & send_bytes = it->second.first;
  //            scounts[iCurDest] = send_bytes;
  //#ifdef _DEBUG_
  //            logfileptr->OFS()<<"P"<<iam<<" ---- "<<send_bytes<<" ---> P"<<iCurDest<<std::endl;
  //#endif
  //          }
  //#else
  //
  //
  //#ifndef BETTER_DISTRIB2
  //          //changing from nrows to pos
  //          if(bufSupDesc.size()>0){
  //            scope_timer_special(a,BETTER_DISTRIB_CNT_TO_POS);
  //            auto prevProcIt = bufSupDesc.begin();
  //            for(Int p = 0 ; p< prevProcIt->first; ++p){
  //              scounts[p] = 0;
  //            }
  //
  //            for(auto procIt = bufSupDesc.begin(); procIt != bufSupDesc.end(); procIt++){
  //              //fill the holes
  //              for(Int p =prevProcIt->first+1; p<procIt->first; ++p){
  //                scounts[p] = 0;
  //              }
  //
  //              Int iCurDest = procIt->first;
  //              int send_bytes = 0;
  //              auto & supMap = procIt->second;
  //              for(auto supIt = supMap.begin() ; supIt != supMap.end(); supIt++){
  //                Int I = supIt->first;
  //                auto & colMap = supIt->second;
  //                for(auto colIt = colMap.begin() ; colIt != colMap.end(); colIt++){
  //                  Idx col = colIt->first;
  //                  //Int count = colIt->second;
  //                  auto & desc = colIt->second;
  //                  Int count = std::get<0>(desc);
  //                  //for each cols, add 1 x Idx, 1 Int,  count x Idx and count x T
  //                  int col_bytes = (count+1)*sizeof(Idx) + sizeof(Int) + count*sizeof(T);
  //                  send_bytes+= col_bytes;
  //
  //                  //update the buf pos (nzval offset will be computed later)
  //                  std::get<1>(desc) = total_send_size;
  //                  total_send_size+= col_bytes;
  //                }
  //              }
  //              scounts[iCurDest] = send_bytes;
  //              prevProcIt = procIt;
  //            }
  //
  //            for(Int p = bufSupDesc.rbegin()->first+1; p<np; ++p){
  //              scounts[p] = 0;
  //            }
  //          }
  //          else{
  //            for(Int p = 0; p<np; ++p){
  //              scounts[p] = 0;
  //            }
  //          }
  //          total_send_size = 0;
  //#else
  //
  //          for(Int p = 0; p<np; ++p){
  //            scounts[p] = 0;
  //          }
  //
  //
  //          for(Int I = 1; I<= NSuper;I++){        
  //            Int iDest = this->Mapping_->Map(I-1,I-1);
  //            auto nnzcnt = bufSupCnt[I-1];
  //            //supId + nnz + nnzcnt triplets
  //            //scounts[iDest]+=nnzcnt*(2*sizeof(Idx)+sizeof(T))+2*sizeof(Int);
  //            scounts[iDest]+=nnzcnt*(2*sizeof(Idx)+sizeof(T));
  //          }
  //
  //
  //          sdispls[0] = 0;
  //          std::partial_sum(scounts.begin(),scounts.end(),&sdispls[1]);
  //          total_send_size = sdispls.back();
  //
  //          {
  //            vector<int> spos(np,0); 
  //            for(Int I = 1; I<= NSuper;I++){        
  //              Int iDest = this->Mapping_->Map(I-1,I-1);
  //              auto nnzcnt = bufSupCnt[I-1];
  //              bufSupPos[I-1] = sdispls[iDest]+spos[iDest];
  //              spos[iDest]+=nnzcnt*(2*sizeof(Idx)+sizeof(T));
  //            }
  //          }
  //
  //
  //#endif
  //
  //          //      {
  //          //        scope_timer(a,NEW_OFFSET);
  //          //        
  //          //        for(Int p = 0; p<np; ++p){
  //          //            scounts[p] = 0;
  //          //        }
  //          //
  //          //        for(auto colIt = colOffset.begin() ; colIt != colOffset.end(); colIt++){
  //          //          Idx col = colIt->first;
  //          //          auto & triplet = colIt->second;
  //          //
  //          //          Int I = SupMembership_[I-1];
  //          //          Int iDest = this->Mapping_->Map(I-1,I-1);
  //          //
  //          //
  //          //          std::get<0>(triplet)+=nrows;
  //          //        }
  //          //      }
  //
  //
  //#endif
  //
  //#ifndef BETTER_DISTRIB2
  //          //compute send displacements
  //          sdispls[0] = 0;
  //          std::partial_sum(scounts.begin(),scounts.end(),&sdispls[1]);
  //          total_send_size = sdispls.back();
  //#endif
  //
  //          auto sposition = sdispls;
  //          Icomm * IsendPtr = new Icomm(total_send_size,MPI_REQUEST_NULL);
  //
  //          SYMPACK_TIMER_SPECIAL_START(serializing);      
  //#ifndef BETTER_DISTRIB
  //          //Fill up the send buffer
  //          for(Int I=1;I<Xsuper_.size();I++){
  //            Idx fc = Xsuper_[I-1];
  //            Idx lc = Xsuper_[I]-1;
  //            Int iWidth = lc-fc+1;
  //            Int iHeight = cc[fc-1];
  //
  //            Icomm & Isend = *IsendPtr;
  //            Int iDest = this->Mapping_->Map(I-1,I-1);
  //
  //#ifdef _DEBUG_
  //            logfileptr->OFS()<<"Supernode "<<I<<" is owned by P"<<iDest<<std::endl;
  //#endif
  //
  //            //look at the owner of the first column of the supernode
  //            Int numColFirst = std::max(1,iSize_ / np);
  //
  //            //Serialize
  //            for(Idx col = fc;col<=lc;col++){
  //              Idx orig_col = Order_.perm[col-1];
  //              Int iOwnerCol = std::min(Int((orig_col-1)/numColFirst),np-1);
  //              if(iam == iOwnerCol){
  //                Int nrows = 0;
  //                Int local_col = (orig_col-(numColFirst)*iOwnerCol);
  //                for(Ptr rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
  //                  Idx orig_row = pMat.Local_.rowind[rowidx-1];
  //                  Idx row = Order_.invp[orig_row-1];
  //
  //
  //                  if(row<col){
  //                    //add the pair (col,row) to processor owning column row
  //                    Int J = SupMembership_[row-1];
  //                    Int iDestJ = this->Mapping_->Map(J-1,J-1);
  //
  //                    T val = pMat.nzvalLocal[rowidx-1];
  //
  //                    //we need to set head to the proper sdispls
  //                    Isend.setHead(sposition[iDestJ]);
  //                    Isend<<row;
  //                    Isend<<1;
  //                    Isend<<col;
  //                    Isend<<val;
  //                    //backup new position for processor iDestJ
  //                    sposition[iDestJ] = Isend.head;
  //                  }
  //                  else{
  //                    nrows++;
  //                  }
  //                }
  //
  //
  //                //restore head
  //                Isend.setHead(sposition[iDest]);
  //                Isend<<col;
  //                Isend<<nrows;
  //
  //                for(Int rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
  //                  Int orig_row = pMat.Local_.rowind[rowidx-1];
  //                  Int row = Order_.invp[orig_row-1];
  //                  if(row>=col){
  //                    Isend<<row;
  //                  }
  //                }
  //
  //                for(Int rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
  //                  Int orig_row = pMat.Local_.rowind[rowidx-1];
  //                  Int row = Order_.invp[orig_row-1];
  //                  if(row>=col){
  //                    Isend<<pMat.nzvalLocal[rowidx-1];
  //                  }
  //                }
  //
  //                //backup new position for processor iDest
  //                sposition[iDest] = Isend.head;
  //
  //
  //              }
  //            }
  //          }
  //#else
  //
  //
  //#ifndef BETTER_DISTRIB2
  //          std::unordered_set<Idx> colInserted;
  //          Idx numColFirst = (std::max(1,iSize_ / np));
  //          for(Idx local_col = 1 ; local_col < pMat.Local_.colptr.size(); local_col++){
  //            Idx orig_col = iam * numColFirst + local_col;
  //            Idx col = Order_.invp[orig_col-1];
  //
  //            Icomm & Isend = *IsendPtr;
  //            Int I = SupMembership_[col-1];
  //            Int iDest = this->Mapping_->Map(I-1,I-1);
  //            auto & posSup=bufSupDesc[iDest][I];
  //
  //            volatile int found = 0;
  //            for(Ptr rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
  //              Idx orig_row = pMat.Local_.rowind[rowidx-1];
  //              Idx row = Order_.invp[orig_row-1];
  //              T & val = pMat.nzvalLocal[rowidx-1];
  //
  //
  //              if(row<col){
  //                //add the pair (col,row) to processor owning column row
  //                Int destCol = row;
  //                Int destRow = col;
  //
  //                Int J = SupMembership_[row-1];
  //                Int iDestJ = this->Mapping_->Map(J-1,J-1);
  //
  //                auto & desc = bufSupDesc[iDestJ][J][destCol];
  //                auto & bufPos = std::get<1>(desc);
  //                auto & nzvalPos = std::get<2>(desc);
  //                int cnt = colInserted.count(destCol);
  //                if(cnt==0){
  //                  auto nrows = std::get<0>(desc);
  //                  Isend.setHead(bufPos);
  //                  //serialize destCol and nrows
  //                  Isend << destCol;
  //                  Isend << nrows;
  //                  bufPos = Isend.head;
  //                  nzvalPos = bufPos+nrows*sizeof(Idx);
  //                  colInserted.insert(destCol);
  //                }
  //
  //                //write the row index
  //                Isend.setHead(bufPos);
  //                Isend<<destRow;
  //                bufPos+=sizeof(Idx);
  //                //write the nzval
  //                Isend.setHead(nzvalPos);
  //                Isend<<val;
  //                nzvalPos+=sizeof(T);
  //
  //
  //              }
  //              else{
  //                found = 1;
  //              }
  //            }
  //
  //            if(found){
  //              auto & desc = posSup[col];
  //              auto nrows = std::get<0>(desc);
  //              auto & bufPos = std::get<1>(desc);
  //              auto & nzvalPos = std::get<2>(desc);
  //              int cnt = colInserted.count(col);
  //              if(cnt==0){
  //                Isend.setHead(bufPos);
  //                //serialize destCol and nrows
  //                Isend << col;
  //                Isend << nrows;
  //                bufPos = Isend.head;
  //                nzvalPos = bufPos+nrows*sizeof(Idx);
  //                colInserted.insert(col);
  //              }
  //
  //              Isend.setHead(bufPos);
  //              for(Ptr rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
  //                Idx orig_row = pMat.Local_.rowind[rowidx-1];
  //                Idx row = Order_.invp[orig_row-1];
  //                if(row>=col){
  //                  Isend<<row;
  //                }
  //              }
  //              bufPos+=nrows*sizeof(Idx);
  //
  //              Isend.setHead(nzvalPos);
  //              for(Ptr rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
  //                Idx orig_row = pMat.Local_.rowind[rowidx-1];
  //                Idx row = Order_.invp[orig_row-1];
  //                if(row>=col){
  //                  T & val = pMat.nzvalLocal[rowidx-1];
  //                  Isend<<val;
  //                }
  //              }
  //              nzvalPos += nrows*sizeof(T);
  //            }
  //
  //
  //          }
  //#else
  //
  //          for(Int I=1;I<Xsuper_.size();I++){
  //            Idx fc = Xsuper_[I-1];
  //            Idx lc = Xsuper_[I]-1;
  //            Int iWidth = lc-fc+1;
  //            Int iHeight = UpdateHeight_[I-1];
  //
  //            Int iDest = this->Mapping_->Map(I-1,I-1);
  //
  //            Icomm & Isend = *IsendPtr;
  //            //look at the owner of the first column of the supernode
  //            Idx numColFirst = std::max(1,iSize_ / np);
  //
  //            //post all the recv and sends
  //            for(Idx col = fc;col<=lc;col++){
  //              //corresponding column in the unsorted matrix A
  //              Idx orig_col = Order_.perm[col-1];
  //              Idx iOwnerCol = std::min((orig_col-1)/numColFirst,(Idx)np-1);
  //              if(iam == iOwnerCol){
  //                Int nrows = 0;
  //                Idx local_col = (orig_col-(numColFirst)*iOwnerCol);
  //                for(Ptr rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
  //                  Idx orig_row = pMat.Local_.rowind[rowidx-1];
  //                  Idx row = Order_.invp[orig_row-1];
  //                  T & val = pMat.nzvalLocal[rowidx-1];
  //
  //                  Int destCol = col;
  //                  Int destRow = row;
  //                  if(row<col){
  //                    destCol = row;
  //                    destRow = col;
  //                  }
  //                  //add the pair (destCol,destRow) to processor owning column destCol
  //                  Int J = SupMembership_[destCol-1];
  //                  Int iDestJ = this->Mapping_->Map(J-1,J-1);
  //
  //                  Isend.setHead(bufSupPos[J-1]);
  //                  //            if(bufSupMarker[J-1]){
  //                  //              bufSupMarker[J-1]=false;
  //                  //              Isend<<J;
  //                  //              Isend<<bufSupCnt[J-1];
  //                  //            }
  //
  //                  Isend<<destRow<<destCol<<val;
  //                  bufSupPos[J-1] = Isend.head;
  //                }
  //              }
  //            }
  //          }
  //
  //
  //
  //
  //
  //
  //#endif
  //
  //#endif
  //          SYMPACK_TIMER_SPECIAL_STOP(serializing);      
  //
  //          //compute the receive structures
  //          size_t total_recv_size = 0;
  //          std::vector<int> rdispls(np+1,0);
  //          std::vector<int> rcounts(np,0);
  //
  //          MPI_Alltoall(&scounts[0],sizeof(int),MPI_BYTE,&rcounts[0],sizeof(int),MPI_BYTE,pMat.comm);
  //
  //          //compute receive displacements
  //          rdispls[0] = 0;
  //          std::partial_sum(rcounts.begin(),rcounts.end(),&rdispls[1]);
  //          total_recv_size = rdispls.back();
  //
  //          Icomm * IrecvPtr = new Icomm(total_recv_size,MPI_REQUEST_NULL);
  //
  //
  //#ifdef _DEBUG_
  //          logfileptr->OFS()<<"scounts: "<<scounts<<std::endl;
  //          logfileptr->OFS()<<"sdispls: "<<sdispls<<std::endl;
  //          logfileptr->OFS()<<"rcounts: "<<rcounts<<std::endl;
  //          logfileptr->OFS()<<"rdispls: "<<rdispls<<std::endl;
  //#endif
  //
  //          MPI_Alltoallv(IsendPtr->front(), &scounts[0], &sdispls[0], MPI_BYTE,
  //              IrecvPtr->front(), &rcounts[0], &rdispls[0], MPI_BYTE,
  //              pMat.comm);
  //
  //          //Need to parse the structure sent from the processor owning the first column of the supernode
  //
  //          SYMPACK_TIMER_SPECIAL_START(deserializing);      
  //#ifndef BETTER_DISTRIB2
  //          for(Int p=0;p<np;++p){
  //            IrecvPtr->setHead(rdispls[p]);
  //            while(IrecvPtr->head < rdispls[p]+rcounts[p]){ 
  //
  //              char * buffer = IrecvPtr->back();
  //              //Deserialize
  //              Idx col = *(Int*)&buffer[0];
  //              Int I = SupMembership_[col-1];
  //              //do a global to local mapping
  //
  //
  //              Int iDest = this->Mapping_->Map(I-1,I-1);
  //              bassert(iam==iDest);
  //
  //              Int fc = Xsuper_[I-1];
  //              Int lc = Xsuper_[I]-1;
  //              Int iWidth = lc-fc+1;
  //
  //              SuperNode<T> * snode = snodeLocal(I);
  //
  //              //nrows of column col sent by processor p
  //              Int nrows = *((Int*)&buffer[sizeof(Int)]);
  //              Idx * rowind = (Idx*)(&buffer[sizeof(Int)+sizeof(Idx)]);
  //              T * nzvalA = (T*)(&buffer[(1+nrows)*sizeof(Idx)+sizeof(Int)]);
  //              //advance in the buffer
  //              IrecvPtr->setHead(IrecvPtr->head + sizeof(Idx) +sizeof(Int) + nrows*(sizeof(Idx)+sizeof(T)));
  //
  //              //Here, do a linear search instead for the blkidx
  //
  //              Ptr colbeg = 1;
  //              Ptr colend = nrows;
  //              //            for(Ptr rowidx = colbeg; rowidx<=colend; ++rowidx){
  //              //              Idx row = rowind[rowidx-1];
  //              //              logfileptr->OFS()<<row<<" ";
  //              //            }
  //              //              logfileptr->OFS()<<std::endl;
  //
  //              if(colbeg<=colend){
  //                //sort rowind and nzvals
  //
  //#if 1 
  //                std::vector<size_t> lperm = sort_permutation(&rowind[colbeg-1],&rowind[colend-1]+1,std::less<Int>());
  //                apply_permutation(&rowind[colbeg-1],&rowind[colend-1]+1,lperm);
  //                apply_permutation(&nzvalA[colbeg-1],&nzvalA[colend-1]+1,lperm);
  //                Idx firstRow = rowind[colbeg-1];
  //                Int blkidx = snode->FindBlockIdx(firstRow);
  //                assert(blkidx!=-1);
  //                Int blk_nrows = snode->NRows(blkidx);    
  //                NZBlockDesc * blk_desc = &snode->GetNZBlockDesc(blkidx);
  //#else
  //                Idx prevRow = 0;
  //                Int blkidx = 0;
  //                NZBlockDesc * blk_desc = &snode->GetNZBlockDesc(blkidx);
  //                Int blk_nrows = snode->NRows(blkidx);
  //#endif
  //
  //                for(Ptr rowidx = colbeg; rowidx<=colend; ++rowidx){
  //                  Idx row = rowind[rowidx-1];
  //
  //#if 1
  //                  while(row>blk_desc->GIndex+blk_nrows-1){
  //                    blkidx++;
  //                    blk_nrows = snode->NRows(blkidx);    
  //                    blk_desc = &snode->GetNZBlockDesc(blkidx);
  //                  }
  //#else
  //                  if(row!=prevRow){
  //                    if(row>blk_desc->GIndex+blk_nrows || row<blk_desc->GIndex){
  //                      blkidx = snode->FindBlockIdx(row);
  //                      blk_desc = &snode->GetNZBlockDesc(blkidx);
  //                      blk_nrows = snode->NRows(blkidx);
  //                    }
  //                    prevRow=row;
  //                  }
  //#endif
  //
  //                  Int local_row = row - blk_desc->GIndex + 1;
  //                  Int local_col = col - fc + 1;
  //                  T * nzval = snode->GetNZval(blk_desc->Offset);
  //                  nzval[(local_row-1)*iWidth+local_col-1] = nzvalA[rowidx-1];
  //                }
  //                //for(Int rowidx = colbeg; rowidx<=colend; ++rowidx){
  //                //  Int row = rowind[rowidx-1];
  //                //  Int blkidx = snode->FindBlockIdx(row);
  //                //  assert(blkidx!=-1);
  //                //  NZBlockDesc & blk_desc = snode->GetNZBlockDesc(blkidx);
  //                //  Int local_row = row - blk_desc.GIndex + 1;
  //                //  Int local_col = col - fc + 1;
  //                //  T * nzval = snode->GetNZval(blk_desc.Offset);
  //                //  nzval[(local_row-1)*iWidth+local_col-1] = nzvalA[rowidx-1];
  //                //}
  //              }
  //            }
  //          }
  //#else
  //#if 0
  //          for(Int p=0;p<np;++p){
  //            IrecvPtr->setHead(rdispls[p]);
  //            while(IrecvPtr->head < rdispls[p]+rcounts[p]){ 
  //
  //              //Deserialize
  //              Int I,nnz;
  //              (*IrecvPtr)>>I>>nnz;
  //
  //              //          Int I = ((Int*)buffer)[0];
  //              //          Int nnz = ((Int*)buffer)[1];
  //              //          char * readBuf = buffer + 2*sizeof(Int);
  //
  //              //prepare for sorting
  //              if(1){
  //                scope_timer(a,SORTING_TRIPLETS);
  //                size_t begin = 0;
  //                size_t end   = nnz;
  //                triplet<T> * buffer = (triplet<T>*)IrecvPtr->back();
  //                //sort 
  //                std::sort(&buffer[begin],&buffer[begin]+end,[](const triplet<T> & a, const triplet<T> & b)->bool{
  //                    bool retval = a.row<b.row;
  //                    if(a.row==b.row){
  //                    retval = a.col<b.col;
  //                    }
  //                    return retval;
  //                    });
  //              }
  //
  //              Int fc = Xsuper_[I-1];
  //              Int lc = Xsuper_[I]-1;
  //              Int iWidth = lc-fc+1;
  //              SuperNode<T> * snode = snodeLocal(I);
  //              Int nnzRead = 0;
  //
  //              Idx prevRow = 0;
  //              Int blkidx = 0;
  //              NZBlockDesc * blk_desc = &snode->GetNZBlockDesc(blkidx);
  //              Int blknrows = snode->NRows(blkidx);
  //              T * nzval = snode->GetNZval(blk_desc->Offset);
  //
  //              while(nnzRead<nnz){
  //                Idx row,col;
  //                T val;
  //                (*IrecvPtr)>>row>>col>>val;
  //
  //                if(row!=prevRow){
  //                  if(row>blk_desc->GIndex+blknrows){
  //                    blkidx = snode->FindBlockIdx(row);
  //                    blk_desc = &snode->GetNZBlockDesc(blkidx);
  //                    blknrows = snode->NRows(blkidx);
  //                    nzval = snode->GetNZval(blk_desc->Offset);
  //                  }
  //                  prevRow=row;
  //                }
  //
  //                Int local_row = row - blk_desc->GIndex + 1;
  //                Int local_col = col - fc + 1;
  //                nzval[(local_row-1)*iWidth+local_col-1] = val;
  //
  //                nnzRead++;
  //              }
  //
  //            }
  //          }
  //#else
  //          IrecvPtr->setHead(0);
  //          //do a big sort
  //
  //          Int prevI = 0;
  //
  //          Idx row,col;
  //          T val;
  //
  //          Int fc = 0;
  //          Int lc = 0;
  //          Int iWidth = 0;
  //          SuperNode<T> * snode = NULL;
  //
  //          Idx prevRow = 0;
  //          Int blkidx = 0;
  //          NZBlockDesc * blk_desc = NULL;
  //          Int blknrows = 0;
  //          T * nzval = NULL;
  //
  //
  //          if(1){
  //            scope_timer(a,SORTING_TRIPLETS);
  //            size_t begin = 0;
  //            size_t end   = total_recv_size/sizeof(triplet<T>);
  //            triplet<T> * buffer = (triplet<T>*)IrecvPtr->back();
  //            //sort 
  //            std::sort(&buffer[begin],&buffer[begin]+end,[&supM = SupMembership_](const triplet<T> & a, const triplet<T> & b)->bool{
  //                Int I = supM[a.col-1];
  //                Int J = supM[b.col-1];
  //
  //                bool retval = I<J;
  //                if(I==J){
  //                retval = a.row<b.row;
  //                if(a.row==b.row){
  //                retval = a.col<b.col;
  //                }
  //                }
  //                return retval;
  //                });
  //          }
  //
  //
  //
  //
  //
  //
  //
  //          triplet<T> * elem = (triplet<T>*)IrecvPtr->back();
  //          Int i = 0;
  //          Int ntrip = total_recv_size/sizeof(triplet<T>);
  //          while(i < ntrip){
  //            triplet<T> * cur_elem = &elem[i];//(triplet<T>*)IrecvPtr->back();
  //            i++;
  //            //IrecvPtr->setHead(IrecvPtr->head+sizeof(triplet<T>));
  //            //(*IrecvPtr)>>row>>col>>val;
  //            Int I = SupMembership_[cur_elem->col-1];
  //
  //            Int iDest = this->Mapping_->Map(I-1,I-1);
  //            bassert(iam==iDest);
  //
  //            if(I!=prevI){
  //              fc = Xsuper_[I-1];
  //              lc = Xsuper_[I]-1;
  //              iWidth = lc-fc+1;
  //              snode = snodeLocal(I);
  //              prevRow = 0;
  //              nzval=NULL;
  //              bassert(cur_elem->row!=prevRow); 
  //              prevI = I;
  //            }
  //
  //            if(cur_elem->row!=prevRow){
  //              if(prevRow==0){
  //                blkidx = snode->FindBlockIdx(cur_elem->row);
  //                blk_desc = &snode->GetNZBlockDesc(blkidx);
  //                blknrows = snode->NRows(blkidx);
  //                nzval = snode->GetNZval(blk_desc->Offset);
  //              }
  //
  //              if(cur_elem->row>blk_desc->GIndex+blknrows){
  //                blkidx = snode->FindBlockIdx(cur_elem->row);
  //                blk_desc = &snode->GetNZBlockDesc(blkidx);
  //                blknrows = snode->NRows(blkidx);
  //                nzval = snode->GetNZval(blk_desc->Offset);
  //              }
  //              prevRow=cur_elem->row;
  //            }
  //            bassert(nzval!=NULL);
  //
  //            Int local_row = cur_elem->row - blk_desc->GIndex + 1;
  //            Int local_col = cur_elem->col - fc + 1;
  //
  //            //logfileptr->OFS()<<I<<" "<<blkidx<<" "<<*blk_desc<<" ("<<row<<" "<<col<<" "<<val<<") "<<local_row<<" "<<local_col<<std::endl;
  //            nzval[(local_row-1)*iWidth+local_col-1] = cur_elem->val;
  //          } 
  //#endif
  //#endif
  //
  //          SYMPACK_TIMER_SPECIAL_STOP(deserializing);      
  //
  //          //after the alltoallv, cleanup
  //          delete IrecvPtr;
  //          delete IsendPtr;
  //
  //
  //
  //
  //
  //
  //#ifdef BETTER_DISTRIB3
  //          if (0){
  //            {
  //              scope_timer(a,NEW_PACKING);
  //              for(Int I=1;I<Xsuper_.size();I++){
  //                Idx fc = Xsuper_[I-1];
  //                Idx lc = Xsuper_[I]-1;
  //                Int iWidth = lc-fc+1;
  //                Int iHeight = UpdateHeight_[I-1];
  //
  //                Int iDest = this->Mapping_->Map(I-1,I-1);
  //
  //                //look at the owner of the first column of the supernode
  //                Idx numColFirst = std::max(1,iSize_ / np);
  //
  //                auto & setI = bufSup[iDest][I];
  //
  //                //post all the recv and sends
  //                for(Idx col = fc;col<=lc;col++){
  //                  //corresponding column in the unsorted matrix A
  //                  Idx orig_col = Order_.perm[col-1];
  //                  Idx iOwnerCol = std::min((orig_col-1)/numColFirst,(Idx)np-1);
  //                  if(iam == iOwnerCol){
  //                    Int nrows = 0;
  //                    Idx local_col = (orig_col-(numColFirst)*iOwnerCol);
  //                    for(Ptr rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
  //                      Idx orig_row = pMat.Local_.rowind[rowidx-1];
  //                      Idx row = Order_.invp[orig_row-1];
  //
  //                      if(row<col){
  //                        //add the pair (col,row) to processor owning column row
  //                        Int J = SupMembership_[row-1];
  //                        Int iDestJ = this->Mapping_->Map(J-1,J-1);
  //                        auto & set = bufSup[iDestJ][J];
  //                        T & val = pMat.nzvalLocal[rowidx-1];
  //                        triplet<T> trip;
  //                        trip.row = col;
  //                        trip.col = row;
  //                        trip.val = val;
  //                        //set.insert(trip);
  //                        set.push(trip);
  //                      }
  //                      else{
  //                        T & val = pMat.nzvalLocal[rowidx-1];
  //                        triplet<T> trip;
  //                        trip.row = row;
  //                        trip.col = col;
  //                        trip.val = val;
  //                        //setI.insert(trip);
  //                        setI.push(trip);
  //                      }
  //                    }
  //                  }
  //                }
  //              }
  //            }
  //
  //            Icomm * IsendPtr;
  //            Icomm * IrecvPtr;
  //            std::vector<int> scounts(np,0);
  //            std::vector<int> sdispls(np+1,0);
  //            size_t total_send_size = 0;
  //
  //
  //            {
  //              scope_timer_special(a,NEW_SERIALIZING);
  //
  //              //map<Int, unordered_map<Int, std::set<triplet<T>, sortTriplet<T> > > > bufSup;
  //              if(bufSup.size()>0){
  //                auto prevProcIt = bufSup.begin();
  //                for(Int p = 0 ; p< prevProcIt->first; ++p){
  //                  scounts[p] = 0;
  //                }
  //
  //                for(auto procIt = bufSup.begin(); procIt != bufSup.end(); procIt++){
  //                  //fill the holes
  //                  for(Int p =prevProcIt->first+1; p<procIt->first; ++p){
  //                    scounts[p] = 0;
  //                  }
  //
  //                  Int iCurDest = procIt->first;
  //                  int send_bytes = 0;
  //                  auto & supMap = procIt->second;
  //                  for(auto supIt = supMap.begin() ; supIt != supMap.end(); supIt++){
  //                    Int I = supIt->first;
  //                    auto & supSet = supIt->second;
  //                    send_bytes += supSet.size()*sizeof(triplet<T>);
  //                    total_send_size+= supSet.size()*sizeof(triplet<T>);
  //                  }
  //                  scounts[iCurDest] = send_bytes;
  //                  prevProcIt = procIt;
  //                }
  //
  //                for(Int p = bufSup.rbegin()->first+1; p<np; ++p){
  //                  scounts[p] = 0;
  //                }
  //              }
  //              else{
  //                for(Int p = 0; p<np; ++p){
  //                  scounts[p] = 0;
  //                }
  //              }
  //              total_send_size = 0;
  //
  //              //compute send displacements
  //              sdispls[0] = 0;
  //              std::partial_sum(scounts.begin(),scounts.end(),&sdispls[1]);
  //              total_send_size = sdispls.back();
  //
  //              IsendPtr = new Icomm(total_send_size,MPI_REQUEST_NULL);
  //              for(auto procIt = bufSup.begin(); procIt != bufSup.end(); procIt++){
  //                Int iCurDest = procIt->first;
  //                auto & supMap = procIt->second;
  //                IsendPtr->setHead(sdispls[iCurDest]);
  //
  //                for(auto supIt = supMap.begin() ; supIt != supMap.end(); supIt++){
  //                  Int I = supIt->first;
  //                  auto & supSet = supIt->second;
  //                  //std::copy(supSet.begin(),supSet.end(),(triplet<T>*)IsendPtr->back());
  //                  //size_t sz = supSet.size()*sizeof(triplet<T>);
  //                  //IsendPtr->setHead(IsendPtr->head + sz );
  //                  scope_timer(b,COPY_PQUEUE);
  //                  while(!supSet.empty()){
  //                    const triplet<T> & trip = supSet.top();
  //                    //            logfileptr->OFS()<<"PACK ("<<trip.row<<" "<<trip.col<<" "<<trip.val<<")"<<std::endl;
  //
  //                    *IsendPtr<<supSet.top();
  //                    supSet.pop();
  //                  }
  //                }
  //              }
  //            }
  //            //compute the receive structures
  //            size_t total_recv_size = 0;
  //            std::vector<int> rdispls(np+1,0);
  //            std::vector<int> rcounts(np,0);
  //
  //            MPI_Alltoall(&scounts[0],sizeof(int),MPI_BYTE,&rcounts[0],sizeof(int),MPI_BYTE,pMat.comm);
  //
  //            //compute receive displacements
  //            rdispls[0] = 0;
  //            std::partial_sum(rcounts.begin(),rcounts.end(),&rdispls[1]);
  //            total_recv_size = rdispls.back();
  //
  //            IrecvPtr = new Icomm(total_recv_size,MPI_REQUEST_NULL);
  //
  //            MPI_Alltoallv(IsendPtr->front(), &scounts[0], &sdispls[0], MPI_BYTE,
  //                IrecvPtr->front(), &rcounts[0], &rdispls[0], MPI_BYTE,
  //                pMat.comm);
  //
  //            //Need to parse the structure sent from the processor owning the first column of the supernode
  //            {
  //              scope_timer_special(a,NEW_DESERIALIZING);
  //              IrecvPtr->setHead(0);
  //              //do a big sort
  //
  //              Int prevI = 0;
  //
  //              Idx row,col;
  //              T val;
  //
  //              Int fc = 0;
  //              Int lc = 0;
  //              Int iWidth = 0;
  //              SuperNode<T> * snode = NULL;
  //
  //              Idx prevRow = 0;
  //              Int blkidx = 0;
  //              NZBlockDesc * blk_desc = NULL;
  //              Int blknrows = 0;
  //              T * nzval = NULL;
  //
  //
  //              if(0){
  //                scope_timer(a,SORTING_TRIPLETS);
  //                size_t begin = 0;
  //                size_t end   = total_recv_size/sizeof(triplet<T>);
  //                triplet<T> * buffer = (triplet<T>*)IrecvPtr->back();
  //                //sort 
  //                std::sort(&buffer[begin],&buffer[begin]+end,[&supM = SupMembership_](const triplet<T> & a, const triplet<T> & b)->bool{
  //                    Int I = supM[a.col-1];
  //                    Int J = supM[b.col-1];
  //
  //                    bool retval = I<J;
  //                    if(I==J){
  //                    retval = a.row<b.row;
  //                    if(a.row==b.row){
  //                    retval = a.col<b.col;
  //                    }
  //                    }
  //                    return retval;
  //                    });
  //              }
  //
  //
  //
  //
  //
  //
  //
  //              triplet<T> * elem = (triplet<T>*)IrecvPtr->back();
  //              Int i = 0;
  //              Int ntrip = total_recv_size/sizeof(triplet<T>);
  //              while(i < ntrip){
  //                triplet<T> * cur_elem = &elem[i];//(triplet<T>*)IrecvPtr->back();
  //                i++;
  //                //IrecvPtr->setHead(IrecvPtr->head+sizeof(triplet<T>));
  //                //(*IrecvPtr)>>row>>col>>val;
  //                Int I = SupMembership_[cur_elem->col-1];
  //
  //                Int iDest = this->Mapping_->Map(I-1,I-1);
  //                bassert(iam==iDest);
  //
  //                if(I!=prevI){
  //                  fc = Xsuper_[I-1];
  //                  lc = Xsuper_[I]-1;
  //                  iWidth = lc-fc+1;
  //                  snode = snodeLocal(I);
  //                  prevRow = 0;
  //                  nzval=NULL;
  //                  bassert(cur_elem->row!=prevRow); 
  //                  prevI = I;
  //                }
  //
  //                if(cur_elem->row!=prevRow){
  //                  if(prevRow==0){
  //                    blkidx = snode->FindBlockIdx(cur_elem->row);
  //                    blk_desc = &snode->GetNZBlockDesc(blkidx);
  //                    blknrows = snode->NRows(blkidx);
  //                    nzval = snode->GetNZval(blk_desc->Offset);
  //                  }
  //
  //                  if(cur_elem->row>blk_desc->GIndex+blknrows){
  //                    blkidx = snode->FindBlockIdx(cur_elem->row);
  //                    blk_desc = &snode->GetNZBlockDesc(blkidx);
  //                    blknrows = snode->NRows(blkidx);
  //                    nzval = snode->GetNZval(blk_desc->Offset);
  //                  }
  //                  prevRow=cur_elem->row;
  //                }
  //                bassert(nzval!=NULL);
  //
  //                Int local_row = cur_elem->row - blk_desc->GIndex + 1;
  //                Int local_col = cur_elem->col - fc + 1;
  //
  //                //            logfileptr->OFS()<<I<<" "<<blkidx<<" "<<*blk_desc<<" ("<<row<<" "<<col<<" "<<val<<") "<<local_row<<" "<<local_col<<std::endl;
  //                //            logfileptr->OFS()<<"("<<cur_elem->row<<" "<<cur_elem->col<<" "<<cur_elem->val<<")"<<std::endl;
  //                nzval[(local_row-1)*iWidth+local_col-1] = cur_elem->val;
  //              }
  //            }
  //
  //
  //
  //            delete IrecvPtr;
  //            delete IsendPtr;
  //
  //          }
  //#endif
  //        }













  template <typename T> void symPACKMatrix<T>::Init(symPACKOptions & options ){
    SYMPACK_TIMER_START(A);
    //scope_timer_special(a,Initialization);
#ifdef _STAT_COMM_
    maxAggreg_ =0;
    sizesAggreg_ =0;
    countAggreg_=0;
    maxFactors_ =0;
    sizesFactors_ =0;
    countFactors_=0;

    maxAggregRecv_ =0;
    sizesAggregRecv_ =0;
    countAggregRecv_=0;
    maxFactorsRecv_ =0;
    sizesFactorsRecv_ =0;
    countFactorsRecv_=0;
#endif

    scope_timer(b,ADVANCE_OUTGOING_COMM);
    //Create the CommEnvironment object if necessary
    if(CommEnv_!=NULL){
      delete CommEnv_;
    }

    //    if(options.commEnv == NULL){
    //      throw std::runtime_error("The communication environment must be initialized in the options");
    //    }


    options_ = options;

    logfileptr->verbose = options_.verbose>0;

    all_np = 0;
    MPI_Comm_size(options_.MPIcomm,&all_np);
    MPI_Comm_rank(options_.MPIcomm,&iam);
    np = options_.used_procs(all_np);


    MPI_Comm workcomm = MPI_COMM_NULL;
    MPI_Comm_split(options_.MPIcomm,iam<np,iam,&workcomm);
    CommEnv_ = new CommEnvironment(workcomm);
    MPI_Comm_free(&workcomm);

    //do another split to contain P0 and all the non working processors
    if(all_np!=np){
      non_workcomm_ = MPI_COMM_NULL;
      MPI_Comm_split(options_.MPIcomm,iam==0||iam>=np,iam==0?0:iam-np+1,&non_workcomm_);
    }

    //SYMPACK_TIMER_START(MAPPING);

    //SYMPACK_TIMER_STOP(MAPPING);

    if(iam<np){
      //Options
      maxIsend_ = options_.maxIsend;
      maxIrecv_ = options_.maxIrecv;
      gMaxIrecv = maxIrecv_+maxIsend_;

      //SYMPACK_TIMER_START(SCHEDULER);
      switch(options_.scheduler){
        case DL:
          scheduler_ = new DLScheduler<std::list<FBTask>::iterator>();
          scheduler2_ = new DLScheduler<FBTask>();
          break;
        case MCT:
          scheduler_ = new MCTScheduler<std::list<FBTask>::iterator>();
          scheduler2_ = new MCTScheduler<FBTask>();
          break;
        case PR:
          scheduler_ = new PRScheduler<std::list<FBTask>::iterator>();
          scheduler2_ = new PRScheduler<FBTask>();
          break;
        case FIFO:
          scheduler_ = new FIFOScheduler<std::list<FBTask>::iterator>();
          scheduler2_ = new FIFOScheduler<FBTask>();
          break;
        default:
          scheduler_ = new DLScheduler<std::list<FBTask>::iterator>();
          scheduler2_ = new DLScheduler<FBTask>();
          break;
      }
      //SYMPACK_TIMER_STOP(SCHEDULER);

      //Local_ = new SparseMatrixStructure();

      //Create an Ordering object to hold the permutation
      //Order_.SetCommEnvironment(CommEnv_);
      //    MPI_Barrier(CommEnv_->MPI_GetComm());
    }

    SYMPACK_TIMER_STOP(A);
    //exit(-1);
  }

  template <typename T> void symPACKMatrix<T>::SymbolicFactorization(DistSparseMatrix<T> & pMat){
    scope_timer(a,SymbolicFactorization);

    if(fullcomm_!=MPI_COMM_NULL){
      MPI_Comm_free(&fullcomm_);
    }

    MPI_Comm_dup(pMat.comm,&fullcomm_);

    iSize_ = pMat.size;
#ifdef EXPLICIT_PERMUTE
    pMat.Localg_.SetSorted(0);
    //this is done on the full communicator
    pMat.ExpandSymmetric();

    graph_ = pMat.GetLocalGraph();
    graph_.SetBaseval(1);
    graph_.SetKeepDiag(0);
    graph_.SetSorted(1);
#else
    //*Local_ = pMat.GetLocalStructure();
    //graph_.SetComm(CommEnv_->MPI_GetComm());
    graph_ = pMat.GetLocalGraph();
    graph_.SetBaseval(1);
    graph_.SetSorted(1);
    //graph_.FromStructure(*Local_);
    graph_.ExpandSymmetric();
#endif    


    logfileptr->OFS()<<"Matrix structure expanded"<<std::endl;

    SparseMatrixGraph * sgraph = NULL;
    {
      {

        double timeSta = get_time();
        SYMPACK_TIMER_START(ORDERING);


        if(options_.orderingStr=="MMD"){
          options_.ordering = symPACK::MMD;
        }
        else if(options_.orderingStr=="AMD"){
          options_.ordering = symPACK::AMD;
        }
        else if(options_.orderingStr=="METIS"){
          options_.ordering = symPACK::METIS;
        }
        else if(options_.orderingStr=="SCOTCH"){
          options_.ordering = symPACK::SCOTCH;
        }
        else if(options_.orderingStr=="PARMETIS"){
          options_.ordering = symPACK::PARMETIS;
        }
        else if(options_.orderingStr=="PTSCOTCH"){
          options_.ordering = symPACK::PTSCOTCH;
        }
        else if(options_.orderingStr=="NATURAL"){
          options_.ordering = symPACK::NATURAL;
        }
        else if(options_.orderingStr=="USER"){
          if(options_.perm ==NULL){
            throw std::logic_error( "When using USER, symPACKOptions.perm must be provided.\n" );
          }
          else{
            options_.ordering = symPACK::USER;
          }
        }
        else{
          throw std::logic_error( "Unknown ordering technique.\n" );
        }



        Order_.NpOrdering = options_.NpOrdering;
        switch(options_.ordering){
          case MMD:
            {
              sgraph = new SparseMatrixGraph();
              graph_.GatherStructure(*sgraph,0);
              sgraph->SetBaseval(1);
              sgraph->SetKeepDiag(0);
              Order_.MMD(*sgraph, graph_.comm);
            }
            break;

          case AMD:
            {
              sgraph = new SparseMatrixGraph();
              graph_.GatherStructure(*sgraph,0);
              sgraph->SetBaseval(1);
              sgraph->SetKeepDiag(0);
              Order_.AMD(*sgraph, graph_.comm);
            }
            break;

          case NDBOX:
            Order_.NDBOX(Size(), graph_.comm);
            break;

          case NDGRID:
            Order_.NDGRID(Size(), graph_.comm);
            break;

#ifdef USE_SCOTCH
          case SCOTCH:
            {

              sgraph = new SparseMatrixGraph();
              graph_.GatherStructure(*sgraph,0);
              sgraph->SetKeepDiag(0);
              Order_.SCOTCH(*sgraph, graph_.comm);
            }
            break;
#endif
#ifdef USE_METIS
          case METIS:
            {
              sgraph = new SparseMatrixGraph();
              graph_.GatherStructure(*sgraph,0);
              sgraph->SetKeepDiag(0);
              Order_.METIS(*sgraph, graph_.comm);
            }
            break;
#endif
#ifdef USE_PARMETIS
          case PARMETIS:
            {
              graph_.SetKeepDiag(0);
              Order_.PARMETIS(graph_);
            }
            break;
#endif
#ifdef USE_PTSCOTCH
          case PTSCOTCH:
            {
              graph_.SetKeepDiag(0);
              Order_.PTSCOTCH(graph_);
            }
            break;
#endif
          case NATURAL:
            Order_.perm.resize(iSize_);
            for(int i=0;i<Order_.perm.size();i++){Order_.perm[i]=i+1;} 
            Order_.invp = Order_.perm;
          default:

            //do nothing: either natural or user provided ordering
            break;
        }
        SYMPACK_TIMER_STOP(ORDERING);

        //The ordering is available on every processor of the full communicator
        double timeStop = get_time();
        if(iam==0){
          std::cout<<"Ordering time: "<<timeStop - timeSta<<std::endl;
        }
        logfileptr->OFS()<<"Ordering done"<<std::endl;
      }
    }

    if(options_.dumpPerm>0){
      logfileptr->OFS()<<"perm = [";
      for (auto i : Order_.perm){ 
        logfileptr->OFS()<<i<<" ";
      }
      logfileptr->OFS()<<"]"<<std::endl;
    }

    std::vector<Int> & cc = cc_;
    std::vector<Int> &rc = rc_;
    //gather the graph if necessary to build the elimination tree
    if(sgraph==NULL){ 
      sgraph = new SparseMatrixGraph();
      graph_.GatherStructure(*sgraph,0);
    }
    //      if(iam<np){
    graph_.SetKeepDiag(1);
    sgraph->SetBaseval(1);
    sgraph->SetKeepDiag(1);

    {

      double timeSta = get_time();
      ETree_.ConstructETree(*sgraph,Order_,fullcomm_);
      ETree_.PostOrderTree(Order_);
      logfileptr->OFS()<<"ETREE done"<<std::endl;
      {
        double timeStart = get_time();
        this->getLColRowCount(*sgraph,cc,rc);
        double timeStop = get_time();
        if(iam==0){
          std::cout<<"Column count (gather + serial + bcast) construction time: "<<timeStop - timeStart<<std::endl;
        }
      }

      if(options_.ordering != NATURAL){
        ETree_.SortChildren(cc,Order_);

        if(options_.dumpPerm>0){
          logfileptr->OFS()<<"perm = [";
          for (auto i : Order_.perm){ 
            logfileptr->OFS()<<i<<" ";
          }
          logfileptr->OFS()<<"]"<<std::endl;
        }


      }

      double timeStop = get_time();
      if(iam==0){
        std::cout<<"Elimination tree construction time: "<<timeStop - timeSta<<std::endl;
      }
    }

    //compute some statistics
    if(iam==0){
      double flops = 0.0;
      int64_t NNZ = 0;
      for(Int i = 0; i<cc.size();++i){
        flops+= (double)pow((double)cc[i],2.0);
        NNZ+=cc[i];
      }
      std::cout<<"Flops: "<<flops<<std::endl;
      std::cout<<"NNZ in L factor: "<<NNZ<<std::endl;
    }
    //      }


    //get rid of the sequential graph
    if(sgraph!=NULL){delete sgraph;}

    { 
      double timeSta = get_time();
      this->findSupernodes(ETree_,Order_,cc,SupMembership_,Xsuper_,options_.relax.maxSize);
      logfileptr->OFS()<<"Supernodes found"<<std::endl;

      if(options_.relax.nrelax0>0){
        this->relaxSupernodes(ETree_, cc,SupMembership_, Xsuper_, options_.relax );
        logfileptr->OFS()<<"Relaxation done"<<std::endl;
      }

      //modify np since it cannot be greater than the number of supernodes
      this->np = std::min(this->np,this->options_.used_procs(Xsuper_.size()-1)); 

      MPI_Comm workcomm = MPI_COMM_NULL;
      MPI_Comm_split(options_.MPIcomm,iam<np,iam,&workcomm);
      bassert(CommEnv_!=NULL);
      delete CommEnv_;
      CommEnv_ = new CommEnvironment(workcomm);
      MPI_Comm_free(&workcomm);

      //do another split to contain P0 and all the non working processors
      if(all_np!=np){
        if(non_workcomm_!=MPI_COMM_NULL){
          MPI_Comm_free(&non_workcomm_);
        }
        non_workcomm_ = MPI_COMM_NULL;
        MPI_Comm_split(options_.MPIcomm,iam==0||iam>=np,iam==0?0:iam-np+1,&non_workcomm_);
      }


      //now create the mapping
      if(this->Mapping_ != NULL){
        delete this->Mapping_;
      }

      Int pmapping = np;
      Int pr = (Int)sqrt((double)pmapping);
      if(options_.mappingTypeStr ==  "MODWRAP2DTREE"){
        this->Mapping_ = new ModWrap2DTreeMapping(pmapping, pr, pr, 1);
      }
      else if(options_.mappingTypeStr ==  "WRAP2D"){
        this->Mapping_ = new Wrap2D(pmapping, pr, pr, 1);
      }
      else if(options_.mappingTypeStr ==  "WRAP2DFORCED"){
        this->Mapping_ = new Wrap2DForced(pmapping, pr, pr, 1);
      }
      else if(options_.mappingTypeStr ==  "MODWRAP2DNS"){
        this->Mapping_ = new Modwrap2DNS(pmapping, pr, pr, 1);
      }
      else if(options_.mappingTypeStr ==  "ROW2D"){
        this->Mapping_ = new Row2D(pmapping, pmapping, pmapping, 1);
      }
      else if(options_.mappingTypeStr ==  "COL2D"){
        this->Mapping_ = new Col2D(pmapping, pmapping, pmapping, 1);
      }
      else{
        this->Mapping_ = new Modwrap2D(pmapping, pr, pr, 1);
      }



      //Compute XsuperDist_
      std::vector<Idx> newVertexDist;
      {
        Idx supPerProc = std::max((size_t)1,(Xsuper_.size()-1) / np);
        XsuperDist_.resize(all_np+1,0);
        XsuperDist_[np] = Xsuper_.size();
        for(int p = 0; p<np; p++){
          XsuperDist_[p]= std::min(XsuperDist_[np],(Int)(p*supPerProc+1));
        }
        for(int p = np+1; p<all_np; p++){
          XsuperDist_[p]= XsuperDist_[p-1];
        }
        XsuperDist_[all_np] = Xsuper_.size();


        newVertexDist.resize(all_np+1,0);
        newVertexDist[all_np] = this->iSize_+1;
        for(int p = 0; p < all_np; p++){
          Int S = XsuperDist_[p];
          newVertexDist[p] = Xsuper_[S-1];
        }
      }


#ifdef EXPLICIT_PERMUTE
      {
        //This should permute the matrix A and move it
        //onto np processor instead of all_np processors
        pMat.Permute(&Order_.invp[0],&newVertexDist[0]);
        graph_ = pMat.GetLocalGraph();
        graph_.SetBaseval(1);
        graph_.SetKeepDiag(1);
        graph_.SetSorted(1);
      }
#endif

      //if(iam<np){
      double timeStaSymb = get_time();
      this->symbolicFactorizationRelaxedDist(cc);

      double timeStopSymb = get_time();
      if(iam==0){
        std::cout<<"Symbolic factorization time: "<<timeStopSymb - timeStaSymb<<std::endl;
      }
      logfileptr->OFS()<<"Symbfact done"<<std::endl;
      //}

      if(options_.order_refinement_str == "SET") {
        //abort();
        //          if(iam<np){
        this->refineSupernodes(2,1,&pMat);

        if(iam==0){
          std::cout<<"Supernode reodering done"<<std::endl;
        }
        //          }

        if(options_.relax.nrelax0>0){
          //            if(iam<np){
          this->relaxSupernodes(ETree_, cc,SupMembership_, Xsuper_, options_.relax );
          logfileptr->OFS()<<"Relaxation done"<<std::endl;
          //            }

          //            //Bcast Xsuper_ again
          //            if( (iam==0 || iam>=np) && all_np!=np){
          //              //MPI_Bcast Xsuper_ to everyone else
          //              size_t nSuper = Xsuper_.size()-1;
          //              MPI_Bcast(&nSuper,sizeof(size_t),MPI_BYTE,0,non_workcomm_);
          //    
          //              Xsuper_.resize(nSuper+1);      
          //              MPI_Bcast(&Xsuper_[0],(nSuper+1)*sizeof(Int),MPI_BYTE,0,non_workcomm_);
          //            }

          //Refresh XsuperDist_
          {
            Idx supPerProc = std::max((size_t)1,(Xsuper_.size()-1) / np);
            XsuperDist_.resize(np+1,0);
            XsuperDist_[np] = Xsuper_.size();
            for(int p = 0; p<np; p++){
              XsuperDist_[p]= std::min(XsuperDist_[np],(Int)(p*supPerProc+1));
            }
            for(int p = np+1; p<all_np; p++){
              XsuperDist_[p]= XsuperDist_[p-1];
            }
            XsuperDist_[all_np] = Xsuper_.size();
          }
        }


        //          if(iam<np){
        {
          double timeSta = get_time();
          this->symbolicFactorizationRelaxedDist(cc);
          double timeStop = get_time();
          if(iam==0){
            std::cout<<"Symbolic factorization time: "<<timeStop - timeSta<<std::endl;
          }
          logfileptr->OFS()<<"Symbfact done"<<std::endl;
        }
        //          }
      }

      double timeStop = get_time();
      if(iam==0){
        std::cout<<"Total symbolic factorization time: "<<timeStop - timeSta<<std::endl;
      }
    }

#ifdef _DEBUG_
    logfileptr->OFS()<<"Membership list is "<<SupMembership_<<std::endl;
    logfileptr->OFS()<<"xsuper "<<Xsuper_<<std::endl;
#endif




    {
      scope_timer(a,LOAD_BALANCE);
      double timeSta = get_time();

      if (this->Balancer_!=NULL){
        delete this->Balancer_;
      }

      std::vector<Int> map;
      if(options_.load_balance_str=="SUBCUBE-FI"){
        if(iam==0){ std::cout<<"Subtree to subcube FI mapping used"<<std::endl;}
        ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
        this->Balancer_ = new SubtreeToSubcube(np,SupETree,Xsuper_,XsuperDist_,SupMembership_,locXlindx_,locLindx_,cc,fullcomm_,true);
      }
      else if(options_.load_balance_str=="SUBCUBE-FO"){
        if(iam==0){ std::cout<<"Subtree to subcube FO mapping used"<<std::endl;}
        ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
        this->Balancer_ = new SubtreeToSubcube(np,SupETree,Xsuper_,XsuperDist_,SupMembership_,locXlindx_,locLindx_,cc,fullcomm_,false);
      }
      else if(options_.load_balance_str=="SUBCUBE-VOLUME-FI"){
        if(iam==0){ std::cout<<"Subtree to subcube volume FI mapping used"<<std::endl;}
        ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
        this->Balancer_ = new SubtreeToSubcubeVolume(np,SupETree,Xsuper_,XsuperDist_,SupMembership_,locXlindx_,locLindx_,cc,fullcomm_,true);
      }
      else if(options_.load_balance_str=="SUBCUBE-VOLUME-FO"){
        if(iam==0){ std::cout<<"Subtree to subcube volume FO mapping used"<<std::endl;}
        ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
        this->Balancer_ = new SubtreeToSubcubeVolume(np,SupETree,Xsuper_,XsuperDist_,SupMembership_,locXlindx_,locLindx_,cc,fullcomm_,false);
      }
      else if(options_.load_balance_str=="NNZ"){
        if(iam==0){ std::cout<<"Load Balancing on NNZ used"<<std::endl;}
        this->Balancer_ = new NNZBalancer(np,Xsuper_,cc);
      }
      else if(options_.load_balance_str=="WORK"){
        if(iam==0){ std::cout<<"Load Balancing on WORK used"<<std::endl;}
        this->Balancer_ = new WorkBalancer(np,Xsuper_,cc);
      }

      if (this->Balancer_!=NULL){
        map = this->Balancer_->GetMap();
        TreeLoadBalancer * testBalancer = dynamic_cast<TreeLoadBalancer*>(this->Balancer_);
        if(testBalancer==NULL){
          this->Mapping_->Update(map);
        }
        else{

          TreeMapping * test = dynamic_cast<TreeMapping*>(this->Mapping_);
          if(test==NULL){
            this->Mapping_->Update(map);
          }
          else{
            //              if(np!=all_np){
            //              abort();
            //              }
            test->Update((TreeLoadBalancer*)this->Balancer_);
          }
        }
      }
      double timeStop = get_time();
      if(iam==0){
        std::cout<<"Load balancing time: "<<timeStop - timeSta<<std::endl;
      }
    }


    //Now split supernodes larger than the maximum size. XlindxLocal_, LindxLocal_, Xsuper_, SupMembership_, Balancer_ and Mapping_ need to up updated




#define _OUTPUT_ETREE_
#ifdef _OUTPUT_ETREE_
    logfileptr->OFS()<<"ETree is "<<ETree_<<std::endl;
    {
      auto supETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
    logfileptr->OFS()<<"Supernodal ETree is "<<supETree<<std::endl;
    logfileptr->OFS()<<"Mapping is ";for(Int i = 0;i<supETree.Size();i++){logfileptr->OFS()<<this->Mapping_->Map(i,i)<<" ";}logfileptr->OFS()<<std::endl;
    logfileptr->OFS()<<"Xsuper is "<<Xsuper_<<std::endl;//;for(Int i = 0;i<supETree.Size();i++){logfileptr->OFS()<<this->Xsuper_->Map(i,i)<<" ";}logfileptr->OFS()<<std::endl;
    } 
#endif

    //starting from now, this only concerns the working processors
    SYMPACK_TIMER_START(Get_UpdateCount);
    GetUpdatingSupernodeCount(UpdateCount_,UpdateWidth_,UpdateHeight_,numBlk_);
    SYMPACK_TIMER_STOP(Get_UpdateCount);

    if(iam<np){

      {
        double timeSta = get_time();
        origTaskLists_.resize(Xsuper_.size(),NULL);
        std::vector<Int> AggregatesToRecv;
        std::vector<Int> LocalAggregates;
        FBGetUpdateCount(UpdatesToDo_,AggregatesToRecv,LocalAggregates);

        generateTaskGraph(taskGraph_, AggregatesToRecv, LocalAggregates);


#define _OUTPUT_TASK_GRAPH_
#ifdef _OUTPUT_TASK_GRAPH_
        logfileptr->OFS()<<"tasks: ";
        for(auto binit = taskGraph_.taskLists_.begin(); binit != taskGraph_.taskLists_.end(); binit++){
          if((*binit)!=NULL){
          for(auto taskit = (*binit)->begin(); taskit!= (*binit)->end(); taskit++){
            logfileptr->OFS()<<"t_"<<taskit->src_snode_id<<"_"<<taskit->tgt_snode_id<<" ";
          }
          }
      
        }
        logfileptr->OFS()<<std::endl; 
        logfileptr->OFS()<<"taskmap: ";
        for(auto binit = taskGraph_.taskLists_.begin(); binit != taskGraph_.taskLists_.end(); binit++){
          if((*binit)!=NULL){
          for(auto taskit = (*binit)->begin(); taskit!= (*binit)->end(); taskit++){
        //for(auto&& bin : taskGraph_.taskLists_)
        //  for(auto&& taskp : *bin)
            logfileptr->OFS()<<iam<<" ";
          }
          }
        }
        logfileptr->OFS()<<std::endl; 
#endif

        double timeStop = get_time();
        if(iam==0){
          std::cout<<"Task graph generation time: "<<timeStop - timeSta<<std::endl;
        }
      }
      MPI_Barrier(CommEnv_->MPI_GetComm());
    }













    Int snodeCount = 0;
    for(Int I=1;I<Xsuper_.size();I++){
      Int iDest = this->Mapping_->Map(I-1,I-1);

      if(iDest==iam){
        ++snodeCount;
      }
    }

    //Resize the local supernodes array
    LocalSupernodes_.reserve(snodeCount);

    //        std::vector<upcxx::global_ptr<SuperNodeDesc > > localFactors;
    //        remoteFactors_.resize(Xsuper_.size()-1);
    //        std::fill((char*)&remoteFactors_[0],(char*)&remoteFactors_[0]+remoteFactors_.size()*sizeof(std::tuple<upcxx::global_ptr<SuperNodeDesc>,Int> ),0);
    //        //logfileptr->OFS()<<"My usable global memory size is: "<<upcxx::my_usable_global_memory_size()<<std::endl;

    SYMPACK_TIMER_START(DISTRIBUTE_CREATE_SNODES);

    for(Int I=1;I<Xsuper_.size();I++){

      Int iDest = this->Mapping_->Map(I-1,I-1);

      //parse the first column to create the supernode structure
      if(iam==iDest){
        Int fc = Xsuper_[I-1];
        Int lc = Xsuper_[I]-1;
        Int iWidth = lc-fc+1;
        Int iHeight = UpdateHeight_[I-1];
        Int nzBlockCnt = numBlk_[I-1];
#ifndef ITREE2
        globToLocSnodes_.push_back(I-1);
#else 
        ITree::Interval snode_inter = { I, I, LocalSupernodes_.size() };
        globToLocSnodes_.Insert(snode_inter);
#endif
        SuperNode<T> * newSnode = CreateSuperNode(options_.decomposition,I,fc,lc,iHeight,iSize_,nzBlockCnt);
        LocalSupernodes_.push_back(newSnode);
      }
    }
    SYMPACK_TIMER_STOP(DISTRIBUTE_CREATE_SNODES);

    //first create the structure of every supernode
    {
      scope_timer(a,DISTRIBUTE_CREATE_REMOTE_NODES);
      std::vector< int > superStructure(all_np);
      std::vector< int > rdisplsStructure;
      std::vector< int > sSuperStructure(all_np);
      std::vector< int > ssizes(all_np,0);
      std::vector< int > sdispls(all_np+1,0);
      std::vector< int > rsizes(all_np,0);

      Int numLocSnode = XsuperDist_[iam+1]-XsuperDist_[iam];
      Int firstSnode = XsuperDist_[iam];

      for(Int locsupno = 1; locsupno<locXlindx_.size(); ++locsupno){
        Idx I = locsupno + firstSnode-1;
        Int iDest = this->Mapping_->Map(I-1,I-1);
        ssizes[iDest] += 1 + 2*numBlk_[I-1]; //1 for supno index + numBlk startrows + numBlk number of rows
      }

      sdispls[0] = 0;
      std::partial_sum(ssizes.begin(),ssizes.end(),&sdispls[1]);

      rdisplsStructure = sdispls;
      sSuperStructure.resize(sdispls.back());
      for(Int locsupno = 1; locsupno<locXlindx_.size(); ++locsupno){
        Idx I = locsupno + firstSnode-1;
        Int iDest = this->Mapping_->Map(I-1,I-1);

        int & tail = rdisplsStructure[iDest];


        Int fc = Xsuper_[I-1];
        Int lc = Xsuper_[I]-1;
        Int iWidth = lc - fc + 1;
        Ptr lfi = locXlindx_[locsupno-1];
        Ptr lli = locXlindx_[locsupno]-1;

        sSuperStructure[tail++] = I;
        //count number of contiguous rows
        for(Ptr sidx = lfi; sidx<=lli;sidx++){
          Idx iStartRow = locLindx_[sidx-1];
          Idx iPrevRow = iStartRow;
          Int iContiguousRows = 1;
          for(Int idx2 = sidx+1; idx2<=lli;idx2++){
            Idx iCurRow = locLindx_[idx2-1];
            if(iStartRow == locLindx_[lfi-1]){
              if(iCurRow>iStartRow+iWidth-1){
                //enforce the first block to be a square diagonal block
                break;
              }
            }

            if(iCurRow==iPrevRow+1){
              sidx++;
              ++iContiguousRows;
              iPrevRow=iCurRow;
            }
            else{
              break;
            }
          }

          sSuperStructure[tail++] = iStartRow;
          sSuperStructure[tail++] = iContiguousRows;
        }
      }

      MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,pMat.comm);
      rdisplsStructure[0] = 0;
      std::partial_sum(rsizes.begin(),rsizes.end(),&rdisplsStructure[1]);
      superStructure.resize(rdisplsStructure.back());

      //turn everything into byte sizes
      for(int p = 0; p<ssizes.size();p++){ ssizes[p]*=sizeof(int); }
      for(int p = 0; p<rsizes.size();p++){ rsizes[p]*=sizeof(int); }
      for(int p = 0; p<sdispls.size();p++){ sdispls[p]*=sizeof(int); }
      for(int p = 0; p<rdisplsStructure.size();p++){ rdisplsStructure[p]*=sizeof(int); }

      //Do the alltoallv to get the structures        
      MPI_Alltoallv(&sSuperStructure[0], &ssizes[0], &sdispls[0], MPI_BYTE,
          &superStructure[0], &rsizes[0], &rdisplsStructure[0], MPI_BYTE,
          pMat.comm);

      //loop through received structure and create supernodes
      for(Int p = 0; p<all_np; p++){
        int pos = rdisplsStructure[p]/sizeof(int);
        int end = rdisplsStructure[p+1]/sizeof(int);
        while(pos<end){
          Int I = superStructure[pos++];
          Int nzBlockCnt = numBlk_[I-1];

          Int fc = Xsuper_[I-1];
          Int lc = Xsuper_[I]-1;
          Int iWidth = lc-fc+1;
          SuperNode<T> * snode = snodeLocal(I);
          if(snode->NZBlockCnt()==0){
            for(Int i = 0; i<nzBlockCnt;i++){
              Int iStartRow = superStructure[pos++];
              Int iContiguousRows = superStructure[pos++];
              snode->AddNZBlock(iContiguousRows , iWidth,iStartRow);
            }
            snode->Shrink();
          }
        }
      } 
    }










    MPI_Barrier(fullcomm_);
  }

  template <typename T> void symPACKMatrix<T>::DistributeMatrix(DistSparseMatrix<T> & pMat){
    scope_timer(a,DistributeA);
#ifdef EXPLICIT_PERMUTE
    pMat.Localg_.SetSorted(0);
    pMat.ExpandSymmetric();
    //matrix should already be distributed ?
    pMat.Permute(&Order_.invp[0]);
#endif
    {
      double timeSta = get_time();
      //Icomm buffer;
      Icomm recv_buffer;


      AsyncComms incomingRecv;
      AsyncComms outgoingSend;
      Int numColFirst = std::max(1,iSize_ / np);

      logfileptr->OFS()<<"Starting Send"<<std::endl;

#ifndef NOTRY
      try
#endif
      {
#ifdef EXPLICIT_PERMUTE
        pMat.ToLowerTriangular();
        pMat.SortGraph();

        //    logfileptr->OFS()<<"*********************************"<<std::endl;
        //    pMat.DumpMatlab();
        //    logfileptr->OFS()<<"*********************************"<<std::endl;

        vector<std::pair<size_t,Icomm *> > send_map(np);

        SYMPACK_TIMER_START(DISTRIBUTE_COUNTING);

        Idx FirstLocalCol = pMat.Localg_.vertexDist[iam];
        Idx LastLocalCol = pMat.Localg_.vertexDist[iam+1];
        //        Int snodeCount = 0;


        //gdb_lock();

        if(FirstLocalCol<LastLocalCol){
          for(Int I=1;I<Xsuper_.size();I++){
            Idx fc = Xsuper_[I-1];
            Idx lc = Xsuper_[I]-1;
            Int iWidth = lc-fc+1;
            Int iHeight = UpdateHeight_[I-1];

            Int iDest = this->Mapping_->Map(I-1,I-1);

            //post all the recv and sends
            for(Idx col = fc;col<=lc;col++){
              //corresponding column in the unsorted matrix A
              if(col>= FirstLocalCol && col < LastLocalCol){
                Idx local_col = col - FirstLocalCol+1;//1 based
                send_map[iDest].first += (pMat.Localg_.colptr[local_col]-pMat.Localg_.colptr[local_col-1])*(sizeof(Idx)+sizeof(T)) + sizeof(Idx) + sizeof(Ptr);
              }
            }
          }
        }
        SYMPACK_TIMER_STOP(DISTRIBUTE_COUNTING);


        //        //Resize the local supernodes array
        //        LocalSupernodes_.reserve(snodeCount);


        //        SYMPACK_TIMER_START(DISTRIBUTE_CREATE_SNODES);
        //        if(FirstLocalCol<LastLocalCol){
        //          for(Int I=1;I<Xsuper_.size();I++){
        //            Int iDest = this->Mapping_->Map(I-1,I-1);
        //            //parse the first column to create the supernode structure
        //            if(iam==iDest){
        //              Int fc = Xsuper_[I-1];
        //              Int lc = Xsuper_[I]-1;
        //              Int iWidth = lc-fc+1;
        //              Int iHeight = UpdateHeight_[I-1];
        //              Int nzBlockCnt = numBlk_[I-1];
        //#ifndef ITREE2
        //              globToLocSnodes_.push_back(I-1);
        //#else 
        //              ITree::Interval snode_inter = { I, I, LocalSupernodes_.size() };
        //              globToLocSnodes_.Insert(snode_inter);
        //#endif
        //              SuperNode<T> * newSnode = CreateSuperNode(options_.decomposition,I,fc,lc,iHeight,iSize_,nzBlockCnt);
        //              LocalSupernodes_.push_back(newSnode);
        //            }
        //          }
        //        }
        //        SYMPACK_TIMER_STOP(DISTRIBUTE_CREATE_SNODES);

        //first create the structure of every supernode
        //        {
        //
        //          scope_timer(a,DISTRIBUTE_CREATE_REMOTE_NODES);
        //          std::vector< int > superStructure(all_np);
        //          std::vector< int > rdisplsStructure;
        //          std::vector< int > sSuperStructure(all_np);
        //          std::vector< int > ssizes(all_np,0);
        //          std::vector< int > sdispls(all_np+1,0);
        //          std::vector< int > rsizes(all_np,0);
        //
        //          Int numLocSnode = XsuperDist_[iam+1]-XsuperDist_[iam];
        //          Int firstSnode = XsuperDist_[iam];
        //
        //          for(Int locsupno = 1; locsupno<locXlindx_.size(); ++locsupno){
        //            Idx I = locsupno + firstSnode-1;
        //            Int iDest = this->Mapping_->Map(I-1,I-1);
        //            ssizes[iDest] += 1 + 2*numBlk_[I-1]; //1 for supno index + numBlk startrows + numBlk number of rows
        //          }
        //
        //          sdispls[0] = 0;
        //          std::partial_sum(ssizes.begin(),ssizes.end(),&sdispls[1]);
        //
        //          rdisplsStructure = sdispls;
        //          sSuperStructure.resize(sdispls.back());
        //          for(Int locsupno = 1; locsupno<locXlindx_.size(); ++locsupno){
        //            Idx I = locsupno + firstSnode-1;
        //            Int fc = Xsuper_[I-1];
        //            Int lc = Xsuper_[I]-1;
        //            Int iWidth = lc - fc + 1;
        //            Ptr lfi = locXlindx_[locsupno-1];
        //            Ptr lli = locXlindx_[locsupno]-1;
        //            Int iDest = this->Mapping_->Map(I-1,I-1);
        //            int & tail = rdisplsStructure[iDest];
        //
        //            sSuperStructure[tail++] = I;
        //            //count number of contiguous rows
        //            for(Ptr sidx = lfi; sidx<=lli;sidx++){
        //              Idx iStartRow = locLindx_[sidx-1];
        //              Idx iPrevRow = iStartRow;
        //              Int iContiguousRows = 1;
        //              for(Int idx2 = sidx+1; idx2<=lli;idx2++){
        //                Idx iCurRow = locLindx_[idx2-1];
        //                if(iStartRow == locLindx_[lfi-1]){
        //                  if(iCurRow>iStartRow+iWidth-1){
        //                    //enforce the first block to be a square diagonal block
        //                    break;
        //                  }
        //                }
        //
        //                if(iCurRow==iPrevRow+1){
        //                  sidx++;
        //                  ++iContiguousRows;
        //                  iPrevRow=iCurRow;
        //                }
        //                else{
        //                  break;
        //                }
        //              }
        //
        //              sSuperStructure[tail++] = iStartRow;
        //              sSuperStructure[tail++] = iContiguousRows;
        //            }
        //          }
        //
        //
        //          MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,pMat.comm);
        //          rdisplsStructure[0] = 0;
        //          std::partial_sum(rsizes.begin(),rsizes.end(),&rdisplsStructure[1]);
        //          superStructure.resize(rdisplsStructure.back());
        //
        //          //turn everything into byte sizes
        //          //      for(int p = 0; p<ssizes.size();p++){ ssizes[p]*=sizeof(int); }
        //          //      for(int p = 0; p<rsizes.size();p++){ rsizes[p]*=sizeof(int); }
        //          //      for(int p = 0; p<sdispls.size();p++){ sdispls[p]*=sizeof(int); }
        //          //      for(int p = 0; p<rdisplsStructure.size();p++){ rdisplsStructure[p]*=sizeof(int); }
        //
        //          //Do the alltoallv to get the structures        
        //          MPI_Alltoallv(&sSuperStructure[0], &ssizes[0], &sdispls[0], MPI_INT,
        //              &superStructure[0], &rsizes[0], &rdisplsStructure[0], MPI_INT,
        //              pMat.comm);
        //          //              CommEnv_->MPI_GetComm());
        //
        //          //loop through received structure and create supernodes
        //          for(Int p = 0; p<all_np; p++){
        //            int pos = rdisplsStructure[p];
        //            int end = rdisplsStructure[p+1];
        //            //        int pos = rdisplsStructure[p]/sizeof(int);
        //            //        int end = rdisplsStructure[p+1]/sizeof(int);
        //            while(pos<end){
        //              Int I = superStructure[pos++];
        //              Int nzBlockCnt = numBlk_[I-1];
        //
        //              Int fc = Xsuper_[I-1];
        //              Int lc = Xsuper_[I]-1;
        //              Int iWidth = lc-fc+1;
        //              SuperNode<T> * snode = snodeLocal(I);
        //              if(snode->NZBlockCnt()==0){
        //                for(Int i = 0; i<nzBlockCnt;i++){
        //                  Int iStartRow = superStructure[pos++];
        //                  Int iContiguousRows = superStructure[pos++];
        //                  snode->AddNZBlock(iContiguousRows , iWidth,iStartRow);
        //                }
        //                snode->Shrink();
        //              }
        //            }
        //          } 
        //        }

        //check if the Alltoallv needs to be performed in multiple steps
        size_t total_send_size = 0;
        std::vector<size_t> stotcounts(all_np,0);
        for(Int i = 0; i< send_map.size(); i++){
          total_send_size += send_map[i].first;
          stotcounts[i] = send_map[i].first;
        }

        Icomm * IsendPtr = new Icomm(total_send_size,MPI_REQUEST_NULL);

        std::vector<size_t> spositions(all_np+1,0);
        spositions[0] = 0;
        std::partial_sum(stotcounts.begin(),stotcounts.end(),&spositions[1]);

        SYMPACK_TIMER_SPECIAL_START(serializing); 
        auto baseval = pMat.Localg_.GetBaseval();
        if(FirstLocalCol<LastLocalCol){
          for(Int I=1;I<Xsuper_.size();I++){
            Idx fc = Xsuper_[I-1];
            Idx lc = Xsuper_[I]-1;
            Int iWidth = lc-fc+1;
            Int iHeight = UpdateHeight_[I-1];

            Int iDest = this->Mapping_->Map(I-1,I-1);

            //post all the recv and sends
            for(Idx col = fc;col<=lc;col++){
              //corresponding column in the unsorted matrix A
              if(col>= FirstLocalCol && col < LastLocalCol){
                Idx local_col = col - FirstLocalCol+1;//1 based
                Ptr colbeg = pMat.Localg_.colptr[local_col-1]-baseval;//0-based
                Ptr colend = pMat.Localg_.colptr[local_col]-baseval;//0-based
                Ptr nrows = colend-colbeg;
                IsendPtr->setHead(spositions[iDest]);
                (*IsendPtr)<<(Idx)col<<nrows;
                Serialize( *IsendPtr,  &pMat.Localg_.rowind[colbeg], nrows);
                Serialize( *IsendPtr,  &pMat.nzvalLocal[colbeg], nrows);
                spositions[iDest] = IsendPtr->head;
              }
            }
          }
        }
        SYMPACK_TIMER_SPECIAL_STOP(serializing);


        Icomm * IrecvPtr = new Icomm(0,MPI_REQUEST_NULL);

        spositions[0] = 0;
        std::partial_sum(stotcounts.begin(),stotcounts.end(),&spositions[1]);
        std::function<void(Icomm&,size_t)> resize_lambda =
          [](Icomm & container, size_t sz){
            container.resize(sz); 
            container.head = 0;
          };

        IsendPtr->setHead(0);
        mpi::Alltoallv((*IsendPtr), &stotcounts[0], &spositions[0], MPI_BYTE,
            (*IrecvPtr),fullcomm_, resize_lambda);
        //(*IrecvPtr),CommEnv_->MPI_GetComm(), resize_lambda);

        SYMPACK_TIMER_SPECIAL_START(deserializing);      
        IrecvPtr->setHead(0);
        while(IrecvPtr->head < IrecvPtr->capacity()){ 
          char * buffer = IrecvPtr->back();
          //Deserialize
          Idx col = *(Int*)&buffer[0];
          Int I = SupMembership_[col-1];

          Int iDest = this->Mapping_->Map(I-1,I-1);
          bassert(iam==iDest);

          Int fc = Xsuper_[I-1];
          Int lc = Xsuper_[I]-1;
          Int iWidth = lc-fc+1;

          SuperNode<T> * snode = snodeLocal(I);

          //nrows of column col sent by processor p
          Ptr nrows = *((Ptr*)&buffer[sizeof(Idx)]);
          Idx * rowind = (Idx*)(&buffer[sizeof(Ptr)+sizeof(Idx)]);
          T * nzvalA = (T*)(&buffer[nrows*sizeof(Idx)+sizeof(Ptr)+sizeof(Idx)]);
          //advance in the buffer
          IrecvPtr->setHead(IrecvPtr->head + sizeof(Idx) +sizeof(Ptr) + nrows*(sizeof(Idx)+sizeof(T)));

          //Here, do a linear search instead for the blkidx
          Ptr colbeg = 1;
          Ptr colend = nrows;

          if(colbeg<=colend){
            //sort rowind and nzvals
            //            std::vector<size_t> lperm = sort_permutation(&rowind[colbeg-1],&rowind[colend-1]+1,std::less<Int>());
            //            apply_permutation(&rowind[colbeg-1],&rowind[colend-1]+1,lperm);
            //            apply_permutation(&nzvalA[colbeg-1],&nzvalA[colend-1]+1,lperm);
            Idx firstRow = rowind[colbeg-1];
            Int blkidx = snode->FindBlockIdx(firstRow);
            assert(blkidx!=-1);
            Int blk_nrows = snode->NRows(blkidx);    
            NZBlockDesc * blk_desc = &snode->GetNZBlockDesc(blkidx);

            for(Ptr rowidx = colbeg; rowidx<=colend; ++rowidx){
              Idx row = rowind[rowidx-1];
              while(row>blk_desc->GIndex+blk_nrows-1){
                blkidx++;
                blk_nrows = snode->NRows(blkidx);    
                blk_desc = &snode->GetNZBlockDesc(blkidx);
              }

              Int local_row = row - blk_desc->GIndex + 1;
              Int local_col = col - fc + 1;
              T * nzval = snode->GetNZval(blk_desc->Offset);
              nzval[(local_row-1)*iWidth+local_col-1] = nzvalA[rowidx-1];
            }
          }
        }
        SYMPACK_TIMER_SPECIAL_STOP(deserializing);

        //after the alltoallv, cleanup
        delete IrecvPtr;
        delete IsendPtr;


#else








        std::map<Int,std::pair<size_t,Icomm *> > send_map;
        for(Int p=0;p<all_np;++p){
          send_map[p].first = 0;
        }
        SYMPACK_TIMER_START(DISTRIBUTE_COUNTING);

        Int baseval = pMat.Localg_.GetBaseval();
        Idx FirstLocalCol = pMat.Localg_.vertexDist[iam] + (1 - baseval); //1-based
        Idx LastLocalCol = pMat.Localg_.vertexDist[iam+1] + (1 - baseval); //1-based

        for(Int I=1;I<Xsuper_.size();I++){
          Idx fc = Xsuper_[I-1];
          Idx lc = Xsuper_[I]-1;
          Int iWidth = lc-fc+1;
          Int iHeight = UpdateHeight_[I-1];

          Int iDest = this->Mapping_->Map(I-1,I-1);

          //post all the recv and sends
          for(Idx col = fc;col<=lc;col++){

            //corresponding column in the unsorted matrix A
            Idx orig_col = Order_.perm[col-1];

            if(orig_col>= FirstLocalCol && orig_col < LastLocalCol){
              Ptr nrows = 0;
              Idx local_col = orig_col - FirstLocalCol+1;//1 based
              for(Ptr rowidx = pMat.Localg_.colptr[local_col-1] + (1-baseval); rowidx<pMat.Localg_.colptr[local_col]+(1-baseval); ++rowidx){
                Idx orig_row = pMat.Localg_.rowind[rowidx-1]+(1-baseval);//1-based
                Idx row = Order_.invp[orig_row-1];

                if(row<col){
                  //add the pair (col,row) to processor owning column row
                  Int J = SupMembership_[row-1];
                  Int iDestJ = this->Mapping_->Map(J-1,J-1);
                  send_map[iDestJ].first += sizeof(Idx)+sizeof(Ptr)+1*(sizeof(Idx)+sizeof(T));
                }
                else{nrows++;}
              }
              send_map[iDest].first += sizeof(Idx)+sizeof(Ptr)+nrows*(sizeof(Idx)+sizeof(T));
            }
          }
        }

        SYMPACK_TIMER_STOP(DISTRIBUTE_COUNTING);



        {
          //allocate one buffer for every remote processor
          //compute the send structures
          size_t total_send_size = 0;
          std::vector<size_t> stotcounts(all_np,0);
          for(auto it = send_map.begin(); it!=send_map.end();it++){
            Int iCurDest = it->first;
            size_t & send_bytes = it->second.first;
            stotcounts[iCurDest] = send_bytes;
          }
          //compute send displacements
          std::vector<size_t> spositions(all_np+1,0);
          spositions[0] = 0;
          std::partial_sum(stotcounts.begin(),stotcounts.end(),&spositions[1]);
          total_send_size = spositions.back();


          Icomm * IsendPtr = new Icomm(total_send_size,MPI_REQUEST_NULL);


          //Fill up the send buffer
          SYMPACK_TIMER_SPECIAL_START(serializing);      
          for(Int I=1;I<Xsuper_.size();I++){
            Idx fc = Xsuper_[I-1];
            Idx lc = Xsuper_[I]-1;
            Int iWidth = lc-fc+1;
            Int iHeight = cc_[fc-1];

            Icomm & Isend = *IsendPtr;
            Int iDest = this->Mapping_->Map(I-1,I-1);

#ifdef _DEBUG_
            logfileptr->OFS()<<"Supernode "<<I<<" is owned by P"<<iDest<<std::endl;
#endif

            //Serialize
            for(Idx col = fc;col<=lc;col++){
              Idx orig_col = Order_.perm[col-1];

              if(orig_col>= FirstLocalCol && orig_col < LastLocalCol){
                Ptr nrows = 0;
                Idx local_col = orig_col - FirstLocalCol+1;//1 based

                for(Ptr rowidx = pMat.Localg_.colptr[local_col-1]+(1-baseval); rowidx<pMat.Localg_.colptr[local_col]+(1-baseval); ++rowidx){
                  Idx orig_row = pMat.Localg_.rowind[rowidx-1]+(1-baseval);
                  Idx row = Order_.invp[orig_row-1];


                  if(row<col){
                    //add the pair (col,row) to processor owning column row
                    Int J = SupMembership_[row-1];
                    Int iDestJ = this->Mapping_->Map(J-1,J-1);

                    T val = pMat.nzvalLocal[rowidx-1];

                    //we need to set head to the proper sdispls
                    Isend.setHead(spositions[iDestJ]);
                    Isend<<row;
                    Isend<<(Ptr)1;
                    Isend<<col;
                    Isend<<val;
                    //backup new position for processor iDestJ
                    spositions[iDestJ] = Isend.head;
                  }
                  else{
                    nrows++;
                  }
                }


                //restore head
                Isend.setHead(spositions[iDest]);
                Isend<<col;
                Isend<<nrows;

                for(Ptr rowidx = pMat.Localg_.colptr[local_col-1]+(1-baseval); rowidx<pMat.Localg_.colptr[local_col]+(1-baseval); ++rowidx){
                  Idx orig_row = pMat.Localg_.rowind[rowidx-1]+(1-baseval);
                  Idx row = Order_.invp[orig_row-1];
                  if(row>=col){
                    Isend<<row;
                  }
                }

                for(Ptr rowidx = pMat.Localg_.colptr[local_col-1]+(1-baseval); rowidx<pMat.Localg_.colptr[local_col]+(1-baseval); ++rowidx){
                  Idx orig_row = pMat.Localg_.rowind[rowidx-1]+(1-baseval);
                  Idx row = Order_.invp[orig_row-1];
                  if(row>=col){
                    Isend<<pMat.nzvalLocal[rowidx-1];
                  }
                }

                //backup new position for processor iDest
                spositions[iDest] = Isend.head;


              }
            }
          }
          SYMPACK_TIMER_SPECIAL_STOP(serializing);      


          Icomm * IrecvPtr = new Icomm(0,MPI_REQUEST_NULL);
          spositions[0] = 0;
          std::partial_sum(stotcounts.begin(),stotcounts.end(),&spositions[1]);
          std::function<void(Icomm&,size_t)> resize_lambda =
            [](Icomm & container, size_t sz){
              container.resize(sz); 
              container.head = 0;
            };

          IsendPtr->setHead(0);
          mpi::Alltoallv((*IsendPtr), &stotcounts[0], &spositions[0], MPI_BYTE,
              (*IrecvPtr),fullcomm_, resize_lambda);


          //Need to parse the structure sent from the processor owning the first column of the supernode

          SYMPACK_TIMER_SPECIAL_START(deserializing);      
          IrecvPtr->setHead(0);
          while(IrecvPtr->head < IrecvPtr->capacity()){ 

            char * buffer = IrecvPtr->back();
            //Deserialize
            Idx col = *(Idx*)&buffer[0];
            //nrows of column col sent by processor p
            Ptr nrows = *((Ptr*)&buffer[sizeof(Idx)]);
            Idx * rowind = (Idx*)(&buffer[sizeof(Ptr)+sizeof(Idx)]);
            T * nzvalA = (T*)(&buffer[(1+nrows)*sizeof(Idx)+sizeof(Ptr)]);

            //advance in the buffer
            IrecvPtr->setHead(IrecvPtr->head + sizeof(Idx) +sizeof(Ptr) + nrows*(sizeof(Idx)+sizeof(T)));


            Int I = SupMembership_[col-1];
            //do a global to local mapping


            Int iDest = this->Mapping_->Map(I-1,I-1);
            bassert(iam==iDest);

            Int fc = Xsuper_[I-1];
            Int lc = Xsuper_[I]-1;
            Int iWidth = lc-fc+1;

            SuperNode<T> * snode = snodeLocal(I);


            //Here, do a linear search instead for the blkidx

            Ptr colbeg = 1;
            Ptr colend = nrows;
            //            for(Ptr rowidx = colbeg; rowidx<=colend; ++rowidx){
            //              Idx row = rowind[rowidx-1];
            //              logfileptr->OFS()<<row<<" ";
            //            }
            //              logfileptr->OFS()<<std::endl;

            if(colbeg<=colend){
              //sort rowind and nzvals

#if 1 
              std::vector<size_t> lperm = sort_permutation(&rowind[colbeg-1],&rowind[colend-1]+1,std::less<Idx>());
              apply_permutation(&rowind[colbeg-1],&rowind[colend-1]+1,lperm);
              apply_permutation(&nzvalA[colbeg-1],&nzvalA[colend-1]+1,lperm);
              Idx firstRow = rowind[colbeg-1];
              Int blkidx = snode->FindBlockIdx(firstRow);
              assert(blkidx!=-1);
              Int blk_nrows = snode->NRows(blkidx);    
              NZBlockDesc * blk_desc = &snode->GetNZBlockDesc(blkidx);
#else
              Idx prevRow = 0;
              Int blkidx = 0;
              NZBlockDesc * blk_desc = &snode->GetNZBlockDesc(blkidx);
              Int blk_nrows = snode->NRows(blkidx);
#endif

              for(Ptr rowidx = colbeg; rowidx<=colend; ++rowidx){
                Idx row = rowind[rowidx-1];

#if 1
                while(row>blk_desc->GIndex+blk_nrows-1){
                  blkidx++;
                  blk_nrows = snode->NRows(blkidx);    
                  blk_desc = &snode->GetNZBlockDesc(blkidx);
                }
#else
                if(row!=prevRow){
                  if(row>blk_desc->GIndex+blk_nrows || row<blk_desc->GIndex){
                    blkidx = snode->FindBlockIdx(row);
                    blk_desc = &snode->GetNZBlockDesc(blkidx);
                    blk_nrows = snode->NRows(blkidx);
                  }
                  prevRow=row;
                }
#endif

                Int local_row = row - blk_desc->GIndex + 1;
                Int local_col = col - fc + 1;
                T * nzval = snode->GetNZval(blk_desc->Offset);
                nzval[(local_row-1)*iWidth+local_col-1] = nzvalA[rowidx-1];
              }
              //for(Int rowidx = colbeg; rowidx<=colend; ++rowidx){
              //  Int row = rowind[rowidx-1];
              //  Int blkidx = snode->FindBlockIdx(row);
              //  assert(blkidx!=-1);
              //  NZBlockDesc & blk_desc = snode->GetNZBlockDesc(blkidx);
              //  Int local_row = row - blk_desc.GIndex + 1;
              //  Int local_col = col - fc + 1;
              //  T * nzval = snode->GetNZval(blk_desc.Offset);
              //  nzval[(local_row-1)*iWidth+local_col-1] = nzvalA[rowidx-1];
              //}
            }
          }
          SYMPACK_TIMER_SPECIAL_STOP(deserializing);      

          //after the alltoallv, cleanup
          delete IrecvPtr;
          delete IsendPtr;
        }
#endif




        //std::vector<upcxx::global_ptr<SuperNodeDesc > > localFactors;
        remoteFactors_.resize(Xsuper_.size()-1);
        std::fill((char*)&remoteFactors_[0],(char*)&remoteFactors_[0]+remoteFactors_.size()*sizeof(std::tuple<upcxx::global_ptr<SuperNodeDesc>,Int> ),0);
        //logfileptr->OFS()<<"My usable global memory size is: "<<upcxx::my_usable_global_memory_size()<<std::endl;

        //      if(iam<np){



        for(Int I=1;I<Xsuper_.size();I++){
          Int iDest = this->Mapping_->Map(I-1,I-1);
          //parse the first column to create the supernode structure
          if(iam==iDest){
            SuperNode<T> * newSnode = snodeLocal(I);
            SuperNodeDesc * meta = newSnode->GetMeta();
            remoteFactors_[I-1] = std::make_tuple( upcxx::global_ptr<SuperNodeDesc>( meta ), meta->blocks_cnt_) ;
          }
        }

        //        MPI_Allreduce( MPI_IN_PLACE, &remoteFactors_[0], remoteFactors_.size()*sizeof( std::tuple<upcxx::global_ptr<SuperNodeDesc>,Int> ), MPI_BYTE, MPI_BOR, CommEnv_->MPI_GetComm());
        MPI_Allreduce( MPI_IN_PLACE, &remoteFactors_[0], remoteFactors_.size()*sizeof( std::tuple<upcxx::global_ptr<SuperNodeDesc>,Int> ), MPI_BYTE, MPI_BOR, fullcomm_);
        //      }
      }
#ifndef NOTRY
      catch(const std::bad_alloc& e){
        std::cout << "Allocation failed: " << e.what() << '\n';
        abort();
      }
#endif
      logfileptr->OFS()<<"Send Done"<<std::endl;

      double timeStop = get_time();
      if(iam==0){
        std::cout<<"Distribution time: "<<timeStop - timeSta<<std::endl;
      }
    }
    //    MPI_Barrier(CommEnv_->MPI_GetComm());
    MPI_Barrier(pMat.comm);
  }




  template <typename T> symPACKMatrix<T>::symPACKMatrix(){
    CommEnv_=NULL;
    Mapping_ = NULL;
    Balancer_ = NULL;
    scheduler_ = NULL;
    scheduler2_ = NULL;
    //Local_=NULL;
    //Global_=NULL;
    //isGlobStructAllocated_ = false;
    non_workcomm_ = MPI_COMM_NULL;
    fullcomm_ = MPI_COMM_NULL;
    if(!upcxx::is_init()){
      upcxx::init(NULL,NULL);
    }

    if(logfileptr==NULL){
      logfileptr = new LogFile(upcxx::myrank(),false);
      logfileptr->OFS()<<"********* LOGFILE OF P"<<upcxx::myrank()<<" *********"<<std::endl;
      logfileptr->OFS()<<"**********************************"<<std::endl;
    }
  }

  template <typename T> symPACKMatrix<T>::symPACKMatrix(DistSparseMatrix<T> & pMat, symPACKOptions & options ):symPACKMatrix(){
    Init(pMat, options);
  }

  template <typename T> symPACKMatrix<T>::~symPACKMatrix(){
    for(Int i=0;i<LocalSupernodes_.size();++i){
      delete LocalSupernodes_[i];
    }

    for(Int i=0;i<Contributions_.size();++i){
      delete Contributions_[i];
    }

#ifdef PREALLOC_IRECV
    //    for(auto it = availRecvBuffers_.begin(); it != availRecvBuffers_.end();++it){
    //      delete *it;
    //    }
#endif

    if(this->scheduler_!=NULL){
      delete this->scheduler_;
    }

    if(this->scheduler2_!=NULL){
      delete this->scheduler2_;
    }

    if(this->Mapping_!=NULL){
      delete this->Mapping_;
    }

    if(this->Balancer_!=NULL){
      delete this->Balancer_;
    }

    if(CommEnv_!=NULL){
      delete CommEnv_;
    }

    //  if(Local_!=NULL){
    //    delete Local_;
    //  }
    //
    //  if(Global_!=NULL){
    //    delete Global_;
    //  }

#ifdef MULTITHREADING
    std::for_each(superLocks_.begin(), superLocks_.end(), [] (omp_lock_t & lock) { omp_destroy_lock(&lock);});
#endif

    if(this->non_workcomm_ != MPI_COMM_NULL){
      MPI_Comm_free(&non_workcomm_);
    }

    if(this->fullcomm_ != MPI_COMM_NULL){
      MPI_Comm_free(&fullcomm_);
    }

  }


  //returns the 1-based index of supernode id global in the local supernode array
  template <typename T> Int symPACKMatrix<T>::snodeLocalIndex(Int global){
#ifndef ITREE2
    auto it = std::lower_bound(globToLocSnodes_.begin(),globToLocSnodes_.end(),global);
    return it - globToLocSnodes_.begin();
#else
    ITree::Interval * ptr = globToLocSnodes_.IntervalSearch(global,global);
    assert(ptr!=NULL);
    return ptr->block_idx;
#endif
  }

  //returns a reference to  a local supernode with id global
  template <typename T> SuperNode<T> * symPACKMatrix<T>::snodeLocal(Int global){
    Int iLocal = snodeLocalIndex(global);
    return LocalSupernodes_[iLocal -1];
  }

  template <typename T> 
    template< class Alloc>
    SuperNode<T,Alloc> * symPACKMatrix<T>::snodeLocal(Int global, std::vector<SuperNode<T,Alloc> *> & snodeColl){
      Int iLocal = snodeLocalIndex(global);
      return snodeColl[iLocal -1];
    }


  template <typename T> 
    void symPACKMatrix<T>::findSupernodes(ETree& tree, Ordering & aOrder, std::vector<Int> & cc,std::vector<Int> & supMembership, std::vector<Int> & xsuper, Int maxSize ){
      SYMPACK_TIMER_START(FindSupernodes);
      Int size = iSize_;
      //TODO: tree order cc supmembership xsuper are all members of the class. no need for argument
      supMembership.resize(size);

      Int nsuper = 1;
      Int supsize = 1;
      supMembership[0] = 1;

      for(Int i =2; i<=size;i++){
        Int prev_parent = tree.PostParent(i-2);
        if(prev_parent == i){
          if(cc[i-2] == cc[i-1]+1 ) {
            if(supsize<maxSize || maxSize==0){
              ++supsize;
              supMembership[i-1] = nsuper;
              continue;
            }
          }
        }

        nsuper++;
        supsize = 1;
        supMembership[i-1] = nsuper;
      }

      xsuper.resize(nsuper+1);
      Int lstsup = nsuper+1;
      for(Int i = size; i>=1;--i){
        Int ksup = supMembership[i-1];
        if(ksup!=lstsup){
          xsuper[lstsup-1] = i + 1; 
        }
        lstsup = ksup;
      }
      xsuper[0]=1;
      SYMPACK_TIMER_STOP(FindSupernodes);
    }



  template <typename T> 
    void symPACKMatrix<T>::getLColRowCount(SparseMatrixGraph & sgraph, std::vector<Int> & cc, std::vector<Int> & rc){
      //The tree need to be postordered
      if(!ETree_.IsPostOrdered()){
        ETree_.PostOrderTree(Order_);
      }

      if(iam == 0 && (!sgraph.IsExpanded() ) ){
        throw std::logic_error( "SparseMatrixGraph must be expanded\n" );
      }

      sgraph.SetBaseval(1);
      sgraph.SetKeepDiag(1);

      SYMPACK_TIMER_START(GetColRowCount_Classic);

      Int size = sgraph.size;

      if(iam==0){
        cc.resize(size);
        rc.resize(size);
        std::vector<Idx> level(size+1);
        std::vector<Int> weight(size+1);
        std::vector<Idx> fdesc(size+1);
        std::vector<Idx> nchild(size+1);
        std::vector<Idx> set(size);
        std::vector<Idx> prvlf(size);
        std::vector<Idx> prvnbr(size);

        Idx xsup = 1;
        level[0] = 0;
        for(Idx k = size; k>=1; --k){
          rc[k-1] = 1;
          cc[k-1] = 0;
          set[k-1] = k;
          prvlf[k-1] = 0;
          prvnbr[k-1] = 0;
          level[k] = level[ETree_.PostParent(k-1)] + 1;
          weight[k] = 1;
          fdesc[k] = k;
          nchild[k] = 0;
        }

        nchild[0] = 0;
        fdesc[0] = 0;
        for(Idx k =1; k<size; ++k){
          Idx parent = ETree_.PostParent(k-1);
          weight[parent] = 0;
          ++nchild[parent];
          Idx ifdesc = fdesc[k];
          if  ( ifdesc < fdesc[parent] ) {
            fdesc[parent] = ifdesc;
          }
        }

        for(Idx lownbr = 1; lownbr<=size; ++lownbr){
          Int lflag = 0;
          Idx ifdesc = fdesc[lownbr];
          Idx oldnbr = Order_.perm[lownbr-1];
          Ptr jstrt = sgraph.colptr[oldnbr-1];
          Ptr jstop = sgraph.colptr[oldnbr] - 1;
          //           -----------------------------------------------
          //           for each ``high neighbor'', hinbr of lownbr ...
          //           -----------------------------------------------
          for(Ptr j = jstrt; j<=jstop;++j){
            Idx hinbr = sgraph.rowind[j-1];
            hinbr = Order_.invp[hinbr-1];
            if  ( hinbr > lownbr )  {
              if  ( ifdesc > prvnbr[hinbr-1] ) {
                //                       -------------------------
                //                       increment weight[lownbr].
                //                       -------------------------
                ++weight[lownbr];
                Idx pleaf = prvlf[hinbr-1];
                //                       -----------------------------------------
                //                       if hinbr has no previous ``low neighbor'' 
                //                       then ...
                //                       -----------------------------------------
                if  ( pleaf == 0 ) {
                  //                           -----------------------------------------
                  //                           ... accumulate lownbr-->hinbr path length 
                  //                               in rowcnt[hinbr].
                  //                           -----------------------------------------
                  rc[hinbr-1] += level[lownbr] - level[hinbr];
                }
                else{
                  //                           -----------------------------------------
                  //                           ... otherwise, lca <-- find[pleaf], which 
                  //                               is the least common ancestor of pleaf 
                  //                               and lownbr.
                  //                               (path halving.)
                  //                           -----------------------------------------
                  Idx last1 = pleaf;
                  Idx last2 = set[last1-1];
                  Idx lca = set[last2-1];
                  while(lca != last2){
                    set[last1-1] = lca;
                    last1 = lca;
                    last2 = set[last1-1];
                    lca = set[last2-1];
                  }
                  //                           -------------------------------------
                  //                           accumulate pleaf-->lca path length in 
                  //                           rowcnt[hinbr].
                  //                           decrement weight(lca).
                  //                           -------------------------------------
                  rc[hinbr-1] += level[lownbr] - level[lca];
                  --weight[lca];
                }
                //                       ----------------------------------------------
                //                       lownbr now becomes ``previous leaf'' of hinbr.
                //                       ----------------------------------------------
                prvlf[hinbr-1] = lownbr;
                lflag = 1;
              }
              //                   --------------------------------------------------
              //                   lownbr now becomes ``previous neighbor'' of hinbr.
              //                   --------------------------------------------------
              prvnbr[hinbr-1] = lownbr;
            }
          }
          //           ----------------------------------------------------
          //           decrement weight ( parent[lownbr] ).
          //           set ( p[lownbr] ) <-- set ( p[lownbr] ) + set[xsup].
          //           ----------------------------------------------------
          Idx parent = ETree_.PostParent(lownbr-1);
          --weight[parent];

          if  ( lflag == 1  || nchild[lownbr] >= 2 ) {
            xsup = lownbr;
          }
          set[xsup-1] = parent;
        }

        for(Int k = 1; k<=size; ++k){
          Int temp = cc[k-1] + weight[k];
          cc[k-1] = temp;
          Int parent = ETree_.PostParent(k-1);
          if  ( parent != 0 ) {
            cc[parent-1] += temp;
          }
        }
      }

      if(iam!=0){
        cc.resize(size);
        rc.resize(size);
      }
      //Broadcast to everyone 
      MPI_Bcast(&cc[0],size*sizeof(Int),MPI_BYTE,0,fullcomm_);
      MPI_Bcast(&rc[0],size*sizeof(Int),MPI_BYTE,0,fullcomm_);


      SYMPACK_TIMER_STOP(GetColRowCount_Classic);
    }


  template <typename T> 
    void symPACKMatrix<T>::getLColRowCount(std::vector<Int> & cc, std::vector<Int> & rc){
      abort();
      //The tree need to be postordered
      if(!ETree_.IsPostOrdered()){
        ETree_.PostOrderTree(Order_);
      }

      graph_.ExpandSymmetric();

      SYMPACK_TIMER_START(GetColRowCount_Classic);

      Int size = graph_.size;
      cc.resize(size);
      rc.resize(size);


      upcxx::shared_array<Idx> rcDist(size);
      upcxx::shared_array<Idx> level(size+1);
      upcxx::shared_array<Idx> nchild(size+1);
      upcxx::shared_array<Int> weight(size+1);
      //ARRAY OF LENGTH NEQNS+1 CONTAINING THE FIRST (I.E., LOWEST-NUMBERED) DESCENDANT.
      upcxx::shared_array<Idx> fdesc(size+1);
      //ARRAY OF LENGTH NEQNS USED TO MAINTAIN THE DISJOINT SETS (I.E., SUBTREES).
      upcxx::shared_array<Idx> set(size);
      //ARRAY OF LENGTH NEQNS USED TO RECORD THE PREVIOUS LEAF OF EACH ROW SUBTREE.
      upcxx::shared_array<Idx> prvlf(size);
      //ARRAY OF LENGTH NEQNS USED TO RECORD THE PREVIOUS ``LOWER NEIGHBOR'' OF EACH NODE.
      upcxx::shared_array<Idx> prvnbr(size);

      upcxx::shared_var<Idx> xsup = 1;
      level[0] = 0;
      for(Idx k = size; k>=1; --k){
        MPI_Barrier(CommEnv_->MPI_GetComm());
        cc[k-1] = 0;
        if( (k-1) % np == iam){
          upcxx::global_ptr<Idx> gptr = &set[k-1]; assert(gptr.where()==iam);
          rcDist[k-1] = 1;
          set[k-1] = k;
          prvlf[k-1] = 0;
          prvnbr[k-1] = 0;
        }
        if( (k) % np == iam){
          level[k] = level[ETree_.PostParent(k-1)] + 1;
          weight[k] = 1;
          fdesc[k] = k;
          nchild[k] = 0;
        }
      }

      nchild[0] = 0;
      fdesc[0] = 0;
      for(Idx k =1; k<size; ++k){
        MPI_Barrier(CommEnv_->MPI_GetComm());
        Idx parent = ETree_.PostParent(k-1);
        if( (parent) % np == iam){
          upcxx::global_ptr<Int> gptr = &weight[parent]; assert(gptr.where()==iam);
          weight[parent] = 0;
          nchild[parent] = nchild[parent] + 1;
        }
        if( (k) % np == iam){
          Idx ifdesc = fdesc[k];
          if  ( ifdesc < fdesc[parent] ) {
            fdesc[parent] = ifdesc;
          }
        }
      }

      MPI_Barrier(CommEnv_->MPI_GetComm());
#if 0
      logfileptr->OFS()<<"nchild: ";
      for(Idx k =0; k<=size; ++k){
        logfileptr->OFS()<<nchild[k]<<" ";
      }
      logfileptr->OFS()<<std::endl;
      logfileptr->OFS()<<"fdesc: ";
      for(Idx k =0; k<=size; ++k){
        logfileptr->OFS()<<fdesc[k]<<" ";
      }
      logfileptr->OFS()<<std::endl;
      logfileptr->OFS()<<"weight: ";
      for(Idx k =0; k<=size; ++k){
        logfileptr->OFS()<<weight[k]<<" ";
      }
      logfileptr->OFS()<<std::endl;
      logfileptr->OFS()<<"level: ";
      for(Idx k =0; k<=size; ++k){
        logfileptr->OFS()<<level[k]<<" ";
      }
      logfileptr->OFS()<<std::endl;
      logfileptr->OFS()<<"prvnbr: ";
      for(Idx k =0; k<size; ++k){
        logfileptr->OFS()<<prvnbr[k]<<" ";
      }
      logfileptr->OFS()<<std::endl;
      logfileptr->OFS()<<"prvlf: ";
      for(Idx k =0; k<size; ++k){
        logfileptr->OFS()<<prvlf[k]<<" ";
      }
      logfileptr->OFS()<<std::endl;
      logfileptr->OFS()<<"set: ";
      for(Idx k =0; k<size; ++k){
        logfileptr->OFS()<<set[k]<<" ";
      }
      logfileptr->OFS()<<std::endl;
      logfileptr->OFS()<<"rcDist: ";
      for(Idx k =0; k<size; ++k){
        logfileptr->OFS()<<rcDist[k]<<" ";
      }
      logfileptr->OFS()<<std::endl;
      MPI_Barrier(CommEnv_->MPI_GetComm());
#endif

      //foreach column in the permuted matrix
      for(Idx lownbr = 1; lownbr<=size; ++lownbr){
        MPI_Barrier(CommEnv_->MPI_GetComm());
        Int lflag = 0;

        // get the original column index
        Idx orig_lownbr = Order_.perm[lownbr-1];

        //if it is a local column
        if(graph_.vertexDist[iam] <= orig_lownbr && orig_lownbr<graph_.vertexDist[iam+1]){
          Idx ifdesc = fdesc[lownbr];
          //get the local column index in the distributed graph
          Idx local_lownbr = orig_lownbr - graph_.LocalFirstVertex() +1;
          Ptr jstrt = graph_.colptr[local_lownbr-1];
          Ptr jstop = graph_.colptr[local_lownbr] - 1;
          //           -----------------------------------------------
          //           for each ``high neighbor'', hinbr of lownbr ...
          //           -----------------------------------------------
          for(Ptr j = jstrt; j<=jstop;++j){
            Idx orig_hinbr = graph_.rowind[j-1];
            Idx hinbr = Order_.invp[orig_hinbr-1];
            if  ( hinbr > lownbr )  {
              if  ( ifdesc > prvnbr[hinbr-1] ) {
                //                       -------------------------
                //                       increment weight[lownbr].
                //                       -------------------------
                weight[lownbr] = weight[lownbr]+1;
                Idx pleaf = prvlf[hinbr-1];
                //                       -----------------------------------------
                //                       if hinbr has no previous ``low neighbor'' 
                //                       then ...
                //                       -----------------------------------------
                if  ( pleaf == 0 ) {
                  //                           -----------------------------------------
                  //                           ... accumulate lownbr-->hinbr path length 
                  //                               in rowcnt[hinbr].
                  //                           -----------------------------------------
                  rcDist[hinbr-1] = rcDist[hinbr-1] +  level[lownbr] - level[hinbr];
                }
                else{
                  //                           -----------------------------------------
                  //                           ... otherwise, lca <-- find[pleaf], which 
                  //                               is the least common ancestor of pleaf 
                  //                               and lownbr.
                  //                               (path halving.)
                  //                           -----------------------------------------

#if 0
                  std::map<Idx,Idx> mySet;

                  Idx last1 = pleaf;
                  Idx last2 = set[last1-1];
                  Idx lca = set[last2-1];
                  mySet[last1-1] = last2;
                  mySet[last2-1] = lca;
                  while(lca != last2){
                    set[last1-1] = lca;
                    mySet[last2-1] = lca;
                    last1 = lca;
                    last2 = set[last1-1];
                    lca = set[last2-1];
                    mySet[last1-1] = last2;
                    mySet[last2-1] = lca;
                  }

                  logfileptr->OFS()<<"local set: ";
                  for(auto it = mySet.begin();it!=mySet.end();it++){
                    logfileptr->OFS()<<it->first<<" -> "<<it->second/*<<" ("<<set[it->first-1]<<")"*/<<" ";
                  }
                  logfileptr->OFS()<<std::endl;


#else
                  Idx last1 = pleaf;
                  Idx last2 = set[last1-1];
                  Idx lca = set[last2-1];
                  while(lca != last2){
                    set[last1-1] = lca;
                    last1 = lca;
                    last2 = set[last1-1];
                    lca = set[last2-1];
                  }
#endif

                  //                           -------------------------------------
                  //                           accumulate pleaf-->lca path length in 
                  //                           rowcnt[hinbr].
                  //                           decrement weight(lca).
                  //                           -------------------------------------
                  rcDist[hinbr-1] = rcDist[hinbr-1] + level[lownbr] - level[lca];
                  weight[lca] = weight[lca] -1;
                }
                //                       ----------------------------------------------
                //                       lownbr now becomes ``previous leaf'' of hinbr.
                //                       ----------------------------------------------
                prvlf[hinbr-1] = lownbr;
                lflag = 1;
              }
              //                   --------------------------------------------------
              //                   lownbr now becomes ``previous neighbor'' of hinbr.
              //                   --------------------------------------------------
              prvnbr[hinbr-1] = lownbr;
            }
          }
          //           ----------------------------------------------------
          //           decrement weight ( parent[lownbr] ).
          //           set ( p[lownbr] ) <-- set ( p[lownbr] ) + set[xsup].
          //           ----------------------------------------------------
          Idx parent = ETree_.PostParent(lownbr-1);
          weight[parent] = weight[parent]-1;


          //merge the sets
          if  ( lflag == 1  || nchild[lownbr] >= 2 ) {
            xsup = lownbr;
          }
          set[xsup-1] = parent;
        }
      }


      MPI_Barrier(CommEnv_->MPI_GetComm());


      for(Int k = 1; k<=size; ++k){
        Int temp = cc[k-1] + weight[k];
        cc[k-1] = temp;
        Int parent = ETree_.PostParent(k-1);
        if  ( parent != 0 ) {
          cc[parent-1] += temp;
        }
      }

      for(Int k = 1; k<=size; ++k){
        rc[k-1] = rcDist[k-1];
      }
      MPI_Barrier(CommEnv_->MPI_GetComm());

      //      logfileptr->OFS()<<"column counts "<<cc<<std::endl;

      SYMPACK_TIMER_STOP(GetColRowCount_Classic);
    }





  template <typename T> 
    void symPACKMatrix<T>::relaxSupernodes(ETree& tree, std::vector<Int> & cc,std::vector<Int> & supMembership, std::vector<Int> & xsuper, RelaxationParameters & params ){
      //todo tree cc supmembership xsuper and relax params are members, no need for arguments
      Int nsuper = xsuper.size()-1;

      DisjointSet sets;
      sets.Initialize(nsuper);
      std::vector<Int> ncols(nsuper);
      std::vector<Int> zeros(nsuper);
      std::vector<Int> newCC(nsuper);
      for(Int ksup=nsuper;ksup>=1;--ksup){
        Int cset = sets.makeSet(ksup);
        sets.Root(cset-1)=ksup;

        Int fstcol = xsuper[ksup-1];
        Int lstcol = xsuper[ksup]-1;
        Int width = lstcol - fstcol +1;
        Int length = cc[fstcol-1];
        ncols[ksup-1] = width;
        zeros[ksup-1] = 0;
        newCC[ksup-1] = length;
      }



      for(Int ksup=nsuper;ksup>=1;--ksup){
        Int fstcol = xsuper[ksup-1];
        Int lstcol = xsuper[ksup]-1;
        Int width = ncols[ksup-1];
        Int length = cc[fstcol-1];

        Int parent_fc = tree.PostParent(lstcol-1);
        if(parent_fc!=0){
          Int parent_snode = supMembership[parent_fc-1];
          Int pset = sets.find(parent_snode);
          parent_snode = sets.Root(pset-1);

          bool merge = (parent_snode == ksup+1);


          if(merge){
            Int parent_width = ncols[parent_snode-1];

            Int parent_fc = xsuper[parent_snode-1];
            Int totzeros = zeros[parent_snode-1];
            Int merged_snode_size = width + parent_width;

            //TODO Shift the loops to make it less obvious that this is coming from cholmod
            //TODO rename nrelax0 as maxSnodeSize
            //TODO rename nrelax1 as ???
            //TODO rename nrelax2 as ???
            //TODO rename zrelax0 as ???
            //TODO rename zrelax1 as ???
            //TODO rename zrelax2 as ???
            merge = false;

            //Supernode is extremely small > merge it
            //TODO TRY THIS
            if(merged_snode_size <=params.maxSize || params.maxSize==0){
              if(merged_snode_size <= params.nrelax0){
                merge = true;
              }
              else if(merged_snode_size <=params.maxSize || params.maxSize==0){
                double nnzchild = cc[fstcol-1];
                double nnzparent = cc[parent_fc-1];
                double xnewzeros = width * (nnzparent + width  - nnzchild);

                //The merge is not creating extra fill=in, proceed safely
                if(xnewzeros == 0){
                  merge = true;
                }
                else{
                  //candidate merged supernode characteristics
                  double xtotzeros = (double)totzeros + xnewzeros;
                  double xmerged_snode_size = (double) merged_snode_size;
                  //new number of nz
                  double xtotsize = (xmerged_snode_size * (xmerged_snode_size+1)/2) + xmerged_snode_size * (nnzparent - parent_width);
                  //percentage of explicit zeros
                  double z = xtotzeros / xtotsize;

                  Int totsize = (merged_snode_size * (merged_snode_size+1)/2) + merged_snode_size * ((Int)nnzparent - parent_width);
                  totzeros += (Int)xnewzeros;

                  //check that we will not have Integer overflow issues with the Ptr type
                  if (xtotsize * sizeof(double)< std::numeric_limits<Ptr>::max()){
                    if (merged_snode_size <= params.nrelax1 && z < params.zrelax0){
                      merge = true;
                    }
                    else if (merged_snode_size <= params.nrelax2 && z < params.zrelax1){
                      merge = true;
                    }
                    else if (z<params.zrelax2){
                      merge = true;
                    }
                  }
                  // merge = ((merged_snode_size <= params.nrelax1 && z < params.zrelax0) 
                  //     || (merged_snode_size <= params.nrelax2 && z < params.zrelax1)
                  //     || (z<params.zrelax2)) &&
                  //   (xtotsize < std::numeric_limits<Int>::max() / sizeof(double));
                }

              }
            }

            //Merge the two supernodes
            if(merge){
              ncols[ksup-1] += ncols[parent_snode-1]; 
              zeros[ksup-1] = totzeros;
              newCC[ksup-1] = width + newCC[parent_snode-1];
              sets.Union(ksup,parent_snode,ksup);
            }
          } 

        }
      }

      std::vector<Int> relXSuper(nsuper+1);
      Int nrSuper = 0;
      for(Int ksup=1;ksup<=nsuper;++ksup){
        Int kset = sets.find(ksup);
        if(ksup == sets.Root(kset-1)){
          Int fstcol = xsuper[ksup-1];
          relXSuper[nrSuper] = fstcol;
          newCC[nrSuper] = newCC[ksup-1];
          ++nrSuper;
        }
      }
      relXSuper[nrSuper] = xsuper[nsuper];
      relXSuper.resize(nrSuper+1);

      for(Int ksup=1;ksup<=nrSuper;++ksup){
        Int fstcol = relXSuper[ksup-1];
        Int lstcol = relXSuper[ksup]-1;
        for(Int col = fstcol; col<=lstcol;++col){
          supMembership[col-1] = ksup;
          cc[col-1] = newCC[ksup-1] - (col-fstcol);
        }
      }


      xsuper = relXSuper;
      ///      //adjust the column counts
      ///      for(Int col=i-2;col>=i-supsize;--col){
      ///        cc[col-1] = cc[col]+1;
      ///      }


    }


  template <typename T> 
    void symPACKMatrix<T>::symbolicFactorizationRelaxedDist(std::vector<Int> & cc){
      scope_timer(a,SymbolicFactorization);
      Int size = iSize_;
      ETree& tree = ETree_;
      Ordering & aOrder = Order_;
      DistSparseMatrixGraph & graph = graph_;
      std::vector<Int> & xsuper = Xsuper_;
      std::vector<Int> & SupMembership = SupMembership_;
      PtrVec & xlindx = locXlindx_;
      IdxVec & lindx = locLindx_;
      MPI_Comm & comm = graph_.comm;//CommEnv_->MPI_GetComm();


      //permute the graph
#ifdef EXPLICIT_PERMUTE
      DistSparseMatrixGraph & pGraph = graph;
#else
      DistSparseMatrixGraph pGraph = graph;
#endif

#if 0
      {
        double tstart = get_time();
        pGraph.Permute(&Order_.invp[0]);
        double tstop = get_time();
        logfileptr->OFS()<<"Permute time: "<<tstop-tstart<<std::endl;
      }



      //recompute xsuper by splitting some supernodes
      {

        //  gdb_lock();

        Idx colPerProc = size / np;
        Int numSplits = 0;
        for(Int snode = 1; snode<=xsuper.size()-1;snode++){
          Idx fstcol = xsuper[snode-1];
          Idx lstcol = xsuper[snode]-1;
          //check if these two columns are on the same processor
          Idx pOwnerFirst = std::min((Idx)np-1, (fstcol-1) / colPerProc);
          Idx pOwnerLast = std::min((Idx)np-1, (lstcol-1) / colPerProc);

          if(pOwnerFirst!=pOwnerLast){
            numSplits += pOwnerLast-pOwnerFirst;
          }
        }

        std::vector<Int> newSnodes(np,0);
        std::vector<Int> newXsuper(xsuper.size()+numSplits);
        Int pos = 0;
        for(Int snode = 1; snode<=xsuper.size()-1;snode++){
          Idx fstcol = xsuper[snode-1];
          Idx lstcol = xsuper[snode]-1;
          //check if these two columns are on the same processor
          Idx pOwnerFirst = std::min((Idx)np-1, (fstcol-1) / colPerProc);
          Idx pOwnerLast = std::min((Idx)np-1, (lstcol-1) / colPerProc);

          newXsuper[pos++] = xsuper[snode-1];
          if(pOwnerFirst!=pOwnerLast){
            for(Idx p = pOwnerFirst; p<pOwnerLast;p++){
              Idx curLstcol = (p+1)*colPerProc+1;//1-based



              //        assert(ETree_.PostParent(curLstcol-2)==curLstcol);
              newXsuper[pos++] = curLstcol;
              newSnodes[p+1]++;
            }
          }
        }
        newXsuper[pos++]=size+1;

        xsuper = newXsuper;
        for(Int snode = 1; snode<=xsuper.size()-1;snode++){
          Int fstcol = xsuper[snode-1];
          Int lstcol = xsuper[snode]-1;
          for(Int col = fstcol; col<=lstcol;++col){
            SupMembership[col-1] = snode;
          }
        }



        //recompute XsuperDist_
        for(int p =1; p<np; p++){
          XsuperDist_[p]= XsuperDist_[p-1] + supPerProc + newSnodes[p-1];
        }
        XsuperDist_[np] = Xsuper_.size();


        //update cc by recomputing the merged structure



      }
#else
      {

#ifndef EXPLICIT_PERMUTE
        {
          double tstart = get_time();

          //recompute vertexDist based on XsuperDist
          //      std::vector<Idx> newVertexDist(np+1);
          //      for(int p = 0; p < np; p++){
          //        newVertexDist[p] = Xsuper_[XsuperDist_[p]-1] + (pGraph.GetBaseval()-1);
          //      }
          //      newVertexDist[np] = pGraph.size+1;

          std::vector<Idx> newVertexDist;
          {
            newVertexDist.resize(all_np+1,0);
            newVertexDist[all_np] = pGraph.size+1;
            for(int p = 0; p < all_np; p++){
              Int S = XsuperDist_[p];
              newVertexDist[p] = Xsuper_[S-1];
            }
          }



          pGraph.Permute(&Order_.invp[0],&newVertexDist[0]);
          double tstop = get_time();


          logfileptr->OFS()<<"Permute time: "<<tstop-tstart<<std::endl;
        }
#endif
      }
#endif

      //logfileptr->OFS()<<Xsuper_<<std::endl;

      Int nsuper = xsuper.size()-1;

      std::list<MPI_Request> mpirequests;


      Ptr nzbeg = 0;
      //nzend points to the last used slot in lindx
      Ptr nzend = 0;

      //tail is the end of list indicator (in rchlnk, not mrglnk)
      Idx tail = size +1;
      Idx head = 0;

      //Array of length nsuper containing the children of 
      //each supernode as a linked list
      std::vector<Idx> mrglnk(nsuper,0);

      //Array of length n+1 containing the current linked list 
      //of merged indices (the "reach" set)
      std::vector<Idx> rchlnk(size+1);

      //Array of length n used to mark indices as they are introduced
      // into each supernode's index set
      std::vector<Int> marker(size,0);

      Int nsuperLocal = XsuperDist_[iam+1]-XsuperDist_[iam];
      Int firstSnode = XsuperDist_[iam];
      Int lastSnode = firstSnode + nsuperLocal-1;
      xlindx.resize(nsuperLocal+1);


      Int firstColumn = Xsuper_[firstSnode-1];
      Int numColumnsLocal = Xsuper_[lastSnode] - firstColumn;


      //Compute the sum of the column count and resize lindx accordingly
      //nofsub will be the local nnz now
      Ptr nofsub = 0;
      for(Int ksup = 1; ksup<=nsuper; ++ksup){
        if(ksup>=firstSnode && ksup<=lastSnode){
          Int fstcol = xsuper[ksup-1];
          nofsub+=cc[fstcol-1];
        }
      }
      lindx.resize(nofsub);

      Ptr point = 1;
      for(Int ksup = 1; ksup<=nsuper; ++ksup){
        if(ksup>=firstSnode && ksup<=lastSnode){
          Int locSnode = ksup - firstSnode +1;
          Int fstcol = xsuper[ksup-1];
          xlindx[locSnode-1] = point;
          point += cc[fstcol-1]; 
        }
      } 
      xlindx[nsuperLocal] = point;
      //logfileptr->OFS()<<xlindx<<std::endl;

      std::vector<Ptr> recvXlindx;
      std::vector<Idx> recvLindx;

      if(iam>0){
        //build the mrglnk array


        for(Int ksup = 1; ksup<firstSnode; ++ksup){
          Int fstcol = xsuper[ksup-1];
          Int lstcol = xsuper[ksup]-1;
          Int width = lstcol - fstcol +1;
          Int length = cc[fstcol-1];

          //if ksup has a parent, insert ksup into its parent's 
          //"merge" list.
          if(length > width){
            Idx pcol = tree.PostParent(fstcol+width-1-1);  
            Int psup = SupMembership[pcol-1];
            mrglnk[ksup-1] = mrglnk[psup-1];
            mrglnk[psup-1] = ksup;
          }
        }

      }


      point = 1;
      for(Int locksup = 1; locksup<=nsuperLocal; ++locksup){
        Int ksup = locksup + firstSnode - 1;
        Int fstcol = xsuper[ksup-1];
        Int lstcol = xsuper[ksup]-1;
        Int width = lstcol - fstcol +1;
        Int length = cc[fstcol-1];
        Ptr knz = 0;
        rchlnk[head] = tail;

        //If ksup has children in the supernodal e-tree
        Int jsup = mrglnk[ksup-1];
        if(jsup>0){
          //copy the indices of the first child jsup into 
          //the linked list, and mark each with the value 
          //ksup.
          Int parentJ = ksup;
          do{
            Int jwidth = xsuper[jsup]-xsuper[jsup-1];

            Ptr * jxlindx = NULL;
            Idx * jlindx = NULL;
            Int locjsup = -1;
            if(jsup>=firstSnode && jsup<=lastSnode){
              locjsup = jsup - firstSnode +1;
              jxlindx = &xlindx[0];
              jlindx = &lindx[0];
            }
            else {
              MPI_Status status;
              recvLindx.resize(size);
              //receive jsup lindx
              Int psrc = 0; for(psrc = 0; psrc<iam;psrc++){ if(XsuperDist_[psrc]<=jsup && jsup<XsuperDist_[psrc+1]){ break; } }
              //logfileptr->OFS()<<"trying to recv "<<jsup<<" max "<<recvLindx.size()*sizeof(Idx)<<" bytes"<<" from P"<<psrc<<std::endl;
              MPI_Recv(&recvLindx[0],recvLindx.size()*sizeof(Idx),MPI_BYTE,psrc,jsup+all_np,comm,&status);
              //get actual number of received elements
              int count = 0;
              MPI_Get_count(&status,MPI_BYTE,&count);
              count/=sizeof(Idx);

              //compute jsup xlindx
              recvXlindx.resize(2);
              recvXlindx[0] = 1;
              recvXlindx[1] = count +1; 

              locjsup = 1;
              jxlindx = &recvXlindx[0];
              jlindx = &recvLindx[0];
            }

            Ptr jnzbeg = jxlindx[locjsup-1] + jwidth;
            Ptr jnzend = jxlindx[locjsup] -1;
            if(parentJ == ksup){
              for(Ptr jptr = jnzend; jptr>=jnzbeg; --jptr){
                Idx newi = jlindx[jptr-1];
                ++knz;
                marker[newi-1] = ksup;
                rchlnk[newi] = rchlnk[head];
                rchlnk[head] = newi;
              }
            }
            else{
              Int nexti = head;
              for(Ptr jptr = jnzbeg; jptr<=jnzend; ++jptr){
                Idx newi = jlindx[jptr-1];
                Idx i;
                do{
                  i = nexti;
                  nexti = rchlnk[i];
                }while(newi > nexti);

                if(newi < nexti){
                  ++knz;
                  rchlnk[i] = newi;
                  rchlnk[newi] = nexti;
                  marker[newi-1] = ksup;
                  nexti = newi;
                }
              }
            }

            parentJ = jsup;
            jsup = mrglnk[jsup-1];
          } while(jsup!=0 && knz < length);


          //TODO do better than this:need to avoid sending unnecessary data
          //receive the speculative sends
          jsup = mrglnk[ksup-1];
          //get the next element of the list
          jsup = mrglnk[jsup-1];
          while(jsup>0){
            Int psrc = 0; for(psrc = 0; psrc<iam;psrc++){ if(XsuperDist_[psrc]<=jsup && jsup<XsuperDist_[psrc+1]){ break; } }
            if(psrc!=iam){
              MPI_Status status;
              MPI_Request request;
              MPI_Irecv(&recvLindx[0],recvLindx.size()*sizeof(Idx),MPI_BYTE,psrc,jsup+all_np,comm,&request);
              MPI_Cancel(&request);
            }
            jsup = mrglnk[jsup-1];
          }
        }

        //structure of a(*,fstcol) has not been examined yet.  
        //"sort" its structure into the linked list,
        //inserting only those indices not already in the
        //list.
        if(knz < length){
          //loop on local columns instead for LOCAL EXPANDED structure
          for(Int row = fstcol; row<=lstcol; ++row){
            Idx newi = row;
            if(newi > fstcol && marker[newi-1] != ksup){
              //position and insert newi in list and
              // mark it with kcol
              Idx nexti = head;
              Idx i;
              do{
                i = nexti;
                nexti = rchlnk[i];
              }while(newi > nexti);
              ++knz;
              rchlnk[i] = newi;
              rchlnk[newi] = nexti;
              marker[newi-1] = ksup;
            }
          }


          for(Int col = fstcol; col<=lstcol; ++col){
            Int local_col = col - firstColumn + 1;
            Ptr knzbeg = pGraph.colptr[local_col-1];
            Ptr knzend = pGraph.colptr[local_col]-1;
            for(Ptr kptr = knzbeg; kptr<=knzend;++kptr){
              Idx newi = pGraph.rowind[kptr-1];

              if(newi > fstcol && marker[newi-1] != ksup){
                //position and insert newi in list and
                // mark it with kcol
                Idx nexti = head;
                Idx i;
                do{
                  i = nexti;
                  nexti = rchlnk[i];
                }while(newi > nexti);
                ++knz;
                rchlnk[i] = newi;
                rchlnk[newi] = nexti;
                marker[newi-1] = ksup;
              }
            }

            if(options_.relax.nrelax0==0){
              break;
            }
          }

        } 

        //if ksup has no children, insert fstcol into the linked list.
        if(rchlnk[head] != fstcol){
          rchlnk[fstcol] = rchlnk[head];
          rchlnk[head] = fstcol;
          ++knz;
        }

        {
          Idx i = head;
          for(Int col = fstcol; col<=lstcol; ++col){
            Int local_col = col - firstColumn + 1;
            Ptr knzbeg = pGraph.colptr[local_col-1];
            Ptr knzend = pGraph.colptr[local_col]-1;
            for(Ptr kptr = knzbeg; kptr<=knzend;++kptr){
              Idx newi = pGraph.rowind[kptr-1];
            }
            if(options_.relax.nrelax0==0){
              break;
            }
          }
        }

        bassert(knz == cc[fstcol-1]);


        //copy indices from linked list into lindx(*).
        nzbeg = nzend+1;
        nzend += knz;
        xlindx[locksup] = nzend+1;
        bassert(nzend+1 == xlindx[locksup]);
        Idx i = head;
        for(Ptr kptr = nzbeg; kptr<=nzend;++kptr){
          i = rchlnk[i];
          lindx[kptr-1] = i;
        } 

        //if ksup has a parent, insert ksup into its parent's 
        //"merge" list.
        if(length > width){
          Idx pcol = tree.PostParent(fstcol+width-1-1);  
          Int psup = SupMembership[pcol-1];
          mrglnk[ksup-1] = mrglnk[psup-1];
          mrglnk[psup-1] = ksup;

          //send L asap
          Int pdest = 0; for(pdest = 0; pdest<all_np;pdest++){ if(XsuperDist_[pdest]<=psup && psup<XsuperDist_[pdest+1]){ break; } }
          //if remote
          if(pdest!=iam){
            mpirequests.push_back(MPI_REQUEST_NULL);
            MPI_Request & request = mpirequests.back();
            MPI_Isend(&lindx[xlindx[locksup-1]-1],length*sizeof(Idx),MPI_BYTE,pdest,ksup+all_np,comm,&request);
          }
        }
      }

      bassert(nzend==0 || iam<np);
      lindx.resize(nzend);
      MPI_Barrier(comm);
    }


  template <typename T> 
    void symPACKMatrix<T>::refineSupernodes(int ordflag,int altflag,DistSparseMatrix<T> * pMat){
      ETree& tree = ETree_;
      Ordering & aOrder = Order_;
      std::vector<Int> & supMembership = SupMembership_; 
      std::vector<Int> & xsuper = Xsuper_; 

      std::vector<int> ixlindx;
      std::vector<int> ilindx;

      std::vector<int>  new_invp;

      //Gather locXlindx_ and locLindx_
      {
        std::vector<Ptr> xlindx;
        std::vector<Idx> lindx;
        //get other proc vertex counts
        Idx localVertexCnt = locXlindx_.size()-1;
        std::vector<Idx> remoteVertexCnt(np,0);
        MPI_Gather(&localVertexCnt,sizeof(localVertexCnt),MPI_BYTE,&remoteVertexCnt[0],sizeof(localVertexCnt),MPI_BYTE,0,CommEnv_->MPI_GetComm());
        Idx totalVertexCnt = std::accumulate(remoteVertexCnt.begin(),remoteVertexCnt.end(),0,std::plus<Idx>());
        if(iam==0){
          xlindx.resize(totalVertexCnt+1);
        }
        //compute receive displacements
        std::vector<int> rsizes(np,0);
        for(int p = 0; p<np;p++){rsizes[p] = (int)remoteVertexCnt[p]*sizeof(Ptr);}
        std::vector<int> rdispls(np+1,0);
        rdispls[0]=0;
        std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
        MPI_Gatherv(&locXlindx_[0],localVertexCnt*sizeof(Ptr),MPI_BYTE,&xlindx[0],&rsizes[0],&rdispls[0],MPI_BYTE,0,CommEnv_->MPI_GetComm());


        Ptr localEdgeCnt = locLindx_.size();
        std::vector<Ptr> remoteEdgeCnt(np,0);
        MPI_Gather(&localEdgeCnt,sizeof(localEdgeCnt),MPI_BYTE,&remoteEdgeCnt[0],sizeof(localEdgeCnt),MPI_BYTE,0,CommEnv_->MPI_GetComm());
        Ptr totalEdgeCnt = std::accumulate(remoteEdgeCnt.begin(),remoteEdgeCnt.end(),0,std::plus<Ptr>());

        //fix xlindx
        if(iam==0){

          Idx pos = remoteVertexCnt[0];
          Ptr offset = 0;
          for(int p=1;p<np;p++){
            offset+=remoteEdgeCnt[p-1]; 
            for(Idx lv = 0; lv < remoteVertexCnt[p]; lv++){
              xlindx[pos++] += offset;//remoteEdgeCnt[p-1];//(vertexDist[p] - baseval);
            }
          }
          xlindx.back()=totalEdgeCnt + 1;

          lindx.resize(totalEdgeCnt);
        }

        //compute receive displacements
        rsizes.assign(np,0);
        for(int p = 0; p<np;p++){rsizes[p] = (int)remoteEdgeCnt[p]*sizeof(Idx);}
        rdispls[0]=0;
        std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
        MPI_Gatherv(&locLindx_[0],localEdgeCnt*sizeof(Idx),MPI_BYTE,&lindx[0],&rsizes[0],&rdispls[0],MPI_BYTE,0,CommEnv_->MPI_GetComm());

        if(iam==0){
          ixlindx.resize(xlindx.size());
          for(int i = 0;i<xlindx.size();i++){
            ixlindx[i] = xlindx[i];
          }
          ilindx.resize(lindx.size());
          for(int i = 0;i<lindx.size();i++){
            ilindx[i] = lindx[i];
          }
        }

        //logfileptr->OFS()<<"locXlindx: "<<locXlindx_<<std::endl;
        //logfileptr->OFS()<<"locLindx: "<<locLindx_<<std::endl;
      }





      if(iam==0){

        //logfileptr->OFS()<<"xlindx: "<<ixlindx<<std::endl;
        //logfileptr->OFS()<<"lindx: "<<ilindx<<std::endl;

        int neqns = iSize_;
        //    int ordflag =2;
        //    int altflag =1;

        int nofsub =ilindx.size();
        int nsuper = xsuper.size()-1;

        std::vector<int>  freeforw(neqns,0);
        std::vector<int>  freeback(neqns,0);
        std::vector<int>  sforw(neqns,0);
        std::vector<int>  sback(neqns,0);
        std::vector<int> setseg_forw(neqns,0);
        std::vector<int>  setseg_back(neqns,0);
        std::vector<int>  nodehead(neqns,0);
        std::vector<int> nodeforw(neqns,0);
        std::vector<int>  nodeback(neqns,0);
        std::vector<int>  setsnode(neqns,0);
        std::vector<int>  supperm(nsuper,0);
        std::vector<int>  mark(neqns+1,0);
        std::vector<int>  set (neqns,0);
        std::vector<int>  compset(neqns,0);
        std::vector<int>  invp2(neqns,0);
        std::vector<int>  heap (2*nsuper,0);

        new_invp.assign(neqns,0);


        std::vector<int>  new_perm(neqns,0);
        std::iota(new_invp.begin(),new_invp.end(),1);
        std::iota(new_perm.begin(),new_perm.end(),1);

        FORTRAN(ordsup) (
            &ordflag, &altflag, &neqns , &nofsub, &nsuper, 
            &xsuper[0], &ixlindx[0], &ilindx[0], &supMembership[0], 
            &new_perm[0], &new_invp[0], 
            //        &aOrder.perm[0], &aOrder.invp[0], 
            &freeforw[0], &freeback[0], &sforw[0], &sback[0], 
            &setseg_forw[0], &setseg_back[0], &nodehead[0], 
            &nodeforw[0], &nodeback[0], 
            &setsnode[0], &supperm[0], &mark[0], &set[0], &compset[0],
            &invp2[0], &heap[0]);

        Order_.Compose(new_invp);

      }
      

#ifdef EXPLICIT_PERMUTE
      
      //Bcast the individual permutation
      new_invp.resize(iSize_);
      MPI_Bcast(&new_invp[0],iSize_*sizeof(int),MPI_BYTE,0,fullcomm_);

      //re permute the matrix
      pMat->Permute(&new_invp[0]);
      graph_ = pMat->GetLocalGraph();
      graph_.SetBaseval(1);
      graph_.SetKeepDiag(1);
      graph_.SetSorted(1);
#endif
      // broadcast invp
      Int N = aOrder.invp.size();
      MPI_Bcast(&aOrder.invp[0],N*sizeof(int),MPI_BYTE,0,fullcomm_);
      MPI_Bcast(&aOrder.perm[0],N*sizeof(int),MPI_BYTE,0,fullcomm_);
    }


  template <typename T> 
    template <class Allocator>
    SuperNode<T,Allocator> * symPACKMatrix<T>::CreateSuperNode(DecompositionType type,Int aiId, Int aiFc, Int aiLc, Int ai_num_rows, Int aiN, Int aiNZBlkCnt){
      SuperNode<T,Allocator> * retval = NULL;
      switch(type){
        case LDL:
          retval = new SuperNodeInd<T,Allocator>( aiId,  aiFc,  aiLc,  ai_num_rows,  aiN,  aiNZBlkCnt);
          break;
        case LL:
          retval = new SuperNode<T,Allocator>( aiId,  aiFc,  aiLc,  ai_num_rows,  aiN,  aiNZBlkCnt);
          break;
        default:
          retval = new SuperNode<T,Allocator>( aiId,  aiFc,  aiLc,  ai_num_rows,  aiN,  aiNZBlkCnt);
          break;
      }
      return retval;
    }

  template <typename T> 
    template <class Allocator>
    SuperNode<T,Allocator> * symPACKMatrix<T>::CreateSuperNode(DecompositionType type,Int aiId, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices){
      SuperNode<T,Allocator> * retval = NULL;
      switch(type){
        case LDL:
          retval = new SuperNodeInd<T,Allocator>( aiId,  aiFc,  aiLc,  aiN,  rowIndices);
          break;
        case LL:
          retval = new SuperNode<T,Allocator>( aiId,  aiFc,  aiLc,  aiN,  rowIndices);
          break;
        default:
          retval = new SuperNode<T,Allocator>( aiId,  aiFc,  aiLc,  aiN,  rowIndices);
          break;
      }
      return retval;
    }



  template <typename T>
    template <class Allocator>
    SuperNode<T,Allocator> * symPACKMatrix<T>::CreateSuperNode(DecompositionType type,char * dataPtr,size_t size, Int firstRow){
      SuperNode<T,Allocator> * retval = NULL;
      switch(type){
        case LDL:
          retval = new SuperNodeInd<T,Allocator>(dataPtr,size,firstRow);
          break;
        case LL:
          retval = new SuperNode<T,Allocator>(dataPtr,size,firstRow);
          break;
        default:
          retval = new SuperNode<T,Allocator>(dataPtr,size,firstRow);
          break;
      }
      return retval;
    }





#include "symPACKMatrix_impl_FB_pull.hpp"


#ifdef _INDEFINITE_
#undef _INDEFINITE_
#endif

}


#endif 
