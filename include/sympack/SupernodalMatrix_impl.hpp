#ifndef _SUPERNODAL_MATRIX_IMPL_HP_
#define _SUPERNODAL_MATRIX_IMPL_HPP_

#include "sympack/SupernodalMatrix.hpp"

#include "sympack/blas.hpp"
#include "sympack/Ordering.hpp"
#include "sympack/LoadBalancer.hpp"
#include "sympack/mpi_interf.hpp"

#include  "sympack/DistSparseMatrixGraph.hpp"

#include <queue>
#include <stdlib.h>
#include <numeric>
#include <algorithm>

//#define _INDEFINITE_

namespace SYMPACK{
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




namespace SYMPACK{

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
    TIMER_START(FB_FIND_TASK);
    //find task corresponding to curUpdate
    auto taskit = taskLists_[tgt-1]->begin();
    for(;
        taskit!=taskLists_[tgt-1]->end();
        taskit++){
      if(taskit->src_snode_id==src && taskit->tgt_snode_id==tgt && taskit->type==type){
        break;
      }
    }
    TIMER_STOP(FB_FIND_TASK);
    return taskit;
  }


  template<typename T> void SupernodalMatrix<T>::generateTaskGraph(supernodalTaskGraph & taskGraph,
      SYMPACK::vector<Int> & AggregatesToRecv,  SYMPACK::vector<Int>& LocalAggregates)
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

    //then do an alltoallv
    //compute send sizes
    vector<int> ssizes(np,0);
    for(auto itp = Updates.begin();itp!=Updates.end();itp++){
      ssizes[itp->first] = itp->second.size()*sizeof(std::pair<Idx,Idx>);
    }

    //compute send displacements
    vector<int> sdispls(np+1,0);
    sdispls[0] = 0;
    std::partial_sum(ssizes.begin(),ssizes.end(),&sdispls[1]);

    //Build the contiguous array of pairs
    vector<std::pair<Idx,Idx> > sendbuf;
    sendbuf.reserve(sdispls.back()/sizeof(std::pair<Idx,Idx>));

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
    vector<std::pair<Idx,Idx> > recvbuf(rdispls.back()/sizeof(std::pair<Idx,Idx>));
    MPI_Alltoallv(&sendbuf[0],&ssizes[0],&sdispls[0],MPI_BYTE,&recvbuf[0],&rsizes[0],&rdispls[0],MPI_BYTE,CommEnv_->MPI_GetComm());

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


    SYMPACK::vector<Int> levels;
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





















  template <typename T> void SupernodalMatrix<T>::Factorize(){
    TIMER_START(FACTORIZATION);
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
    TIMER_STOP(FACTORIZATION);

#ifdef PROFILE_COMM
    logfileptr->OFS()<<"Local volume of communication: "<<gVolComm<<endl;
    logfileptr->OFS()<<"Local number of messages: "<<gNumMsg<<endl;

    size_t totalVolComm = 0;
    team_->reduce(&gVolComm,&totalVolComm,1,0,UPCXX_SUM);
    size_t totalNumMsg = 0;
    team_->reduce(&gNumMsg,&totalNumMsg,1,0,UPCXX_SUM);

    if(iam==0){
      cout<<"Total volume of communication: "<<totalVolComm<<endl;
      cout<<"Total number of messages: "<<totalNumMsg<<endl;
    }

    gVolComm=0;
    gNumMsg=0;
#endif
  }


  //Solve related routines
  template <typename T> void SupernodalMatrix<T>::Solve(T * RHS, int nrhs,  T * Xptr) {
    TIMER_START(SPARSE_SOLVE);

    Int n = iSize_;
    Int iam = CommEnv_->MPI_Rank();
    Int np  = CommEnv_->MPI_Size();

    SYMPACK::vector<Int> children(Xsuper_.size());
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
    logfileptr->OFS()<<"Children SYMPACK::vector is"<<children<<std::endl;
#endif


    SYMPACK::vector<Int> UpdatesToDo = children;

    Contributions_.resize(LocalSupernodes_.size());
    SYMPACK::vector<std::stack<Int> > LocalUpdates(LocalSupernodes_.size());

    AsyncComms outgoingSend;

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
//#ifdef _INDEFINITE_
//        SuperNode<T,MallocAllocator> * contrib = new SuperNodeInd<T,MallocAllocator>(I,1,nrhs, cur_snode->NRowsBelowBlock(0) ,iSize_);
//#else
//        SuperNode<T,MallocAllocator> * contrib = new SuperNode<T,MallocAllocator>(I,1,nrhs, cur_snode->NRowsBelowBlock(0) ,iSize_);
//#endif
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


    SYMPACK::vector<char> src_blocks;


    //forward-substitution phase
    //Sending contrib up the tree
    //Start from the leaves of the tree
    TIMER_START(SPARSE_FWD_SUBST);

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
          contrib->forward_update(dist_contrib,iOwner);
          //this->forward_update(dist_contrib,contrib);
          //delete contributions[(contrib_snode_id-1) / np];
          --UpdatesToDo[I-1];
        }

        //do remote updates
        size_t max_bytes;
        Int nz_cnt;
        while(UpdatesToDo[I-1]>0){
          //receive children contrib
#ifdef _DEBUG_
          logfileptr->OFS()<<UpdatesToDo[I-1]<<" contribs left"<<endl;
#endif


          MPI_Status recv_status;
          int bytes_received = 0;

          TIMER_START(RECV_MPI);
          MPI_Probe(MPI_ANY_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
          src_blocks.resize(bytes_received);
          MPI_Recv(&src_blocks[0],bytes_received,MPI_BYTE,recv_status.MPI_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
          TIMER_STOP(RECV_MPI);
//#ifdef _INDEFINITE_
//          SuperNode<T,MallocAllocator> * dist_contrib = new SuperNodeInd<T,MallocAllocator>(&src_blocks[0],bytes_received);
//#else
//          SuperNode<T,MallocAllocator> * dist_contrib = new SuperNode<T,MallocAllocator>(&src_blocks[0],bytes_received);
//#endif
          SuperNode<T,MallocAllocator> * dist_contrib = CreateSuperNode<MallocAllocator>(options_.decomposition,&src_blocks[0],bytes_received);
          //TODO Put this back
          //Deserialize(&src_blocks[0],dist_contrib);
#ifdef _DEBUG_
          logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib.Id()<<std::endl;
#endif

          Int iOwner = recv_status.MPI_SOURCE;
          contrib->forward_update(dist_contrib,iOwner);
          //this->forward_update(dist_contrib,contrib);

          --UpdatesToDo[I-1];
          delete dist_contrib;
        }

        assert(UpdatesToDo[I-1]==0);

        if(UpdatesToDo[I-1]==0){

          //now compute MY contribution
          contrib->forward_update_contrib(RHS,cur_snode,this->Order_.perm);

#if 0
          //This corresponds to the i loop in dtrsm
          for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){

            NZBlockDesc & cur_desc = contrib->GetNZBlockDesc(blkidx);
            NZBlockDesc & chol_desc = cur_snode->GetNZBlockDesc(blkidx);
            NZBlockDesc & diag_desc = contrib->GetNZBlockDesc(0);

            Int cur_nrows = contrib->NRows(blkidx);
            Int chol_nrows = cur_snode->NRows(blkidx);
            Int diag_nrows = contrib->NRows(0);

            T * cur_nzval = contrib->GetNZval(cur_desc.Offset);
            T * chol_nzval = cur_snode->GetNZval(chol_desc.Offset);
            T * diag_nzval = contrib->GetNZval(diag_desc.Offset);

            //compute my contribution
            //Handle the diagonal block
            if(blkidx==0){
              //TODO That's where we can use the selective inversion
              //if we are processing the "pivot" block
              for(Int kk = 0; kk<cur_snode->Size(); ++kk){
                for(Int j = 0; j<nrhs;++j){
                  Int srcRow = Order_.perm[diag_desc.GIndex+kk-1];
                  //diag_nzval[kk*nrhs+j] = (B(srcRow-1,j) + diag_nzval[kk*nrhs+j]) / chol_nzval[kk*cur_snode->Size()+kk];
                  diag_nzval[kk*nrhs+j] = (RHS[srcRow-1 + j*n] + diag_nzval[kk*nrhs+j]) / chol_nzval[kk*cur_snode->Size()+kk];
                  for(Int i = kk+1; i<cur_nrows;++i){
                    diag_nzval[i*nrhs+j] += -diag_nzval[kk*nrhs+j]*chol_nzval[i*cur_snode->Size()+kk];
                  }
                }
              }
            }
            else{
              for(Int kk = 0; kk<cur_snode->Size(); ++kk){
                for(Int j = 0; j<nrhs;++j){
                  for(Int i = 0; i<cur_nrows;++i){
                    cur_nzval[i*nrhs+j] += -diag_nzval[kk*nrhs+j]*chol_nzval[i*cur_snode->Size()+kk];
                  }
                }
              }
            }
          }
#endif

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
              if(next_local_contrib< parent_snode_id){
                //need to push the prev src_last_row
                ContribsToSend.push(DelayedComm(contrib->Id(),parent_snode_id,1,src_first_row));
#ifdef _DEBUG_DELAY_
                cout<<"P"<<iam<<" has delayed update from Contrib "<<I<<" to "<<parent_snode_id<<" from row "<<src_first_row<<endl;
#endif
                isSkipped= true;
              }

              if(!isSkipped){
                //Create a new Icomm buffer, serialize the contribution
                // in it and add it to the outgoing comm list

                Icomm * send_buffer = new Icomm();
                contrib->Serialize(*send_buffer,src_nzblk_idx,src_first_row);
                AddOutgoingComm(outgoingSend,send_buffer);

                if( outgoingSend.size() > maxIsend_){
                  TIMER_START(SEND_MPI);
                  MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,parent_snode_id,CommEnv_->MPI_GetComm());
                  TIMER_STOP(SEND_MPI);
                  outgoingSend.pop_back();
                }
                else{
                  TIMER_START(SEND_MPI);
                  MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,parent_snode_id,CommEnv_->MPI_GetComm(),&outgoingSend.back()->Request);
                  TIMER_STOP(SEND_MPI);
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

#if 0
#ifdef _CHECK_RESULT_SEQ_
      //if(iLocalI>0 && iLocalI<=LocalSupernodes_.size())
      {
        MPI_Barrier(CommEnv_->MPI_GetComm());
        NumMat<T> tmp = B;
        GetSolution(tmp);


        Int nrows = 0;
        for(Int ii=1; ii<=I;++ii){ nrows+= Xsuper_[ii] - Xsuper_[ii-1];}

        NumMat<T> tmp3(tmp.size(),tmp.n());
        NumMat<T> tmpFwd(tmp.size(),tmp.n());
        for(Int i = 0; i<tmp.size();++i){
          for(Int j = 0; j<tmp.n();++j){
            //          tmp3(Order_.perm[i]-1,j) = tmp(i,j);
            //          tmpFwd(Order_.perm[i]-1,j) = forwardSol(i,j);

            tmp3(i,j) = tmp(Order_.perm[i]-1,j);
            tmpFwd(i,j) = forwardSol(Order_.perm[i]-1,j);
          }
        }

        NumMat<T> tmp2 = tmp3;

        blas::Axpy(tmp.size()*tmp.n(),-1.0,&tmpFwd(0,0),1,&tmp2(0,0),1);
        double norm = lapack::Lange('F',nrows,tmp.n(),&tmp2(0,0),tmp.size());
        logfileptr->OFS()<<"Norm after SuperNode2 "<<I<<" is "<<norm<<std::endl; 

        if(abs(norm)>=1e-1){
          for(Int i = 0;i<nrows/*tmp.size()*/;++i){
            logfileptr->OFS()<<tmpFwd(i,0)<<"       "<<tmp3(i,0)<<std::endl;
          }
        }
      }
#endif
#endif
      ++iLocalI;
    }

    while(!outgoingSend.empty()){
      AdvanceOutgoing(outgoingSend);
    } 
    MPI_Barrier(CommEnv_->MPI_GetComm());

    TIMER_STOP(SPARSE_FWD_SUBST);

    //DumpContrib();

    //Back-substitution phase
    TIMER_START(SPARSE_BACK_SUBST);

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

        SYMPACK::vector<Int> isBlockUpdated(contrib->NZBlockCnt(),0);

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

//#ifdef _INDEFINITE_
//            dist_contrib = new SuperNodeInd<T,MallocAllocator>(&src_blocks[0],bytes_received);
//#else
//            dist_contrib = new SuperNode<T,MallocAllocator>(&src_blocks[0],bytes_received);
//#endif
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
#if 0
#ifdef _DEBUG_
        logfileptr->OFS()<<"BACK Processing contrib "<<I<<std::endl;
#endif

        NZBlockDesc & diag_desc = cur_snode->GetNZBlockDesc(0);
        NZBlockDesc & tgt_desc = contrib->GetNZBlockDesc(0);

        T* diag_nzval = cur_snode->GetNZval(diag_desc.Offset);
        T* tgt_nzval = contrib->GetNZval(tgt_desc.Offset);

        for(Int j = 0; j<nrhs;++j){
          for(Int ii = cur_snode->Size()-1; ii>=0; --ii){
            T temp = tgt_nzval[ii*nrhs+j];

            //This corresponds to the k loop in dtrsm
            for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){
              NZBlockDesc & chol_desc = cur_snode->GetNZBlockDesc(blkidx);
              Int chol_nrows = cur_snode->NRows(blkidx);

              Int src_blkidx = contrib->FindBlockIdx(chol_desc.GIndex);
              NZBlockDesc & cur_desc = contrib->GetNZBlockDesc(src_blkidx);
              Int cur_nrows = contrib->NRows(src_blkidx);

              T* chol_nzval = cur_snode->GetNZval(chol_desc.Offset);
              T* cur_nzval = contrib->GetNZval(cur_desc.Offset);

              for(Int kk = 0; kk< chol_nrows; ++kk){
                if(chol_desc.GIndex+kk>cur_snode->FirstCol()+ii){
                  Int src_row = chol_desc.GIndex - cur_desc.GIndex +kk;
                  if(src_row< cur_nrows){
                    temp += -chol_nzval[kk*cur_snode->Size()+ii]*cur_nzval[src_row*nrhs+j];
                  }
                }
              }
            }

            temp = temp / diag_nzval[ii*cur_snode->Size()+ii];
            tgt_nzval[ii*nrhs+j] = temp;
          }
        }
#endif

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
                  if(next_local_contrib > child_snode_id){
                    //need to push the prev src_last_row
                    //Send
                    Int src_nzblk_idx = 0;
                    NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);
                    Int src_first_row = pivot_desc.GIndex;
                    ContribsToSendDown.push(DelayedComm(contrib->Id(),child_snode_id,0,src_first_row));
#ifdef _DEBUG_DELAY_
                    cout<<"P"<<iam<<" has delayed update from Contrib "<<I<<" to "<<child_snode_id<<" from row "<<src_first_row<<endl;
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
                      TIMER_START(SEND_MPI);
                      MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm());
                      TIMER_STOP(SEND_MPI);
                      outgoingSend.pop_back();
                    }
                    else{
                      TIMER_START(SEND_MPI);
                      MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm(),&outgoingSend.back()->Request);
                      TIMER_STOP(SEND_MPI);
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
    TIMER_STOP(SPARSE_BACK_SUBST);

    MPI_Barrier(CommEnv_->MPI_GetComm());

    //DumpContrib();
    TIMER_STOP(SPARSE_SOLVE);

  }

  template<typename T> void SupernodalMatrix<T>::GetSolution(T * B, int nrhs){
    Int n = iSize_;
    Int iam = CommEnv_->MPI_Rank();
    Int np  = CommEnv_->MPI_Size();
    //    Int nrhs = B.n();
    //Gather B from everybody and put it in the original matrix order
    SYMPACK::vector<T> tmp_nzval;
    for(Int I=1; I<Xsuper_.size();++I){
      Int iOwner = this->Mapping_->Map(I-1,I-1);
      T * data;
      Int snode_size = Xsuper_[I] - Xsuper_[I-1];
      Int nzcnt = snode_size * nrhs;
      tmp_nzval.resize(nzcnt);

      if( iOwner == iam ){
        SuperNode<T,MallocAllocator> * contrib = snodeLocal(I,Contributions_);
        data = contrib->GetNZval(0);
      }
      else{
        data = &tmp_nzval[0];
      }

      MPI_Bcast(data,nzcnt*sizeof(T),MPI_BYTE,iOwner,CommEnv_->MPI_GetComm());

      for(Int i = 0; i<snode_size; ++i){ 
        for(Int j = 0; j<nrhs; ++j){
          Int destRow = Xsuper_[I-1] + i;
          destRow = Order_.perm[destRow - 1];
          //B(destRow-1, j) = data[i*nrhs + j];
          B[destRow-1 + j*n] = data[i*nrhs + j];
        }
      }
    }
    MPI_Barrier(CommEnv_->MPI_GetComm());
  }


  template <typename T>
    template< class Alloc>
    void SupernodalMatrix<T>::SendDelayedMessagesUp(Int iLocalI, CommList & MsgToSend, AsyncComms & OutgoingSend, SYMPACK::vector<SuperNode<T,Alloc> *> & snodeColl){
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
        logfileptr->OFS()<<"Picked { "<<src_snode_id<<" -> "<<tgt_snode_id<<" }"<<endl;
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
    void SupernodalMatrix<T>::SendDelayedMessagesDown(Int iLocalI, DownCommList & MsgToSend, AsyncComms & OutgoingSend, SYMPACK::vector<SuperNode<T,Alloc> *> & snodeColl){
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
        logfileptr->OFS()<<"Picked { "<<src_snode_id<<" -> "<<tgt_snode_id<<" }"<<endl;
#endif

        if(tgt_snode_id>next_snode_id || is_last /*|| OutgoingSend.size() <= maxIsend_*/){

          SuperNode<T,Alloc> * prev_src_snode = snodeLocal(src_snode_id,snodeColl);
          //this can be sent now

          Int iTarget = this->Mapping_->Map(tgt_snode_id-1,tgt_snode_id-1);
          if(iTarget != iam){
#ifdef _DEBUG_DELAY_
            logfileptr->OFS()<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode->Id()<<" to Supernode "<<tgt_snode_id<<endl;
            cout<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode->Id()<<" to Supernode "<<tgt_snode_id<<endl;
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
              TIMER_START(SEND_MPI);
              MPI_Send(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm());
              TIMER_STOP(SEND_MPI);
              OutgoingSend.pop_back();
            }
            else{
              TIMER_START(SEND_MPI);
              MPI_Isend(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm(),&OutgoingSend.back()->Request);
              TIMER_STOP(SEND_MPI);
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
    void SupernodalMatrix<T>::SendMessage(const DelayedComm & comm, AsyncComms & OutgoingSend, SYMPACK::vector<SuperNode<T,Alloc> *> & snodeColl){
      Int src_snode_id = comm.src_snode_id;
      Int tgt_snode_id = comm.tgt_snode_id;
      Int src_nzblk_idx = comm.src_nzblk_idx;
      Int src_first_row = comm.src_first_row;


      SuperNode<T,Alloc> * prev_src_snode = snodeLocal(src_snode_id,snodeColl);

      //this can be sent now
      Int iTarget = this->Mapping_->Map(tgt_snode_id-1,tgt_snode_id-1);
      if(iTarget != iam){
#ifdef _DEBUG_DELAY_
        logfileptr->OFS()<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode->Id()<<" to Supernode "<<tgt_snode_id<<endl;
        cout<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode->Id()<<" to Supernode "<<tgt_snode_id<<endl;
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
          TIMER_START(SEND_MPI);
          MPI_Send(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm());
          TIMER_STOP(SEND_MPI);
          OutgoingSend.pop_back();
        }
        else{
          TIMER_START(SEND_MPI);
          MPI_Isend(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm(),&OutgoingSend.back()->Request);
          TIMER_STOP(SEND_MPI);
        }

#ifdef _STAT_COMM_
        maxFactors_ = max(maxFactors_,send_buffer->size());
        sizesFactors_ += send_buffer->size();
        countFactors_++;
#endif


      }
    }




  template <typename T> void SupernodalMatrix<T>::AdvanceOutgoing(AsyncComms & outgoingSend){
    TIMER_START(ADVANCE_OUTGOING_COMM);
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
    TIMER_STOP(ADVANCE_OUTGOING_COMM);
  }


  template<typename T> void SupernodalMatrix<T>::AddOutgoingComm(AsyncComms & outgoingSend, Icomm * send_buffer){
    outgoingSend.push_back(send_buffer);
  }


  template <typename T> void SupernodalMatrix<T>::GetUpdatingSupernodeCount(SYMPACK::vector<Int> & sc,SYMPACK::vector<Int> & mw, SYMPACK::vector<Int> & mh, SYMPACK::vector<Int> & numBlk){
    sc.resize(Xsuper_.size(),I_ZERO);
    SYMPACK::vector<Int> marker(Xsuper_.size(),I_ZERO);
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
            mw[supno-1] = max(mw[supno-1],last_col - first_col+1);
          }
        }
        prevSnode = supno;
      }
      numBlk[s-1] = nzBlockCnt;
    }

    //do an allreduce on sc, mw and mh
    MPI_Allreduce(MPI_IN_PLACE,&sc[0],sc.size(),MPI_INT,MPI_SUM,CommEnv_->MPI_GetComm());
    MPI_Allreduce(MPI_IN_PLACE,&mw[0],mw.size(),MPI_INT,MPI_MAX,CommEnv_->MPI_GetComm());
    MPI_Allreduce(MPI_IN_PLACE,&mh[0],mh.size(),MPI_INT,MPI_MAX,CommEnv_->MPI_GetComm());
    MPI_Allreduce(MPI_IN_PLACE,&numBlk[0],numBlk.size(),MPI_INT,MPI_SUM,CommEnv_->MPI_GetComm());
  }


  template <typename T> SparseMatrixStructure SupernodalMatrix<T>::GetLocalStructure() const {
    return *Local_;
  }

  template <typename T> SparseMatrixStructure SupernodalMatrix<T>::GetGlobalStructure(){
    if(isGlobStructAllocated_){
      Global_=new SparseMatrixStructure();
      Local_->ToGlobal(*Global_,CommEnv_->MPI_GetComm());
      isGlobStructAllocated_= true;
    }
    return *Global_;
  }


  template<typename T>
    void SupernodalMatrix<T>::Dump(){
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
    void SupernodalMatrix<T>::DumpContrib(){
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



  template <typename T> void SupernodalMatrix<T>::Init(DistSparseMatrix<T> & pMat, NGCholOptions & options ){


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

    TIMER_START(SUPERMATRIX_INIT);
    //Create the CommEnvironment object if necessary
    if(CommEnv_!=NULL){
      delete CommEnv_;
    }

    if(options.commEnv == NULL){
      throw std::runtime_error("The communication environment must be initialized in the options");
    }

    options_ = options;

    CommEnv_ = new CommEnvironment(*options.commEnv);

    TIMER_START(MAPPING);

    //create the mapping
    Int pmapping = options.used_procs(np);
    Int pr = (Int)sqrt((double)pmapping);
    if(options.mappingTypeStr ==  "MODWRAP2DTREE"){
      this->Mapping_ = new ModWrap2DTreeMapping(pmapping, pr, pr, 1);
    }
    else if(options.mappingTypeStr ==  "WRAP2D"){
      this->Mapping_ = new Wrap2D(pmapping, pr, pr, 1);
    }
    else if(options.mappingTypeStr ==  "WRAP2DFORCED"){
      this->Mapping_ = new Wrap2DForced(pmapping, pr, pr, 1);
    }
    else if(options.mappingTypeStr ==  "MODWRAP2DNS"){
      this->Mapping_ = new Modwrap2DNS(pmapping, pr, pr, 1);
    }
    else if(options.mappingTypeStr ==  "ROW2D"){
      this->Mapping_ = new Row2D(pmapping, pmapping, pmapping, 1);
    }
    else if(options.mappingTypeStr ==  "COL2D"){
      this->Mapping_ = new Col2D(pmapping, pmapping, pmapping, 1);
    }
    else{
      this->Mapping_ = new Modwrap2D(pmapping, pr, pr, 1);
    }

    TIMER_STOP(MAPPING);

    //Options
    maxIsend_ = options.maxIsend;
    maxIrecv_ = options.maxIrecv;

    gMaxIrecv = maxIrecv_+maxIsend_;

    TIMER_START(SCHEDULER);
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
    TIMER_STOP(SCHEDULER);

    iSize_ = pMat.size;
    Local_ = new SparseMatrixStructure();
    *Local_ = pMat.GetLocalStructure();

    graph_.SetComm(CommEnv_->MPI_GetComm());
    graph_.SetBaseval(1);
    graph_.SetKeepDiag(1);
    graph_.SetSorted(1);
    graph_.FromStructure(*Local_);
    graph_.ExpandSymmetric();
    logfileptr->OFS()<<"Matrix structure expanded"<<endl;

    SparseMatrixGraph * sgraph = NULL;
    {

#if 0
      Global_ = new SparseMatrixStructure();
      Local_->ToGlobal(*Global_,CommEnv_->MPI_GetComm());
      Global_->ExpandSymmetric();
      isGlobStructAllocated_ = true;
#endif

      {
        //Create an Ordering object to hold the permutation
        Order_.SetCommEnvironment(CommEnv_);
#if 0
        Order_.SetStructure(*Global_);
        logfileptr->OFS()<<"Structure set"<<endl;
#endif

        double timeSta = get_time();
        TIMER_START(ORDERING);
        switch(options_.ordering){
          case MMD:
            {
              sgraph = new SparseMatrixGraph();
              graph_.GatherStructure(*sgraph,0);
              sgraph->SetBaseval(1);
              sgraph->SetKeepDiag(0);
              //Reoder the matrix with MMD
              Order_.MMD(*sgraph);
            }
              break;
          case AMD:
            {
              sgraph = new SparseMatrixGraph();
              graph_.GatherStructure(*sgraph,0);
              sgraph->SetBaseval(1);
              sgraph->SetKeepDiag(0);
            Order_.AMD(*sgraph);
            }
            break;
            
          case NDBOX:
            Order_.NDBOX(Size());
            break;
          case NDGRID:
            Order_.NDGRID(Size());
            break;
#ifdef USE_SCOTCH
          case SCOTCH:
            {

              sgraph = new SparseMatrixGraph();
              graph_.GatherStructure(*sgraph,0);
              sgraph->SetKeepDiag(0);
              Order_.SCOTCH(*sgraph);
            }
            break;
#endif
#ifdef USE_METIS
          case METIS:
            {
              sgraph = new SparseMatrixGraph();
              graph_.GatherStructure(*sgraph,0);
              sgraph->SetKeepDiag(0);
              Order_.METIS(*sgraph);
            }
            break;
#endif
#ifdef USE_PARMETIS
          case PARMETIS:
            {
              DistSparseMatrixGraph graph;
              graph.SetComm(CommEnv_->MPI_GetComm());
              graph.SetBaseval(0);
              graph.SetKeepDiag(0);
              graph.SetSorted(1);
              graph.FromStructure(*Local_);
              graph.ExpandSymmetric();
              Order_.PARMETIS(graph);
            }
            break;
#endif
#ifdef USE_PTSCOTCH
          case PTSCOTCH:
            {
              DistSparseMatrixGraph graph;
              graph.SetComm(CommEnv_->MPI_GetComm());
              graph.SetBaseval(0);
              graph.SetKeepDiag(0);
              graph.SetSorted(1);
              graph.FromStructure(*Local_);
              graph.ExpandSymmetric();
              Order_.PTSCOTCH(graph);
            }
            break;
#endif
          default:
            //do nothing: either natural or user provided ordering
            break;
        }
        TIMER_STOP(ORDERING);

        double timeStop = get_time();
        if(iam==0){
          cout<<"Ordering time: "<<timeStop - timeSta<<endl;
        }
        logfileptr->OFS()<<"Ordering done"<<endl;
      }
    }

#if 0
    if(Global_==NULL){
      //TODO THIS WILL HAVE TO BE REMOVED AT SOME POINT
      Global_ = new SparseMatrixStructure();
      Local_->ToGlobal(*Global_,CommEnv_->MPI_GetComm());
      Global_->ExpandSymmetric();
      isGlobStructAllocated_ = true;
      logfileptr->OFS()<<"Matrix structure expanded"<<endl;
    }
#endif

    SYMPACK::vector<Int> cc,rc;
    {

      double timeSta = get_time();
      //ETree_.ConstructETree(graph_,Order_);
      {
        if(sgraph==NULL){ 
          sgraph = new SparseMatrixGraph();
          graph_.GatherStructure(*sgraph,0);
        }
        sgraph->SetKeepDiag(1);
        sgraph->SetBaseval(1);
        ETree_.ConstructETree(*sgraph,Order_,team_);
      }

#if 0
      Ordering newOrder;
      newOrder.SetCommEnvironment(CommEnv_);
      newOrder.SetStructure(*Global_);
      newOrder.perm = Order_.perm;
      newOrder.invp = Order_.invp;
#endif


      ETree_.PostOrderTree(Order_);
      logfileptr->OFS()<<"ETREE done"<<endl;
#if 0
      {
        double timeStart = get_time();
        this->getLColRowCount(cc,rc);
        double timeStop = get_time();
        if(iam==0){
          cout<<"Column count (upcxx) construction time: "<<timeStop - timeStart<<endl;
        }
      }
#endif
      {
        double timeStart = get_time();
        if(sgraph==NULL){ 
          sgraph = new SparseMatrixGraph();
          graph_.GatherStructure(*sgraph,0);
        }
        sgraph->SetKeepDiag(1);
        sgraph->SetBaseval(1);
        //Global_->GetLColRowCount(tree,*sgraph,cc,rc);
        this->getLColRowCount(*sgraph,cc,rc);
        //bcast rc, cc
        double timeStop = get_time();
        if(iam==0){
          cout<<"Column count (gather + serial + bcast) construction time: "<<timeStop - timeStart<<endl;
        }
      }
#if 0
      { 
        SYMPACK::vector<Int> cc2,rc2;
        ETree tree;
        tree.ConstructETree(*Global_,newOrder);
        tree.PostOrderTree(newOrder);
        double timeStart = get_time();
        Global_->GetLColRowCount(tree,newOrder,cc2,rc2);
        double timeStop = get_time();
        if(iam==0){
          cout<<"Column count (serial) construction time: "<<timeStop - timeStart<<endl;
        }
        for(int i = 0; i<cc.size(); i++){ bassert(cc[i]==cc2[i]);}
        for(int i = 0; i<rc.size(); i++){ bassert(rc[i]==rc2[i]);}
      }
#endif

      ETree_.SortChildren(cc,Order_);
      double timeStop = get_time();
      if(iam==0){
        cout<<"Elimination tree construction time: "<<timeStop - timeSta<<endl;
      }
      //logfileptr->OFS()<<"perm: "<<Order_.perm<<endl;
    }

    if(sgraph!=NULL){delete sgraph;}


#ifdef _DEBUG_
    logfileptr->OFS()<<"colcnt "<<cc<<std::endl;
    logfileptr->OFS()<<"rowcnt "<<rc<<std::endl;
#endif
    double flops = 0.0;
    for(Int i = 0; i<cc.size();++i){
      flops+= (double)pow((double)cc[i],2.0);
    }

    if(iam==0){
      cout<<"Flops: "<<flops<<endl;
    }

#if 1
    int64_t NNZ = 0;
    for(Int i = 0; i<cc.size();++i){
      NNZ+=cc[i];
    }
    if(iam==0){
      cout<<"NNZ in L factor: "<<NNZ<<endl;
    }
#endif

#ifndef _NO_COMPUTATION_
    {
      double timeSta = get_time();
      this->findSupernodes(ETree_,Order_,cc,SupMembership_,Xsuper_,options.maxSnode);
      //logfileptr->OFS()<<"xsuper: "<<Xsuper_<<endl;

      //Compute XsuperDist_
      {
        Idx supPerProc = (Xsuper_.size()-1) / np;
        XsuperDist_.resize(np+1,0);
        for(int p =0; p<np; p++){
          XsuperDist_[p]= p*supPerProc+1;
        }
        XsuperDist_[np] = Xsuper_.size();
      }

      logfileptr->OFS()<<"Supernodes found"<<endl;

      if(options_.relax.nrelax0>0){
        this->relaxSupernodes(ETree_, cc,SupMembership_, Xsuper_, options_.relax );

        logfileptr->OFS()<<"Relaxation done"<<endl;
        //Refresh XsuperDist_
        {
          Idx supPerProc = (Xsuper_.size()-1) / np;
          XsuperDist_.resize(np+1,0);
          for(int p =0; p<np; p++){
            XsuperDist_[p]= p*supPerProc+1;
          }
          XsuperDist_[np] = Xsuper_.size();
        }
      }

      {
        double timeSta = get_time();
        this->symbolicFactorizationRelaxedDist(cc);

        double timeStop = get_time();
        if(iam==0){
          cout<<"Symbolic factorization time: "<<timeStop - timeSta<<endl;
        }
      }

      logfileptr->OFS()<<"Symbfact done"<<endl;
      if(options.order_refinement_str == "SET")
      {
        this->refineSupernodes(2,1);
        if(iam==0){
          cout<<"Supernode reodering done"<<endl;
        }

        if(options_.relax.nrelax0>0){
          this->relaxSupernodes(ETree_, cc,SupMembership_, Xsuper_, options_.relax );
          logfileptr->OFS()<<"Relaxation done"<<endl;
          //Refresh XsuperDist_
          {
            Idx supPerProc = (Xsuper_.size()-1) / np;
            XsuperDist_.resize(np+1,0);
            for(int p =0; p<np; p++){
              XsuperDist_[p]= p*supPerProc+1;
            }
            XsuperDist_[np] = Xsuper_.size();
          }
        }


        {

          double timeSta = get_time();
          this->symbolicFactorizationRelaxedDist(cc);
          double timeStop = get_time();
          if(iam==0){
            cout<<"Symbolic factorization time: "<<timeStop - timeSta<<endl;
          }
        }
        logfileptr->OFS()<<"Symbfact done"<<endl;
      }

  double timeStop = get_time();
  if(iam==0){
    cout<<"Total symbolic factorization time: "<<timeStop - timeSta<<endl;
  }
}

#ifdef _DEBUG_
logfileptr->OFS()<<"Membership list is "<<SupMembership_<<std::endl;
logfileptr->OFS()<<"xsuper "<<Xsuper_<<std::endl;
#endif


#ifdef _OUTPUT_ETREE_
logfileptr->OFS()<<"ETree is "<<ETree_<<std::endl;
logfileptr->OFS()<<"Supernodal ETree is "<<ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_)<<std::endl;
#endif


{
  double timeSta = get_time();

  TIMER_START(LOAD_BALANCE);
  if(options.load_balance_str=="SUBCUBE-FI"){
    if(iam==0){ cout<<"Subtree to subcube FI mapping used"<<endl;}
    ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
    this->Balancer_ = new SubtreeToSubcube(np,SupETree,Xsuper_,XsuperDist_,SupMembership_,locXlindx_,locLindx_,cc,CommEnv_,true);
#ifdef _DEBUG_MAPPING_
    logfileptr->OFS()<<"Proc Mapping: "<<this->Balancer_->GetMap()<<endl;
    logfileptr->OFS()<<"SupETree: "<<SupETree<<endl;
#endif


    TreeMapping * test = dynamic_cast<TreeMapping*>(this->Mapping_);
    if(test==NULL){
      this->Mapping_->Update(this->Balancer_->GetMap());
    }
    else{
      test->Update((TreeLoadBalancer*)this->Balancer_);
    }

  }
  else if(options.load_balance_str=="SUBCUBE-FO"){
    if(iam==0){ cout<<"Subtree to subcube FO mapping used"<<endl;}
    ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
    this->Balancer_ = new SubtreeToSubcube(np,SupETree,Xsuper_,XsuperDist_,SupMembership_,locXlindx_,locLindx_,cc,CommEnv_,false);
#ifdef _DEBUG_MAPPING_
    logfileptr->OFS()<<"Proc Mapping: "<<this->Balancer_->GetMap()<<endl;
    logfileptr->OFS()<<"SupETree: "<<SupETree<<endl;
#endif

    TreeMapping * test = dynamic_cast<TreeMapping*>(this->Mapping_);
    if(test==NULL){
      this->Mapping_->Update(this->Balancer_->GetMap());
    }
    else{
      test->Update((TreeLoadBalancer*)this->Balancer_);
    }


  }
  else if(options.load_balance_str=="SUBCUBE-VOLUME-FI"){
    if(iam==0){ cout<<"Subtree to subcube volume FI mapping used"<<endl;}
    ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
    this->Balancer_ = new SubtreeToSubcubeVolume(np,SupETree,Xsuper_,XsuperDist_,SupMembership_,locXlindx_,locLindx_,cc,CommEnv_,true);
#ifdef _DEBUG_MAPPING_
    logfileptr->OFS()<<"Proc Mapping: "<<this->Balancer_->GetMap()<<endl;
    logfileptr->OFS()<<"SupETree: "<<SupETree<<endl;
#endif

    TreeMapping * test = dynamic_cast<TreeMapping*>(this->Mapping_);
    if(test==NULL){
      this->Mapping_->Update(this->Balancer_->GetMap());
    }
    else{
      test->Update((TreeLoadBalancer*)this->Balancer_);
    }


  }
  else if(options.load_balance_str=="SUBCUBE-VOLUME-FO"){
    if(iam==0){ cout<<"Subtree to subcube volume FO mapping used"<<endl;}
    ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
    this->Balancer_ = new SubtreeToSubcubeVolume(np,SupETree,Xsuper_,XsuperDist_,SupMembership_,locXlindx_,locLindx_,cc,CommEnv_,false);
#ifdef _DEBUG_MAPPING_
    logfileptr->OFS()<<"Proc Mapping: "<<this->Balancer_->GetMap()<<endl;
    logfileptr->OFS()<<"SupETree: "<<SupETree<<endl;
#endif

    TreeMapping * test = dynamic_cast<TreeMapping*>(this->Mapping_);
    if(test==NULL){
      this->Mapping_->Update(this->Balancer_->GetMap());
    }
    else{
      test->Update((TreeLoadBalancer*)this->Balancer_);
    }


  }
  else if(options.load_balance_str=="NNZ"){
    if(iam==0){ cout<<"Load Balancing on NNZ used"<<endl;}
    this->Balancer_ = new NNZBalancer(np,Xsuper_,cc);
#ifdef _DEBUG_MAPPING_
    logfileptr->OFS()<<"Proc Mapping: "<<this->Balancer_->GetMap()<<endl;
#endif
    this->Mapping_->Update(this->Balancer_->GetMap());
  }
  else if(options.load_balance_str=="OLDSUBCUBE"){
    if(iam==0){ cout<<"OLD Subtree to subcube mapping used"<<endl;}
    ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
    //compute number of children and load
    SYMPACK::vector<double> SubTreeLoad(SupETree.Size()+1,0.0);
    SYMPACK::vector<Int> children(SupETree.Size()+1,0);
    for(Int I=1;I<=SupETree.Size();I++){
      Int parent = SupETree.Parent(I-1);
      ++children[parent];
      Int fc = Xsuper_[I-1];
      Int width = Xsuper_[I] - Xsuper_[I-1];
      Int height = cc[fc-1];
      double local_load = width*height*height;
      SubTreeLoad[I]+=local_load;
      SubTreeLoad[parent]+=SubTreeLoad[I];
    }

    logfileptr->OFS()<<"SubTreeLoad is "<<SubTreeLoad<<endl;


    //procmaps[0]/pstart[0] represents the complete list
    SYMPACK::vector<SYMPACK::vector<Int> * > procmaps(SupETree.Size()+1);
    for(Int i = 0; i<procmaps.size();++i){ procmaps[i] = new SYMPACK::vector<Int>();}
    SYMPACK::vector<Int> pstart(SupETree.Size()+1,0);
    procmaps[0]->reserve(np);
    for(Int p = 0;p<np;++p){ procmaps[0]->push_back(p);}


    for(Int I=SupETree.Size(); I>= 1;I--){

      Int parent = SupETree.Parent(I-1);

      //split parent's proc list
      double parent_load = 0.0;

      if(parent!=0){
        Int fc = Xsuper_[parent-1];
        Int width = Xsuper_[parent] - Xsuper_[parent-1];
        Int height = cc[fc-1];
        parent_load = width*height*height;
      }


      double proportion = SubTreeLoad[I]/(SubTreeLoad[parent]-parent_load);
      Int npParent = procmaps[parent]->size();
      Int pFirstIdx = min(pstart[parent],npParent-1);
      //          Int numProcs = max(1,children[parent]==1?npParent-pFirstIdx:min((Int)std::floor(npParent*proportion),npParent-pFirstIdx));
      Int npIdeal =(Int)std::round(npParent*proportion);
      Int numProcs = max(1,min(npParent-pFirstIdx,npIdeal));
      Int pFirst = procmaps[parent]->at(pFirstIdx);

      logfileptr->OFS()<<I<<" npParent = "<<npParent<<" pstartParent = "<<pstart[parent]<<" childrenParent = "<<children[parent]<<" pFirst = "<<pFirst<<" numProcs = "<<numProcs<<" proportion = "<<proportion<<endl; 
      pstart[parent]+= numProcs;//= min(pstart[parent] + numProcs,npParent-1);
      //          children[parent]--;

      //          Int pFirstIdx = pstart[parent]*npParent/children[parent];
      //          Int pFirst = procmaps[parent]->at(pFirstIdx);
      //          Int numProcs = max(pstart[parent]==children[parent]-1?npParent-pFirstIdx:npParent/children[parent],1);



      procmaps[I]->reserve(numProcs);
      for(Int p = pFirst; p<pFirst+numProcs;++p){ procmaps[I]->push_back(p);}

      logfileptr->OFS()<<I<<": "; 
      for(Int i = 0; i<procmaps[I]->size();++i){logfileptr->OFS()<<procmaps[I]->at(i)<<" ";}
      logfileptr->OFS()<<endl;
      //
      //          if(iam==0){
      //            cout<<I<<": "; 
      //            for(Int i = 0; i<procmaps[I]->size();++i){cout<<procmaps[I]->at(i)<<" ";}
      //            cout<<endl;
      //          }
      pstart[parent]++;
    }


    //now choose which processor to get
    SYMPACK::vector<Int> procMap(SupETree.Size());
    SYMPACK::vector<double> load(np,0.0);
    for(Int I=1;I<=SupETree.Size();I++){
      Int minLoadP= -1;
      double minLoad = -1;
      for(Int i = 0; i<procmaps[I]->size();++i){
        Int proc = procmaps[I]->at(i);
        if(load[proc]<minLoad || minLoad==-1){
          minLoad = load[proc];
          minLoadP = proc;
        }
      }

      procMap[I-1] = minLoadP;


      Int fc = Xsuper_[I-1];
      Int width = Xsuper_[I] - Xsuper_[I-1];
      Int height = cc[fc-1];
      double local_load = width*height*height;

      load[minLoadP]+=local_load;
    }


    //for(Int i = 0; i<procMap.size();++i){ logfileptr->OFS()<<i+1<<" is on "<<procMap[i]<<endl;}
    logfileptr->OFS()<<"Proc load: "<<load<<endl;
    //Update the mapping
    logfileptr->OFS()<<"Proc Mapping: "<<procMap<<endl;
    this->Mapping_->Update(procMap);

    for(Int i = 0; i<procmaps.size();++i){ delete procmaps[i];}
  } 
  TIMER_STOP(LOAD_BALANCE);

  double timeStop = get_time();
  if(iam==0){
    cout<<"Load balancing time: "<<timeStop - timeSta<<endl;
  }
}



////    switch(options.load_balance){
////      case SUBCUBE:
////        {
////          LoadBalancer * balancer;
////          if(iam==0){ cout<<"Subtree to subcube mapping used"<<endl;}
////          ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
////          //balancer = new SubtreeToSubcube(np,SupETree,Xsuper_,SupMembership_,xlindx_,lindx_,cc);
////          balancer = new SubtreeToSubcubeVolume(np,SupETree,Xsuper_,SupMembership_,xlindx_,lindx_,cc);
////////
////////
////////          //Int np = 4;
////////
////////          //compute children array and subtree costs
////////
////////          //        if(iam == 0){
////////          //          cout<<"Supernodal ETree is "<<SupETree<<std::endl;
////////          //        }
////////
////////          //compute number of children and load
////////          SYMPACK::vector<double> SubTreeLoad(SupETree.Size()+1,0.0);
////////          SYMPACK::vector<Int> children(SupETree.Size()+1,0);
////////          for(Int I=1;I<=SupETree.Size();I++){
////////            Int parent = SupETree.Parent(I-1);
////////            ++children[parent];
////////            Int fc = Xsuper_[I-1];
////////            Int width = Xsuper_[I] - Xsuper_[I-1];
////////            Int height = cc[fc-1];
////////            double local_load = width*height*height;
////////            SubTreeLoad[I]+=local_load;
////////            SubTreeLoad[parent]+=SubTreeLoad[I];
////////          }
////////
////////          logfileptr->OFS()<<"SubTreeLoad is "<<SubTreeLoad<<endl;
////////
////////
////////          //procmaps[0]/pstart[0] represents the complete list
////////          SYMPACK::vector<SYMPACK::vector<Int> * > procmaps(SupETree.Size()+1);
////////          for(Int i = 0; i<procmaps.size();++i){ procmaps[i] = new SYMPACK::vector<Int>();}
////////          SYMPACK::vector<Int> pstart(SupETree.Size()+1,0);
////////          procmaps[0]->reserve(np);
////////          for(Int p = 0;p<np;++p){ procmaps[0]->push_back(p);}
////////
////////
////////          for(Int I=SupETree.Size(); I>= 1;I--){
////////
////////            Int parent = SupETree.Parent(I-1);
////////
////////            //split parent's proc list
////////            double parent_load = 0.0;
////////
////////            if(parent!=0){
////////              Int fc = Xsuper_[parent-1];
////////              Int width = Xsuper_[parent] - Xsuper_[parent-1];
////////              Int height = cc[fc-1];
////////              parent_load = width*height*height;
////////            }
////////
////////
////////            double proportion = SubTreeLoad[I]/(SubTreeLoad[parent]-parent_load);
////////            Int npParent = procmaps[parent]->size();
////////            Int pFirstIdx = min(pstart[parent],npParent-1);
////////            //          Int numProcs = max(1,children[parent]==1?npParent-pFirstIdx:min((Int)std::floor(npParent*proportion),npParent-pFirstIdx));
////////            Int npIdeal =(Int)std::round(npParent*proportion);
////////            Int numProcs = max(1,min(npParent-pFirstIdx,npIdeal));
////////            Int pFirst = procmaps[parent]->at(pFirstIdx);
////////
////////            logfileptr->OFS()<<I<<" npParent = "<<npParent<<" pstartParent = "<<pstart[parent]<<" childrenParent = "<<children[parent]<<" pFirst = "<<pFirst<<" numProcs = "<<numProcs<<" proportion = "<<proportion<<endl; 
////////            pstart[parent]+= numProcs;//= min(pstart[parent] + numProcs,npParent-1);
////////            //          children[parent]--;
////////
////////            //          Int pFirstIdx = pstart[parent]*npParent/children[parent];
////////            //          Int pFirst = procmaps[parent]->at(pFirstIdx);
////////            //          Int numProcs = max(pstart[parent]==children[parent]-1?npParent-pFirstIdx:npParent/children[parent],1);
////////
////////
////////
////////            procmaps[I]->reserve(numProcs);
////////            for(Int p = pFirst; p<pFirst+numProcs;++p){ procmaps[I]->push_back(p);}
////////
////////            logfileptr->OFS()<<I<<": "; 
////////            for(Int i = 0; i<procmaps[I]->size();++i){logfileptr->OFS()<<procmaps[I]->at(i)<<" ";}
////////            logfileptr->OFS()<<endl;
////////            //
////////            //          if(iam==0){
////////            //            cout<<I<<": "; 
////////            //            for(Int i = 0; i<procmaps[I]->size();++i){cout<<procmaps[I]->at(i)<<" ";}
////////            //            cout<<endl;
////////            //          }
////////            pstart[parent]++;
////////          }
////////
////////
////////          //now choose which processor to get
////////          SYMPACK::vector<Int> procMap(SupETree.Size());
////////          SYMPACK::vector<double> load(np,0.0);
////////          for(Int I=1;I<=SupETree.Size();I++){
////////            Int minLoadP= -1;
////////            double minLoad = -1;
////////            for(Int i = 0; i<procmaps[I]->size();++i){
////////              Int proc = procmaps[I]->at(i);
////////              if(load[proc]<minLoad || minLoad==-1){
////////                minLoad = load[proc];
////////                minLoadP = proc;
////////              }
////////            }
////////
////////            procMap[I-1] = minLoadP;
////////
////////
////////            Int fc = Xsuper_[I-1];
////////            Int width = Xsuper_[I] - Xsuper_[I-1];
////////            Int height = cc[fc-1];
////////            double local_load = width*height*height;
////////
////////            load[minLoadP]+=local_load;
////////          }
////////
////////
////////          //for(Int i = 0; i<procMap.size();++i){ logfileptr->OFS()<<i+1<<" is on "<<procMap[i]<<endl;}
////////          logfileptr->OFS()<<"Proc load: "<<load<<endl;
////          //Update the mapping
////          logfileptr->OFS()<<"Proc Mapping: "<<balancer->GetMap()<<endl;
////          this->Mapping_->Update(balancer->GetMap());
////
////////          for(Int i = 0; i<procmaps.size();++i){ delete procmaps[i];}
////delete balancer;
////        }       
////        break;
////      case SUBCUBE_NNZ:
////        {
////          if(iam==0){ cout<<"Subtree to subcube mapping used"<<endl;}
////
////          //Int np = 4;
////
////          //compute children array and subtree costs
////          ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
////
////          //        if(iam == 0){
////          //          cout<<"Supernodal ETree is "<<SupETree<<std::endl;
////          //        }
////
////          //compute number of children and load
////          SYMPACK::vector<double> SubTreeLoad(SupETree.Size()+1,0.0);
////          SYMPACK::vector<Int> children(SupETree.Size()+1,0);
////          for(Int I=1;I<=SupETree.Size();I++){
////            Int parent = SupETree.Parent(I-1);
////            ++children[parent];
////            Int fc = Xsuper_[I-1];
////            Int width = Xsuper_[I] - Xsuper_[I-1];
////            Int height = cc[fc-1];
////            double local_load = width*height;
////            SubTreeLoad[I]+=local_load;
////            SubTreeLoad[parent]+=SubTreeLoad[I];
////          }
////
////          logfileptr->OFS()<<"SubTreeLoad is "<<SubTreeLoad<<endl;
////
////
////          //procmaps[0]/pstart[0] represents the complete list
////          SYMPACK::vector<SYMPACK::vector<Int> * > procmaps(SupETree.Size()+1);
////          for(Int i = 0; i<procmaps.size();++i){ procmaps[i] = new SYMPACK::vector<Int>();}
////          SYMPACK::vector<Int> pstart(SupETree.Size()+1,0);
////          procmaps[0]->reserve(np);
////          for(Int p = 0;p<np;++p){ procmaps[0]->push_back(p);}
////
////
////          for(Int I=SupETree.Size(); I>= 1;I--){
////
////            Int parent = SupETree.Parent(I-1);
////
////            //split parent's proc list
////            double parent_load = 0.0;
////
////            if(parent!=0){
////              Int fc = Xsuper_[parent-1];
////              Int width = Xsuper_[parent] - Xsuper_[parent-1];
////              Int height = cc[fc-1];
////              parent_load = width*height;
////            }
////
////
////            double proportion = SubTreeLoad[I]/(SubTreeLoad[parent]-parent_load);
////            Int npParent = procmaps[parent]->size();
////            Int pFirstIdx = min(pstart[parent],npParent-1);
////            //          Int numProcs = max(1,children[parent]==1?npParent-pFirstIdx:min((Int)std::floor(npParent*proportion),npParent-pFirstIdx));
////            Int npIdeal =(Int)std::round(npParent*proportion);
////            Int numProcs = max(1,min(npParent-pFirstIdx,npIdeal));
////            Int pFirst = procmaps[parent]->at(pFirstIdx);
////
////            logfileptr->OFS()<<I<<" npParent = "<<npParent<<" pstartParent = "<<pstart[parent]<<" childrenParent = "<<children[parent]<<" pFirst = "<<pFirst<<" numProcs = "<<numProcs<<" proportion = "<<proportion<<endl; 
////            pstart[parent]+= numProcs;//= min(pstart[parent] + numProcs,npParent-1);
////            //          children[parent]--;
////
////            //          Int pFirstIdx = pstart[parent]*npParent/children[parent];
////            //          Int pFirst = procmaps[parent]->at(pFirstIdx);
////            //          Int numProcs = max(pstart[parent]==children[parent]-1?npParent-pFirstIdx:npParent/children[parent],1);
////
////
////
////            procmaps[I]->reserve(numProcs);
////            for(Int p = pFirst; p<pFirst+numProcs;++p){ procmaps[I]->push_back(p);}
////
////            logfileptr->OFS()<<I<<": "; 
////            for(Int i = 0; i<procmaps[I]->size();++i){logfileptr->OFS()<<procmaps[I]->at(i)<<" ";}
////            logfileptr->OFS()<<endl;
////            //
////            //          if(iam==0){
////            //            cout<<I<<": "; 
////            //            for(Int i = 0; i<procmaps[I]->size();++i){cout<<procmaps[I]->at(i)<<" ";}
////            //            cout<<endl;
////            //          }
////            pstart[parent]++;
////          }
////
////
////          //now choose which processor to get
////          SYMPACK::vector<Int> procMap(SupETree.Size());
////          SYMPACK::vector<double> load(np,0.0);
////          for(Int I=1;I<=SupETree.Size();I++){
////            Int minLoadP= -1;
////            double minLoad = -1;
////            for(Int i = 0; i<procmaps[I]->size();++i){
////              Int proc = procmaps[I]->at(i);
////              if(load[proc]<minLoad || minLoad==-1){
////                minLoad = load[proc];
////                minLoadP = proc;
////              }
////            }
////
////            procMap[I-1] = minLoadP;
////
////
////            Int fc = Xsuper_[I-1];
////            Int width = Xsuper_[I] - Xsuper_[I-1];
////            Int height = cc[fc-1];
////            double local_load = width*height;
////
////            load[minLoadP]+=local_load;
////          }
////
////
////          //for(Int i = 0; i<procMap.size();++i){ logfileptr->OFS()<<i+1<<" is on "<<procMap[i]<<endl;}
////          logfileptr->OFS()<<"Proc load: "<<load<<endl;
////          logfileptr->OFS()<<"Proc Mapping: "<<procMap<<endl;
////
////          //Update the mapping
////          this->Mapping_->Update(procMap);
////
////          for(Int i = 0; i<procmaps.size();++i){ delete procmaps[i];}
////        }       
////        break;
////
////      case NNZ:
////        {
////          LoadBalancer * balancer;
////          if(iam==0){ cout<<"Load Balancing on NNZ used"<<endl;}
////          balancer = new NNZBalancer(np,Xsuper_,cc);
////
//////          //Do a greedy load balancing heuristic
//////          SYMPACK::vector<Int> procMap(Xsuper_.size()-1);
//////          //Do a greedy heuristic to balance the number of nnz ?
//////          SYMPACK::vector<double> load(np,0.0);
//////
//////          for(Int i = 1; i< Xsuper_.size();  ++i){
//////            //find least loaded processor
//////            SYMPACK::vector<double>::iterator it = std::min_element(load.begin(),load.end());
//////            Int proc = (Int)(it - load.begin());
//////            Int width = Xsuper_[i] - Xsuper_[i-1];
//////            Int height = cc[i-1];
//////            *it += (double)(width*height);
//////            procMap[i-1] = proc;
//////          } 
//////
//////          logfileptr->OFS()<<"Proc load: "<<load<<endl;
//////          logfileptr->OFS()<<"Proc Mapping: "<<procMap<<endl;
////
////          //Update the mapping
//////          this->Mapping_->Update(procMap);
////
////          logfileptr->OFS()<<"Proc Mapping: "<<balancer->GetMap()<<endl;
////          this->Mapping_->Update(balancer->GetMap());
////
////////          for(Int i = 0; i<procmaps.size();++i){ delete procmaps[i];}
////delete balancer;
////        }
////        break;
////      case FLOPS:
////        {
////          if(iam==0){ cout<<"Load Balancing on FLOPS used"<<endl;}
////          //Do a greedy load balancing heuristic
////          SYMPACK::vector<Int> procMap(Xsuper_.size()-1);
////          //Do a greedy heuristic to balance the number of nnz ?
////          SYMPACK::vector<double> load(np,0.0);
////
////          for(Int i = 1; i< Xsuper_.size();  ++i){
////            //find least loaded processor
////            SYMPACK::vector<double>::iterator it = std::min_element(load.begin(),load.end());
////            Int proc = (Int)(it - load.begin());
////            Int width = Xsuper_[i] - Xsuper_[i-1];
////            Int height = cc[i-1];
////            *it += (double)(height*height);
////            procMap[i-1] = proc;
////          } 
////
////          logfileptr->OFS()<<"Proc load: "<<load<<endl;
////          logfileptr->OFS()<<"Proc Mapping: "<<procMap<<endl;
////
////          //Update the mapping
////          this->Mapping_->Update(procMap);
////        }
////        break;
////
////      default:
////        {
////        }
////        break;
////    }
////


#ifdef _DEBUG_MAPPING_
//    this->Mapping_->Dump(2*np);
#endif

SYMPACK::vector<Int> numBlk;
TIMER_START(Get_UpdateCount);
GetUpdatingSupernodeCount(UpdateCount_,UpdateWidth_,UpdateHeight_,numBlk);
TIMER_STOP(Get_UpdateCount);




Int iam = CommEnv_->MPI_Rank();
Int np  = CommEnv_->MPI_Size();


//copy

{
  double timeSta = get_time();
  TIMER_START(DISTRIBUTING_MATRIX);
  //Icomm buffer;
  Icomm recv_buffer;


  AsyncComms incomingRecv;
  AsyncComms outgoingSend;
  Int numColFirst = std::max(1,iSize_ / np);

  logfileptr->OFS()<<"Starting Send"<<endl;
  try{
    TIMER_START(DISTRIBUTE_COUNTING);
    //first, count 
    map<Int,pair<size_t,Icomm *> > send_map;

    //make a copy of colptr to store the current position in each column

    Idx firstCol = iam*(iSize_/np);

    for(Int p=0;p<np;++p){
      send_map[p].first = 0;
    }

    Int snodeCount = 0;
    for(Int I=1;I<Xsuper_.size();I++){
      Int fc = Xsuper_[I-1];
      Int lc = Xsuper_[I]-1;
      Int iWidth = lc-fc+1;
      Int iHeight = UpdateHeight_[I-1];

      Int iDest = this->Mapping_->Map(I-1,I-1);

      if(iDest==iam){
        ++snodeCount;
      }

      //look at the owner of the first column of the supernode
      Int numColFirst = std::max(1,iSize_ / np);

      //post all the recv and sends
      for(Int col = fc;col<=lc;col++){
        //corresponding column in the unsorted matrix A
        Int orig_col = Order_.perm[col-1];
        Int iOwnerCol = std::min((orig_col-1)/numColFirst,np-1);
        size_t & send_bytes = send_map[iDest].first;
        if(iam == iOwnerCol){
          Int nrows = 0;
          Int local_col = (orig_col-(numColFirst)*iOwnerCol);
          for(Int rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
            Int orig_row = pMat.Local_.rowind[rowidx-1];
            Int row = Order_.invp[orig_row-1];

            if(row<col){
              //add the pair (col,row) to processor owning column row
              Int J = SupMembership_[row-1];
              Int iDestJ = this->Mapping_->Map(J-1,J-1);
              send_map[iDestJ].first += sizeof(Int)+sizeof(Int)+1*(sizeof(Int)+sizeof(T));
            }

            if(row>=col){
              //add the pair (row,col) to iDest
              nrows++;
            }
          }
          send_bytes += sizeof(Int)+sizeof(Int)+nrows*(sizeof(Int)+sizeof(T));
        }
      }
    }

    TIMER_STOP(DISTRIBUTE_COUNTING);

    //Resize the local supernodes array
    LocalSupernodes_.reserve(snodeCount);


    TIMER_START(DISTRIBUTE_CREATE_SNODES);

    for(Int I=1;I<Xsuper_.size();I++){

      Int iDest = this->Mapping_->Map(I-1,I-1);

      //parse the first column to create the supernode structure
      if(iam==iDest){
        Int fc = Xsuper_[I-1];
        Int lc = Xsuper_[I]-1;
        Int iWidth = lc-fc+1;
        Int iHeight = UpdateHeight_[I-1];
        Int nzBlockCnt = numBlk[I-1];
#ifndef ITREE2
        globToLocSnodes_.push_back(I-1);
#else 
        ITree::Interval snode_inter = { I, I, LocalSupernodes_.size() };
        globToLocSnodes_.Insert(snode_inter);
#endif

//#ifdef _INDEFINITE_
//        LocalSupernodes_.push_back( new SuperNodeInd<T>(I,fc,lc,iHeight,iSize_,nzBlockCnt));
//#else
//        LocalSupernodes_.push_back( new SuperNode<T>(I,fc,lc,iHeight,iSize_,nzBlockCnt));
//#endif
        LocalSupernodes_.push_back( CreateSuperNode(options_.decomposition,I,fc,lc,iHeight,iSize_,nzBlockCnt));
      }
    }

    TIMER_STOP(DISTRIBUTE_CREATE_SNODES);

    //first create the structure of every supernode
    {
      std::vector< int > superStructure(np);
      std::vector< int > rdisplsStructure;
      std::vector< int > sSuperStructure(np);
      std::vector< int > ssizes(np,0);
      std::vector< int > sdispls(np+1,0);
      std::vector< int > rsizes(np,0);

      //      Int numLocSnode = ( (Xsuper_.size()-1) / np);
      //      Int firstSnode = iam*numLocSnode + 1;

      Int numLocSnode = XsuperDist_[iam+1]-XsuperDist_[iam];
      Int firstSnode = XsuperDist_[iam];


      for(Int locsupno = 1; locsupno<locXlindx_.size(); ++locsupno){
        Idx I = locsupno + firstSnode-1;
        Int iDest = this->Mapping_->Map(I-1,I-1);
        ssizes[iDest] += 1 + 2*numBlk[I-1]; //1 for supno index + numBlk startrows + numBlk number of rows
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

      MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,CommEnv_->MPI_GetComm());
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
          CommEnv_->MPI_GetComm());

      //loop through received structure and create supernodes
      for(Int p = 0; p<np; p++){
        int pos = rdisplsStructure[p]/sizeof(int);
        int end = rdisplsStructure[p+1]/sizeof(int);
        while(pos<end){
          Int I = superStructure[pos++];
          Int nzBlockCnt = numBlk[I-1];

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



    {
      //allocate one buffer for every remote processor
      //compute the send structures
      size_t total_send_size = 0;
      SYMPACK::vector<int> sdispls(np+1,0);
      SYMPACK::vector<int> scounts(np,0);
      for(auto it = send_map.begin(); it!=send_map.end();it++){
        Int iCurDest = it->first;
        size_t & send_bytes = it->second.first;
        scounts[iCurDest] = send_bytes;
#ifdef _DEBUG_
        logfileptr->OFS()<<"P"<<iam<<" ---- "<<send_bytes<<" ---> P"<<iCurDest<<endl;
#endif
      }


      //compute send displacements
      sdispls[0] = 0;
      std::partial_sum(scounts.begin(),scounts.end(),&sdispls[1]);
      total_send_size = sdispls.back();

      auto sposition = sdispls;
      Icomm * IsendPtr = new Icomm(total_send_size,MPI_REQUEST_NULL);

      //Fill up the send buffer
      for(Int I=1;I<Xsuper_.size();I++){
        Int fc = Xsuper_[I-1];
        Int lc = Xsuper_[I]-1;
        Int iWidth = lc-fc+1;
        //  Ptr fi = xlindx_[I-1];
        //  Ptr li = xlindx_[I]-1;
        Int iHeight = cc[fc-1];//li-fi+1;

        Icomm & Isend = *IsendPtr;
        Int iDest = this->Mapping_->Map(I-1,I-1);

#ifdef _DEBUG_
        logfileptr->OFS()<<"Supernode "<<I<<" is owned by P"<<iDest<<std::endl;
#endif

        //look at the owner of the first column of the supernode
        Int numColFirst = std::max(1,iSize_ / np);

        //Serialize
        for(Int col = fc;col<=lc;col++){
          Int orig_col = Order_.perm[col-1];
          Int iOwnerCol = std::min((orig_col-1)/numColFirst,np-1);
          if(iam == iOwnerCol){
            Int nrows = 0;
            Int local_col = (orig_col-(numColFirst)*iOwnerCol);
            for(Int rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
              Int orig_row = pMat.Local_.rowind[rowidx-1];
              Int row = Order_.invp[orig_row-1];


              if(row<col){
                //add the pair (col,row) to processor owning column row
                Int J = SupMembership_[row-1];
                Int iDestJ = this->Mapping_->Map(J-1,J-1);

                T val = pMat.nzvalLocal[rowidx-1];

                //we need to set head to the proper sdispls
                Isend.setHead(sposition[iDestJ]);
                Isend<<row;
                Isend<<1;
                Isend<<col;
                Isend<<val;
                //backup new position for processor iDestJ
                sposition[iDestJ] = Isend.head;
              }


              if(row>=col){
                nrows++;
              }
            }


            //restore head
            Isend.setHead(sposition[iDest]);
            Isend<<col;
            Isend<<nrows;

            for(Int rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
              Int orig_row = pMat.Local_.rowind[rowidx-1];
              Int row = Order_.invp[orig_row-1];
              if(row>=col){
                Isend<<row;
              }
            }

            for(Int rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
              Int orig_row = pMat.Local_.rowind[rowidx-1];
              Int row = Order_.invp[orig_row-1];
              if(row>=col){
                Isend<<pMat.nzvalLocal[rowidx-1];
              }
            }

            //backup new position for processor iDest
            sposition[iDest] = Isend.head;


          }
        }
      }

      //compute the receive structures
      size_t total_recv_size = 0;
      SYMPACK::vector<int> rdispls(np+1,0);
      SYMPACK::vector<int> rcounts(np,0);

      MPI_Alltoall(&scounts[0],sizeof(int),MPI_BYTE,&rcounts[0],sizeof(int),MPI_BYTE,CommEnv_->MPI_GetComm());

      //compute receive displacements
      rdispls[0] = 0;
      std::partial_sum(rcounts.begin(),rcounts.end(),&rdispls[1]);
      total_recv_size = rdispls.back();

      Icomm * IrecvPtr = new Icomm(total_recv_size,MPI_REQUEST_NULL);


#ifdef _DEBUG_
      logfileptr->OFS()<<"scounts: "<<scounts<<endl;
      logfileptr->OFS()<<"sdispls: "<<sdispls<<endl;
      logfileptr->OFS()<<"rcounts: "<<rcounts<<endl;
      logfileptr->OFS()<<"rdispls: "<<rdispls<<endl;
#endif

      MPI_Alltoallv(IsendPtr->front(), &scounts[0], &sdispls[0], MPI_BYTE,
          IrecvPtr->front(), &rcounts[0], &rdispls[0], MPI_BYTE,
          CommEnv_->MPI_GetComm());

      //Need to parse the structure sent from the processor owning the first column of the supernode

      for(Int p=0;p<np;++p){
        IrecvPtr->setHead(rdispls[p]);
        while(IrecvPtr->head < rdispls[p]+rcounts[p]){ 
          char * buffer = IrecvPtr->back();

          //Deserialize
          Int col = *(Int*)&buffer[0];
          Int I = SupMembership_[col-1];

          Int iDest = this->Mapping_->Map(I-1,I-1);
          if(iam != iDest){gdb_lock();}

          assert(iam==iDest);

          Int fc = Xsuper_[I-1];
          Int lc = Xsuper_[I]-1;
          Int iWidth = lc-fc+1;

          SuperNode<T> * snode = snodeLocal(I);

          //nrows of column col sent by processor p
          Int nrows = *((Int*)&buffer[sizeof(Int)]);
          Int * rowind = (Int*)(&buffer[2*sizeof(Int)]);
          T * nzvalA = (T*)(&buffer[(2+nrows)*sizeof(Int)]);
          //advance in the buffer
          IrecvPtr->setHead(IrecvPtr->head + 2*sizeof(Int) + nrows*(sizeof(Int)+sizeof(T)));

          Int colbeg = 1;
          Int colend = nrows;
          for(Int rowidx = colbeg; rowidx<=colend; ++rowidx){
            Int row = rowind[rowidx-1];
            //logfileptr->OFS()<<nzvalA[rowidx-1]<<" ";
            //Int orig_row = rowind[rowidx-1];
            //Int row = Order_.invp[orig_row-1];
            //if(row>=col){
            Int blkidx = snode->FindBlockIdx(row);
            assert(blkidx!=-1);
            NZBlockDesc & blk_desc = snode->GetNZBlockDesc(blkidx);
            Int local_row = row - blk_desc.GIndex + 1;
            Int local_col = col - fc + 1;
            T * nzval = snode->GetNZval(blk_desc.Offset);
            nzval[(local_row-1)*iWidth+local_col-1] = nzvalA[rowidx-1];
            //}
          }
        }
      }


      //after the alltoallv, cleanup
      delete IrecvPtr;
      delete IsendPtr;
    }

  }
  catch(const std::bad_alloc& e){
    std::cout << "Allocation failed: " << e.what() << '\n';
    abort();
  }
  logfileptr->OFS()<<"Send Done"<<endl;

  TIMER_STOP(DISTRIBUTING_MATRIX);
  double timeStop = get_time();
  if(iam==0){
    cout<<"Distribution time: "<<timeStop - timeSta<<endl;
  }
}

#ifdef _DEBUG_
for(Int I=1;I<Xsuper_.size();I++){
  Int src_first_col = Xsuper_[I-1];
  Int src_last_col = Xsuper_[I]-1;
  Int iOwner = this->Mapping_->Map(I-1,I-1);
  //If I own the column, factor it
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
#endif
//Dump();

#if 0

//allocate backup buffer
{
  auto maxwit = max_element(&UpdateWidth_[0],&UpdateWidth_[0]+UpdateWidth_.size());
  auto maxhit = max_element(&UpdateHeight_[0],&UpdateHeight_[0]+UpdateHeight_.size());

  int maxh = *maxhit;
  int maxw = *maxwit;
  for(Int I=1;I<Xsuper_.size();++I){
    int width = Xsuper_[I] - Xsuper_[I-1];
    maxw = max(maxw,width);
  }


  Int max_bytes = 6*sizeof(Int); 
  Int nrows = maxh;
  Int ncols = maxw;
  Int nz_cnt = nrows * ncols;
  Int nblocks = nrows;
  max_bytes += (nblocks)*sizeof(NZBlockDesc);
  max_bytes += nz_cnt*sizeof(T);

  bool success = backupBuffer_.AllocLocal(max_bytes);
  assert(success);

}

#endif


#ifdef PREALLOC_IRECV
#endif



//delete Local_;
//delete Global_;
#endif




origTaskLists_.resize(Xsuper_.size(),NULL);
SYMPACK::vector<Int> AggregatesToRecv;
SYMPACK::vector<Int> LocalAggregates;
FBGetUpdateCount(UpdatesToDo_,AggregatesToRecv,LocalAggregates);

generateTaskGraph(taskGraph_, AggregatesToRecv, LocalAggregates);

//    generateTaskGraph(localTaskCount_, origTaskLists_);
//    for(int i = 0; i<origTaskLists_.size(); ++i){
//      if(origTaskLists_[i] != NULL){
//        for(auto taskit = origTaskLists_[i]->begin(); taskit!=origTaskLists_[i]->end();taskit++){
//          FBTask & curUpdate = *taskit;
//          if(curUpdate.type==FACTOR){
//            //      curUpdate.rank = levels[curUpdate.src_snode_id-1];
//            curUpdate.remote_deps = AggregatesToRecv[curUpdate.tgt_snode_id-1];
//            curUpdate.local_deps = LocalAggregates[curUpdate.tgt_snode_id-1];
//          }
//          else if(curUpdate.type==UPDATE){
//            Int iOwner = Mapping_->Map(curUpdate.src_snode_id-1,curUpdate.src_snode_id-1);
//            //        curUpdate.rank = levels[curUpdate.src_snode_id-1];
//            //If I own the factor, it is a local dependency
//            if(iam==iOwner){
//              curUpdate.remote_deps = 0;
//              curUpdate.local_deps = 1;
//            }
//            else{
//              curUpdate.remote_deps = 1;
//              curUpdate.local_deps = 0;
//            }
//          }
//        }
//      }
//    }
//
//    SYMPACK::vector<Int> levels;
//    levels.resize(Xsuper_.size());
//    levels[0]=0;
//    Int numLevel = 0; 
//    for(Int i=Xsuper_.size()-1-1; i>=0; i-- ){     
//      Int fcol = Xsuper_[i];
//      bassert(fcol-1>=0);
//      Int pcol =ETree_.PostParent(fcol-1);
//      Int supno = pcol>0?SupMembership_[pcol-1]:0;
//      levels[i] = levels[supno]+1;
//    }
//
//    //sort the task queues
//    {
//      for(int i = 0; i<origTaskLists_.size(); ++i){
//        if(origTaskLists_[i] != NULL){
//          for(auto taskit = origTaskLists_[i]->begin(); taskit!=origTaskLists_[i]->end();taskit++){
//            if(taskit->remote_deps==0 && taskit->local_deps==0){
//              taskit->rank = levels[taskit->src_snode_id];
//              scheduler_->push(taskit);
//            }
//          }
//        }
//      }
//    }







MPI_Barrier(CommEnv_->MPI_GetComm());


TIMER_STOP(SUPERMATRIX_INIT);

}




template <typename T> SupernodalMatrix<T>::SupernodalMatrix(){
  CommEnv_=NULL;
  Mapping_ = NULL;
  Balancer_ = NULL;
  scheduler_ = NULL;
  scheduler2_ = NULL;
  Local_=NULL;
  Global_=NULL;
  isGlobStructAllocated_ = false;
}

template <typename T> SupernodalMatrix<T>::SupernodalMatrix(DistSparseMatrix<T> & pMat, NGCholOptions & options ){
  CommEnv_ = NULL;
  Mapping_ = NULL;
  Balancer_ = NULL;
  scheduler_ = NULL;
  scheduler2_ = NULL;
  Local_=NULL;
  Global_=NULL;
  isGlobStructAllocated_ = false;
  Init(pMat, options);
}

template <typename T> SupernodalMatrix<T>::~SupernodalMatrix(){
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

  if(Local_!=NULL){
    delete Local_;
  }

  if(Global_!=NULL){
    delete Global_;
  }

#ifdef MULTITHREADING
  std::for_each(superLocks_.begin(), superLocks_.end(), [] (omp_lock_t & lock) { omp_destroy_lock(&lock);});
#endif



}


//returns the 1-based index of supernode id global in the local supernode array
template <typename T> Int SupernodalMatrix<T>::snodeLocalIndex(Int global){
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
template <typename T> SuperNode<T> * SupernodalMatrix<T>::snodeLocal(Int global){
  Int iLocal = snodeLocalIndex(global);
  return LocalSupernodes_[iLocal -1];
}

template <typename T> 
template< class Alloc>
SuperNode<T,Alloc> * SupernodalMatrix<T>::snodeLocal(Int global, SYMPACK::vector<SuperNode<T,Alloc> *> & snodeColl){
  Int iLocal = snodeLocalIndex(global);
  return snodeColl[iLocal -1];
}


template <typename T> 
void SupernodalMatrix<T>::findSupernodes(ETree& tree, Ordering & aOrder, SYMPACK::vector<Int> & cc,SYMPACK::vector<Int> & supMembership, SYMPACK::vector<Int> & xsuper, Int maxSize ){
  TIMER_START(FindSupernodes);
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
        if(supsize<maxSize || maxSize==-1){
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
  TIMER_STOP(FindSupernodes);
}



template <typename T> 
void SupernodalMatrix<T>::getLColRowCount(SparseMatrixGraph & sgraph, SYMPACK::vector<Int> & cc, SYMPACK::vector<Int> & rc){
  //The tree need to be postordered
  if(!ETree_.IsPostOrdered()){
    ETree_.PostOrderTree(Order_);
  }

  if(iam == 0 && (!sgraph.IsExpanded() ) ){
    throw std::logic_error( "SparseMatrixGraph must be expanded\n" );
  }

  sgraph.SetBaseval(1);
  sgraph.SetKeepDiag(1);

  TIMER_START(GetColRowCount_Classic);

  Int size = sgraph.size;

  if(iam==0){
    cc.resize(size);
    rc.resize(size);
    SYMPACK::vector<Idx> level(size+1);
    SYMPACK::vector<Int> weight(size+1);
    SYMPACK::vector<Idx> fdesc(size+1);
    SYMPACK::vector<Idx> nchild(size+1);
    SYMPACK::vector<Idx> set(size);
    SYMPACK::vector<Idx> prvlf(size);
    SYMPACK::vector<Idx> prvnbr(size);

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

  team_->barrier();
  if(iam!=0){
    cc.resize(size);
    rc.resize(size);
  }
  //Broadcast  
  team_->bcast(&cc[0], &cc[0], size*sizeof(Int), 0);
  team_->bcast(&rc[0], &rc[0], size*sizeof(Int), 0);


  TIMER_STOP(GetColRowCount_Classic);
}


template <typename T> 
void SupernodalMatrix<T>::getLColRowCount(SYMPACK::vector<Int> & cc, SYMPACK::vector<Int> & rc){
  //The tree need to be postordered
  if(!ETree_.IsPostOrdered()){
    ETree_.PostOrderTree(Order_);
  }

  graph_.ExpandSymmetric();

  TIMER_START(GetColRowCount_Classic);

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
    team_->barrier();
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
    team_->barrier();
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

  team_->barrier();
#if 0
  logfileptr->OFS()<<"nchild: ";
  for(Idx k =0; k<=size; ++k){
    logfileptr->OFS()<<nchild[k]<<" ";
  }
  logfileptr->OFS()<<endl;
  logfileptr->OFS()<<"fdesc: ";
  for(Idx k =0; k<=size; ++k){
    logfileptr->OFS()<<fdesc[k]<<" ";
  }
  logfileptr->OFS()<<endl;
  logfileptr->OFS()<<"weight: ";
  for(Idx k =0; k<=size; ++k){
    logfileptr->OFS()<<weight[k]<<" ";
  }
  logfileptr->OFS()<<endl;
  logfileptr->OFS()<<"level: ";
  for(Idx k =0; k<=size; ++k){
    logfileptr->OFS()<<level[k]<<" ";
  }
  logfileptr->OFS()<<endl;
  logfileptr->OFS()<<"prvnbr: ";
  for(Idx k =0; k<size; ++k){
    logfileptr->OFS()<<prvnbr[k]<<" ";
  }
  logfileptr->OFS()<<endl;
  logfileptr->OFS()<<"prvlf: ";
  for(Idx k =0; k<size; ++k){
    logfileptr->OFS()<<prvlf[k]<<" ";
  }
  logfileptr->OFS()<<endl;
  logfileptr->OFS()<<"set: ";
  for(Idx k =0; k<size; ++k){
    logfileptr->OFS()<<set[k]<<" ";
  }
  logfileptr->OFS()<<endl;
  logfileptr->OFS()<<"rcDist: ";
  for(Idx k =0; k<size; ++k){
    logfileptr->OFS()<<rcDist[k]<<" ";
  }
  logfileptr->OFS()<<endl;
  team_->barrier();
#endif

  //foreach column in the permuted matrix
  for(Idx lownbr = 1; lownbr<=size; ++lownbr){
    team_->barrier();
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
              logfileptr->OFS()<<endl;


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


  team_->barrier();


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
  team_->barrier();

  //      logfileptr->OFS()<<"column counts "<<cc<<std::endl;

  TIMER_STOP(GetColRowCount_Classic);
}





template <typename T> 
void SupernodalMatrix<T>::relaxSupernodes(ETree& tree, SYMPACK::vector<Int> & cc,SYMPACK::vector<Int> & supMembership, SYMPACK::vector<Int> & xsuper, RelaxationParameters & params ){
  //todo tree cc supmembership xsuper and relax params are members, no need for arguments
  Int nsuper = xsuper.size()-1;

  DisjointSet sets;
  sets.Initialize(nsuper);
  SYMPACK::vector<Int> ncols(nsuper);
  SYMPACK::vector<Int> zeros(nsuper);
  SYMPACK::vector<Int> newCC(nsuper);
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

    Int parent_fstcol = tree.PostParent(lstcol-1);
    if(parent_fstcol!=0){
      Int parent_snode = supMembership[parent_fstcol-1];
      Int pset = sets.find(parent_snode);
      parent_snode = sets.Root(pset-1);

      bool merge = (parent_snode == ksup+1);


      if(merge){
        Int parent_width = ncols[parent_snode-1];

        Int parent_fstcol = xsuper[parent_snode-1];
        Int parent_lstcol = xsuper[parent_snode]-1;
        Int totzeros = zeros[parent_snode-1];
        Int fused_cols = width + parent_width;

        merge = false;
        if(fused_cols <= params.nrelax0){
          merge = true;
        }
        else if(fused_cols <=params.maxSize || params.maxSize==-1){
          double child_lnz = cc[fstcol-1];
          double parent_lnz = cc[parent_fstcol-1];
          double xnewzeros = width * (parent_lnz + width  - child_lnz);

          if(xnewzeros == 0){
            merge = true;
          }
          else{
            //all these values are the values corresponding to the merged snode
            double xtotzeros = (double)totzeros + xnewzeros;
            double xfused_cols = (double) fused_cols;
            //new number of nz
            double xtotsize = (xfused_cols * (xfused_cols+1)/2) + xfused_cols * (parent_lnz - parent_width);
            //percentage of explicit zeros
            double z = xtotzeros / xtotsize;

            Int totsize = (fused_cols * (fused_cols+1)/2) + fused_cols * ((Int)parent_lnz - parent_width);
            totzeros += (Int)xnewzeros;

            merge = ((fused_cols <= params.nrelax1 && z < params.zrelax0) 
                || (fused_cols <= params.nrelax2 && z < params.zrelax1)
                || (z<params.zrelax2)) &&
              (xtotsize < std::numeric_limits<Int>::max() / sizeof(double));
          }

        }

        if(merge){
          //              std::cout<<"merge "<<ksup<<" and "<<parent_snode<<std::endl;
          ncols[ksup-1] += ncols[parent_snode-1]; 
          zeros[ksup-1] = totzeros;
          newCC[ksup-1] = width + newCC[parent_snode-1];
          sets.Union(ksup,parent_snode,ksup);
        }
      } 

    }
  }

  SYMPACK::vector<Int> relXSuper(nsuper+1);
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
void SupernodalMatrix<T>::symbolicFactorizationRelaxedDist(SYMPACK::vector<Int> & cc){
  TIMER_START(SymbolicFactorization);
  Int size = iSize_;
  ETree& tree = ETree_;
  Ordering & aOrder = Order_;
  DistSparseMatrixGraph & graph = graph_;
  SYMPACK::vector<Int> & xsuper = Xsuper_;
  SYMPACK::vector<Int> & SupMembership = SupMembership_;
  PtrVec & xlindx = locXlindx_;
  IdxVec & lindx = locLindx_;
  MPI_Comm & comm = CommEnv_->MPI_GetComm();



  //permute the graph
  DistSparseMatrixGraph pGraph = graph;

#if 0
  {
    double tstart = get_time();
    pGraph.Permute(&Order_.invp[0]);
    double tstop = get_time();
    logfileptr->OFS()<<"Permute time: "<<tstop-tstart<<endl;
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
      Idx pOwnerFirst = min((Idx)np-1, (fstcol-1) / colPerProc);
      Idx pOwnerLast = min((Idx)np-1, (lstcol-1) / colPerProc);

      if(pOwnerFirst!=pOwnerLast){
        numSplits += pOwnerLast-pOwnerFirst;
      }
    }

    SYMPACK::vector<Int> newSnodes(np,0);
    SYMPACK::vector<Int> newXsuper(xsuper.size()+numSplits);
    Int pos = 0;
    for(Int snode = 1; snode<=xsuper.size()-1;snode++){
      Idx fstcol = xsuper[snode-1];
      Idx lstcol = xsuper[snode]-1;
      //check if these two columns are on the same processor
      Idx pOwnerFirst = min((Idx)np-1, (fstcol-1) / colPerProc);
      Idx pOwnerLast = min((Idx)np-1, (lstcol-1) / colPerProc);

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
#if 1
    {
      double tstart = get_time();

      //recompute vertexDist based on XsuperDist
      std::vector<Idx> newVertexDist(np+1);
      for(int p = 0; p < np; p++){
        newVertexDist[p] = Xsuper_[XsuperDist_[p]-1] + (pGraph.GetBaseval()-1);
      }
      newVertexDist[np] = pGraph.size+1;


      pGraph.Permute(&Order_.invp[0],&newVertexDist[0]);
      double tstop = get_time();
      logfileptr->OFS()<<"Permute time: "<<tstop-tstart<<endl;
    }
#else
    {
      double tstart = get_time();
      pGraph.Permute(&Order_.invp[0]);
      double tstop = get_time();
      logfileptr->OFS()<<"Permute time: "<<tstop-tstart<<endl;
    }
    double tstart = get_time();
    //redistribute the graph according to the supernodal partition
    pGraph.RedistributeSupernodal(xsuper.size()-1,&xsuper[0],&XsuperDist_[0],&SupMembership[0]);
    double tstop = get_time();
    logfileptr->OFS()<<"Redistribute time: "<<tstop-tstart<<endl;
#endif
  }
#endif

  //logfileptr->OFS()<<Xsuper_<<endl;

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
  SYMPACK::vector<Idx> mrglnk(nsuper,0);

  //Array of length n+1 containing the current linked list 
  //of merged indices (the "reach" set)
  SYMPACK::vector<Idx> rchlnk(size+1);

  //Array of length n used to mark indices as they are introduced
  // into each supernode's index set
  SYMPACK::vector<Int> marker(size,0);

  Int nsuperLocal = XsuperDist_[iam+1]-XsuperDist_[iam];//(iam!=np-1)?nsuper/np:nsuper-iam*(nsuper/np);
  Int firstSnode = XsuperDist_[iam];//iam*(nsuper/np)+1;
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
  //logfileptr->OFS()<<xlindx<<endl;

  SYMPACK::vector<Ptr> recvXlindx;
  SYMPACK::vector<Idx> recvLindx;

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

#if 0
  Idx colPerProc = size / np;
  for(Int locksup = 1; locksup<=nsuperLocal; ++locksup){
    Int ksup = locksup + firstSnode - 1;
    Int fstcol = xsuper[ksup-1];
    Int lstcol = xsuper[ksup]-1;

    Idx pOwnerFirst = min((Idx)np-1, (fstcol-1) / colPerProc);
    Idx pOwnerLast = min((Idx)np-1, (lstcol-1) / colPerProc);
    bassert(pOwnerFirst==pOwnerLast && iam==pOwnerFirst);
  }
#endif

  //logfileptr->OFS()<<"P"<<iam<<" is active"<<endl;
  //logfileptr->OFS()<<mrglnk<<endl;
  //logfileptr->OFS()<<rchlnk<<endl;


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

        //logfileptr->OFS()<<jsup<<endl;
        Int jwidth = xsuper[jsup]-xsuper[jsup-1];

        Ptr * jxlindx = NULL;
        Idx * jlindx = NULL;
        Int locjsup = -1;
        if(jsup>=firstSnode && jsup<=lastSnode){
          locjsup = jsup - firstSnode +1;
          jxlindx = &xlindx[0];
          jlindx = &lindx[0];
        }
        else{
          MPI_Status status;
          recvLindx.resize(size);
          //receive jsup lindx
          //Int psrc = min( (Int)np-1, (jsup-1) / (nsuper/np) );
          Int psrc = 0; for(psrc = 0; psrc<iam;psrc++){ if(XsuperDist_[psrc]<=jsup && jsup<XsuperDist_[psrc+1]){ break; } }
          //logfileptr->OFS()<<"trying to recv "<<jsup<<" max "<<recvLindx.size()*sizeof(Idx)<<" bytes"<<" from P"<<psrc<<endl;
          MPI_Recv(&recvLindx[0],recvLindx.size()*sizeof(Idx),MPI_BYTE,psrc,jsup+np,comm,&status);
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
          MPI_Irecv(&recvLindx[0],recvLindx.size()*sizeof(Idx),MPI_BYTE,psrc,jsup+np,comm,&request);
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
      //      logfileptr->OFS()<<"col "<<fstcol<<" : ";
      //      for(Ptr kptr = nzend+1; kptr<=nzend+knz;++kptr){
      //        i = rchlnk[i];
      //        logfileptr->OFS()<<i<<" ";
      //      } 
      //      logfileptr->OFS()<<endl;


      for(Int col = fstcol; col<=lstcol; ++col){
        Int local_col = col - firstColumn + 1;

        //        logfileptr->OFS()<<"     col "<<col<<" : ";
        Ptr knzbeg = pGraph.colptr[local_col-1];
        Ptr knzend = pGraph.colptr[local_col]-1;
        for(Ptr kptr = knzbeg; kptr<=knzend;++kptr){
          Idx newi = pGraph.rowind[kptr-1];
          //          if(newi > fstcol){
          //            logfileptr->OFS()<<newi<<" ";
          //          }
        }
        //        logfileptr->OFS()<<endl;
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
      //        Idx pcol = lindx[xlindx[locksup-1] + width -1];
      //logfileptr->OFS()<<pcol2<<" vs "<<pcol<<endl;
      //assert(pcol==pcol2);
      Int psup = SupMembership[pcol-1];
      mrglnk[ksup-1] = mrglnk[psup-1];
      mrglnk[psup-1] = ksup;


      //send L asap
      //Int pdest = min( (Int)np-1, (psup-1) / (nsuper/np) );

      Int pdest = 0; for(pdest = 0; pdest<np;pdest++){ if(XsuperDist_[pdest]<=psup && psup<XsuperDist_[pdest+1]){ break; } }
      //if remote
      if(pdest!=iam){
        mpirequests.push_back(MPI_REQUEST_NULL);
        MPI_Request & request = mpirequests.back();
        //logfileptr->OFS()<<"sending "<<length*sizeof(Idx)<<" bytes of "<<ksup<<" to P"<<pdest<<" for "<<ksup<<endl;
        MPI_Isend(&lindx[xlindx[locksup-1]-1],length*sizeof(Idx),MPI_BYTE,pdest,ksup+np,comm,&request);
      }
    }

  }



  //lindx.resize(nzend+1);
  lindx.resize(nzend);
  MPI_Barrier(comm);
  TIMER_STOP(SymbolicFactorization);
}


template <typename T> 
void SupernodalMatrix<T>::refineSupernodes(int ordflag,int altflag){
  ETree& tree = ETree_;
  Ordering & aOrder = Order_;
  SYMPACK::vector<Int> & supMembership = SupMembership_; 
  SYMPACK::vector<Int> & xsuper = Xsuper_; 

  SYMPACK::vector<int> ixlindx;
  SYMPACK::vector<int> ilindx;

  //Gather locXlindx_ and locLindx_
  {
    SYMPACK::vector<Ptr> xlindx;
    SYMPACK::vector<Idx> lindx;
    //get other proc vertex counts
    Idx localVertexCnt = locXlindx_.size()-1;
    SYMPACK::vector<Idx> remoteVertexCnt(np,0);
    MPI_Gather(&localVertexCnt,sizeof(localVertexCnt),MPI_BYTE,&remoteVertexCnt[0],sizeof(localVertexCnt),MPI_BYTE,0,CommEnv_->MPI_GetComm());
    Idx totalVertexCnt = std::accumulate(remoteVertexCnt.begin(),remoteVertexCnt.end(),0,std::plus<Idx>());
    if(iam==0){
      xlindx.resize(totalVertexCnt+1);
    }
    //compute receive displacements
    SYMPACK::vector<int> rsizes(np,0);
    for(int p = 0; p<np;p++){rsizes[p] = (int)remoteVertexCnt[p]*sizeof(Ptr);}
    SYMPACK::vector<int> rdispls(np+1,0);
    rdispls[0]=0;
    std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
    MPI_Gatherv(&locXlindx_[0],localVertexCnt*sizeof(Ptr),MPI_BYTE,&xlindx[0],&rsizes[0],&rdispls[0],MPI_BYTE,0,CommEnv_->MPI_GetComm());


    Ptr localEdgeCnt = locLindx_.size();
    SYMPACK::vector<Ptr> remoteEdgeCnt(np,0);
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

    //logfileptr->OFS()<<"locXlindx: "<<locXlindx_<<endl;
    //logfileptr->OFS()<<"locLindx: "<<locLindx_<<endl;
  }





  if(iam==0){

    //logfileptr->OFS()<<"xlindx: "<<ixlindx<<endl;
    //logfileptr->OFS()<<"lindx: "<<ilindx<<endl;

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


    FORTRAN(ordsup) (
        &ordflag, &altflag, &neqns , &nofsub, &nsuper, 
        &xsuper[0], &ixlindx[0], &ilindx[0], &supMembership[0], &aOrder.perm[0],
        &aOrder.invp[0], 
        &freeforw[0], &freeback[0], &sforw[0], &sback[0], 
        &setseg_forw[0], &setseg_back[0], &nodehead[0], 
        &nodeforw[0], &nodeback[0], 
        &setsnode[0], &supperm[0], &mark[0], &set[0], &compset[0],
        &invp2[0], &heap[0]);

  }
  //broadcast perm and invp

  // broadcast invp
  Int N = aOrder.invp.size();
  MPI_Bcast(&aOrder.invp[0],N*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
  MPI_Bcast(&aOrder.perm[0],N*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
}


template <typename T> 
template <class Allocator>
SuperNode<T,Allocator> * SupernodalMatrix<T>::CreateSuperNode(DecompositionType type,Int aiId, Int aiFc, Int aiLc, Int ai_num_rows, Int aiN, Int aiNZBlkCnt){
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
 SuperNode<T,Allocator> * SupernodalMatrix<T>::CreateSuperNode(DecompositionType type,char * dataPtr,size_t size, Idx firstRow){
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





#include "SupernodalMatrix_impl_FB_pull.hpp"


#ifdef _INDEFINITE_
#undef _INDEFINITE_
#endif

}


#endif 
