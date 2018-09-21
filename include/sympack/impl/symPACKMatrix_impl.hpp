/*
   Copyright (c) 2016 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.  

Author: Mathias Jacquelin

This file is part of symPACK. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to Lawrence Berkeley National
Laboratory, without imposing a separate written license agreement for such
Enhancements, then you hereby grant the following license: a non-exclusive,
royalty-free perpetual license to install, use, modify, prepare derivative
works, incorporate into other computer software, distribute, and sublicense
such enhancements or derivative works thereof, in binary and source code form.
*/

#ifndef _SYMPACK_MATRIX_IMPL_HPP_
#define _SYMPACK_MATRIX_IMPL_HPP_

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
#include <type_traits>

namespace symPACK{
  extern "C" {
    void FORTRAN(ordsup) (int * ordflag, int *  altflag, int *  NEQNS, int *  nofsub, int *  nsuper, 
        int * xsuper, int *  xlindx, int *  lindx , int *  snode , int *  perm  , 
        int * invp  , int *  freeforw, int *  freeback, int *  sforw, int *  sback, 
        int * setseg_forw, int *  setseg_back, int *  nodehead, 
        int * nodeforw, int *  nodeback, 
        int *  setsnode, int *  supperm, int *  mark, int *  set  , int *  compset,
        int *  invp2 , int *  heap                             );

    void FORTRAN(ordsup_ind_tsp_paths2)
      (  int * nadj  , int * neqns , int * nofsub, int * nsuper, int * supsiz,
         int * xsuper, int * xlindx, int * lindx , int * snode , int * xadj  , 
         int * adjncy, int * etpar , int * perm  , int * invp  , int * iflag , 
         int * xskadj, int * sklenf, int * sklenb, int * skadj , int * invp2 , 
         int * link  , int * fstloc, int * sperm , int * fstloc2, int * dist1, 
         int * suppar, int * iwsiz , int * iwork , int * rep             );

    void FORTRAN(ordsup_ind_tsp_paths)
      (  int * nadj  , int * neqns , int * nofsub, int * nsuper, int * supsiz,
         int * xsuper, int * xlindx, int * lindx , int * snode , int * xadj  , 
         int * adjncy, int * etpar , int * perm  , int * invp  , int * iflag , 
         int * xskadj, int * sklenf, int * sklenb, int * skadj , int * invp2 , 
         int * link  , int * fstloc, int * sperm , int * fstloc2, int * dist1, 
         int * suppar, int * iwsiz , int * iwork              );


  }
}

#define NEW_SOLVE


#define SPLIT_AT_BOUNDARY

namespace symPACK{


  template<typename T> 
  inline void symPACKMatrix<T>::generateTaskGraph(supernodalTaskGraph<FBTask> & taskGraph,
      std::vector<Int> & AggregatesToRecv,  std::vector<Int>& LocalAggregates)
  {
    //we will need to communicate if only partial xlindx_, lindx_
    //idea: build tasklist per processor and then exchange
    std::map<Idx, std::list<std::pair<Idx,Idx> >  > Updates;
    std::vector<int> marker(this->np,0);
    Int numLocSnode = this->XsuperDist_[this->iam+1]-this->XsuperDist_[this->iam];
    Int firstSnode = this->XsuperDist_[this->iam];
    for(Int locsupno = 1; locsupno<this->locXlindx_.size(); ++locsupno){
      Idx I = locsupno + firstSnode-1;
      Int iOwner = Mapping_->Map(I-1,I-1);

      Int first_col = this->Xsuper_[I-1];
      Int last_col = this->Xsuper_[I]-1;

      //Create the factor task on the owner
      Updates[iOwner].push_back(std::make_pair(I,I));

      Int J = -1; 
      Ptr lfi = this->locXlindx_[locsupno-1];
      Ptr lli = this->locXlindx_[locsupno]-1;
      Idx prevSnode = -1;
      for(Ptr sidx = lfi; sidx<=lli;sidx++){
        Idx row = this->locLindx_[sidx-1];
        J = this->SupMembership_[row-1];
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
    vector<int> ssizes(this->np,0);
    for(auto itp = Updates.begin();itp!=Updates.end();itp++){
      ssizes[itp->first] = itp->second.size();//*sizeof(std::pair<Idx,Idx>);
    }

    //compute send displacements
    vector<int> sdispls(this->np+1,0);
    sdispls[0] = 0;
    std::partial_sum(ssizes.begin(),ssizes.end(),&sdispls[1]);

    //Build the contiguous array of pairs
    vector<std::pair<Idx,Idx> > sendbuf;
    sendbuf.reserve(sdispls.back());///sizeof(std::pair<Idx,Idx>));

    for(auto itp = Updates.begin();itp!=Updates.end();itp++){
      sendbuf.insert(sendbuf.end(),itp->second.begin(),itp->second.end());
    }

    //gather receive sizes
    vector<int> rsizes(this->np,0);
    MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,CommEnv_->MPI_GetComm());

    //compute receive displacements
    vector<int> rdispls(this->np+1,0);
    rdispls[0] = 0;
    std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);


    //Now do the alltoallv
    vector<std::pair<Idx,Idx> > recvbuf(rdispls.back());///sizeof(std::pair<Idx,Idx>));
    MPI_Alltoallv(&sendbuf[0],&ssizes[0],&sdispls[0],type,&recvbuf[0],&rsizes[0],&rdispls[0],type,CommEnv_->MPI_GetComm());

    MPI_Type_free(&type);

    //now process recvbuf and add the tasks to my tasklists
    taskGraph.taskLists_.resize(this->TotalSupernodeCnt(),NULL);
    //Int localTaskCount = 0;
    for(auto it = recvbuf.begin();it!=recvbuf.end();it++){
      FBTask curUpdate;
      curUpdate.src_snode_id = it->first;
      curUpdate.tgt_snode_id = it->second;
      curUpdate.type=(it->first==it->second)?FACTOR:UPDATE;

      Int J = curUpdate.tgt_snode_id;

      taskGraph.addTask(curUpdate);

      ////create the list if needed
      //if(taskGraph.taskLists_[J-1] == NULL){
      //  taskGraph.taskLists_[J-1]=new std::list<FBTask>();
      //}
      //taskGraph.taskLists_[J-1]->push_back(curUpdate);

      //TODO
      //localTaskCount++;
    }
    //taskGraph.localTaskCount_=localTaskCount;



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
            if(this->iam==iOwner){
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
    levels.resize(this->Xsuper_.size());
    levels[0]=0;
    Int numLevel = 0; 
    for(Int i=this->Xsuper_.size()-1-1; i>=0; i-- ){     
      Int fcol = this->Xsuper_[i];
      bassert(fcol-1>=0);
      Int pcol =this->ETree_.PostParent(fcol-1);
      Int supno = pcol>0?this->SupMembership_[pcol-1]:0;
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


#if 1
  template<typename T> 
    inline void symPACKMatrix<T>::generateTaskGraph_New(taskGraph & graph,
      std::vector<Int> & AggregatesToRecv,  std::vector<Int>& LocalAggregates, std::vector<Int> & mw, std::vector<Int> & mh)
    {


      //we will need to communicate if only partial xlindx_, lindx_
      //idea: build tasklist per processor and then exchange
      //tuple is: src_snode,tgt_snode,src_first_row,op_type,
      using upd_tuple_t = std::tuple<Idx,Idx,Idx,Factorization::op_type>;
      //using upd_:std::tuple<Idx,Idx,Idx,Factorization::op_type> ; 
      std::map<Idx, std::list< upd_tuple_t> > Updates;
      std::vector<int> marker(this->np,0);

      std::vector< std::map< Idx, std::pair<Idx,std::set<Idx> > > > updCnt(this->Xsuper_.size() -1 );
      
        Int numLocSnode = this->XsuperDist_[this->iam+1]-this->XsuperDist_[this->iam];
        Int firstSnode = this->XsuperDist_[this->iam];
        for(Int locsupno = 1; locsupno<this->locXlindx_.size(); ++locsupno){
          Idx I = locsupno + firstSnode-1;
          Int first_col = this->Xsuper_[I-1];
          Int last_col = this->Xsuper_[I]-1;
          Int iOwner = Mapping_->Map(I-1,I-1);

          //Create the factor task on the owner
          Updates[iOwner].push_back(std::make_tuple(I,I,first_col,Factorization::op_type::FACTOR));

          Int J = -1; 
          Ptr lfi = this->locXlindx_[locsupno-1];
          Ptr lli = this->locXlindx_[locsupno]-1;

          Idx blockIdx = -1;
          Idx prevRow = -1;
          Idx prevSnode = -1;

          for(Ptr sidx = lfi; sidx<=lli;sidx++){
            Idx row = this->locLindx_[sidx-1];
            J = this->SupMembership_[row-1];

            //Split at boundary or after diagonal block
            if(J!=prevSnode){

              Int iUpdater = Mapping_->Map(J-1,I-1);

              if(J>I){


            //if(J==26){gdb_lock();}

                if(marker[iUpdater]!=I){
                  //create the task on iUpdater
                  Updates[iUpdater].push_back(std::make_tuple(I,J,row,Factorization::op_type::UPDATE));
                }

#ifdef _SEQ_SPECIAL_CASE_
                if(Multithreading::NumThread==1){
                  //create only one task if I don't own the factor
                  if(iUpdater!=iOwner){
                    marker[iUpdater]=I;
                  }
                }
#endif
              }
            }
            prevSnode = J;
          }
        }

      MPI_Datatype type;
      MPI_Type_contiguous( sizeof(upd_tuple_t), MPI_BYTE, &type );
      MPI_Type_commit(&type);

      //then do an alltoallv
      //compute send sizes
      vector<int> ssizes(this->np,0);
      for(auto itp = Updates.begin();itp!=Updates.end();itp++){
        ssizes[itp->first] = itp->second.size();//*sizeof(std::pair<Idx,Idx>);
      }

      //compute send displacements
      vector<int> sdispls(this->np+1,0);
      sdispls[0] = 0;
      std::partial_sum(ssizes.begin(),ssizes.end(),&sdispls[1]);

      //Build the contiguous array of pairs
      vector<upd_tuple_t> sendbuf;
      sendbuf.reserve(sdispls.back());

      for(auto itp = Updates.begin();itp!=Updates.end();itp++){
        sendbuf.insert(sendbuf.end(),itp->second.begin(),itp->second.end());
      }

      //gather receive sizes
      vector<int> rsizes(this->np,0);
      MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,CommEnv_->MPI_GetComm());

      //compute receive displacements
      vector<int> rdispls(this->np+1,0);
      rdispls[0] = 0;
      std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);


      //Now do the alltoallv
      vector<upd_tuple_t> recvbuf(rdispls.back());///sizeof(std::pair<Idx,Idx>));
      MPI_Alltoallv(&sendbuf[0],&ssizes[0],&sdispls[0],type,&recvbuf[0],&rsizes[0],&rdispls[0],type,CommEnv_->MPI_GetComm());

      MPI_Type_free(&type);

      std::hash<std::string> hash_fn;


      //now process recvbuf and add the tasks to the taskGraph
      size_t max_mem_req = 0;
      for(auto it = recvbuf.begin();it!=recvbuf.end();it++){
        std::shared_ptr<GenericTask> pTask(new SparseTask);
        SparseTask & Task = *(SparseTask*)pTask.get();

        Task.meta.resize(3*sizeof(Int)+sizeof(Factorization::op_type));
        Int * meta = reinterpret_cast<Int*>(Task.meta.data());
        meta[0] = std::get<0>(*it);
        meta[1] = std::get<1>(*it);
        meta[2] = std::get<2>(*it);
        Factorization::op_type & type = *reinterpret_cast<Factorization::op_type*>(&meta[3]);
        type = std::get<3>(*it);//(it->first==it->second)?Factorization::op_type::FACTOR:Factorization::op_type::UPDATE;

        Int src = meta[0];
        Int tgt = meta[1];
        size_t mem_req = 0;
        switch(type){
          case Factorization::op_type::FACTOR:
            {
              Task.remote_deps = AggregatesToRecv[tgt-1];
              Task.local_deps = LocalAggregates[tgt-1];
              //memory req is: aggregate to recv * mw[tgt] * mh[tgt]
              mem_req = Task.remote_deps * mw[tgt] * mh[tgt] * sizeof(T);

              //extra work storage for potrf
              if(this->options_.decomposition == DecompositionType::LDL){
                Int NB = lapack::Ilaenv( 1, "DPOTRF", "U", mw[src], -1, -1, -1 );
                mem_req += NB * mw[src] * sizeof(T); 
              }
            }
            break;
          case Factorization::op_type::UPDATE:
            {
                Int iOwner = Mapping_->Map(src-1,src-1);
                //If I own the factor, it is a local dependency
                if(this->iam==iOwner){
                  Task.remote_deps = 0;
                  Task.local_deps = 1;
                }
                else{
                  Task.remote_deps = 1;
                  Task.local_deps = 0;                    
                }

                Int iOwnerTgt = Mapping_->Map(tgt-1,tgt-1);
                if(this->iam!=iOwnerTgt){
                  //memory req is: mw[tgt] * mh[tgt]
                  mem_req += mw[tgt] * mh[tgt] * sizeof(T);
                }
                if(this->iam!=iOwner){
                  //memory req is: mw[src]*mh[src]
                  mem_req += mw[src] * mh[src] * sizeof(T);
                }

                if(this->options_.decomposition == DecompositionType::LL){
                  mem_req += mh[src]*mw[tgt]*sizeof(T) + mw[tgt]*sizeof(Idx) + mh[src]*sizeof(Ptr);
                }
                else{
                  mem_req += mh[src]*mw[tgt]*sizeof(T) + mw[tgt]*sizeof(Idx) + mh[src]*sizeof(Ptr);
                  mem_req += mw[src]*mw[tgt]*sizeof(T);
                }
            }
            break;
        }


        max_mem_req = std::max(mem_req,max_mem_req);



        Task.getHash = [&](char * ameta=nullptr)->GenericTask::id_type{
          char * pmeta = ameta;
          //if(pmeta==nullptr){ pmeta = Task.meta.data(); } 
          std::stringstream sstr;
          sstr<<meta[0]<<"_"<<meta[1]<<"_"<<0<<"_"<<(Int)(*reinterpret_cast<Factorization::op_type*>(&meta[3]));
          return hash_fn(sstr.str());
        };
        Task.id = Task.getHash(Task.meta.data());

        graph.addTask(pTask);


      }
      logfileptr->OFS()<<"Maximum single task memory requirement is: "<<max_mem_req<<std::endl;
    }
#endif









  template <typename T> 
   inline  void symPACKMatrix<T>::Factorize(){
    SYMPACK_TIMER_START(NUMERICAL_FACT);
    if(this->iam<this->np){
      switch(this->options_.factorization){
        case FANBOTH:
          FanBoth_New();
          break;
        case FANOUT:
          FanBoth();
          break;
        default:
          FanBoth();
          break;
      }
      SYMPACK_TIMER_STOP(NUMERICAL_FACT);

#ifdef _SYMPACK_PROFILE_COMM_
      //TODO this needs debugging
      if(this->iam<this->np){
        logfileptr->OFS()<<"Local volume of communication: "<<gVolComm<<std::endl;
        logfileptr->OFS()<<"Local number of messages: "<<gNumMsg<<std::endl;

        size_t totalVolComm = 0;
        MPI_Reduce(&gVolComm,&totalVolComm,1,MPI_UINT64_T,MPI_SUM,0,CommEnv_->MPI_GetComm());
        size_t totalNumMsg = 0;
        MPI_Reduce(&gNumMsg,&totalNumMsg,1,MPI_UINT64_T,MPI_SUM,0,CommEnv_->MPI_GetComm());

        if(this->iam==0){
          std::cout<<"Total volume of communication: "<<totalVolComm<<std::endl;
          std::cout<<"Total number of messages: "<<totalNumMsg<<std::endl;
        }

        gVolComm=0;
        gNumMsg=0;
      }
#endif
    }
  }


  template <typename T> inline void symPACKMatrix<T>::Solve(T * RHS, int nrhs,  T * Xptr) {
    scope_timer(a,SPARSE_SOLVE);

    if (this->options_.iterRefinement){
      abort();  
      this->solveNew_(RHS,nrhs,Xptr);
      //do{
      //  this->solve_(RHS,nrhs,Xptr);

      //  //update

      //  //compute residual
      //}
      //while();
    }
    else{

      for(auto ptr: Contributions_){ delete ptr; }
      Contributions_.clear();
      Contributions2_.clear();

#ifndef NEW_SOLVE
      this->solve_(RHS,nrhs,Xptr);
#else
//      if(this->iam<this->np){
        switch(this->options_.factorization){
          case FANBOTH:
            this->solveNew2_(RHS,nrhs,Xptr);
            break;
          case FANOUT:
            this->solveNew2_(RHS,nrhs,Xptr);
            break;
          default:
            this->solveNew2_(RHS,nrhs,Xptr);
            break;
        }
//      }

      //      this->solveNew_(RHS,nrhs,Xptr);
#endif
    }

  }

  //Solve related routines
  template <typename T> inline void symPACKMatrix<T>::solve_(T * RHS, int nrhs,  T * Xptr) {
    scope_timer(a,SPARSE_SOLVE_INTERNAL);

    Int n = this->iSize_;
    //Int this->iam = CommEnv_->MPI_Rank();
    //Int this->np  = CommEnv_->MPI_Size();

    if(this->iam<this->np){
      std::vector<Int> children(this->Xsuper_.size());
      SetValue(children,0);

      for(Int I=1;I<this->Xsuper_.size()-1;I++){
        Int fc = this->Xsuper_[I-1];
        Int lc = this->Xsuper_[I]-1;

        Int parent = this->ETree_.PostParent(lc-1);
        if(parent!=0){
          ++children[this->SupMembership_[parent-1]-1];
        }
      }

#ifdef _DEBUG_
      logfileptr->OFS()<<"Children std::vector is"<<children<<std::endl;
#endif


      std::vector<Int> UpdatesToDo = children;

      Contributions_.resize(LocalSupernodes_.size());
      std::vector<std::stack<Int> > LocalUpdates(LocalSupernodes_.size());

      AsyncComms outgoingSend;

      //This corresponds to the k loop in dtrsm
      for(Int I=1;I<this->Xsuper_.size();I++){
        Int iOwner = this->Mapping_->Map(I-1,I-1);
        //If I own the column, factor it
        if( iOwner == this->iam ){
          //Create the blocks of my contrib with the same nz structure as L
          //and the same width as the final solution
          //MEMORY CONSUMPTION TOO HIGH ?
          Int iLocalI = snodeLocalIndex(I);
          SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];
          SuperNode<T,MallocAllocator> * contrib = CreateSuperNode<MallocAllocator>(this->options_.decomposition,I,cur_snode->FirstRow(),1,nrhs, cur_snode->NRowsBelowBlock(0) ,this->iSize_,this->options_.panel);

          Contributions_[iLocalI-1] = contrib;


          for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){
            NZBlockDesc & cur_desc = cur_snode->GetNZBlockDesc(blkidx);
            contrib->AddNZBlock(cur_snode->NRows(blkidx), cur_desc.GIndex);
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
          Int parent = this->ETree_.PostParent(cur_snode->LastCol()-1);
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

            Int iOwner = this->iam;
            contrib->forward_update(dist_contrib,iOwner,this->iam);
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
            SuperNode<T,MallocAllocator> * dist_contrib = CreateSuperNode<MallocAllocator>(this->options_.decomposition,&src_blocks[0],bytes_received);
            //TODO Put this back
            //Deserialize(&src_blocks[0],dist_contrib);
#ifdef _DEBUG_
            logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib.Id()<<std::endl;
#endif

            Int iOwner = recv_status.MPI_SOURCE;
            contrib->forward_update(dist_contrib,iOwner,this->iam);

            --UpdatesToDo[I-1];
            delete dist_contrib;
          }

          assert(UpdatesToDo[I-1]==0);

          if(UpdatesToDo[I-1]==0){

            //now compute MY contribution
            contrib->forward_update_contrib(RHS,cur_snode,this->Order_.perm);

            //send to my parent
            if(parent!=0){
              Int parent_snode_id = this->SupMembership_[parent-1];

              Int iTarget = this->Mapping_->Map(parent_snode_id-1,parent_snode_id-1);

              if(iTarget!=this->iam){
#ifdef _DEBUG_
                logfileptr->OFS()<<"Remote Supernode "<<parent_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif

                Int tgt_first_col = this->Xsuper_[parent_snode_id-1];
                Int tgt_last_col = this->Xsuper_[parent_snode_id]-1;
                Int src_nzblk_idx = 1;
                NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);

                Int src_first_row = pivot_desc.GIndex;

                bool isSkipped= false;

                Int next_local_contrib = (iLocalI < Contributions_.size())?Contributions_[iLocalI]->Id():this->Xsuper_.size();
                if(next_local_contrib< parent_snode_id || true){
                  //need to push the prev src_last_row
                  ContribsToSend.push(DelayedComm(contrib->Id(),parent_snode_id,1,src_first_row));
#ifdef _DEBUG_DELAY_
                  std::cout<<"P"<<this->iam<<" has delayed update from Contrib "<<I<<" to "<<parent_snode_id<<" from row "<<src_first_row<<std::endl;
#endif
                  isSkipped= true;
                }

                if(!isSkipped){
                  //Create a new Icomm buffer, serialize the contribution
                  // in it and add it to the outgoing comm list

                  Icomm * send_buffer = new Icomm();
                  contrib->Serialize(*send_buffer,src_nzblk_idx,src_first_row);
                  AddOutgoingComm(outgoingSend,send_buffer);

                  if( outgoingSend.size() > this->options_.maxIsend){
                    MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,parent_snode_id,CommEnv_->MPI_GetComm());
                    //that deletes the Icomm
                    outgoingSend.pop_back();
                  }
                  else{
                    MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,parent_snode_id,CommEnv_->MPI_GetComm(),&outgoingSend.back()->Request);
                  }
#ifdef _DEBUG_            
                  logfileptr->OFS()<<"     Send contribution "<<I<<" to Supernode "<<parent_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
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

          Int parent = this->ETree_.PostParent(cur_snode->LastCol()-1);

          std::vector<Int> isBlockUpdated(contrib->NZBlockCnt(),0);

          if(parent!=0){
            Int parent_snode_id = this->SupMembership_[parent-1];
            Int iTarget = this->Mapping_->Map(parent_snode_id-1,parent_snode_id-1);
            //Do all my updates (Local and remote)
            //Local updates
            SuperNode<T,MallocAllocator> * dist_contrib;
            if(!LocalUpdates[iLocalI-1].empty()){
              Int contrib_snode_id = LocalUpdates[iLocalI-1].top();
              LocalUpdates[iLocalI-1].pop();
              dist_contrib = snodeLocal(contrib_snode_id,Contributions_);
              contrib->back_update(dist_contrib);
            }
            else{
              //Receive parent contrib
              MPI_Status recv_status;
              Int bytes_received = 0;
              MPI_Probe(iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);
              MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
              src_blocks.resize(bytes_received);

              MPI_Recv(&src_blocks[0],bytes_received,MPI_BYTE,iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);

              dist_contrib = CreateSuperNode<MallocAllocator>(this->options_.decomposition,&src_blocks[0],bytes_received);
              //TODO Replace this
              //Deserialize(&src_blocks[0],*dist_contrib); 
              dist_contrib->InitIdxToBlk();


#ifdef _DEBUG_
              logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib->Id()<<std::endl;
#endif
              contrib->back_update(dist_contrib);
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
              Int child_snode_id = this->SupMembership_[colIdx-1];

              Int parent = this->ETree_.PostParent(colIdx-1);
              if(parent!=0){
                if(this->SupMembership_[parent-1]==cur_snode->Id()){
                  Int iTarget = this->Mapping_->Map(child_snode_id-1,child_snode_id-1);

                  if(iTarget!=this->iam){

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
                      std::cout<<"P"<<this->iam<<" has delayed update from Contrib "<<I<<" to "<<child_snode_id<<" from row "<<src_first_row<<std::endl;
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


                      if( outgoingSend.size() > this->options_.maxIsend){
                        MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm());
                        outgoingSend.pop_back();
                      }
                      else{
                        MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,child_snode_id,CommEnv_->MPI_GetComm(),&outgoingSend.back()->Request);
                      }

#ifdef _DEBUG_            
                      logfileptr->OFS()<<"     Send contribution "<<I<<" to Supernode "<<child_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
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
              colIdx = this->Xsuper_[child_snode_id-1]-1;
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

  }

  template<typename T> inline void symPACKMatrix<T>::GetSolution(T * B, int nrhs){
    Int n = this->iSize_;
    //Int this->iam = CommEnv_->MPI_Rank();
    //Int this->np  = CommEnv_->MPI_Size();

    //if(this->iam<this->np)
    {
      std::fill(B,B+n*nrhs,T(0.0));

      //    Int nrhs = B.n();
      //Gather B from everybody and put it in the original matrix order
      std::vector<T> tmp_nzval;
      for(Int I=1; I<this->Xsuper_.size();++I){
        Int iOwner = this->Mapping_->Map(I-1,I-1);
        T * data;
        Int snode_size = this->Xsuper_[I] - this->Xsuper_[I-1];
        Int nzcnt = snode_size * nrhs;

        if( iOwner == this->iam ){
#ifndef NEW_SOLVE
          SuperNode<T,MallocAllocator> * contrib = snodeLocal(I,Contributions_);
#else
          Int Ilocal = snodeLocalIndex(I); 
          auto contrib = std::dynamic_pointer_cast< SuperNode<T,UpcxxAllocator> >(Contributions2_[Ilocal-1]);
#endif
          //logfileptr->OFS()<<*contrib<<std::endl;
          data = contrib->GetNZval(0);
        }
        else{
          tmp_nzval.resize(nzcnt);
          data = &tmp_nzval[0];
        }

        //MPI_Bcast(data,nzcnt*sizeof(T),MPI_BYTE,iOwner,CommEnv_->MPI_GetComm());

        if(this->iam==iOwner){
          for(Int i = 0; i<snode_size; ++i){ 
            Int destRow = this->Xsuper_[I-1] + i;
            destRow = this->Order_.perm[destRow - 1];
            auto * Bptr = &B[destRow-1];
            auto * dataPtr = &data[i*nrhs];
            for(Int j = 0; j<nrhs; ++j){
              Bptr[j*n] = dataPtr[j];
            }
          }
        }
      }

      mpi::Allreduce((T*)MPI_IN_PLACE,&B[0],n*nrhs,MPI_SUM,this->fullcomm_);
      MPI_Barrier(CommEnv_->MPI_GetComm());
    }
  }


  template <typename T>
    template< class Alloc>
    inline void symPACKMatrix<T>::SendDelayedMessagesUp(Int iLocalI, CommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode<T,Alloc> *> & snodeColl){
      if(snodeColl.empty() || MsgToSend.empty()) { return;}

      //Index of the last global snode to do
      Int last_snode_id = this->Xsuper_.size()-1;
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

        if(tgt_snode_id < next_snode_id || is_last ){
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
    inline void symPACKMatrix<T>::SendDelayedMessagesDown(Int iLocalI, DownCommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode<T,Alloc> *> & snodeColl){
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

        if(tgt_snode_id>next_snode_id || is_last ){

          SuperNode<T,Alloc> * prev_src_snode = snodeLocal(src_snode_id,snodeColl);
          //this can be sent now

          Int iTarget = this->Mapping_->Map(tgt_snode_id-1,tgt_snode_id-1);
          if(iTarget != this->iam){
#ifdef _DEBUG_DELAY_
            logfileptr->OFS()<<"P"<<this->iam<<" has sent update from Supernode "<<prev_src_snode->Id()<<" to Supernode "<<tgt_snode_id<<std::endl;
            std::cout<<"P"<<this->iam<<" has sent update from Supernode "<<prev_src_snode->Id()<<" to Supernode "<<tgt_snode_id<<std::endl;
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


            if( OutgoingSend.size() > this->options_.maxIsend){
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
    inline void symPACKMatrix<T>::SendMessage(const DelayedComm & comm, AsyncComms & OutgoingSend, std::vector<SuperNode<T,Alloc> *> & snodeColl){
      Int src_snode_id = comm.src_snode_id;
      Int tgt_snode_id = comm.tgt_snode_id;
      Int src_nzblk_idx = comm.src_nzblk_idx;
      Int src_first_row = comm.src_first_row;


      SuperNode<T,Alloc> * prev_src_snode = snodeLocal(src_snode_id,snodeColl);

      //this can be sent now
      Int iTarget = this->Mapping_->Map(tgt_snode_id-1,tgt_snode_id-1);
      if(iTarget != this->iam){
#ifdef _DEBUG_DELAY_
        logfileptr->OFS()<<"P"<<this->iam<<" has sent update from Supernode "<<prev_src_snode->Id()<<" to Supernode "<<tgt_snode_id<<std::endl;
        std::cout<<"P"<<this->iam<<" has sent update from Supernode "<<prev_src_snode->Id()<<" to Supernode "<<tgt_snode_id<<std::endl;
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


        if( OutgoingSend.size() > this->options_.maxIsend){
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




  template <typename T> inline void symPACKMatrix<T>::AdvanceOutgoing(AsyncComms & outgoingSend){
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


  template<typename T> inline void symPACKMatrix<T>::AddOutgoingComm(AsyncComms & outgoingSend, Icomm * send_buffer){
    outgoingSend.push_back(send_buffer);
  }


  template <typename T> inline void symPACKMatrix<T>::GetUpdatingSupernodeCount(std::vector<Int> & sc,std::vector<Int> & mw, std::vector<Int> & mh, std::vector<Int> & numBlk){
    sc.resize(this->Xsuper_.size(),I_ZERO);
    std::vector<Int> marker(this->Xsuper_.size(),I_ZERO);
    mw.resize(this->Xsuper_.size(),I_ZERO);
    mh.resize(this->Xsuper_.size(),I_ZERO);
    numBlk.resize(this->Xsuper_.size(),I_ZERO);

    //Int numLocSnode = ( (this->Xsuper_.size()-1) / this->np);
    //Int firstSnode = this->iam*numLocSnode + 1;

    Int numLocSnode = this->XsuperDist_[this->iam+1]-this->XsuperDist_[this->iam];
    Int firstSnode = this->XsuperDist_[this->iam];

    for(Int locsupno = 1; locsupno<this->locXlindx_.size(); ++locsupno){
      Idx s = locsupno + firstSnode-1;

      Int first_col = this->Xsuper_[s-1];
      Int last_col = this->Xsuper_[s]-1;

      Ptr lfi = this->locXlindx_[locsupno-1];
      Ptr lli = this->locXlindx_[locsupno]-1;
      mh[s-1] = lli-lfi+1;



      //count number of contiguous blocks (at least one diagonal block)
      Idx iPrevRow = this->locLindx_[lfi-1]-1;
      Idx iFirstRow = this->locLindx_[lfi-1];

      Idx width = mh[s-1]; 

#ifdef SPLIT_AT_BOUNDARY
      Int nextSup = s+1;
      Idx next_fc = first_col;
      Idx next_lc = last_col;
#endif

      Idx nzBlockCnt = 1;
      Idx prevSnode = -1;
      for(Ptr sidx = lfi; sidx<=lli;sidx++){

        Idx row = this->locLindx_[sidx-1];


        //enforce the first block to be a square diagonal block
        if(nzBlockCnt==1 && row>last_col){
          nzBlockCnt++;
        }
        else{
#ifdef SPLIT_AT_BOUNDARY
          if( nextSup<=this->Xsuper_.size()-1){
            next_fc = this->Xsuper_[nextSup-1];
            next_lc = this->Xsuper_[nextSup]-1;
          }
          while(row>=next_lc){
            nextSup++;
            if( nextSup<=this->Xsuper_.size()-1){
              next_fc = this->Xsuper_[nextSup-1];
              next_lc = this->Xsuper_[nextSup]-1;
              //             logfileptr->OFS()<<"P Moving on to next supernode"<<std::endl;
            }
            else{
              break;
            }
          } 

          if(row==iPrevRow+1){
            if( nextSup<this->Xsuper_.size()-0 && row==next_fc){
              //              logfileptr->OFS()<<"P "<<row<<" Next fc and lc are: "<<next_fc<<" "<<next_lc<<std::endl;
              //              logfileptr->OFS()<<"Crossed boundary"<<std::endl;
              nzBlockCnt++;
            }
          }
#endif


          if(row!=iPrevRow+1){
            nzBlockCnt++;
          }
        }
        iPrevRow=row;

        Int supno = this->SupMembership_[row-1];

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
    MPI_Allreduce(MPI_IN_PLACE,&sc[0],sc.size(),MPI_INT,MPI_SUM,this->fullcomm_);
    MPI_Allreduce(MPI_IN_PLACE,&mw[0],mw.size(),MPI_INT,MPI_MAX,this->fullcomm_);
    MPI_Allreduce(MPI_IN_PLACE,&mh[0],mh.size(),MPI_INT,MPI_MAX,this->fullcomm_);
    MPI_Allreduce(MPI_IN_PLACE,&numBlk[0],numBlk.size(),MPI_INT,MPI_SUM,this->fullcomm_);
  }

  template<typename T>
    inline void symPACKMatrix<T>::DumpMatlab(){
      logfileptr->OFS()<<"+sparse([";
      for(Int I=1;I<this->Xsuper_.size();I++){
        Int src_first_col = this->Xsuper_[I-1];
        Int src_last_col = this->Xsuper_[I]-1;
        Int iOwner = this->Mapping_->Map(I-1,I-1);
        if( iOwner == this->iam ){
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
      for(Int I=1;I<this->Xsuper_.size();I++){
        Int src_first_col = this->Xsuper_[I-1];
        Int src_last_col = this->Xsuper_[I]-1;
        Int iOwner = this->Mapping_->Map(I-1,I-1);
        if( iOwner == this->iam ){
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
      for(Int I=1;I<this->Xsuper_.size();I++){
        Int src_first_col = this->Xsuper_[I-1];
        Int src_last_col = this->Xsuper_[I]-1;
        Int iOwner = this->Mapping_->Map(I-1,I-1);
        if( iOwner == this->iam ){
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
      logfileptr->OFS()<<"],"<<this->iSize_<<","<<this->iSize_<<")"<<std::endl;



    }

  template<typename T>
    inline void symPACKMatrix<T>::Dump(){
      for(Int I=1;I<this->Xsuper_.size();I++){
        Int src_first_col = this->Xsuper_[I-1];
        Int src_last_col = this->Xsuper_[I]-1;
        Int iOwner = this->Mapping_->Map(I-1,I-1);
        if( iOwner == this->iam ){
          SuperNode<T> * src_snode = snodeLocal(I);


          logfileptr->OFS()<<"+++++++++++++"<<I<<"++++++++"<<std::endl;

          logfileptr->OFS()<<"cols: ";
          for(Int i = 0; i< src_snode->Size(); ++i){
            logfileptr->OFS()<<" "<<this->Order_.perm[src_first_col+i-1];
          }
          logfileptr->OFS()<<std::endl;
          logfileptr->OFS()<<std::endl;
          for(int blkidx=0;blkidx<src_snode->NZBlockCnt();++blkidx){

            NZBlockDesc & desc = src_snode->GetNZBlockDesc(blkidx);
            T * val = src_snode->GetNZval(desc.Offset);
            Int nRows = src_snode->NRows(blkidx);

            Int row = desc.GIndex;
            for(Int i = 0; i< nRows; ++i){
              logfileptr->OFS()<<row+i<<" | "<<this->Order_.perm[row+i-1]<<":   ";
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
    inline void symPACKMatrix<T>::DumpContrib(){
      for(Int I=1;I<this->Xsuper_.size();I++){
        Int src_first_col = this->Xsuper_[I-1];
        Int src_last_col = this->Xsuper_[I]-1;
        Int iOwner = this->Mapping_->Map(I-1,I-1);
        if( iOwner == this->iam ){
          SuperNode<T,MallocAllocator> * src_snode = snodeLocal(I,Contributions_);


          logfileptr->OFS()<<"+++++++++++++"<<I<<"++++++++"<<std::endl;

          logfileptr->OFS()<<"cols: ";
          for(Int i = 0; i< src_snode->Size(); ++i){
            logfileptr->OFS()<<" "<<this->Order_.perm[src_first_col+i-1];
          }
          logfileptr->OFS()<<std::endl;
          logfileptr->OFS()<<std::endl;
          for(int blkidx=0;blkidx<src_snode->NZBlockCnt();++blkidx){

            NZBlockDesc & desc = src_snode->GetNZBlockDesc(blkidx);
            T * val = src_snode->GetNZval(desc.Offset);
            Int nRows = src_snode->NRows(blkidx);

            Int row = desc.GIndex;
            for(Int i = 0; i< nRows; ++i){
              logfileptr->OFS()<<row+i<<" | "<<this->Order_.perm[row+i-1]<<":   ";
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
  //        Int NSuper = this->Xsuper_.size()-1;
  //        vector<Int> bufSupCnt(NSuper);
  //        vector<bool> bufSupMarker(NSuper,false);
  //        vector<size_t> bufSupPos(NSuper,0);
  //        //vector<Idx> lastLocalSup(this->np,0);
  //#endif
  //#endif
  //
  //
  //#ifndef BETTER_DISTRIB
  //        for(Int p=0;p<this->np;++p){
  //          send_map[p].first = 0;
  //        }
  //#endif
  //
  //        SYMPACK_TIMER_START(DISTRIBUTE_COUNTING);
  //
  //        Int snodeCount = 0;
  //        for(Int I=1;I<this->Xsuper_.size();I++){
  //          Idx fc = this->Xsuper_[I-1];
  //          Idx lc = this->Xsuper_[I]-1;
  //          Int iWidth = lc-fc+1;
  //          Int iHeight = UpdateHeight_[I-1];
  //
  //          Int iDest = this->Mapping_->Map(I-1,I-1);
  //
  //          if(iDest==this->iam){
  //            ++snodeCount;
  //          }
  //
  //          //look at the owner of the first column of the supernode
  //          Idx numColFirst = std::max(1,this->iSize_ / this->np);
  //
  //          //post all the recv and sends
  //          for(Idx col = fc;col<=lc;col++){
  //            //corresponding column in the unsorted matrix A
  //            Idx orig_col = this->Order_.perm[col-1];
  //            Idx iOwnerCol = std::min((orig_col-1)/numColFirst,(Idx)this->np-1);
  //#ifndef BETTER_DISTRIB
  //            size_t & send_bytes = send_map[iDest].first;
  //#endif
  //            if(this->iam == iOwnerCol){
  //              Int nrows = 0;
  //              Idx local_col = (orig_col-(numColFirst)*iOwnerCol);
  //              for(Ptr rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
  //                Idx orig_row = pMat.Local_.rowind[rowidx-1];
  //                Idx row = this->Order_.invp[orig_row-1];
  //
  //                if(row<col){
  //                  //add the pair (col,row) to processor owning column row
  //                  Int J = this->SupMembership_[row-1];
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
  //        remoteFactors_.resize(this->Xsuper_.size()-1);
  //        std::fill((char*)&remoteFactors_[0],(char*)&remoteFactors_[0]+remoteFactors_.size()*sizeof(std::tuple<upcxx::global_ptr<SuperNodeDesc>,Int> ),0);
  //
  //        //logfileptr->OFS()<<"My usable global memory size is: "<<upcxx::my_usable_global_memory_size()<<std::endl;
  //
  //        SYMPACK_TIMER_START(DISTRIBUTE_CREATE_SNODES);
  //
  //        for(Int I=1;I<this->Xsuper_.size();I++){
  //
  //          Int iDest = this->Mapping_->Map(I-1,I-1);
  //
  //          //parse the first column to create the supernode structure
  //          if(this->iam==iDest){
  //            Int fc = this->Xsuper_[I-1];
  //            Int lc = this->Xsuper_[I]-1;
  //            Int iWidth = lc-fc+1;
  //            Int iHeight = UpdateHeight_[I-1];
  //            Int nzBlockCnt = numBlk_[I-1];
  //#ifndef ITREE2
  //            globToLocSnodes_.push_back(I-1);
  //#else 
  //            ITree::Interval snode_inter = { I, I, LocalSupernodes_.size() };
  //            globToLocSnodes_.Insert(snode_inter);
  //#endif
  //            SuperNode<T> * newSnode = CreateSuperNode(this->options_.decomposition,I,fc,lc,iHeight,this->iSize_,nzBlockCnt);
  //            LocalSupernodes_.push_back(newSnode);
  //          }
  //        }
  //
  //
  //        SYMPACK_TIMER_STOP(DISTRIBUTE_CREATE_SNODES);
  //
  //        //first create the structure of every supernode
  //        {
  //          std::vector< int > superStructure(this->np);
  //          std::vector< int > rdisplsStructure;
  //          std::vector< int > sSuperStructure(this->np);
  //          std::vector< int > ssizes(this->np,0);
  //          std::vector< int > sdispls(this->np+1,0);
  //          std::vector< int > rsizes(this->np,0);
  //
  //          Int numLocSnode = this->XsuperDist_[this->iam+1]-this->XsuperDist_[this->iam];
  //          Int firstSnode = this->XsuperDist_[this->iam];
  //
  //          for(Int locsupno = 1; locsupno<this->locXlindx_.size(); ++locsupno){
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
  //          for(Int locsupno = 1; locsupno<this->locXlindx_.size(); ++locsupno){
  //            Idx I = locsupno + firstSnode-1;
  //            Int iDest = this->Mapping_->Map(I-1,I-1);
  //
  //            int & tail = rdisplsStructure[iDest];
  //
  //
  //            Int fc = this->Xsuper_[I-1];
  //            Int lc = this->Xsuper_[I]-1;
  //            Int iWidth = lc - fc + 1;
  //            Ptr lfi = this->locXlindx_[locsupno-1];
  //            Ptr lli = this->locXlindx_[locsupno]-1;
  //
  //            sSuperStructure[tail++] = I;
  //            //count number of contiguous rows
  //            for(Ptr sidx = lfi; sidx<=lli;sidx++){
  //              Idx iStartRow = this->locLindx_[sidx-1];
  //              Idx iPrevRow = iStartRow;
  //              Int iContiguousRows = 1;
  //              for(Int idx2 = sidx+1; idx2<=lli;idx2++){
  //                Idx iCurRow = this->locLindx_[idx2-1];
  //                if(iStartRow == this->locLindx_[lfi-1]){
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
  //          for(Int p = 0; p<this->np; p++){
  //            int pos = rdisplsStructure[p]/sizeof(int);
  //            int end = rdisplsStructure[p+1]/sizeof(int);
  //            while(pos<end){
  //              Int I = superStructure[pos++];
  //              Int nzBlockCnt = numBlk_[I-1];
  //
  //              Int fc = this->Xsuper_[I-1];
  //              Int lc = this->Xsuper_[I]-1;
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
  //          std::vector<int> sdispls(this->np+1,0);
  //          std::vector<int> scounts(this->np,0);
  //#ifndef BETTER_DISTRIB
  //          for(auto it = send_map.begin(); it!=send_map.end();it++){
  //            Int iCurDest = it->first;
  //            size_t & send_bytes = it->second.first;
  //            scounts[iCurDest] = send_bytes;
  //#ifdef _DEBUG_
  //            logfileptr->OFS()<<"P"<<this->iam<<" ---- "<<send_bytes<<" ---> P"<<iCurDest<<std::endl;
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
  //            for(Int p = bufSupDesc.rbegin()->first+1; p<this->np; ++p){
  //              scounts[p] = 0;
  //            }
  //          }
  //          else{
  //            for(Int p = 0; p<this->np; ++p){
  //              scounts[p] = 0;
  //            }
  //          }
  //          total_send_size = 0;
  //#else
  //
  //          for(Int p = 0; p<this->np; ++p){
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
  //            vector<int> spos(this->np,0); 
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
  //          //        for(Int p = 0; p<this->np; ++p){
  //          //            scounts[p] = 0;
  //          //        }
  //          //
  //          //        for(auto colIt = colOffset.begin() ; colIt != colOffset.end(); colIt++){
  //          //          Idx col = colIt->first;
  //          //          auto & triplet = colIt->second;
  //          //
  //          //          Int I = this->SupMembership_[I-1];
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
  //          for(Int I=1;I<this->Xsuper_.size();I++){
  //            Idx fc = this->Xsuper_[I-1];
  //            Idx lc = this->Xsuper_[I]-1;
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
  //            Int numColFirst = std::max(1,this->iSize_ / this->np);
  //
  //            //Serialize
  //            for(Idx col = fc;col<=lc;col++){
  //              Idx orig_col = this->Order_.perm[col-1];
  //              Int iOwnerCol = std::min(Int((orig_col-1)/numColFirst),this->np-1);
  //              if(this->iam == iOwnerCol){
  //                Int nrows = 0;
  //                Int local_col = (orig_col-(numColFirst)*iOwnerCol);
  //                for(Ptr rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
  //                  Idx orig_row = pMat.Local_.rowind[rowidx-1];
  //                  Idx row = this->Order_.invp[orig_row-1];
  //
  //
  //                  if(row<col){
  //                    //add the pair (col,row) to processor owning column row
  //                    Int J = this->SupMembership_[row-1];
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
  //                  Int row = this->Order_.invp[orig_row-1];
  //                  if(row>=col){
  //                    Isend<<row;
  //                  }
  //                }
  //
  //                for(Int rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
  //                  Int orig_row = pMat.Local_.rowind[rowidx-1];
  //                  Int row = this->Order_.invp[orig_row-1];
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
  //          Idx numColFirst = (std::max(1,this->iSize_ / this->np));
  //          for(Idx local_col = 1 ; local_col < pMat.Local_.colptr.size(); local_col++){
  //            Idx orig_col = this->iam * numColFirst + local_col;
  //            Idx col = this->Order_.invp[orig_col-1];
  //
  //            Icomm & Isend = *IsendPtr;
  //            Int I = this->SupMembership_[col-1];
  //            Int iDest = this->Mapping_->Map(I-1,I-1);
  //            auto & posSup=bufSupDesc[iDest][I];
  //
  //            volatile int found = 0;
  //            for(Ptr rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
  //              Idx orig_row = pMat.Local_.rowind[rowidx-1];
  //              Idx row = this->Order_.invp[orig_row-1];
  //              T & val = pMat.nzvalLocal[rowidx-1];
  //
  //
  //              if(row<col){
  //                //add the pair (col,row) to processor owning column row
  //                Int destCol = row;
  //                Int destRow = col;
  //
  //                Int J = this->SupMembership_[row-1];
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
  //                Idx row = this->Order_.invp[orig_row-1];
  //                if(row>=col){
  //                  Isend<<row;
  //                }
  //              }
  //              bufPos+=nrows*sizeof(Idx);
  //
  //              Isend.setHead(nzvalPos);
  //              for(Ptr rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
  //                Idx orig_row = pMat.Local_.rowind[rowidx-1];
  //                Idx row = this->Order_.invp[orig_row-1];
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
  //          for(Int I=1;I<this->Xsuper_.size();I++){
  //            Idx fc = this->Xsuper_[I-1];
  //            Idx lc = this->Xsuper_[I]-1;
  //            Int iWidth = lc-fc+1;
  //            Int iHeight = UpdateHeight_[I-1];
  //
  //            Int iDest = this->Mapping_->Map(I-1,I-1);
  //
  //            Icomm & Isend = *IsendPtr;
  //            //look at the owner of the first column of the supernode
  //            Idx numColFirst = std::max(1,this->iSize_ / this->np);
  //
  //            //post all the recv and sends
  //            for(Idx col = fc;col<=lc;col++){
  //              //corresponding column in the unsorted matrix A
  //              Idx orig_col = this->Order_.perm[col-1];
  //              Idx iOwnerCol = std::min((orig_col-1)/numColFirst,(Idx)this->np-1);
  //              if(this->iam == iOwnerCol){
  //                Int nrows = 0;
  //                Idx local_col = (orig_col-(numColFirst)*iOwnerCol);
  //                for(Ptr rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
  //                  Idx orig_row = pMat.Local_.rowind[rowidx-1];
  //                  Idx row = this->Order_.invp[orig_row-1];
  //                  T & val = pMat.nzvalLocal[rowidx-1];
  //
  //                  Int destCol = col;
  //                  Int destRow = row;
  //                  if(row<col){
  //                    destCol = row;
  //                    destRow = col;
  //                  }
  //                  //add the pair (destCol,destRow) to processor owning column destCol
  //                  Int J = this->SupMembership_[destCol-1];
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
  //          std::vector<int> rdispls(this->np+1,0);
  //          std::vector<int> rcounts(this->np,0);
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
  //          for(Int p=0;p<this->np;++p){
  //            IrecvPtr->setHead(rdispls[p]);
  //            while(IrecvPtr->head < rdispls[p]+rcounts[p]){ 
  //
  //              char * buffer = IrecvPtr->back();
  //              //Deserialize
  //              Idx col = *(Int*)&buffer[0];
  //              Int I = this->SupMembership_[col-1];
  //              //do a global to local mapping
  //
  //
  //              Int iDest = this->Mapping_->Map(I-1,I-1);
  //              bassert(this->iam==iDest);
  //
  //              Int fc = this->Xsuper_[I-1];
  //              Int lc = this->Xsuper_[I]-1;
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
  //          for(Int p=0;p<this->np;++p){
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
  //              Int fc = this->Xsuper_[I-1];
  //              Int lc = this->Xsuper_[I]-1;
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
  //            std::sort(&buffer[begin],&buffer[begin]+end,[&supM = this->SupMembership_](const triplet<T> & a, const triplet<T> & b)->bool{
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
  //            Int I = this->SupMembership_[cur_elem->col-1];
  //
  //            Int iDest = this->Mapping_->Map(I-1,I-1);
  //            bassert(this->iam==iDest);
  //
  //            if(I!=prevI){
  //              fc = this->Xsuper_[I-1];
  //              lc = this->Xsuper_[I]-1;
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
  //              for(Int I=1;I<this->Xsuper_.size();I++){
  //                Idx fc = this->Xsuper_[I-1];
  //                Idx lc = this->Xsuper_[I]-1;
  //                Int iWidth = lc-fc+1;
  //                Int iHeight = UpdateHeight_[I-1];
  //
  //                Int iDest = this->Mapping_->Map(I-1,I-1);
  //
  //                //look at the owner of the first column of the supernode
  //                Idx numColFirst = std::max(1,this->iSize_ / this->np);
  //
  //                auto & setI = bufSup[iDest][I];
  //
  //                //post all the recv and sends
  //                for(Idx col = fc;col<=lc;col++){
  //                  //corresponding column in the unsorted matrix A
  //                  Idx orig_col = this->Order_.perm[col-1];
  //                  Idx iOwnerCol = std::min((orig_col-1)/numColFirst,(Idx)this->np-1);
  //                  if(this->iam == iOwnerCol){
  //                    Int nrows = 0;
  //                    Idx local_col = (orig_col-(numColFirst)*iOwnerCol);
  //                    for(Ptr rowidx = pMat.Local_.colptr[local_col-1]; rowidx<pMat.Local_.colptr[local_col]; ++rowidx){
  //                      Idx orig_row = pMat.Local_.rowind[rowidx-1];
  //                      Idx row = this->Order_.invp[orig_row-1];
  //
  //                      if(row<col){
  //                        //add the pair (col,row) to processor owning column row
  //                        Int J = this->SupMembership_[row-1];
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
  //            std::vector<int> scounts(this->np,0);
  //            std::vector<int> sdispls(this->np+1,0);
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
  //                for(Int p = bufSup.rbegin()->first+1; p<this->np; ++p){
  //                  scounts[p] = 0;
  //                }
  //              }
  //              else{
  //                for(Int p = 0; p<this->np; ++p){
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
  //            std::vector<int> rdispls(this->np+1,0);
  //            std::vector<int> rcounts(this->np,0);
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
  //                std::sort(&buffer[begin],&buffer[begin]+end,[&supM = this->SupMembership_](const triplet<T> & a, const triplet<T> & b)->bool{
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
  //                Int I = this->SupMembership_[cur_elem->col-1];
  //
  //                Int iDest = this->Mapping_->Map(I-1,I-1);
  //                bassert(this->iam==iDest);
  //
  //                if(I!=prevI){
  //                  fc = this->Xsuper_[I-1];
  //                  lc = this->Xsuper_[I]-1;
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













  template <typename T> inline void symPACKMatrix<T>::Init(symPACKOptions & options ){
    scope_timer(a,symPACKMatrix::Init);
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

    //Create the CommEnvironment object if necessary
    if(CommEnv_!=NULL){
      delete CommEnv_;
    }

    //    if(options.commEnv == NULL){
    //      throw std::runtime_error("The communication environment must be initialized in the options");
    //    }


    this->options_ = options;
    logfileptr->verbose = this->options_.verbose>0;

    this->all_np = 0;
    MPI_Comm_size(this->options_.MPIcomm,&this->all_np);
    MPI_Comm_rank(this->options_.MPIcomm,&this->iam);

#ifdef NEW_UPCXX
    this->iam = upcxx::rank_me();
    this->all_np = upcxx::rank_n();
#endif

    this->np = this->options_.used_procs(this->all_np);


    MPI_Comm_split(this->options_.MPIcomm,this->iam<this->np,this->iam,&this->workcomm_);
    CommEnv_ = new CommEnvironment(this->workcomm_);

//    Int new_rank = (this->iam<this->np)?this->iam:this->iam-this->np;
//    upcxx::team_all.split(this->iam<this->np,new_rank, this->team_);



    //do another split to contain P0 and all the non working processors
    if(this->all_np!=this->np){
      this->non_workcomm_ = MPI_COMM_NULL;
      MPI_Comm_split(this->options_.MPIcomm,this->iam==0||this->iam>=this->np,this->iam==0?0:this->iam-this->np+1,&this->non_workcomm_);
    }

    //SYMPACK_TIMER_START(MAPPING);

    //SYMPACK_TIMER_STOP(MAPPING);

    if(this->iam<this->np){
      //Options
      if(this->options_.maxIrecv==-1 && this->options_.maxIsend==-1){
        gMaxIrecv = -1;
      }
      else{
        gMaxIrecv = this->options_.maxIrecv + this->options_.maxIsend;
      }

      //SYMPACK_TIMER_START(SCHEDULER);
      switch(this->options_.scheduler){
        case DL:
          scheduler_ = new DLScheduler<std::list<FBTask>::iterator>();
          scheduler2_ = new DLScheduler<FBTask>();
          scheduler_new_ = std::make_shared<DLScheduler< std::shared_ptr<GenericTask> > >( );
          break;
        case MCT:
          scheduler_ = new MCTScheduler<std::list<FBTask>::iterator>();
          scheduler2_ = new MCTScheduler<FBTask>();
          //scheduler_new_.reset(new MCTScheduler< std::shared_ptr<GenericTask> >( ));
          break;
        case PR:
          scheduler_ = new PRScheduler<std::list<FBTask>::iterator>();
          scheduler2_ = new PRScheduler<FBTask>();
          //scheduler_new_.reset(new PRScheduler< std::shared_ptr<GenericTask> >( ));
          break;
        case FIFO:
          scheduler_ = new FIFOScheduler<std::list<FBTask>::iterator>();
          scheduler2_ = new FIFOScheduler<FBTask>();
          scheduler_new_ = std::make_shared<FIFOScheduler< std::shared_ptr<GenericTask> > >( );
          break;
        default:
          scheduler_ = new DLScheduler<std::list<FBTask>::iterator>();
          scheduler2_ = new DLScheduler<FBTask>();
          scheduler_new_ = std::make_shared<DLScheduler< std::shared_ptr<GenericTask> > >( );
          break;
      }
      //SYMPACK_TIMER_STOP(SCHEDULER);

      //Create an Ordering object to hold the permutation
      //this->Order_.SetCommEnvironment(CommEnv_);
      //    MPI_Barrier(CommEnv_->MPI_GetComm());
    }

  }

  template <typename T> inline void symPACKMatrix<T>::SymbolicFactorization(DistSparseMatrix<T> & pMat){
    scope_timer(a,SymbolicFactorization);

    //This has to be declared here to be able to debug ...
    std::vector<int, Mallocator<int> > xadj;
    std::vector<int, Mallocator<int> > adj;
    Idx row;
    int fc,lc,colbeg,colend,col;
    Ptr supbeg,supend,rowidx;
    Int I;

#ifndef NO_MPI
    if(this->fullcomm_!=MPI_COMM_NULL){
      MPI_Comm_free(&this->fullcomm_);
    }
    MPI_Comm_dup(pMat.comm,&this->fullcomm_);
#endif

    this->iSize_ = pMat.size;
#ifdef EXPLICIT_PERMUTE
    pMat.Localg_.SetSorted(0);
    //this is done on the full communicator
    pMat.ExpandSymmetric();

    this->graph_ = pMat.GetLocalGraph();
    this->graph_.SetBaseval(1);
    this->graph_.SetKeepDiag(0);
    this->graph_.SetSorted(1);
#else
    this->graph_ = pMat.GetLocalGraph();
    this->graph_.SetBaseval(1);
    this->graph_.SetSorted(1);
    this->graph_.ExpandSymmetric();
#endif    


    logfileptr->OFS()<<"Matrix structure expanded"<<std::endl;

    SparseMatrixGraph * sgraph = NULL;
    {
      {

        double timeSta = get_time();
        SYMPACK_TIMER_START(ORDERING);

        if(this->options_.orderingStr=="MMD"){
          this->options_.ordering = symPACK::MMD;
        }
        else if(this->options_.orderingStr=="RCM"){
          this->options_.ordering = symPACK::RCM;
        }
        else if(this->options_.orderingStr=="AMD"){
          this->options_.ordering = symPACK::AMD;
        }
#ifdef USE_METIS
        else if(this->options_.orderingStr=="METIS"){
          this->options_.ordering = symPACK::METIS;
        }
#endif
#ifdef USE_SCOTCH
        else if(this->options_.orderingStr=="SCOTCH"){
          this->options_.ordering = symPACK::SCOTCH;
        }
#endif
#ifdef USE_PARMETIS
        else if(this->options_.orderingStr=="PARMETIS"){
          this->options_.ordering = symPACK::PARMETIS;
        }
#endif
#ifdef USE_PTSCOTCH
        else if(this->options_.orderingStr=="PTSCOTCH"){
          this->options_.ordering = symPACK::PTSCOTCH;
        }
#endif
        else if(this->options_.orderingStr=="NATURAL"){
          this->options_.ordering = symPACK::NATURAL;
        }
        else if(this->options_.orderingStr=="USER"){
          if(this->options_.perm ==NULL){
            throw std::logic_error( "When using USER, symPACKOptions.perm must be provided.\n" );
          }
          else{
            this->options_.ordering = symPACK::USER;
          }
        }
        else{
          std::stringstream sstr;
          sstr<<"This ordering method is not supported by symPACK. Valid options are:";
          sstr<<"NATURAL MMD AMD RCM NDBOX NDGRID USER ";
#ifdef USE_SCOTCH
          sstr<<"SCOTCH ";
#endif
#ifdef USE_PTSCOTCH
          sstr<<"PTSCOTCH ";
#endif
#ifdef USE_METIS
          sstr<<"METIS ";
#endif
#ifdef USE_PARMETIS
          sstr<<"PARMETIS ";
#endif
          sstr<<std::endl;
          throw std::logic_error( sstr.str()  );

        }

        this->Order_.NpOrdering = this->options_.NpOrdering;
        switch(this->options_.ordering){
          case MMD:
            {
              sgraph = new SparseMatrixGraph();
              this->graph_.GatherStructure(*sgraph,0);
              sgraph->SetBaseval(1);
              sgraph->SetKeepDiag(0);
              this->Order_.MMD(*sgraph, this->graph_.comm);
            }
            break;

          case RCM:
            {
              sgraph = new SparseMatrixGraph();
              this->graph_.GatherStructure(*sgraph,0);
              sgraph->SetBaseval(1);
              sgraph->SetKeepDiag(0);
              this->Order_.RCM(*sgraph, this->graph_.comm);
            }
            break;

          case AMD:
            {
              sgraph = new SparseMatrixGraph();
              this->graph_.GatherStructure(*sgraph,0);
              sgraph->SetBaseval(1);
              sgraph->SetKeepDiag(0);
              this->Order_.AMD(*sgraph, this->graph_.comm);
            }
            break;

          case USER:
            {
              this->Order_.perm.resize(this->iSize_);
              //find baseval
              auto baseval = std::min_element(&this->options_.perm[0],&this->options_.perm[0]+this->iSize_);
              //now compute everything in 1 based 
              for(int i=0;i<this->Order_.perm.size();i++){this->Order_.perm[i] = this->options_.perm[i] - *baseval + 1; /*1 based*/}
              this->Order_.invp.resize(this->iSize_);
              for(int i=0;i<this->Order_.perm.size();i++){this->Order_.invp[this->Order_.perm[i]-1] = i;} 
            }
            break;

          case NDBOX:
            this->Order_.NDBOX(this->Size(), this->graph_.comm);
            break;

          case NDGRID:
            this->Order_.NDGRID(this->Size(), this->graph_.comm);
            break;

#ifdef USE_SCOTCH
          case SCOTCH:
            {

              sgraph = new SparseMatrixGraph();
              this->graph_.GatherStructure(*sgraph,0);
              sgraph->SetKeepDiag(0);
              this->Order_.SCOTCH(*sgraph, this->graph_.comm);
            }
            break;
#endif
#ifdef USE_METIS
          case METIS:
            {
              sgraph = new SparseMatrixGraph();
              this->graph_.GatherStructure(*sgraph,0);
              sgraph->SetKeepDiag(0);
              this->Order_.METIS(*sgraph, this->graph_.comm);
            }
            break;
#endif
#ifdef USE_PARMETIS
          case PARMETIS:
            {
              this->graph_.SetKeepDiag(0);
              this->Order_.PARMETIS(this->graph_);
            }
            break;
#endif
#ifdef USE_PTSCOTCH
          case PTSCOTCH:
            {
              this->graph_.SetKeepDiag(0);
              this->Order_.PTSCOTCH(this->graph_);
            }
            break;
#endif
          case NATURAL:
            this->Order_.perm.resize(this->iSize_);
            for(int i=0;i<this->Order_.perm.size();i++){this->Order_.perm[i]=i+1;} 
            this->Order_.invp = this->Order_.perm;
            break;
          default:
            {
              std::stringstream sstr;
              sstr<<"This ordering method is not supported by symPACK. Valid options are:";
              sstr<<"MMD AMD RCM NDBOX NDGRID USER ";
#ifdef USE_SCOTCH
              sstr<<"SCOTCH ";
#endif
#ifdef USE_PTSCOTCH
              sstr<<"PTSCOTCH ";
#endif
#ifdef USE_METIS
              sstr<<"METIS ";
#endif
#ifdef USE_PARMETIS
              sstr<<"PARMETIS ";
#endif
              sstr<<std::endl;
              throw std::logic_error( sstr.str()  );
              //do nothing: either natural or user provided ordering
            }
            break;
        }
        SYMPACK_TIMER_STOP(ORDERING);

        //The ordering is available on every processor of the full communicator
        double timeStop = get_time();
        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Ordering time: "<<timeStop - timeSta<<std::endl;
        }
        logfileptr->OFS()<<"Ordering done"<<std::endl;
      }
    }

    if(this->options_.dumpPerm>0){
      logfileptr->OFS()<<"perm = [";
      for (auto i : this->Order_.perm){ 
        logfileptr->OFS()<<i<<" ";
      }
      logfileptr->OFS()<<"]"<<std::endl;
    }

    std::vector<Int> & cc = cc_;
    std::vector<Int> &rc = rc_;

    if(sgraph==NULL){ 
      logfileptr->OFS()<<"copying graph"<<std::endl;
      DistSparseMatrixGraph graph = this->graph_;
      double timeSta = get_time();

      graph.SetKeepDiag(1);
      graph.SetBaseval(1);
      logfileptr->OFS()<<"graph copied"<<std::endl;

      int expandedGraph = graph.expanded;
      //Expand to unsymmetric storage
      graph.ExpandSymmetric();

      logfileptr->OFS()<<"graph expanded"<<std::endl;
      //for(Idx i = 1; i<=graph.LocalVertexCount(); ++i){
      //  Int node = this->Order_.perm[i-1];
      //  logfileptr->OFS()<<node<<": ";

      //  Ptr jstrt = graph.colptr[node-1];
      //  Ptr jstop = graph.colptr[node] - 1;
      //  for(Ptr j = jstrt; j<=jstop; ++j){
      //    Idx nbr = graph.rowind[j-1];
      //    nbr = this->Order_.invp[nbr-1];
      //    logfileptr->OFS()<<nbr<<" ";
      //  }
      //  logfileptr->OFS()<<std::endl;
      //}




      auto backinvp = this->Order_.invp;
      //auto backperm = this->Order_.perm;



      //assume matrix is not permuted yet
      graph.Permute(this->Order_.invp.data());
//      graph.Permute(this->Order_.perm.data());

      logfileptr->OFS()<<"graph permuted 1/2"<<std::endl;


      this->ETree_.ConstructETree(graph,this->Order_);
      this->ETree_.PostOrderTree(this->Order_);
      logfileptr->OFS()<<"ETree computed and postordered"<<std::endl;

      std::vector<Int> relinvp;
      this->Order_.GetRelativeInvp(backinvp,relinvp);
//      logfileptr->OFS()<<"relinvp "<<relinvp<<std::endl;

      graph.Permute(relinvp.data());


      logfileptr->OFS()<<"graph permuted 2/2"<<std::endl;


      //  logfileptr->OFS()<<std::endl;
      //  logfileptr->OFS()<<std::endl;
      //for(Idx i = 1; i<=graph.LocalVertexCount(); ++i){
      //  Int node = this->Order_.perm[i-1];
      //  logfileptr->OFS()<<node<<": ";

      //  Ptr jstrt = graph.colptr[i-1];
      //  Ptr jstop = graph.colptr[i] - 1;
      //  for(Ptr j = jstrt; j<=jstop; ++j){
      //    Idx nbr = graph.rowind[j-1];
      //    logfileptr->OFS()<<nbr<<" ";
      //  }
      //  logfileptr->OFS()<<std::endl;
      //}



      //logfileptr->OFS()<<this->ETree_<<std::endl;
      //logfileptr->OFS()<<"ETREE done"<<std::endl;

      {
        double timeStart = get_time();
        this->getLColRowCount(graph,cc,rc);
        double timeStop = get_time();
        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Column count (distributed) construction time: "<<timeStop - timeStart<<std::endl;
        }
      }

//      graph.Permute(this->Order_.perm.data());
//      //graph.Permute(this->Order_.invp.data());
//      //If graph was not expanded, restore to symmetric storage
//      if(!expandedGraph){
//        graph.ToSymmetric();
//      }

      if(this->options_.ordering != NATURAL){
        this->ETree_.SortChildren(cc,this->Order_);

        if(this->options_.dumpPerm>0){
          logfileptr->OFS()<<"perm = [";
          for (auto i : this->Order_.perm){ 
            logfileptr->OFS()<<i<<" ";
          }
          logfileptr->OFS()<<"]"<<std::endl;
        }
      }

      double timeStop = get_time();
      if(this->iam==0 && this->options_.verbose){
        std::cout<<"Elimination tree construction time: "<<timeStop - timeSta<<std::endl;
      }

      //logfileptr->OFS()<<this->ETree_<<std::endl;
      //logfileptr->OFS()<<"------------------------------------------"<<std::endl;
      //logfileptr->OFS()<<cc<<std::endl;
      //logfileptr->OFS()<<rc<<std::endl;
      //logfileptr->OFS()<<"------------------------------------------"<<std::endl;

      //this->Order_.invp = backinvp;
      //this->Order_.perm = backperm;


    }
    else
    {
      //gather the graph if necessary to build the elimination tree
      if(sgraph==NULL){ 
        sgraph = new SparseMatrixGraph();
        this->graph_.GatherStructure(*sgraph,0);
      }
      this->graph_.SetKeepDiag(1);
      sgraph->SetBaseval(1);
      sgraph->SetKeepDiag(1);



      {

        double timeSta = get_time();
        this->ETree_.ConstructETree(*sgraph,this->Order_,this->fullcomm_);
        this->ETree_.PostOrderTree(this->Order_);

      //for(Idx i = 1; i<=sgraph->LocalVertexCount(); ++i){
      //  Int node = this->Order_.perm[i-1];
      //  logfileptr->OFS()<<node<<": ";

      //  Ptr jstrt = sgraph->colptr[node-1];
      //  Ptr jstop = sgraph->colptr[node] - 1;
      //  for(Ptr j = jstrt; j<=jstop; ++j){
      //    Idx nbr = sgraph->rowind[j-1];
      //    nbr = this->Order_.invp[nbr-1];
      //    logfileptr->OFS()<<nbr<<" ";
      //  }
      //  logfileptr->OFS()<<std::endl;
      //}


        //logfileptr->OFS()<<this->ETree_<<std::endl;
        //logfileptr->OFS()<<"ETREE done"<<std::endl;



        {
          double timeStart = get_time();
          this->getLColRowCount(*sgraph,cc,rc);
          double timeStop = get_time();
          if(this->iam==0 && this->options_.verbose){
            std::cout<<"Column count (gather + serial + bcast) construction time: "<<timeStop - timeStart<<std::endl;
          }
        }

        if(this->options_.ordering != NATURAL){
          this->ETree_.SortChildren(cc,this->Order_);
          if(this->options_.dumpPerm>0){
            logfileptr->OFS()<<"perm = [";
            for (auto i : this->Order_.perm){ 
              logfileptr->OFS()<<i<<" ";
            }
            logfileptr->OFS()<<"]"<<std::endl;
          }
        }

        double timeStop = get_time();
        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Elimination tree construction time: "<<timeStop - timeSta<<std::endl;
        }
      }

      //logfileptr->OFS()<<"------------------------------------------"<<std::endl;
      //logfileptr->OFS()<<cc<<std::endl;
      //logfileptr->OFS()<<rc<<std::endl;
      //logfileptr->OFS()<<"------------------------------------------"<<std::endl;

    }

    //compute some statistics
    if(this->iam==0 && this->options_.verbose){
      double flops = 0.0;
      int64_t NNZ = 0;
      for(Int i = 0; i<cc.size();++i){
        flops+= (double)pow((double)cc[i],2.0);
        NNZ+=cc[i];
      }
      std::cout<<"Flops: "<<flops<<std::endl;
      std::cout<<"NNZ in L factor: "<<NNZ<<std::endl;
    }


    //get rid of the sequential graph
    if(sgraph!=NULL && this->options_.order_refinement_str.substr(0,4) != "TSPB"){delete sgraph;}

    { 
      double timeSta = get_time();
      this->findSupernodes(this->ETree_,this->Order_,cc,this->SupMembership_,this->Xsuper_,this->options_.relax.maxSize);
      logfileptr->OFS()<<"Supernodes found"<<std::endl;

      if(this->options_.relax.nrelax0>0)
      {
        this->relaxSupernodes(this->ETree_, cc,this->SupMembership_, this->Xsuper_, this->options_.relax );
        logfileptr->OFS()<<"Relaxation done"<<std::endl;
      }

      //modify this->np since it cannot be greater than the number of supernodes
      this->np = std::min(this->np,this->options_.used_procs(this->Xsuper_.size()-1)); 

      if(this->workcomm_!=MPI_COMM_NULL){
        MPI_Comm_free(&this->workcomm_);
      }
      MPI_Comm_split(this->options_.MPIcomm,this->iam<this->np,this->iam,&this->workcomm_);
      bassert(CommEnv_!=NULL);
      delete CommEnv_;
      CommEnv_ = new CommEnvironment(this->workcomm_);
      this->group_.reset( new RankGroup( this->workcomm_ ) ); 

      //do another split to contain P0 and all the non working processors
      if(this->all_np!=this->np){
        if(this->non_workcomm_!=MPI_COMM_NULL){
          MPI_Comm_free(&this->non_workcomm_);
        }
        this->non_workcomm_ = MPI_COMM_NULL;
        MPI_Comm_split(this->options_.MPIcomm,this->iam==0||this->iam>=this->np,this->iam==0?0:this->iam-this->np+1,&this->non_workcomm_);
      }

      //now create the mapping
      if(this->Mapping_ != NULL){
        delete this->Mapping_;
      }

      Int pmapping = this->np;
      Int pr = (Int)sqrt((double)pmapping);
      if(this->options_.mappingTypeStr ==  "MODWRAP2DTREE"){
        this->Mapping_ = new ModWrap2DTreeMapping(pmapping, pr, pr, 1);
      }
      else if(this->options_.mappingTypeStr ==  "WRAP2D"){
        this->Mapping_ = new Wrap2D(pmapping, pr, pr, 1);
      }
      else if(this->options_.mappingTypeStr ==  "WRAP2DFORCED"){
        this->Mapping_ = new Wrap2DForced(pmapping, pr, pr, 1);
      }
      else if(this->options_.mappingTypeStr ==  "MODWRAP2DNS"){
        this->Mapping_ = new Modwrap2DNS(pmapping, pr, pr, 1);
      }
      else if(this->options_.mappingTypeStr ==  "ROW2D"){
        this->Mapping_ = new Row2D(pmapping, pmapping, pmapping, 1);
      }
      else if(this->options_.mappingTypeStr ==  "COL2D"){
        this->Mapping_ = new Col2D(pmapping, pmapping, pmapping, 1);
      }
      else{
        this->Mapping_ = new Modwrap2D(pmapping, pr, pr, 1);
      }

      //Compute this->XsuperDist_
      std::vector<Idx> newVertexDist;
      {
        Idx supPerProc = std::max((size_t)1,(this->Xsuper_.size()-1) / this->np);
        this->XsuperDist_.resize(this->all_np+1,0);
        this->XsuperDist_[this->np] = this->Xsuper_.size();
        for(int p = 0; p<this->np; p++){
          this->XsuperDist_[p]= std::min(this->XsuperDist_[this->np],(Int)(p*supPerProc+1));
        }
        for(int p = this->np+1; p<this->all_np; p++){
          this->XsuperDist_[p]= this->XsuperDist_[p-1];
        }
        this->XsuperDist_[this->all_np] = this->Xsuper_.size();


        newVertexDist.resize(this->all_np+1,0);
        newVertexDist[this->all_np] = this->iSize_+1;
        for(int p = 0; p < this->all_np; p++){
          Int S = this->XsuperDist_[p];
          newVertexDist[p] = this->Xsuper_[S-1];
        }
      }


#ifdef EXPLICIT_PERMUTE
      {
        //This should permute the matrix A and move it
        //onto this->np processor instead of this->all_np processors
        pMat.Permute(&this->Order_.invp[0],&newVertexDist[0]);
        this->graph_ = pMat.GetLocalGraph();
        this->graph_.SetBaseval(1);
        this->graph_.SetKeepDiag(1);
        this->graph_.SetSorted(1);
      }
#endif

      double timeStaSymb = get_time();
      this->symbolicFactorizationRelaxedDist(cc);

      double timeStopSymb = get_time();
      if(this->iam==0 && this->options_.verbose){
        std::cout<<"Symbolic factorization time: "<<timeStopSymb - timeStaSymb<<std::endl;
      }
      logfileptr->OFS()<<"Symbfact done"<<std::endl;


      if(this->options_.order_refinement_str != "NO") {
        double timeSta = get_time();
        if(this->options_.order_refinement_str == "SET"){ 
          this->refineSupernodes(3,1,&pMat);
        }
        else if(this->options_.order_refinement_str == "SET10"){ 
          this->refineSupernodes(1,0,&pMat);
        }
        else if(this->options_.order_refinement_str == "SET11"){ 
          this->refineSupernodes(1,1,&pMat);
        }
        else if(this->options_.order_refinement_str == "SET20"){ 
          this->refineSupernodes(2,0,&pMat);
        }
        else if(this->options_.order_refinement_str == "SET21"){ 
          this->refineSupernodes(2,1,&pMat);
        }
        else if(this->options_.order_refinement_str == "SET30"){ 
          this->refineSupernodes(3,0,&pMat);
        }
        else if(this->options_.order_refinement_str == "SET31"){ 
          this->refineSupernodes(3,1,&pMat);
        }
        else if(this->options_.order_refinement_str == "TSP"){
          auto SupETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);

          std::vector<Ptr> xlindx;
          std::vector<Idx> lindx;
          this->gatherLStructure(xlindx, lindx);

          if(this->iam==0){
            TSP::SymbolMatrix * symbmtx = TSP::GetPastixSymbolMatrix(this->Xsuper_,this->SupMembership_, xlindx, lindx);
            TSP::Order * psorder = TSP::GetPastixOrder(symbmtx,this->Xsuper_, SupETree, &this->Order_.perm[0], &this->Order_.invp[0]);

            double timeSta = get_time();
            TSP::symbolReordering( symbmtx, psorder, 0, std::numeric_limits<int>::max(), 0 );
            double timeStop = get_time();
            if(this->iam==0 && this->options_.verbose){
              std::cout<<"TSP reordering done in "<<timeStop-timeSta<<std::endl;
            }

            //overwrite order
            for(int i = 0; i < this->Order_.perm.size(); ++i){
              this->Order_.perm[i] = psorder->peritab[i]+1;
              this->Order_.invp[i] = psorder->permtab[i]+1;
            }
          }

          MPI_Bcast(this->Order_.perm.data(),this->Order_.perm.size()*sizeof(Int),MPI_BYTE,0,this->fullcomm_);
          MPI_Bcast(this->Order_.invp.data(),this->Order_.invp.size()*sizeof(Int),MPI_BYTE,0,this->fullcomm_);
        } 
        else if(this->options_.order_refinement_str.substr(0,4) == "TSPB"){

          if(sgraph==NULL){ 
            sgraph = new SparseMatrixGraph();
            this->graph_.GatherStructure(*sgraph,0);
          }

          ETree& tree = this->ETree_;
          Ordering & aOrder = this->Order_;
          std::vector<Int> & supMembership = this->SupMembership_; 
          std::vector<Int> & xsuper = this->Xsuper_; 

          std::vector<int> ixlindx;
          std::vector<int> ilindx;

          std::vector<int>  new_invp;

          //Gather this->locXlindx_ and this->locLindx_
          {
            std::vector<Ptr> xlindx;
            std::vector<Idx> lindx;
            this->gatherLStructure(xlindx, lindx);

            if(this->iam==0){
              ixlindx.resize(xlindx.size());
              for(int i = 0;i<xlindx.size();i++){
                ixlindx[i] = xlindx[i];
              }
              ilindx.resize(lindx.size());
              for(int i = 0;i<lindx.size();i++){
                ilindx[i] = lindx[i];
              }
            }
          }

          if(this->iam==0){
            bassert(sgraph!=nullptr);
            int neqns = this->iSize_;

            int nofsub =ilindx.size();
            int nsuper = xsuper.size()-1;


            new_invp.assign(neqns,0);


            std::vector<int>  new_perm(neqns,0);

            for(size_t i =0; i<new_perm.size(); i++){ new_invp[i] = this->Order_.invp[i];}
            for(size_t i =0; i<new_perm.size(); i++){ new_perm[i] = this->Order_.perm[i];}

            int supsiz = 0;
            for(I = 1; I <= nsuper; I++){
              Int fc = this->Xsuper_[I-1];
              Int lc = this->Xsuper_[I]-1;
              supsiz = std::max(supsiz,lc-fc+1);
            }

            sgraph->SetKeepDiag(0);
            xadj.resize(sgraph->colptr.size());
            adj.resize(sgraph->rowind.size());
            int nadj = adj.size();
            for(size_t i = 0; i< sgraph->colptr.size(); i++){ xadj[i] = int(sgraph->colptr[i]); }
            for(size_t i = 0; i< sgraph->rowind.size(); i++){ adj[i] = int(sgraph->rowind[i]); }

            if(sgraph!=NULL){delete sgraph;}

            std::vector<int, Mallocator<int> > etpar(neqns);
            for(int i = 0; i<neqns; i++){ etpar[i] = tree.PostParent(i); }


            std::vector<int, Mallocator<int> > xskadj(neqns+1);
            std::vector<int, Mallocator<int> > sklenf(neqns);
            std::vector<int, Mallocator<int> > sklenb(neqns);
            std::vector<int, Mallocator<int> > skadj(nadj);
            std::vector<int, Mallocator<int> > invp2(neqns);
            std::vector<int, Mallocator<int> > link(neqns);
            std::vector<int, Mallocator<int> > fstloc(nsuper);
            std::vector<int, Mallocator<int> > sperm(nsuper);
            std::vector<int, Mallocator<int> > fstloc2(neqns);
            std::vector<int, Mallocator<int> > dist1((supsiz+1)*(supsiz+1));
            std::vector<int, Mallocator<int> > suppar(nsuper);
            int iwsiz = 8*neqns+3;
            std::vector<int, Mallocator<int> > iwork(iwsiz);

            int iflag = 0;
            double timeSta = get_time();
            if(this->options_.order_refinement_str == "TSPB"){
              std::vector<int, Mallocator<int> > rep(neqns);
              FORTRAN(ordsup_ind_tsp_paths2)
                ( 
                 &nadj, &neqns , &nofsub, &nsuper, &supsiz, 
                 xsuper.data(), ixlindx.data(), ilindx.data(), supMembership.data(), 
                 xadj.data(), adj.data(), 
                 etpar.data(),
                 new_perm.data(), new_invp.data(), 
                 &iflag , 
                 xskadj.data(), sklenf.data(), sklenb.data(), 
                 skadj.data() , invp2.data() , link.data()  , 
                 fstloc.data(), sperm.data() , fstloc2.data(), dist1.data(), 
                 suppar.data(), &iwsiz , iwork.data() , rep.data()             );
            }
            else{
              FORTRAN(ordsup_ind_tsp_paths)
                ( 
                 &nadj, &neqns , &nofsub, &nsuper, &supsiz, 
                 xsuper.data(), ixlindx.data(), ilindx.data(), supMembership.data(), 
                 xadj.data(), adj.data(), 
                 etpar.data(),
                 new_perm.data(), new_invp.data(), 
                 &iflag , 
                 xskadj.data(), sklenf.data(), sklenb.data(), 
                 skadj.data() , invp2.data() , link.data()  , 
                 fstloc.data(), sperm.data() , fstloc2.data(), dist1.data(), 
                 suppar.data(), &iwsiz , iwork.data() );
            }

            double timeStop = get_time();
            if(this->iam==0 && this->options_.verbose){
              std::cout<<"TSPB reordering done in "<<timeStop-timeSta<<std::endl;
            }


            //this->Order_.Compose(new_invp);
            for(size_t i =0; i<new_perm.size(); i++){ this->Order_.invp[i] = new_invp[i];}
            for(size_t i =0; i<new_perm.size(); i++){ this->Order_.perm[i] = new_perm[i];}
          }

          // broadcast invp
          Int N = aOrder.invp.size();
          MPI_Bcast(&aOrder.invp[0],N*sizeof(Int),MPI_BYTE,0,this->fullcomm_);
          MPI_Bcast(&aOrder.perm[0],N*sizeof(Int),MPI_BYTE,0,this->fullcomm_);

        }

        double timeStop = get_time();

        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Supernode reordering done in "<<timeStop-timeSta<<std::endl;
        }

        {
          double timeSta = get_time();
          this->symbolicFactorizationRelaxedDist(cc);
          double timeStop = get_time();
          if(this->iam==0 && this->options_.verbose){
            std::cout<<"Symbolic factorization time: "<<timeStop - timeSta<<std::endl;
          }
          logfileptr->OFS()<<"Symbfact done"<<std::endl;
        }
      }

      double timeStop = get_time();
      if(this->iam==0 && this->options_.verbose){
        std::cout<<"Total symbolic factorization time: "<<timeStop - timeSta<<std::endl;
      }

      //Print statistics
      if(this->options_.print_stats){
        OrderStats stats;
        stats.get(this->Xsuper_, this->XsuperDist_, this->locXlindx_, this->locLindx_, this->fullcomm_);
        if (this->iam==0){
          stats.print();
        }
      }
    }

#ifdef _DEBUG_
    logfileptr->OFS()<<"Membership list is "<<this->SupMembership_<<std::endl;
    logfileptr->OFS()<<"xsuper "<<this->Xsuper_<<std::endl;
#endif

    {
      scope_timer(a,LOAD_BALANCE);
      double timeSta = get_time();

      if (this->Balancer_!=NULL){
        delete this->Balancer_;
      }

      std::vector<Int> map;
      if(this->options_.load_balance_str=="SUBCUBE-FI"){
        if(this->iam==0 && this->options_.verbose){ std::cout<<"Subtree to subcube FI mapping used"<<std::endl;}
        ETree SupETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);
        this->Balancer_ = new SubtreeToSubcube(this->np,SupETree,this->Xsuper_,this->XsuperDist_,this->SupMembership_,this->locXlindx_,this->locLindx_,cc,this->fullcomm_,true);
      }
      else if(this->options_.load_balance_str=="SUBCUBE-FO"){
        if(this->iam==0 && this->options_.verbose){ std::cout<<"Subtree to subcube FO mapping used"<<std::endl;}
        ETree SupETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);
        this->Balancer_ = new SubtreeToSubcube(this->np,SupETree,this->Xsuper_,this->XsuperDist_,this->SupMembership_,this->locXlindx_,this->locLindx_,cc,this->fullcomm_,false);
      }
      else if(this->options_.load_balance_str=="SUBCUBE-VOLUME-FI"){
        if(this->iam==0 && this->options_.verbose){ std::cout<<"Subtree to subcube volume FI mapping used"<<std::endl;}
        ETree SupETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);
        this->Balancer_ = new SubtreeToSubcubeVolume(this->np,SupETree,this->Xsuper_,this->XsuperDist_,this->SupMembership_,this->locXlindx_,this->locLindx_,cc,this->fullcomm_,true);
      }
      else if(this->options_.load_balance_str=="SUBCUBE-VOLUME-FO"){
        if(this->iam==0 && this->options_.verbose){ std::cout<<"Subtree to subcube volume FO mapping used"<<std::endl;}
        ETree SupETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);
        this->Balancer_ = new SubtreeToSubcubeVolume(this->np,SupETree,this->Xsuper_,this->XsuperDist_,this->SupMembership_,this->locXlindx_,this->locLindx_,cc,this->fullcomm_,false);
      }
      else if(this->options_.load_balance_str=="NNZ"){
        if(this->iam==0 && this->options_.verbose){ std::cout<<"Load Balancing on NNZ used"<<std::endl;}
        this->Balancer_ = new NNZBalancer(this->np,this->Xsuper_,cc);
      }
      else if(this->options_.load_balance_str=="WORK"){
        if(this->iam==0 && this->options_.verbose){ std::cout<<"Load Balancing on WORK used"<<std::endl;}
        this->Balancer_ = new WorkBalancer(this->np,this->Xsuper_,cc);
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
            test->Update((TreeLoadBalancer*)this->Balancer_);
          }
        }
      }
      double timeStop = get_time();
      if(this->iam==0 && this->options_.verbose){
        std::cout<<"Load balancing time: "<<timeStop - timeSta<<std::endl;
      }
    }


    //Now split supernodes larger than the maximum size. XlindxLocal_, LindxLocal_, this->Xsuper_, this->SupMembership_, Balancer_ and Mapping_ need to up updated

#ifdef _OUTPUT_ETREE_
    logfileptr->OFS()<<"ETree is "<<this->ETree_<<std::endl;
    {
      auto supETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);
      logfileptr->OFS()<<"Supernodal ETree is "<<supETree<<std::endl;
      logfileptr->OFS()<<"Mapping is ";for(Int i = 0;i<supETree.Size();i++){logfileptr->OFS()<<this->Mapping_->Map(i,i)<<" ";}logfileptr->OFS()<<std::endl;
      logfileptr->OFS()<<"Xsuper is "<<this->Xsuper_<<std::endl;//;for(Int i = 0;i<supETree.Size();i++){logfileptr->OFS()<<this->Xsuper_->Map(i,i)<<" ";}logfileptr->OFS()<<std::endl;
    } 
#endif

#ifndef _MAP_DEBUG_
    //starting from now, this only concerns the working processors
    SYMPACK_TIMER_START(Get_UpdateCount);
    GetUpdatingSupernodeCount(UpdateCount_,UpdateWidth_,UpdateHeight_,numBlk_);
    SYMPACK_TIMER_STOP(Get_UpdateCount);

    if(this->iam<this->np){

      {
        double timeSta = get_time();
        std::vector<Int> AggregatesToRecv;
        std::vector<Int> LocalAggregates;


        FBGetUpdateCount(UpdatesToDo_,AggregatesToRecv,LocalAggregates);
        generateTaskGraph(taskGraph_, AggregatesToRecv, LocalAggregates);
          {
          utility::scope_memprofiler m("symPACKMatrix_task_graph");
        generateTaskGraph_New(taskGraph_New_, AggregatesToRecv, LocalAggregates,UpdateWidth_,UpdateHeight_);
          }


        //#define _OUTPUT_TASK_GRAPH_
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
              logfileptr->OFS()<<this->iam<<" ";
            }
          }
        }
        logfileptr->OFS()<<std::endl; 
#endif

        double timeStop = get_time();
        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Task graph generation time: "<<timeStop - timeSta<<std::endl;
        }
      }
      MPI_Barrier(CommEnv_->MPI_GetComm());
    }

    SYMPACK_TIMER_START(DISTRIBUTE_CREATE_SNODES);


    Int snodeCount = 0;
    for(Int I=1;I<this->Xsuper_.size();I++){
      Int iDest = this->Mapping_->Map(I-1,I-1);
      if(iDest==this->iam){
        ++snodeCount;
      }
    }
    //Resize the local supernodes array
    for(auto ptr : LocalSupernodes_){ delete ptr; }
    LocalSupernodes_.clear();
    globToLocSnodes_.clear();
    LocalSupernodes_.reserve(snodeCount);

    for(Int I=1;I<this->Xsuper_.size();I++){

      Int iDest = this->Mapping_->Map(I-1,I-1);

      //parse the first column to create the supernode structure
      if(this->iam==iDest){
        Int fc = this->Xsuper_[I-1];
        Int lc = this->Xsuper_[I]-1;
        Int iWidth = lc-fc+1;
        Int iHeight = UpdateHeight_[I-1];
        Int nzBlockCnt = numBlk_[I-1];
#ifndef ITREE2
        globToLocSnodes_.push_back(I-1);
#else 
        ITree::Interval snode_inter = { I, I, LocalSupernodes_.size() };
        globToLocSnodes_.Insert(snode_inter);
#endif
        SuperNode<T> * newSnode = NULL;
        try{
          newSnode = CreateSuperNode(this->options_.decomposition,I,fc,fc,lc,iHeight,this->iSize_,nzBlockCnt,this->options_.panel);
        }
        catch(const MemoryAllocationException & e){
          std::stringstream sstr;
          sstr<<"There is not enough memory on the system to allocate the factors. Try to increase GASNET_MAX_SEGSIZE value."<<std::endl;
          logfileptr->OFS()<<sstr.str();
          std::cerr<<sstr.str();
          throw(e);
        }

        LocalSupernodes_.push_back(newSnode);
      }
    }
    SYMPACK_TIMER_STOP(DISTRIBUTE_CREATE_SNODES);

    //first create the structure of every supernode

    {
      scope_timer(a,DISTRIBUTE_CREATE_REMOTE_NODES);
      std::vector< int > superStructure(this->all_np);
      std::vector< int > rdisplsStructure;
      std::vector< int > sSuperStructure(this->all_np);
      std::vector< int > ssizes(this->all_np,0);
      std::vector< int > sdispls(this->all_np+1,0);
      std::vector< int > rsizes(this->all_np,0);

      Int numLocSnode = this->XsuperDist_[this->iam+1]-this->XsuperDist_[this->iam];
      Int firstSnode = this->XsuperDist_[this->iam];

      for(Int locsupno = 1; locsupno<this->locXlindx_.size(); ++locsupno){
        Idx I = locsupno + firstSnode-1;
        Int iDest = this->Mapping_->Map(I-1,I-1);
        ssizes[iDest] += 1 + 2*numBlk_[I-1]; //1 for supno index + numBlk startrows + numBlk number of rows
      }

      sdispls[0] = 0;
      std::partial_sum(ssizes.begin(),ssizes.end(),&sdispls[1]);

      rdisplsStructure = sdispls;
      sSuperStructure.resize(sdispls.back());
      for(Int locsupno = 1; locsupno<this->locXlindx_.size(); ++locsupno){
        Idx I = locsupno + firstSnode-1;
        Int iDest = this->Mapping_->Map(I-1,I-1);
        int & tail = rdisplsStructure[iDest];

        Idx fc = this->Xsuper_[I-1];
        Idx lc = this->Xsuper_[I]-1;
        Int iWidth = lc - fc + 1;
        Ptr lfi = this->locXlindx_[locsupno-1];
        Ptr lli = this->locXlindx_[locsupno]-1;

#ifdef SPLIT_AT_BOUNDARY
        Int nextSup = I+1;
        Idx next_fc = fc;
        Idx next_lc = lc;
#endif

        sSuperStructure[tail++] = I;
        //count number of contiguous rows
        for(Ptr sidx = lfi; sidx<=lli;sidx++){
          Idx iStartRow = this->locLindx_[sidx-1];
          Idx iPrevRow = iStartRow;
          Int iContiguousRows = 1;
          for(Int idx2 = sidx+1; idx2<=lli;idx2++){
            Idx iCurRow = this->locLindx_[idx2-1];
            if(iStartRow == this->locLindx_[lfi-1]){
              if(iCurRow>iStartRow+iWidth-1){
                //enforce the first block to be a square diagonal block
                break;
              }
            }

#ifdef SPLIT_AT_BOUNDARY
            if( nextSup<=this->Xsuper_.size()-1){
              next_fc = this->Xsuper_[nextSup-1];
              next_lc = this->Xsuper_[nextSup]-1;
            }
            while(iCurRow>=next_lc){
              nextSup++;
              if( nextSup<=this->Xsuper_.size()-1){
                next_fc = this->Xsuper_[nextSup-1];
                next_lc = this->Xsuper_[nextSup]-1;
              }
              else{
                break;
              }
            } 

            if( nextSup<this->Xsuper_.size()-0 &&   iCurRow==next_fc){
              break;
            }
#endif

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

      MPI_Datatype type;
      MPI_Type_contiguous( sizeof(int), MPI_BYTE, &type );
      MPI_Type_commit(&type);

      //Do the alltoallv to get the structures        
      MPI_Alltoallv(&sSuperStructure[0], &ssizes[0], &sdispls[0], type,
          &superStructure[0], &rsizes[0], &rdisplsStructure[0], type,
          pMat.comm);

      //loop through received structure and create supernodes
      for(Int p = 0; p<this->all_np; p++){
        int pos = rdisplsStructure[p];
        int end = rdisplsStructure[p+1];
        while(pos<end){
          Int I = superStructure[pos++];
          Int nzBlockCnt = numBlk_[I-1];

          Int fc = this->Xsuper_[I-1];
          Int lc = this->Xsuper_[I]-1;
          Int iWidth = lc-fc+1;
          SuperNode<T> * snode = snodeLocal(I);
          if(snode->NZBlockCnt()==0){
            for(Int i = 0; i<nzBlockCnt;i++){
              Int iStartRow = superStructure[pos++];
              Int iContiguousRows = superStructure[pos++];
              snode->AddNZBlock(iContiguousRows , iStartRow);
            }
            snode->Shrink();
          }
        }
      } 
      MPI_Type_free(&type);
    }
#endif
    MPI_Barrier(this->fullcomm_);
  }

  template <typename T> inline void symPACKMatrix<T>::DistributeMatrix(DistSparseMatrix<T> & pMat){
    scope_timer(a,DistributeA);
#ifdef EXPLICIT_PERMUTE
    pMat.Localg_.SetSorted(0);
    pMat.ExpandSymmetric();
    //matrix should already be distributed ?
    pMat.Permute(&this->Order_.invp[0]);
#endif
    {
      double timeSta = get_time();
      //Icomm buffer;
      Icomm recv_buffer;


      AsyncComms incomingRecv;
      AsyncComms outgoingSend;
      Int numColFirst = std::max((Int)1,this->iSize_ / this->np);

      logfileptr->OFS()<<"Starting Send"<<std::endl;

#ifndef NOTRY
      try
#endif
      {
#ifdef EXPLICIT_PERMUTE
        pMat.ToLowerTriangular();
        pMat.SortGraph();

        vector<std::pair<size_t,Icomm *> > send_map(this->np);

        SYMPACK_TIMER_START(DISTRIBUTE_COUNTING);

        Idx FirstLocalCol = pMat.Localg_.vertexDist[this->iam];
        Idx LastLocalCol = pMat.Localg_.vertexDist[this->iam+1];
        if(FirstLocalCol<LastLocalCol){
          for(Int I=1;I<this->Xsuper_.size();I++){
            Idx fc = this->Xsuper_[I-1];
            Idx lc = this->Xsuper_[I]-1;
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


        size_t total_send_size = 0;
        std::vector<size_t> stotcounts(this->all_np,0);
        for(Int i = 0; i< send_map.size(); i++){
          total_send_size += send_map[i].first;
          stotcounts[i] = send_map[i].first;
        }

        Icomm * IsendPtr = new Icomm(total_send_size,MPI_REQUEST_NULL);

        std::vector<size_t> spositions(this->all_np+1,0);
        spositions[0] = 0;
        std::partial_sum(stotcounts.begin(),stotcounts.end(),&spositions[1]);

        SYMPACK_TIMER_SPECIAL_START(serializing); 
        auto baseval = pMat.Localg_.GetBaseval();
        if(FirstLocalCol<LastLocalCol){
          for(Int I=1;I<this->Xsuper_.size();I++){
            Idx fc = this->Xsuper_[I-1];
            Idx lc = this->Xsuper_[I]-1;
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
            (*IrecvPtr),this->fullcomm_, resize_lambda);

        //restore zero values in the LocalSupernodes_
        for(auto snode : LocalSupernodes_){
          snode->clear();
        }

        SYMPACK_TIMER_SPECIAL_START(deserializing);      
        IrecvPtr->setHead(0);
        while(IrecvPtr->head < IrecvPtr->capacity()){ 
          char * buffer = IrecvPtr->back();
          //Deserialize
          Idx col = *(Int*)&buffer[0];
          Int I = this->SupMembership_[col-1];

          Int iDest = this->Mapping_->Map(I-1,I-1);
          bassert(this->iam==iDest);

          Int fc = this->Xsuper_[I-1];
          Int lc = this->Xsuper_[I]-1;
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
        for(Int p=0;p<this->all_np;++p){
          send_map[p].first = 0;
        }


        typedef std::conditional< sizeof(Idx) < sizeof(Ptr), Idx, Ptr>::type minTypeIdxPtr;
        using minType = typename std::conditional< sizeof(minTypeIdxPtr) < sizeof(T), minTypeIdxPtr, T >::type;

        size_t minSize = sizeof(minType);
        size_t IdxToMin = sizeof(Idx) / minSize;
        size_t PtrToMin = sizeof(Ptr) / minSize;
        size_t TToMin = sizeof(T) / minSize;

        Int baseval = pMat.Localg_.GetBaseval();
        Idx FirstLocalCol = pMat.Localg_.vertexDist[this->iam] + (1 - baseval); //1-based
        Idx LastLocalCol = pMat.Localg_.vertexDist[this->iam+1] + (1 - baseval); //1-based

        SYMPACK_TIMER_START(DISTRIBUTE_COUNTING);
        for(Int I=1;I<this->Xsuper_.size();I++){
          Idx fc = this->Xsuper_[I-1];
          Idx lc = this->Xsuper_[I]-1;
          Int iWidth = lc-fc+1;
          Int iHeight = UpdateHeight_[I-1];

          Int iDest = this->Mapping_->Map(I-1,I-1);

          //post all the recv and sends
          for(Idx col = fc;col<=lc;col++){
            //corresponding column in the unsorted matrix A
            Idx orig_col = this->Order_.perm[col-1];
            if(orig_col>= FirstLocalCol && orig_col < LastLocalCol){
              Ptr nrows = 0;
              Idx local_col = orig_col - FirstLocalCol+1;//1 based
              for(Ptr rowidx = pMat.Localg_.colptr[local_col-1] + (1-baseval); rowidx<pMat.Localg_.colptr[local_col]+(1-baseval); ++rowidx){
                Idx orig_row = pMat.Localg_.rowind[rowidx-1]+(1-baseval);//1-based
                Idx row = this->Order_.invp[orig_row-1];

                if(row<col){
                  //add the pair (col,row) to processor owning column row
                  Int J = this->SupMembership_[row-1];
                  Int iDestJ = this->Mapping_->Map(J-1,J-1);
                  send_map[iDestJ].first += IdxToMin+PtrToMin+ (IdxToMin + TToMin);
                }
                else{nrows++;}
              }
              send_map[iDest].first += IdxToMin+PtrToMin+nrows*(IdxToMin + TToMin);
            }
          }
        }
        SYMPACK_TIMER_STOP(DISTRIBUTE_COUNTING);



        {
          //allocate one buffer for every remote processor
          //compute the send structures
          size_t total_send_size = 0;
          std::vector<size_t> stotcounts(this->all_np,0);
          for(auto it = send_map.begin(); it!=send_map.end();it++){
            Int iCurDest = it->first;
            size_t & send_bytes = it->second.first;
            stotcounts[iCurDest] = send_bytes;
          }
          //compute send displacements
          std::vector<size_t> spositions(this->all_np+1,0);
          spositions[0] = 0;
          std::partial_sum(stotcounts.begin(),stotcounts.end(),&spositions[1]);
          total_send_size = spositions.back();

          std::vector<minType, Mallocator<minType> > sendBuffer(total_send_size);

          //Fill up the send buffer
          {
          utility::scope_memprofiler m("symPACKMatrix::DistributeMatrix::Serializing");
          SYMPACK_TIMER_SPECIAL_START(serializing);      
          for(Int I=1;I<this->Xsuper_.size();I++){
            Idx fc = this->Xsuper_[I-1];
            Idx lc = this->Xsuper_[I]-1;
            Int iWidth = lc-fc+1;
            Int iHeight = cc_[fc-1];

            Int iDest = this->Mapping_->Map(I-1,I-1);
#ifdef _DEBUG_
            logfileptr->OFS()<<"Supernode "<<I<<" is owned by P"<<iDest<<std::endl;
#endif

            //Serialize
            for(Idx col = fc;col<=lc;col++){
              Idx orig_col = this->Order_.perm[col-1];

              if(orig_col>= FirstLocalCol && orig_col < LastLocalCol){
                Ptr nrows = 0;
                Idx local_col = orig_col - FirstLocalCol+1;//1 based

                for(Ptr rowidx = pMat.Localg_.colptr[local_col-1]+(1-baseval); rowidx<pMat.Localg_.colptr[local_col]+(1-baseval); ++rowidx){
                  Idx orig_row = pMat.Localg_.rowind[rowidx-1]+(1-baseval);
                  Idx row = this->Order_.invp[orig_row-1];

                  if(row<col){
                    //add the pair (col,row) to processor owning column row
                    Int J = this->SupMembership_[row-1];
                    Int iDestJ = this->Mapping_->Map(J-1,J-1);

                    T val = pMat.nzvalLocal[rowidx-1];

                    *((Idx*)&sendBuffer[spositions[iDestJ]]) = row;
                    spositions[iDestJ]+=IdxToMin;
                    *((Ptr*)&sendBuffer[spositions[iDestJ]]) = 1;
                    spositions[iDestJ]+=PtrToMin;
                    *((Idx*)&sendBuffer[spositions[iDestJ]]) = col;
                    spositions[iDestJ]+=IdxToMin;
                    *((T*)&sendBuffer[spositions[iDestJ]]) = val;
                    spositions[iDestJ]+=TToMin;

                  }
                  else{
                    nrows++;
                  }
                }

                *((Idx*)&sendBuffer[spositions[iDest]]) = col;
                spositions[iDest]+=IdxToMin;
                *((Ptr*)&sendBuffer[spositions[iDest]]) = nrows;
                spositions[iDest]+=PtrToMin;

                for(Ptr rowidx = pMat.Localg_.colptr[local_col-1]+(1-baseval); rowidx<pMat.Localg_.colptr[local_col]+(1-baseval); ++rowidx){
                  Idx orig_row = pMat.Localg_.rowind[rowidx-1]+(1-baseval);
                  Idx row = this->Order_.invp[orig_row-1];
                  if(row>=col){
                    *((Idx*)&sendBuffer[spositions[iDest]]) = row;
                    spositions[iDest]+=IdxToMin;
                  }
                }

                for(Ptr rowidx = pMat.Localg_.colptr[local_col-1]+(1-baseval); rowidx<pMat.Localg_.colptr[local_col]+(1-baseval); ++rowidx){
                  Idx orig_row = pMat.Localg_.rowind[rowidx-1]+(1-baseval);
                  Idx row = this->Order_.invp[orig_row-1];
                  if(row>=col){
                    auto & val = pMat.nzvalLocal[rowidx-1];
                    *((T*)&sendBuffer[spositions[iDest]]) = val;
                    spositions[iDest]+=TToMin;
                  }
                }

              }
            }
          }
          SYMPACK_TIMER_SPECIAL_STOP(serializing);      
          }


          spositions[0] = 0;
          std::partial_sum(stotcounts.begin(),stotcounts.end(),&spositions[1]);

          size_t total_recv_size = 0;
          std::vector<minType, Mallocator<minType> > recvBuffer;
          std::function<void(std::vector<minType, Mallocator<minType> >&,size_t)> resize_lambda =
            [](std::vector<minType, Mallocator<minType> >& container, size_t sz){
              container.resize(sz);
            };


          MPI_Datatype type;
          MPI_Type_contiguous( sizeof(minType), MPI_BYTE, &type );
          MPI_Type_commit(&type);

          mpi::Alltoallv(sendBuffer, &stotcounts[0], &spositions[0], type ,
              recvBuffer,this->fullcomm_, resize_lambda);

          total_recv_size = recvBuffer.size();

          MPI_Type_free(&type);
          //Need to parse the structure sent from the processor owning the first column of the supernode

          {
          utility::scope_memprofiler m("symPACKMatrix::DistributeMatrix::Deserializing");
          SYMPACK_TIMER_SPECIAL_START(deserializing);      
          size_t head = 0;

          while(head < total_recv_size)
          { 
            //Deserialize
            Idx col = *((Idx*)&recvBuffer[head]);
            head+=IdxToMin;
            //nrows of column col sent by processor p
            Ptr nrows = *((Ptr*)&recvBuffer[head]);
            head+=PtrToMin;
            Idx * rowind = (Idx*)(&recvBuffer[head]);
            head+=nrows*IdxToMin;
            T * nzvalA = (T*)(&recvBuffer[head]);
            head+=nrows*TToMin;

            Int I = this->SupMembership_[col-1];
            //do a global to local mapping

            Int iDest = this->Mapping_->Map(I-1,I-1);
            bassert(this->iam==iDest);

            Int fc = this->Xsuper_[I-1];
            Int lc = this->Xsuper_[I]-1;
            Int iWidth = lc-fc+1;
            SuperNode<T> * snode = snodeLocal(I);

            //Here, do a linear search instead for the blkidx
            Ptr colbeg = 1;
            Ptr colend = nrows;
            if(colbeg<=colend){
              //sort rowind and nzvals

              std::vector<size_t> lperm = sort_permutation(&rowind[colbeg-1],&rowind[colend-1]+1,std::less<Idx>());
              apply_permutation(&rowind[colbeg-1],&rowind[colend-1]+1,lperm);
              apply_permutation(&nzvalA[colbeg-1],&nzvalA[colend-1]+1,lperm);
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
          }
        }
#endif

        remoteFactors_.resize(this->Xsuper_.size()-1);
#ifdef NEW_UPCXX
        std::fill((char*)&remoteFactors_[0],(char*)&remoteFactors_[0]+remoteFactors_.size()*sizeof(std::tuple<upcxx::global_ptr<char>,Int> ),0);
#else
        std::fill((char*)&remoteFactors_[0],(char*)&remoteFactors_[0]+remoteFactors_.size()*sizeof(std::tuple<upcxx::global_ptr<SuperNodeDesc>,Int> ),0);
#endif

        for(Int I=1;I<this->Xsuper_.size();I++){
          Int iDest = this->Mapping_->Map(I-1,I-1);
          //parse the first column to create the supernode structure
          if(this->iam==iDest){
            SuperNode<T> * newSnode = snodeLocal(I);
            SuperNodeDesc * meta = newSnode->GetMeta();
#ifdef NEW_UPCXX
            remoteFactors_[I-1] = std::make_tuple( upcxx::to_global_ptr<char>( (char*)meta ), meta->blocks_cnt_) ;
#else
            remoteFactors_[I-1] = std::make_tuple( upcxx::global_ptr<SuperNodeDesc>( meta ), meta->blocks_cnt_) ;
#endif
          }
        }

        auto bor_op = []( void *in, void *inout, int *len, MPI_Datatype *dptr ){ 
          size_t i; 
          using F = std::tuple<upcxx::global_ptr<SuperNodeDesc>,Int>; 
          char * pinout = (char*)inout;
          char * pin = (char*)in;
#pragma unroll
          for (i=0; i< (*len)*sizeof(F); ++i) { 
            pinout[i] |= pin[i];
          } 
        };

        MPI_Op MPI_SYMPACK_BOR; 
        MPI_Op_create( bor_op, true, &MPI_SYMPACK_BOR ); 

        MPI_Datatype type;
        MPI_Type_contiguous( sizeof(std::tuple<upcxx::global_ptr<SuperNodeDesc>,Int> ), MPI_BYTE, &type );
        MPI_Type_commit(&type);
        MPI_Allreduce( MPI_IN_PLACE, &remoteFactors_[0], remoteFactors_.size(), type, MPI_SYMPACK_BOR, this->fullcomm_);
        MPI_Type_free(&type);
        MPI_Op_free(&MPI_SYMPACK_BOR);
      }
#ifndef NOTRY
      catch(const std::bad_alloc& e){
        std::cout << "Allocation failed: " << e.what() << '\n';
        abort();
      }
#endif
      logfileptr->OFS()<<"Send Done"<<std::endl;

      double timeStop = get_time();
      if(this->iam==0 && this->options_.verbose){
        std::cout<<"Distribution time: "<<timeStop - timeSta<<std::endl;
      }
    }
    MPI_Barrier(pMat.comm);
  }




  template <typename T> inline symPACKMatrix<T>::symPACKMatrix():
    symPACKMatrixMeta<T>(){
    CommEnv_=NULL;
    //team_=nullptr;
    Mapping_ = NULL;
    Balancer_ = NULL;
    scheduler_ = NULL;
    scheduler2_ = NULL;
    //Local_=NULL;
    //Global_=NULL;
    //isGlobStructAllocated_ = false;
//    if(!upcxx::is_init()){
//      upcxx::init(NULL,NULL);
//    }

    if(logfileptr==NULL){
#ifdef NEW_UPCXX
      logfileptr = new LogFile(upcxx::rank_me(),false);
      logfileptr->OFS()<<"********* LOGFILE OF P"<<upcxx::rank_me()<<" *********"<<std::endl;
#else
      logfileptr = new LogFile(upcxx::myrank(),false);
      logfileptr->OFS()<<"********* LOGFILE OF P"<<upcxx::myrank()<<" *********"<<std::endl;
#endif
      logfileptr->OFS()<<"**********************************"<<std::endl;
    }

    logfileptr->OFS()<<"Shared node size "<<shmNode_.shmsize<<", rank "<<shmNode_.shmrank<<std::endl;


  }

  template <typename T> inline symPACKMatrix<T>::symPACKMatrix(DistSparseMatrix<T> & pMat, symPACKOptions & options ):symPACKMatrix(){
    Init(pMat, options);
  }

  template <typename T> inline symPACKMatrix<T>::~symPACKMatrix(){
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

    //if(team_!=nullptr){
    //  delete team_;
    //}
    //  if(Local_!=NULL){
    //    delete Local_;
    //  }
    //
    //  if(Global_!=NULL){
    //    delete Global_;
    //  }



  }


  //returns the 1-based index of supernode id global in the local supernode array
  template <typename T> inline Int symPACKMatrix<T>::snodeLocalIndex(Int global){
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
  template <typename T> inline SuperNode<T> * symPACKMatrix<T>::snodeLocal(Int global){
    Int iLocal = snodeLocalIndex(global);
    return LocalSupernodes_[iLocal -1];
  }

  template <typename T> 
    template< class Alloc>
    inline SuperNode<T,Alloc> * symPACKMatrix<T>::snodeLocal(Int global, std::vector<SuperNode<T,Alloc> *> & snodeColl){
      Int iLocal = snodeLocalIndex(global);
      return snodeColl[iLocal -1];
    }


  template <typename T> 
    inline void symPACKMatrixMeta<T>::findSupernodes(ETree& tree, Ordering & aOrder, std::vector<Int> & cc,std::vector<Int> & supMembership, std::vector<Int> & xsuper, Int maxSize ){
      SYMPACK_TIMER_START(FindSupernodes);
      Int size = this->iSize_;
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
    inline void symPACKMatrixMeta<T>::getLColRowCount(DistSparseMatrixGraph & dgraph, std::vector<Int> & cc, std::vector<Int> & rc){
      scope_timer(q,GetColRowCount_Classic);

      //The tree need to be postordered
      if(!this->ETree_.IsPostOrdered()){
        this->ETree_.PostOrderTree(this->Order_);
      }

      if(this->iam == 0 && (!dgraph.IsExpanded() ) ){
        throw std::logic_error( "DistSparseMatrixGraph must be expanded and permuted\n" );
      }

      dgraph.SetBaseval(1);
      dgraph.SetKeepDiag(1);


      int mpisize;
      MPI_Comm_size(dgraph.GetComm(),&mpisize);

      int mpirank;
      MPI_Comm_rank(dgraph.GetComm(),&mpirank);

      Int size = dgraph.size;

      MPI_Datatype Idxtype;
      MPI_Type_contiguous( sizeof(Idx), MPI_BYTE, &Idxtype );
      MPI_Type_commit(&Idxtype);

      MPI_Datatype Inttype;
      MPI_Type_contiguous( sizeof(Int), MPI_BYTE, &Inttype );
      MPI_Type_commit(&Inttype);

      {
        std::vector<Idx> level(size+1);
        std::vector<Idx> nchild(size+1,0);
        std::vector<Idx> fdesc(size+1);

        std::vector<Idx> storage;
        Int * weight = nullptr;
        Idx * set = nullptr;
        Idx * prvlf = nullptr;
        Idx * prvnbr = nullptr;
        Int * drc = nullptr;
        Idx * pxsup = nullptr;
        //Idx * fdesc = nullptr;

        if(mpirank==0){
          storage.resize(5*size+1+1);
          weight = (Int*)&storage[0];
          set = &storage[size+1];
          prvlf = &storage[2*size+1];
          prvnbr = &storage[3*size+1];
          drc = (Int*)&storage[4*size+1];
          //fdesc = &storage[5*size+1];
          pxsup = &storage[5*size+1];
          *pxsup = 1;
        }

        level[0] = 0;
        for(Idx k = size; k>=1; --k){
          if(mpirank==0){
            drc[k-1] = 1;
            set[k-1] = k;
            prvlf[k-1] = 0;
            prvnbr[k-1] = 0;
            weight[k] = 1;
          }
          fdesc[k] = k;
          level[k] = level[this->ETree_.PostParent(k-1)] + 1;
          nchild[k] = 0;
        }

        nchild[0] = 0;
        fdesc[0] = 0;
        for(Idx k =1; k<size; ++k){
          Idx parent = this->ETree_.PostParent(k-1);

          ++nchild[parent];
          if(mpirank==0){
            weight[parent] = 0;
          }

            Idx ifdesc = fdesc[k];
            if  ( ifdesc < fdesc[parent] ) {
              fdesc[parent] = ifdesc;
            }
        }




        Idx firstLocCol = dgraph.LocalFirstVertex()-1;

        if(mpirank>0){
          //If something is coming, we can resize          
          MPI_Probe(mpirank-1,mpirank-1,this->graph_.GetComm(),MPI_STATUS_IGNORE);

          storage.resize(5*size+1+1);
          weight = (Int*)&storage[0];
          set = &storage[size+1];
          prvlf = &storage[2*size+1];
          prvnbr = &storage[3*size+1];
          drc = (Int*)&storage[4*size+1];
          pxsup = &storage[5*size+1];
          MPI_Recv(storage.data(),storage.size(),Idxtype,mpirank-1,mpirank-1,this->graph_.GetComm(),MPI_STATUS_IGNORE);
        }

        //logfileptr->OFS()<<"weight: "; for(Int k = 0; k<=size; ++k){ logfileptr->OFS()<<weight[k]<<" "; } logfileptr->OFS()<<std::endl;
        //logfileptr->OFS()<<"set: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<set[k]<<" "; } logfileptr->OFS()<<std::endl;
        //logfileptr->OFS()<<"prvlf: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<prvlf[k]<<" "; } logfileptr->OFS()<<std::endl;
        //logfileptr->OFS()<<"prvnbr: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<prvnbr[k]<<" "; } logfileptr->OFS()<<std::endl;
        //logfileptr->OFS()<<"fdesc: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<fdesc[k]<<" "; } logfileptr->OFS()<<std::endl;
        //logfileptr->OFS()<<"level: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<level[k]<<" "; } logfileptr->OFS()<<std::endl;
        //logfileptr->OFS()<<"nchild: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<nchild[k]<<" "; } logfileptr->OFS()<<std::endl;
        //logfileptr->OFS()<<"drc: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<drc[k]<<" "; } logfileptr->OFS()<<std::endl;

        for(Idx loclownbr = 1; loclownbr<=dgraph.LocalVertexCount(); ++loclownbr){
          Idx lownbr = firstLocCol + loclownbr;

          //logfileptr->OFS()<<lownbr<<" weight: "; for(Int k = 0; k<=size; ++k){ logfileptr->OFS()<<weight[k]<<" "; } logfileptr->OFS()<<std::endl;
          //logfileptr->OFS()<<lownbr<<" set: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<set[k]<<" "; } logfileptr->OFS()<<std::endl;
          //logfileptr->OFS()<<lownbr<<" prvlf: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<prvlf[k]<<" "; } logfileptr->OFS()<<std::endl;
          //logfileptr->OFS()<<lownbr<<" prvnbr: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<prvnbr[k]<<" "; } logfileptr->OFS()<<std::endl;
          //logfileptr->OFS()<<lownbr<<" rc: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<drc[k]<<" "; } logfileptr->OFS()<<std::endl;

          Int lflag = 0;
          Idx ifdesc = fdesc[lownbr];
          Ptr jstrt = dgraph.colptr[loclownbr-1];
          Ptr jstop = dgraph.colptr[loclownbr] - 1;
          //           -----------------------------------------------
          //           for each ``high neighbor'', hinbr of lownbr ...
          //           -----------------------------------------------
          //std::vector<Idx> hack(&dgraph.rowind[jstrt-1],&dgraph.rowind[jstop-1]+1);
          //std::sort(hack.begin(),hack.end(),[this]( const Idx & a, const Idx &b){ return this->Order_.perm[a-1]<this->Order_.perm[b-1]; });
          //logfileptr->OFS()<<lownbr<<" adj: "<<hack<<std::endl;

          for(Ptr j = jstrt; j<=jstop;++j){
            Idx hinbr = dgraph.rowind[j-1];
            //Idx hinbr = hack[j-jstrt];
            //bassert(hinbr >= lownbr);
            if  ( hinbr > lownbr )  {
              //logfileptr->OFS()<<hinbr<<" vs "<<lownbr<<std::endl;
              if  ( ifdesc > prvnbr[hinbr-1] ) {
                //logfileptr->OFS()<<ifdesc<<" vs "<<prvnbr[hinbr-1]<<std::endl;
                //                       -------------------------
                //                       increment weight[lownbr].
                //                       -------------------------
                //logfileptr->OFS()<<"lownbr "<<lownbr<<std::endl;
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
                  drc[hinbr-1] += level[lownbr] - level[hinbr];
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
                  drc[hinbr-1] += level[lownbr] - level[lca];
                  --weight[lca];
                  //logfileptr->OFS()<<"lca "<<lca<<std::endl;
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
          Idx parent = this->ETree_.PostParent(lownbr-1);
          //logfileptr->OFS()<<"parent "<<parent<<std::endl;
          --weight[parent];

          if  ( lflag == 1  || nchild[lownbr] >= 2 ) {
            *pxsup = lownbr;
          }
          //logfileptr->OFS()<<"xsup "<<*pxsup<<std::endl;
          set[*pxsup-1] = parent;
        }

        if(mpirank<mpisize-1){
          MPI_Send(storage.data(),storage.size(),Idxtype,mpirank+1,mpirank,this->graph_.GetComm());
        }
        else{
          //logfileptr->OFS()<<"weight: "; for(Int k = 0; k<=size; ++k){ logfileptr->OFS()<<weight[k]<<" "; } logfileptr->OFS()<<std::endl;

          cc.assign(size,0);
          for(Int k = 1; k<=size; ++k){
            Int temp = cc[k-1] + weight[k];
            cc[k-1] = temp;
            Int parent = this->ETree_.PostParent(k-1);
            if  ( parent != 0 ) {
              cc[parent-1] += temp;
            }
          }

          rc.resize(size);
          for(size_t i=0;i<size;i++){ rc[i] = drc[i];}
        }
      }

      //Broadcast to everyone

      if (mpirank<mpisize-1){
        cc.resize(size);
      }

      MPI_Bcast(&cc[0],size,Inttype,mpisize-1,this->fullcomm_);

      rc.resize(size);
      MPI_Bcast(&rc[0],size,Inttype,mpisize-1,this->fullcomm_);



      MPI_Type_free(&Inttype);
      MPI_Type_free(&Idxtype);

    }




  template <typename T> 
    inline void symPACKMatrixMeta<T>::getLColRowCount(SparseMatrixGraph & sgraph, std::vector<Int> & cc, std::vector<Int> & rc){
      scope_timer(q,GetColRowCount_Classic);
      //The tree need to be postordered
      if(!this->ETree_.IsPostOrdered()){
        this->ETree_.PostOrderTree(this->Order_);
      }

      if(this->iam == 0 && (!sgraph.IsExpanded() ) ){
        throw std::logic_error( "SparseMatrixGraph must be expanded\n" );
      }

      sgraph.SetBaseval(1);
      sgraph.SetKeepDiag(1);
      //TODO EXPERIMENTAL
      //sgraph.SortEdges();

      MPI_Datatype Inttype;
      MPI_Type_contiguous( sizeof(Int), MPI_BYTE, &Inttype );
      MPI_Type_commit(&Inttype);

      Int size = sgraph.size;

      if(this->iam==0){
        cc.resize(size);
        rc.resize(size);
        std::vector<Idx> level(size+1);
        std::vector<Int> weight(size+1,0);
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
          level[k] = level[this->ETree_.PostParent(k-1)] + 1;
          weight[k] = 1;
          fdesc[k] = k;
          nchild[k] = 0;
        }

        nchild[0] = 0;
        fdesc[0] = 0;
        for(Idx k =1; k<size; ++k){
          Idx parent = this->ETree_.PostParent(k-1);
          weight[parent] = 0;
          ++nchild[parent];
          Idx ifdesc = fdesc[k];
          if  ( ifdesc < fdesc[parent] ) {
            fdesc[parent] = ifdesc;
          }
        }


          //logfileptr->OFS()<<"weight: "; for(Int k = 0; k<=size; ++k){ logfileptr->OFS()<<weight[k]<<" "; } logfileptr->OFS()<<std::endl;
          //logfileptr->OFS()<<"set: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<set[k]<<" "; } logfileptr->OFS()<<std::endl;
          //logfileptr->OFS()<<"prvlf: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<prvlf[k]<<" "; } logfileptr->OFS()<<std::endl;
          //logfileptr->OFS()<<"prvnbr: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<prvnbr[k]<<" "; } logfileptr->OFS()<<std::endl;
          //logfileptr->OFS()<<"fdesc: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<fdesc[k]<<" "; } logfileptr->OFS()<<std::endl;
          //logfileptr->OFS()<<"level: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<level[k]<<" "; } logfileptr->OFS()<<std::endl;
          //logfileptr->OFS()<<"nchild: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<nchild[k]<<" "; } logfileptr->OFS()<<std::endl;
          //logfileptr->OFS()<<"rc: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<rc[k]<<" "; } logfileptr->OFS()<<std::endl;




        for(Idx lownbr = 1; lownbr<=size; ++lownbr){

          //logfileptr->OFS()<<lownbr<<" weight: "; for(Int k = 0; k<=size; ++k){ logfileptr->OFS()<<weight[k]<<" "; } logfileptr->OFS()<<std::endl;
          //logfileptr->OFS()<<lownbr<<" set: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<set[k]<<" "; } logfileptr->OFS()<<std::endl;
          //logfileptr->OFS()<<lownbr<<" prvlf: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<prvlf[k]<<" "; } logfileptr->OFS()<<std::endl;
          //logfileptr->OFS()<<lownbr<<" prvnbr: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<prvnbr[k]<<" "; } logfileptr->OFS()<<std::endl;
          //logfileptr->OFS()<<lownbr<<" rc: "; for(Int k = 0; k<size; ++k){ logfileptr->OFS()<<rc[k]<<" "; } logfileptr->OFS()<<std::endl;

          Int lflag = 0;
          Idx ifdesc = fdesc[lownbr];
          Idx oldnbr = this->Order_.perm[lownbr-1];
          Ptr jstrt = sgraph.colptr[oldnbr-1];
          Ptr jstop = sgraph.colptr[oldnbr] - 1;


          //std::vector<Idx> hack(&sgraph.rowind[jstrt-1],&sgraph.rowind[jstop-1]+1);
          //std::for_each(hack.begin(),hack.end(),[this]( Idx & a){ a = this->Order_.invp[a-1]; });
          //logfileptr->OFS()<<lownbr<<" adj: "<<hack<<std::endl;


          //           -----------------------------------------------
          //           for each ``high neighbor'', hinbr of lownbr ...
          //           -----------------------------------------------
          for(Ptr j = jstrt; j<=jstop;++j){
            Idx hinbr = sgraph.rowind[j-1];
            hinbr = this->Order_.invp[hinbr-1];
            if  ( hinbr > lownbr )  {
              //logfileptr->OFS()<<hinbr<<" vs "<<lownbr<<std::endl;
              if  ( ifdesc > prvnbr[hinbr-1] ) {
                //logfileptr->OFS()<<ifdesc<<" vs "<<prvnbr[hinbr-1]<<std::endl;
                //                       -------------------------
                //                       increment weight[lownbr].
                //                       -------------------------
                //logfileptr->OFS()<<"lownbr "<<lownbr<<std::endl;
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
                  //logfileptr->OFS()<<"lca "<<lca<<std::endl;
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
          Idx parent = this->ETree_.PostParent(lownbr-1);
          //logfileptr->OFS()<<"parent "<<parent<<std::endl;
          --weight[parent];

          if  ( lflag == 1  || nchild[lownbr] >= 2 ) {
            xsup = lownbr;
          }
          //logfileptr->OFS()<<"xsup "<<xsup<<std::endl;
          set[xsup-1] = parent;
        }

          //logfileptr->OFS()<<"weight: "; for(Int k = 0; k<=size; ++k){ logfileptr->OFS()<<weight[k]<<" "; } logfileptr->OFS()<<std::endl;

        for(Int k = 1; k<=size; ++k){
          Int temp = cc[k-1] + weight[k];
          cc[k-1] = temp;
          Int parent = this->ETree_.PostParent(k-1);
          if  ( parent != 0 ) {
            cc[parent-1] += temp;
          }
        }
      }

      if(this->iam!=0){
        cc.resize(size);
        rc.resize(size);
      }
      //Broadcast to everyone 
      MPI_Bcast(&cc[0],size,Inttype,0,this->fullcomm_);
      MPI_Bcast(&rc[0],size,Inttype,0,this->fullcomm_);

      MPI_Type_free(&Inttype);
    }





  template <typename T> 
    inline void symPACKMatrixMeta<T>::relaxSupernodes(ETree& tree, std::vector<Int> & cc,std::vector<Int> & supMembership, std::vector<Int> & xsuper, RelaxationParameters & params ){
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

            //TODO rename nrelax0 as maxSnodeSize
            //TODO rename nrelax1 as ???
            //TODO rename nrelax2 as ???
            //TODO rename zrelax0 as ???
            //TODO rename zrelax1 as ???
            //TODO rename zrelax2 as ???
            merge = false;

            //Supernode is extremely small > merge it
            //TODO TRY THIS
            if( (merged_snode_size <=params.maxSize || params.maxSize==0 ) && params.nrelax0>0){
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
    inline void symPACKMatrixMeta<T>::symbolicFactorizationRelaxedDist(std::vector<Int> & cc){
      scope_timer(a,SymbolicFactorization);
      Int size = this->iSize_;
      ETree& tree = this->ETree_;
      Ordering & aOrder = this->Order_;
      DistSparseMatrixGraph & graph = this->graph_;
      std::vector<Int> & xsuper = this->Xsuper_;
      std::vector<Int> & SupMembership = this->SupMembership_;
      PtrVec & xlindx = this->locXlindx_;
      IdxVec & lindx = this->locLindx_;
      MPI_Comm & comm = this->graph_.comm;//CommEnv_->MPI_GetComm();


      //permute the graph
#ifdef EXPLICIT_PERMUTE
      DistSparseMatrixGraph & pGraph = graph;
#else
      DistSparseMatrixGraph pGraph = graph;
#endif

#if 0
      {
        double tstart = get_time();
        pGraph.Permute(&this->Order_.invp[0]);
        double tstop = get_time();
        logfileptr->OFS()<<"Permute time: "<<tstop-tstart<<std::endl;
      }



      //recompute xsuper by splitting some supernodes
      {

        //  gdb_lock();

        Idx colPerProc = size / this->np;
        Int numSplits = 0;
        for(Int snode = 1; snode<=xsuper.size()-1;snode++){
          Idx fstcol = xsuper[snode-1];
          Idx lstcol = xsuper[snode]-1;
          //check if these two columns are on the same processor
          Idx pOwnerFirst = std::min((Idx)this->np-1, (fstcol-1) / colPerProc);
          Idx pOwnerLast = std::min((Idx)this->np-1, (lstcol-1) / colPerProc);

          if(pOwnerFirst!=pOwnerLast){
            numSplits += pOwnerLast-pOwnerFirst;
          }
        }

        std::vector<Int> newSnodes(this->np,0);
        std::vector<Int> newXsuper(xsuper.size()+numSplits);
        Int pos = 0;
        for(Int snode = 1; snode<=xsuper.size()-1;snode++){
          Idx fstcol = xsuper[snode-1];
          Idx lstcol = xsuper[snode]-1;
          //check if these two columns are on the same processor
          Idx pOwnerFirst = std::min((Idx)this->np-1, (fstcol-1) / colPerProc);
          Idx pOwnerLast = std::min((Idx)this->np-1, (lstcol-1) / colPerProc);

          newXsuper[pos++] = xsuper[snode-1];
          if(pOwnerFirst!=pOwnerLast){
            for(Idx p = pOwnerFirst; p<pOwnerLast;p++){
              Idx curLstcol = (p+1)*colPerProc+1;//1-based



              //        assert(this->ETree_.PostParent(curLstcol-2)==curLstcol);
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



        //recompute this->XsuperDist_
        for(int p =1; p<this->np; p++){
          this->XsuperDist_[p]= this->XsuperDist_[p-1] + supPerProc + newSnodes[p-1];
        }
        this->XsuperDist_[this->np] = this->Xsuper_.size();


        //update cc by recomputing the merged structure



      }
#else
      {

#ifndef EXPLICIT_PERMUTE
        {
          double tstart = get_time();

          //recompute vertexDist based on XsuperDist
          //      std::vector<Idx> newVertexDist(this->np+1);
          //      for(int p = 0; p < this->np; p++){
          //        newVertexDist[p] = this->Xsuper_[this->XsuperDist_[p]-1] + (pGraph.GetBaseval()-1);
          //      }
          //      newVertexDist[this->np] = pGraph.size+1;

          std::vector<Idx> newVertexDist;
          {
            newVertexDist.resize(this->all_np+1,0);
            newVertexDist[this->all_np] = pGraph.size+1;
            for(int p = 0; p < this->all_np; p++){
              Int S = this->XsuperDist_[p];
              newVertexDist[p] = this->Xsuper_[S-1];
            }
          }



          pGraph.Permute(&this->Order_.invp[0],&newVertexDist[0]);
          double tstop = get_time();


          logfileptr->OFS()<<"Permute time: "<<tstop-tstart<<std::endl;
        }
#endif
      }
#endif

      //logfileptr->OFS()<<this->Xsuper_<<std::endl;

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

      Int nsuperLocal = this->XsuperDist_[this->iam+1]-this->XsuperDist_[this->iam];
      Int firstSnode = this->XsuperDist_[this->iam];
      Int lastSnode = firstSnode + nsuperLocal-1;
      xlindx.resize(nsuperLocal+1);


      Int firstColumn = this->Xsuper_[firstSnode-1];
      Int numColumnsLocal = this->Xsuper_[lastSnode] - firstColumn;


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

      if(this->iam>0){
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
              Int psrc = 0; for(psrc = 0; psrc<this->iam;psrc++){ if(this->XsuperDist_[psrc]<=jsup && jsup<this->XsuperDist_[psrc+1]){ break; } }
              //logfileptr->OFS()<<"trying to recv "<<jsup<<" max "<<recvLindx.size()*sizeof(Idx)<<" bytes"<<" from P"<<psrc<<std::endl;
              MPI_Recv(&recvLindx[0],recvLindx.size()*sizeof(Idx),MPI_BYTE,psrc,jsup+this->all_np,comm,&status);
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


          //TODO do better than this:need to ainline void sending unnecessary data
          //receive the speculative sends
          jsup = mrglnk[ksup-1];
          //get the next element of the list
          jsup = mrglnk[jsup-1];
          while(jsup>0){
            Int psrc = 0; for(psrc = 0; psrc<this->iam;psrc++){ if(this->XsuperDist_[psrc]<=jsup && jsup<this->XsuperDist_[psrc+1]){ break; } }
            if(psrc!=this->iam){
              MPI_Status status;
              MPI_Request request;
              MPI_Irecv(&recvLindx[0],recvLindx.size()*sizeof(Idx),MPI_BYTE,psrc,jsup+this->all_np,comm,&request);
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

              if(this->options_.relax.nrelax0==0 && this->options_.order_refinement_str == "NO") {
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
              if(this->options_.relax.nrelax0==0 && this->options_.order_refinement_str == "NO") {
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
          Int pdest = 0; for(pdest = 0; pdest<this->all_np;pdest++){ if(this->XsuperDist_[pdest]<=psup && psup<this->XsuperDist_[pdest+1]){ break; } }
          //if remote
          if(pdest!=this->iam){
            mpirequests.push_back(MPI_REQUEST_NULL);
            MPI_Request & request = mpirequests.back();
            MPI_Isend(&lindx[xlindx[locksup-1]-1],length*sizeof(Idx),MPI_BYTE,pdest,ksup+this->all_np,comm,&request);
          }
        }
      }

      bassert(nzend==0 || this->iam<this->np);
      lindx.resize(nzend);
      MPI_Barrier(comm);
    }




  template <typename T> 
    inline void symPACKMatrixMeta<T>::gatherLStructure(std::vector<Ptr>& xlindx, std::vector<Idx> & lindx){
      //Gather this->locXlindx_ and this->locLindx_
      //get other proc vertex counts
      Idx localVertexCnt = this->locXlindx_.size()-1;
      std::vector<Idx> remoteVertexCnt(this->np,0);
      MPI_Gather(&localVertexCnt,sizeof(localVertexCnt),MPI_BYTE,&remoteVertexCnt[0],sizeof(localVertexCnt),MPI_BYTE,0,this->workcomm_);
      Idx totalVertexCnt = std::accumulate(remoteVertexCnt.begin(),remoteVertexCnt.end(),0,std::plus<Idx>());
      if(this->iam==0){
        xlindx.resize(totalVertexCnt+1);
      }
      //compute receive displacements

    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(Ptr), MPI_BYTE, &type );
    MPI_Type_commit(&type);
    MPI_Datatype typeIdx;
    MPI_Type_contiguous( sizeof(Idx), MPI_BYTE, &typeIdx );
    MPI_Type_commit(&typeIdx);

      std::vector<int> rsizes(this->np,0);
      for(int p = 0; p<this->np;p++){rsizes[p] = (int)remoteVertexCnt[p];}
      std::vector<int> rdispls(this->np+1,0);
      rdispls[0]=0;
      std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
      MPI_Gatherv(&this->locXlindx_[0],localVertexCnt,type,&xlindx[0],&rsizes[0],&rdispls[0],type,0,this->workcomm_);


      Ptr localEdgeCnt = this->locLindx_.size();
      std::vector<Ptr> remoteEdgeCnt(this->np,0);
      MPI_Gather(&localEdgeCnt,sizeof(localEdgeCnt),MPI_BYTE,&remoteEdgeCnt[0],sizeof(localEdgeCnt),MPI_BYTE,0,this->workcomm_);
      Ptr totalEdgeCnt = std::accumulate(remoteEdgeCnt.begin(),remoteEdgeCnt.end(),0,std::plus<Ptr>());

      //fix xlindx
      if(this->iam==0){

        Idx pos = remoteVertexCnt[0];
        Ptr offset = 0;
        for(int p=1;p<this->np;p++){
          offset+=remoteEdgeCnt[p-1]; 
          for(Idx lv = 0; lv < remoteVertexCnt[p]; lv++){
            xlindx[pos++] += offset;//remoteEdgeCnt[p-1];//(vertexDist[p] - baseval);
          }
        }
        xlindx.back()=totalEdgeCnt + 1;

        lindx.resize(totalEdgeCnt);
      }

      //compute receive displacements
      rsizes.assign(this->np,0);
      for(int p = 0; p<this->np;p++){rsizes[p] = (int)remoteEdgeCnt[p];}
      rdispls[0]=0;
      std::partial_sum(rsizes.begin(),rsizes.end(),rdispls.begin()+1);
      MPI_Gatherv(&this->locLindx_[0],localEdgeCnt,typeIdx,&lindx[0],&rsizes[0],&rdispls[0],typeIdx,0,this->workcomm_);
    MPI_Type_free(&typeIdx);
    MPI_Type_free(&type);
    }




  template <typename T> 
    inline void symPACKMatrixMeta<T>::refineSupernodes(int ordflag,int altflag,DistSparseMatrix<T> * pMat){
      ETree& tree = this->ETree_;
      Ordering & aOrder = this->Order_;
      std::vector<Int> & supMembership = this->SupMembership_; 
      std::vector<Int> & xsuper = this->Xsuper_; 

      std::vector<int> ixlindx;
      std::vector<int> ilindx;

      std::vector<int>  new_invp;

      //Gather this->locXlindx_ and this->locLindx_
      {
        std::vector<Ptr> xlindx;
        std::vector<Idx> lindx;
        this->gatherLStructure(xlindx, lindx);

        if(this->iam==0){
          ixlindx.resize(xlindx.size());
          for(int i = 0;i<xlindx.size();i++){
            ixlindx[i] = xlindx[i];
          }
          ilindx.resize(lindx.size());
          for(int i = 0;i<lindx.size();i++){
            ilindx[i] = lindx[i];
          }
        }
      }

      if(this->iam==0){

        //logfileptr->OFS()<<"xlindx: "<<ixlindx<<std::endl;
        //logfileptr->OFS()<<"lindx: "<<ilindx<<std::endl;

        int neqns = this->iSize_;
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

        this->Order_.Compose(new_invp);

      }


#ifdef EXPLICIT_PERMUTE

      //Bcast the individual permutation
      new_invp.resize(this->iSize_);
      MPI_Bcast(&new_invp[0],this->iSize_*sizeof(int),MPI_BYTE,0,this->fullcomm_);

      //re permute the matrix
      pMat->Permute(&new_invp[0]);
      this->graph_ = pMat->GetLocalGraph();
      this->graph_.SetBaseval(1);
      this->graph_.SetKeepDiag(1);
      this->graph_.SetSorted(1);
#endif
      // broadcast invp
      Int N = aOrder.invp.size();
      MPI_Bcast(&aOrder.invp[0],N*sizeof(int),MPI_BYTE,0,this->fullcomm_);
      MPI_Bcast(&aOrder.perm[0],N*sizeof(int),MPI_BYTE,0,this->fullcomm_);
    }


  template <typename T> 
    template <class Allocator>
    inline SuperNode<T,Allocator> * symPACKMatrix<T>::CreateSuperNode(DecompositionType type,Int aiId, Int aiFr, Int aiFc, Int aiLc, Int ai_num_rows, Int aiN, Int aiNZBlkCnt, Int panel){
      SuperNode<T,Allocator> * retval = NULL;
        try{
      switch(type){
        case DecompositionType::LDL:
          retval = new SuperNodeInd<T,Allocator>( aiId,  aiFr, aiFc,  aiLc,  ai_num_rows,  aiN,  aiNZBlkCnt, panel);
          break;
        case DecompositionType::LL:
          retval = new SuperNode<T,Allocator>( aiId, aiFr,  aiFc,  aiLc,  ai_num_rows,  aiN,  aiNZBlkCnt, panel);
          break;
        default:
          retval = new SuperNode<T,Allocator>( aiId, aiFr,  aiFc,  aiLc,  ai_num_rows,  aiN,  aiNZBlkCnt, panel);
          break;
      }
        }
        catch(const MemoryAllocationException & e){
          throw;
        }
      return retval;
    }

  template <typename T> 
    template <class Allocator>
    inline SuperNode<T,Allocator> * symPACKMatrix<T>::CreateSuperNode(DecompositionType type,Int aiId, Int aiFr, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices, Int panel){
      SuperNode<T,Allocator> * retval = NULL;
      switch(type){
        case DecompositionType::LDL:
          retval = new SuperNodeInd<T,Allocator>( aiId, aiFr,  aiFc,  aiLc,  aiN,  rowIndices, panel);
          break;
        case DecompositionType::LL:
          retval = new SuperNode<T,Allocator>( aiId, aiFr,  aiFc,  aiLc,  aiN,  rowIndices, panel);
          break;
        default:
          retval = new SuperNode<T,Allocator>( aiId, aiFr,  aiFc,  aiLc,  aiN,  rowIndices, panel);
          break;
      }
      return retval;
    }

  template <typename T> 
    template <class Allocator>
    inline SuperNode<T,Allocator> * symPACKMatrix<T>::CreateSuperNode(DecompositionType type){
      SuperNode<T,Allocator> * retval = NULL;
        switch(type){
          case DecompositionType::LDL:
            retval = new SuperNodeInd<T,Allocator>();
            break;
          case DecompositionType::LL:
            retval = new SuperNode<T,Allocator>();
            break;
          default:
            retval = new SuperNode<T,Allocator>();
            break;
        }
      return retval;
    }


  template <typename T>
    template <class Allocator>
    inline SuperNode<T,Allocator> * symPACKMatrix<T>::CreateSuperNode(DecompositionType type,char * dataPtr,size_t size, Int firstRow){
      SuperNode<T,Allocator> * retval = NULL;
      switch(type){
        case DecompositionType::LDL:
          retval = new SuperNodeInd<T,Allocator>(dataPtr,size,firstRow);
          break;
        case DecompositionType::LL:
          retval = new SuperNode<T,Allocator>(dataPtr,size,firstRow);
          break;
        default:
          retval = new SuperNode<T,Allocator>(dataPtr,size,firstRow);
          break;
      }
      return retval;
    }





#include "symPACKMatrix_impl_FB_pull.hpp"
#include "symPACKMatrix_impl_solve.hpp"


#ifdef _INDEFINITE_
#undef _INDEFINITE_
#endif

  namespace TSP{

    inline void countBlock(const vector<Int> & Xsuper, const vector<Ptr> & Xlindx, const vector<Idx> & Lindx, const vector<Int> & supMembership, vector<int> & blkCount){
      blkCount.assign(Xlindx.size(),0);

      int totalBlocks = 0;
      for(Int I = 1; I<Xlindx.size();++I){
        //count number of contiguous blocks (at least one diagonal block)
        Idx fc = Xsuper[I-1];
        Idx lc = Xsuper[I]-1;
        Ptr fi = Xlindx[I-1];
        Ptr li = Xlindx[I]-1;
        Idx iPrevRow = Lindx[fi-1]-1;
        Idx iFirstRow = Lindx[fi-1];

        Int width = lc - fc + 1; 

        //only go through last column
        for(Idx col = lc; col<=lc;col++){ 
          //1 to count the diagonal block, 0 to skip it
          int nzBlockCnt = 1;
          Int prevFacing = I;
          for(Ptr idx = fi; idx<=li;idx++){
            Idx iRow = Lindx[idx-1];
            Int facing = supMembership[iRow-1];
            //enforce the first block to be a square diagonal block
            if(nzBlockCnt==1 && iRow>col){
              nzBlockCnt++;
              iFirstRow=iRow;
            }
            else if(iRow!=iPrevRow+1 || prevFacing != facing){
              nzBlockCnt++;
              iFirstRow=iRow;
            }
            iPrevRow=iRow;
            prevFacing = facing;
          }

          //totalBlocks+=nzBlockCnt;
          if(col==lc){
            blkCount[I-1] = nzBlockCnt;
            totalBlocks+=nzBlockCnt;
          }
        }
      }
      blkCount.back() = totalBlocks;
    }


    inline void symbolBuildRowtab(SymbolMatrix *symbptr) {
      SymbolCblk *cblk;
      SymbolBlok *blok;
      int *innbr, *intmp, *browtab;
      int  itercblk;
      int  cblknbr;
      int  edgenbr = symbptr->bloknbr - symbptr->cblknbr;

      cblknbr = symbptr->cblknbr;

      innbr = new int[cblknbr];
      std::fill(innbr,innbr+cblknbr,0);

      /* Count the number of ithis->nput edge per cblk */
      cblk = symbptr->cblktab;
      blok = symbptr->bloktab;
      for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
      {
#ifdef VERBOSE1
        cout<<itercblk<<"/"<<cblknbr<<std::endl;
        cout<<"["<<cblk[0].fcolnum<<".."<<cblk[0].lcolnum<<"] "<<cblk[0].bloknum<<std::endl;
        cout<<"["<<cblk[1].fcolnum<<".."<<cblk[1].lcolnum<<"] "<<cblk[1].bloknum<<std::endl;
#endif

        int iterblok = cblk[0].bloknum + 1;
        int lbloknum = cblk[1].bloknum;

#ifdef VERBOSE1
        cout<<"Looping from "<<iterblok<<" to "<<lbloknum<<std::endl;
#endif

        /* Skip diagonal block */
        blok++;

        /* Off-diagonal blocks */
        for( ; iterblok < lbloknum; iterblok++, blok++)
        {
#ifdef VERBOSE1
          cout<<iterblok<<" ["<< blok->frownum <<".."<<blok->lrownum<< "] facing "<<blok->fcblknm<<std::endl;
#endif
          innbr[ blok->fcblknm ]++;
        }
      }

#ifdef VERBOSE1
      for(int i =0;i<cblknbr;++i){
        cout<<i+1<<": "<<innbr[i]<<std::endl;
      }
#endif

      /* Initialize the brownum fields */
      cblk = symbptr->cblktab;
      cblk->brownum = 0;
      for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
      {
        cblk[1].brownum = cblk[0].brownum + innbr[ itercblk ];
        innbr[itercblk] = cblk[0].brownum;
      }


#ifdef VERBOSE1
      for(int i =0;i<=cblknbr;++i){
        cout<<i+1<<": "<<symbptr->cblktab[i].brownum<<std::endl;
      }
#endif

      assert( cblk[0].brownum == edgenbr );

      /* Initialize the browtab */
      browtab = new int[edgenbr];

      cblk = symbptr->cblktab;
      blok = symbptr->bloktab;
      for(itercblk=0; itercblk<symbptr->cblknbr; itercblk++, cblk++)
      {
        int iterblok = cblk[0].bloknum + 1;
        int lbloknum = cblk[1].bloknum;

        /* Skip diagonal block */
        blok++;

        /* Off-diagonal blocks */
        for( ; iterblok < lbloknum; iterblok++, blok++)
        {
          intmp = innbr + blok->fcblknm;
          browtab[ *intmp ] = iterblok;
          (*intmp)++;
        }
      }
      //if (symbptr->browtab == NULL) {
      //   delete [] symbptr->browtab;
      //}
      symbptr->browtab = browtab;

      delete [] innbr;
      return;
    }



    inline Order * GetPastixOrder(SymbolMatrix * symbptr,const vector<Int> & xsuper, const ETree & etree, const Int * perm, const Int * iperm){

      Order * order = new Order();
      order->baseval = 0;
      order->vertnbr = symbptr->nodenbr;
      order->cblknbr = symbptr->cblknbr;

      order->permtab = new int[order->vertnbr];
      for(int i =0; i<order->vertnbr;++i){order->permtab[i] = iperm[i]-1;}

      order->peritab = new int[order->vertnbr];
      for(int i =0; i<order->vertnbr;++i){order->peritab[i] = perm[i]-1;}

      order->rangtab = new int[xsuper.size()];
      for(int i = 0;i<xsuper.size();++i){order->rangtab[i] = xsuper[i]-1;} 

      assert(xsuper.size()-1 == etree.Size());
      order->treetab = new int[etree.Size()];
      for(int i = 0;i<etree.Size();++i){order->treetab[i] = etree.Parent(i)-1;} 

      return order;
    };


    inline SymbolMatrix *  GetPastixSymbolMatrix(const vector<Int> & xsuper,const vector<Int> & supMembership, vector<Ptr> & xlindx, vector<Idx> & lindx){
      //count the number of blocks and blocks per supernode
      vector<int> blkCount;
      countBlock(xsuper, xlindx, lindx, supMembership, blkCount);


      SymbolMatrix * symbmtx = new SymbolMatrix;

      symbmtx->baseval = 0;
      //ignore dof
      symbmtx->cblknbr = xsuper.size()-1;
      symbmtx->bloknbr = blkCount.back();
      //nodenbr should be n
      symbmtx->nodenbr = xsuper.back()-1;





      symbmtx->cblktab = new SymbolCblk[symbmtx->cblknbr+1];
      symbmtx->bloktab = new SymbolBlok[symbmtx->bloknbr];
      Int blockIdx = 1;
      for(Int I = 1; I<xlindx.size();++I){
        //count number of contiguous blocks (at least one diagonal block)

        Int fc = xsuper[I-1];
        Int lc = xsuper[I]-1;
        Ptr fi = xlindx[I-1];
        Ptr li = xlindx[I]-1;
        Idx iPrevRow = lindx[fi-1]-1;
        Idx iFirstRow = lindx[fi-1];

        Int width = lc - fc + 1; 

        //only go through last column
        for(Int col = lc; col<=lc;col++){ 
          //1 to count the diagonal block, 0 to skip it
          Int nzBlockCnt = 1;
          Int prevFacing = I;
          for(Ptr idx = fi; idx<=li;idx++){
            Idx iRow = lindx[idx-1];
            //enforce the first block to be a square diagonal block
            Int facing = supMembership[iRow-1];
            if(nzBlockCnt==1 && iRow>col){

              symbmtx->bloktab[blockIdx-1].frownum = fc-1;  /*< First row index            */
              symbmtx->bloktab[blockIdx-1].lrownum = lc-1;  /*< Last row index (inclusive) */
              symbmtx->bloktab[blockIdx-1].lcblknm = I-1;  /*< Local column block         */
              symbmtx->bloktab[blockIdx-1].fcblknm = I-1;  /*< Facing column block        */ //>> CBLK que tu update (ou cblk couran pr le bloc diagonal)

              symbmtx->cblktab[I-1].fcolnum = fc-1;
              symbmtx->cblktab[I-1].lcolnum = lc-1;
              symbmtx->cblktab[I-1].bloknum = blockIdx-1;  /*< First block in column (diagonal) */

              blockIdx++;
              nzBlockCnt++;
              iFirstRow=iRow;
            }
            else if(iRow!=iPrevRow+1 || prevFacing != facing){

              Int facingSnode = supMembership[iFirstRow-1];
              symbmtx->bloktab[blockIdx-1].frownum = iFirstRow-1;  /*< First row index            */
              symbmtx->bloktab[blockIdx-1].lrownum = iPrevRow-1;  /*< Last row index (inclusive) */
              symbmtx->bloktab[blockIdx-1].lcblknm = I-1;  /*< Local column block         */
              symbmtx->bloktab[blockIdx-1].fcblknm = facingSnode-1;  /*< Facing column block        */ //>> CBLK que tu update (ou cblk couran pr le bloc diagonal)

              blockIdx++;
              nzBlockCnt++;
              iFirstRow=iRow;
            }

            iPrevRow=iRow;
            prevFacing = facing;
          }

          if(col==lc){

            if(nzBlockCnt==1){
              symbmtx->bloktab[blockIdx-1].frownum = fc-1;  /*< First row index            */
              symbmtx->bloktab[blockIdx-1].lrownum = lc-1;  /*< Last row index (inclusive) */
              symbmtx->bloktab[blockIdx-1].lcblknm = I-1;  /*< Local column block         */
              symbmtx->bloktab[blockIdx-1].fcblknm = I-1;  /*< Facing column block        */ //>> CBLK que tu update (ou cblk couran pr le bloc diagonal)

              symbmtx->cblktab[I-1].fcolnum = fc-1;
              symbmtx->cblktab[I-1].lcolnum = lc-1;
              symbmtx->cblktab[I-1].bloknum = blockIdx-1;  /*< First block in column (diagonal) */
            }
            else{

              Int facingSnode = supMembership[iFirstRow-1];
              symbmtx->bloktab[blockIdx-1].frownum = iFirstRow-1;  /*< First row index            */
              symbmtx->bloktab[blockIdx-1].lrownum = iPrevRow-1;  /*< Last row index (inclusive) */
              symbmtx->bloktab[blockIdx-1].lcblknm = I-1;  /*< Local column block         */
              symbmtx->bloktab[blockIdx-1].fcblknm = facingSnode-1;  /*< Facing column block        */ //>> CBLK que tu update (ou cblk couran pr le bloc diagonal)

            }

            blockIdx++;
          }

        }
      }

      symbmtx->cblktab[xlindx.size()-1].fcolnum = xsuper.back()-1;
      symbmtx->cblktab[xlindx.size()-1].lcolnum = xsuper.back()-1;
      symbmtx->cblktab[xlindx.size()-1].bloknum = blockIdx-1;  /*< First block in column (diagonal) */

      symbolBuildRowtab(symbmtx);

      return symbmtx;
    }




    ///**
    // *
    // * @file symbol_reordering.c
    // *
    // *  PaStiX symbol structure routines
    // *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
    // *  LaBRI, University of Bordeaux 1 and IPB.
    // *
    // * @version 5.1.0
    // * @author Gregoire Pichon
    // * @author Mathieu Faverge
    // * @author Pierre Ramet
    // * @date 2015-04
    // *
    // **/
    //#include "common.h"
    //#include "symbol.h"
    //#include "order.h"

    static inline int
      compute_cblklevel( const int *treetab,
          int *levels,
          int  cblknum )
      {
        /* cblknum level has already been computed */
        if ( levels[cblknum] != 0 ) {
          return levels[cblknum];
        }
        else {
          int father = treetab[cblknum];

          if ( father == -1 ) {
            return 1;
          }
          else {
            return compute_cblklevel( treetab, levels, father ) + 1;
          }
        }
      }

    static inline int
      hamming_distance(int **vectors,
          int  *vectors_size,
          int   xi,
          int   xj,
          int   stop)
      {
        /* For the fictive vertex */
        if (xi == -1){
          return vectors_size[xj];
        }
        if (xj == -1){
          return vectors_size[xi];
        }

        int sum = 0;
        int *set1 = vectors[xi];
        int *set2 = vectors[xj];
        int *end1 = vectors[xi] + vectors_size[xi];
        int *end2 = vectors[xj] + vectors_size[xj];

        if (vectors_size[xi] - vectors_size[xj] >= stop){
          return stop;
        }
        if (vectors_size[xj] - vectors_size[xi] >= stop){
          return stop;
        }

        while((set1 < end1) && (set2 < end2))
        {
          if( *set1 == *set2)
          {
            set1++;
            set2++;
          }
          else if( *set1 < *set2 )
          {
            while (( set1 < end1 ) && ( *set1 < *set2 ))
            {
              sum ++;
              set1++;
            }
          }
          else if( *set1 > *set2 )
          {
            while (( set2 < end2 ) && ( *set1 > *set2 ))
            {
              sum ++;
              set2++;
            }
          }
          else
          {
            assert(0);
          }

          /* The computation is stopped if sum overlapped a given limit */
          if (sum >= stop){
            return stop;
          }
        }

        sum += end1-set1;
        sum += end2-set2;

        if (sum >= stop){
          return stop;
        }

        return sum;
      }


    static inline void
      symbol_reorder_tsp(int size, Order *order, int sn_id,
          int **lw_vectors, int *lw_vectors_size,
          int **up_vectors, int *up_vectors_size,
          int stop_criteria, int stop_when_fitting)
      {

        if ( size < 3 ) {
          return;
        }

        int  i, j, k, l, elected;
        int *tmpinvp;
        int *tmplen;
        int distance;

        tmpinvp = (int*)malloc((size+1)*sizeof(int));
        tmplen = (int*)malloc((size+1)*sizeof(int));
        memset(tmplen, 0, (size+1)*sizeof(int));

        tmpinvp[0] = -1;
        tmpinvp[1] = 0;

        distance = hamming_distance(lw_vectors, lw_vectors_size, 0, -1, stop_criteria);

        tmplen[0] = distance;
        tmplen[1] = distance;

        int min_cut = -1;
        for(i=1; i<size; i++) {
          int first_pos;
          int last_pos;

          int lw_before_pos;
          int lw_after_pos;

          int up_before_pos;
          int up_after_pos;

          /* Start by adding the row in first position */
          lw_before_pos = hamming_distance(lw_vectors, lw_vectors_size, i, tmpinvp[0], stop_criteria);
          lw_after_pos  = hamming_distance(lw_vectors, lw_vectors_size, i, tmpinvp[1], stop_criteria);
          up_after_pos  = hamming_distance(up_vectors, up_vectors_size, i, tmpinvp[1], 1);

          int minl = lw_before_pos + lw_after_pos - tmplen[0];
          int mpos = 1;
          int min_cut = -1;

          for(j=1; j<i; j++ ){
            up_before_pos = up_after_pos;
            up_after_pos  = hamming_distance(up_vectors, up_vectors_size, i, tmpinvp[j+1], 1);

            if ( up_before_pos < 1 ||
                up_after_pos  < 1 )
            {

              /* If split was used previously, this first distance may not be already computed */
              if (lw_after_pos == -1)
                lw_before_pos = hamming_distance(lw_vectors, lw_vectors_size, i, tmpinvp[j], stop_criteria);
              else
                lw_before_pos = lw_after_pos;


              lw_after_pos = hamming_distance(lw_vectors, lw_vectors_size, i, tmpinvp[j+1], stop_criteria);

              l = lw_before_pos + lw_after_pos - tmplen[j];


              /* Minimize the cut between two lines, for the same TSP result */
              if ( l == minl ) {
                if (lw_before_pos < min_cut){
                  min_cut = lw_before_pos;
                  minl = l; mpos = j+1;
                }
                if (lw_after_pos < min_cut){
                  min_cut = lw_after_pos;
                  minl = l; mpos = j+1;
                }
              }

              /* Position that minimizes TSP */
              if ( l < minl ) {
                minl = l; mpos = j+1;

                min_cut = lw_before_pos;
                if (lw_after_pos < min_cut){
                  min_cut = lw_after_pos;
                }
              }

              if ( l < minl ) {
                minl = l; mpos = j+1;
                min_cut = lw_before_pos;
                if (lw_after_pos < min_cut){
                  min_cut = lw_after_pos;
                }
              }


              /* Stop if two lines are equal (already done tmpinvp[j]) */
              if (lw_after_pos == 0){
                min_cut = 0;
                minl = l; mpos = j+1;
                j = i;
              }
            }
            else{
              lw_after_pos = -1;
            }

          }

          /* Test between last and first */
          first_pos = hamming_distance(lw_vectors, lw_vectors_size, i, tmpinvp[0], stop_criteria);
          last_pos  = hamming_distance(lw_vectors, lw_vectors_size, i, tmpinvp[i], stop_criteria);

          lw_before_pos = hamming_distance(lw_vectors, lw_vectors_size, i, tmpinvp[mpos-1], stop_criteria);
          lw_after_pos  = hamming_distance(lw_vectors, lw_vectors_size, i, tmpinvp[mpos  ], stop_criteria);

          l = first_pos + last_pos - tmplen[i];
          if ( l < minl ) {
            minl = l; mpos = i+1;
          }

          if (mpos > 0){
            tmplen[mpos-1] = lw_before_pos;
          }

          if (mpos < (i+1))
          {
            int tmpi, tmpl;
            k = i;
            l = lw_after_pos;

            /* Insert the line in the tmpinvp/tmplen arrays */
            for(j=mpos; j<i+2; j++ )
            {
              tmpi = tmpinvp[j];
              tmpl = tmplen[j];

              tmpinvp[j] = k;
              tmplen[j]  = l;

              k = tmpi;
              l = tmpl;
            }
          }
          else {
            tmpinvp[i+1] = i;
            tmplen[i+1]  = first_pos;
          }
        }

        elected = 0;
        for (i=0; i<size; i++)
        {
          if (tmpinvp[i] == -1){
            elected = i;
          }
        }

        int *sn_connected;
        sn_connected = (int*)malloc(size*sizeof(int));
        {
          //TODO
          int *peritab = order->peritab + order->rangtab[sn_id];
          for (i=0; i<size; i++)
          {
            sn_connected[i] = peritab[ tmpinvp[(i + 1 + elected)%(size+1)] ];
          }
          memcpy( peritab, sn_connected, size * sizeof(int) );
        }

        free(sn_connected);
        free(tmpinvp);
        free(tmplen);
      }

    static inline void
      symbol_reorder_cblk( const SymbolMatrix *symbptr,
          const SymbolCblk   *cblk,
          Order              *order,
          const int *levels,
          int        cblklvl,
          int       *depthweight,
          int        depthmax,
          int        split_level,
          int                 stop_criteria,
          int                 stop_when_fitting,
          double             *time_compute_vectors,
          double             *time_update_perm)
      {
        SymbolBlok *blok;
        int **up_vectors, *up_vectors_size;
        int **lw_vectors, *lw_vectors_size;
        int size = cblk->lcolnum - cblk->fcolnum + 1;
        int local_split_level = split_level;
        int i, iterblok;
        int *brow = symbptr->browtab;
        double timer;

        /**
         * Compute hamming vectors in two subsets:
         *   - The upper subset contains the cblk with level higher than the split_level
         *     in the elimination tree, (or depth lower than levels[cblk])
         *   - The lower subset contains the cblk with level lower than the split_level
         *     in the elimination tree, (or depth higher than levels[cblk])
         *
         * The delimitation between the lower and upper levels is made such that
         * the upper level represents 17% to 25% of the total number of cblk.
         */
        //  clockStart(timer);
        {
          int weight = 0;

          /* Compute the weigth of each level */
          //MATHIAS: this is a loop through blocks in current column 
          for(iterblok=cblk[0].brownum; iterblok<cblk[1].brownum; iterblok++)
          {
            int blokweight;
            blok = symbptr->bloktab + brow[iterblok];
            blokweight = blok->lrownum - blok->frownum + 1;
            depthweight[ levels[ blok->lcblknm ] - 1 ] += blokweight;
            weight += blokweight;
          }

          /**
           * Compute the split_level:
           *    We start with the given split_level parameter
           *    and we try to correct it to minimize the following iterative process
           */
          {
            /* Current for each line within the current cblk the number of contributions */
            int up_total = 0;
            int lw_total = 0;
            int sign = 0;

split:
            up_total = 0;
            lw_total = 0;

            for(i=0; i<local_split_level; i++)
            {
              up_total += depthweight[i];
            }
            for(; i<depthmax; i++)
            {
              lw_total += depthweight[i];
            }

            /* If there are too many upper bloks */
            if ( (lw_total < (5 * up_total)) &&
                (lw_total > 10) && (up_total > 10) && (sign <= 0))
            {
              local_split_level--;
              sign--;
              goto split;
            }

            /* If there are too many lower bloks */
            if ( (lw_total > (3 * up_total)) &&
                (lw_total > 10) && (up_total > 10) && (sign >= 0) )
            {
              local_split_level++;
              sign++;
              goto split;
            }

            /* Adjust to depth of the level array */
            /* symbol_reorder_cblk( symbptr, cblk, order, */
            /*                      levels, levels[itercblk], */
            /*                      depthweight + levels[itercblk], maxdepth-levels[itercblk], */
            /*                      split_level, stop_criteria, stop_when_fitting, */
            /*                      &time_compute_vectors, &time_update_perm); */
            /* local_split_level += cblklvl; */
            /* for(i=0; (i<local_split_level) && (depthweight[i] != 0); i++) */
            /* for(; (i<depthmax) && (depthweight[i] != 0); i++) */
          }

          /* Compute the Hamming vector size for each row of the cblk */
          up_vectors_size = (int*)malloc(size*sizeof(int));
          memset(up_vectors_size, 0, size * sizeof(int));
          lw_vectors_size = (int*)malloc(size*sizeof(int));
          memset(lw_vectors_size, 0, size * sizeof(int));

          for(iterblok=cblk[0].brownum; iterblok<cblk[1].brownum; iterblok++)
          {
            blok = symbptr->bloktab + brow[iterblok];

            /* For upper levels in nested dissection */
            if (levels[blok->lcblknm] <= local_split_level){
              for (i=blok->frownum; i<=blok->lrownum; i++){
                int index = i - cblk->fcolnum;
                up_vectors_size[index]++;
              }
            }
            else{
              for (i=blok->frownum; i<=blok->lrownum; i++){
                int index = i - cblk->fcolnum;
                lw_vectors_size[index]++;
              }
            }
          }

          /* Initiate Hamming vectors structure */
          lw_vectors = (int**)malloc(size*sizeof(int*));
          up_vectors = (int**)malloc(size*sizeof(int*));
          for (i=0; i<size; i++) {
            lw_vectors[i] = (int*)malloc(lw_vectors_size[i]*sizeof(int));
            up_vectors[i]=(int*)malloc(up_vectors_size[i]*sizeof(int));
            memset(lw_vectors[i], 0, lw_vectors_size[i] * sizeof(int));
            memset(up_vectors[i], 0, up_vectors_size[i] * sizeof(int));
          }
          memset(lw_vectors_size, 0, size * sizeof(int));
          memset(up_vectors_size, 0, size * sizeof(int));

          /* Fill-in vectors structure with contributing cblks */
          for(iterblok=cblk[0].brownum; iterblok<cblk[1].brownum; iterblok++)
          {
            blok = symbptr->bloktab + brow[iterblok];

            /* For upper levels in nested dissection */
            if (levels[blok->lcblknm] <= local_split_level) {
              for (i=blok->frownum; i<=blok->lrownum; i++){
                int index = i - cblk->fcolnum;
                up_vectors[index][up_vectors_size[index]] = blok->lcblknm;
                up_vectors_size[index]++;
              }
            }
            else{
              for (i=blok->frownum; i<=blok->lrownum; i++){
                int index = i - cblk->fcolnum;
                lw_vectors[index][lw_vectors_size[index]] = blok->lcblknm;
                lw_vectors_size[index]++;
              }
            }
          }
        }

        //  clockStop(timer);
        //  *time_compute_vectors += clockVal(timer);

        //  clockStart(timer);
        {
          /* Apply the pseudo-TSP algorithm to the rows in the current supernode */
          symbol_reorder_tsp(size, order, cblk - symbptr->cblktab,
              lw_vectors, lw_vectors_size,
              up_vectors, up_vectors_size,
              stop_criteria, stop_when_fitting);
        }
        //  clockStop(timer);
        //  *time_update_perm += clockVal(timer);

        for (i=0; i<size; i++){
          free(lw_vectors[i]);
          free(up_vectors[i]);
        }

        free(lw_vectors);
        free(up_vectors);
        free(lw_vectors_size);
        free(up_vectors_size);
      }

    /* For split_level parameter */
    /* The chosen level to reduce computational cost: no effects if set to 0 */
    /* A first comparison is computed according to upper levels */
    /* If hamming distances are equal, the computation goes through lower levels */

    /* For stop_criteria parameter */
    /* Criteria to limit the number of comparisons when computing hamming distances */

    /* For stop_when_fitting parameter */
    /* Criteria to insert a line when no extra-blok is created */
    /* If set to 0, the algorithm will minimize the cut between two lines */

    inline void
      symbolReordering( const SymbolMatrix *symbptr,
          Order *order,
          int split_level,
          int stop_criteria,
          int stop_when_fitting )
      {
        SymbolCblk  *cblk;
        int itercblk;
        int cblknbr = symbptr->cblknbr;

        double time_compute_vectors = 0.;
        double time_update_perm     = 0.;

        int i, maxdepth;
        int *levels, *depthweight;

        /* Create the level array to compute the depth of each cblk and the maximum depth */
        {
          maxdepth = 0;
          levels = new int[cblknbr];
          for(int i = 0; i<cblknbr; ++i){levels[i] = 0;}

          for (i=0; i<cblknbr; i++) {
            levels[i] = compute_cblklevel( order->treetab, levels, i );
            maxdepth = std::max( maxdepth, levels[i] );
          }
        }

        /**
         * Solves the Traveler Salesman Problem on each cblk to minimize the number
         * of off-diagonal blocks per row
         */
        cblk = symbptr->cblktab;
        depthweight = new int[maxdepth];
        for (itercblk=0; itercblk<cblknbr; itercblk++, cblk++) {

          memset( depthweight, 0, maxdepth * sizeof(int) );

          symbol_reorder_cblk( symbptr, cblk, order,
              levels, levels[itercblk],
              depthweight, maxdepth,
              split_level, stop_criteria, stop_when_fitting,
              &time_compute_vectors, &time_update_perm);
        }

        printf("Time to compute vectors  %lf s\n", time_compute_vectors);
        printf("Time to update  perm     %lf s\n", time_update_perm);

        /* Update the permutation */
        for (i=0; i<symbptr->nodenbr; i++) {
          order->permtab[ order->peritab[i] ] = i;
        }
        delete [] levels; levels = NULL;
        delete [] depthweight; depthweight = NULL;
      }


  } //namespace TSP



}


#endif //_SYMPACK_MATRIX_IMPL_HP_

