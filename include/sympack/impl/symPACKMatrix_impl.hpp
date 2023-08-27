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

#define SPLIT_AT_BOUNDARY

namespace symPACK{


  template<typename T> 
    inline void symPACKMatrix<T>::generateTaskGraph(taskGraph & graph,
      std::vector<Int> & AggregatesToRecv,  std::vector<Int>& LocalAggregates, std::vector<Int> & mw, std::vector<Int> & mh)
    {
      //we will need to communicate if only partial xlindx_, lindx_
      //idea: build tasklist per processor and then exchange
      //tuple is: src_snode,tgt_snode,src_first_row,op_type,
      using upd_tuple_t = std::tuple<Idx,Idx,Idx,Factorization::op_type>;
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

          Idx blockIdx  = (Idx)-1;
          Idx prevRow   = (Idx)-1;
          Idx prevSnode = (Idx)-1;

          for(Ptr sidx = lfi; sidx<=lli;sidx++){
            Idx row = this->locLindx_[sidx-1];
            J = this->SupMembership_[row-1];

            //Split at boundary or after diagonal block
            if(J!=prevSnode){

              Int iUpdater = Mapping_->Map(J-1,I-1);

              if(J>I){
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
        ssizes[itp->first] = itp->second.size();
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
      vector<upd_tuple_t> recvbuf(rdispls.back());
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
        type = std::get<3>(*it);

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
          default: break; // silence -Wswitch warnings from clang
        }


        max_mem_req = std::max(mem_req,max_mem_req);



        Task.getHash = [&](char * ameta)->GenericTask::id_type{
          char * pmeta = ameta;
          std::stringstream sstr;
          sstr<<meta[0]<<"_"<<meta[1]<<"_"<<0<<"_"<<(Int)(*reinterpret_cast<Factorization::op_type*>(&meta[3]));
          return hash_fn(sstr.str());
        };
        Task.id = Task.getHash(Task.meta.data());

        graph.addTask(pTask);


      }
      logfileptr->OFS()<<"Maximum single task memory requirement is: "<<max_mem_req<<std::endl;
    }









  template <typename T> 
   inline  void symPACKMatrix<T>::Factorize(){
    SYMPACK_TIMER_START(NUMERICAL_FACT);
    if(this->iam<this->np){
      switch(this->options_.factorization){
        case FANBOTH:
          FanBoth();
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
          symPACKOS<<"Total volume of communication: "<<totalVolComm<<std::endl;
          symPACKOS<<"Total number of messages: "<<totalNumMsg<<std::endl;
        }

        gVolComm=0;
        gNumMsg=0;
      }
#endif
    }
  }


  template <typename T> inline void symPACKMatrix<T>::Solve(T * RHS, int nrhs, int rhs_size,  T * Xptr) {
    scope_timer(a,SPARSE_SOLVE);
    if (this->options_.iterRefinement) {
      throw std::runtime_error("Iterative refinement is not yet implemented.");
      this->solve_(RHS,nrhs,Xptr);
    }
    else{
        Contributions2_.clear();
        switch(this->options_.factorization){
          case FANBOTH:
            this->solve_(RHS,nrhs,Xptr);
            break;
          case FANOUT:
            this->solve_(RHS,nrhs,Xptr);
            break;
          default:
            this->solve_(RHS,nrhs,Xptr);
            break;
        }
    }
  }

  template<typename T> inline void symPACKMatrix<T>::GetSolution(T * B, int nrhs){
    Int n = this->iSize_;

    {
      std::fill(B,B+n*nrhs,T(0.0));

      //Gather B from everybody and put it in the original matrix order
      std::vector<T> tmp_nzval;
      for(Int I=1; I<this->Xsuper_.size();++I){
        Int iOwner = this->Mapping_->Map(I-1,I-1);
        T * data;
        Int snode_size = this->Xsuper_[I] - this->Xsuper_[I-1];
        Int nzcnt = snode_size * nrhs;

        if( iOwner == this->iam ){
          Int Ilocal = snodeLocalIndex(I); 
          auto contrib = std::dynamic_pointer_cast< SuperNode<T,UpcxxAllocator> >(Contributions2_[Ilocal-1]);
          data = contrib->GetNZval(0);
        }
        else{
          tmp_nzval.resize(nzcnt);
          data = &tmp_nzval[0];
        }

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
            symPACKOS<<"P"<<this->iam<<" has sent update from Supernode "<<prev_src_snode->Id()<<" to Supernode "<<tgt_snode_id<<std::endl;
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
        symPACKOS<<"P"<<this->iam<<" has sent update from Supernode "<<prev_src_snode->Id()<<" to Supernode "<<tgt_snode_id<<std::endl;
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
      Idx prevSnode = (Idx)-1;
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
            }
            else{
              break;
            }
          } 

          if(row==iPrevRow+1){
            if( nextSup<this->Xsuper_.size()-0 && row==next_fc){
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
                logfileptr->OFS()<<std::scientific<<"("<<row+i<<","<<src_first_col+j<<") "<<ToMatlabScalar(val[i*src_snode->Size()+j])<<" "<<std::endl;
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
    if(CommEnv_!=nullptr){
      delete CommEnv_;
    }

    this->options_ = options;
    logfileptr->verbose = this->options_.verbose>0;
    if(this->options_.verbose==0){
      symPACKOS.rdbuf(nullptr);
    }

#ifndef NO_MPI
    if(this->fullcomm_!=MPI_COMM_NULL){
      MPI_Comm_free(&this->fullcomm_);
    }
    MPI_Comm_dup(this->options_.MPIcomm,&this->fullcomm_);
#endif


    if(this->workcomm_!=MPI_COMM_NULL){
      MPI_Comm_free(&this->workcomm_);
    }


    this->all_np = 0;
    MPI_Comm_size(this->options_.MPIcomm,&this->all_np);
    MPI_Comm_rank(this->options_.MPIcomm,&this->iam);

    this->iam = upcxx::rank_me();
    this->all_np = upcxx::rank_n();

    this->np = this->options_.used_procs(this->all_np);


    MPI_Comm_split(this->options_.MPIcomm,this->iam<this->np,this->iam,&this->workcomm_);
    if (this->workteam_!=nullptr)  { this->workteam_->destroy();}
    this->workteam_.reset(new upcxx::team(upcxx::world().split(this->iam<this->np,this->iam)));
    CommEnv_ = new CommEnvironment(this->workcomm_);

    //do another split to contain P0 and all the non working processors
    if(this->all_np!=this->np){
      this->non_workcomm_ = MPI_COMM_NULL;
      MPI_Comm_split(this->options_.MPIcomm,this->iam==0||this->iam>=this->np,this->iam==0?0:this->iam-this->np+1,&this->non_workcomm_);
    }

    if(this->iam<this->np){
      //Options
      if(this->options_.maxIrecv==-1 && this->options_.maxIsend==-1){
        gMaxIrecv = -1;
      }
      else{
        gMaxIrecv = this->options_.maxIrecv + this->options_.maxIsend;
      }

      
      switch(this->options_.scheduler){
        case DL:
          scheduler_new_ = std::make_shared<DLScheduler< std::shared_ptr<GenericTask> > >( );
          break;
        case MCT:
          //scheduler2_ = new MCTScheduler<FBTask>();
          break;
        case PR:
          //scheduler2_ = new PRScheduler<FBTask>();
          break;
        case FIFO:
          scheduler_new_ = std::make_shared<FIFOScheduler< std::shared_ptr<GenericTask> > >( );
          break;
        default:
          scheduler_new_ = std::make_shared<DLScheduler< std::shared_ptr<GenericTask> > >( );
          break;
      }
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

    this->iSize_ = pMat.size;
    this->graph_ = pMat.GetLocalGraph();
    this->graph_.SetBaseval(1);
    this->graph_.SetSorted(1);
    this->graph_.ExpandSymmetric();

    logfileptr->OFS()<<"Matrix structure expanded"<<std::endl;

    SparseMatrixGraph * sgraph = nullptr;
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
          if(this->options_.perm ==nullptr){
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
          symPACKOS<<"Ordering time: "<<timeStop - timeSta<<std::endl;
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

    if(sgraph==nullptr){ 
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

      auto backinvp = this->Order_.invp;

      //assume matrix is not permuted yet
      graph.Permute(this->Order_.invp.data());

      logfileptr->OFS()<<"graph permuted 1/2"<<std::endl;


      this->ETree_.ConstructETree(graph,this->Order_);
      this->ETree_.PostOrderTree(this->Order_);
      logfileptr->OFS()<<"ETree computed and postordered"<<std::endl;

      std::vector<Int> relinvp;
      this->Order_.GetRelativeInvp(backinvp,relinvp);

      graph.Permute(relinvp.data());


      logfileptr->OFS()<<"graph permuted 2/2"<<std::endl;

      {
        double timeStart = get_time();
        this->getLColRowCount(graph,cc,rc);
        double timeStop = get_time();
        if(this->iam==0 && this->options_.verbose){
          symPACKOS<<"Column count (distributed) construction time: "<<timeStop - timeStart<<std::endl;
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
        symPACKOS<<"Elimination tree construction time: "<<timeStop - timeSta<<std::endl;
      }
    }
    else
    {
      //gather the graph if necessary to build the elimination tree
      if(sgraph==nullptr){ 
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

        {
          double timeStart = get_time();
          this->getLColRowCount(*sgraph,cc,rc);
          double timeStop = get_time();
          if(this->iam==0 && this->options_.verbose){
            symPACKOS<<"Column count (gather + serial + bcast) construction time: "<<timeStop - timeStart<<std::endl;
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
          symPACKOS<<"Elimination tree construction time: "<<timeStop - timeSta<<std::endl;
        }
      }
    }

    //compute some statistics
    if(this->iam==0 && this->options_.verbose){
      double flops = 0.0;
      int64_t NNZ = 0;
      for(Int i = 0; i<cc.size();++i){
        flops+= (double)pow((double)cc[i],2.0);
        NNZ+=cc[i];
      }
      symPACKOS<<"Flops: "<<flops<<std::endl;
      symPACKOS<<"NNZ in L factor: "<<NNZ<<std::endl;
    }


    //get rid of the sequential graph
    if(sgraph!=nullptr && this->options_.order_refinement_str.substr(0,4) != "TSPB"){delete sgraph;}

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

      this->workteam_->destroy();
      this->workteam_.reset(new upcxx::team(upcxx::world().split(this->iam<this->np,this->iam)));

      bassert(CommEnv_!=nullptr);
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
      if(this->Mapping_ != nullptr){
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
      else if(this->options_.mappingTypeStr ==  "FANOUT"){
        this->Mapping_ = new Row2D(pmapping, pmapping, pmapping, 1);
      }
      else if(this->options_.mappingTypeStr ==  "FANIN"){
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

      double timeStaSymb = get_time();
      this->symbolicFactorizationRelaxedDist(cc);

      double timeStopSymb = get_time();
      if(this->iam==0 && this->options_.verbose){
        symPACKOS<<"Symbolic factorization time: "<<timeStopSymb - timeStaSymb<<std::endl;
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
        else if(this->options_.order_refinement_str.substr(0,4) == "TSPB"){

          if(sgraph==nullptr){ 
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

            if(sgraph!=nullptr){delete sgraph;}

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
              symPACKOS<<"TSPB reordering done in "<<timeStop-timeSta<<std::endl;
            }


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
          symPACKOS<<"Supernode reordering done in "<<timeStop-timeSta<<std::endl;
        }

        {
          double timeSta = get_time();
          this->symbolicFactorizationRelaxedDist(cc);
          double timeStop = get_time();
          if(this->iam==0 && this->options_.verbose){
            symPACKOS<<"Symbolic factorization time: "<<timeStop - timeSta<<std::endl;
          }
          logfileptr->OFS()<<"Symbfact done"<<std::endl;
        }
      }

      double timeStop = get_time();
      if(this->iam==0 && this->options_.verbose){
        symPACKOS<<"Total symbolic factorization time: "<<timeStop - timeSta<<std::endl;
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

      if (this->Balancer_!=nullptr){
        delete this->Balancer_;
      }

      std::vector<Int> map;
      if(this->options_.load_balance_str=="SUBCUBE-FI"){
        if(this->iam==0 && this->options_.verbose){ symPACKOS<<"Subtree to subcube FI mapping used"<<std::endl;}
        ETree SupETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);
        this->Balancer_ = new SubtreeToSubcube(this->np,SupETree,this->Xsuper_,this->XsuperDist_,this->SupMembership_,this->locXlindx_,this->locLindx_,cc,this->fullcomm_,true);
      }
      else if(this->options_.load_balance_str=="SUBCUBE-FO"){
        if(this->iam==0 && this->options_.verbose){ symPACKOS<<"Subtree to subcube FO mapping used"<<std::endl;}
        ETree SupETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);
        this->Balancer_ = new SubtreeToSubcube(this->np,SupETree,this->Xsuper_,this->XsuperDist_,this->SupMembership_,this->locXlindx_,this->locLindx_,cc,this->fullcomm_,false);
      }
      else if(this->options_.load_balance_str=="SUBCUBE-VOLUME-FI"){
        if(this->iam==0 && this->options_.verbose){ symPACKOS<<"Subtree to subcube volume FI mapping used"<<std::endl;}
        ETree SupETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);
        this->Balancer_ = new SubtreeToSubcubeVolume(this->np,SupETree,this->Xsuper_,this->XsuperDist_,this->SupMembership_,this->locXlindx_,this->locLindx_,cc,this->fullcomm_,true);
      }
      else if(this->options_.load_balance_str=="SUBCUBE-VOLUME-FO"){
        if(this->iam==0 && this->options_.verbose){ symPACKOS<<"Subtree to subcube volume FO mapping used"<<std::endl;}
        ETree SupETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);
        this->Balancer_ = new SubtreeToSubcubeVolume(this->np,SupETree,this->Xsuper_,this->XsuperDist_,this->SupMembership_,this->locXlindx_,this->locLindx_,cc,this->fullcomm_,false);
      }
      else if(this->options_.load_balance_str=="NNZ"){
        if(this->iam==0 && this->options_.verbose){ symPACKOS<<"Load Balancing on NNZ used"<<std::endl;}
        this->Balancer_ = new NNZBalancer(this->np,this->Xsuper_,cc);
      }
      else if(this->options_.load_balance_str=="WORK"){
        if(this->iam==0 && this->options_.verbose){ symPACKOS<<"Load Balancing on WORK used"<<std::endl;}
        this->Balancer_ = new WorkBalancer(this->np,this->Xsuper_,cc);
      }


      if (this->Balancer_!=nullptr){
        map = this->Balancer_->GetMap();
        TreeLoadBalancer * testBalancer = dynamic_cast<TreeLoadBalancer*>(this->Balancer_);
        if(testBalancer==nullptr){
          this->Mapping_->Update(map);
        }
        else{

          TreeMapping * test = dynamic_cast<TreeMapping*>(this->Mapping_);
          if(test==nullptr){
            this->Mapping_->Update(map);
          }
          else{
            test->Update((TreeLoadBalancer*)this->Balancer_);
          }
        }
      }
      double timeStop = get_time();
      if(this->iam==0 && this->options_.verbose){
        symPACKOS<<"Load balancing time: "<<timeStop - timeSta<<std::endl;
      }
    }


    //Now split supernodes larger than the maximum size. XlindxLocal_, LindxLocal_, this->Xsuper_, this->SupMembership_, Balancer_ and Mapping_ need to up updated

#ifdef _OUTPUT_ETREE_
    logfileptr->OFS()<<"ETree is "<<this->ETree_<<std::endl;
    {
      auto supETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);
      logfileptr->OFS()<<"Supernodal ETree is "<<supETree<<std::endl;
      logfileptr->OFS()<<"Mapping is ";for(Int i = 0;i<supETree.Size();i++){logfileptr->OFS()<<this->Mapping_->Map(i,i)<<" ";}logfileptr->OFS()<<std::endl;
      logfileptr->OFS()<<"Xsuper is "<<this->Xsuper_<<std::endl;
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
        {
          generateTaskGraph(taskGraph_, AggregatesToRecv, LocalAggregates,UpdateWidth_,UpdateHeight_);
        }

        double timeStop = get_time();
        if(this->iam==0 && this->options_.verbose){
          symPACKOS<<"Task graph generation time: "<<timeStop - timeSta<<std::endl;
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
        SuperNode<T> * newSnode = nullptr;
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
          bassert(this->Mapping_->Map(I-1,I-1)==this->iam);
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
    {
      double timeSta = get_time();
      Icomm recv_buffer;


      AsyncComms incomingRecv;
      AsyncComms outgoingSend;
      Int numColFirst = std::max((Int)1,this->iSize_ / this->np);

      logfileptr->OFS()<<"Starting Send"<<std::endl;

#ifndef NOTRY
      try
#endif
      {
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

        remoteFactors_.resize(this->Xsuper_.size()-1);
        std::fill((char*)&remoteFactors_[0],(char*)&remoteFactors_[0]+remoteFactors_.size()*sizeof(std::tuple<upcxx::global_ptr<char>,Int> ),0);

        for(Int I=1;I<this->Xsuper_.size();I++){
          Int iDest = this->Mapping_->Map(I-1,I-1);
          //parse the first column to create the supernode structure
          if(this->iam==iDest){
            SuperNode<T> * newSnode = snodeLocal(I);
            SuperNodeDesc * meta = newSnode->GetMeta();
            remoteFactors_[I-1] = std::make_tuple( upcxx::to_global_ptr<char>( (char*)meta ), meta->blocks_cnt_) ;
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
        symPACKOS << "Allocation failed: " << e.what() << '\n';
        abort();
      }
#endif
      logfileptr->OFS()<<"Send Done"<<std::endl;

      double timeStop = get_time();
      if(this->iam==0 && this->options_.verbose){
        symPACKOS<<"Distribution time: "<<timeStop - timeSta<<std::endl;
      }
    }
    MPI_Barrier(pMat.comm);
  }




  template <typename T> inline symPACKMatrix<T>::symPACKMatrix():
    symPACKMatrixMeta<T>(){
    CommEnv_=nullptr;
    Mapping_ = nullptr;
    Balancer_ = nullptr;

    if(logfileptr==nullptr){
      logfileptr = new LogFile(upcxx::rank_me(),false);
      logfileptr->OFS()<<"********* LOGFILE OF P"<<upcxx::rank_me()<<" *********"<<std::endl;
      logfileptr->OFS()<<"**********************************"<<std::endl;
    }

    logfileptr->OFS()<<"Shared node size "<<shmNode_.shmsize<<", rank "<<shmNode_.shmrank<<std::endl;


  }

  template <typename T> inline symPACKMatrix<T>::symPACKMatrix(DistSparseMatrix<T> & pMat, symPACKOptions & options ):symPACKMatrix(){
    Init(options);
    DistributeMatrix(pMat);
  }

  template <typename T> inline symPACKMatrix<T>::~symPACKMatrix(){
    for(Int i=0;i<LocalSupernodes_.size();++i){
      delete LocalSupernodes_[i];
    }

    for(Int i=0;i<Contributions_.size();++i){
      delete Contributions_[i];
    }


    if(this->Mapping_!=nullptr){
      delete this->Mapping_;
    }

    if(this->Balancer_!=nullptr){
      delete this->Balancer_;
    }

    if(CommEnv_!=nullptr){
      delete CommEnv_;
    }

    if(this->fullcomm_!=MPI_COMM_NULL){
      MPI_Comm_free(&this->fullcomm_);
    }

    if(this->workcomm_!=MPI_COMM_NULL){
      MPI_Comm_free(&this->workcomm_);
    }


  }


  //returns the 1-based index of supernode id global in the local supernode array
  template <typename T> inline Int symPACKMatrix<T>::snodeLocalIndex(Int global){
#ifndef ITREE2
    bassert(global<=globToLocSnodes_.back()+1);
    auto it = std::lower_bound(globToLocSnodes_.begin(),globToLocSnodes_.end(),global);
    return it - globToLocSnodes_.begin();
#else
    ITree::Interval * ptr = globToLocSnodes_.IntervalSearch(global,global);
    assert(ptr!=nullptr);
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
    template <class Allocator>
    inline SuperNode<T,Allocator> * symPACKMatrix<T>::CreateSuperNode(DecompositionType type,Int aiId, Int aiFr, Int aiFc, Int aiLc, Int ai_num_rows, Int aiN, Int aiNZBlkCnt, Int panel){
      SuperNode<T,Allocator> * retval = nullptr;
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
      SuperNode<T,Allocator> * retval = nullptr;
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
      SuperNode<T,Allocator> * retval = nullptr;
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
      SuperNode<T,Allocator> * retval = nullptr;
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



}


#endif //_SYMPACK_MATRIX_IMPL_HP_

