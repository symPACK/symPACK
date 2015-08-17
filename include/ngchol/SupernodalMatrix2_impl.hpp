#ifndef _SUPERNODAL_MATRIX_IMPL_HP_
#define _SUPERNODAL_MATRIX_IMPL_HPP_

#include "ngchol/SupernodalMatrix2.hpp"

#include "ngchol/blas.hpp"
#include "ngchol/Ordering.hpp"
#include "ngchol/LoadBalancer.hpp"
#include "ngchol/mpi_interf.hpp"

#include <queue>
#include <stdlib.h>

#ifdef _DEADLOCK_
#define TOP() front()
#else
#define TOP() top()
#endif

namespace LIBCHOLESKY{

template <typename T> void SupernodalMatrix2<T>::Factorize(){
  TIMER_START(FACTORIZATION);
  switch(options_.factorization){
    case FANBOTH:
      FanBoth();
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

  template <typename T> void SupernodalMatrix2<T>::forward_update(SuperNode2<T> * src_contrib,SuperNode2<T> * tgt_contrib){

    Int iam = CommEnv_->MPI_Rank();
    Int np  = CommEnv_->MPI_Size();

    Int iOwner = this->Mapping_->Map(src_contrib->Id()-1,src_contrib->Id()-1);
    Int src_ncols = src_contrib->Size();
    Int tgt_ncols = tgt_contrib->Size();

    Int startBlock = (iam==iOwner)?1:0;
    Int tgt_blkidx = -1;

    Int src_blkidx = startBlock;
    while(src_blkidx<src_contrib->NZBlockCnt()){
      NZBlockDesc2 & src_desc = src_contrib->GetNZBlockDesc(src_blkidx);
      Int src_nrows = src_contrib->NRows(src_blkidx);
  
      if(tgt_blkidx<1){
          tgt_blkidx = tgt_contrib->FindBlockIdx(src_desc.GIndex);
      }

      NZBlockDesc2 & tgt_desc = tgt_contrib->GetNZBlockDesc(tgt_blkidx);
      Int tgt_nrows = tgt_contrib->NRows(tgt_blkidx);

      Int src_local_fr = max(tgt_desc.GIndex - src_desc.GIndex,0);
      Int src_lr = src_desc.GIndex+src_nrows-1;

      Int tgt_local_fr = max(src_desc.GIndex - tgt_desc.GIndex,0);
      Int tgt_lr = tgt_desc.GIndex+tgt_nrows-1;
      Int tgt_local_lr = min(src_lr,tgt_lr) - tgt_desc.GIndex;

      T * src = &src_contrib->GetNZval(src_desc.Offset)[src_local_fr*src_ncols];
      T * tgt = &tgt_contrib->GetNZval(tgt_desc.Offset)[tgt_local_fr*tgt_ncols];


//TODO understand why this axpy makes the code crash
//      for(Int i=0; i<(tgt_local_lr - tgt_local_fr +1)*src_ncols;++i,++tgt,++src){
//        *tgt+=*src;
//      }
      for(Int i=0; i<(tgt_local_lr - tgt_local_fr +1)*src_ncols;++i){
        tgt[i]+=src[i];
      }
//      blas::Axpy((tgt_local_lr - tgt_local_fr +1)*src_ncols,
//          ONE<T>(),src,1,tgt,1);

      if(src_lr>tgt_lr){
        //the src block hasn't been completely used and is
        // updating some lines in the nz block just below the diagonal block
        //this case happens only locally for the diagonal block
//        assert(tgt_blkidx==0);
        //skip to the next tgt nz block
        ++tgt_blkidx;
      }
      else{
        //skip to the next src nz block
        ++src_blkidx;
        if(src_blkidx<src_contrib->NZBlockCnt()){
          NZBlockDesc2 & next_src_desc = src_contrib->GetNZBlockDesc(src_blkidx);
          tgt_blkidx = tgt_contrib->FindBlockIdx(next_src_desc.GIndex);
        }
      }
    }
  }

  template <typename T> void SupernodalMatrix2<T>::back_update(SuperNode2<T> * src_contrib,SuperNode2<T> * tgt_contrib){
    Int nrhs = tgt_contrib->Size();
    for(Int blkidx = 1; blkidx<tgt_contrib->NZBlockCnt();++blkidx){
      NZBlockDesc2 & tgt_desc = tgt_contrib->GetNZBlockDesc(blkidx);
      Int tgt_nrows = tgt_contrib->NRows(blkidx);

      Int src_nzblk_idx = src_contrib->FindBlockIdx(tgt_desc.GIndex);
      Int tgt_lr = tgt_desc.GIndex+tgt_nrows-1;
      Int src_lr; 
      do
      {
        NZBlockDesc2 & src_desc = src_contrib->GetNZBlockDesc(src_nzblk_idx);
        Int src_nrows = src_contrib->NRows(src_nzblk_idx);
        src_lr = src_desc.GIndex+src_nrows-1;

        Int src_local_fr = max(tgt_desc.GIndex - src_desc.GIndex,0);

        Int tgt_local_fr = max(src_desc.GIndex - tgt_desc.GIndex,0);
        Int tgt_local_lr = min(src_lr,tgt_lr) - tgt_desc.GIndex;

        T * src = &src_contrib->GetNZval(src_desc.Offset)[src_local_fr*nrhs];
        T * tgt = &tgt_contrib->GetNZval(tgt_desc.Offset)[tgt_local_fr*nrhs];

        std::copy(src,src+(tgt_local_lr - tgt_local_fr +1)*nrhs,tgt);
        //              lapack::Lacpy('N',nrhs,(tgt_local_lr - tgt_local_fr +1),
        //                  src,nrhs,  
        //                  tgt,nrhs);

        //do other block
        if(tgt_lr>src_lr){
//          assert(src_nzblk_idx==0);
          src_nzblk_idx++;
        }

      } while(tgt_lr>src_lr);
    }
  }

  template <typename T> void SupernodalMatrix2<T>::Solve(NumMat<T> * RHS,  NumMat<T> * Xptr) {
    TIMER_START(SPARSE_SOLVE);

    NumMat<T> & B = *RHS;
    Int nrhs = B.n();

    Int iam = CommEnv_->MPI_Rank();
    Int np  = CommEnv_->MPI_Size();

    IntNumVec children(Xsuper_.m());
    SetValue(children,0);

    for(Int I=1;I<Xsuper_.m()-1;I++){
      Int fc = Xsuper_(I-1);
      Int lc = Xsuper_(I)-1;

      Int parent = ETree_.PostParent(lc-1);
      if(parent!=0){
        ++children(SupMembership_(parent-1)-1);
      }
    }

#ifdef _DEBUG_
    logfileptr->OFS()<<"Children vector is"<<children<<std::endl;
#endif


    IntNumVec UpdatesToDo = children;

    Contributions_.resize(LocalSupernodes_.size());
    std::vector<std::stack<Int> > LocalUpdates(LocalSupernodes_.size());

    AsyncComms outgoingSend;

    //This corresponds to the k loop in dtrsm
    for(Int I=1;I<Xsuper_.m();I++){
      Int iOwner = this->Mapping_->Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam ){
        //Create the blocks of my contrib with the same nz structure as L
        //and the same width as the final solution
        //Int iLocalI = (I-1) / np +1 ;
        Int iLocalI = snodeLocalIndex(I);
        SuperNode2<T> * cur_snode = LocalSupernodes_[iLocalI-1];
        SuperNode2<T> * contrib = new SuperNode2<T>(I,1,nrhs, cur_snode->NRowsBelowBlock(0) ,iSize_);
        Contributions_[iLocalI-1] = contrib;


        for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){
          NZBlockDesc2 & cur_desc = cur_snode->GetNZBlockDesc(blkidx);
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
    TIMER_START(SPARSE_FWD_SUBST);

    Int I =1;
    Int iLocalI =1;
    while(iLocalI<=LocalSupernodes_.size() || !ContribsToSend.empty() || !outgoingSend.empty()){
      //Check for completion of outgoing communication
      AdvanceOutgoing(outgoingSend);

      //process some of the delayed send

      if(iLocalI>0 && iLocalI<=LocalSupernodes_.size()){

      //If I own the column, factor it
        SuperNode2<T> * cur_snode = LocalSupernodes_[iLocalI-1];
        I = cur_snode->Id();
        Int parent = ETree_.PostParent(cur_snode->LastCol()-1);
        SuperNode2<T> * contrib = Contributions_[iLocalI-1];


        //Do all my updates (Local and remote)
        //Local updates
        while(!LocalUpdates[iLocalI-1].empty()){
          Int contrib_snode_id = LocalUpdates[iLocalI-1].top();
          LocalUpdates[iLocalI-1].pop();

          SuperNode2<T> * dist_contrib = snodeLocal(contrib_snode_id,Contributions_);

#ifdef _DEBUG_
          logfileptr->OFS()<<"LOCAL Supernode "<<I<<" is updated by contrib of Supernode "<<contrib_snode_id<<std::endl;
#endif

          forward_update(dist_contrib,contrib);
          //delete contributions[(contrib_snode_id-1) / np];
          --UpdatesToDo(I-1);
        }

        //do remote updates
        size_t max_bytes;
        Int nz_cnt;
        while(UpdatesToDo(I-1)>0){
          //receive children contrib
#ifdef _DEBUG_
          logfileptr->OFS()<<UpdatesToDo(I-1)<<" contribs left"<<endl;
#endif


          MPI_Status recv_status;
          int bytes_received = 0;

          TIMER_START(RECV_MPI);
          MPI_Probe(MPI_ANY_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
          src_blocks.resize(bytes_received);
          MPI_Recv(&src_blocks[0],bytes_received,MPI_BYTE,recv_status.MPI_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
          TIMER_STOP(RECV_MPI);
          SuperNode2<T> dist_contrib(&src_blocks[0],bytes_received);
          //TODO Put this back
          //Deserialize(&src_blocks[0],dist_contrib);
#ifdef _DEBUG_
          logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib.Id()<<std::endl;
#endif

          forward_update(&dist_contrib,contrib);

          --UpdatesToDo(I-1);

        }

        assert(UpdatesToDo(I-1)==0);

        if(UpdatesToDo(I-1)==0){

          //This corresponds to the i loop in dtrsm
          for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){

            NZBlockDesc2 & cur_desc = contrib->GetNZBlockDesc(blkidx);
            NZBlockDesc2 & chol_desc = cur_snode->GetNZBlockDesc(blkidx);
            NZBlockDesc2 & diag_desc = contrib->GetNZBlockDesc(0);

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
                  //logfileptr->OFS()<<"                      B = "<<B(Order_.perm[diag_desc.GIndex+kk-1]-1,j)<<std::endl;
                  Int srcRow = Order_.perm[diag_desc.GIndex+kk-1];
                  diag_nzval[kk*nrhs+j] = (B(srcRow-1,j) + diag_nzval[kk*nrhs+j]) / chol_nzval[kk*cur_snode->Size()+kk];
                  //diag_nzval[kk*nrhs+j] = (B(diag_desc.GIndex-1+kk,j) + diag_nzval[kk*nrhs+j]) / chol_nzval[kk*cur_snode->Size()+kk];
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


          //send to my parent
          if(parent!=0){
            Int parent_snode_id = SupMembership_[parent-1];

            Int iTarget = this->Mapping_->Map(parent_snode_id-1,parent_snode_id-1);

            if(iTarget!=iam){
#ifdef _DEBUG_
              logfileptr->OFS()<<"Remote Supernode "<<parent_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif

              Int tgt_first_col = Xsuper_(parent_snode_id-1);
              Int tgt_last_col = Xsuper_(parent_snode_id)-1;
              Int src_nzblk_idx = 1;
              NZBlockDesc2 & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);

              Int src_first_row = pivot_desc.GIndex;

              bool isSkipped= false;

              Int next_local_contrib = (iLocalI < Contributions_.size())?Contributions_[iLocalI]->Id():Xsuper_.m();
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
                    Serialize(*send_buffer,*contrib,src_nzblk_idx,src_first_row);
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

      NumMat<T> tmp3(tmp.m(),tmp.n());
      NumMat<T> tmpFwd(tmp.m(),tmp.n());
      for(Int i = 0; i<tmp.m();++i){
        for(Int j = 0; j<tmp.n();++j){
//          tmp3(Order_.perm[i]-1,j) = tmp(i,j);
//          tmpFwd(Order_.perm[i]-1,j) = forwardSol(i,j);

          tmp3(i,j) = tmp(Order_.perm[i]-1,j);
          tmpFwd(i,j) = forwardSol(Order_.perm[i]-1,j);
        }
      }

      NumMat<T> tmp2 = tmp3;
      
      blas::Axpy(tmp.m()*tmp.n(),-1.0,&tmpFwd(0,0),1,&tmp2(0,0),1);
      double norm = lapack::Lange('F',nrows,tmp.n(),&tmp2(0,0),tmp.m());
      logfileptr->OFS()<<"Norm after SuperNode2 "<<I<<" is "<<norm<<std::endl; 

        if(abs(norm)>=1e-1){
          for(Int i = 0;i<nrows/*tmp.m()*/;++i){
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


    //Back-substitution phase
    TIMER_START(SPARSE_BACK_SUBST);

    //start from the root of the tree
    iLocalI = LocalSupernodes_.size() ;
    while(iLocalI>0|| !ContribsToSend.empty() || !outgoingSend.empty()){

      //Check for completion of outgoing communication
      AdvanceOutgoing(outgoingSend);

      if(iLocalI>0 && iLocalI<=LocalSupernodes_.size()){
        SuperNode2<T> * cur_snode = LocalSupernodes_[iLocalI-1];
        I = cur_snode->Id();

        SuperNode2<T> * contrib = Contributions_[iLocalI-1];

        Int parent = ETree_.PostParent(cur_snode->LastCol()-1);

        std::vector<Int> isBlockUpdated(contrib->NZBlockCnt(),0);

        if(parent!=0){
          Int parent_snode_id = SupMembership_[parent-1];
          Int iTarget = this->Mapping_->Map(parent_snode_id-1,parent_snode_id-1);
          //Do all my updates (Local and remote)
          //Local updates
          SuperNode2<T> * dist_contrib;
          if(!LocalUpdates[iLocalI-1].empty()){
            Int contrib_snode_id = LocalUpdates[iLocalI-1].top();
            LocalUpdates[iLocalI-1].pop();
            dist_contrib = snodeLocal(contrib_snode_id,Contributions_);
            back_update(dist_contrib,contrib);
          }
          else{
            //Receive parent contrib
            MPI_Status recv_status;
            Int bytes_received = 0;
            MPI_Probe(iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);
            MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
            src_blocks.resize(bytes_received);

            MPI_Recv(&src_blocks[0],bytes_received,MPI_BYTE,iTarget,I,CommEnv_->MPI_GetComm(),&recv_status);

            dist_contrib = new SuperNode2<T>(&src_blocks[0],bytes_received);
            //TODO Replace this
            //Deserialize(&src_blocks[0],*dist_contrib); 
            dist_contrib->InitIdxToBlk();


            #ifdef _DEBUG_
            logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib->Id()<<std::endl;
            #endif

            back_update(dist_contrib,contrib);
            delete dist_contrib;
          }
        }

        //now compute MY contribution
#ifdef _DEBUG_
        logfileptr->OFS()<<"BACK Processing contrib "<<I<<std::endl;
#endif

        NZBlockDesc2 & diag_desc = cur_snode->GetNZBlockDesc(0);
        NZBlockDesc2 & tgt_desc = contrib->GetNZBlockDesc(0);

        T* diag_nzval = cur_snode->GetNZval(diag_desc.Offset);
        T* tgt_nzval = contrib->GetNZval(tgt_desc.Offset);

        for(Int j = 0; j<nrhs;++j){
          for(Int ii = cur_snode->Size()-1; ii>=0; --ii){
            T temp = tgt_nzval[ii*nrhs+j];

            //This corresponds to the k loop in dtrsm
            for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){
              NZBlockDesc2 & chol_desc = cur_snode->GetNZBlockDesc(blkidx);
              Int chol_nrows = cur_snode->NRows(blkidx);

              Int src_blkidx = contrib->FindBlockIdx(chol_desc.GIndex);
              NZBlockDesc2 & cur_desc = contrib->GetNZBlockDesc(src_blkidx);
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


        //send to my children
        Int colIdx = cur_snode->FirstCol()-1;
        if(colIdx>0){
          Int children_found = 0;
          while(children_found<children(I-1)){
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
                    NZBlockDesc2 & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);
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
                    NZBlockDesc2 & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);
                    Int src_first_row = pivot_desc.GIndex;
                    //Create a new Icomm buffer, serialize the contribution
                    // in it and add it to the outgoing comm list
                    Icomm * send_buffer = new Icomm();
                    //TODO replace this
                    Serialize(*send_buffer,*contrib,src_nzblk_idx,src_first_row);
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
    TIMER_STOP(SPARSE_SOLVE);

}

template<typename T> void SupernodalMatrix2<T>::GetSolution(NumMat<T> & B){
    Int iam = CommEnv_->MPI_Rank();
    Int np  = CommEnv_->MPI_Size();
    Int nrhs = B.n();
    //Gather B from everybody and put it in the original matrix order
    std::vector<T> tmp_nzval;
    for(Int I=1; I<Xsuper_.m();++I){
      Int iOwner = this->Mapping_->Map(I-1,I-1);
      T * data;
      Int snode_size = Xsuper_[I] - Xsuper_[I-1];
      Int nzcnt = snode_size * nrhs;
      tmp_nzval.resize(nzcnt);

      if( iOwner == iam ){
        SuperNode2<T> * contrib = snodeLocal(I,Contributions_);
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
          B(destRow-1, j) = data[i*nrhs + j];
        }
      }
    }
      MPI_Barrier(CommEnv_->MPI_GetComm());
}


template <typename T> void SupernodalMatrix2<T>::SendDelayedMessagesUp(Int iLocalI, CommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode2<T> *> & snodeColl){
  if(snodeColl.empty() || MsgToSend.empty()) { return;}

  //Index of the last global snode to do
  Int last_snode_id = Xsuper_.m()-1;
  //Index of the last local supernode
  Int last_local_id = snodeColl.back()->Id();
  //Index of the last PROCESSED supernode
  Int prev_snode_id = iLocalI<=snodeColl.size()?snodeColl[iLocalI-1]->Id():last_local_id;
  //Index of the next local supernode
  Int next_snode_id = prev_snode_id>=last_local_id?last_snode_id+1:snodeColl[iLocalI]->Id();

  bool is_last = prev_snode_id>=last_local_id;

  while( MsgToSend.size()>0){
    //Pull the highest priority message
#ifdef _DEADLOCK_
    const DelayedComm & comm = MsgToSend.front();
#else
    const DelayedComm & comm = MsgToSend.top();
#endif
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


template <typename T> void SupernodalMatrix2<T>::SendDelayedMessagesDown(Int iLocalI, DownCommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode2<T> *> & snodeColl){
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
    const DelayedComm & comm = MsgToSend.TOP();

    Int src_snode_id = comm.src_snode_id;
    Int tgt_snode_id = comm.tgt_snode_id;
    Int src_nzblk_idx = comm.src_nzblk_idx;
    Int src_first_row = comm.src_first_row;

#ifdef _DEBUG_DELAY_
    logfileptr->OFS()<<"Picked { "<<src_snode_id<<" -> "<<tgt_snode_id<<" }"<<endl;
#endif

    if(tgt_snode_id>next_snode_id || is_last /*|| OutgoingSend.size() <= maxIsend_*/){

      SuperNode2<T> & prev_src_snode = *snodeLocal(src_snode_id,snodeColl);
      //this can be sent now

      Int iTarget = this->Mapping_->Map(tgt_snode_id-1,tgt_snode_id-1);
      if(iTarget != iam){
#ifdef _DEBUG_DELAY_
        logfileptr->OFS()<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode.Id()<<" to Supernode "<<tgt_snode_id<<endl;
        cout<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode.Id()<<" to Supernode "<<tgt_snode_id<<endl;
#endif

#ifdef _DEBUG_
        logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated by Supernode "<<prev_src_snode.Id()<<" rows "<<src_first_row/*<<" to "<<src_last_row*/<<" "<<src_nzblk_idx<<std::endl;
#endif

        NZBlockDesc2 & pivot_desc = prev_src_snode.GetNZBlockDesc(src_nzblk_idx);
        //Create a new Icomm buffer, serialize the contribution
        // in it and add it to the outgoing comm list
        Icomm * send_buffer = new Icomm();
        //TODO replace this
        Serialize(*send_buffer,prev_src_snode,src_nzblk_idx,src_first_row);
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

template <typename T> void SupernodalMatrix2<T>::SendMessage(const DelayedComm & comm, AsyncComms & OutgoingSend, std::vector<SuperNode2<T> *> & snodeColl){
    Int src_snode_id = comm.src_snode_id;
    Int tgt_snode_id = comm.tgt_snode_id;
    Int src_nzblk_idx = comm.src_nzblk_idx;
    Int src_first_row = comm.src_first_row;


      SuperNode2<T> & prev_src_snode = *snodeLocal(src_snode_id,snodeColl);

      //this can be sent now
      Int iTarget = this->Mapping_->Map(tgt_snode_id-1,tgt_snode_id-1);
      if(iTarget != iam){
#ifdef _DEBUG_DELAY_
        logfileptr->OFS()<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode.Id()<<" to Supernode "<<tgt_snode_id<<endl;
        cout<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode.Id()<<" to Supernode "<<tgt_snode_id<<endl;
#endif

#ifdef _DEBUG_
        logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated by Supernode "<<prev_src_snode.Id()<<" rows "<<src_first_row/*<<" to "<<src_last_row*/<<" "<<src_nzblk_idx<<std::endl;
#endif
        NZBlockDesc2 & pivot_desc = prev_src_snode.GetNZBlockDesc(src_nzblk_idx);
        //Create a new Icomm buffer, serialize the contribution
        // in it and add it to the outgoing comm list
        Icomm * send_buffer = new Icomm();
        //TODO replace this
        Serialize(*send_buffer,prev_src_snode,src_nzblk_idx,src_first_row);
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




template <typename T> void SupernodalMatrix2<T>::AdvanceOutgoing(AsyncComms & outgoingSend){
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


  template<typename T> void SupernodalMatrix2<T>::AddOutgoingComm(AsyncComms & outgoingSend, Icomm * send_buffer){
    outgoingSend.push_back(send_buffer);
  }


  template <typename T> void SupernodalMatrix2<T>::GetUpdatingSupernodeCount(std::vector<Int> & sc,std::vector<Int> & mw, std::vector<Int> & mh){
    sc.resize(Xsuper_.m(),I_ZERO);
    std::vector<Int> marker(Xsuper_.m(),I_ZERO);
    mw.resize(Xsuper_.m(),I_ZERO);
    mh.resize(Xsuper_.m(),I_ZERO);

    for(Int s = 1; s<Xsuper_.m(); ++s){
      Int first_col = Xsuper_[s-1];
      Int last_col = Xsuper_[s]-1;

      Ptr fi = xlindx_[s-1];
      Ptr li = xlindx_[s]-1;

      mh[s-1] = li-fi+1;

      Int iOwner = this->Mapping_->Map(s-1,s-1);
#ifdef _DEBUG_UPDATES_
      logfileptr->OFS()<<"Supernode "<<s<<" on P"<<iOwner<<" updates: ";
#endif

      for(Ptr row_idx = fi; row_idx<=li;++row_idx){
        Idx row = lindx_[row_idx-1];
        Int supno = SupMembership_(row-1);

        if(marker[supno-1]!=s && supno!=s){

#ifdef _DEBUG_UPDATES_
          logfileptr->OFS()<<supno<<" ";
#endif
          ++sc[supno-1];
          marker[supno-1] = s;

          mw[supno-1] = max(mw[supno-1],last_col - first_col+1);

        }
      }

#ifdef _DEBUG_UPDATES_
      logfileptr->OFS()<<std::endl;
#endif
    }
  }


  template <typename T> SparseMatrixStructure SupernodalMatrix2<T>::GetLocalStructure() const {
    return *Local_;
  }

  template <typename T> SparseMatrixStructure SupernodalMatrix2<T>::GetGlobalStructure(){
    if(isGlobStructAllocated_){
      Global_=new SparseMatrixStructure();
      Local_->ToGlobal(*Global_,CommEnv_->MPI_GetComm());
      isGlobStructAllocated_= true;
    }
    return *Global_;
  }


template<typename T>
void SupernodalMatrix2<T>::Dump(){
    for(Int I=1;I<Xsuper_.m();I++){
      Int src_first_col = Xsuper_(I-1);
      Int src_last_col = Xsuper_(I)-1;
      Int iOwner = this->Mapping_->Map(I-1,I-1);
      if( iOwner == iam ){
        SuperNode2<T> & src_snode = *snodeLocal(I);
          

        logfileptr->OFS()<<"+++++++++++++"<<I<<"++++++++"<<std::endl;

          logfileptr->OFS()<<"cols: ";
          for(Int i = 0; i< src_snode.Size(); ++i){
              logfileptr->OFS()<<" "<<Order_.perm[src_first_col+i-1];
          }
            logfileptr->OFS()<<std::endl;
            logfileptr->OFS()<<std::endl;
        for(int blkidx=0;blkidx<src_snode.NZBlockCnt();++blkidx){

          NZBlockDesc2 & desc = src_snode.GetNZBlockDesc(blkidx);
          T * val = src_snode.GetNZval(desc.Offset);
          Int nRows = src_snode.NRows(blkidx);

          Int row = desc.GIndex;
          for(Int i = 0; i< nRows; ++i){
              logfileptr->OFS()<<row+i<<" | "<<Order_.perm[row+i-1]<<":   ";
            for(Int j = 0; j< src_snode.Size(); ++j){
              logfileptr->OFS()<<val[i*src_snode.Size()+j]<<" ";
            }
            logfileptr->OFS()<<std::endl;
          }

        logfileptr->OFS()<<"_______________________________"<<std::endl;
        }
      }
    }

}


  template <typename T> void SupernodalMatrix2<T>::Init(const DistSparseMatrix<T> & pMat, NGCholOptions & options ){


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


    //create the mapping
    Int pmapping = options.used_procs(np);
    Int pr = (Int)sqrt((double)pmapping);
    switch(options.mappingType){
      case WRAP2D:
        this->Mapping_ = new Wrap2D(pmapping, pr, pr, 1);
        break;
      case WRAP2DFORCED:
        this->Mapping_ = new Wrap2DForced(pmapping, pr, pr, 1);
        break;
      case MODWRAP2DNS:
        this->Mapping_ = new Modwrap2DNS(pmapping, pr, pr, 1);
        break;
      case ROW2D:
        this->Mapping_ = new Row2D(pmapping, pmapping, pmapping, 1);
        break;
      case COL2D:
        this->Mapping_ = new Col2D(pmapping, pmapping, pmapping, 1);
        break;
      case MODWRAP2D: default:
        this->Mapping_ = new Modwrap2D(pmapping, pr, pr, 1);
        break;
    }

    //Options
    maxIsend_ = options.maxIsend;
    maxIrecv_ = options.maxIrecv;

    gMaxIrecv = maxIrecv_;


    switch(options_.scheduler){
      case DL:
        scheduler_ = new DLScheduler<std::list<FBTask>::iterator>();
        break;
      case MCT:
        scheduler_ = new MCTScheduler<std::list<FBTask>::iterator>();
        break;
      case PR:
        scheduler_ = new PRScheduler<std::list<FBTask>::iterator>();
        break;
      default:
        scheduler_ = new DLScheduler<std::list<FBTask>::iterator>();
        break;
    }



    iSize_ = pMat.size;
    Local_ = new SparseMatrixStructure();
    Global_ = new SparseMatrixStructure();
    *Local_ = pMat.GetLocalStructure();

    Local_->ToGlobal(*Global_,CommEnv_->MPI_GetComm());
    Global_->ExpandSymmetric();
    isGlobStructAllocated_ = true;

logfileptr->OFS()<<"Matrix expanded"<<endl;

    //Create an Ordering object to hold the permutation
    Order_.SetStructure(*Global_);
logfileptr->OFS()<<"Structure set"<<endl;

    switch(options_.ordering){
      case MMD:
        //Reoder the matrix with MMD
        Order_.MMD();
        break;
      case AMD:
        Order_.AMD();
        break;
      default:
        Order_.MMD();
        break;
    }
logfileptr->OFS()<<"Ordering done"<<endl;

    ETree_.ConstructETree(*Global_,Order_);
    ETree_.PostOrderTree(Order_);
logfileptr->OFS()<<"ETREE done"<<endl;
    IntNumVec cc,rc;
    Global_->GetLColRowCount(ETree_,Order_,cc,rc);
    ETree_.SortChildren(cc,Order_);

#ifdef _DEBUG_
    logfileptr->OFS()<<"colcnt "<<cc<<std::endl;
    logfileptr->OFS()<<"rowcnt "<<rc<<std::endl;
#endif
    double flops = 0.0;
    for(Int i = 0; i<cc.m();++i){
      flops+= (double)pow((double)cc[i],2.0);
    }

    if(iam==0){
      cout<<"Flops: "<<flops<<endl;
    }

    Global_->FindSupernodes(ETree_,Order_,cc,SupMembership_,Xsuper_,options.maxSnode);

logfileptr->OFS()<<"Supernodes found"<<endl;

#ifdef RELAXED_SNODE
    Global_->RelaxSupernodes(ETree_, cc,SupMembership_, Xsuper_, options.maxSnode );
logfileptr->OFS()<<"Relaxation done"<<endl;
    Global_->SymbolicFactorizationRelaxed(ETree_,Order_,cc,Xsuper_,SupMembership_,xlindx_,lindx_);

#else
    Global_->SymbolicFactorization(ETree_,Order_,cc,Xsuper_,SupMembership_,xlindx_,lindx_);
#endif

logfileptr->OFS()<<"Symbfact done"<<endl;

#ifdef REFINED_SNODE
    IntNumVec permRefined;
    //    IntNumVec newPerm(Size());
    Global_->RefineSupernodes(ETree_,Order_, SupMembership_, Xsuper_, xlindx_, lindx_, permRefined);

    //      Perm_.Resize(Size());
    //      for(Int i =0; i<Perm_.m();++i){
    //        Perm_[permRefined[i]-1] = permChild[i];
    //      }

    //    Perm_ = permRefined;
#else
    //Combine permMMD with Postorder
    //    Perm_.Resize(Size());
    //    for(Int i =0; i<Perm_.m();++i){
    //      Perm_[i] = ETree_.FromPostOrder(i+1);
    //    }
#endif

#ifdef _DEBUG_
    logfileptr->OFS()<<"Membership list is "<<SupMembership_<<std::endl;
    logfileptr->OFS()<<"xsuper "<<Xsuper_<<std::endl;
#endif


#ifdef _OUTPUT_ETREE_
    logfileptr->OFS()<<"ETree is "<<ETree_<<std::endl;
    logfileptr->OFS()<<"Supernodal ETree is "<<ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_)<<std::endl;
#endif


isXlindxAllocated_=true;
isLindxAllocated_=true;

    switch(options.load_balance){
      case SUBCUBE:
        {
    LoadBalancer * balancer;
          if(iam==0){ cout<<"Subtree to subcube mapping used"<<endl;}
          ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
          balancer = new SubtreeToSubcube(np,SupETree,Xsuper_,SupMembership_,xlindx_,lindx_,cc);
////
////
////          //Int np = 4;
////
////          //compute children array and subtree costs
////
////          //        if(iam == 0){
////          //          cout<<"Supernodal ETree is "<<SupETree<<std::endl;
////          //        }
////
////          //compute number of children and load
////          vector<double> SubTreeLoad(SupETree.Size()+1,0.0);
////          vector<Int> children(SupETree.Size()+1,0);
////          for(Int I=1;I<=SupETree.Size();I++){
////            Int parent = SupETree.Parent(I-1);
////            ++children[parent];
////            Int fc = Xsuper_[I-1];
////            Int width = Xsuper_[I] - Xsuper_[I-1];
////            Int height = cc[fc-1];
////            double local_load = width*height*height;
////            SubTreeLoad[I]+=local_load;
////            SubTreeLoad[parent]+=SubTreeLoad[I];
////          }
////
////          logfileptr->OFS()<<"SubTreeLoad is "<<SubTreeLoad<<endl;
////
////
////          //procmaps[0]/pstart[0] represents the complete list
////          vector<vector<Int> * > procmaps(SupETree.Size()+1);
////          for(Int i = 0; i<procmaps.size();++i){ procmaps[i] = new vector<Int>();}
////          vector<Int> pstart(SupETree.Size()+1,0);
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
////              parent_load = width*height*height;
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
////          std::vector<Int> procMap(SupETree.Size());
////          std::vector<double> load(np,0.0);
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
////            double local_load = width*height*height;
////
////            load[minLoadP]+=local_load;
////          }
////
////
////          //for(Int i = 0; i<procMap.size();++i){ logfileptr->OFS()<<i+1<<" is on "<<procMap[i]<<endl;}
////          logfileptr->OFS()<<"Proc load: "<<load<<endl;
          //Update the mapping
          logfileptr->OFS()<<"Proc Mapping: "<<balancer->GetMap()<<endl;
          this->Mapping_->Update(balancer->GetMap());

////          for(Int i = 0; i<procmaps.size();++i){ delete procmaps[i];}
delete balancer;
        }       
        break;
      case SUBCUBE_NNZ:
        {
          if(iam==0){ cout<<"Subtree to subcube mapping used"<<endl;}

          //Int np = 4;

          //compute children array and subtree costs
          ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);

          //        if(iam == 0){
          //          cout<<"Supernodal ETree is "<<SupETree<<std::endl;
          //        }

          //compute number of children and load
          vector<double> SubTreeLoad(SupETree.Size()+1,0.0);
          vector<Int> children(SupETree.Size()+1,0);
          for(Int I=1;I<=SupETree.Size();I++){
            Int parent = SupETree.Parent(I-1);
            ++children[parent];
            Int fc = Xsuper_[I-1];
            Int width = Xsuper_[I] - Xsuper_[I-1];
            Int height = cc[fc-1];
            double local_load = width*height;
            SubTreeLoad[I]+=local_load;
            SubTreeLoad[parent]+=SubTreeLoad[I];
          }

          logfileptr->OFS()<<"SubTreeLoad is "<<SubTreeLoad<<endl;


          //procmaps[0]/pstart[0] represents the complete list
          vector<vector<Int> * > procmaps(SupETree.Size()+1);
          for(Int i = 0; i<procmaps.size();++i){ procmaps[i] = new vector<Int>();}
          vector<Int> pstart(SupETree.Size()+1,0);
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
              parent_load = width*height;
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
          std::vector<Int> procMap(SupETree.Size());
          std::vector<double> load(np,0.0);
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
            double local_load = width*height;

            load[minLoadP]+=local_load;
          }


          //for(Int i = 0; i<procMap.size();++i){ logfileptr->OFS()<<i+1<<" is on "<<procMap[i]<<endl;}
          logfileptr->OFS()<<"Proc load: "<<load<<endl;
          logfileptr->OFS()<<"Proc Mapping: "<<procMap<<endl;

          //Update the mapping
          this->Mapping_->Update(procMap);

          for(Int i = 0; i<procmaps.size();++i){ delete procmaps[i];}
        }       
        break;

      case NNZ:
        {
          LoadBalancer * balancer;
          if(iam==0){ cout<<"Load Balancing on NNZ used"<<endl;}
          balancer = new NNZBalancer(np,Xsuper_,cc);

//          //Do a greedy load balancing heuristic
//          std::vector<Int> procMap(Xsuper_.m()-1);
//          //Do a greedy heuristic to balance the number of nnz ?
//          std::vector<double> load(np,0.0);
//
//          for(Int i = 1; i< Xsuper_.m();  ++i){
//            //find least loaded processor
//            vector<double>::iterator it = std::min_element(load.begin(),load.end());
//            Int proc = (Int)(it - load.begin());
//            Int width = Xsuper_[i] - Xsuper_[i-1];
//            Int height = cc[i-1];
//            *it += (double)(width*height);
//            procMap[i-1] = proc;
//          } 
//
//          logfileptr->OFS()<<"Proc load: "<<load<<endl;
//          logfileptr->OFS()<<"Proc Mapping: "<<procMap<<endl;

          //Update the mapping
//          this->Mapping_->Update(procMap);

          logfileptr->OFS()<<"Proc Mapping: "<<balancer->GetMap()<<endl;
          this->Mapping_->Update(balancer->GetMap());

////          for(Int i = 0; i<procmaps.size();++i){ delete procmaps[i];}
delete balancer;
        }
        break;
      case FLOPS:
        {
          if(iam==0){ cout<<"Load Balancing on FLOPS used"<<endl;}
          //Do a greedy load balancing heuristic
          std::vector<Int> procMap(Xsuper_.m()-1);
          //Do a greedy heuristic to balance the number of nnz ?
          std::vector<double> load(np,0.0);

          for(Int i = 1; i< Xsuper_.m();  ++i){
            //find least loaded processor
            vector<double>::iterator it = std::min_element(load.begin(),load.end());
            Int proc = (Int)(it - load.begin());
            Int width = Xsuper_[i] - Xsuper_[i-1];
            Int height = cc[i-1];
            *it += (double)(height*height);
            procMap[i-1] = proc;
          } 

          logfileptr->OFS()<<"Proc load: "<<load<<endl;
          logfileptr->OFS()<<"Proc Mapping: "<<procMap<<endl;

          //Update the mapping
          this->Mapping_->Update(procMap);
        }
        break;

      default:
        {
        }
        break;
    }

#ifdef _DEBUG_MAPPING_
    this->Mapping_->Dump(2*np);
#endif

    GetUpdatingSupernodeCount(UpdateCount_,UpdateWidth_,UpdateHeight_);




    Int iam = CommEnv_->MPI_Rank();
    Int np  = CommEnv_->MPI_Size();

    DistSparseMatrix<T> ExpA(CommEnv_->MPI_GetComm());
    {
      double timeSta = get_time();
      ExpA.size = pMat.size;
      ExpA.nnz = 2*pMat.nnz-pMat.size;
      Int numColFirst = std::max(1,ExpA.size / np);
      std::vector<Int> numColLocalVec(np,numColFirst);
      numColLocalVec[np-1] = ExpA.size - numColFirst * (np-1);  // Modify the last entry	
      Int numColLocal = numColLocalVec[iam];
      //Expand A to symmetric storage
      Int localFirstCol = iam*numColFirst+1;
      Int localLastCol = localFirstCol+numColLocal-1;

      //Initialize the Local structure
      ExpA.Local_.size = ExpA.size;
      ExpA.Local_.colptr.Resize(numColLocal+1);

      for( Int i = 0; i < numColLocal + 1; i++ ){
        ExpA.Local_.colptr[i] = Global_->expColptr[iam * numColFirst+i] - Global_->expColptr[iam * numColFirst] + 1;
      }

      ExpA.Local_.nnz = ExpA.Local_.colptr[numColLocal] - ExpA.Local_.colptr[0];

      ExpA.Local_.rowind.Resize(ExpA.Local_.nnz);

      Int globalColBeg = Global_->expColptr[iam*numColFirst];
      std::copy(&Global_->expRowind[globalColBeg-1],&Global_->expRowind[globalColBeg-1]+ExpA.Local_.nnz,ExpA.Local_.rowind.Data());
      ExpA.nzvalLocal.Resize(ExpA.Local_.nnz);

#ifdef _DEBUG_
      logfileptr->OFS()<<"Global_.colptr: "<<Global_->colptr<<endl;
      logfileptr->OFS()<<"Global_.rowind: "<<Global_->rowind<<endl;
      logfileptr->OFS()<<"Global_.expColptr: "<<Global_->expColptr<<endl;
      logfileptr->OFS()<<"Global_.expRowind: "<<Global_->expRowind<<endl;
      logfileptr->OFS()<<"ExpA.colptr: "<<ExpA.Local_.colptr<<endl;
      logfileptr->OFS()<<"ExpA.rowind: "<<ExpA.Local_.rowind<<endl;
      logfileptr->OFS()<<"pMat.colptr: "<<pMat.Local_.colptr<<endl;
      logfileptr->OFS()<<"pMat.rowind: "<<pMat.Local_.rowind<<endl;
#endif

      IntNumVec localColHead = pMat.Local_.colptr;
      IntNumVec localDestColHead = ExpA.Local_.colptr;

      NumVec<T> recvNzval;
      IntNumVec recvColptr;
      IntNumVec recvRowind;
      for(Int proc = 0; proc<np; ++proc){
        //communicate colptr
        //Broadcast size
        Int size = iam==proc?pMat.Local_.colptr.m():0;
        MPI_Bcast(&size,sizeof(size),MPI_BYTE,proc,ExpA.comm);
        //Broadcast colptr
        if(iam==proc){
          MPI_Bcast(pMat.Local_.colptr.Data(),size*sizeof(Int),MPI_BYTE,proc,ExpA.comm);
        }
        else{
          recvColptr.Resize(size);
          MPI_Bcast(recvColptr.Data(),size*sizeof(Int),MPI_BYTE,proc,ExpA.comm);
        }

        //communicate rowind
        size = iam==proc?pMat.Local_.rowind.m():0;
        MPI_Bcast(&size,sizeof(size),MPI_BYTE,proc,ExpA.comm);
        //Broadcast rowind
        if(iam==proc){
          MPI_Bcast(pMat.Local_.rowind.Data(),size*sizeof(Int),MPI_BYTE,proc,ExpA.comm);
        }
        else{
          recvRowind.Resize(size);
          MPI_Bcast(recvRowind.Data(),size*sizeof(Int),MPI_BYTE,proc,ExpA.comm);
        }


        //communicate nzvalLocal
        //Broadcast size
        size = iam==proc?pMat.nzvalLocal.m():0;
        MPI_Bcast(&size,sizeof(size),MPI_BYTE,proc,ExpA.comm);
        //Broadcast nzvalLocal
        if(iam==proc){
          MPI_Bcast(pMat.nzvalLocal.Data(),size*sizeof(T),MPI_BYTE,proc,ExpA.comm);
        }
        else{
          recvNzval.Resize(size);
          MPI_Bcast(recvNzval.Data(),size*sizeof(T),MPI_BYTE,proc,ExpA.comm);
        }

        int colptrSize = iam==proc?pMat.Local_.colptr.m()-1:recvColptr.m()-1;
        int * colptr = iam==proc?pMat.Local_.colptr.Data():recvColptr.Data();
        int * rowind = iam==proc?pMat.Local_.rowind.Data():recvRowind.Data();
        T * nzval = iam==proc?pMat.nzvalLocal.Data():recvNzval.Data();
        //Parse the received data and copy it in the ExpA.nzvalLocal
        Int recvFirstCol = proc*numColFirst+1;
        Int recvLastCol = recvFirstCol+numColLocalVec[proc]-1;
        for(Int colIdx = 1; colIdx<=colptrSize;++colIdx){
          Int col = recvFirstCol + colIdx-1;
          Int colBeg = colptr[colIdx-1];
          Int colEnd = colptr[colIdx]-1;
          for(Int rowIdx = colBeg; rowIdx<=colEnd;++rowIdx){
            Int row = rowind[rowIdx-1];
            if(row>=localFirstCol && row<=localLastCol){
              //copy only upper triangular part
              if(col<row){
                Int localCol = row - localFirstCol +1;
                Int & localDestColIdx = localDestColHead[localCol-1];
                ExpA.nzvalLocal[localDestColIdx-1] = nzval[rowIdx-1];
                localDestColIdx++;
              }
            }
          }
        }
      }
      //copy upper triangular part
      for(Int colIdx = 1; colIdx<pMat.Local_.colptr.m();++colIdx){
        Int & localDestColIdx = localDestColHead[colIdx-1];
        Int colBeg = pMat.Local_.colptr[colIdx-1];
        Int colEnd = pMat.Local_.colptr[colIdx]-1;
        for(Int rowIdx = colBeg; rowIdx<=colEnd;++rowIdx){
          Int row = pMat.Local_.rowind[rowIdx-1];
          ExpA.nzvalLocal[localDestColIdx-1] = pMat.nzvalLocal[rowIdx-1];
          localDestColIdx++;
        }
      }
      double timeEnd = get_time();


#ifdef _DEBUG_
      logfileptr->OFS()<<"ExpA.nzvalLocal: "<<ExpA.nzvalLocal<<endl;
      logfileptr->OFS()<<"pMat.nzvalLocal: "<<pMat.nzvalLocal<<endl;
#endif

      //now we can do the copy by doing sends of whole columns to the dest processor
      if(iam==0){
        cout<<"Time for expanding matrix to asymmetric: "<<timeEnd-timeSta<<endl;
      }
    }

    //copy

    TIMER_START(DISTRIBUTING_MATRIX);
    Icomm buffer;
    std::vector<T> denseA;
    Icomm recv_buffer;


    AsyncComms incomingRecv;
    AsyncComms outgoingSend;
      Int numColFirst = std::max(1,iSize_ / np);

      logfileptr->OFS()<<"Starting Send"<<endl;
      try{
        //first, count 
        map<Int,pair<size_t,Icomm *> > recv_map;
        map<Int,pair<size_t,Icomm *> > send_map;
        Int snodeCount = 0;
        for(Int I=1;I<Xsuper_.m();I++){
          Int fc = Xsuper_(I-1);
          Int lc = Xsuper_(I)-1;
          Int iWidth = lc-fc+1;
          Ptr fi = xlindx_(I-1);
          Ptr li = xlindx_(I)-1;
          Int iHeight = li-fi+1;

          Int iDest = this->Mapping_->Map(I-1,I-1);

          if(iDest==iam){
            ++snodeCount;
          }

          //look at the owner of the first column of the supernode
          Int numColFirst = std::max(1,iSize_ / np);

          //post all the recv and sends
          //      logfileptr->OFS()<<"Sending columns of SuperNode2 "<<I<<endl;
          for(Int col = fc;col<=lc;col++){
            //corresponding column in the unsorted matrix A
            Int orig_col = Order_.perm[col-1];
            Int iOwnerCol = std::min((orig_col-1)/numColFirst,np-1);

            if(iam == iDest){
              if(iam != iOwnerCol){
                Int nrows = Global_->expColptr[orig_col] - Global_->expColptr[orig_col-1];
                size_t & recv_bytes = recv_map[iOwnerCol].first;
                recv_bytes += sizeof(Int)+sizeof(Int)+nrows*(sizeof(Int) + sizeof(T));
              }
            }
            else if(iam == iOwnerCol){
              Int local_col = (orig_col-(numColFirst)*iOwnerCol);
              Int nrows = ExpA.Local_.colptr[local_col]-ExpA.Local_.colptr[local_col-1];
              size_t & send_bytes = send_map[iDest].first;
              send_bytes += sizeof(Int)+sizeof(Int)+nrows*(sizeof(Int)+sizeof(T));
            }
          }
        }


        //Resize the local supernodes array
        LocalSupernodes_.reserve(snodeCount);


        
        for(Int I=1;I<Xsuper_.m();I++){

          Int iDest = this->Mapping_->Map(I-1,I-1);

          //parse the first column to create the supernode structure
          if(iam==iDest){
            Int fc = Xsuper_(I-1);
            Int lc = Xsuper_(I)-1;
            Int iWidth = lc-fc+1;
            Ptr fi = xlindx_(I-1);
            Ptr li = xlindx_(I)-1;
            Int iHeight = li-fi+1;

#ifndef ITREE2
            globToLocSnodes_.push_back(I-1);
#else 
            ITree::Interval snode_inter = { I, I, LocalSupernodes_.size() };
            globToLocSnodes_.Insert(snode_inter);
#endif

                  //count number of contiguous blocks (at least one diagonal block)
                  Int nzBlockCnt = 1;
                  Idx iPrevRow = lindx_(fi-1)-1;
                  for(Ptr idx = fi; idx<=li;idx++){
                    Idx iRow = lindx_(idx-1);
    
                    //enforce the first block to be a square diagonal block
                    if(nzBlockCnt==1 && iRow>fc/*lindx_(fi-1)*/+iWidth-1){
                      nzBlockCnt++;
                    }
                    else if(iRow!=iPrevRow+1){
                      nzBlockCnt++;
                    }

                    iPrevRow=iRow;
                  }

            LocalSupernodes_.push_back( new SuperNode2<T>(I,fc,lc,iHeight,iSize_,nzBlockCnt));
          }
        }



        //Create the remote ones when we receive something
        for(Int p=0;p<np;++p){
          if(iam==p){
            //ONE BY ONE            
            //second, pack and alloc send/recv buffer
            for(auto it = send_map.begin(); it!=send_map.end();it++){
              Int iCurDest = it->first;
              size_t & send_bytes = it->second.first;


              Icomm * IsendPtr = new Icomm(send_bytes,MPI_REQUEST_NULL);
              //get the Isend
              Icomm & Isend = *IsendPtr;

              for(Int I=1;I<Xsuper_.m();I++){
                Int fc = Xsuper_(I-1);
                Int lc = Xsuper_(I)-1;
                Int iWidth = lc-fc+1;
                Ptr fi = xlindx_(I-1);
                Ptr li = xlindx_(I)-1;
                Int iHeight = li-fi+1;

                Int iDest = this->Mapping_->Map(I-1,I-1);

                if(iDest == iCurDest){
#ifdef _DEBUG_
                  logfileptr->OFS()<<"Supernode "<<I<<" is owned by P"<<iDest<<std::endl;
#endif

                  //look at the owner of the first column of the supernode
                  Int numColFirst = std::max(1,iSize_ / np);

                  //post all the recv and sends
                  //      logfileptr->OFS()<<"Sending columns of SuperNode2 "<<I<<endl;
                  bool isSend = false;
                  for(Int col = fc;col<=lc;col++){
                    //corresponding column in the unsorted matrix A
                    Int orig_col = Order_.perm[col-1];
                    Int iOwnerCol = std::min((orig_col-1)/numColFirst,np-1);

                    if(iam == iOwnerCol){
                      isSend = true;
                    }
                  }

                  if(isSend){

                    //Serialize
                    //logfileptr->OFS()<<"Sending SuperNode2 "<<I<<" cols {";
                    for(Int col = fc;col<=lc;col++){
                      Int orig_col = Order_.perm[col-1];
                      Int iOwnerCol = std::min((orig_col-1)/numColFirst,np-1);
                      if(iam==iOwnerCol && iDest!=iam){

                        Int local_col = (orig_col-(numColFirst)*iOwnerCol);
                        Int nrows = ExpA.Local_.colptr[local_col]-ExpA.Local_.colptr[local_col-1];
                        Isend<<col;
                        Isend<<nrows;
                        Serialize(Isend,(char*)&ExpA.Local_.rowind[ExpA.Local_.colptr[local_col-1]-1],nrows*sizeof(Int));
                        Serialize(Isend,(char*)&ExpA.nzvalLocal[ExpA.Local_.colptr[local_col-1]-1],nrows*sizeof(T));
                      }
                    }
                  }
                }
              }


              Int tag = p;
              MPI_Send(Isend.front(),Isend.size(),MPI_BYTE,iCurDest,tag,CommEnv_->MPI_GetComm()/*,&Isend.Request*/);



              delete IsendPtr;
            }
          }
          else{
            auto it = recv_map.find(p);
            if(it!=recv_map.end()){

              Int iOwnerCol = it->first;
              size_t & exp_recv_bytes = it->second.first;

              Icomm * Irecv = new Icomm(exp_recv_bytes,MPI_REQUEST_NULL);
              Int tag = p;
              assert(p==iOwnerCol);
              MPI_Status recv_status;
              MPI_Recv(Irecv->front(),Irecv->capacity(),MPI_BYTE,iOwnerCol,tag,CommEnv_->MPI_GetComm(),&recv_status/*,&Irecv.Request*/);

              int recv_bytes;
              MPI_Get_count(&recv_status,MPI_BYTE,&recv_bytes);
              Irecv->setHead(0);

              //logfileptr->OFS()<<"Receiving SuperNode2 "<<I<<" from P"<<recv_status.MPI_SOURCE<<" ("<<recv_bytes<<") cols {";
              while(Irecv->head <Irecv->capacity()){ 
                char * buffer = Irecv->back();

                //Deserialize
                Int col = *(Int*)&buffer[0];

                Int I = SupMembership_[col-1];

                Int iDest = this->Mapping_->Map(I-1,I-1);
                assert(iam==iDest);

                Int fc = Xsuper_(I-1);
                Int lc = Xsuper_(I)-1;
                Int iWidth = lc-fc+1;
                Ptr fi = xlindx_(I-1);
                Ptr li = xlindx_(I)-1;
                Int iHeight = li-fi+1;

                SuperNode2<T> & snode = *snodeLocal(I);

                if(snode.NZBlockCnt()==0){
                  for(Ptr idx = fi; idx<=li;idx++){
                    Idx iStartRow = lindx_(idx-1);
                    Idx iPrevRow = iStartRow;
                    Int iContiguousRows = 1;
                    for(Int idx2 = idx+1; idx2<=li;idx2++){
                      Idx iCurRow = lindx_(idx2-1);
                      if(iStartRow == lindx_(fi-1)){
                        if(iCurRow>iStartRow+iWidth-1){
                          //enforce the first block to be a square diagonal block
                          break;
                        }
                      }

                      if(iCurRow==iPrevRow+1){
                        idx++;
                        ++iContiguousRows;
                        iPrevRow=iCurRow;
                      }
                      else{
                        break;
                      }
                    }

                    Int iCurNZcnt = iContiguousRows * iWidth;
                    snode.AddNZBlock(iContiguousRows , iWidth,iStartRow);
                  }

                  snode.Shrink();
                }

                Int nrows = *((Int*)&buffer[sizeof(Int)]);
                Int * rowind = (Int*)(&buffer[2*sizeof(Int)]);
                T * nzvalA = (T*)(&buffer[(2+nrows)*sizeof(Int)]);
                //advance in the buffer
                Irecv->setHead(Irecv->head + 2*sizeof(Int) + nrows*(sizeof(Int)+sizeof(T)));

                //logfileptr->OFS()<<col<<" ";

                Int orig_col = Order_.perm[col-1];
                Int iOwnerCol = std::min((orig_col-1)/numColFirst,np-1);

                Int colbeg = 1;
                Int colend = nrows;

                for(Int rowidx = colbeg; rowidx<=colend; ++rowidx){
                  Int orig_row = rowind[rowidx-1];
                  Int row = Order_.invp[orig_row-1];

                  if(row>=col){
                    Int blkidx = snode.FindBlockIdx(row);
                    NZBlockDesc2 & blk_desc = snode.GetNZBlockDesc(blkidx);
                    Int local_row = row - blk_desc.GIndex + 1;
                    Int local_col = col - fc + 1;
                    T * nzval = snode.GetNZval(blk_desc.Offset);
                    nzval[(local_row-1)*iWidth+local_col-1] = nzvalA[rowidx-1];
                  }
                }
              }
              //logfileptr->OFS()<<"}"<<endl;



              delete Irecv;


            }

          }
        }


        //do the local and resize my SuperNode2 structure
        for(Int I=1;I<Xsuper_.m();I++){

          Int iDest = this->Mapping_->Map(I-1,I-1);

#ifdef _DEBUG_
          logfileptr->OFS()<<"Supernode "<<I<<" is owned by P"<<iDest<<std::endl;
#endif

          //parse the first column to create the supernode structure
          if(iam==iDest){
            Int fc = Xsuper_(I-1);
            Int lc = Xsuper_(I)-1;
            Int iWidth = lc-fc+1;
            Ptr fi = xlindx_(I-1);
            Ptr li = xlindx_(I)-1;
            Int iHeight = li-fi+1;

//#ifndef ITREE2
//            globToLocSnodes_.push_back(I-1);
//#else 
//            ITree::Interval snode_inter = { I, I, LocalSupernodes_.size() };
//            globToLocSnodes_.Insert(snode_inter);
//#endif
//
//            LocalSupernodes_.push_back( new SuperNode2<T>(I,fc,lc,iHeight,iSize_));
            SuperNode2<T> & snode = *snodeLocal(I);

            if(snode.NZBlockCnt()==0){
              for(Ptr idx = fi; idx<=li;idx++){
                Idx iStartRow = lindx_(idx-1);
                Idx iPrevRow = iStartRow;
                Int iContiguousRows = 1;
                for(Int idx2 = idx+1; idx2<=li;idx2++){
                  Idx iCurRow = lindx_(idx2-1);
                  if(iStartRow == lindx_(fi-1)){
                    if(iCurRow>iStartRow+iWidth-1){
                      //enforce the first block to be a square diagonal block
                      break;
                    }
                  }

                  if(iCurRow==iPrevRow+1){
                    idx++;
                    ++iContiguousRows;
                    iPrevRow=iCurRow;
                  }
                  else{
                    break;
                  }
                }

                Int iCurNZcnt = iContiguousRows * iWidth;
                snode.AddNZBlock(iContiguousRows , iWidth,iStartRow);
              }

              snode.Shrink();
            }






            //do the local
            for(Int col = fc;col<=lc;col++){
              //corresponding column in the unsorted matrix A
              Int orig_col = Order_.perm[col-1];
              Int iOwnerCol = std::min((orig_col-1)/numColFirst,np-1);

              if(iam == iOwnerCol){
                int * colptr = ExpA.Local_.colptr.Data();
                int * rowind = ExpA.Local_.rowind.Data();
                T * nzvalA = ExpA.nzvalLocal.Data();

                Int local_col = (orig_col-(numColFirst)*iOwnerCol);

                Int colbeg = colptr[local_col-1];
                Int colend = colptr[local_col]-1;
                for(Int rowidx = colbeg; rowidx<=colend; ++rowidx){
                  Int orig_row = rowind[rowidx-1];
                  Int row = Order_.invp[orig_row-1];

                  if(row>=col){
                    Int blkidx = snode.FindBlockIdx(row);
                    NZBlockDesc2 & blk_desc = snode.GetNZBlockDesc(blkidx);
                    Int local_row = row - blk_desc.GIndex + 1;
                    Int local_col = col - fc + 1;
                    T * nzval = snode.GetNZval(blk_desc.Offset);
                    nzval[(local_row-1)*iWidth+local_col-1] = nzvalA[rowidx-1];
                  }
                }
              }
            }
          }
        }
      }
      catch(const std::bad_alloc& e){
        std::cout << "Allocation failed: " << e.what() << '\n';
        abort();
      }
      logfileptr->OFS()<<"Send Done"<<endl;

    TIMER_STOP(DISTRIBUTING_MATRIX);

#ifdef _DEBUG_
    for(Int I=1;I<Xsuper_.m();I++){
      Int src_first_col = Xsuper_(I-1);
      Int src_last_col = Xsuper_(I)-1;
      Int iOwner = this->Mapping_->Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam ){
        SuperNode2<T> & src_snode = *snodeLocal(I);


        logfileptr->OFS()<<"+++++++++++++"<<I<<"++++++++"<<std::endl;

        logfileptr->OFS()<<"cols: ";
        for(Int i = 0; i< src_snode.Size(); ++i){
          logfileptr->OFS()<<" "<<Order_.perm[src_first_col+i-1];
        }
        logfileptr->OFS()<<std::endl;
        logfileptr->OFS()<<std::endl;
        for(int blkidx=0;blkidx<src_snode.NZBlockCnt();++blkidx){

          NZBlockDesc2 & desc = src_snode.GetNZBlockDesc(blkidx);
          T * val = src_snode.GetNZval(desc.Offset);
          Int nRows = src_snode.NRows(blkidx);

          Int row = desc.GIndex;
          for(Int i = 0; i< nRows; ++i){
            logfileptr->OFS()<<row+i<<" | "<<Order_.perm[row+i-1]<<":   ";
            for(Int j = 0; j< src_snode.Size(); ++j){
              logfileptr->OFS()<<val[i*src_snode.Size()+j]<<" ";
            }
            logfileptr->OFS()<<std::endl;
          }

          logfileptr->OFS()<<"_______________________________"<<std::endl;
        }
      }
    }
#endif
    Dump();


#ifdef PREALLOC_IRECV
//    TIMER_START(ALLOC_IRECV_BUFFER);
//    //Init asyncRecv buffers
//    auto maxwit = max_element(&UpdateWidth_[0],&UpdateWidth_[0]+UpdateWidth_.m());
//    auto maxhit = max_element(&UpdateHeight_[0],&UpdateHeight_[0]+UpdateHeight_.m());
//
//    int maxh = *maxhit;
//    int maxw = *maxwit;
//    for(Int I=1;I<Xsuper_.m();++I){
//      int width = Xsuper_(I) - Xsuper_(I-1);
//      maxw = max(maxw,width);
//    }
//    
//
//      Int max_bytes = 6*sizeof(Int); 
//      Int nrows = maxh;
//      Int ncols = maxw;
//      Int nz_cnt = nrows * ncols;
//      Int nblocks = nrows;
//      max_bytes += (nblocks)*sizeof(NZBlockDesc2);
//      max_bytes += nz_cnt*sizeof(T);
//
//
//
//    for(int i = 0; i<maxIrecv_;++i){
//      availRecvBuffers_.push_back( new Icomm(max_bytes,MPI_REQUEST_NULL) );
//    }
//    TIMER_STOP(ALLOC_IRECV_BUFFER);
#endif



    //xlindx_.Clear();
    //lindx_.Clear();
    //delete Local_;
    //delete Global_;

    MPI_Barrier(CommEnv_->MPI_GetComm());


    TIMER_STOP(SUPERMATRIX_INIT);

  }




  template <typename T> SupernodalMatrix2<T>::SupernodalMatrix2(){
    CommEnv_=NULL;
    Mapping_ = NULL;
    scheduler_ = NULL;
    Local_=NULL;
    Global_=NULL;
    isGlobStructAllocated_ = false;
    isXlindxAllocated_=false;
    isLindxAllocated_=false;
  }

  template <typename T> SupernodalMatrix2<T>::SupernodalMatrix2(const DistSparseMatrix<T> & pMat, NGCholOptions & options ){
    CommEnv_ = NULL;
    Mapping_ = NULL;
    scheduler_ = NULL;
    Local_=NULL;
    Global_=NULL;
    isGlobStructAllocated_ = false;
    isXlindxAllocated_=false;
    isLindxAllocated_=false;
    Init(pMat, options);
  }

  template <typename T> SupernodalMatrix2<T>::~SupernodalMatrix2(){
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

    if(this->Mapping_!=NULL){
      delete this->Mapping_;
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
  }


//returns the 1-based index of supernode id global in the local supernode array
template <typename T> Int SupernodalMatrix2<T>::snodeLocalIndex(Int global){
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
template <typename T> SuperNode2<T> * SupernodalMatrix2<T>::snodeLocal(Int global){
      Int iLocal = snodeLocalIndex(global);
      return LocalSupernodes_[iLocal -1];
}

template <typename T> SuperNode2<T> * SupernodalMatrix2<T>::snodeLocal(Int global, std::vector<SuperNode2<T> *> & snodeColl){
      Int iLocal = snodeLocalIndex(global);
      return snodeColl[iLocal -1];
}


#include "SupernodalMatrix2_impl_FB_pull.hpp"


}


#endif 
