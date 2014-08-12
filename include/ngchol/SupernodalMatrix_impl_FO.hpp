#ifndef _SUPERNODAL_MATRIX_IMPL_FO_HPP_
#define _SUPERNODAL_MATRIX_IMPL_FO_HPP_

  struct FOUpdate{
    Int tgt_snode_id;
    Int src_first_row;
    Int src_last_row;
    Int src_nzblk_idx;
    Int src_next_nzblk_idx;
    std::vector<bool> is_factor_sent;

    FOUpdate(Int np){
      tgt_snode_id = 0;
      src_first_row = 0;
      src_last_row = 0;
      src_nzblk_idx = 0;
      src_next_nzblk_idx = 0;
      is_factor_sent.assign(np,false);
    }
  };




template <typename T> void SupernodalMatrix<T>::SendDelayedMessages(Int iLocalI, CommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode<T> *> & snodeColl, bool reverse){
//  typedef volatile LIBCHOLESKY::Int Int;
 
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

  if(!MsgToSend.empty()){

#ifdef _DEBUG_DELAY_
    {
      CommList tmp = MsgToSend;
      logfileptr->OFS()<<"Queue : ";
      while( tmp.size()>0){
        //Pull the highest priority message
        const DelayedComm & comm = tmp.top();
        Int src_snode_id = comm.src_snode_id;
        Int tgt_id = comm.tgt_snode_id;
        Int src_nzblk_idx = comm.src_nzblk_idx;
        Int src_first_row = comm.src_first_row;
        logfileptr->OFS()<<" { "<<src_snode_id<<" -> "<<tgt_id<<" }";
        tmp.pop();
      }
      logfileptr->OFS()<<endl;


    }
#endif
    while( MsgToSend.size()>0){


#ifdef _DEBUG_DELAY_
      logfileptr->OFS()<<"Async comms: "<<OutgoingSend<<endl;
#endif
      //Check for completion of outgoing communication
//      AdvanceOutgoing(OutgoingSend);


      //Pull the highest priority message
      const DelayedComm & comm = MsgToSend.top();
      Int src_snode_id = comm.src_snode_id;
      Int tgt_id = comm.tgt_snode_id;
      Int src_nzblk_idx = comm.src_nzblk_idx;
      Int src_first_row = comm.src_first_row;

#ifdef _DEBUG_DELAY_
      logfileptr->OFS()<<"Picked { "<<src_snode_id<<" -> "<<tgt_id<<" }"<<endl;
#endif

      Int iLocalSrc = (src_snode_id-1) / np +1 ;
      SuperNode<T> & prev_src_snode = *snodeColl[iLocalSrc -1];

      //        assert(prev_src_snode.Id()==src_snode_id);

      bool is_less = tgt_id < next_snode_id;
//      bool is_sendable = ((!is_last_local && is_less)|| ( is_after_last_local) || (!is_last_snode && is_last_local));
      bool is_sendable = (is_less || is_last);

      if(is_sendable){

        //this can be sent now
        Int tgt_snode_id = tgt_id;

        Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
        if(iTarget != iam){
#ifdef _DEBUG_
          logfileptr->OFS()<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode.Id()<<" to Supernode "<<tgt_snode_id<<endl;
          cout<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode.Id()<<" to Supernode "<<tgt_snode_id<<endl;
#endif

#ifdef _DEBUG_
          logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated by Supernode "<<prev_src_snode.Id()<<" rows "<<src_first_row/*<<" to "<<src_last_row*/<<" "<<src_nzblk_idx<<std::endl;
#endif

          //Send
          Int tgt_first_col = Xsuper_(tgt_snode_id-1);
          Int tgt_last_col = Xsuper_(tgt_snode_id)-1;
          NZBlockDesc & pivot_desc = prev_src_snode.GetNZBlockDesc(src_nzblk_idx);
          Int local_first_row = src_first_row-pivot_desc.GIndex;
          Int nzblk_cnt = prev_src_snode.NZBlockCnt()-src_nzblk_idx;
          Int nz_cnt = (prev_src_snode.NRowsBelowBlock(src_nzblk_idx)
              - local_first_row )*prev_src_snode.Size();
          //            assert(nz_cnt>0);

          T * nzval_ptr = prev_src_snode.GetNZval(pivot_desc.Offset
              +local_first_row*prev_src_snode.Size());

          TIMER_START(SEND_MALLOC);
          AddOutgoingComm(OutgoingSend, prev_src_snode.Id(), prev_src_snode.Size(), src_first_row, pivot_desc, nzblk_cnt, nzval_ptr, nz_cnt);
          TIMER_STOP(SEND_MALLOC);

          if(OutgoingSend.size() > maxIsend_){
            TIMER_START(SEND_MPI);
            MPI_Send(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm());
            TIMER_STOP(SEND_MPI);

#ifdef _DEBUG_DELAY_
            logfileptr->OFS()<<"DELAYED Sending "<<nz_cnt*sizeof(T)<<" bytes to P"<<iTarget<<std::endl;
            logfileptr->OFS()<<"DELAYED     Send factor "<<prev_src_snode.Id()<<" to node"<<tgt_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
            logfileptr->OFS()<<"DELAYED Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
            logfileptr->OFS()<<"DELAYED Sending "<<OutgoingSend.back()->size()<<" bytes to P"<<iTarget<<std::endl;
#endif
            OutgoingSend.pop_back();

          }
          else{
            TIMER_START(SEND_MPI);
            MPI_Isend(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm(),&OutgoingSend.back()->Request);
            TIMER_STOP(SEND_MPI);

#ifdef _DEBUG_DELAY_
            logfileptr->OFS()<<"DELAYED Sending "<<nz_cnt*sizeof(T)<<" bytes to P"<<iTarget<<std::endl;
            logfileptr->OFS()<<"DELAYED     Send factor "<<prev_src_snode.Id()<<" to node"<<tgt_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
            logfileptr->OFS()<<"DELAYED Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
            logfileptr->OFS()<<"DELAYED Sending "<<OutgoingSend.back()->size()<<" bytes to P"<<iTarget<<std::endl;
#endif
          }






        }
      }
      else{
//        gdb_lock();
      }

      if(is_sendable){
        //remove from the list
        MsgToSend.pop();
        //it = MsgToSend.erase(it);
      }
      else{
        break;
        //it++;
      }
    }


  }
  }

  template<typename T> void SupernodalMatrix<T>::AsyncRecvFactors(Int iLocalI, std::vector<AsyncComms> & incomingRecvArr,IntNumVec & FactorsToRecv,IntNumVec & UpdatesToDo){

    for(Int nextLocalI = iLocalI;nextLocalI<=LocalSupernodes_.size();++nextLocalI){
      SuperNode<T> * next_src_snode = LocalSupernodes_[nextLocalI-1];
      AsyncComms & incomingRecv = incomingRecvArr[nextLocalI-1]; 


      Int IrecvCnt = 0; 
      Int maxRecvCnt = min(UpdatesToDo[nextLocalI-1],FactorsToRecv[nextLocalI-1]);

      for(Int idx =0; idx<maxRecvCnt && incomingRecvCnt_ + IrecvCnt < maxIrecv_;
          ++idx){
        Int max_bytes = 5*sizeof(Int); 
        //The upper bound must be of the width of the "largest" child
        Int nrows = next_src_snode->NRowsBelowBlock(0);
        Int ncols = UpdateWidth_(next_src_snode->Id()-1);
        Int nz_cnt = nrows * ncols;

        Int nblocks = nrows;//std::max((Int)ceil(nrows/2)+1,next_src_snode.NZBlockCnt());
        max_bytes += (nblocks)*sizeof(NZBlockDesc);
        max_bytes += nz_cnt*sizeof(T);

        incomingRecv.push_back(new Icomm(max_bytes,MPI_REQUEST_NULL));
        Icomm & Irecv = *incomingRecv.back();
        MPI_Irecv(Irecv.front(),Irecv.capacity(),MPI_BYTE,MPI_ANY_SOURCE,next_src_snode->Id(),CommEnv_->MPI_GetComm(),&Irecv.Request);
        ++IrecvCnt;
      }

      //  logfileptr->OFS()<<"Posting "<<IrecvCnt<<" IRECV for Supernode "<<next_src_snode->Id()<<std::endl;
      incomingRecvCnt_+=IrecvCnt;
      FactorsToRecv[nextLocalI-1] = maxRecvCnt - IrecvCnt;

      //          assert(FactorsToRecv[nextLocalI-1]>=0);

      if( incomingRecvCnt_ >= maxIrecv_){
        break;
      }
    }
#ifdef _DEBUG_
    logfileptr->OFS()<<"FactorsToRecv: "<<FactorsToRecv<<std::endl;
#endif

  }

  template <typename T> inline AsyncComms::iterator SupernodalMatrix<T>::WaitIncomingFactors(AsyncComms & cur_incomingRecv, MPI_Status & recv_status, AsyncComms & outgoingSend) {
    //        vector<MPI_Request> reqs;
    //        reqs.reserve(cur_incomingRecv.size());
    //      for(AsyncComms::iterator it = cur_incomingRecv.begin(); it!=cur_incomingRecv.end();++it){
    //        reqs.push_back((*it)->Request);
    //      }
    //
    //          TIMER_START(IRECV_MPI2);
    //
    //          AsyncComms::iterator it = cur_incomingRecv.begin();
    //          if(cur_incomingRecv.size()==0){
    //            it = cur_incomingRecv.end();
    //          }
    //          else{
    //          int index=-1;
    //          MPI_Waitany(reqs.size(),&reqs[0],&index,&recv_status);
    //          advance(it,index);
    //          }
    //          TIMER_STOP(IRECV_MPI2);
    //          return it;
    //



    TIMER_START(IRECV_MPI);
    if(cur_incomingRecv.size()==0){
      TIMER_STOP(IRECV_MPI);
      return cur_incomingRecv.end();
    }
    else{
      Int done = 0;
      while(cur_incomingRecv.size()>0){
        Int index = 0;
        for(AsyncComms::iterator it = cur_incomingRecv.begin(); it!=cur_incomingRecv.end();++it, ++index){
          Icomm * curComm = *it;
          //MPI_Request req = (curComm->Request);
          //MPI_Test(&(curComm->Request),&done,&recv_status);
 //     set_mpi_handler(MPI_COMM_WORLD);
      int error_code = MPI_Test(&(curComm->Request),&done,&recv_status);
//      check_mpi_error(error_code, MPI_COMM_WORLD, true);

          //Test if comm is done
          if(done==1){

            Int bytes_received = 0;
            MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
            curComm->setHead(bytes_received);

            TIMER_STOP(IRECV_MPI);
            return it;
          }
        }

        //   AdvanceOutgoing(outgoingSend);
      }
    }
  }








  /// Routines used to perform a FanOut update
#ifdef SINGLE_BLAS
  template <typename T> inline void SupernodalMatrix<T>::UpdateSuperNode(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode, Int &pivot_idx, NumMat<T> & tmpBuf,IntNumVec & src_colindx, IntNumVec & src_rowindx, IntNumVec & src_to_tgt_offset
      , Int  pivot_fr)
  {

    TIMER_START(UPDATE_SNODE);

    TIMER_START(UPDATE_SNODE_FIND_INDEX);
    Int first_pivot_idx = -1;
    Int tgt_fc = pivot_fr;
    if(tgt_fc ==I_ZERO ){
      tgt_fc = tgt_snode.FirstCol();
      //find the pivot idx
      do {first_pivot_idx = src_snode.FindBlockIdx(tgt_fc); tgt_fc++;}
      while(first_pivot_idx<0 && tgt_fc<=tgt_snode.LastCol());
      tgt_fc--;
    }
    else{
      first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
    }
    NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);

    Int tgt_lc = tgt_snode.LastCol();
    Int last_pivot_idx = -1;
    //find the pivot idx
    do {last_pivot_idx = src_snode.FindBlockIdx(tgt_lc); tgt_lc--;}
    while(last_pivot_idx<0 && tgt_lc>=tgt_fc);
    tgt_lc++;
    NZBlockDesc & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);
    TIMER_STOP(UPDATE_SNODE_FIND_INDEX);







    //determine the first column that will be updated in the target supernode
    Int tgt_local_fc =  tgt_fc - tgt_snode.FirstCol();
    Int tgt_local_lc =  tgt_lc - tgt_snode.FirstCol();

    Int tgt_nrows = tgt_snode.NRowsBelowBlock(0);
    Int src_nrows = src_snode.NRowsBelowBlock(first_pivot_idx) - (tgt_fc - first_pivot_desc.GIndex);
    Int src_lr = tgt_fc+src_nrows-1;
    src_nrows = src_lr - tgt_fc + 1;

    Int tgt_width = src_snode.NRowsBelowBlock(first_pivot_idx) - (tgt_fc - first_pivot_desc.GIndex) - src_snode.NRowsBelowBlock(last_pivot_idx) + (tgt_lc - last_pivot_desc.GIndex)+1;



    T * pivot = &(src_snode.GetNZval(first_pivot_desc.Offset)[(tgt_fc-first_pivot_desc.GIndex)*src_snode.Size()]);
    T * tgt = tgt_snode.GetNZval(0);

    //If the target supernode has the same structure,
    //The GEMM is directly done in place
    if(src_nrows == tgt_nrows){
      //if(src_snode.NRowsBelowBlock(last_pivot_idx+1) == tgt_snode.NRowsBelowBlock(1) && first_pivot_idx == last_pivot_idx)
      Int tgt_offset = (tgt_fc - tgt_snode.FirstCol());
      tgt = &tgt[tgt_offset];

      //everything is in row-major
      TIMER_START(UPDATE_SNODE_GEMM);
      blas::Gemm('T','N',tgt_width, src_nrows,src_snode.Size(),MINUS_ONE<T>(),pivot,src_snode.Size(),pivot,src_snode.Size(),ONE<T>(),tgt,tgt_width);
      TIMER_STOP(UPDATE_SNODE_GEMM);
    }
    else{
      //Compute the update in a temporary buffer
#ifdef _DEBUG_
      tmpBuf.Resize(tgt_width,src_nrows);
#endif

      T * buf = tmpBuf.Data();

      //everything is in row-major
      TIMER_START(UPDATE_SNODE_GEMM);
      blas::Gemm('T','N',tgt_width, src_nrows,src_snode.Size(),MINUS_ONE<T>(),pivot,src_snode.Size(),pivot,src_snode.Size(),ZERO<T>(),buf,tgt_width);
      TIMER_STOP(UPDATE_SNODE_GEMM);

#ifdef _DEBUG_
      logfileptr->OFS()<<"tmpBuf is "<<tmpBuf<<std::endl;
#endif

      //now add the update to the target supernode
      TIMER_START(UPDATE_SNODE_INDEX_MAP);
      Int src_snode_size = src_snode.Size();
      Int tgt_snode_size = tgt_snode.Size();



      if(tgt_snode_size==1){

        Int rowidx = 0;
        for(Int blkidx = first_pivot_idx; blkidx < src_snode.NZBlockCnt(); ++blkidx){
          NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
          Int cur_src_nrows = src_snode.NRows(blkidx);
          Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
          Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);

          //NOTE: Except for the last pivot block which MIGHT 
          //      be splitted onto multiple blocks in the target
          //      The others MUST reside into single target block
          Int row = cur_src_fr;
          while(row<=cur_src_lr){
            Int tgt_blk_idx = tgt_snode.FindBlockIdx(row);
assert(tgt_blk_idx!=-1);
            NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
            Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + tgt_snode.NRows(tgt_blk_idx)-1);
            Int tgtOffset = cur_tgt_desc.Offset + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
            for(Int cr = row ;cr<=lr;++cr){
              tgt[tgtOffset + (cr - row)*tgt_snode_size] += buf[rowidx]; 
              rowidx++;
            }
            row += (lr-row+1);
          }
        }

      }
      else{
        src_colindx.Resize(tgt_width);
        src_to_tgt_offset.Resize(src_nrows);

        Int colidx = 0;
        Int rowidx = 0;
        Int offset = 0;


        for(Int blkidx = first_pivot_idx; blkidx < src_snode.NZBlockCnt(); ++blkidx){
          NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
          Int cur_src_nrows = src_snode.NRows(blkidx);
          Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
          Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
          cur_src_nrows = cur_src_lr - cur_src_fr +1;

          //Except for the last pivot block which MIGHT be splitted onto multiple blocks in the target
          //The other one MUST reside into a single block in the target
          Int row = cur_src_fr;
          while(row<=cur_src_lr){
            Int tgt_blk_idx = tgt_snode.FindBlockIdx(row);
            NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
            Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + tgt_snode.NRows(tgt_blk_idx)-1);
            Int tgtOffset = cur_tgt_desc.Offset + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
            for(Int cr = row ;cr<=lr;++cr){
              if(cr<=tgt_lc){
                src_colindx[colidx++] = cr;
              }
              offset+=tgt_width;
              src_to_tgt_offset[rowidx] = tgtOffset + (cr - row)*tgt_snode_size;
              rowidx++;
            }
            row += (lr-row+1);
          }
        }


        //Multiple cases to consider
        // same structure between src and tgt : src_nrows == tgt_nrows
        // tgt has only one column 
        // single pivot block first_pivot idx == last_pivot_idx updating contiguous columns
        // full sparse case (done right now)
        TIMER_STOP(UPDATE_SNODE_INDEX_MAP);
        /////    TIMER_START(UPDATE_SNODE_ASSEMBLY);

        if(first_pivot_idx==last_pivot_idx){
          Int tgt_offset = (tgt_fc - tgt_snode.FirstCol());
          for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
            blas::Axpy(tgt_width,ONE<T>(),&buf[rowidx*tgt_width],1,&tgt[src_to_tgt_offset[rowidx] + tgt_offset],1);
          }
        }
        else{
          for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
            for(Int colidx = 0; colidx< src_colindx.m();++colidx){
              Int col = src_colindx[colidx];
              Int tgt_colidx = col - tgt_snode.FirstCol();
              tgt[src_to_tgt_offset[rowidx] + tgt_colidx] += buf[rowidx*tgt_width+colidx]; 
            }
          }
        }

        /////    TIMER_STOP(UPDATE_SNODE_ASSEMBLY);
      }
    }
    TIMER_STOP(UPDATE_SNODE);
  }






#else
  template <typename T> inline void SupernodalMatrix<T>::UpdateSuperNode(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode, Int &pivot_idx, Int  pivot_fr) {

    TIMER_START(UPDATE_SNODE);

    NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(pivot_idx);
    Int first_pivot_fr = pivot_fr;
    if(first_pivot_fr ==I_ZERO ){
      first_pivot_fr = first_pivot_desc.GIndex;
    }

    //start with the first pivot
    for(int cur_piv_idx=pivot_idx;cur_piv_idx<src_snode.NZBlockCnt();++cur_piv_idx){


      NZBlockDesc & pivot_desc = src_snode.GetNZBlockDesc(cur_piv_idx);

      if(pivot_fr ==I_ZERO || cur_piv_idx != pivot_idx){
        pivot_fr = pivot_desc.GIndex;
      }

#ifdef _DEBUG_
      assert(pivot_fr >= pivot_desc.GIndex);
#endif

      if(pivot_fr>tgt_snode.LastCol()){
        break;
      }

      Int pivot_nrows = src_snode.NRows(cur_piv_idx);
      Int pivot_lr = min(pivot_desc.GIndex + pivot_nrows -1, tgt_snode.LastCol());
      T * pivot = &(src_snode.GetNZval(pivot_desc.Offset)[(pivot_fr-pivot_desc.GIndex)*src_snode.Size()]);

      //determine the first column that will be updated in the target supernode
      Int tgt_updated_fc =  pivot_fr - tgt_snode.FirstCol();
      Int tgt_updated_lc =  pivot_lr - tgt_snode.FirstCol();

      for(int src_idx=pivot_idx;src_idx<src_snode.NZBlockCnt();++src_idx){

        NZBlockDesc & src_desc = src_snode.GetNZBlockDesc(src_idx);

        Int src_fr = max(first_pivot_fr, src_desc.GIndex) ;
        //Int src_fr = src_desc.GIndex;
        Int src_nrows = src_snode.NRows(src_idx);
        Int src_lr = src_desc.GIndex+src_nrows-1;

        do{
          //TODO Need to be replaced by GlobToLoc index
          Int tgt_idx = tgt_snode.FindBlockIdx(src_fr);

          if(tgt_idx<0){
            break;
          }

#ifdef _DEBUG_
          assert(tgt_idx!=-1 && tgt_idx<tgt_snode.NZBlockCnt());
#endif

          NZBlockDesc & tgt_desc = tgt_snode.GetNZBlockDesc(tgt_idx);
          Int tgt_nrows = tgt_snode.NRows(tgt_idx);
          Int tgt_fr = tgt_desc.GIndex;
          Int tgt_lr = tgt_desc.GIndex+tgt_nrows-1;

          Int update_fr = max(tgt_fr, src_fr);
          Int update_lr = min(tgt_lr, src_lr);

#ifdef _DEBUG_
          assert(update_fr >= tgt_fr); 
          assert(update_lr >= update_fr); 
          assert(update_lr <= tgt_lr); 
#endif


#ifdef _DEBUG_
          logfileptr->OFS()<<"L("<<update_fr<<".."<<update_lr<<","<<tgt_snode.FirstCol()+tgt_updated_fc<<".."<< tgt_snode.FirstCol()+tgt_updated_lc <<") -= L("<<update_fr<<".."<<update_lr<<",:) * L("<<pivot_fr<<".."<<pivot_lr<<",:)'"<<endl;
#endif



          //Update tgt_nzblk with src_nzblk

          T * src = &(src_snode.GetNZval(src_desc.Offset)[(update_fr - src_desc.GIndex)*src_snode.Size()]);
          T * tgt = &(tgt_snode.GetNZval(tgt_desc.Offset)[(update_fr - tgt_desc.GIndex)*tgt_snode.Size()+tgt_updated_fc]);

          //everything is in row-major
          TIMER_START(UPDATE_SNODE_GEMM);
          blas::Gemm('T','N',pivot_lr-pivot_fr+1, update_lr - update_fr + 1,src_snode.Size(),MINUS_ONE<T>(),pivot,src_snode.Size(),src,src_snode.Size(),ONE<T>(),tgt,tgt_snode.Size());
          TIMER_STOP(UPDATE_SNODE_GEMM);

          if(tgt_idx+1<tgt_snode.NZBlockCnt()){
            NZBlockDesc & next_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_idx+1);
            src_fr = next_tgt_desc.GIndex;
          }
          else{
            break;
          }
        }while(src_fr<=src_lr);
      }
    } //end for pivots

    TIMER_STOP(UPDATE_SNODE);
  }
#endif





template <typename T> void SupernodalMatrix<T>::FanOut( ){

  TIMER_START(FACTORIZATION_FO);

  Real timeSta, timeEnd;
  timeSta =  get_time( );


  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();
  IntNumVec UpdatesToDo = UpdateCount_;

  CommList FactorsToSend; 
  AsyncComms outgoingSend;

  std::vector<AsyncComms> incomingRecvArr(LocalSupernodes_.size());
  incomingRecvCnt_ = 0;
  IntNumVec FactorsToRecv(LocalSupernodes_.size());

  std::vector<T> src_nzval;
  std::vector<char> src_blocks;

  std::vector<std::queue<SnodeUpdate> > LocalUpdates(LocalSupernodes_.size());

  Int maxwidth = 0;
  for(Int i = 1; i<Xsuper_.m(); ++i){
    Int width =Xsuper_(i) - Xsuper_(i-1);
    if(width>=maxwidth){
      maxwidth = width;
    }
  }
  TempUpdateBuffers<T> tmpBufs(Size(),maxwidth);



  for(Int i = LocalSupernodes_.size()-1; i>=0; --i){
    SuperNode<T> & cur_snode = *LocalSupernodes_[i];

    SnodeUpdate curUpdate;
    while(cur_snode.FindNextUpdate(curUpdate,Xsuper_,SupMembership_)){ 
      Int iTarget = Mapping_.Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
      if(iTarget == iam){
        Int iLocalJ = (curUpdate.tgt_snode_id-1) / np +1 ;
        LocalUpdates[iLocalJ-1].push((SnodeUpdate)curUpdate);

#ifdef _DEBUG_
        logfileptr->OFS()<<"FUTURE LOCAL Supernode "<<curUpdate.tgt_snode_id<<" is going to be updated by Supernode "<<cur_snode.Id()<<" from row "<<curUpdate.src_first_row<<" "<<curUpdate.blkidx<<std::endl;
#endif
      }
    }
    FactorsToRecv[i] = UpdatesToDo(cur_snode.Id()-1) % np;//UpdatesToDo(cur_snode.Id()-1) - LocalUpdates[i].size();
  }


#ifdef _DEBUG_
  logfileptr->OFS()<<"FactorsToRecv: "<<FactorsToRecv<<std::endl;
#endif



  //FanOut looking cholesky factorization
  Int I =1;
  Int iLocalI=1;
  while(iLocalI<=LocalSupernodes_.size() || !FactorsToSend.empty() || !outgoingSend.empty()){

    //Check for completion of outgoing communication
    AdvanceOutgoing(outgoingSend);

    if(iLocalI>0 && iLocalI<=LocalSupernodes_.size()){
      SuperNode<T> & src_snode = *LocalSupernodes_[iLocalI -1];
      I = src_snode.Id();
      Int src_first_col = src_snode.FirstCol();
      Int src_last_col = src_snode.LastCol();

#ifdef _DEBUG_
      logfileptr->OFS()<<"Processing Supernode "<<I<<std::endl;
#endif

      //Launch Irecv for subsequent local supernodes if I can
      AsyncRecvFactors(iLocalI,incomingRecvArr,FactorsToRecv,UpdatesToDo);

      //Do all my updates (Local and remote)
      //Local updates
      while(!LocalUpdates[iLocalI-1].empty()){
        SnodeUpdate & curUpdate = LocalUpdates[iLocalI-1].front();
        SuperNode<T> & local_src_snode = *LocalSupernodes_[(curUpdate.src_snode_id-1) / np];

#ifdef TRACK_PROGRESS
  Real timeStaTask =  get_time( );
#endif
        src_snode.Update(local_src_snode,curUpdate, tmpBufs);
#ifdef TRACK_PROGRESS
  timeEnd =  get_time( );
  progressptr->OFS()<<curUpdate.tgt_snode_id<<" "<<timeStaTask - timeSta<<" "<<timeEnd - timeSta<<" U"<<curUpdate.src_snode_id<<"-"<<curUpdate.tgt_snode_id<<" "<<iam<< endl;
#endif


        LocalUpdates[iLocalI-1].pop();
        --UpdatesToDo(I-1);
#ifdef _DEBUG_
        logfileptr->OFS()<<"LOCAL Supernode "<<src_snode.Id()<<" is updated by Supernode "<<curUpdate.src_snode_id<<" from row "<<curUpdate.src_first_row<<" "<<curUpdate.blkidx<<std::endl;
        logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left"<<endl;
#endif
      }

      //Remote updates

      //first wait for the Irecv
      AsyncComms & cur_incomingRecv = incomingRecvArr[iLocalI-1];
      MPI_Status recv_status;

      AsyncComms::iterator it = WaitIncomingFactors(cur_incomingRecv, recv_status,outgoingSend);
      while( it != cur_incomingRecv.end() ){
        Icomm * curComm = *it;

        SuperNode<T> dist_src_snode;
        Deserialize(curComm->front(),dist_src_snode);
#ifdef _DEBUG_
        logfileptr->OFS()<<"IRECV Supernode "<<dist_src_snode.Id()<<std::endl;
#endif

        //Update everything I own with that factor
        //Update the ancestors

        TIMER_START(UPDATE_ANCESTORS);
        SnodeUpdate curUpdate;
        while(dist_src_snode.FindNextUpdate(curUpdate,Xsuper_,SupMembership_,false)){ 
          Int iTarget = Mapping_.Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
          if(iTarget == iam){
            Int iLocalJ = (curUpdate.tgt_snode_id-1) / np +1 ;
            SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];

#ifdef TRACK_PROGRESS
  Real timeStaTask =  get_time( );
#endif
            tgt_snode.Update(dist_src_snode,curUpdate, tmpBufs);
#ifdef TRACK_PROGRESS
  timeEnd =  get_time( );
  progressptr->OFS()<<curUpdate.tgt_snode_id<<" "<<timeStaTask - timeSta<<" "<<timeEnd - timeSta<<" U"<<curUpdate.src_snode_id<<"-"<<curUpdate.tgt_snode_id<<" "<<iam<< endl;
#endif



            --UpdatesToDo(curUpdate.tgt_snode_id-1);

#ifdef _DEBUG_
            logfileptr->OFS()<<"IRECV Supernode "<<curUpdate.tgt_snode_id<<" is updated by Supernode "<<dist_src_snode.Id()<<" rows "<<curUpdate.src_first_row<<" "<<curUpdate.blkidx<<std::endl;
            logfileptr->OFS()<<UpdatesToDo(curUpdate.tgt_snode_id-1)<<" updates left for Supernode "<<curUpdate.tgt_snode_id<<endl;
#endif

          }
        }
        TIMER_STOP(UPDATE_ANCESTORS);

        //delete the request from the list
        cur_incomingRecv.erase(it);
        --incomingRecvCnt_;

#ifdef _DEBUG_
        logfileptr->OFS()<<cur_incomingRecv.size()<<" async recv to do for Supernode "<<I<<endl;
        logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left for Supernode "<<I<<endl;
#endif
        if(UpdatesToDo(src_snode.Id()-1)==0){
          //cancel all requests
          incomingRecvCnt_-=cur_incomingRecv.size(); 
          cur_incomingRecv.clear();
        }
        it = WaitIncomingFactors(cur_incomingRecv,recv_status,outgoingSend);
      }




      src_blocks.resize(0);

      Int nz_cnt;
      Int max_bytes;
      while(UpdatesToDo(I-1)>0){
#ifdef _DEBUG_
        logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left"<<endl;
#endif


        TIMER_START(RECV_MALLOC);
        if(src_blocks.size()==0){
          max_bytes = 5*sizeof(Int); 
          //The upper bound must be of the width of the "largest" child
#ifdef _DEBUG_
          logfileptr->OFS()<<"Maximum width is "<<UpdateWidth_(I-1)<<std::endl;
#endif

          Int nrows = src_snode.NRowsBelowBlock(0);
          Int ncols = UpdateWidth_(I-1);
          nz_cnt = nrows * ncols;

          Int nblocks = nrows;//std::max((Int)ceil(nrows/2)+1,src_snode.NZBlockCnt());
          max_bytes += (nblocks)*sizeof(NZBlockDesc);
          max_bytes += nz_cnt*sizeof(T); 

          src_blocks.resize(max_bytes);
        }
        TIMER_STOP(RECV_MALLOC);

        TIMER_START(RECV_MPI);
        MPI_Status recv_status;
        int bytes_received = 0;

#ifdef PROBE_FIRST
        MPI_Probe(MPI_ANY_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
        MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);

        bool doabort = false;
        int prev_size = 0;
        if(src_blocks.size()<bytes_received){
          prev_size = src_blocks.size();
          doabort = true;
          //receive anyway
          src_blocks.resize(bytes_received);
        }
#endif


#ifdef PROBE_FIRST
        MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,recv_status.MPI_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
#else
        //receive the index array
        MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,MPI_ANY_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
#endif

        SuperNode<T> dist_src_snode;
        Deserialize(&src_blocks[0],dist_src_snode);
#ifdef PROBE_FIRST
        if(doabort){
          cout<<"We have a problem !!!! on P"<<iam<<"\n";
          gdb_lock();

          abort();
        }
#endif
#ifdef _DEBUG_
        logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode.Id()<<std::endl;
#endif



        SnodeUpdate curUpdate;
        TIMER_START(UPDATE_ANCESTORS);
        while(dist_src_snode.FindNextUpdate(curUpdate,Xsuper_,SupMembership_,false)){ 
          Int iTarget = Mapping_.Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
          if(iTarget == iam){
            Int iLocalJ = (curUpdate.tgt_snode_id-1) / np +1 ;
            SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];
#ifdef TRACK_PROGRESS
  Real timeStaTask =  get_time( );
#endif
            tgt_snode.Update(dist_src_snode,curUpdate, tmpBufs);
#ifdef TRACK_PROGRESS
  timeEnd =  get_time( );
  progressptr->OFS()<<curUpdate.tgt_snode_id<<" "<<timeStaTask - timeSta<<" "<<timeEnd - timeSta<<" U"<<curUpdate.src_snode_id<<"-"<<curUpdate.tgt_snode_id<<" "<<iam<< endl;
#endif


            --UpdatesToDo(curUpdate.tgt_snode_id-1);

#if defined(_DEBUG_) || defined(_DEBUG_DELAY_)
            logfileptr->OFS()<<"RECV Supernode "<<curUpdate.tgt_snode_id<<" is updated by Supernode "<<dist_src_snode.Id()<<" rows "<<curUpdate.src_first_row<<" "<<curUpdate.blkidx<<std::endl;
            logfileptr->OFS()<<UpdatesToDo(curUpdate.tgt_snode_id-1)<<" updates left for Supernode "<<curUpdate.tgt_snode_id<<endl;
#endif
          }
        }
        TIMER_STOP(UPDATE_ANCESTORS);
      }
      //clear the buffer
      //        { vector<char>().swap(src_blocks);  }
      //        { vector<T>().swap(src_nzval);  }


      timeEnd =  get_time( );
#ifdef _DEBUG_
      if(UpdatesToDo(src_snode.Id()-1)!=0){gdb_lock();}
      assert(UpdatesToDo(src_snode.Id()-1)==0);
      logfileptr->OFS()<<"  Factoring Supernode "<<I<<" at "<<timeEnd-timeSta<<" s"<<std::endl;
#endif

#ifdef TRACK_PROGRESS
  Real timeStaTask =  get_time( );
#endif
      TIMER_START(FACTOR_PANEL);
      src_snode.Factorize();
      TIMER_STOP(FACTOR_PANEL);

#ifdef TRACK_PROGRESS
  timeEnd =  get_time( );
  progressptr->OFS()<<src_snode.Id()<<" "<<timeStaTask - timeSta<<" "<<timeEnd - timeSta<<" F"<<src_snode.Id()<<" "<<iam<<endl;
#endif





      //Send my factor to my ancestors. 
      BolNumVec is_factor_sent(np);
      SetValue(is_factor_sent,false);

      BolNumVec is_skipped(np);
      SetValue(is_skipped,false);

      SnodeUpdate curUpdate;
      TIMER_START(FIND_UPDATED_ANCESTORS);
      while(src_snode.FindNextUpdate(curUpdate,Xsuper_,SupMembership_)){ 
        Int iTarget = Mapping_.Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);

        if(iTarget != iam){

#ifdef _DEBUG_
          logfileptr->OFS()<<"Remote Supernode "<<curUpdate.tgt_snode_id<<" is updated by Supernode "<<I<<" rows "<<curUpdate.src_first_row<<" "<<curUpdate.blkidx<<std::endl;

          if(is_factor_sent[iTarget]){
            logfileptr->OFS()<<"Remote Supernode "<<curUpdate.tgt_snode_id<<" is updated Supernode "<<I<<" in factor "<<is_factor_sent[iTarget]<<std::endl;
          }
#endif

          if(!is_factor_sent[iTarget] && !is_skipped[iTarget] ){

            //need a std::unordered_set to check whether 
            Int next_local_snode = (iLocalI < LocalSupernodes_.size())?LocalSupernodes_[iLocalI]->Id():Xsuper_.m();
            if( next_local_snode< curUpdate.tgt_snode_id){
              //need to push the prev src_last_row
              FactorsToSend.push(DelayedComm(src_snode.Id(),curUpdate.tgt_snode_id,curUpdate.blkidx,curUpdate.src_first_row));
#ifdef _DEBUG_DELAY_
              logfileptr->OFS()<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<curUpdate.tgt_snode_id<<" from row "<<curUpdate.src_first_row<<endl;
              cout<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<curUpdate.tgt_snode_id<<" from row "<<curUpdate.src_first_row<<endl;
#endif
              is_skipped[iTarget] = true;
            }
            if(!is_skipped[iTarget]){
              is_factor_sent[iTarget] = true;

              //Send
              NZBlockDesc & pivot_desc = src_snode.GetNZBlockDesc(curUpdate.blkidx);

              Icomm * send_buffer = new Icomm();
              Serialize(*send_buffer,src_snode,curUpdate.blkidx,curUpdate.src_first_row);
              AddOutgoingComm(outgoingSend,send_buffer);

              if( outgoingSend.size() > maxIsend_){
                TIMER_START(SEND_MPI);
                MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,curUpdate.tgt_snode_id,CommEnv_->MPI_GetComm());
                TIMER_STOP(SEND_MPI);
                outgoingSend.pop_back();
              }
              else{
                TIMER_START(SEND_MPI);
                MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,curUpdate.tgt_snode_id,CommEnv_->MPI_GetComm(),&outgoingSend.back()->Request);
                TIMER_STOP(SEND_MPI);
              }
            }

          }
        }
      }
      TIMER_STOP(FIND_UPDATED_ANCESTORS);
    }


    //process some of the delayed send
    SendDelayedMessagesUp(iLocalI,FactorsToSend,outgoingSend,LocalSupernodes_);
    iLocalI++;
  }


  for(Int idx = 0; idx <incomingRecvArr.size();++idx){
    assert(incomingRecvArr[idx].size()==0);
  } 

  MPI_Barrier(CommEnv_->MPI_GetComm());

  TIMER_STOP(FACTORIZATION_FO);
}


  template <typename T> void SupernodalMatrix<T>::FanOut2( ){

    TIMER_START(FACTORIZATION_FO);

    Real timeSta, timeEnd;
    timeSta =  get_time( );


    Int iam = CommEnv_->MPI_Rank();
    Int np  = CommEnv_->MPI_Size();


    IntNumVec UpdatesToDo = UpdateCount_;
    AsyncComms outgoingSend;

    IntNumVec FactorsToRecv(LocalSupernodes_.size());
    incomingRecvCnt_ = 0;


    std::vector<FOUpdate> nextUpdate(LocalSupernodes_.size(),FOUpdate(np));
    std::vector<AsyncComms> incomingRecvArr(LocalSupernodes_.size());

    std::vector<T> src_nzval;
    std::vector<char> src_blocks;


//    std::vector<std::queue<LocalUpdate> > LocalUpdates(LocalSupernodes_.size());

    for(Int i = LocalSupernodes_.size()-1; i>=0; --i){
      SuperNode<T> & cur_snode = *LocalSupernodes_[i];

      Int tgt_snode_id = 0;
      Int src_first_row = 0;
      Int src_last_row = 0;
      Int src_nzblk_idx = 0;
      Int src_next_nzblk_idx = 0;

      while(FindNextUpdate(cur_snode, tgt_snode_id, src_first_row, src_nzblk_idx, src_last_row, src_next_nzblk_idx)){ 
        Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
        if(iTarget == iam){
          Int iLocalJ = (tgt_snode_id-1) / np +1 ;
//          LocalUpdates[iLocalJ-1].push(LocalUpdate(cur_snode.Id(),src_nzblk_idx,src_first_row));
//
//#ifdef _DEBUG_
//          logfileptr->OFS()<<"FUTURE LOCAL Supernode "<<tgt_snode_id<<" is going to be updated by Supernode "<<cur_snode.Id()<<" from row "<<src_first_row<<" "<<src_nzblk_idx<<std::endl;
//#endif
          //        FactorsToRecv[iLocalJ-1]-=UpdatesToDo(cur_snode.Id()-1)+1;
        }
      }

      FactorsToRecv[i] = UpdatesToDo(cur_snode.Id()-1) % np;//UpdatesToDo(cur_snode.Id()-1) - LocalUpdates[i].size();
    }


#ifdef _DEBUG_
    logfileptr->OFS()<<"FactorsToRecv: "<<FactorsToRecv<<std::endl;
#endif




    Int maxwidth = 0;
    for(Int i = 1; i<Xsuper_.m(); ++i){
      Int width =Xsuper_(i) - Xsuper_(i-1);
      if(width>=maxwidth){
        maxwidth = width;
      }
    }

    NumMat<T> tmpBuf(iSize_,maxwidth);
    IntNumVec src_colindx(maxwidth);
    IntNumVec src_rowindx(iSize_);
    IntNumVec src_to_tgt_offset(iSize_);


    //dummy right looking cholesky factorization
    Int I =1;
    Int iLocalI = 1;
    bool UpdatesAllSent = false;
    while(I<Xsuper_.m() || /*!FactorsToSend.empty() ||*/ !outgoingSend.empty()){

      //Check for completion of outgoing communication
      AdvanceOutgoing(outgoingSend);

      //process some of the delayed send
//      SendDelayedMessages(I,FactorsToSend,outgoingSend);


      
      if(I<Xsuper_.m()){
        Int src_first_col = Xsuper_(I-1);
        Int src_last_col = Xsuper_(I)-1;
        Int iOwner = Mapping_.Map(I-1,I-1);
        //If I own the column, factor it
        if( iOwner == iam ){

          for(Int iLocalJ=1;iLocalJ<iLocalI;++iLocalJ){

            SuperNode<T> & src_snode = *LocalSupernodes_[iLocalJ -1];
            FOUpdate & curUpdate = nextUpdate[iLocalJ-1];


            //Send my factor to my ancestors. 

            Int tgt_snode_id = curUpdate.tgt_snode_id;
            Int src_nzblk_idx = curUpdate.src_nzblk_idx;
            Int src_first_row = curUpdate.src_first_row;
            Int src_last_row = curUpdate.src_last_row;
            Int src_next_nzblk_idx = curUpdate.src_next_nzblk_idx;

            TIMER_START(FIND_UPDATED_ANCESTORS);
            while(FindNextUpdate(src_snode, tgt_snode_id, src_first_row, src_nzblk_idx, src_last_row, src_next_nzblk_idx)){ 


              if(tgt_snode_id>I){
                 break;
              }

              curUpdate.tgt_snode_id = tgt_snode_id;
              curUpdate.src_nzblk_idx = src_nzblk_idx;
              curUpdate.src_first_row = src_first_row;
              curUpdate.src_last_row = src_last_row;
              curUpdate.src_next_nzblk_idx = src_next_nzblk_idx;


              Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
              if(iTarget != iam){
#ifdef _DEBUG_
                logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated by Supernode "<<src_snode.Id()<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;
                if(curUpdate.is_factor_sent[iTarget]){
                  logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated Supernode "<<src_snode.Id()<<" in factor "<<curUpdate.is_factor_sent[iTarget]<<std::endl;
                }
#endif

                if(!curUpdate.is_factor_sent[iTarget]){
                  curUpdate.is_factor_sent[iTarget] = true;

                  Int tgt_first_col = Xsuper_(tgt_snode_id-1);
                  Int tgt_last_col = Xsuper_(tgt_snode_id)-1;

                  //Send
                  NZBlockDesc & pivot_desc = src_snode.GetNZBlockDesc(src_nzblk_idx);

                  Int local_first_row = src_first_row-pivot_desc.GIndex;
                  Int nzblk_cnt = src_snode.NZBlockCnt()-src_nzblk_idx;

                  Int nz_cnt = (src_snode.NRowsBelowBlock(src_nzblk_idx) - local_first_row )*src_snode.Size();
                  //              assert(nz_cnt>0);

                  T * nzval_ptr = src_snode.GetNZval(pivot_desc.Offset
                      +local_first_row*src_snode.Size());


                  TIMER_START(SEND_MALLOC);
                  AddOutgoingComm(outgoingSend, src_snode.Id(), src_snode.Size(), src_first_row, pivot_desc, nzblk_cnt, nzval_ptr, nz_cnt);
                  TIMER_STOP(SEND_MALLOC);

                  if(outgoingSend.size() > maxIsend_){
                    TIMER_START(SEND_MPI);
                    MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm());
                    TIMER_STOP(SEND_MPI);

#ifdef _DEBUG_            
                    logfileptr->OFS()<<"Sending "<<nz_cnt*sizeof(T)<<" bytes to P"<<iTarget<<std::endl;
                    logfileptr->OFS()<<"     Sent factor "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
                    logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
                    logfileptr->OFS()<<"Sending "<<outgoingSend.back()->size()<<" bytes to P"<<iTarget<<std::endl;
#endif
                    outgoingSend.pop_back();
                  }
                  else{
                    TIMER_START(SEND_MPI);
                    MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm(),&outgoingSend.back()->Request);
                    TIMER_STOP(SEND_MPI);

#ifdef _DEBUG_            
                    logfileptr->OFS()<<"Sending "<<nz_cnt*sizeof(T)<<" bytes to P"<<iTarget<<std::endl;
                    logfileptr->OFS()<<"     Sent factor "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
                    logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
                    logfileptr->OFS()<<"Sending "<<outgoingSend.back()->size()<<" bytes to P"<<iTarget<<std::endl;
#endif
                  }
                }
              }
              else{
                //do we do the update now ?

                Int iLocalTgt = (tgt_snode_id-1) / np +1 ;
                SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalTgt -1];
//                assert(tgt_snode_id == tgt_snode.Id());
#ifdef _DEBUG_
                logfileptr->OFS()<<"LOCAL Supernode "<<tgt_snode_id<<" is updated by Supernode "<<src_snode.Id()<<" from row "<<src_first_row<<" "<<src_nzblk_idx<<std::endl;
#endif

#ifdef SINGLE_BLAS
                UpdateSuperNode(src_snode,tgt_snode,src_nzblk_idx,tmpBuf,src_colindx,src_rowindx,src_to_tgt_offset, src_first_row);
#else
                UpdateSuperNode(src_snode,tgt_snode,src_nzblk_idx, src_first_row);
#endif
                --UpdatesToDo(tgt_snode_id-1);
#ifdef _DEBUG_
                logfileptr->OFS()<<UpdatesToDo(tgt_snode_id-1)<<" updates left"<<endl;
#endif
              }





            }
            TIMER_STOP(FIND_UPDATED_ANCESTORS);
          } 
















#ifdef _DEBUG_
          logfileptr->OFS()<<"Processing Supernode "<<I<<std::endl;
#endif

          SuperNode<T> & src_snode = *LocalSupernodes_[iLocalI -1];



          //Launch Irecv for subsequent local supernodes if I can
          AsyncRecvFactors(iLocalI,incomingRecvArr,FactorsToRecv,UpdatesToDo);


          //Remote updates

          //first wait for the Irecv
          AsyncComms & cur_incomingRecv = incomingRecvArr[iLocalI-1];
          MPI_Status recv_status;

          AsyncComms::iterator it = WaitIncomingFactors(cur_incomingRecv, recv_status,outgoingSend);
          while( it != cur_incomingRecv.end() ){
            Icomm * curComm = *it;


            std::vector<char> & src_blocks = *curComm->pSrcBlocks;
            Int src_snode_id = *(Int*)&src_blocks[0];
            Int src_nzblk_cnt = *(((Int*)&src_blocks[0])+1);
            NZBlockDesc * src_blocks_ptr = 
              reinterpret_cast<NZBlockDesc*>(&src_blocks[2*sizeof(Int)]);
            Int src_nzval_cnt = *(Int*)(src_blocks_ptr + src_nzblk_cnt);
            T * src_nzval_ptr = (T*)((Int*)(src_blocks_ptr + src_nzblk_cnt)+1);

            //Create the dummy supernode for that data
            SuperNode<T> dist_src_snode(src_snode_id,Xsuper_[src_snode_id-1],Xsuper_[src_snode_id]-1, src_blocks_ptr, src_nzblk_cnt, src_nzval_ptr, src_nzval_cnt);


            //              logfileptr->OFS()<<"IRECV Supernode "<<dist_src_snode.Id()<<std::endl;
#ifdef _DEBUG_
            logfileptr->OFS()<<"IRECV Supernode "<<dist_src_snode.Id()<<std::endl;
#endif

            //Update everything I own with that factor
            //Update the ancestors
            Int tgt_snode_id = 0;
            Int src_first_row = 0;
            Int src_last_row = 0;
            Int src_nzblk_idx = 0;

            TIMER_START(UPDATE_ANCESTORS);

            Int src_next_nzblk_idx = 0;
            while(FindNextUpdate(dist_src_snode, tgt_snode_id, src_first_row, src_nzblk_idx, src_last_row, src_next_nzblk_idx)){ 
              Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
              if(iTarget == iam){
#ifdef _DEBUG_
                logfileptr->OFS()<<"IRECV Supernode "<<tgt_snode_id<<" is updated by Supernode "<<dist_src_snode.Id()<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;
#endif

                Int iLocalJ = (tgt_snode_id-1) / np +1 ;
                SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];



#ifdef SINGLE_BLAS
                UpdateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, tmpBuf,src_colindx,src_rowindx,src_to_tgt_offset,src_first_row);
#else
                UpdateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, src_first_row);
#endif

                --UpdatesToDo(tgt_snode_id-1);



#ifdef _DEBUG_
                logfileptr->OFS()<<UpdatesToDo(tgt_snode_id-1)<<" updates left for Supernode "<<tgt_snode_id<<endl;
#endif

              }
            }
            TIMER_STOP(UPDATE_ANCESTORS);





            //delete the request from the list
//            assert(distance(cur_incomingRecv.begin(),it)<cur_incomingRecv.size());

            cur_incomingRecv.erase(it);
            --incomingRecvCnt_;

#ifdef _DEBUG_
            logfileptr->OFS()<<cur_incomingRecv.size()<<" async recv to do for Supernode "<<I<<endl;
            logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left for Supernode "<<I<<endl;
#endif
            if(UpdatesToDo(src_snode.Id()-1)==0){
              //cancel all requests
              incomingRecvCnt_-=cur_incomingRecv.size(); 
              cur_incomingRecv.clear();
            }


            it = WaitIncomingFactors(cur_incomingRecv,recv_status,outgoingSend);
          }


          src_blocks.resize(0);

          Int nz_cnt;
          Int max_bytes;
  //        std::vector<bool> received(Xsuper_.m(),false);
          while(UpdatesToDo(I-1)>0){
#ifdef _DEBUG_
            logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left"<<endl;
#endif


            TIMER_START(RECV_MALLOC);
            if(src_blocks.size()==0){
              max_bytes = 3*sizeof(Int); 
              //The upper bound must be of the width of the "largest" child
#ifdef _DEBUG_
              logfileptr->OFS()<<"Maximum width is "<<UpdateWidth_(I-1)<<std::endl;
#endif

              Int nrows = src_snode.NRowsBelowBlock(0);
              Int ncols = UpdateWidth_(I-1);
              nz_cnt = nrows * ncols;

              Int nblocks = nrows;//std::max((Int)ceil(nrows/2)+1,src_snode.NZBlockCnt());
              max_bytes += (nblocks)*sizeof(NZBlockDesc);
              max_bytes += nz_cnt*sizeof(T); 

              src_blocks.resize(max_bytes);
            }
            TIMER_STOP(RECV_MALLOC);

            TIMER_START(RECV_MPI);
            MPI_Status recv_status;
            int bytes_received = 0;

#ifdef PROBE_FIRST
            MPI_Probe(MPI_ANY_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
            MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);

            bool doabort = false;
            int prev_size = 0;
            if(src_blocks.size()<bytes_received){
              prev_size = src_blocks.size();
              doabort = true;
              //receive anyway
              src_blocks.resize(bytes_received);
            }
#endif


#ifdef PROBE_FIRST
            MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,recv_status.MPI_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
#else
            //receive the index array
            MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,MPI_ANY_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
#endif
            MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);

            Int src_snode_id = *(Int*)&src_blocks[0];
            Int src_nzblk_cnt = *(((Int*)&src_blocks[0])+1);
            NZBlockDesc * src_blocks_ptr = 
              reinterpret_cast<NZBlockDesc*>(&src_blocks[2*sizeof(Int)]);
            Int src_nzval_cnt = *(Int*)(src_blocks_ptr + src_nzblk_cnt);
            T * src_nzval_ptr = (T*)((Int*)(src_blocks_ptr + src_nzblk_cnt)+1);
            TIMER_STOP(RECV_MPI);
            //Create the dummy supernode for that data
            SuperNode<T> dist_src_snode(src_snode_id,Xsuper_[src_snode_id-1],Xsuper_[src_snode_id]-1, src_blocks_ptr, src_nzblk_cnt, src_nzval_ptr, src_nzval_cnt);


//            assert(!received[src_snode_id-1]);
//            received[src_snode_id-1] = true;

#ifdef PROBE_FIRST
            if(doabort){
              cout<<"We have a problem !!!! on P"<<iam<<"\n";
              gdb_lock();

              abort();
            }
#endif
#ifdef _DEBUG_
            logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode.Id()<<std::endl;
#endif

            //Update everything I own with that factor
            //Update the ancestors
            Int tgt_snode_id = 0;
            Int src_first_row = 0;
            Int src_last_row = 0;
            Int src_nzblk_idx = 0;


            TIMER_START(UPDATE_ANCESTORS);
            Int src_next_nzblk_idx = 0;
            while(FindNextUpdate(dist_src_snode, tgt_snode_id, src_first_row, src_nzblk_idx, src_last_row, src_next_nzblk_idx)){ 
              Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
              if(iTarget == iam){
#ifdef _DEBUG_
                logfileptr->OFS()<<"RECV Supernode "<<tgt_snode_id<<" is updated by Supernode "<<dist_src_snode.Id()<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;
#endif

                Int iLocalJ = (tgt_snode_id-1) / np +1 ;
                SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];



#ifdef SINGLE_BLAS
                UpdateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, tmpBuf,src_colindx,src_rowindx,src_to_tgt_offset,src_first_row);
#else
                UpdateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, src_first_row);
#endif

                --UpdatesToDo(tgt_snode_id-1);
#ifdef _DEBUG_
                logfileptr->OFS()<<UpdatesToDo(tgt_snode_id-1)<<" updates left for Supernode "<<tgt_snode_id<<endl;
#endif



              }
            }
            TIMER_STOP(UPDATE_ANCESTORS);
          }


 //                 timeEnd =  get_time( );

//          assert(UpdatesToDo(src_snode.Id()-1)==0);
#ifdef _DEBUG_
          if(UpdatesToDo(src_snode.Id()-1)!=0){gdb_lock();}
          assert(UpdatesToDo(src_snode.Id()-1)==0);
          logfileptr->OFS()<<"  Factoring Supernode "<<I<<" at "<<timeEnd-timeSta<<" s"<<std::endl;
#endif
          //if((I*100/(Xsuper_.m()-1)) % 20 == 0){
          //        logfileptr->OFS()<<"  Factoring Supernode "<<I<<" at "<<timeEnd-timeSta<<" s"<<std::endl;
          //}

          TIMER_START(FACTOR_PANEL);
          //Factorize Diagonal block
          NZBlockDesc & diag_desc = src_snode.GetNZBlockDesc(0);
          for(Int col = 0; col<src_snode.Size();col+=BLOCKSIZE){
            Int bw = min(BLOCKSIZE,src_snode.Size()-col);
            T * diag_nzval = &src_snode.GetNZval(diag_desc.Offset)[col+col*src_snode.Size()];
            lapack::Potrf( 'U', bw, diag_nzval, src_snode.Size());
            T * nzblk_nzval = &src_snode.GetNZval(diag_desc.Offset)[col+(col+bw)*src_snode.Size()];
            blas::Trsm('L','U','T','N',bw, src_snode.NRowsBelowBlock(0)-(col+bw), ONE<T>(),  diag_nzval, src_snode.Size(), nzblk_nzval, src_snode.Size());

            //update the rest !!! (next blocks columns)
            T * tgt_nzval = &src_snode.GetNZval(diag_desc.Offset)[col+bw+(col+bw)*src_snode.Size()];
            blas::Gemm('T','N',src_snode.Size()-(col+bw), src_snode.NRowsBelowBlock(0)-(col+bw),bw,MINUS_ONE<T>(),nzblk_nzval,src_snode.Size(),nzblk_nzval,src_snode.Size(),ONE<T>(),tgt_nzval,src_snode.Size());
          }


          //          T * diag_nzval = src_snode.GetNZval(diag_desc.Offset);
          //          lapack::Potrf( 'U', src_snode.Size(), diag_nzval, src_snode.Size());
          //        if(src_snode.NZBlockCnt()>1){
          //          NZBlockDesc & nzblk_desc = src_snode.GetNZBlockDesc(1);
          //          T * nzblk_nzval = src_snode.GetNZval(nzblk_desc.Offset);
          //          blas::Trsm('L','U','T','N',src_snode.Size(), src_snode.NRowsBelowBlock(1), ONE<T>(),  diag_nzval, src_snode.Size(), nzblk_nzval, src_snode.Size());
          //        }

#ifdef _DEBUG_
          //        logfileptr->OFS()<<src_snode<<std::endl;
#endif

          TIMER_STOP(FACTOR_PANEL);

          //if this is the last local supernode, send everything
          if(iLocalI == LocalSupernodes_.size()){

            for(Int iLocalJ=1;iLocalJ<=iLocalI;++iLocalJ){

              SuperNode<T> & src_snode = *LocalSupernodes_[iLocalJ -1];
              FOUpdate & curUpdate = nextUpdate[iLocalJ-1];


              //Send my factor to my ancestors. 

              Int tgt_snode_id = curUpdate.tgt_snode_id;
              Int src_nzblk_idx = curUpdate.src_nzblk_idx;
              Int src_first_row = curUpdate.src_first_row;
              Int src_last_row = curUpdate.src_last_row;
              Int src_next_nzblk_idx = curUpdate.src_next_nzblk_idx;

              TIMER_START(FIND_UPDATED_ANCESTORS);
              while(FindNextUpdate(src_snode, tgt_snode_id, src_first_row, src_nzblk_idx, src_last_row, src_next_nzblk_idx)){ 
                Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
                if(iTarget != iam){
                  #ifdef _DEBUG_
                  logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated by Supernode "<<src_snode.Id()<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;
                  if(curUpdate.is_factor_sent[iTarget]){
                    logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated Supernode "<<src_snode.Id()<<" in factor "<<curUpdate.is_factor_sent[iTarget]<<std::endl;
                  }
                  #endif

                  if(!curUpdate.is_factor_sent[iTarget]){

                    curUpdate.is_factor_sent[iTarget] = true;

                    Int tgt_first_col = Xsuper_(tgt_snode_id-1);
                    Int tgt_last_col = Xsuper_(tgt_snode_id)-1;

                    //Send
                    NZBlockDesc & pivot_desc = src_snode.GetNZBlockDesc(src_nzblk_idx);

                    Int local_first_row = src_first_row-pivot_desc.GIndex;
                    Int nzblk_cnt = src_snode.NZBlockCnt()-src_nzblk_idx;

                    Int nz_cnt = (src_snode.NRowsBelowBlock(src_nzblk_idx) - local_first_row )*src_snode.Size();
                    //              assert(nz_cnt>0);

                    T * nzval_ptr = src_snode.GetNZval(pivot_desc.Offset
                        +local_first_row*src_snode.Size());


                    TIMER_START(SEND_MALLOC);
                    AddOutgoingComm(outgoingSend, src_snode.Id(), src_snode.Size(), src_first_row, pivot_desc, nzblk_cnt, nzval_ptr, nz_cnt);
                    TIMER_STOP(SEND_MALLOC);

                    if(outgoingSend.size() > maxIsend_){
                      TIMER_START(SEND_MPI);
                      MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm());
                      TIMER_STOP(SEND_MPI);

#ifdef _DEBUG_            
                      logfileptr->OFS()<<"Sending "<<nz_cnt*sizeof(T)<<" bytes to P"<<iTarget<<std::endl;
                      logfileptr->OFS()<<"     Sent factor "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
                      logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
                      logfileptr->OFS()<<"Sending "<<outgoingSend.back()->size()<<" bytes to P"<<iTarget<<std::endl;
#endif
                      outgoingSend.pop_back();
                    }
                    else{
                      TIMER_START(SEND_MPI);
                      MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,CommEnv_->MPI_GetComm(),&outgoingSend.back()->Request);
                      TIMER_STOP(SEND_MPI);

#ifdef _DEBUG_            
                      logfileptr->OFS()<<"Sending "<<nz_cnt*sizeof(T)<<" bytes to P"<<iTarget<<std::endl;
                      logfileptr->OFS()<<"     Sent factor "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
                      logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
                      logfileptr->OFS()<<"Sending "<<outgoingSend.back()->size()<<" bytes to P"<<iTarget<<std::endl;
#endif
                    }
                  }
                }
                else{
                  //do we do the update now ?

                  Int iLocalTgt = (tgt_snode_id-1) / np +1 ;
                  SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalTgt -1];
#ifdef _DEBUG_
                  logfileptr->OFS()<<"LOCAL Supernode "<<tgt_snode.Id()<<" is updated by Supernode "<<src_snode.Id()<<" from row "<<src_first_row<<" "<<src_nzblk_idx<<std::endl;
#endif

#ifdef SINGLE_BLAS
                  UpdateSuperNode(src_snode,tgt_snode,src_nzblk_idx,tmpBuf,src_colindx,src_rowindx,src_to_tgt_offset, src_first_row);
#else
                  UpdateSuperNode(src_snode,tgt_snode,src_nzblk_idx, src_first_row);
#endif
                  --UpdatesToDo(tgt_snode_id-1);
#ifdef _DEBUG_
                  logfileptr->OFS()<<UpdatesToDo(tgt_snode_id-1)<<" updates left"<<endl;
#endif
                }


                curUpdate.tgt_snode_id = tgt_snode_id;
                curUpdate.src_nzblk_idx = src_nzblk_idx;
                curUpdate.src_first_row = src_first_row;
                curUpdate.src_last_row = src_last_row;
                curUpdate.src_next_nzblk_idx = src_next_nzblk_idx;





              }
              TIMER_STOP(FIND_UPDATED_ANCESTORS);
            } 



          }


          iLocalI++;
        }
        //      MPI_Barrier(CommEnv_->MPI_GetComm());


        //      {
        //      NumMat<T> tmp;
        //      GetFullFactors(tmp);
        //      }


      }
        I++;
    }


    for(Int idx = 0; idx <incomingRecvArr.size();++idx){
      assert(incomingRecvArr[idx].size()==0);
    } 

    MPI_Barrier(CommEnv_->MPI_GetComm());

    TIMER_STOP(FACTORIZATION_FO);
  }


#endif // _SUPERNODAL_MATRIX_IMPL_FO_HPP_
