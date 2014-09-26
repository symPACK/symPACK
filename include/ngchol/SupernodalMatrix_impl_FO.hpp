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
#ifdef _DEADLOCK_
        const DelayedComm & comm = tmp.front();
#else
        const DelayedComm & comm = tmp.top();
#endif
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
#ifdef _DEADLOCK_
      const DelayedComm & comm = MsgToSend.front();
#else
      const DelayedComm & comm = MsgToSend.top();
#endif
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

        Int iTarget = Mapping_->Map(tgt_snode_id-1,tgt_snode_id-1);
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
          int error_code = MPI_Test(&(curComm->Request),&done,&recv_status);

          //Test if comm is done
          if(done==1){

            Int bytes_received = 0;
            MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
            curComm->setHead(bytes_received);

            TIMER_STOP(IRECV_MPI);
            return it;
          }
        }
      }
    }
  }




template <typename T> void SupernodalMatrix<T>::FanOut( ){

  TIMER_START(FACTORIZATION_FO);

  xlindx_.Clear();
  lindx_.Clear();

  Real timeSta, timeEnd;
  timeSta =  get_time( );


  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();

  IntNumVec UpdatesToDo = UpdateCount_;





  CommList FactorsToSend; 
  //FBCommList FactorsToSend; 

  std::vector<AsyncComms> incomingRecvArr(LocalSupernodes_.size());
  incomingRecvCnt_ = 0;
  IntNumVec FactorsToRecv(LocalSupernodes_.size());

  std::vector<char> src_blocks;

  std::vector<std::queue<SnodeUpdate> > LocalUpdates(LocalSupernodes_.size());

  Int maxwidth = 0;
  for(Int i = 1; i<Xsuper_.m(); ++i){
    Int width =Xsuper_(i) - Xsuper_(i-1);
    if(width>=maxwidth){
      maxwidth = width;
    }
  }
  
  tmpBufs.Resize(Size(),maxwidth);



  for(Int i = LocalSupernodes_.size()-1; i>=0; --i){
    SuperNode<T> & cur_snode = *LocalSupernodes_[i];

    SnodeUpdate curUpdate;
    while(cur_snode.FindNextUpdate(curUpdate,Xsuper_,SupMembership_)){ 
      Int iTarget = Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
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

#ifdef _DEBUG_PROGRESS_
      logfileptr->OFS()<<"Processing Supernode "<<I<<std::endl;
#endif

      //Launch Irecv for subsequent local supernodes if I can
      AsyncRecvFactors(iLocalI,incomingRecvArr,FactorsToRecv,UpdatesToDo);

      //Do all my updates (Local and remote)
      //Local updates
      while(!LocalUpdates[iLocalI-1].empty()){
        SnodeUpdate & curUpdate = LocalUpdates[iLocalI-1].front();
        SuperNode<T> & local_src_snode = *LocalSupernodes_[(curUpdate.src_snode_id-1) / np];


#ifdef _DEBUG_PROGRESS_
logfileptr->OFS()<<"Processing update from Supernode "<<curUpdate.src_snode_id<<" to Supernode "<<curUpdate.tgt_snode_id<<endl;
#endif

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
#ifndef _LINEAR_SEARCH_FCLC_
        dist_src_snode.InitIdxToBlk();
#endif

#ifdef _DEBUG_
        logfileptr->OFS()<<"IRECV Supernode "<<dist_src_snode.Id()<<std::endl;
#endif

        //Update everything I own with that factor
        //Update the ancestors

        TIMER_START(UPDATE_ANCESTORS);
        SnodeUpdate curUpdate;
        while(dist_src_snode.FindNextUpdate(curUpdate,Xsuper_,SupMembership_,false)){ 
          Int iTarget = Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
          if(iTarget == iam){
            Int iLocalJ = (curUpdate.tgt_snode_id-1) / np +1 ;
            SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];

#ifdef _DEBUG_PROGRESS_
logfileptr->OFS()<<"Processing update from Supernode "<<curUpdate.src_snode_id<<" to Supernode "<<curUpdate.tgt_snode_id<<endl;
#endif

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
        TIMER_STOP(RECV_MPI);

        SuperNode<T> dist_src_snode;
        Deserialize(&src_blocks[0],dist_src_snode);
#ifndef _LINEAR_SEARCH_FCLC_
        dist_src_snode.InitIdxToBlk();
#endif

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
          Int iTarget = Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
          if(iTarget == iam){
            Int iLocalJ = (curUpdate.tgt_snode_id-1) / np +1 ;
            SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];
#ifdef _DEBUG_PROGRESS_
logfileptr->OFS()<<"Processing update from Supernode "<<curUpdate.src_snode_id<<" to Supernode "<<curUpdate.tgt_snode_id<<endl;
#endif
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


      timeEnd =  get_time( );
#ifdef _DEBUG_
      if(UpdatesToDo(src_snode.Id()-1)!=0){gdb_lock();}
      assert(UpdatesToDo(src_snode.Id()-1)==0);
      logfileptr->OFS()<<"  Factoring Supernode "<<I<<" at "<<timeEnd-timeSta<<" s"<<std::endl;
#endif
#ifdef _DEBUG_PROGRESS_
      logfileptr->OFS()<<"  Factoring Supernode "<<I<<std::endl;
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
        Int iTarget = Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);

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
#ifndef _DEADLOCK_
//            if( next_local_snode< curUpdate.tgt_snode_id){
              //need to push the prev src_last_row
              FactorsToSend.push(DelayedComm(src_snode.Id(),curUpdate.tgt_snode_id,curUpdate.blkidx,curUpdate.src_first_row));
//              Int tag =curUpdate.tgt_snode_id;
//              FactorsToSend.push(FBDelayedComm(FACTOR,(void*)&src_snode,curUpdate.src_snode_id,curUpdate.tgt_snode_id,curUpdate.blkidx,curUpdate.src_first_row,iTarget,tag));
#ifdef _DEBUG_DELAY_
              logfileptr->OFS()<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<curUpdate.tgt_snode_id<<" from row "<<curUpdate.src_first_row<<endl;
              cout<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<curUpdate.tgt_snode_id<<" from row "<<curUpdate.src_first_row<<endl;
#endif
              is_skipped[iTarget] = true;


              SendDelayedMessagesUp(iLocalI,FactorsToSend,outgoingSend,LocalSupernodes_);


//            }
#else
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
#endif

          }
        }
      }
      TIMER_STOP(FIND_UPDATED_ANCESTORS);
    }


    //process some of the delayed send
    SendDelayedMessagesUp(iLocalI,FactorsToSend,outgoingSend,LocalSupernodes_);
//    SendDelayedMessagesUp(FactorsToSend, outgoingSend, FBTasks & taskList, NULL);
    iLocalI++;
  }


  for(Int idx = 0; idx <incomingRecvArr.size();++idx){
    assert(incomingRecvArr[idx].size()==0);
  } 

  MPI_Barrier(CommEnv_->MPI_GetComm());

  tmpBufs.Clear();

  TIMER_STOP(FACTORIZATION_FO);
}



#endif // _SUPERNODAL_MATRIX_IMPL_FO_HPP_
