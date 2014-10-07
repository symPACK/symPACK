#ifndef _SUPERNODAL_MATRIX_IMPL_FB_HPP_
#define _SUPERNODAL_MATRIX_IMPL_FB_HPP_

template<typename T> void SupernodalMatrix<T>::FBAsyncRecv(Int iLocalI, std::vector<AsyncComms> & incomingRecvAggArr, std::vector<AsyncComms * > & incomingRecvFactArr, IntNumVec & AggregatesToRecv, IntNumVec & FactorsToRecv){

  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();
  std::vector<SnodeUpdateFB> tmpTasks;
  Int nextLocalI = iLocalI;
  while(!LocalTasks.empty()){ 
    SnodeUpdateFB curTask = LocalTasks.top();

    Int IrecvCnt = 0; 
    if(curTask.type == FACTOR){
      AsyncComms & incomingRecvAgg = incomingRecvAggArr[nextLocalI-1]; 
      nextLocalI++;
      //this is a factorization task: we have to receive aggregates
      Int maxRecvCnt = AggregatesToRecv[curTask.tgt_snode_id-1];
      Int tag = AGG_TAG(curTask.src_snode_id,curTask.tgt_snode_id);

      for(Int idx =0; idx<maxRecvCnt && incomingRecvCnt_ + IrecvCnt < maxIrecv_;
          ++idx){
        Int max_bytes = getAggBufSize<T>(curTask, Xsuper_, UpdateHeight_);

        TIMER_START(ICOMM_MALLOC_AGG);
        incomingRecvAgg.push_back(new Icomm(max_bytes,MPI_REQUEST_NULL));
        TIMER_STOP(ICOMM_MALLOC_AGG);
        Icomm & Irecv = *incomingRecvAgg.back();
#ifdef _SEPARATE_COMM_
        MPI_Irecv(Irecv.front(),Irecv.capacity(),MPI_BYTE,MPI_ANY_SOURCE,tag,FBAggCommEnv_->MPI_GetComm(),&Irecv.Request);
#else
        MPI_Irecv(Irecv.front(),Irecv.capacity(),MPI_BYTE,MPI_ANY_SOURCE,tag,CommEnv_->MPI_GetComm(),&Irecv.Request);
#endif
        ++IrecvCnt;
      }

      AggregatesToRecv[curTask.tgt_snode_id-1] -= IrecvCnt;
    }
    else{

      //this is an update task: we potentially have to receive factors
      if(incomingRecvFactArr[curTask.tgt_snode_id-1] == NULL){
        incomingRecvFactArr[curTask.tgt_snode_id-1] = new AsyncComms();
      }
      AsyncComms * incomingRecvFact = incomingRecvFactArr[curTask.tgt_snode_id-1];

      Int maxRecvCnt = FactorsToRecv[curTask.tgt_snode_id-1];
      Int tag = FACT_TAG(curTask.src_snode_id,curTask.tgt_snode_id);

      for(Int idx =0; idx<maxRecvCnt && incomingRecvCnt_ + IrecvCnt < maxIrecv_;
          ++idx){
        Int max_bytes = getFactBufSize<T>(curTask, UpdateWidth_, UpdateHeight_);

        TIMER_START(ICOMM_MALLOC_FACT);
        incomingRecvFact->push_back(new Icomm(max_bytes,MPI_REQUEST_NULL));
        TIMER_STOP(ICOMM_MALLOC_FACT);
        Icomm & Irecv = *incomingRecvFact->back();
        MPI_Irecv(Irecv.front(),Irecv.capacity(),MPI_BYTE,MPI_ANY_SOURCE,tag,CommEnv_->MPI_GetComm(),&Irecv.Request);
        ++IrecvCnt;
      }
      FactorsToRecv[curTask.tgt_snode_id-1] -= IrecvCnt;
      assert(FactorsToRecv[curTask.tgt_snode_id-1]>=0);
    }

    incomingRecvCnt_+=IrecvCnt;

    if( incomingRecvCnt_ >= maxIrecv_){
      break;
    }

    LocalTasks.pop();
    tmpTasks.push_back(curTask);
  }

        TIMER_START(RESTORE_TASK_QUEUE);
  //put the tasks back into the task queue
  for(auto it = tmpTasks.begin(); it != tmpTasks.end(); it++){
    LocalTasks.push(*it);
  }
        TIMER_STOP(RESTORE_TASK_QUEUE);


}
















template <typename T> void SupernodalMatrix<T>::FBFactorizationTask(SnodeUpdateFB & curTask, Int iLocalI, IntNumVec & AggregatesDone, IntNumVec & AggregatesToRecv, std::vector<char> & src_blocks,std::vector<AsyncComms> & incomingRecvAggArr)
{

      Int src_snode_id = curTask.src_snode_id;
      Int tgt_snode_id = curTask.tgt_snode_id;
      Int I = src_snode_id;
      SuperNode<T> & src_snode = *LocalSupernodes_[iLocalI -1];
        Int src_first_col = src_snode.FirstCol();
        Int src_last_col = src_snode.LastCol();

#ifdef _DEBUG_
        logfileptr->OFS()<<"Processing Supernode "<<I<<std::endl;
#endif

#ifdef _DEBUG_PROGRESS_
        logfileptr->OFS()<<"Processing Supernode "<<I<<std::endl;
#endif


        //Receiving aggregates

        TIMER_START(RECV_AGGREGATES);
        src_blocks.resize(0);

        Int nz_cnt;
        
        Int max_bytes;
#ifdef _SEPARATE_COMM_
        max_bytes = getMaxBufSize<T>( UpdateWidth_, UpdateHeight_);
#endif




      //first wait for the Irecv
      AsyncComms & cur_incomingRecv = incomingRecvAggArr[iLocalI-1];
      MPI_Status recv_status;

        TIMER_START(IRECV_MPI_AGG);
      AsyncComms::iterator it = WaitIncomingFactors(cur_incomingRecv, recv_status,outgoingSend);
      while( it != cur_incomingRecv.end() ){
        Icomm * curComm = *it;

          SuperNode<T> dist_src_snode;
          size_t read_bytes = Deserialize(curComm->front(),dist_src_snode);
          //Deserialize the number of aggregates
          //Int * aggregatesCnt = (Int *)(curComm->front()+read_bytes);


#ifdef _DEBUG_
          logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode.Id()<<std::endl;
#endif

          src_snode.Aggregate(dist_src_snode);

        //delete the request from the list
        cur_incomingRecv.erase(it);
        --incomingRecvCnt_;

        it = WaitIncomingFactors(cur_incomingRecv,recv_status,outgoingSend);
      }
        TIMER_STOP(IRECV_MPI_AGG);







        //Wait for all the aggregates BUT receive from any
        while(AggregatesToRecv(tgt_snode_id-1)>0){
#ifdef _DEBUG_
          logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left"<<endl;
#endif

          TIMER_START(RECV_MPI);
          MPI_Status recv_status;

#ifdef _SEPARATE_COMM_
          TIMER_START(RECV_MALLOC);
          src_blocks.resize(max_bytes);
          TIMER_STOP(RECV_MALLOC);

          //Do an IProbe on the TAG we want to receive 
          Int tag = AGG_TAG(curTask.src_snode_id,curTask.tgt_snode_id);
          Int flag = 0;
          MPI_Iprobe(MPI_ANY_SOURCE,tag,FBAggCommEnv_->MPI_GetComm(), &flag,&recv_status);
          if(flag){
            MPI_Recv(&src_blocks[0],src_blocks.size(),MPI_BYTE,recv_status.MPI_SOURCE,tag,FBAggCommEnv_->MPI_GetComm(),&recv_status);
          }
          //if nothing available, receive from any tag
          else{
            MPI_Recv(&src_blocks[0],src_blocks.size(),MPI_BYTE,MPI_ANY_SOURCE,MPI_ANY_TAG,FBAggCommEnv_->MPI_GetComm(),&recv_status);
            tag = recv_status.MPI_TAG;
          }

          Int tgt_id = AGG_TAG_TO_ID(tag);

          Int iLocalJ = globToLocSnodes_.IntervalSearch(tgt_id,tgt_id)->block_idx;
          SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ-1];

          SuperNode<T> dist_src_snode;
          Deserialize(&src_blocks[0],dist_src_snode);

#ifdef _DEBUG_
          logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode.Id()<<" for Supernode "<<tgt_snode.Id()<<std::endl;
#endif


          tgt_snode.Aggregate(dist_src_snode);
          AggregatesToRecv[tgt_snode.Id()-1] -=1;


#else
          Int max_bytes = getAggBufSize<T>(curTask, Xsuper_, UpdateHeight_);

          TIMER_START(RECV_MALLOC);
          src_blocks.resize(max_bytes);
          TIMER_STOP(RECV_MALLOC);

          Int tag = AGG_TAG(curTask.src_snode_id,curTask.tgt_snode_id);
          MPI_Recv(&src_blocks[0],src_blocks.size(),MPI_BYTE,MPI_ANY_SOURCE,tag,CommEnv_->MPI_GetComm(),&recv_status);
          TIMER_STOP(RECV_MPI);

          SuperNode<T> dist_src_snode;
          size_t read_bytes = Deserialize(&src_blocks[0],dist_src_snode);
          //Deserialize the number of aggregates
          //Int * aggregatesCnt = (Int *)(&src_blocks[0]+read_bytes);

#ifdef _DEBUG_
          logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode.Id()<<std::endl;
#endif

          src_snode.Aggregate(dist_src_snode);
          AggregatesToRecv[src_snode.Id()-1] -=1;
#endif

        }
        //clear the buffer
        //{ vector<char>().swap(src_blocks);  }

        TIMER_STOP(RECV_AGGREGATES);

          assert(AggregatesToRecv[src_snode.Id()-1]==0);

#ifdef TRACK_PROGRESS
  Real timeStaTask =  get_time( );
#endif

#ifdef _DEBUG_PROGRESS_
        logfileptr->OFS()<<"Factoring Supernode "<<I<<std::endl;
#endif

        TIMER_START(FACTOR_PANEL);
        src_snode.Factorize();
        TIMER_STOP(FACTOR_PANEL);

#ifdef TRACK_PROGRESS
  timeEnd =  get_time( );
  *progstr<<src_snode.Id()<<" "<<timeStaTask - timeSta<<" "<<timeEnd - timeSta<<" F"<<src_snode.Id()<<" "<<iam<<endl;
#endif

        //Sending factors

        //Send my factor to my ancestors. 
        BolNumVec is_factor_sent(np);
        SetValue(is_factor_sent,false);
        BolNumVec is_skipped(np);
        SetValue(is_skipped,false);

        SnodeUpdate curUpdate;
        TIMER_START(FIND_UPDATED_ANCESTORS);
        while(src_snode.FindNextUpdate(curUpdate,Xsuper_,SupMembership_)){ 
          Int iTarget = Mapping_->Map(curUpdate.tgt_snode_id-1,src_snode.Id()-1);

          if(iTarget != iam){
            if(!is_factor_sent[iTarget]){


                Int tag = FACT_TAG(curUpdate.src_snode_id,curUpdate.tgt_snode_id);
                FBDelayedComm comm(FACTOR,(void*)&src_snode,curUpdate.src_snode_id,curUpdate.tgt_snode_id,curUpdate.blkidx,curUpdate.src_first_row,iTarget,tag);
                //Try to do the async send first
                if(outgoingSend.size() < maxIsend_){
                  SendMessage(comm, outgoingSend);
                }
                else{
TIMER_START(PUSH_MSG);
                  MsgToSend.push(comm);
TIMER_STOP(PUSH_MSG);
                }
#ifdef _DEBUG_DELAY_
                logfileptr->OFS()<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<curUpdate.tgt_snode_id<<" from row "<<curUpdate.src_first_row<<endl;
                cout<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<curUpdate.tgt_snode_id<<" from row "<<curUpdate.src_first_row<<endl;
#endif  
                is_factor_sent[iTarget] = true;

            }

          }
        }
        TIMER_STOP(FIND_UPDATED_ANCESTORS);

       //process some of the delayed send
//       AdvanceOutgoing(outgoingSend);
//       SendDelayedMessagesUp(MsgToSend,outgoingSend,LocalTasks);




}



template <typename T> void SupernodalMatrix<T>::FBUpdateTask(SnodeUpdateFB & curTask, IntNumVec & UpdatesToDo, IntNumVec & AggregatesDone,std::vector< SuperNode<T> * > & aggVectors, std::vector<char> & src_blocks, std::vector<AsyncComms * > & incomingRecvFactArr, IntNumVec & FactorsToRecv)
{
  Int src_snode_id = curTask.src_snode_id;
  Int tgt_snode_id = curTask.tgt_snode_id;

  src_snode_id = abs(src_snode_id);
  bool is_first_local = curTask.src_snode_id <0;
  curTask.src_snode_id = src_snode_id;

  SuperNode<T> * cur_src_snode; 

  src_blocks.resize(0);
  AsyncComms::iterator it;
  AsyncComms * cur_incomingRecv = incomingRecvFactArr[tgt_snode_id-1];
  cur_src_snode = FBRecvFactor(curTask,src_blocks,cur_incomingRecv,it,FactorsToRecv);

  //need to be updated because we might have received from someone else
  src_snode_id = cur_src_snode->Id();
#ifdef _DEBUG_PROGRESS_
  if(abs(curTask.src_snode_id) != src_snode_id){ 
    cout<<"YOUHOU WE HAVE ASYNC HERE !!! expected: "<< abs(curTask.src_snode_id) << " vs received: "<<src_snode_id<<endl;
    logfileptr->OFS()<<"YOUHOU WE HAVE ASYNC HERE !!! expected: "<< abs(curTask.src_snode_id) << " vs received: "<<src_snode_id<<endl;
  }
#endif
  Int iSrcOwner = Mapping_->Map(src_snode_id-1,src_snode_id-1);


  //Update everything src_snode_id own with that factor
  //Update the ancestors
  SnodeUpdate curUpdate;
  if(is_first_local){

    curUpdate.tgt_snode_id  = curTask.tgt_snode_id;
    curUpdate.src_first_row = Xsuper_[curTask.tgt_snode_id-1];
    while( (curUpdate.next_blkidx = cur_src_snode->FindBlockIdx(curUpdate.src_first_row)) == -1){curUpdate.src_first_row++;} 
    curUpdate.src_next_row  = curUpdate.src_first_row;
    curUpdate.blkidx        = curUpdate.next_blkidx;
  }

  TIMER_START(UPDATE_ANCESTORS);
  while(cur_src_snode->FindNextUpdate(curUpdate,Xsuper_,SupMembership_,iam==iSrcOwner)){
    Int iUpdater = Mapping_->Map(curUpdate.tgt_snode_id-1,cur_src_snode->Id()-1);
    if(iUpdater == iam){

#ifdef _DEBUG_PROGRESS_
      logfileptr->OFS()<<"implicit Task: {"<<curUpdate.src_snode_id<<" -> "<<curUpdate.tgt_snode_id<<"}"<<std::endl;
      logfileptr->OFS()<<"Processing update from Supernode "<<curUpdate.src_snode_id<<" to Supernode "<<curUpdate.tgt_snode_id<<endl;
#endif
      SuperNode<T> * tgt_aggreg;

      Int iTarget = Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
      if(iTarget == iam){
        //the aggregate vector is directly the target snode
        //Int iLocalJ = (curUpdate.tgt_snode_id-1) / np +1 ;
        Int iLocalJ = globToLocSnodes_.IntervalSearch(curUpdate.tgt_snode_id,curUpdate.tgt_snode_id)->block_idx;
        tgt_aggreg = LocalSupernodes_[iLocalJ -1];

        assert(curUpdate.tgt_snode_id == tgt_aggreg->Id());

      }
      else{
        //Check if src_snode_id already have an aggregate vector
        //                if(aggVectors[curUpdate.tgt_snode_id-1]==NULL)
        if(AggregatesDone[curUpdate.tgt_snode_id-1]==0){
#ifdef COMPACT_AGGREGATES
          aggVectors[curUpdate.tgt_snode_id-1] = new SuperNode<T>(curUpdate.tgt_snode_id, Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id]-1);
#else
          aggVectors[curUpdate.tgt_snode_id-1] = new SuperNode<T>(curUpdate.tgt_snode_id, Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id]-1,  xlindx_, lindx_);
#endif
        }
        tgt_aggreg = aggVectors[curUpdate.tgt_snode_id-1];



      }

#ifdef _DEBUG_
      logfileptr->OFS()<<"RECV Supernode "<<curUpdate.tgt_snode_id<<" is updated by Supernode "<<cur_src_snode->Id()<<" rows "<<curUpdate.src_first_row<<" "<<curUpdate.blkidx<<std::endl;
#endif


      //Update the aggregate

#ifdef TRACK_PROGRESS
      Real timeStaTask =  get_time( );
#endif

      tgt_aggreg->UpdateAggregate(*cur_src_snode,curUpdate,tmpBufs,iTarget);

#ifdef TRACK_PROGRESS
      timeEnd =  get_time( );
      *progstr<<curUpdate.tgt_snode_id<<" "<<timeStaTask - timeSta<<" "<<timeEnd - timeSta<<" U"<<curUpdate.src_snode_id<<"-"<<curUpdate.tgt_snode_id<<" "<<iam<< endl;
#endif

      --UpdatesToDo[curUpdate.tgt_snode_id-1];
      ++AggregatesDone[curUpdate.tgt_snode_id-1];
#ifdef _DEBUG_
      logfileptr->OFS()<<UpdatesToDo(curUpdate.tgt_snode_id-1)<<" updates left for Supernode "<<curUpdate.tgt_snode_id<<endl;
#endif

      //Send the aggregate if it's the last
      //If this is my last update sent it to curUpdate.tgt_snode_id
      if(UpdatesToDo[curUpdate.tgt_snode_id-1]==0){
        if(iTarget != iam){
#ifdef _DEBUG_
          logfileptr->OFS()<<"Remote Supernode "<<curUpdate.tgt_snode_id<<" is updated by Supernode "<<src_snode_id<<std::endl;
#endif



          //if(curUpdate.tgt_snode_id==14){logfileptr->OFS()<<*tgt_aggreg<<endl;}

          //Push the comm in the comm list
          Int tag = AGG_TAG(curUpdate.src_snode_id,curUpdate.tgt_snode_id);
          NZBlockDesc & pivot_desc = tgt_aggreg->GetNZBlockDesc(0);

          FBDelayedComm comm(AGGREGATE,(void*)tgt_aggreg,curUpdate.src_snode_id,curUpdate.tgt_snode_id,0,pivot_desc.GIndex,iTarget,tag,AggregatesDone[curUpdate.tgt_snode_id-1]);

            if(outgoingSend.size() < maxIsend_){
              SendMessage(comm, outgoingSend);
            }
            else{
TIMER_START(PUSH_MSG);
              MsgToSend.push(comm);
TIMER_STOP(PUSH_MSG);
            }




        }
      }

      //if local update, push a new task in the queue and stop the while loop
      if(iam==iSrcOwner ){
        break;
      }
    }
  }
  TIMER_STOP(UPDATE_ANCESTORS);

  //process some of the delayed send
//  AdvanceOutgoing(outgoingSend);
//  SendDelayedMessagesUp(MsgToSend,outgoingSend,LocalTasks);



  if(iSrcOwner!=iam){
    //if the recv was from an async recv, delete the request
    if(cur_incomingRecv != NULL){
      if(it != cur_incomingRecv->end()){
        //delete the request from the list
        cur_incomingRecv->erase(it);
        --incomingRecvCnt_;
      }
      if(cur_incomingRecv->empty() && FactorsToRecv[tgt_snode_id-1]==0){
        delete incomingRecvFactArr[tgt_snode_id-1];
        incomingRecvFactArr[tgt_snode_id-1] = NULL;
      }
    }
    delete cur_src_snode;
    //clear the buffer
    //{ vector<char>().swap(src_blocks);  }
  }

}



template <typename T> void SupernodalMatrix<T>::FanBoth()
{
  TIMER_START(FACTORIZATION_FB);

  Real timeSta, timeEnd;


  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();

  IntNumVec UpdatesToDo;
  IntNumVec AggregatesToRecv;

  std::vector<AsyncComms> incomingRecvAggArr(LocalSupernodes_.size());
  std::vector<AsyncComms * > incomingRecvFactArr(Xsuper_.m(),NULL);
  incomingRecvCnt_ = 0;
  
  double timesta2 = get_time(); 
  FBGetUpdateCount(UpdatesToDo,AggregatesToRecv);
#ifdef COMPACT_AGGREGATES
  xlindx_.Clear();
  lindx_.Clear();
#endif


#ifdef _SEPARATE_COMM_
  FBAggCommEnv_ = new CommEnvironment(*CommEnv_);
#endif

  double timeend2 = get_time(); 
  logfileptr->OFS()<<"Update count time: "<<timeend2-timesta2<<endl;


  IntNumVec AggregatesDone(Xsuper_.m());
  SetValue(AggregatesDone,0);


  //Array for Irecv of factors
  std::vector<AsyncComms> incomingRecvArr(LocalSupernodes_.size());
  incomingRecvCnt_ = 0;

  std::vector<T> src_nzval;
  std::vector<char> src_blocks;


  Int maxwidth = 0;
  for(Int i = 1; i<Xsuper_.m(); ++i){
    Int width =Xsuper_(i) - Xsuper_(i-1);
    if(width>=maxwidth){
      maxwidth = width;
    }
  }
  tmpBufs.Resize(Size(),maxwidth);
  std::vector< SuperNode<T> * > aggVectors(Xsuper_.m()-1,NULL);


  timeSta =  get_time( );

  TIMER_START(BUILD_TASK_LIST);
  //IntNumVec FactorsToRecv = UpdatesToDo;
  IntNumVec FactorsToRecv(Xsuper_.m());
  SetValue(FactorsToRecv,I_ZERO);
  for(Int I = 1; I<Xsuper_.m(); ++I){
    Int iOwner = Mapping_->Map(I-1,I-1);
    if(iam==iOwner){
      SnodeUpdateFB curUpdate;
      curUpdate.type=FACTOR;
      curUpdate.src_snode_id = I;
      curUpdate.tgt_snode_id = I;
      LocalTasks.push(curUpdate);
      Int J = -1;
     
      bool is_first = true; 
      while( (J=FBUpdate(I,J))!= -1 ){
        SnodeUpdateFB curUpdate;
        curUpdate.type=AGGREGATE;
        //set it to be negative as it is a local update
//        FactorsToRecv[J-1]--;

        curUpdate.src_snode_id = is_first?I:-I;
        curUpdate.tgt_snode_id = J;
        LocalTasks.push(curUpdate);
        is_first=false;
      }
    }
    else{
      Int J = FBUpdate(I); 
      if(J!=-1){
        SnodeUpdateFB curUpdate;
        curUpdate.type=AGGREGATE;
        curUpdate.src_snode_id = I;
        curUpdate.tgt_snode_id = J;
        LocalTasks.push(curUpdate);
        FactorsToRecv[J-1]++;
      }
    }
  }
  TIMER_STOP(BUILD_TASK_LIST);

#ifdef _DEBUG_UPDATES_
  logfileptr->OFS()<<"UpdatesToDo: "<<UpdatesToDo<<endl;
  logfileptr->OFS()<<"AggregatesToRecv: "<<AggregatesToRecv<<endl;
  logfileptr->OFS()<<"FactorsToRecv: "<<FactorsToRecv<<endl;
#endif

  Int iLocalI = 1;
  while(!LocalTasks.empty() || !MsgToSend.empty() || !outgoingSend.empty()){

    //Check for completion of outgoing communication
    AdvanceOutgoing(outgoingSend);
    //We have to have a list of remote supernodes for which we are computing updates

    if(!LocalTasks.empty()){
      //make a copy because we need to pop it
#ifdef _DEADLOCK_
      SnodeUpdateFB curTask = LocalTasks.front();
#else
      SnodeUpdateFB curTask = LocalTasks.top();
#endif


#ifdef _DEBUG_DELAY_
  {
    FBTasks tmp = LocalTasks;
    logfileptr->OFS()<<"Task Queue : ";
    while( tmp.size()>0){
      //Pull the highest priority message
#ifdef _DEADLOCK_
      const SnodeUpdateFB & comm = tmp.front();
#else
      const SnodeUpdateFB & comm = tmp.top();
#endif
      Int src_snode_id = comm.src_snode_id;
      Int tgt_id = comm.tgt_snode_id;
      logfileptr->OFS()<<" { "<<src_snode_id<<" -> "<<tgt_id<<" }";
      tmp.pop();
    }
    logfileptr->OFS()<<endl;
  }
#endif



#ifdef _DEBUG_PROGRESS_
        logfileptr->OFS()<<"picked Task: {"<<curTask.src_snode_id<<" -> "<<curTask.tgt_snode_id<<"}"<<std::endl;
#endif

      //Launch Irecv for subsequent local supernodes if I can
      FBAsyncRecv(iLocalI,incomingRecvAggArr,incomingRecvFactArr,AggregatesToRecv, FactorsToRecv);

      //If it is a factorization
      if(curTask.type == FACTOR){
        FBFactorizationTask(curTask,iLocalI,AggregatesDone,AggregatesToRecv, src_blocks,incomingRecvAggArr );
        ++iLocalI; 
      }
      else{
        FBUpdateTask(curTask, UpdatesToDo, AggregatesDone, aggVectors, src_blocks,incomingRecvFactArr,FactorsToRecv);
      }
      LocalTasks.pop();
    }

    //process some of the delayed send
    AdvanceOutgoing(outgoingSend);
    SendDelayedMessagesUp(MsgToSend,outgoingSend,LocalTasks);

  }


  for(Int idx = 0; idx <incomingRecvAggArr.size();++idx){
    assert(incomingRecvAggArr[idx].size()==0);
  } 

  for(Int idx = 0; idx <incomingRecvFactArr.size();++idx){
    if(incomingRecvFactArr[idx]!=NULL){
      if(incomingRecvFactArr[idx]->size()!=0){gdb_lock();}
      assert(incomingRecvFactArr[idx]->size()==0);
    }
  } 

  MPI_Barrier(CommEnv_->MPI_GetComm());

  tmpBufs.Clear();

#ifdef _SEPARATE_COMM_
  delete FBAggCommEnv_;
#endif

  TIMER_STOP(FACTORIZATION_FB);
}


template <typename T> void SupernodalMatrix<T>::FBGetUpdateCount(IntNumVec & sc, IntNumVec & atr){
  sc.Resize(Xsuper_.m());
  SetValue(sc,I_ZERO);

  atr.Resize(Xsuper_.m());
  SetValue(atr,I_ZERO);


  IntNumVec marker(Xsuper_.m());
  SetValue(marker,I_ZERO);

  std::vector<bool>isSent(Xsuper_.m()*np,false);

  for(Int s = 1; s<Xsuper_.m(); ++s){
    Int first_col = Xsuper_(s-1);
    Int last_col = Xsuper_(s)-1;

    Int fi = xlindx_(s-1);
    Int li = xlindx_(s)-1;


#ifdef _DEBUG_UPDATES_
    logfileptr->OFS()<<"Supernode "<<s<<" updates: ";
#endif


    for(Int row_idx = fi; row_idx<=li;++row_idx){
      Int row = lindx_(row_idx-1);
      Int supno = SupMembership_(row-1);

      if(marker(supno-1)!=s && supno!=s){


        Int iFactorizer = Mapping_->Map(supno-1,supno-1);
        Int iUpdater = Mapping_->Map(supno-1,s-1);

    if(!isSent[(supno-1)*np+iUpdater] && iUpdater!=iFactorizer){
      atr[supno-1]++;
      isSent[(supno-1)*np+iUpdater]=true;
    }


        if(iam == iUpdater){

#ifdef _DEBUG_UPDATES_
          logfileptr->OFS()<<supno<<" ";
#endif
          ++sc[supno-1];
//          lu[supno-1] = s;
        }
        else{

#ifdef _DEBUG_UPDATES_
          logfileptr->OFS()<<supno<<" [P"<<iUpdater<<"] ";
#endif
        }

        marker(supno-1) = s;
      }
    }

#ifdef _DEBUG_UPDATES_
    logfileptr->OFS()<<std::endl;
#endif

  }
}


template<typename T> SuperNode<T> * SupernodalMatrix<T>::FBRecvFactor(const SnodeUpdateFB & curTask, std::vector<char> & src_blocks,AsyncComms * cur_incomingRecv,AsyncComms::iterator & it, IntNumVec & FactorsToRecv){
  TIMER_START(RECV_FACTORS);
  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();
  SuperNode<T> * cur_src_snode = NULL;
  Int iSrcOwner = Mapping_->Map(curTask.src_snode_id-1,curTask.src_snode_id-1);
  

  //This is a local update, iSrcOwner matters
  if(iSrcOwner==iam){
    //Int iLocalI = (curTask.src_snode_id-1) / np +1 ;
    Int iLocalI = globToLocSnodes_.IntervalSearch(curTask.src_snode_id,curTask.src_snode_id)->block_idx;
    cur_src_snode = LocalSupernodes_[iLocalI -1];
  }
  else{

    Int tag = FACT_TAG(curTask.src_snode_id,curTask.tgt_snode_id);
#if 0
    TIMER_START(RECV_MPI);
    MPI_Status recv_status;
    int bytes_received = 0;

    //This is a remote update, iSrcOwner doesn't matter
    MPI_Probe(MPI_ANY_SOURCE,tag,CommEnv_->MPI_GetComm(),&recv_status);
    MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);

    TIMER_START(RECV_MALLOC);
    src_blocks.resize(bytes_received);
    TIMER_STOP(RECV_MALLOC);

//    logfileptr->OFS()<<"RECV Bfore"<<endl;
    MPI_Recv(&src_blocks[0],bytes_received,MPI_BYTE,recv_status.MPI_SOURCE,recv_status.MPI_TAG,CommEnv_->MPI_GetComm(),&recv_status);
    MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
//    logfileptr->OFS()<<"RECV After"<<endl;


    TIMER_STOP(RECV_MPI);
#else
    TIMER_START(RECV_MPI);
      SuperNode<T> * dist_src_snode = new SuperNode<T>();
      //Receive the factor
      //first wait for the Irecv
      bool do_blocking_recv = true;
      if(cur_incomingRecv != NULL){
        TIMER_START(IRECV_MPI_FACT);
        MPI_Status recv_status;
        it = WaitIncomingFactors(*cur_incomingRecv, recv_status,outgoingSend);
        if( it != cur_incomingRecv->end() ){
          Icomm * curComm = *it;
          Deserialize(curComm->front(),*dist_src_snode);

          do_blocking_recv = false;
        }
        TIMER_STOP(IRECV_MPI_FACT);
      }

      //do a blocking recv otherwise
      if( do_blocking_recv){

if(FactorsToRecv[curTask.tgt_snode_id-1]<=0){gdb_lock();}
assert(FactorsToRecv[curTask.tgt_snode_id-1]>0);
        MPI_Status recv_status;
        int bytes_received = 0;

        Int max_bytes = getFactBufSize<T>(curTask, UpdateWidth_, UpdateHeight_);

        TIMER_START(RECV_MALLOC);
        src_blocks.resize(max_bytes);
        TIMER_STOP(RECV_MALLOC);

        MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,MPI_ANY_SOURCE,tag,CommEnv_->MPI_GetComm(),&recv_status);
        //    logfileptr->OFS()<<"RECV After"<<endl;
#endif
        FactorsToRecv[curTask.tgt_snode_id-1]--;

        Deserialize(&src_blocks[0],*dist_src_snode);
        //#if not(defined(_LINEAR_SEARCH_FCLC_) || defined(_LAZY_INIT) )
      }

    dist_src_snode->InitIdxToBlk();
//#endif

    TIMER_STOP(RECV_MPI);




#ifdef _DEBUG_PROGRESS_
    logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode->Id()<<std::endl;
#endif

#ifdef _DEBUG_
    logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode->Id()<<std::endl;
#endif
    cur_src_snode = dist_src_snode;

  }
  TIMER_STOP(RECV_FACTORS);
  return cur_src_snode;
}

template<typename T> Int SupernodalMatrix<T>::FBUpdate(Int I,Int prevJ){
  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();
  //Check if I have anything to update with that supernode
  //look at lindx_
  Int fi = xlindx_(I-1);
  Int li = xlindx_(I)-1;

  Int iOwner = Mapping_->Map(I-1,I-1);


  Int firstUpdate = -1; 
  Int J = -1; 
  for(Int idx = fi; idx<=li;++idx){
    Int row = lindx_[idx-1];
    J = SupMembership_[row-1];
    Int iUpdater = Mapping_->Map(J-1,I-1);


    if(iUpdater == iam && J>I){
      if(iUpdater==iOwner){
        if(J>prevJ){
          firstUpdate = J;
          break;
        }
      }
      else{
        if(J>prevJ){
          firstUpdate = J;
          break;
        }
      }
    }
  }

  return firstUpdate;
}






#endif
