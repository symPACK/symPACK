#ifndef _SUPERNODAL_MATRIX_IMPL_FB_HPP_
#define _SUPERNODAL_MATRIX_IMPL_FB_HPP_

template<typename T> Int SupernodalMatrix<T>::FBTaskAsyncRecv(Int iLocalI, SnodeUpdateFB & curTask, std::vector<AsyncComms> & incomingRecvAggArr, std::vector<AsyncComms * > & incomingRecvFactArr, IntNumVec & AggregatesToRecv, IntNumVec & FactorsToRecv){
    Int IrecvCnt = 0; 
    if(curTask.type == FACTOR){

      assert(curTask.remote_deps == AggregatesToRecv[curTask.tgt_snode_id-1]);
      AsyncComms & incomingRecvAgg = incomingRecvAggArr[iLocalI-1]; 
      //this is a factorization task: we have to receive aggregates
      Int maxRecvCnt = AggregatesToRecv[curTask.tgt_snode_id-1];
      Int tag = AGG_TAG(curTask.src_snode_id,curTask.tgt_snode_id);

//if(maxRecvCnt>0){
//gdb_lock();
//}
      
      while(maxRecvCnt>0 && incomingRecvCnt_ + IrecvCnt < maxIrecv_)
//      for(Int idx =0; idx<maxRecvCnt && incomingRecvCnt_ + IrecvCnt < maxIrecv_;
//         ++idx)
      {
        Int max_bytes = getAggBufSize<T>(curTask, Xsuper_, UpdateHeight_);


#ifdef PREALLOC_IRECV
        list<Icomm*>::iterator pos = availRecvBuffers_.begin();
        usedRecvBuffers_.splice(usedRecvBuffers_.begin(),availRecvBuffers_,pos);

        Icomm * bufPtr = *usedRecvBuffers_.begin();
        bufPtr->SetPos(usedRecvBuffers_.begin());
        bufPtr->reset();

        incomingRecvAgg.push_back(bufPtr);
#else
        incomingRecvAgg.push_back(new Icomm(max_bytes,MPI_REQUEST_NULL));
#endif

        Icomm & Irecv = *incomingRecvAgg.back();
        MPI_Irecv(Irecv.front(),Irecv.capacity(),MPI_BYTE,MPI_ANY_SOURCE,tag,CommEnv_->MPI_GetComm(),&Irecv.Request);
        ++IrecvCnt;
        maxRecvCnt--;
      }

//if(IrecvCnt>0){gdb_lock();}

      AggregatesToRecv[curTask.tgt_snode_id-1] -= IrecvCnt;
      assert(AggregatesToRecv[curTask.tgt_snode_id-1]>=0);
    }
    else{
      assert(curTask.remote_deps <= FactorsToRecv[curTask.tgt_snode_id-1]);
      Int maxRecvCnt = FactorsToRecv[curTask.tgt_snode_id-1];
//      if(maxRecvCnt==0){ return 0;}

      //this is an update task: we potentially have to receive factors
      if(incomingRecvFactArr[curTask.tgt_snode_id-1] == NULL){
        incomingRecvFactArr[curTask.tgt_snode_id-1] = new AsyncComms();
      }
      AsyncComms * incomingRecvFact = incomingRecvFactArr[curTask.tgt_snode_id-1];

      Int tag = FACT_TAG(curTask.src_snode_id,curTask.tgt_snode_id);

      while(maxRecvCnt>0 && incomingRecvCnt_ + IrecvCnt < maxIrecv_)
      //for(Int idx =0; idx<maxRecvCnt && incomingRecvCnt_ + IrecvCnt < maxIrecv_;
      //    ++idx)
      {
        Int max_bytes = getFactBufSize<T>(curTask, UpdateWidth_, UpdateHeight_);


#ifdef PREALLOC_IRECV
        list<Icomm*>::iterator pos = availRecvBuffers_.begin();
        usedRecvBuffers_.splice(usedRecvBuffers_.begin(),availRecvBuffers_,pos);

        Icomm * bufPtr = *usedRecvBuffers_.begin();
        bufPtr->SetPos(usedRecvBuffers_.begin());
        bufPtr->reset();

        incomingRecvFact->push_back(bufPtr);
#else
        incomingRecvFact->push_back(new Icomm(max_bytes,MPI_REQUEST_NULL));
#endif
        Icomm & Irecv = *incomingRecvFact->back();
        MPI_Irecv(Irecv.front(),Irecv.capacity(),MPI_BYTE,MPI_ANY_SOURCE,tag,CommEnv_->MPI_GetComm(),&Irecv.Request);
        ++IrecvCnt;
        maxRecvCnt--;
      }
      FactorsToRecv[curTask.tgt_snode_id-1] -= IrecvCnt;

//      if(FactorsToRecv[curTask.tgt_snode_id-1]<0){gdb_lock();}

      assert(FactorsToRecv[curTask.tgt_snode_id-1]>=0);
    }

    return IrecvCnt;
}

template<typename T> void SupernodalMatrix<T>::FBAsyncRecv(Int iLocalI, std::vector<AsyncComms> & incomingRecvAggArr, std::vector<AsyncComms * > & incomingRecvFactArr, IntNumVec & AggregatesToRecv, IntNumVec & FactorsToRecv){

      TIMER_START(ASYNC_RECV);
  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();
  std::vector<SnodeUpdateFB> tmpTasks;
  Int nextLocalI = iLocalI;

  Int iterCount = 0;
  volatile Int IrecvCnt = 0; 

  if(LocalTasks.size()>1){


  while(!LocalTasks.empty()){ 
      iterCount++;

if(iterCount>10){break;}

//      TIMER_START(ASYNC_RECV_TOP);
#ifdef _TASKLIST_
    SnodeUpdateFB curTask = LocalTasks.front();
#else
    SnodeUpdateFB curTask = LocalTasks.TOP();
#endif
//      TIMER_STOP(ASYNC_RECV_TOP);



    IrecvCnt = FBTaskAsyncRecv(nextLocalI, curTask, incomingRecvAggArr, incomingRecvFactArr, AggregatesToRecv, FactorsToRecv);
    if(curTask.type == FACTOR){
      nextLocalI++;
    }



//
//    IrecvCnt = 0; 
//    if(curTask.type == FACTOR){
//
//      AsyncComms & incomingRecvAgg = incomingRecvAggArr[nextLocalI-1]; 
//      nextLocalI++;
//      //this is a factorization task: we have to receive aggregates
//      Int maxRecvCnt = AggregatesToRecv[curTask.tgt_snode_id-1];
//      Int tag = AGG_TAG(curTask.src_snode_id,curTask.tgt_snode_id);
//
////if(maxRecvCnt>0){
////gdb_lock();
////}
//      for(Int idx =0; idx<maxRecvCnt && incomingRecvCnt_ + IrecvCnt < maxIrecv_;
//          ++idx){
//        Int max_bytes = getAggBufSize<T>(curTask, Xsuper_, UpdateHeight_);
//
//
//#ifdef PREALLOC_IRECV
//        list<Icomm*>::iterator pos = availRecvBuffers_.begin();
//        usedRecvBuffers_.splice(usedRecvBuffers_.begin(),availRecvBuffers_,pos);
//
//        Icomm * bufPtr = *usedRecvBuffers_.begin();
//        bufPtr->SetPos(usedRecvBuffers_.begin());
//        bufPtr->reset();
//
//        incomingRecvAgg.push_back(bufPtr);
//#else
//        incomingRecvAgg.push_back(new Icomm(max_bytes,MPI_REQUEST_NULL));
//#endif
//
//        Icomm & Irecv = *incomingRecvAgg.back();
//#ifdef _SEPARATE_COMM_
//        MPI_Irecv(Irecv.front(),Irecv.capacity(),MPI_BYTE,MPI_ANY_SOURCE,tag,FBAggCommEnv_->MPI_GetComm(),&Irecv.Request);
//#else
//        MPI_Irecv(Irecv.front(),Irecv.capacity(),MPI_BYTE,MPI_ANY_SOURCE,tag,CommEnv_->MPI_GetComm(),&Irecv.Request);
//#endif
//        ++IrecvCnt;
//      }
//
////if(IrecvCnt>0){gdb_lock();}
//
//      AggregatesToRecv[curTask.tgt_snode_id-1] -= IrecvCnt;
//    }
//    else{
//
//      //this is an update task: we potentially have to receive factors
//      if(incomingRecvFactArr[curTask.tgt_snode_id-1] == NULL){
//        incomingRecvFactArr[curTask.tgt_snode_id-1] = new AsyncComms();
//      }
//      AsyncComms * incomingRecvFact = incomingRecvFactArr[curTask.tgt_snode_id-1];
//
//      Int maxRecvCnt = FactorsToRecv[curTask.tgt_snode_id-1];
//      Int tag = FACT_TAG(curTask.src_snode_id,curTask.tgt_snode_id);
//
//      for(Int idx =0; idx<maxRecvCnt && incomingRecvCnt_ + IrecvCnt < maxIrecv_;
//          ++idx){
//        Int max_bytes = getFactBufSize<T>(curTask, UpdateWidth_, UpdateHeight_);
//
//
//#ifdef PREALLOC_IRECV
//        list<Icomm*>::iterator pos = availRecvBuffers_.begin();
//        usedRecvBuffers_.splice(usedRecvBuffers_.begin(),availRecvBuffers_,pos);
//
//        Icomm * bufPtr = *usedRecvBuffers_.begin();
//        bufPtr->SetPos(usedRecvBuffers_.begin());
//        bufPtr->reset();
//
//        incomingRecvFact->push_back(bufPtr);
//#else
//        incomingRecvFact->push_back(new Icomm(max_bytes,MPI_REQUEST_NULL));
//#endif
//        Icomm & Irecv = *incomingRecvFact->back();
//        MPI_Irecv(Irecv.front(),Irecv.capacity(),MPI_BYTE,MPI_ANY_SOURCE,tag,CommEnv_->MPI_GetComm(),&Irecv.Request);
//        ++IrecvCnt;
//      }
//      FactorsToRecv[curTask.tgt_snode_id-1] -= IrecvCnt;
//
//      if(FactorsToRecv[curTask.tgt_snode_id-1]<0){gdb_lock();}
//
//      assert(FactorsToRecv[curTask.tgt_snode_id-1]>=0);
//    }

    incomingRecvCnt_+=IrecvCnt;

if(IrecvCnt>0){
logfileptr->OFS()<<"incomingRecvCnt_: "<<incomingRecvCnt_<<endl;
logfileptr->OFS()<<"IrecvCnt: "<<IrecvCnt<<endl;
logfileptr->OFS()<<"maxIrecv_: "<<maxIrecv_<<endl;
}

    if( incomingRecvCnt_ >= maxIrecv_){
      break;
    }

//      TIMER_START(ASYNC_RECV_POP);
#ifdef _TASKLIST_
    LocalTasks.pop_front();
#else
    LocalTasks.pop();
#endif
//      TIMER_STOP(ASYNC_RECV_POP);
//        TIMER_START(RESTORE_TASK_QUEUE);
    tmpTasks.push_back(curTask);
//        TIMER_STOP(RESTORE_TASK_QUEUE);
  }
//if(iterCount>0){
//logfileptr->OFS()<<"ITER COUNT: "<<iterCount<<endl;
//}
        TIMER_START(RESTORE_TASK_QUEUE);
  //put the tasks back into the task queue
  for(auto it = tmpTasks.begin(); it != tmpTasks.end(); it++){
    LocalTasks.push(*it);
  }
        TIMER_STOP(RESTORE_TASK_QUEUE);
  }

      TIMER_STOP(ASYNC_RECV);
}
















template <typename T> void SupernodalMatrix<T>::FBFactorizationTask(SnodeUpdateFB & curTask, Int iLocalI, IntNumVec & AggregatesDone,  IntNumVec & FactorsToRecv, IntNumVec & AggregatesToRecv, std::vector<char> & src_blocks,std::vector<AsyncComms> & incomingRecvAggArr, std::vector<AsyncComms * > & incomingRecvFactArr)
{
  scope_timer(FB_FACTORIZATION_TASK);

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

  
  //if possible issue new Irecv and do the updates in any order

  if( incomingRecvCnt_ < maxIrecv_){
    Int IrecvCnt = FBTaskAsyncRecv(iLocalI, curTask, incomingRecvAggArr, incomingRecvFactArr, AggregatesToRecv, FactorsToRecv);
    incomingRecvCnt_+=IrecvCnt;
#ifdef _DEBUG_
    if(IrecvCnt>0){
      logfileptr->OFS()<<"FACT Posting "<<IrecvCnt<<" additionnal Irecv"<<std::endl;
    }
#endif
  }

  //first wait for the Irecv that have been issued previously
  AsyncComms & cur_incomingRecv = incomingRecvAggArr[iLocalI-1];
  MPI_Status recv_status;




  TIMER_START(IRECV_MPI_AGG);
  AsyncComms::iterator it = WaitIncomingFactors(cur_incomingRecv, recv_status,outgoingSend);
  while( it != cur_incomingRecv.end() ){
    Icomm * curComm = *it;
    SuperNode<T> dist_src_snode;
    size_t read_bytes = Deserialize(curComm->front(),dist_src_snode);

#ifdef _STAT_COMM_
    maxAggregRecv_ = max(maxAggregRecv_,read_bytes);
    sizesAggregRecv_ += read_bytes;
    countAggregRecv_++;
#endif

    //Deserialize the number of aggregates
    //Int * aggregatesCnt = (Int *)(curComm->front()+read_bytes);


//#ifdef _DEBUG_
    logfileptr->OFS()<<"IRECV Supernode "<<dist_src_snode.Id()<<std::endl;
//#endif

    src_snode.Aggregate(dist_src_snode);

    //delete the request from the list
#ifdef PREALLOC_IRECV
    //put the icomm back into the avail list
    availRecvBuffers_.splice(availRecvBuffers_.begin(),usedRecvBuffers_,(*it)->GetPos());
    Icomm * bufPtr = *availRecvBuffers_.begin();
    bufPtr->SetPos(availRecvBuffers_.begin());
#endif
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

#ifdef _STAT_COMM_
    maxAggregRecv_ = max(maxAggregRecv_,read_bytes);
    sizesAggregRecv_ += read_bytes;
    countAggregRecv_++;
#endif

#ifdef _DEBUG_
    logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode.Id()<<std::endl;
#endif

    src_snode.Aggregate(dist_src_snode);
    AggregatesToRecv[src_snode.Id()-1] -=1;

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

  if(np>1){
    //Send my factor to my ancestors. 
    BolNumVec is_factor_sent(np);
    SetValue(is_factor_sent,false);
    BolNumVec is_skipped(np);
    SetValue(is_skipped,false);

    SnodeUpdate curUpdate;
    TIMER_START(FIND_UPDATED_ANCESTORS);
    while(src_snode.FindNextUpdate(curUpdate,Xsuper_,SupMembership_)){ 
      Int iTarget = this->Mapping_->Map(curUpdate.tgt_snode_id-1,src_snode.Id()-1);

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

#ifndef _LOAD_BALANCE
          AdvanceOutgoing(outgoingSend);
          SendDelayedMessagesUp(MsgToSend,outgoingSend,LocalTasks);
#endif

        }

      }
    }
    TIMER_STOP(FIND_UPDATED_ANCESTORS);
  }
  //process some of the delayed send
  //       AdvanceOutgoing(outgoingSend);
  //       SendDelayedMessagesUp(MsgToSend,outgoingSend,LocalTasks);
}


template <typename T> void SupernodalMatrix<T>::FBUpdateTask(SnodeUpdateFB & curTask, IntNumVec & UpdatesToDo, IntNumVec & AggregatesDone,std::vector< SuperNode<T> * > & aggVectors, std::vector<char> & src_blocks,std::vector<AsyncComms> & incomingRecvAggArr, std::vector<AsyncComms * > & incomingRecvFactArr,  IntNumVec & FactorsToRecv, IntNumVec & AggregatesToRecv)
{
  scope_timer(FB_UPDATE_TASK);
  Int src_snode_id = curTask.src_snode_id;
  Int tgt_snode_id = curTask.tgt_snode_id;

  src_snode_id = abs(src_snode_id);
  bool is_first_local = curTask.src_snode_id <0;
  curTask.src_snode_id = src_snode_id;

  SuperNode<T> * cur_src_snode; 

  src_blocks.resize(0);

  Int iExpectedSrcOwner = this->Mapping_->Map(abs(curTask.src_snode_id)-1,abs(curTask.src_snode_id)-1);


  {

    //launch Irecv ?
  if( incomingRecvCnt_ < maxIrecv_){
    Int IrecvCnt = FBTaskAsyncRecv(-1, curTask, incomingRecvAggArr, incomingRecvFactArr, AggregatesToRecv, FactorsToRecv);
    incomingRecvCnt_+=IrecvCnt;
#ifdef _DEBUG_
    if(IrecvCnt>0){
      logfileptr->OFS()<<"AGG Posting "<<IrecvCnt<<" additionnal Irecv"<<std::endl;
    }
#endif
  }

  AsyncComms::iterator it;
  AsyncComms * cur_incomingRecv = incomingRecvFactArr[tgt_snode_id-1];

      cur_src_snode = FBRecvFactor(curTask,src_blocks,cur_incomingRecv,it,FactorsToRecv);

      /*if(cur_src_snode!=NULL)*/{

        //need to be updated because we might have received from someone else
#ifdef _DEBUG_PROGRESS_
        if(abs(curTask.src_snode_id) != src_snode_id){ 
          cout<<"YOUHOU WE HAVE ASYNC HERE !!! expected: "<< abs(curTask.src_snode_id) << " vs received: "<<cur_src_snode->Id()<<endl;
          logfileptr->OFS()<<"YOUHOU WE HAVE ASYNC HERE !!! expected: "<< abs(curTask.src_snode_id) << " vs received: "<<cur_src_snode->Id()<<endl;
        }
#endif
        Int iSrcOwner = this->Mapping_->Map(cur_src_snode->Id()-1,cur_src_snode->Id()-1);


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
          Int iUpdater = this->Mapping_->Map(curUpdate.tgt_snode_id-1,cur_src_snode->Id()-1);
          if(iUpdater == iam){

#ifdef _DEBUG_PROGRESS_
            logfileptr->OFS()<<"implicit Task: {"<<curUpdate.src_snode_id<<" -> "<<curUpdate.tgt_snode_id<<"}"<<std::endl;
            logfileptr->OFS()<<"Processing update from Supernode "<<curUpdate.src_snode_id<<" to Supernode "<<curUpdate.tgt_snode_id<<endl;
#endif
            SuperNode<T> * tgt_aggreg;

            Int iTarget = this->Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
            if(iTarget == iam){
              //the aggregate vector is directly the target snode
              tgt_aggreg = snodeLocal(curUpdate.tgt_snode_id);
              assert(curUpdate.tgt_snode_id == tgt_aggreg->Id());
            }
            else{
              //Check if src_snode_id already have an aggregate vector
              if(AggregatesDone[curUpdate.tgt_snode_id-1]==0){
                aggVectors[curUpdate.tgt_snode_id-1] = new SuperNode<T>(curUpdate.tgt_snode_id, Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id]-1,iSize_);
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
                logfileptr->OFS()<<"Remote Supernode "<<curUpdate.tgt_snode_id<<" is updated by Supernode "<<cur_src_snode->Id()<<std::endl;
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


#ifndef _LOAD_BALANCE
                AdvanceOutgoing(outgoingSend);
                SendDelayedMessagesUp(MsgToSend,outgoingSend,LocalTasks);
#endif


              }
            }

            //if local update, push a new task in the queue and stop the while loop
            if(iam==iSrcOwner ){
              break;
            }
          }
        }
        TIMER_STOP(UPDATE_ANCESTORS);
      }


    //process some of the delayed send
    //  AdvanceOutgoing(outgoingSend);
    //  SendDelayedMessagesUp(MsgToSend,outgoingSend,LocalTasks);



    if(iExpectedSrcOwner!=iam){
      //if the recv was from an async recv, delete the request
      if(cur_incomingRecv != NULL){
        if(it != cur_incomingRecv->end()){
          //delete the request from the list
#ifdef PREALLOC_IRECV
        //put the icomm back into the avail list
        availRecvBuffers_.splice(availRecvBuffers_.begin(),usedRecvBuffers_,(*it)->GetPos());
        Icomm * bufPtr = *availRecvBuffers_.begin();
        bufPtr->SetPos(availRecvBuffers_.begin());
#endif
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
}


template <typename T> void SupernodalMatrix<T>::FanBoth()
{
  scope_timer(FACTORIZATION_FB);

  TIMER_START(FB_INIT);
  Real timeSta, timeEnd;

  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();

  IntNumVec UpdatesToDo;
  IntNumVec AggregatesToRecv;

  std::vector<AsyncComms> incomingRecvAggArr(LocalSupernodes_.size());
  std::vector<AsyncComms * > incomingRecvFactArr(Xsuper_.m(),NULL);
  incomingRecvCnt_ = 0;
 
   
  FBGetUpdateCount(UpdatesToDo,AggregatesToRecv);
  xlindx_.Clear();
  lindx_.Clear();



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


    Int iOwner = this->Mapping_->Map(I-1,I-1);
    if(iam==iOwner){
      SnodeUpdateFB curUpdate;
      curUpdate.type=FACTOR;
      curUpdate.src_snode_id = I;
      curUpdate.tgt_snode_id = I;
      curUpdate.remote_deps = AggregatesToRecv[curUpdate.tgt_snode_id-1];
#ifdef _TASKLIST_
      LocalTasks.push_back(curUpdate);
#else
      LocalTasks.push(curUpdate);
#endif

      //TODO have to replace this
      Int J = -1;
     
      bool is_first = true; 
      while( (J=FBUpdate(I,J))!= -1 ){
        SnodeUpdateFB curUpdate;
        curUpdate.type=AGGREGATE;
        //set it to be negative as it is a local update
//        FactorsToRecv[J-1]--;

        curUpdate.src_snode_id = is_first?I:-I;
        curUpdate.tgt_snode_id = J;
#ifdef _TASKLIST_
        LocalTasks.push_back(curUpdate);
#else
        LocalTasks.push(curUpdate);
#endif
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
#ifdef _TASKLIST_
        LocalTasks.push_back(curUpdate);
#else
        LocalTasks.push(curUpdate);
#endif

        curUpdate.remote_deps = 1;
        FactorsToRecv[J-1]++;
      }
    }
  }

#ifdef _TASKLIST_
  //sort the task queue
  {
    SnodeUpdateFBCompare comp;
    //LocalTasks.sort(CompareSnodeUpdateFB);
    LocalTasks.sort(comp);
  }
#endif

  TIMER_STOP(BUILD_TASK_LIST);

#ifdef _DEBUG_UPDATES_
  logfileptr->OFS()<<"UpdatesToDo: "<<UpdatesToDo<<endl;
  logfileptr->OFS()<<"AggregatesToRecv: "<<AggregatesToRecv<<endl;
  logfileptr->OFS()<<"FactorsToRecv: "<<FactorsToRecv<<endl;
#endif
  TIMER_STOP(FB_INIT);

  Int iLocalI = 1;
  while(!LocalTasks.empty() || !MsgToSend.empty() || !outgoingSend.empty()){

    //Check for completion of outgoing communication
    AdvanceOutgoing(outgoingSend);
    //We have to have a list of remote supernodes for which we are computing updates

    if(!LocalTasks.empty()){
      //make a copy because we need to pop it
      TIMER_START(POP_TASK);
#ifdef _TASKLIST_
      SnodeUpdateFB curTask = LocalTasks.front();
#else
      SnodeUpdateFB curTask = LocalTasks.TOP();
#endif
      TIMER_STOP(POP_TASK);

      DUMP_TASK_LIST();
#ifdef _DEBUG_PROGRESS_
        logfileptr->OFS()<<"picked Task: {"<<curTask.src_snode_id<<" -> "<<curTask.tgt_snode_id<<"}. Still "<<LocalTasks.size()<<" to do"<<std::endl;
#endif

      //Launch Irecv for subsequent local supernodes if I can
      //FBAsyncRecv(iLocalI,incomingRecvAggArr,incomingRecvFactArr,AggregatesToRecv, FactorsToRecv);

      //If it is a factorization
      if(curTask.type == FACTOR){
        FBFactorizationTask(curTask,iLocalI,AggregatesDone,FactorsToRecv,AggregatesToRecv, src_blocks,incomingRecvAggArr,incomingRecvFactArr);
        ++iLocalI; 
      }
      else{
        FBUpdateTask(curTask, UpdatesToDo, AggregatesDone, aggVectors, src_blocks,incomingRecvAggArr,incomingRecvFactArr,FactorsToRecv,AggregatesToRecv);
      }
      TIMER_START(POP_TASK);
#ifdef _TASKLIST_
      LocalTasks.pop_front();
#else
      LocalTasks.pop();
#endif
      TIMER_STOP(POP_TASK);
    }

    //process some of the delayed send
    AdvanceOutgoing(outgoingSend);
    SendDelayedMessagesUp(MsgToSend,outgoingSend,LocalTasks);

  }

  TIMER_START(SAFETY_ASSERT);
  for(Int idx = 0; idx <incomingRecvAggArr.size();++idx){
    assert(incomingRecvAggArr[idx].size()==0);
  } 

  bool babort=false;
  logfileptr->OFS()<<"CHECK: ";
  for(Int idx = 0; idx <incomingRecvFactArr.size();++idx){
    if(incomingRecvFactArr[idx]!=NULL){
      if(incomingRecvFactArr[idx]->size()!=0){logfileptr->OFS()<<idx<<" "; babort=true;}
//      if(incomingRecvFactArr[idx]->size()!=0){gdb_lock();}
//      assert(incomingRecvFactArr[idx]->size()==0);
    }
  }
  logfileptr->OFS()<<endl;
  if(babort){abort();}
 
  TIMER_STOP(SAFETY_ASSERT);

  TIMER_START(FINAL_BARRIER);
  MPI_Barrier(CommEnv_->MPI_GetComm());
  TIMER_STOP(FINAL_BARRIER);

  tmpBufs.Clear();

}


template <typename T> void SupernodalMatrix<T>::FBGetUpdateCount(IntNumVec & sc, IntNumVec & atr){
  scope_timer(FB_GET_UPDATE_COUNT);
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

    Idx64 fi = xlindx_[s-1];
    Idx64 li = xlindx_[s]-1;


//#ifdef _DEBUG_UPDATES_
//    logfileptr->OFS()<<"Supernode "<<s<<" updates: ";
//#endif


    for(Idx64 row_idx = fi; row_idx<=li;++row_idx){
      Idx32 row = lindx_[row_idx-1];
      Int supno = SupMembership_(row-1);

      if(marker(supno-1)!=s && supno!=s){


        Int iFactorizer = this->Mapping_->Map(supno-1,supno-1);
        Int iUpdater = this->Mapping_->Map(supno-1,s-1);

#ifdef _DEBUG_UPDATES_
          logfileptr->OFS()<<iFactorizer<<" "<<iUpdater<<endl;
#endif

        
        if(!isSent[(supno-1)*np+iUpdater] && iUpdater!=iFactorizer){
          atr[supno-1]++;
          isSent[(supno-1)*np+iUpdater]=true;
        }


        if(iam == iUpdater){

//#ifdef _DEBUG_UPDATES_
//          logfileptr->OFS()<<supno<<" ";
//#endif
          ++sc[supno-1];
        }
        else{

//#ifdef _DEBUG_UPDATES_
//          logfileptr->OFS()<<supno<<" [P"<<iUpdater<<"] ";
//#endif
        }

        marker(supno-1) = s;
      }
    }

//#ifdef _DEBUG_UPDATES_
//    logfileptr->OFS()<<std::endl;
//#endif

  }
}


template<typename T> SuperNode<T> * SupernodalMatrix<T>::FBRecvFactor(const SnodeUpdateFB & curTask, std::vector<char> & src_blocks,AsyncComms * cur_incomingRecv,AsyncComms::iterator & it, IntNumVec & FactorsToRecv)
{
  scope_timer(RECV_FACTORS);
  //TIMER_START(RECV_FACTORS);
  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();
  SuperNode<T> * cur_src_snode = NULL;
  Int iSrcOwner = this->Mapping_->Map(curTask.src_snode_id-1,curTask.src_snode_id-1);
  

  //This is a local update, iSrcOwner matters
  if(iSrcOwner==iam){
    cur_src_snode = snodeLocal(curTask.src_snode_id);
  }
  else /*if(FactorsToRecv[curTask.tgt_snode_id-1]>0)*/
  {

    Int tag = FACT_TAG(curTask.src_snode_id,curTask.tgt_snode_id);

    TIMER_START(RECV_MPI);
      SuperNode<T> * dist_src_snode = new SuperNode<T>();
      //Receive the factor
      //first wait for the Irecv
      bool do_blocking_recv = true;
#if 1
//      if(curTask.src_snode_id==1 && curTask.tgt_snode_id==2){gdb_lock();}
      if(cur_incomingRecv != NULL){
        TIMER_START(IRECV_MPI_FACT);
        MPI_Status recv_status;
        it = WaitIncomingFactors(*cur_incomingRecv, recv_status,outgoingSend);
        if( it != cur_incomingRecv->end() ){
          Icomm * curComm = *it;
          size_t read_bytes = Deserialize(curComm->front(),*dist_src_snode);

#ifdef _STAT_COMM_
          maxFactorsRecv_ = max(maxFactorsRecv_,read_bytes);
          sizesFactorsRecv_ += read_bytes;
          countFactorsRecv_++;
#endif
          do_blocking_recv = false;
        }

//      if( do_blocking_recv){
//        if(FactorsToRecv[curTask.tgt_snode_id-1]==0){
//assert(cur_incomingRecv->size()==0);
//        }
//      }

        TIMER_STOP(IRECV_MPI_FACT);
      }
#endif

      //do a blocking recv otherwise
      if( do_blocking_recv){

//if(FactorsToRecv[curTask.tgt_snode_id-1]<=0){gdb_lock();}
//logfileptr->OFS()<<FactorsToRecv[curTask.tgt_snode_id-1]<<endl;

assert(FactorsToRecv[curTask.tgt_snode_id-1]>0);
        MPI_Status recv_status;
        int bytes_received = 0;

        Int max_bytes = getFactBufSize<T>(curTask, UpdateWidth_, UpdateHeight_);

        TIMER_START(RECV_MALLOC);
        src_blocks.resize(max_bytes);
        TIMER_STOP(RECV_MALLOC);

        MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,MPI_ANY_SOURCE,tag,CommEnv_->MPI_GetComm(),&recv_status);
        //    logfileptr->OFS()<<"RECV After"<<endl;
        FactorsToRecv[curTask.tgt_snode_id-1]--;

        size_t read_bytes = Deserialize(&src_blocks[0],*dist_src_snode);

#ifdef _STAT_COMM_
          maxFactorsRecv_ = max(maxFactorsRecv_,read_bytes);
          sizesFactorsRecv_ += read_bytes;
          countFactorsRecv_++;
#endif

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
//  TIMER_STOP(RECV_FACTORS);
  return cur_src_snode;
}

template<typename T> Int SupernodalMatrix<T>::FBUpdate(Int I,Int prevJ){
  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();
  //Check if I have anything to update with that supernode
  //look at lindx_
  Idx64 fi = xlindx_(I-1);
  Idx64 li = xlindx_(I)-1;

  Int iOwner = this->Mapping_->Map(I-1,I-1);


  Int firstUpdate = -1; 
  Int J = -1; 
  for(Idx64 idx = fi; idx<=li;++idx){
    Idx32 row = lindx_[idx-1];
    J = SupMembership_[row-1];
    Int iUpdater = this->Mapping_->Map(J-1,I-1);


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
