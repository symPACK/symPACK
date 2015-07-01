#ifndef _SUPERNODAL_MATRIX_IMPL_FB_PULL_HPP_
#define _SUPERNODAL_MATRIX_IMPL_FB_PULL_HPP_


template <typename T> void SupernodalMatrix2<T>::FanBoth() {
  scope_timer(a,FACTORIZATION_FB);

  TIMER_START(FB_INIT);
  Real timeSta, timeEnd;

  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();

  incomingRecvCnt_ = 0;

   
  std::vector<Int> UpdatesToDo;
  std::vector<Int> AggregatesToRecv;
  std::vector<Int> LocalAggregates;
  FBGetUpdateCount(UpdatesToDo,AggregatesToRecv,LocalAggregates);
  xlindx_.Clear();
  lindx_.Clear();
  std::vector<Int> AggregatesDone(Xsuper_.m(),I_ZERO);

  //tmp buffer space
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
  std::vector< SuperNode2<T> * > aggVectors(Xsuper_.m()-1,NULL);



  timeSta =  get_time( );
  TIMER_START(BUILD_TASK_LIST);
  //resize taskLists_
  taskLists_.resize(Xsuper_.m(),NULL);

  Int localTaskCount =0;
  for(Int I = 1; I<Xsuper_.m(); ++I){
    Int iOwner = Mapping_->Map(I-1,I-1);


    if(iam==iOwner){
      FBTask curUpdate;
      curUpdate.type=FACTOR;
      curUpdate.src_snode_id = I;
      curUpdate.tgt_snode_id = I;
      curUpdate.remote_deps = AggregatesToRecv[curUpdate.tgt_snode_id-1];
      curUpdate.local_deps = LocalAggregates[curUpdate.tgt_snode_id-1];

      //create the list if needed
      if(taskLists_[I-1] == NULL){
        taskLists_[I-1]=new std::list<FBTask>();
      }

      taskLists_[I-1]->push_back(curUpdate);
      localTaskCount++;

//      //task might be ready
//      auto taskit = --taskLists_[I-1]->end();
//      if(taskit->remote_deps==0 && taskit->local_deps==0){
//        readyTasks_.push_back(taskit);    
//      }

    }


      //TODO Update tasks receiving REMOTE factors need only the FIRST task
      Int J = -1;
      while( (J=FBUpdate(I,J))!= -1 ){
        FBTask curUpdate;
        curUpdate.type=AGGREGATE;
        curUpdate.src_snode_id = I;
        curUpdate.tgt_snode_id = J;

        //If I own the factor, it is a local dependency
        if(iam==iOwner){
          curUpdate.remote_deps = 0;
          curUpdate.local_deps = 1;
        }
        else{
          curUpdate.remote_deps = 1;
          curUpdate.local_deps = 0;
        }


        //create the list if needed
        if(taskLists_[J-1] == NULL){
          taskLists_[J-1]=new std::list<FBTask>();
        }
        taskLists_[J-1]->push_back(curUpdate);
        localTaskCount++;

        //create only one task if I don't own the factor
        if(iam!=iOwner){
          break;
        }
      }
  }

  //sort the task queues
  {
    FBTaskCompare comp;
    for(int i = 0; i<taskLists_.size(); ++i){
      if(taskLists_[i] != NULL){
        taskLists_[i]->sort(comp);
        //now see if there are some ready tasks
        for(auto taskit = taskLists_[i]->begin(); taskit!=taskLists_[i]->end();taskit++){
          if(taskit->remote_deps==0 && taskit->local_deps==0){
            readyTasks_.push_back(taskit);    
          }
        }
      }
    }
  }

  TIMER_STOP(BUILD_TASK_LIST);

#ifdef _DEBUG_UPDATES_
  logfileptr->OFS()<<"UpdatesToDo: "<<UpdatesToDo<<endl;
  logfileptr->OFS()<<"AggregatesToRecv: "<<AggregatesToRecv<<endl;
  logfileptr->OFS()<<"LocalAggregates: "<<LocalAggregates<<endl;


  Int cnt = 0;
  logfileptr->OFS()<<"All Tasks: "<<endl;
  for(Int I = 1; I<Xsuper_.m(); ++I){
    logfileptr->OFS()<<I<<": "<<endl;
    if(taskLists_[I-1]!=NULL){
    for(auto taskit = taskLists_[I-1]->begin();
        taskit!=taskLists_[I-1]->end();
        taskit++){
      logfileptr->OFS()<<"   T("<<taskit->src_snode_id<<","<<taskit->tgt_snode_id<<") ";
      cnt++;
    }
    }
    logfileptr->OFS()<<endl;
  }
  assert(localTaskCount == cnt);

  logfileptr->OFS()<<"Ready Tasks: "<<endl;
  for(auto rtaskit = readyTasks_.begin(); rtaskit!=readyTasks_.end();rtaskit++){
    auto taskit = *rtaskit;
    
    logfileptr->OFS()<<"   T("<<taskit->src_snode_id<<","<<taskit->tgt_snode_id<<") ";
  }
  logfileptr->OFS()<<endl;

#endif
  TIMER_STOP(FB_INIT);

#if 1

bool doPrint = true;
Int prevCnt = -1;
while(localTaskCount>0){
  CheckIncomingMessages();

  if(!readyTasks_.empty()){
    //Pick a ready task
    auto taskit = readyTasks_.front();
    auto rtaskit = readyTasks_.begin();
    //auto taskit = *rtaskit;
    //process task
    FBTask & curTask = *taskit;
#ifdef _DEBUG_PROGRESS_
    logfileptr->OFS()<<"Processing T("<<taskit->src_snode_id<<","<<taskit->tgt_snode_id<<") "<<endl;
#endif
    Int iLocalTGT = snodeLocalIndex(curTask.tgt_snode_id);
    switch(curTask.type){
      case FACTOR:
        {
          FBFactorizationTask(curTask, iLocalTGT);
        }
        break;
      case AGGREGATE:
        {
          FBUpdateTask(curTask, UpdatesToDo, AggregatesDone,aggVectors, AggregatesToRecv, AggregatesToRecv,localTaskCount);
        }
        break;
defaut:
        {
          abort();
        }
        break;
    }

    //remove task
    taskLists_[curTask.tgt_snode_id-1]->erase(taskit);
    readyTasks_.erase(rtaskit);
    localTaskCount--;
  }

#if 0
  //dump what still needs to be done
  if(localTaskCount==prevCnt){
    if(doPrint){
    logfileptr->OFS()<<"=================================="<<endl;
    logfileptr->OFS()<<"Still to do: "<<endl;
    Int cnt = 0;
    for(Int I = 1; I<Xsuper_.m(); ++I){
      logfileptr->OFS()<<I<<": "<<endl;
      if(taskLists_[I-1]!=NULL){
        for(auto taskit = taskLists_[I-1]->begin();
            taskit!=taskLists_[I-1]->end();
            taskit++){
          logfileptr->OFS()<<"   T("<<taskit->src_snode_id<<","<<taskit->tgt_snode_id<<") ";
          cnt++;
        }
      }
      logfileptr->OFS()<<endl;
    }
    logfileptr->OFS()<<"=================================="<<endl;
    if(cnt==0){
      logfileptr->OFS()<<localTaskCount<<endl;
      break;
    }
    doPrint = false;
    }
  }
  else{
    doPrint = true;
  }
#endif
  prevCnt=localTaskCount;
}
#endif

upcxx::async_wait();
upcxx::barrier();

////  Int iLocalI = 1;
////  while(!LocalTasks.empty() || !MsgToSend.empty() || !outgoingSend.empty()){
////
////     if(!LocalTasks.empty()){
////      //Peek task
////      TIMER_START(POP_TASK);
////#ifdef _TASKLIST_
////      SnodeUpdateFB curTask = LocalTasks.front();
////#else
////      SnodeUpdateFB curTask = LocalTasks.TOP();
////#endif
////      TIMER_STOP(POP_TASK);
////
////      //check if all the incoming communications are present
////
////
////     }
////   
////    //Check for completion of incoming communication
////    upcxx::advance();
////    bool comm_found = false;
////    if(!gIncomingRecvAsync.empty()){
////
////      //find if there is some finished async comm
////      auto it = TestAsyncIncomingMessage();
////      if(it!=gIncomingRecvAsync.end()){
////        comm_found = true;
////
////        IncomingMessage * msg = *it;
////        msg->Wait(); 
////        logfileptr->OFS()<<"Received async msg from "<<msg->Sender()<<endl;
////        int * ptr = (int*)msg->GetLocalPtr();
////        logfileptr->OFS()<<"Message is {"<<ptr[0]<<"..."<<ptr[numint-1]<<"} "<<endl;
////
////        remote_delete(msg->GetRemotePtr());
////
////        delete msg;
////
////        gIncomingRecvAsync.erase(it);
////        numRecv--;
////      }
////    }
////
////    //We have to have a list of remote supernodes for which we are computing updates
////
////    if(!LocalTasks.empty()){
////      //make a copy because we need to pop it
////      TIMER_START(POP_TASK);
////#ifdef _TASKLIST_
////      SnodeUpdateFB curTask = LocalTasks.front();
////#else
////      SnodeUpdateFB curTask = LocalTasks.TOP();
////#endif
////      TIMER_STOP(POP_TASK);
////
////      DUMP_TASK_LIST();
////#ifdef _DEBUG_PROGRESS_
////        logfileptr->OFS()<<"picked Task: {"<<curTask.src_snode_id<<" -> "<<curTask.tgt_snode_id<<"}. Still "<<LocalTasks.size()<<" to do"<<std::endl;
////#endif
////
////      //Launch Irecv for subsequent local supernodes if I can
////      //FBAsyncRecv(iLocalI,incomingRecvAggArr,incomingRecvFactArr,AggregatesToRecv, FactorsToRecv);
////
////      //If it is a factorization
////      if(curTask.type == FACTOR){
////        FBFactorizationTask(curTask,iLocalI,AggregatesDone,FactorsToRecv,AggregatesToRecv, src_blocks,incomingRecvAggArr,incomingRecvFactArr);
////        ++iLocalI; 
////      }
////      else{
////        FBUpdateTask(curTask, UpdatesToDo, AggregatesDone, aggVectors, src_blocks,incomingRecvAggArr,incomingRecvFactArr,FactorsToRecv,AggregatesToRecv);
////      }
////      TIMER_START(POP_TASK);
////#ifdef _TASKLIST_
////      LocalTasks.pop_front();
////#else
////      LocalTasks.pop();
////#endif
////      TIMER_STOP(POP_TASK);
////    }
////
////    //process some of the delayed send
////    AdvanceOutgoing(outgoingSend);
////    SendDelayedMessagesUp(MsgToSend,outgoingSend,LocalTasks);
////
////  }
////
////  TIMER_START(SAFETY_ASSERT);
////  for(Int idx = 0; idx <incomingRecvAggArr.size();++idx){
////    assert(incomingRecvAggArr[idx].size()==0);
////  } 
////
////  bool babort=false;
////  logfileptr->OFS()<<"CHECK: ";
////  for(Int idx = 0; idx <incomingRecvFactArr.size();++idx){
////    if(incomingRecvFactArr[idx]!=NULL){
////      if(incomingRecvFactArr[idx]->size()!=0){logfileptr->OFS()<<idx<<" "; babort=true;}
//////      if(incomingRecvFactArr[idx]->size()!=0){gdb_lock();}
//////      assert(incomingRecvFactArr[idx]->size()==0);
////    }
////  }
////  logfileptr->OFS()<<endl;
////  if(babort){abort();}
//// 
////  TIMER_STOP(SAFETY_ASSERT);
////
////  TIMER_START(FINAL_BARRIER);
////  MPI_Barrier(CommEnv_->MPI_GetComm());
////  TIMER_STOP(FINAL_BARRIER);

  upcxx::barrier();
  tmpBufs.Clear();

}


template <typename T> void SupernodalMatrix2<T>::FBGetUpdateCount(std::vector<Int> & UpdatesToDo, std::vector<Int> & AggregatesToRecv,std::vector<Int> & LocalAggregates){
  scope_timer(a,FB_GET_UPDATE_COUNT);
  UpdatesToDo.resize(Xsuper_.m(),I_ZERO);
  AggregatesToRecv.resize(Xsuper_.m(),I_ZERO);
  LocalAggregates.resize(Xsuper_.m(),I_ZERO);


  std::vector<Int> marker(Xsuper_.m(),I_ZERO);
  std::vector<bool>isSent(Xsuper_.m()*np,false);

  for(Int s = 1; s<Xsuper_.m(); ++s){
    Int first_col = Xsuper_[s-1];
    Int last_col = Xsuper_[s]-1;

    Idx64 fi = xlindx_[s-1];
    Idx64 li = xlindx_[s]-1;


    for(Idx64 row_idx = fi; row_idx<=li;++row_idx){
      Idx32 row = lindx_[row_idx-1];
      Int supno = SupMembership_[row-1];

      if(marker[supno-1]!=s && supno!=s){


        Int iFactorizer = Mapping_->Map(supno-1,supno-1);
        Int iUpdater = Mapping_->Map(supno-1,s-1);

        if( iUpdater==iFactorizer){
          LocalAggregates[supno-1]++;
        }
        else{
          if(!isSent[(supno-1)*np+iUpdater]){
            AggregatesToRecv[supno-1]++;
            isSent[(supno-1)*np+iUpdater]=true;
          }
        }


        if(iam == iUpdater){
          ++UpdatesToDo[supno-1];
        }

        marker[supno-1] = s;
      }
    }
  }
}

template<typename T> Int SupernodalMatrix2<T>::FBUpdate(Int I,Int prevJ){
  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();
  //Check if I have anything to update with that supernode
  //look at lindx_
  Idx64 fi = xlindx_(I-1);
  Idx64 li = xlindx_(I)-1;

  Int iOwner = Mapping_->Map(I-1,I-1);


  Int firstUpdate = -1; 
  Int J = -1; 
  for(Idx64 idx = fi; idx<=li;++idx){
    Idx32 row = lindx_[idx-1];
    J = SupMembership_[row-1];
    Int iUpdater = Mapping_->Map(J-1,I-1);


    if(iUpdater == iam && J>I){
      //if(iUpdater==iOwner){
        if(J>prevJ){
          firstUpdate = J;
          break;
        }
      //}
      //else{
      //  if(J>prevJ){
      //    firstUpdate = J;
      //    break;
      //  }
      //}
    }
  }

  return firstUpdate;
}



template <typename T> void SupernodalMatrix2<T>::FBFactorizationTask(FBTask & curTask, Int iLocalI){
  scope_timer(a,FB_FACTORIZATION_TASK);

  Int src_snode_id = curTask.src_snode_id;
  Int tgt_snode_id = curTask.tgt_snode_id;
  Int I = src_snode_id;
  SuperNode2<T> & src_snode = *LocalSupernodes_[iLocalI -1];
  Int src_first_col = src_snode.FirstCol();
  Int src_last_col = src_snode.LastCol();

#ifdef _DEBUG_PROGRESS_
  logfileptr->OFS()<<"Processing Supernode "<<I<<std::endl;
#endif


  //Applying aggregates
  TIMER_START(APPLY_AGGREGATES);
  //TODO we might have to create a third AGGREGATE type of task
  for(auto msgit = curTask.data.begin();msgit!=curTask.data.end();msgit++){
    IncomingMessage * msgPtr = *msgit;
    assert(msgPtr->IsDone());
    char* dataPtr = msgPtr->GetLocalPtr();

//if( msgPtr->meta.src == 84  && msgPtr->meta.tgt == 86 ){gdb_lock();}

    SuperNode2<T> dist_src_snode(dataPtr,msgPtr->Size());
    dist_src_snode.InitIdxToBlk();

    src_snode.Aggregate(dist_src_snode);

    //ask for remote delete
    remote_delete(msgPtr->GetRemotePtr());

    delete msgPtr;
  }
  TIMER_STOP(APPLY_AGGREGATES);

#ifdef _DEBUG_PROGRESS_
  logfileptr->OFS()<<"Factoring Supernode "<<I<<std::endl;
#endif

  TIMER_START(FACTOR_PANEL);
  src_snode.Factorize();
  TIMER_STOP(FACTOR_PANEL);

  //Sending factors and update local tasks
    //Send my factor to my ancestors. 
    BolNumVec is_factor_sent(np);
    SetValue(is_factor_sent,false);

    SnodeUpdate curUpdate;
    TIMER_START(FIND_UPDATED_ANCESTORS);
    while(src_snode.FindNextUpdate(curUpdate,Xsuper_,SupMembership_)){ 
      Int iTarget = this->Mapping_->Map(curUpdate.tgt_snode_id-1,src_snode.Id()-1);

      if(iTarget != iam){
        if(!is_factor_sent[iTarget]){
          MsgMetadata meta;
         
//if(curUpdate.tgt_snode_id==29 && curUpdate.src_snode_id==28){gdb_lock(0);}
          //TODO Replace all this by a Serialize function
          NZBlockDesc2 & nzblk_desc = src_snode.GetNZBlockDesc(curUpdate.blkidx);
          Int local_first_row = curUpdate.src_first_row - nzblk_desc.GIndex;
    Int nzblk_cnt = src_snode.NZBlockCnt() - curUpdate.blkidx;
    Int nzval_cnt_ = src_snode.Size()*(src_snode.NRowsBelowBlock(curUpdate.blkidx)-local_first_row);
    T* nzval_ptr = src_snode.GetNZval(nzblk_desc.Offset) + local_first_row*src_snode.Size();

          upcxx::global_ptr<char> sendPtr((char*)nzval_ptr);
          //the size of the message is the number of bytes between sendPtr and the address of nzblk_desc
          
          meta.src = curUpdate.src_snode_id;
          meta.tgt = curUpdate.tgt_snode_id;
          meta.GIndex = curUpdate.src_first_row;

          char * last_byte_ptr = (char*)&nzblk_desc + sizeof(NZBlockDesc2);
          size_t msgSize = last_byte_ptr - (char*)nzval_ptr;
          signal_data(sendPtr, msgSize, iTarget, meta);
          is_factor_sent[iTarget] = true;
        }
      }
      else{
        //Update local tasks
        //find task corresponding to curUpdate

        auto taskit = find_task(curUpdate.src_snode_id,curUpdate.tgt_snode_id);
        taskit->local_deps--;
        if(taskit->remote_deps==0 && taskit->local_deps==0){
          readyTasks_.push_back(taskit);    
        }
      }
    }
    TIMER_STOP(FIND_UPDATED_ANCESTORS);



}


template <typename T> std::list<FBTask>::iterator SupernodalMatrix2<T>::find_task(Int src, Int tgt){
        //find task corresponding to curUpdate
        for(auto taskit = taskLists_[tgt-1]->begin();
            taskit!=taskLists_[tgt-1]->end();
            taskit++){
          if(taskit->src_snode_id==src && taskit->tgt_snode_id==tgt){
            return taskit;
          }
        }
abort();
        return taskLists_[tgt-1]->end();
}

template <typename T> void SupernodalMatrix2<T>::FBUpdateTask(FBTask & curTask, std::vector<Int> & UpdatesToDo, std::vector<Int> & AggregatesDone,std::vector< SuperNode2<T> * > & aggVectors,  std::vector<Int> & FactorsToRecv, std::vector<Int> & AggregatesToRecv,Int & localTaskCount)
{
  scope_timer(a,FB_UPDATE_TASK);
  Int src_snode_id = curTask.src_snode_id;
  Int tgt_snode_id = curTask.tgt_snode_id;

  src_snode_id = abs(src_snode_id);
  bool is_first_local = curTask.src_snode_id <0;
  curTask.src_snode_id = src_snode_id;

  SuperNode2<T> * cur_src_snode; 


  Int iSrcOwner = this->Mapping_->Map(abs(curTask.src_snode_id)-1,abs(curTask.src_snode_id)-1);
  //if(iSrcOwner!=iam){gdb_lock();}

  {
    //AsyncComms::iterator it;
    //AsyncComms * cur_incomingRecv = incomingRecvFactArr[tgt_snode_id-1];


    IncomingMessage * msgPtr = NULL;
    //Local or remote factor
    //we have only one local or one remote incoming aggregate
    if(curTask.data.size()==0){
      cur_src_snode = snodeLocal(curTask.src_snode_id);
    }
    else{
      auto msgit = curTask.data.begin();
      msgPtr = *msgit;
      assert(msgPtr->IsDone());
      char* dataPtr = msgPtr->GetLocalPtr();

      cur_src_snode = new SuperNode2<T>(dataPtr,msgPtr->Size(),msgPtr->meta.GIndex);
      cur_src_snode->InitIdxToBlk();
      //TODO this has to be written again for the new supernode type
      //size_t read_bytes = Deserialize(dataPtr,*cur_src_snode);
    }





//    cur_src_snode = FBRecvFactor(curTask,src_blocks,cur_incomingRecv,it,FactorsToRecv);



    //Update everything src_snode_id own with that factor
    //Update the ancestors
    SnodeUpdate curUpdate;
    TIMER_START(UPDATE_ANCESTORS);
    while(cur_src_snode->FindNextUpdate(curUpdate,Xsuper_,SupMembership_,iam==iSrcOwner)){

      //skip if this update is "lower"
      if(curUpdate.tgt_snode_id<curTask.tgt_snode_id){continue;}

      Int iUpdater = this->Mapping_->Map(curUpdate.tgt_snode_id-1,cur_src_snode->Id()-1);
      if(iUpdater == iam){
#ifdef _DEBUG_PROGRESS_
        logfileptr->OFS()<<"implicit Task: {"<<curUpdate.src_snode_id<<" -> "<<curUpdate.tgt_snode_id<<"}"<<std::endl;
        logfileptr->OFS()<<"Processing update from Supernode "<<curUpdate.src_snode_id<<" to Supernode "<<curUpdate.tgt_snode_id<<endl;
#endif

        SuperNode2<T> * tgt_aggreg;

        Int iTarget = this->Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
        if(iTarget == iam){
          //the aggregate vector is directly the target snode
          tgt_aggreg = snodeLocal(curUpdate.tgt_snode_id);
          assert(curUpdate.tgt_snode_id == tgt_aggreg->Id());
        }
        else{
          //Check if src_snode_id already have an aggregate vector
          if(AggregatesDone[curUpdate.tgt_snode_id-1]==0){
            //use number of rows below factor as initializer
            Int iWidth =Xsuper_[curUpdate.tgt_snode_id] - Xsuper_[curUpdate.tgt_snode_id-1]; 
//            Int blkidx = cur_src_snode->FindBlockIdx(Xsuper_[curUpdate.tgt_snode_id-1]);
//            if(blkidx < cur_src_snode->NZBlockCnt()){
//
//              Int blkidx2 = cur_src_snode->FindBlockIdx(Xsuper_[curUpdate.tgt_snode_id]-1);
//              if(blkidx!=blkidx2){
//              }
//
//            }


            aggVectors[curUpdate.tgt_snode_id-1] = new SuperNode2<T>(curUpdate.tgt_snode_id, Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id]-1, iWidth, iSize_);
          }
          tgt_aggreg = aggVectors[curUpdate.tgt_snode_id-1];
        }

#ifdef _DEBUG_
        logfileptr->OFS()<<"RECV Supernode "<<curUpdate.tgt_snode_id<<" is updated by Supernode "<<cur_src_snode->Id()<<" rows "<<curUpdate.src_first_row<<" "<<curUpdate.blkidx<<std::endl;
#endif


        //Update the aggregate
        tgt_aggreg->UpdateAggregate(*cur_src_snode,curUpdate,tmpBufs,iTarget);


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

            tgt_aggreg->Shrink();

            MsgMetadata meta;


//if(curUpdate.tgt_snode_id==86 && curUpdate.src_snode_id==84){gdb_lock();}

            NZBlockDesc2 & nzblk_desc = tgt_aggreg->GetNZBlockDesc(0);
            T* nzval_ptr = tgt_aggreg->GetNZval(0);

            //this is an aggregate
            meta.src = curUpdate.src_snode_id;
            meta.tgt = curUpdate.tgt_snode_id;
            meta.GIndex = nzblk_desc.GIndex;


            upcxx::global_ptr<char> sendPtr = tgt_aggreg->GetGlobalPtr(meta.GIndex);
            //the size of the message is the number of bytes between sendPtr and the address of nzblk_desc
            size_t msgSize = tgt_aggreg->StorageSize();
            signal_data(sendPtr, msgSize, iTarget, meta);

          }
        }

#if 0
        //IF IT IS AN IMPLICIT TASK, MAYBE IT WAS NOT WORTH CREATING IT IN THE FIRST PLACE
        if(curTask.tgt_snode_id!=curUpdate.tgt_snode_id){
          auto taskit = find_task(curUpdate.src_snode_id, curUpdate.tgt_snode_id);
          taskLists_[curUpdate.tgt_snode_id-1]->erase(taskit);
          localTaskCount--;

//          for(auto taskit = taskLists_[curUpdate.tgt_snode_id-1]->begin();
//              taskit!=taskLists_[curUpdate.tgt_snode_id-1]->end();
//              taskit++){
//            if(taskit->src_snode_id==curUpdate.src_snode_id
//                && taskit->tgt_snode_id==curUpdate.tgt_snode_id){
//              taskLists_[curUpdate.tgt_snode_id-1]->erase(taskit); 
//              localTaskCount--; 
//              break;
//            }
//          }

        }
#endif

        if(iTarget == iam)
        {
          //update the dependency of the factorization
          auto taskit = find_task(curUpdate.tgt_snode_id,curUpdate.tgt_snode_id);
          taskit->local_deps--;
          if(taskit->remote_deps==0 && taskit->local_deps==0){
            readyTasks_.push_back(taskit);    
          }
        }

        //if local update, push a new task in the queue and stop the while loop
        if(iam==iSrcOwner ){
          break;
        }

      }
    }
    TIMER_STOP(UPDATE_ANCESTORS);



    if(curTask.data.size()>0){
      delete cur_src_snode;
      auto msgit = curTask.data.begin();
      IncomingMessage * msgPtr = *msgit;


      delete msgPtr;
    }





  }
}


template <typename T> void SupernodalMatrix2<T>::CheckIncomingMessages(){
//return;

  //call advance
    upcxx::advance();
    bool comm_found = false;
    
    IncomingMessage * msg = NULL;

      //find if there is some finished async comm
      auto it = TestAsyncIncomingMessage();
      if(it!=gIncomingRecvAsync.end()){
        msg = *it;
        gIncomingRecvAsync.erase(it);
      }
      else if(!gIncomingRecv.empty()){
        auto it = gIncomingRecv.begin();
        msg = *it;
        gIncomingRecv.erase(it);
      }

      if(msg!=NULL){
        msg->Wait(); 

        //int * ptr = (int*)msg->GetLocalPtr();
        //update task dependency
        //find task corresponding to curUpdate
//if(msg->meta.src==28 && msg->meta.tgt==29){gdb_lock(1);}

        Int iOwner = Mapping_->Map(msg->meta.tgt-1,msg->meta.tgt-1);
        Int iUpdater = Mapping_->Map(msg->meta.tgt-1,msg->meta.src-1);
        //this is an aggregate
        std::list<FBTask>::iterator taskit;
        if(iOwner==iam && iUpdater!=iam){
          //need to update dependencies of the factorization task
          taskit = find_task(msg->meta.tgt,msg->meta.tgt);
        }
        else{
          taskit = find_task(msg->meta.src,msg->meta.tgt);
        } 

        taskit->remote_deps--;
        taskit->data.push_back(msg);
        if(taskit->remote_deps==0 && taskit->local_deps==0){
          readyTasks_.push_back(taskit);    
        }
      }
       

}


#endif
