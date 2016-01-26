#ifndef _SUPERNODAL_MATRIX_IMPL_FB_PULL_HPP_
#define _SUPERNODAL_MATRIX_IMPL_FB_PULL_HPP_

template <typename T> void SupernodalMatrix2<T>::FanBoth_Static() {
  TIMER_START(FACTORIZATION_FB);

  TIMER_START(FB_INIT);
  Real timeSta, timeEnd;

  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();

  incomingRecvCnt_ = 0;


  std::vector<Int> UpdatesToDo;
  std::vector<Int> AggregatesToRecv;
  std::vector<Int> LocalAggregates;
  FBGetUpdateCount(UpdatesToDo,AggregatesToRecv,LocalAggregates);
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


  //build level structure
  std::vector<Int> levels(Xsuper_.m());
  levels[Xsuper_.m()-1]=-1;
  Int numLevel = 0; 
  for(Int i=Xsuper_.m()-1-1; i>=0; i-- ){ 
    levels[i] = levels[ETree_.PostParent(i)-1]+1;
  }




  localTaskCount_ =0;
  for(Int I = 1; I<Xsuper_.m(); ++I){
    Int iOwner = Mapping_->Map(I-1,I-1);


    if(iam==iOwner){
      FBTask curUpdate;
      curUpdate.type=FACTOR;
      curUpdate.src_snode_id = I;
      curUpdate.tgt_snode_id = I;
      curUpdate.remote_deps = AggregatesToRecv[curUpdate.tgt_snode_id-1];
      curUpdate.local_deps = LocalAggregates[curUpdate.tgt_snode_id-1];

      //{
      //  SuperNode2<T> & src_snode = *snodeLocal(curUpdate.src_snode_id);
      //  curUpdate.rank = 2.0*pow((double)src_snode.Size(),3.0)/3.0 + src_snode.Size()*(src_snode.Size()+1.0)*(src_snode.NRowsBelowBlock(0)-src_snode.Size())/2.0;
      //}
      curUpdate.rank = levels[curUpdate.src_snode_id-1];

      //create the list if needed
      if(taskLists_[I-1] == NULL){
        taskLists_[I-1]=new std::list<FBTask>();
      }

      taskLists_[I-1]->push_back(curUpdate);
      localTaskCount_++;

      //      //task might be ready
      //      auto taskit = --taskLists_[I-1]->end();
      //      if(taskit->remote_deps==0 && taskit->local_deps==0){
      //        readyTasks_.push_back(taskit);    
      //      }

    }

    //TODO create tasks only if local factor
    //TODO Update tasks receiving REMOTE factors need only the FIRST task
#ifdef DYN_TASK_CREATE
    if(iam==iOwner){
#endif
      Int J = -1;
      while( (J=FBUpdate(I,J))!= -1 ){
        FBTask curUpdate;
        curUpdate.type=UPDATE;
        curUpdate.src_snode_id = I;
        curUpdate.tgt_snode_id = J;
        curUpdate.rank = levels[curUpdate.src_snode_id-1];

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
        localTaskCount_++;

        //create only one task if I don't own the factor
        if(iam!=iOwner){
          break;
        }
      }
#ifdef DYN_TASK_CREATE
    }
#endif
  }

  xlindx_.Clear();
  lindx_.Clear();

  //sort the task queues
  {
    for(int i = 0; i<taskLists_.size(); ++i){
      if(taskLists_[i] != NULL){
        //        taskLists_[i]->sort(comp);
        //now see if there are some ready tasks
        for(auto taskit = taskLists_[i]->begin(); taskit!=taskLists_[i]->end();taskit++){


            //compute cost 
            taskit->rank = levels[taskit->src_snode_id-1];
            //            taskit->update_rank();
            scheduler_->push(taskit);

        }
      }
    }
  }

  TIMER_STOP(BUILD_TASK_LIST);

#ifdef _DEBUG_UPDATES_
  //<<<<<<< HEAD
  //  logfileptr->OFS()<<"UpdatesToDo: "<<UpdatesToDo<<endl;
  //  logfileptr->OFS()<<"AggregatesToRecv: "<<AggregatesToRecv<<endl;
  //  logfileptr->OFS()<<"LocalAggregates: "<<LocalAggregates<<endl;
  //
  //
  //  Int cnt = 0;
  //  logfileptr->OFS()<<"All Tasks: "<<endl;
  //  for(Int I = 1; I<Xsuper_.m(); ++I){
  //    logfileptr->OFS()<<I<<": "<<endl;
  //    if(taskLists_[I-1]!=NULL){
  //      for(auto taskit = taskLists_[I-1]->begin();
  //          taskit!=taskLists_[I-1]->end();
  //          taskit++){
  //        logfileptr->OFS()<<"   T("<<taskit->src_snode_id<<","<<taskit->tgt_snode_id<<") ";
  //        cnt++;
  //      }
  //    }
  //    logfileptr->OFS()<<endl;
  //  }
  //  assert(localTaskCount_ == cnt);
  //
  //  logfileptr->OFS()<<"Ready Tasks: "<<endl;
  //  auto tmp = scheduler_->GetQueue();
  //  //std::priority_queue<std::list<FBTask>::iterator, vector<std::list<FBTask>::iterator>, FBTaskCompare > tmp = scheduler_->GetQueue();
  //  while(!tmp.empty()){
  //    auto taskit = tmp.top();
  //    logfileptr->OFS()<<"   T("<<taskit->src_snode_id<<","<<taskit->tgt_snode_id<<") ";
  //    tmp.pop();
  //  }
  //  logfileptr->OFS()<<endl;
  //
  //=======
  ////  logfileptr->OFS()<<"UpdatesToDo: "<<UpdatesToDo<<endl;
  ////  logfileptr->OFS()<<"AggregatesToRecv: "<<AggregatesToRecv<<endl;
  ////  logfileptr->OFS()<<"LocalAggregates: "<<LocalAggregates<<endl;
  ////
  ////
  ////  Int cnt = 0;
  ////  logfileptr->OFS()<<"All Tasks: "<<endl;
  ////  for(Int I = 1; I<Xsuper_.m(); ++I){
  ////    logfileptr->OFS()<<I<<": "<<endl;
  ////    if(taskLists_[I-1]!=NULL){
  ////      for(auto taskit = taskLists_[I-1]->begin();
  ////          taskit!=taskLists_[I-1]->end();
  ////          taskit++){
  ////        logfileptr->OFS()<<"   T("<<taskit->src_snode_id<<","<<taskit->tgt_snode_id<<") ";
  ////        cnt++;
  ////      }
  ////    }
  ////    logfileptr->OFS()<<endl;
  ////  }
  ////  assert(localTaskCount_ == cnt);
  ////
  ////  logfileptr->OFS()<<"Ready Tasks: "<<endl;
  ////  std::priority_queue<std::list<FBTask>::iterator, vector<std::list<FBTask>::iterator>, FBTaskCompare > tmp = scheduler_->GetQueue();
  ////  while(!tmp.empty()){
  ////    auto taskit = tmp.top();
  ////    logfileptr->OFS()<<"   T("<<taskit->src_snode_id<<","<<taskit->tgt_snode_id<<") ";
  ////    tmp.pop();
  ////  }
  ////  logfileptr->OFS()<<endl;
  ////
  //>>>>>>> 1a872570a77db18b529e4cf8848b080f3f96dd77
#endif
  TIMER_STOP(FB_INIT);

#if 1

  bool doPrint = true;
  Int prevCnt = -1;
  while(localTaskCount_>0){
    //  if(scheduler_->done()){
    //CheckIncomingMessages();
    // }

    //logfileptr->OFS()<<"scheduler size: "<<scheduler_->size()<<" vs "<<localTaskCount_<<endl;
    if(!scheduler_->done()){
      //  TIMER_START(SORT_TASK);
      //    readyTasks_.sort(comp);
      //  TIMER_STOP(SORT_TASK);

      //Pick a ready task
      auto taskit = scheduler_->top();
      //auto taskit = *rtaskit;
      //process task
      FBTask & curTask = *taskit;


      while(taskit->remote_deps>0){
        CheckIncomingMessages(true);
      }

      scheduler_->pop();
      assert(taskit->local_deps==0);


      //    assert(find(taskLists_[curTask.tgt_snode_id-1]->begin(),taskLists_[curTask.tgt_snode_id-1]->end(),curTask)!=taskLists_[curTask.tgt_snode_id-1]->end());
#ifdef _DEBUG_PROGRESS_
      logfileptr->OFS()<<"Processing T("<<taskit->src_snode_id<<","<<taskit->tgt_snode_id<<") "<<endl;
#endif
      Int iLocalTGT = snodeLocalIndex(curTask.tgt_snode_id);
      switch(curTask.type){
        case FACTOR:
          {
            FBFactorizationTask(curTask, iLocalTGT,true);
          }
          break;
        case AGGREGATE:
          {
            FBAggregationTask(curTask, iLocalTGT,true);
          }
          break;
        case UPDATE:
          {
            FBUpdateTask(curTask, UpdatesToDo, AggregatesDone,aggVectors, AggregatesToRecv, AggregatesToRecv,localTaskCount_,true);
          }
          break;
defaut:
          {
            abort();
          }
          break;
      }

      TIMER_START(REMOVE_TASK);
      //remove task
      taskLists_[curTask.tgt_snode_id-1]->erase(taskit);
      localTaskCount_--;
      TIMER_STOP(REMOVE_TASK);
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
    prevCnt=localTaskCount_;
  }
#endif

  TIMER_START(BARRIER);
  upcxx::async_wait();
  team_->barrier();
  TIMER_STOP(BARRIER);

  tmpBufs.Clear();

  TIMER_STOP(FACTORIZATION_FB);
}



template <typename T> void SupernodalMatrix2<T>::FanBoth() {
  TIMER_START(FACTORIZATION_FB);

  TIMER_START(FB_INIT);
  Real timeSta, timeEnd;

  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();

  incomingRecvCnt_ = 0;


  std::vector<Int> UpdatesToDo;
  std::vector<Int> AggregatesToRecv;
  std::vector<Int> LocalAggregates;
  FBGetUpdateCount(UpdatesToDo,AggregatesToRecv,LocalAggregates);
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


  //build level structure
  std::vector<Int> levels(Xsuper_.m());
  levels[Xsuper_.m()-1]=-1;
  Int numLevel = 0; 
  for(Int i=Xsuper_.m()-1-1; i>=0; i-- ){ 
    levels[i] = levels[ETree_.PostParent(i)-1]+1;
  }




  localTaskCount_ =0;
  for(Int I = 1; I<Xsuper_.m(); ++I){
    Int iOwner = Mapping_->Map(I-1,I-1);


    if(iam==iOwner){
      FBTask curUpdate;
      curUpdate.type=FACTOR;
      curUpdate.src_snode_id = I;
      curUpdate.tgt_snode_id = I;
      curUpdate.remote_deps = AggregatesToRecv[curUpdate.tgt_snode_id-1];
      curUpdate.local_deps = LocalAggregates[curUpdate.tgt_snode_id-1];

      //{
      //  SuperNode2<T> & src_snode = *snodeLocal(curUpdate.src_snode_id);
      //  curUpdate.rank = 2.0*pow((double)src_snode.Size(),3.0)/3.0 + src_snode.Size()*(src_snode.Size()+1.0)*(src_snode.NRowsBelowBlock(0)-src_snode.Size())/2.0;
      //}
      curUpdate.rank = levels[curUpdate.src_snode_id-1];

      //create the list if needed
      if(taskLists_[I-1] == NULL){
        taskLists_[I-1]=new std::list<FBTask>();
      }

      taskLists_[I-1]->push_back(curUpdate);
      localTaskCount_++;

      //      //task might be ready
      //      auto taskit = --taskLists_[I-1]->end();
      //      if(taskit->remote_deps==0 && taskit->local_deps==0){
      //        readyTasks_.push_back(taskit);    
      //      }

    }

    //TODO create tasks only if local factor
    //TODO Update tasks receiving REMOTE factors need only the FIRST task
#ifdef DYN_TASK_CREATE
    if(iam==iOwner){
#endif
      Int J = -1;
      while( (J=FBUpdate(I,J))!= -1 ){
        FBTask curUpdate;
        curUpdate.type=UPDATE;
        curUpdate.src_snode_id = I;
        curUpdate.tgt_snode_id = J;
        curUpdate.rank = levels[curUpdate.src_snode_id-1];

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
        localTaskCount_++;

        //create only one task if I don't own the factor
        if(iam!=iOwner){
          break;
        }
      }
#ifdef DYN_TASK_CREATE
    }
#endif
  }

  xlindx_.Clear();
  lindx_.Clear();

  //sort the task queues
  {
    for(int i = 0; i<taskLists_.size(); ++i){
      if(taskLists_[i] != NULL){
        //        taskLists_[i]->sort(comp);
        //now see if there are some ready tasks
        for(auto taskit = taskLists_[i]->begin(); taskit!=taskLists_[i]->end();taskit++){
          if(taskit->remote_deps==0 && taskit->local_deps==0){
          //compute cost 
          taskit->rank = levels[taskit->src_snode_id-1];
          //            taskit->update_rank();
          scheduler_->push(taskit);
          }
        }
      }
    }
  }

  TIMER_STOP(BUILD_TASK_LIST);

#ifdef _DEBUG_UPDATES_
  //<<<<<<< HEAD
  //  logfileptr->OFS()<<"UpdatesToDo: "<<UpdatesToDo<<endl;
  //  logfileptr->OFS()<<"AggregatesToRecv: "<<AggregatesToRecv<<endl;
  //  logfileptr->OFS()<<"LocalAggregates: "<<LocalAggregates<<endl;
  //
  //
  //  Int cnt = 0;
  //  logfileptr->OFS()<<"All Tasks: "<<endl;
  //  for(Int I = 1; I<Xsuper_.m(); ++I){
  //    logfileptr->OFS()<<I<<": "<<endl;
  //    if(taskLists_[I-1]!=NULL){
  //      for(auto taskit = taskLists_[I-1]->begin();
  //          taskit!=taskLists_[I-1]->end();
  //          taskit++){
  //        logfileptr->OFS()<<"   T("<<taskit->src_snode_id<<","<<taskit->tgt_snode_id<<") ";
  //        cnt++;
  //      }
  //    }
  //    logfileptr->OFS()<<endl;
  //  }
  //  assert(localTaskCount_ == cnt);
  //
  //  logfileptr->OFS()<<"Ready Tasks: "<<endl;
  //  auto tmp = scheduler_->GetQueue();
  //  //std::priority_queue<std::list<FBTask>::iterator, vector<std::list<FBTask>::iterator>, FBTaskCompare > tmp = scheduler_->GetQueue();
  //  while(!tmp.empty()){
  //    auto taskit = tmp.top();
  //    logfileptr->OFS()<<"   T("<<taskit->src_snode_id<<","<<taskit->tgt_snode_id<<") ";
  //    tmp.pop();
  //  }
  //  logfileptr->OFS()<<endl;
  //
  //=======
  ////  logfileptr->OFS()<<"UpdatesToDo: "<<UpdatesToDo<<endl;
  ////  logfileptr->OFS()<<"AggregatesToRecv: "<<AggregatesToRecv<<endl;
  ////  logfileptr->OFS()<<"LocalAggregates: "<<LocalAggregates<<endl;
  ////
  ////
  ////  Int cnt = 0;
  ////  logfileptr->OFS()<<"All Tasks: "<<endl;
  ////  for(Int I = 1; I<Xsuper_.m(); ++I){
  ////    logfileptr->OFS()<<I<<": "<<endl;
  ////    if(taskLists_[I-1]!=NULL){
  ////      for(auto taskit = taskLists_[I-1]->begin();
  ////          taskit!=taskLists_[I-1]->end();
  ////          taskit++){
  ////        logfileptr->OFS()<<"   T("<<taskit->src_snode_id<<","<<taskit->tgt_snode_id<<") ";
  ////        cnt++;
  ////      }
  ////    }
  ////    logfileptr->OFS()<<endl;
  ////  }
  ////  assert(localTaskCount_ == cnt);
  ////
  ////  logfileptr->OFS()<<"Ready Tasks: "<<endl;
  ////  std::priority_queue<std::list<FBTask>::iterator, vector<std::list<FBTask>::iterator>, FBTaskCompare > tmp = scheduler_->GetQueue();
  ////  while(!tmp.empty()){
  ////    auto taskit = tmp.top();
  ////    logfileptr->OFS()<<"   T("<<taskit->src_snode_id<<","<<taskit->tgt_snode_id<<") ";
  ////    tmp.pop();
  ////  }
  ////  logfileptr->OFS()<<endl;
  ////
  //>>>>>>> 1a872570a77db18b529e4cf8848b080f3f96dd77
#endif
  TIMER_STOP(FB_INIT);

#if 1

  bool doPrint = true;
  Int prevCnt = -1;
  while(localTaskCount_>0){
    //  if(scheduler_->done()){
    CheckIncomingMessages();
    //  }

    if(!scheduler_->done()){
      //  TIMER_START(SORT_TASK);
      //    readyTasks_.sort(comp);
      //  TIMER_STOP(SORT_TASK);

      //Pick a ready task
      auto taskit = scheduler_->top();
      scheduler_->pop();
      //auto taskit = *rtaskit;
      //process task
      FBTask & curTask = *taskit;




      //    assert(find(taskLists_[curTask.tgt_snode_id-1]->begin(),taskLists_[curTask.tgt_snode_id-1]->end(),curTask)!=taskLists_[curTask.tgt_snode_id-1]->end());
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
            FBAggregationTask(curTask, iLocalTGT);
          }
          break;
        case UPDATE:
          {
            FBUpdateTask(curTask, UpdatesToDo, AggregatesDone,aggVectors, AggregatesToRecv, AggregatesToRecv,localTaskCount_);
          }
          break;
defaut:
          {
            abort();
          }
          break;
      }

      TIMER_START(REMOVE_TASK);
      //remove task
      taskLists_[curTask.tgt_snode_id-1]->erase(taskit);
      localTaskCount_--;
      TIMER_STOP(REMOVE_TASK);
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
    prevCnt=localTaskCount_;
  }
#endif

  TIMER_START(BARRIER);
  upcxx::async_wait();
  team_->barrier();
  TIMER_STOP(BARRIER);

  tmpBufs.Clear();

  TIMER_STOP(FACTORIZATION_FB);
}


template <typename T> void SupernodalMatrix2<T>::FBGetUpdateCount(std::vector<Int> & UpdatesToDo, std::vector<Int> & AggregatesToRecv,std::vector<Int> & LocalAggregates){
  TIMER_START(FB_GET_UPDATE_COUNT);
  UpdatesToDo.resize(Xsuper_.m(),I_ZERO);
  AggregatesToRecv.resize(Xsuper_.m(),I_ZERO);
  LocalAggregates.resize(Xsuper_.m(),I_ZERO);


  std::vector<Int> marker(Xsuper_.m(),I_ZERO);
  std::vector<bool>isSent(Xsuper_.m()*np,false);

  for(Int s = 1; s<Xsuper_.m(); ++s){
    Int first_col = Xsuper_[s-1];
    Int last_col = Xsuper_[s]-1;

    Ptr fi = xlindx_[s-1];
    Ptr li = xlindx_[s]-1;


    for(Ptr row_idx = fi; row_idx<=li;++row_idx){
      Idx row = lindx_[row_idx-1];
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
  TIMER_STOP(FB_GET_UPDATE_COUNT);
}

template<typename T> Int SupernodalMatrix2<T>::FBUpdate(Int I,Int prevJ){
  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();
  //Check if I have anything to update with that supernode
  //look at lindx_
  Ptr fi = xlindx_(I-1);
  Ptr li = xlindx_(I)-1;

  Int iOwner = Mapping_->Map(I-1,I-1);


  Int firstUpdate = -1; 
  Int J = -1; 
  for(Ptr idx = fi; idx<=li;++idx){
    Idx row = lindx_[idx-1];
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


template <typename T> void SupernodalMatrix2<T>::FBAggregationTask(FBTask & curTask, Int iLocalI, bool is_static){
  TIMER_START(FB_AGGREGATION_TASK);

#ifdef _DEBUG_PROGRESS_
  logfileptr->OFS()<<"Processing T_AGGREG("<<curTask.src_snode_id<<","<<curTask.tgt_snode_id<<") "<<endl;
#endif
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
  //TODO we might have to create a third UPDATE type of task
  //process them one by one
  assert(curTask.data.size()==1);
  for(auto msgit = curTask.data.begin();msgit!=curTask.data.end();msgit++){
    IncomingMessage * msgPtr = *msgit;
    assert(msgPtr->IsDone());
    char* dataPtr = msgPtr->GetLocalPtr();

    SuperNode2<T> dist_src_snode(dataPtr,msgPtr->Size());
    dist_src_snode.InitIdxToBlk();

    src_snode.Aggregate(dist_src_snode);

    if(msgPtr->IsLocal()){
      msgPtr->DeallocRemote();
    }

    delete msgPtr;
  }


  //need to update dependencies of the factorization task
  auto taskit = find_task(tgt_snode_id,tgt_snode_id,FACTOR);
  taskit->remote_deps--;

  if(!is_static){
    if(taskit->remote_deps==0 && taskit->local_deps==0){

      //    //compute cost 
      //    if(taskit->type==FACTOR){
      //      //TODO this is the factorization cost only. There is also the aggregation cost to take into account
      //      SuperNode2<T> & src_snode = *snodeLocal(taskit->src_snode_id);
      //      taskit->rank = 2.0*pow((double)src_snode.Size(),3.0)/3.0 + src_snode.Size()*(src_snode.Size()+1.0)*(src_snode.NRowsBelowBlock(0)-src_snode.Size())/2.0;
      //    }
      //    else if(taskit->type == UPDATE){
      //      //TODO loop through the source snode and compute cost for ONE target or ALL targets (if source is remote)
      //      taskit->rank = 0.0;
      //    }
      //compute cost 
      //taskit->update_rank();
      scheduler_->push(taskit);    
    }
  }
  TIMER_STOP(FB_AGGREGATION_TASK);
}



template <typename T> void SupernodalMatrix2<T>::FBFactorizationTask(FBTask & curTask, Int iLocalI, bool is_static){
  TIMER_START(FB_FACTORIZATION_TASK);

  Int src_snode_id = curTask.src_snode_id;
  Int tgt_snode_id = curTask.tgt_snode_id;
  Int I = src_snode_id;
  SuperNode2<T> & src_snode = *LocalSupernodes_[iLocalI -1];
  Int src_first_col = src_snode.FirstCol();
  Int src_last_col = src_snode.LastCol();

#ifdef _DEBUG_PROGRESS_
  logfileptr->OFS()<<"Processing Supernode "<<I<<std::endl;
#endif

#ifndef SEPARATE_AGGREG
  //Applying aggregates
  TIMER_START(APPLY_AGGREGATES);
  //TODO we might have to create a third UPDATE type of task
  for(auto msgit = curTask.data.begin();msgit!=curTask.data.end();msgit++){
    IncomingMessage * msgPtr = *msgit;
    assert(msgPtr->IsDone());
    char* dataPtr = msgPtr->GetLocalPtr();

    //if( msgPtr->meta.src == 84  && msgPtr->meta.tgt == 86 ){gdb_lock();}

    SuperNode2<T> dist_src_snode(dataPtr,msgPtr->Size());
    dist_src_snode.InitIdxToBlk();

    src_snode.Aggregate(dist_src_snode);


    delete msgPtr;
  }
  TIMER_STOP(APPLY_AGGREGATES);
#endif

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

        //Send factor 
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

      auto taskit = find_task(curUpdate.src_snode_id,curUpdate.tgt_snode_id,UPDATE);
      taskit->local_deps--;
      if(!is_static){
        if(taskit->remote_deps==0 && taskit->local_deps==0){

          //            //compute cost 
          //            if(taskit->type==FACTOR){
          //              //TODO this is the factorization cost only. There is also the aggregation cost to take into account
          //              SuperNode2<T> & src_snode = *snodeLocal(taskit->src_snode_id);
          //              taskit->rank = 2.0*pow((double)src_snode.Size(),3.0)/3.0 + src_snode.Size()*(src_snode.Size()+1.0)*(src_snode.NRowsBelowBlock(0)-src_snode.Size())/2.0;
          //            }
          //            else if(taskit->type == UPDATE){
          //              //TODO loop through the source snode and compute cost for ONE target or ALL targets (if source is remote)
          //              taskit->rank = 0.0;
          //            }

          //compute cost 
          //taskit->update_rank();
          scheduler_->push(taskit);    
        }
      }
    }
  }
  TIMER_STOP(FIND_UPDATED_ANCESTORS);



  TIMER_STOP(FB_FACTORIZATION_TASK);
}


template <typename T> std::list<FBTask>::iterator SupernodalMatrix2<T>::find_task(Int src, Int tgt, TaskType type ){
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

template <typename T> void SupernodalMatrix2<T>::FBUpdateTask(FBTask & curTask, std::vector<Int> & UpdatesToDo, std::vector<Int> & AggregatesDone,std::vector< SuperNode2<T> * > & aggVectors,  std::vector<Int> & FactorsToRecv, std::vector<Int> & AggregatesToRecv,Int & localTaskCount, bool is_static)
{
  TIMER_START(FB_UPDATE_TASK);
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
          if(/*AggregatesDone[curUpdate.tgt_snode_id-1]==0*/aggVectors[curUpdate.tgt_snode_id-1]==NULL){
            //use number of rows below factor as initializer
            Int iWidth =Xsuper_[curUpdate.tgt_snode_id] - Xsuper_[curUpdate.tgt_snode_id-1]; 

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
        //++AggregatesDone[curUpdate.tgt_snode_id-1];
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

            //            logfileptr->OFS()<<"Remote Supernode "<<curUpdate.tgt_snode_id<<" on P"<<iTarget<<" is updated by Supernode "<<cur_src_snode->Id()<<std::endl;

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
          auto taskit = find_task(curUpdate.tgt_snode_id,curUpdate.tgt_snode_id,FACTOR);
          taskit->local_deps--;
          if(!is_static){
            if(taskit->remote_deps==0 && taskit->local_deps==0){
              //            //compute cost 
              //            if(taskit->type==FACTOR){
              //              //TODO this is the factorization cost only. There is also the aggregation cost to take into account
              //              SuperNode2<T> & src_snode = *snodeLocal(taskit->src_snode_id);
              //              taskit->rank = 2.0*pow((double)src_snode.Size(),3.0)/3.0 + src_snode.Size()*(src_snode.Size()+1.0)*(src_snode.NRowsBelowBlock(0)-src_snode.Size())/2.0;
              //            }
              //            else if(taskit->type == UPDATE){
              //              //TODO loop through the source snode and compute cost for ONE target or ALL targets (if source is remote)
              //              taskit->rank = 0.0;
              //            }
              //compute cost 
              //taskit->update_rank();
              scheduler_->push(taskit);    
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



    if(curTask.data.size()>0){
      delete cur_src_snode;
      auto msgit = curTask.data.begin();
      IncomingMessage * msgPtr = *msgit;


      delete msgPtr;
    }





  }
  TIMER_STOP(FB_UPDATE_TASK);
}


template <typename T> void SupernodalMatrix2<T>::CheckIncomingMessages(bool is_static){
  TIMER_START(CHECK_MESSAGE);
  //return;

  //if(1 || upcxx::peek())
  {
    //call advance
    TIMER_START(UPCXX_ADVANCE);
    upcxx::advance();
    TIMER_STOP(UPCXX_ADVANCE);
    bool comm_found = false;

    IncomingMessage * msg = NULL;

    FBTask * curTask = NULL;
    if(is_static){
      curTask = &*scheduler_->top();
    }


    do{
      msg=NULL;
      //    if(!gIncomingRecv.empty()){
      TIMER_START(MV_MSG_SYNC);

      //if we have some room, turn blocking comms into async comms
      if(gIncomingRecvAsync.size() < gMaxIrecv || gMaxIrecv==-1){
        while((gIncomingRecvAsync.size() < gMaxIrecv || gMaxIrecv==-1) && !gIncomingRecv.empty()){
          if(!is_static){
          //auto it = gIncomingRecv.top();
          auto it = gIncomingRecv.front();
          (it)->AllocLocal();
          (it)->AsyncGet();
          gIncomingRecvAsync.push_back(it);
          //gIncomingRecv.pop();
          gIncomingRecv.pop_front();

          //auto it = gIncomingRecv.begin();
          //(*it)->AllocLocal();
          //(*it)->AsyncGet();
          //gIncomingRecvAsync.splice(gIncomingRecvAsync.end(),gIncomingRecv,it);
#ifdef _DEBUG_PROGRESS_
          logfileptr->OFS()<<"TRANSFERRED TO ASYNC COMM"<<endl;
#endif
          }
          else{


//              Int iOwner = Mapping_->Map(msg->meta.tgt-1,msg->meta.tgt-1);
//              Int iUpdater = Mapping_->Map(msg->meta.tgt-1,msg->meta.src-1);
//              if(curTask->type == FACTOR){
//                if(iOwner==iam && iUpdater!=iam){
//                  //this is an aggregate
//                  taskit = find_task(msg->meta.tgt,msg->meta.tgt,FACTOR);
//                }
//              }
//              else if(curTask->type == UPDATE){
//                if(iOwner!=iam || iUpdater==iam){
//                  //this is a factor
//                  taskit = find_task(msg->meta.src,msg->meta.tgt,UPDATE);
//                }
//              }




            for(auto it = gIncomingRecv.begin();it!=gIncomingRecv.end();it++){
              //find a message corresponding to current task
              IncomingMessage * curMsg = *it;
              if(curMsg->meta.tgt==curTask->tgt_snode_id){
                (curMsg)->AllocLocal();
                (curMsg)->AsyncGet();
                gIncomingRecvAsync.push_back(curMsg);
                //gIncomingRecv.pop();
                gIncomingRecv.erase(it);
                break;
              }
            }
          }
        }
      }
      TIMER_STOP(MV_MSG_SYNC);
      //    }

#ifdef HANDLE_LOCAL_POINTER
      //find if there is some local (SHMEM) comms
      if(!gIncomingRecvLocal.empty()){
        auto it = gIncomingRecvLocal.begin();
        msg = *it;
        gIncomingRecvLocal.erase(it);

#ifdef _DEBUG_PROGRESS_
        logfileptr->OFS()<<"COMM LOCAL: AVAILABLE MSG("<<msg->meta.src<<","<<msg->meta.tgt<<") from P"<<msg->remote_ptr.where()<<endl;
#endif
      }
      else
#endif
      {

        //find if there is some finished async comm
        auto it = TestAsyncIncomingMessage();
        if(it!=gIncomingRecvAsync.end()){
          TIMER_START(RM_MSG_ASYNC);
          msg = *it;
          gIncomingRecvAsync.erase(it);
          TIMER_STOP(RM_MSG_ASYNC);
#ifdef _DEBUG_PROGRESS_
          logfileptr->OFS()<<"COMM ASYNC: AVAILABLE MSG("<<msg->meta.src<<","<<msg->meta.tgt<<") from P"<<msg->remote_ptr.where()<<endl;
#endif
        }
        else if((is_static || scheduler_->done()) && !gIncomingRecv.empty()){
          //find a "high priority" task

#if 0
          {
            auto tmp = gIncomingRecv;
            while(!tmp.empty()){
              auto it = tmp.top();
              msg = it;
              tmp.pop();
              logfileptr->OFS()<<"COMM: AVAILABLE MSG("<<msg->meta.src<<","<<msg->meta.tgt<<") from P"<<msg->remote_ptr.where()<<endl;
            }
          }
#endif

          TIMER_START(RM_MSG_SYNC);
if(!is_static){
          //auto it = gIncomingRecv.top();
          auto it = gIncomingRecv.front();
          msg = it;
          //gIncomingRecv.pop();
          gIncomingRecv.pop_front();
          //auto it = gIncomingRecv.begin();
          //msg = *it;
          //gIncomingRecv.erase(it);
}
else{
            for(auto it = gIncomingRecv.begin();it!=gIncomingRecv.end();it++){
              //find a message corresponding to current task
              IncomingMessage * curMsg = *it;
              if(curMsg->meta.tgt==curTask->tgt_snode_id){
                msg = curMsg;
                gIncomingRecv.erase(it);
                break;
              }
            }
}
          TIMER_STOP(RM_MSG_SYNC);
        }

      }

      if(msg!=NULL){
        scope_timer(a,WAIT_AND_UPDATE_DEPS);
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
#ifdef _DEBUG_PROGRESS_
          logfileptr->OFS()<<"COMM: FETCHED AGGREG MSG("<<msg->meta.src<<","<<msg->meta.tgt<<") from P"<<iUpdater<<endl;
#endif

#ifdef SEPARATE_AGGREG
          //insert an aggregation task
          Int I = msg->meta.src;
          Int J = msg->meta.tgt;
          FBTask curUpdate;
          curUpdate.type=AGGREGATE;
          curUpdate.src_snode_id = I;
          curUpdate.tgt_snode_id = J;

          curUpdate.remote_deps = 1;
          curUpdate.local_deps = 0;

          if(taskLists_[J-1] == NULL){
            taskLists_[J-1]=new std::list<FBTask>();
          }
          taskLists_[J-1]->push_back(curUpdate);
          localTaskCount_++;
          taskit = --taskLists_[J-1]->end();

          {
            SuperNode2<T> & tgt_snode = *snodeLocal(taskit->tgt_snode_id);
            char* dataPtr = msg->GetLocalPtr();
            SuperNode2<T> dist_src_snode(dataPtr,msg->Size());
            taskit->rank = dist_src_snode.NRowsBelowBlock(0)*tgt_snode.Size();
          }
#else
          //need to update dependencies of the factorization task
          taskit = find_task(msg->meta.tgt,msg->meta.tgt,FACTOR);

#endif
        }
        else{
          Int iOwner = Mapping_->Map(msg->meta.src-1,msg->meta.src-1);
#ifdef _DEBUG_PROGRESS_
          logfileptr->OFS()<<"COMM: FETCHED FACTOR MSG("<<msg->meta.src<<","<<msg->meta.tgt<<") from P"<<iOwner<<endl;
#endif

#ifdef DYN_TASK_CREATE
          {
            FBTask curTask;
            curTask.type=UPDATE;
            curTask.src_snode_id = msg->meta.src;
            curTask.tgt_snode_id = msg->meta.tgt;

            //If I own the factor, it is a local dependency
            curTask.remote_deps = 1;
            curTask.local_deps = 0;

            //create the list if needed
            Int J = msg->meta.tgt;
            if(taskLists_[J-1] == NULL){
              taskLists_[J-1]=new std::list<FBTask>();
            }
            taskLists_[J-1]->push_back(curTask);
            localTaskCount_++;
            taskit = --taskLists_[J-1]->end();

            //compute cost
            if(1){
              char* dataPtr = msg->GetLocalPtr();
              SuperNode2<T> dist_src_snode(dataPtr,msg->Size(),msg->meta.GIndex);
              dist_src_snode.InitIdxToBlk();
              SnodeUpdate curUpdate;
              taskit->rank = 0.0;
              while(dist_src_snode.FindNextUpdate(curUpdate,Xsuper_,SupMembership_,false)){
                //skip if this update is "lower"
                if(curUpdate.tgt_snode_id<curTask.tgt_snode_id){continue;}
                Int iUpdater = this->Mapping_->Map(curUpdate.tgt_snode_id-1,curTask.src_snode_id-1);
                if(iUpdater == iam){
                  Int tgt_snode_width = Xsuper_[curUpdate.tgt_snode_id] - Xsuper_[curUpdate.tgt_snode_id-1];
                  taskit->rank += dist_src_snode.NRowsBelowBlock(0)*pow(tgt_snode_width,2.0);
                }
              }
              taskit->rank*=-1;
            }

          }
#else
          taskit = find_task(msg->meta.src,msg->meta.tgt,UPDATE);
#endif
        } 

        taskit->remote_deps--;
        taskit->data.push_back(msg);


        //if this is a factor task, then we should delete the aggregate
        if(!msg->IsLocal() && taskit->type==AGGREGATE){
          msg->DeallocRemote();
          //ask for remote delete
          //          remote_delete(msg->GetRemotePtr());
        }

        if(taskit->type==AGGREGATE){
          assert(taskit->remote_deps==0 && taskit->local_deps==0);
        }

        if(!is_static){
          if(taskit->remote_deps==0 && taskit->local_deps==0){

            //        //compute cost 
            //        if(taskit->type==FACTOR){
            //          //TODO this is the factorization cost only. There is also the aggregation cost to take into account
            //          SuperNode2<T> & src_snode = *snodeLocal(taskit->src_snode_id);
            //          taskit->rank = 2.0*pow((double)src_snode.Size(),3.0)/3.0 + src_snode.Size()*(src_snode.Size()+1.0)*(src_snode.NRowsBelowBlock(0)-src_snode.Size())/2.0;
            //        }
            //        else if(taskit->type == UPDATE){
            //          //TODO loop through the source snode and compute cost for ONE target or ALL targets (if source is remote)
            //          taskit->rank = 0.0;
            //        }
            //compute cost 
            //taskit->update_rank();
            scheduler_->push(taskit);    
          }
        }
      }


    }while(msg!=NULL);
  } 
  TIMER_STOP(CHECK_MESSAGE);
}


#endif
