#ifndef _SUPERNODAL_MATRIX_IMPL_FB_HPP_
#define _SUPERNODAL_MATRIX_IMPL_FB_HPP_

//void MyProbe(std::function<bool(Int)> tag_validator){
//    bool flag = false;
//    while(!flag){
//          MPI_Status recv_status;
//          int bytes_received = 0;
//
//          MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,CommEnv_->MPI_GetComm(),&recv_status);
//          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
//
//          Int tag = recv_status.MPI_TAG;
//        if
//    }
// }


template <typename T> void SupernodalMatrix<T>::FBFactorizationTask(SnodeUpdateFB & curTask, Int iLocalI, IntNumVec & AggregatesDone, IntNumVec & AggregatesToRecv, std::vector<char> & src_blocks)
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
        //Wait for all the aggregates BUT receive from any
        while(AggregatesToRecv(tgt_snode_id-1)>0){
#ifdef _DEBUG_
          logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left"<<endl;
#endif




            //does not work because we are receiving factors as well. Should work if I create a duplicate comm for aggregates
//          TIMER_START(RECV_MPI);
//          MPI_Status recv_status;
//          int bytes_received = 0;
//
//          MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,CommEnv_->MPI_GetComm(),&recv_status);
//          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
//
//          Int tag = recv_status.MPI_TAG;
//          
//          TIMER_START(RECV_MALLOC);
//          src_blocks.resize(bytes_received);
//          TIMER_STOP(RECV_MALLOC);
//          MPI_Recv(&src_blocks[0],src_blocks.size(),MPI_BYTE,recv_status.MPI_SOURCE,tag,CommEnv_->MPI_GetComm(),&recv_status);
//          TIMER_STOP(RECV_MPI);
//
//          SuperNode<T> dist_src_snode;
//          size_t read_bytes = Deserialize(&src_blocks[0],dist_src_snode);
//          //Deserialize the number of aggregates
//          Int * aggregatesCnt = (Int *)(&src_blocks[0]+read_bytes);
//
//
//#ifdef _DEBUG_
//          logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode.Id()<<std::endl;
//#endif
//          logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode.Id()<<std::endl;
//
//          Int local_tgt_id = (dist_src_snode.Id()-1)/np +1;
//          SuperNode<T> & tgt_snode = *LocalSupernodes_[local_tgt_id -1];
//          logfileptr->OFS()<<"TGT Supernode "<<tgt_snode.Id()<<std::endl;
//
//assert(dist_src_snode.Id()==tgt_snode.Id());
//
//          tgt_snode.Aggregate(dist_src_snode);
//          AggregatesToRecv(tgt_snode.Id()-1) -= *aggregatesCnt;
















          TIMER_START(RECV_MPI);
          MPI_Status recv_status;
          int bytes_received = 0;

          Int tag = AGG_TAG(src_snode_id,tgt_snode_id);
          MPI_Probe(MPI_ANY_SOURCE,tag,CommEnv_->MPI_GetComm(),&recv_status);
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
          TIMER_START(RECV_MALLOC);
          src_blocks.resize(bytes_received);
          TIMER_STOP(RECV_MALLOC);
          MPI_Recv(&src_blocks[0],src_blocks.size(),MPI_BYTE,recv_status.MPI_SOURCE,tag,CommEnv_->MPI_GetComm(),&recv_status);
          TIMER_STOP(RECV_MPI);

          SuperNode<T> dist_src_snode;
          size_t read_bytes = Deserialize(&src_blocks[0],dist_src_snode);
          //Deserialize the number of aggregates
          Int * aggregatesCnt = (Int *)(&src_blocks[0]+read_bytes);


#ifdef _DEBUG_
          logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode.Id()<<std::endl;
#endif

          src_snode.Aggregate(dist_src_snode);
          AggregatesToRecv(src_snode.Id()-1) -= *aggregatesCnt;

        }
        //clear the buffer
        //{ vector<char>().swap(src_blocks);  }


        TIMER_STOP(RECV_AGGREGATES);


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
#ifndef _USE_TAU_
                MsgToSend.push(FBDelayedComm<T>(FACTOR,&src_snode,curUpdate.src_snode_id,curUpdate.tgt_snode_id,curUpdate.blkidx,curUpdate.src_first_row,iTarget,tag));
#else
                MsgToSend.push(FBDelayedComm(FACTOR,(void*)&src_snode,curUpdate.src_snode_id,curUpdate.tgt_snode_id,curUpdate.blkidx,curUpdate.src_first_row,iTarget,tag));
#endif
#ifdef _DEBUG_DELAY_
                logfileptr->OFS()<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<curUpdate.tgt_snode_id<<" from row "<<curUpdate.src_first_row<<endl;
                cout<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<curUpdate.tgt_snode_id<<" from row "<<curUpdate.src_first_row<<endl;
#endif  
                is_factor_sent[iTarget] = true;

              //process some of the delayed send
              SendDelayedMessagesUp(MsgToSend,outgoingSend,LocalTasks,&AggregatesDone[0]);
            }

          }
        }
        TIMER_STOP(FIND_UPDATED_ANCESTORS);





}



template <typename T> void SupernodalMatrix<T>::FBUpdateTask(SnodeUpdateFB & curTask, IntNumVec & UpdatesToDo, IntNumVec & AggregatesDone, IntNumVec & AggregatesToRecv,std::vector< SuperNode<T> * > & aggVectors, std::vector<char> & src_blocks)
{

      Int src_snode_id = curTask.src_snode_id;
      Int tgt_snode_id = curTask.tgt_snode_id;

          //Receive the factor
          src_blocks.resize(0);
          src_snode_id = abs(src_snode_id);
          SuperNode<T> * cur_src_snode = FBRecvFactor(src_snode_id,tgt_snode_id,src_blocks);

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
          if(curTask.src_snode_id<0){
//gdb_lock();

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
              //if(cur_src_snode->Id()==2){gdb_lock();}
              SuperNode<T> * tgt_aggreg;

              Int iTarget = Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
              if(iTarget == iam){
                //the aggregate vector is directly the target snode
                Int iLocalJ = (curUpdate.tgt_snode_id-1) / np +1 ;
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
                if(iTarget == iam){
                  AggregatesToRecv[curUpdate.tgt_snode_id-1]-=AggregatesDone[curUpdate.tgt_snode_id-1];
#ifdef _DEBUG_
                  logfileptr->OFS()<<AggregatesToRecv(curUpdate.tgt_snode_id-1)<<" aggregates left for Supernode "<<curUpdate.tgt_snode_id<<endl;
#endif
                }
                else{
#ifdef _DEBUG_
                  logfileptr->OFS()<<"Remote Supernode "<<curUpdate.tgt_snode_id<<" is updated by Supernode "<<src_snode_id<<std::endl;
#endif
                  //Push the comm in the comm list
                    Int tag = AGG_TAG(curUpdate.src_snode_id,curUpdate.tgt_snode_id);
                    NZBlockDesc & pivot_desc = tgt_aggreg->GetNZBlockDesc(0);
#ifndef _USE_TAU_
                    MsgToSend.push(FBDelayedComm<T>(AGGREGATE,tgt_aggreg,curUpdate.src_snode_id,curUpdate.tgt_snode_id,0,pivot_desc.GIndex,iTarget,tag,AggregatesDone[curUpdate.tgt_snode_id-1]));
#else
                    MsgToSend.push(FBDelayedComm(AGGREGATE,(void*)tgt_aggreg,curUpdate.src_snode_id,curUpdate.tgt_snode_id,0,pivot_desc.GIndex,iTarget,tag,AggregatesDone[curUpdate.tgt_snode_id-1]));
//                    MsgToSend.push(FBDelayedComm(AGGREGATE,(void*)tgt_aggreg,curUpdate.tgt_snode_id,curUpdate.src_snode_id,0,pivot_desc.GIndex,iTarget,tag,AggregatesDone[curUpdate.tgt_snode_id-1]));
//                    MsgToSend.push(FBDelayedComm(AGGREGATE,(void*)tgt_aggreg,curUpdate.tgt_snode_id,curUpdate.tgt_snode_id,0,pivot_desc.GIndex,iTarget,tag,AggregatesDone[curUpdate.tgt_snode_id-1]));
#endif
                    //MsgToSend.emplace(AGGREGATE,tgt_aggreg,curUpdate.src_snode_id,curUpdate.tgt_snode_id,0,pivot_desc.GIndex,iTarget,tag,AggregatesDone[curUpdate.tgt_snode_id-1]);
#ifdef _DEBUG_DELAY_
                    //                            logfileptr->OFS()<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<curUpdate.tgt_snode_id<<" from row "<<curUpdate.src_first_row<<endl;
                    //                            cout<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<curUpdate.tgt_snode_id<<" from row "<<curUpdate.src_first_row<<endl;
#endif

                  //process some of the delayed send
                  SendDelayedMessagesUp(MsgToSend,outgoingSend,LocalTasks,&AggregatesDone[0]);

                }
              }

//if local update, push a new task in the queue and stop the while loop
if(iam==iSrcOwner ){
break;

//  bool next_local = false;
//
//  while(cur_src_snode->FindNextUpdate(curUpdate,Xsuper_,SupMembership_,iam==iSrcOwner)){
//    Int iUpdater = Mapping_->Map(curUpdate.tgt_snode_id-1,cur_src_snode->Id()-1);
//    if(iUpdater == iam){
//      next_local = true;
//      break;
//    }
//  }
//
//  if(next_local){
////    gdb_lock();
//    SnodeUpdateFB newTask;
//    newTask.src_snode_id = -curUpdate.src_snode_id;
//    newTask.tgt_snode_id = curUpdate.tgt_snode_id;
//    LocalTasks.push(newTask);
//    break;
//  }
}


            }



          }
          TIMER_STOP(UPDATE_ANCESTORS);

          if(iSrcOwner!=iam){
            delete cur_src_snode;
            //clear the buffer
            //{ vector<char>().swap(src_blocks);  }

            { vector<char>().swap(src_blocks);  }
          }

}



template <typename T> void SupernodalMatrix<T>::FanBoth()
{
  TIMER_START(FACTORIZATION_FB);

  Real timeSta, timeEnd;


  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();

  IntNumVec UpdatesToDo;
  IntNumVec LastUpdate;
  IntNumVec AggregatesToRecv;

  UpdatesToDo = UpdateCount_;
  AggregatesToRecv = UpdateCount_;
  
  double timesta2 = get_time(); 
  FBGetUpdateCount(UpdatesToDo,LastUpdate);
  double timeend2 = get_time(); 
  logfileptr->OFS()<<"Update count time: "<<timeend2-timesta2<<endl;

#ifdef _DEBUG_UPDATES_
  logfileptr->OFS()<<"UpdatesToDo: "<<UpdatesToDo<<endl;
  logfileptr->OFS()<<"AggregatesToRecv: "<<AggregatesToRecv<<endl;
#endif

  IntNumVec AggregatesDone(Xsuper_.m());
  SetValue(AggregatesDone,0);


  //Array for Irecv of factors
  std::vector<AsyncComms> incomingRecvArr(LocalSupernodes_.size());
  incomingRecvCnt_ = 0;
  IntNumVec FactorsToRecv(LocalSupernodes_.size());

  std::vector<T> src_nzval;
  std::vector<char> src_blocks;

  //std::vector<std::queue<SnodeUpdate> > LocalUpdates(LocalSupernodes_.size());

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

  for(Int I = 1; I<Xsuper_.m(); ++I){
    Int iOwner = Mapping_->Map(I-1,I-1);
    if(iam==iOwner){
      SnodeUpdateFB curUpdate;
      curUpdate.src_snode_id = I;
      curUpdate.tgt_snode_id = I;
      LocalTasks.push(curUpdate);
      Int J = -1;
     
      bool is_first = true; 
      while( (J=FBUpdate(I,J))!= -1 ){
        SnodeUpdateFB curUpdate;
        //set it to be negative as it is a local update
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
        curUpdate.src_snode_id = I;
        curUpdate.tgt_snode_id = J;
        LocalTasks.push(curUpdate);
      }
    }
  }


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

#ifdef _DEADLOCK_
      const SnodeUpdateFB * nextTask = (!LocalTasks.empty())?&LocalTasks.front():NULL;
#else
      const SnodeUpdateFB * nextTask = (!LocalTasks.empty())?&LocalTasks.top():NULL;
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

      Int src_snode_id = curTask.src_snode_id;
      Int tgt_snode_id = curTask.tgt_snode_id;
      //If it is a factorization
      if(src_snode_id == tgt_snode_id){
        FBFactorizationTask(curTask,iLocalI,AggregatesDone,AggregatesToRecv,src_blocks);
        ++iLocalI; 
      }
      else{
        FBUpdateTask(curTask, UpdatesToDo, AggregatesDone, AggregatesToRecv,aggVectors, src_blocks);
      }
      LocalTasks.pop();
    }

    //process some of the delayed send
    SendDelayedMessagesUp(MsgToSend,outgoingSend,LocalTasks,&AggregatesDone[0]);

  }



  MPI_Barrier(CommEnv_->MPI_GetComm());

  tmpBufs.Clear();
  TIMER_STOP(FACTORIZATION_FB);
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
//if(prevJ!=-1){ gdb_lock();}
          firstUpdate = J;
          break;
        }
      }
      else{
        firstUpdate = J;
        break;
      }
    }
  }

  return firstUpdate;
}


template <typename T> void SupernodalMatrix<T>::FBGetUpdateCount(IntNumVec & sc,IntNumVec & lu){
  sc.Resize(Xsuper_.m());
  SetValue(sc,I_ZERO);

//  lu.Resize(Xsuper_.m());
//  SetValue(lu,I_ZERO);

  IntNumVec marker(Xsuper_.m());
  SetValue(marker,I_ZERO);

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



        Int iUpdater = Mapping_->Map(supno-1,s-1);

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


template<typename T> SuperNode<T> * SupernodalMatrix<T>::FBRecvFactor(Int src_snode_id,Int tgt_snode_id, std::vector<char> & src_blocks){
  TIMER_START(RECV_FACTORS);
  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();
  SuperNode<T> * cur_src_snode = NULL;
  Int iSrcOwner = Mapping_->Map(src_snode_id-1,src_snode_id-1);
  

  //This is a local update, iSrcOwner matters
  if(iSrcOwner==iam){
    Int iLocalI = (src_snode_id-1) / np +1 ;
    cur_src_snode = LocalSupernodes_[iLocalI -1];
  }
  else{
    Int nz_cnt;

    TIMER_START(RECV_MPI);
    MPI_Status recv_status;
    int bytes_received = 0;

    //gdb_lock();
    //This is a remote update, iSrcOwner doesn't matter
    Int tag = FACT_TAG(src_snode_id,tgt_snode_id);
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


    SuperNode<T> * dist_src_snode = new SuperNode<T>();
    Deserialize(&src_blocks[0],*dist_src_snode);
#if not(defined(_LINEAR_SEARCH_FCLC_) || defined(_LAZY_INIT))
    dist_src_snode->InitIdxToBlk();
#endif

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















#endif
