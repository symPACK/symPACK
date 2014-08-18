#ifndef _SUPERNODAL_MATRIX_IMPL_FB_HPP_
#define _SUPERNODAL_MATRIX_IMPL_FB_HPP_



template <typename T> void SupernodalMatrix<T>::FanBoth(){

  TIMER_START(FACTORIZATION_FB);








  Real timeSta, timeEnd;


  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();

  IntNumVec UpdatesToDo;
  IntNumVec LastUpdate;
  IntNumVec AggregatesToRecv;
  //gdb_lock(0);
  UpdatesToDo = UpdateCount_;
  AggregatesToRecv = UpdateCount_;
  //FBGetUpdateCount(AggregatesToRecv,LastUpdate);
  
  double timesta2 = get_time(); 
  FBGetUpdateCount(UpdatesToDo,LastUpdate);
  double timeend2 = get_time(); 
  logfileptr->OFS()<<"Update count time: "<<timeend2-timesta2<<endl;

#ifdef _DEBUG_UPDATES_
  //logfileptr->OFS()<<"LastUpdate: "<<LastUpdate<<endl;
  logfileptr->OFS()<<"UpdatesToDo: "<<UpdatesToDo<<endl;
  logfileptr->OFS()<<"AggregatesToRecv: "<<AggregatesToRecv<<endl;
#endif

  IntNumVec AggregatesDone(Xsuper_.m());
  SetValue(AggregatesDone,0);


#ifndef _USE_TAU_
  FBCommList<T> MsgToSend; 
#else
  FBCommList MsgToSend; 
#endif
  AsyncComms outgoingSend;

  //CommList UpdatesToSend; 


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
  TempUpdateBuffers<T> tmpBufs(Size(),maxwidth);




  std::vector< SuperNode<T> * > aggVectors(Xsuper_.m()-1,NULL);





  //
  //  //Initialize the list of local updates as well as the factors to receive
  //  for(Int i = LocalSupernodes_.size()-1; i>=0; --i){
  //    SuperNode<T> & cur_snode = *LocalSupernodes_[i];
  //
  //    SnodeUpdate curUpdate;
  //    while(cur_snode.FindNextUpdate(curUpdate,Xsuper_,SupMembership_)){ 
  //      Int iTarget = Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
  //      if(iTarget == iam){
  //        Int iLocalJ = (curUpdate.tgt_snode_id-1) / np +1 ;
  //        LocalUpdates[iLocalJ-1].push((SnodeUpdate)curUpdate);
  //
  //#ifdef _DEBUG_
  //        logfileptr->OFS()<<"FUTURE LOCAL Supernode "<<curUpdate.tgt_snode_id<<" is going to be updated by Supernode "<<cur_snode.Id()<<" from row "<<curUpdate.src_first_row<<" "<<curUpdate.blkidx<<std::endl;
  //#endif
  //      }
  //    }
  //    FactorsToRecv[i] = UpdatesToDo(cur_snode.Id()-1) % np;//UpdatesToDo(cur_snode.Id()-1) - LocalUpdates[i].size();
  //  }
  //
  //
  //#ifdef _DEBUG_
  //  logfileptr->OFS()<<"FactorsToRecv: "<<FactorsToRecv<<std::endl;
  //#endif




  timeSta =  get_time( );

  FBTasks LocalTasks;

  std::vector<SnodeUpdate> LocalUpdates;

  std::vector< Int > ProcessedSnodes;
  for(Int I = 1; I<Xsuper_.m(); ++I){
    Int iOwner = Mapping_->Map(I-1,I-1);
    Int J = FBUpdate(I); 
    if(iam==iOwner){
      SnodeUpdateFB curUpdate;
      curUpdate.src_snode_id = I;
      curUpdate.tgt_snode_id = I;
      LocalTasks.push(curUpdate);
      
      if(J!=-1){
      SnodeUpdateFB curUpdate;
      curUpdate.src_snode_id = I;
      curUpdate.tgt_snode_id = J;
      LocalTasks.push(curUpdate);
      }

      ProcessedSnodes.push_back(I);
      //      if(J!=-1){
      //I may have to do local updates
      ProcessedSnodes.push_back(J);
      //      }
    }
    else{

      if(J!=-1){
      SnodeUpdateFB curUpdate;
      curUpdate.src_snode_id = I;
      curUpdate.tgt_snode_id = J;
      LocalTasks.push(curUpdate);
      }
      //TODO replace this by a PLINDX ? 
      //am I updating the snode ?
      if(J!=-1){
        ProcessedSnodes.push_back(I);
        ProcessedSnodes.push_back(J);
      }
    }


    if(J!=-1){
      SnodeUpdate curUpdate;
      curUpdate.src_snode_id = I;
      curUpdate.tgt_snode_id = J;
      LocalUpdates.push_back(curUpdate);
    }



  }


  Int I =1;
  Int iLocalI = 1;
  while(!LocalTasks.empty() || !MsgToSend.empty() || !outgoingSend.empty()){

    //Check for completion of outgoing communication
    AdvanceOutgoing(outgoingSend);
    //We have to have a list of remote supernodes for which we are computing updates

    if(!LocalTasks.empty()){
      //make a copy because we need to pop it
      SnodeUpdateFB curTask = LocalTasks.top();
      LocalTasks.pop();

      const SnodeUpdateFB * nextTask = (!LocalTasks.empty())?&LocalTasks.top():NULL;

#ifdef _DEBUG_DELAY_
  {
    FBTasks tmp = LocalTasks;
    logfileptr->OFS()<<"Task Queue : ";
    while( tmp.size()>0){
      //Pull the highest priority message
      const SnodeUpdateFB & comm = tmp.top();
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
      I = src_snode_id;
      //If it is a factorization
      if(src_snode_id == tgt_snode_id){
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
        src_blocks.resize(0);

        Int nz_cnt;
        Int max_bytes;
        while(AggregatesToRecv(I-1)>0){
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
            Int ncols = src_snode.Size();
            nz_cnt = nrows * ncols;

            max_bytes += nrows*sizeof(NZBlockDesc);
            max_bytes += nz_cnt*sizeof(T); 


            max_bytes += sizeof(Int); 

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
            gdb_lock();
            prev_size = src_blocks.size();
            doabort = true;
            //receive anyway
            src_blocks.resize(bytes_received);
          }
#endif


//          logfileptr->OFS()<<"RECV2 Bfore"<<endl;
#ifdef PROBE_FIRST
          MPI_Recv(&src_blocks[0],src_blocks.size(),MPI_BYTE,recv_status.MPI_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
#else
          //receive the index array
          MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,MPI_ANY_SOURCE,I,CommEnv_->MPI_GetComm(),&recv_status);
#endif
//          logfileptr->OFS()<<"RECV2 after"<<endl;

          SuperNode<T> dist_src_snode;
          size_t read_bytes = Deserialize(&src_blocks[0],dist_src_snode);

          //Deserialize the number of aggregates
          Int * aggregatesCnt = (Int *)(&src_blocks[0]+read_bytes);


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
          dist_src_snode.FindNextUpdate(curUpdate,Xsuper_,SupMembership_,false);
          FBAggregateSuperNode(dist_src_snode,src_snode,curUpdate.blkidx, curUpdate.src_first_row);
          AggregatesToRecv(curUpdate.tgt_snode_id-1) -= *aggregatesCnt;

        }
        //clear the buffer
        { vector<char>().swap(src_blocks);  }




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
            if(!is_factor_sent[iTarget] /*&& !is_skipped[iTarget]*/ ){


                Int tag = FACT_TAG(curUpdate.src_snode_id,curUpdate.tgt_snode_id);
#ifndef _USE_TAU_
                MsgToSend.push(FBDelayedComm<T>(FACTOR,&src_snode,curUpdate.src_snode_id,curUpdate.tgt_snode_id,curUpdate.blkidx,curUpdate.src_first_row,iTarget,tag));
#else
                MsgToSend.push(FBDelayedComm(FACTOR,(void*)&src_snode,curUpdate.src_snode_id,curUpdate.tgt_snode_id,curUpdate.blkidx,curUpdate.src_first_row,iTarget,tag));
#endif
//                MsgToSend.emplace(FACTOR,&src_snode,curUpdate.src_snode_id,curUpdate.tgt_snode_id,curUpdate.blkidx,curUpdate.src_first_row,iTarget,tag,0);
#ifdef _DEBUG_DELAY_
                logfileptr->OFS()<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<curUpdate.tgt_snode_id<<" from row "<<curUpdate.src_first_row<<endl;
                cout<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<curUpdate.tgt_snode_id<<" from row "<<curUpdate.src_first_row<<endl;
#endif  
//                is_skipped[iTarget] = true;
                is_factor_sent[iTarget] = true;
            }

          }
        }
        TIMER_STOP(FIND_UPDATED_ANCESTORS);



        ++iLocalI; 
      }
      else{
          //Receive the factor
          src_blocks.resize(0);
          SuperNode<T> * cur_src_snode = FBRecvFactor(src_snode_id,tgt_snode_id,src_blocks);


          Int iSrcOwner = Mapping_->Map(src_snode_id-1,src_snode_id-1);


          //Update everything src_snode_id own with that factor
          //Update the ancestors
          SnodeUpdate curUpdate;

          TIMER_START(UPDATE_ANCESTORS);

          while(cur_src_snode->FindNextUpdate(curUpdate,Xsuper_,SupMembership_,iam==iSrcOwner)){

#ifdef _DEBUG_PROGRESS_
        logfileptr->OFS()<<"Implicit Task: {"<<curUpdate.src_snode_id<<" -> "<<curUpdate.tgt_snode_id<<"}"<<std::endl;
logfileptr->OFS()<<"Processing update from Supernode "<<curUpdate.src_snode_id<<" to Supernode "<<curUpdate.tgt_snode_id<<endl;
#endif
            Int iUpdater = Mapping_->Map(curUpdate.tgt_snode_id-1,cur_src_snode->Id()-1);
            if(iUpdater == iam){
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
                  aggVectors[curUpdate.tgt_snode_id-1] = new SuperNode<T>(curUpdate.tgt_snode_id, Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id]-1,  xlindx_, lindx_);
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


              tgt_aggreg->Update(*cur_src_snode,curUpdate,tmpBufs);

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
#endif
                    //MsgToSend.emplace(AGGREGATE,tgt_aggreg,curUpdate.src_snode_id,curUpdate.tgt_snode_id,0,pivot_desc.GIndex,iTarget,tag,AggregatesDone[curUpdate.tgt_snode_id-1]);
#ifdef _DEBUG_DELAY_
                    //                            logfileptr->OFS()<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<curUpdate.tgt_snode_id<<" from row "<<curUpdate.src_first_row<<endl;
                    //                            cout<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<curUpdate.tgt_snode_id<<" from row "<<curUpdate.src_first_row<<endl;
#endif


                }
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
    }

    //process some of the delayed send
    SendDelayedMessagesUp(MsgToSend,outgoingSend,LocalTasks,&AggregatesDone[0]);

  }



  MPI_Barrier(CommEnv_->MPI_GetComm());

  TIMER_STOP(FACTORIZATION_FO);
}

template<typename T> Int SupernodalMatrix<T>::FBUpdate(Int I){
  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();
  //Check if I have anything to update with that supernode
  //look at lindx_
  Int fi = xlindx_(I-1);
  Int li = xlindx_(I)-1;


  Int firstUpdate = -1; 
  Int J = -1; 
  for(Int idx = fi; idx<=li;++idx){
    Int row = lindx_[idx-1];
    J = SupMembership_[row-1];
    Int iUpdater = Mapping_->Map(J-1,I-1);
    if(iUpdater == iam && J>I){
      firstUpdate = J;
      break;
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
  Int iam = CommEnv_->MPI_Rank();
  Int np  = CommEnv_->MPI_Size();
  SuperNode<T> * cur_src_snode = NULL;
  Int iSrcOwner = Mapping_->Map(src_snode_id-1,src_snode_id-1);

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
    MPI_Probe(iSrcOwner,FACT_TAG(src_snode_id,tgt_snode_id),CommEnv_->MPI_GetComm(),&recv_status);
    MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);

    TIMER_START(RECV_MALLOC);
    src_blocks.resize(bytes_received);
    TIMER_STOP(RECV_MALLOC);

//    logfileptr->OFS()<<"RECV Bfore"<<endl;
    MPI_Recv(&src_blocks[0],bytes_received,MPI_BYTE,recv_status.MPI_SOURCE,recv_status.MPI_TAG,CommEnv_->MPI_GetComm(),&recv_status);
    MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
//    logfileptr->OFS()<<"RECV After"<<endl;




    SuperNode<T> * dist_src_snode = new SuperNode<T>();
    Deserialize(&src_blocks[0],*dist_src_snode);

#ifdef _DEBUG_
    logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode->Id()<<std::endl;
#endif
    cur_src_snode = dist_src_snode;

  }
  return cur_src_snode;
}


template <typename T> inline void SupernodalMatrix<T>::FBAggregateSuperNode(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode, Int &pivot_idx, Int  pivot_fr)
{

  TIMER_START(AGGREGATE_SNODE);

#ifdef SINGLE_BLAS

  TIMER_START(AGGREGATE_SNODE_FIND_INDEX);
  Int first_pivot_idx = -1;
  Int tgt_fc = pivot_fr;
  if(tgt_fc ==I_ZERO ){
#ifdef FAST_INDEX_SEARCH
    Int tgt_fc = tgt_snode.FirstCol();
    first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
    if(first_pivot_idx<0){
      tgt_fc = -first_pivot_idx;
    }
#else
    tgt_fc = tgt_snode.FirstCol();
    //find the pivot idx
    do {first_pivot_idx = src_snode.FindBlockIdx(tgt_fc); tgt_fc++;}
    while(first_pivot_idx<0 && tgt_fc<=tgt_snode.LastCol());
    tgt_fc--;
#endif
  }
  else{
    first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
  }
  NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);

#ifdef FAST_INDEX_SEARCH
  //    TIMER_START(UPDATE_SNODE_FIND_INDEX_LAST2);
  Int tgt_lc = tgt_snode.LastCol();
  Int last_pivot_idx = src_snode.FindBlockIdx(tgt_lc);
  if(last_pivot_idx<0){
    if(last_pivot_idx == -(iSize_+1)){
      last_pivot_idx = src_snode.NZBlockCnt()-1;
    }
    else{
      last_pivot_idx = src_snode.FindBlockIdx(-last_pivot_idx)-1;
    }
    assert(last_pivot_idx>=0);
  }
  NZBlockDesc & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);
  tgt_lc = min(tgt_lc,last_pivot_desc.GIndex + src_snode.NRows(last_pivot_idx)-1);
  //    TIMER_STOP(UPDATE_SNODE_FIND_INDEX_LAST2);
#else
  //    TIMER_START(UPDATE_SNODE_FIND_INDEX_LAST);
  Int tgt_lc = tgt_snode.LastCol();
  Int last_pivot_idx = -1;
  //find the pivot idx
  do {last_pivot_idx = src_snode.FindBlockIdx(tgt_lc); tgt_lc--;}
  while(last_pivot_idx<0 && tgt_lc>=tgt_fc);
  tgt_lc++;
  NZBlockDesc & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);
  //    TIMER_STOP(UPDATE_SNODE_FIND_INDEX_LAST);
#endif

  TIMER_STOP(AGGREGATE_SNODE_FIND_INDEX);

  //determine the first column that will be updated in the target supernode
  Int tgt_local_fc =  tgt_fc - tgt_snode.FirstCol();
  Int tgt_local_lc =  tgt_lc - tgt_snode.FirstCol();

  Int src_nrows = src_snode.NRowsBelowBlock(first_pivot_idx) - (tgt_fc - first_pivot_desc.GIndex);
  Int src_lr = tgt_fc+src_nrows-1;
  src_nrows = src_lr - tgt_fc + 1;

  Int tgt_width = src_snode.NRowsBelowBlock(first_pivot_idx) - (tgt_fc - first_pivot_desc.GIndex) - src_snode.NRowsBelowBlock(last_pivot_idx) + (tgt_lc - last_pivot_desc.GIndex)+1;
  //tmpBuf.Resize(tgt_width,src_nrows);

  T * pivot = &(src_snode.GetNZval(first_pivot_desc.Offset)[(tgt_fc-first_pivot_desc.GIndex)*src_snode.Size()]);

  if(1){

    TIMER_START(AGGREGATE_SNODE_INDEX_MAP);
    Int src_snode_size = src_snode.Size();
    Int tgt_snode_size = tgt_snode.Size();
    IntNumVec src_colindx(tgt_width);
    IntNumVec src_rowindx(src_nrows);
    IntNumVec src_to_tgt_offset(src_nrows);
    IntNumVec src_offset(src_nrows);
    SetValue(src_offset,-1);
    SetValue(src_to_tgt_offset,-1);

    Int colidx = 0;
    Int rowidx = 0;
    Int offset = 0;
    for(Int blkidx = first_pivot_idx; blkidx < src_snode.NZBlockCnt(); ++blkidx){
      NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
      Int cur_src_nrows = src_snode.NRows(blkidx);
      Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
      Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
      cur_src_nrows = cur_src_lr - cur_src_fr +1;

      for(Int row = cur_src_fr; row<= cur_src_lr;++row){
        if(row<=tgt_lc){
          src_colindx[colidx++] = row;
        }
        src_rowindx[rowidx] = row;
        src_offset[rowidx] = offset;
        //cur_block_desc.Offset - (first_pivot_desc.Offset + (tgt_fc - first_pivot_desc.GIndex)*src_snode_size ) + (row - cur_block_desc.GIndex)*src_snode_size;
        offset+=src_snode_size;

        Int tgt_blk_idx = tgt_snode.FindBlockIdx(row);
        NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
        src_to_tgt_offset[rowidx] = cur_tgt_desc.Offset + (row - cur_tgt_desc.GIndex)*tgt_snode_size; 
        rowidx++;
      }
    }
    TIMER_STOP(AGGREGATE_SNODE_INDEX_MAP);

#ifdef _DEBUG_ 
    logfileptr->OFS()<<"src_rowindx :"<<src_rowindx<<std::endl;
    logfileptr->OFS()<<"src_colindx :"<<src_colindx<<std::endl;
    logfileptr->OFS()<<"Index map tgt :"<<src_to_tgt_offset<<std::endl;
    logfileptr->OFS()<<"Index map src :"<<src_offset<<std::endl;
#endif

    TIMER_START(AGGREGATE_SNODE_ASSEMBLY);
    T* tgt = tgt_snode.GetNZval(0);
    for(Int rowidx = 0; rowidx < src_rowindx.m(); ++rowidx){
      Int row = src_rowindx[rowidx];
      for(Int colidx = 0; colidx< src_colindx.m();++colidx){
        Int col = src_colindx[colidx];
        Int tgt_colidx = col - tgt_snode.FirstCol();
        Int src_colidx = col - src_snode.FirstCol();
        tgt[src_to_tgt_offset[rowidx] + tgt_colidx] += pivot[src_offset[rowidx]+src_colidx]; 
      }
    }
    TIMER_STOP(AGGREGATE_SNODE_ASSEMBLY);
    //logfileptr->OFS()<<"After "<<std::endl<<tgt_snode<<std::endl;



  }








#else
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

#endif
  TIMER_STOP(AGGREGATE_SNODE);
}








#ifndef _USE_TAU_
template <typename T> void SupernodalMatrix<T>::SendDelayedMessagesUp(FBCommList<T> & MsgToSend, AsyncComms & OutgoingSend, FBTasks & taskList, Int * AggregatesDone)
#else
template <typename T> void SupernodalMatrix<T>::SendDelayedMessagesUp(FBCommList & MsgToSend, AsyncComms & OutgoingSend, FBTasks & taskList, Int * AggregatesDone)
#endif
{
  if(MsgToSend.empty()) { return;}
  bool is_last = taskList.empty();

  //Index of the last global snode to do
  const SnodeUpdateFB * nextTask = (!taskList.empty())?&taskList.top():NULL;
  FBTaskType nextType = nextTask!=NULL?(nextTask->src_snode_id==nextTask->tgt_snode_id?FACTOR:AGGREGATE):FACTOR;

#ifndef _USE_TAU_
#ifdef _DEBUG_DELAY_
  {
    FBCommList<T> tmp = MsgToSend;
    logfileptr->OFS()<<"Comm "<<"Queue : ";
    while( tmp.size()>0){
      //Pull the highest priority message
      const FBDelayedComm<T> & comm = tmp.top();
      if(comm.type==FACTOR){
        logfileptr->OFS()<<" { F "<<comm.src_data->Id()<<" -> "<<comm.tgt_snode_id<<" }";
      }
      else{
        logfileptr->OFS()<<" { A "<<comm.src_data->Id()<<" -> "<<comm.tgt_snode_id<<" }";
      }
      tmp.pop();
    }
    logfileptr->OFS()<<endl;
  }
#endif
#endif

#ifndef _USE_TAU_
  FBDelayedCommCompare<T> comparator;
#else
  FBDelayedCommCompare comparator;
#endif
//comparator.unitTest();

  while( MsgToSend.size()>0){
    //Pull the highest priority message
#ifndef _USE_TAU_
    const FBDelayedComm<T> & comm = MsgToSend.top();
    const Int & tgt_snode_id = comm.tgt_snode_id;
    const Int & src_nzblk_idx = comm.src_nzblk_idx;
    const Int & src_first_row = comm.src_first_row;
    const FBTaskType & type = comm.type;
    SuperNode<T> * src_data = comm.src_data;
    Int src_snode_id = src_data->Id();
#else
    const FBDelayedComm & comm = MsgToSend.top();
    const Int & tgt_snode_id = comm.tgt_snode_id;
    const Int & src_nzblk_idx = comm.src_nzblk_idx;
    const Int & src_first_row = comm.src_first_row;
    const FBTaskType & type = comm.type;
    SuperNode<T> * src_data = (SuperNode<T> *)comm.src_data;
    Int src_snode_id = comm.src_snode_id;
#endif
    
    bool doSend = true;
    if(nextTask!=NULL){
      //gdb_lock(3);
      doSend = comparator.compare(nextTask->src_snode_id,nextTask->tgt_snode_id,
                                                      nextType,src_snode_id,tgt_snode_id,type);
    }
    //if we still have async send, we can do the send anyway
    doSend = doSend ||  OutgoingSend.size() <= maxIsend_ ;
    //if it is the last thing we have to do, do it anyway
    doSend = doSend || is_last;

#ifdef _DEBUG_PROGRESS_
      if(nextTask!=NULL){
        logfileptr->OFS()<<"Comm "<<(type==FACTOR?"F":"A")<<" {"<<src_snode_id<<" -> "<<tgt_snode_id<<"} vs Task "
                <<(nextType==FACTOR?"F":"A")<<" {"<<nextTask->src_snode_id<<" -> "<<nextTask->tgt_snode_id<<"}"<<std::endl;
      }
#endif



    if(doSend){
      SuperNode<T> & src_snode = *src_data;

#ifdef _DEBUG_DELAY_
      if(type==FACTOR){
        logfileptr->OFS()<<"Picked Comm { F "<<src_snode_id<<" -> "<<tgt_snode_id<<" }"<<endl;
      }
      else{
        logfileptr->OFS()<<"Picked Comm { A "<<src_snode_id<<" -> "<<tgt_snode_id<<" }"<<endl;
      }
#endif
      //this can be sent now
      Int iTarget, tag;
      if(type==FACTOR){
        iTarget = FACT_TARGET(Mapping_,src_snode_id,tgt_snode_id);
        tag = FACT_TAG(src_snode_id,tgt_snode_id);
      }
      else{
        iTarget = AGG_TARGET(Mapping_,src_snode_id,tgt_snode_id);
        tag = AGG_TAG(src_snode_id,tgt_snode_id);
      }
//      iTarget = comm.target;


      if(iTarget != iam){

#ifdef _DEBUG_DELAY_
      if(type==FACTOR){
        logfileptr->OFS()<<"P"<<iam<<" is sending Factor from Supernode "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<endl;
        cout<<"P"<<iam<<" is sending Factor from Supernode "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<endl;
      }
      else{
        logfileptr->OFS()<<"P"<<iam<<" is sending Aggregate from Supernode "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<endl;
        cout<<"P"<<iam<<" is sending Aggregate from Supernode "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<endl;
      }
#endif
        NZBlockDesc & pivot_desc = src_snode.GetNZBlockDesc(src_nzblk_idx);
        //Create a new Icomm buffer, serialize the contribution
        // in it and add it to the outgoing comm list
        Icomm * send_buffer = new Icomm();
        //if this is an aggregate add the count to the end
        if(type==AGGREGATE){
#ifdef _DEBUG_
          assert(AggregatesDone != NULL);
#endif
          Serialize(*send_buffer,src_snode,src_nzblk_idx,src_first_row,sizeof(Int));
          *send_buffer<<AggregatesDone[tgt_snode_id-1];
          //set it to 0 agains to indicate that it has been freed
          AggregatesDone[tgt_snode_id-1] = 0;
        }
        else{
          Serialize(*send_buffer,src_snode,src_nzblk_idx,src_first_row);
        }
        AddOutgoingComm(OutgoingSend,send_buffer);

        if( OutgoingSend.size() > maxIsend_){
          TIMER_START(SEND_MPI);
          MPI_Send(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tag,CommEnv_->MPI_GetComm());
          TIMER_STOP(SEND_MPI);
          OutgoingSend.pop_back();
        }
        else{
          TIMER_START(SEND_MPI);
          MPI_Isend(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tag,CommEnv_->MPI_GetComm(),&OutgoingSend.back()->Request);
          TIMER_STOP(SEND_MPI);
        }
#ifdef _DEBUG_DELAY_
      if(type==FACTOR){
        logfileptr->OFS()<<"P"<<iam<<" has sent Factor from Supernode "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<endl;
        cout<<"P"<<iam<<" has sent Factor from Supernode "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<endl;
      }
      else{
        logfileptr->OFS()<<"P"<<iam<<" has sent Aggregate from Supernode "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<endl;
        cout<<"P"<<iam<<" has sent Aggregate from Supernode "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<endl;
      }
#endif


      }
      //remove from the list
      //This will also free the aggregate

        if(comm.type==AGGREGATE){
          delete src_data;
        }
      MsgToSend.pop();
    }
    else{
      break;
    }
  }
}














#endif
