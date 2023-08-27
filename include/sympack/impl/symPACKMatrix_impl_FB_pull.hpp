#ifndef _SYMPACK_MATRIX_IMPL_FB_PULL_HPP_
#define _SYMPACK_MATRIX_IMPL_FB_PULL_HPP_

#include "sympack/symPACKMatrix.hpp"

//#define _NO_COMP_
//#define SP_THREADS




template <typename T> inline void symPACKMatrix<T>::dfs_traversal(std::vector<std::list<Int> > & tree,int node,std::list<Int> & frontier) {
  for (std::list<Int>::iterator it = tree[node].begin(); it!=tree[node].end(); it++) {
    Int I = *it;

    Int iOwner = Mapping_->Map(I-1,I-1);
    if (iOwner == this->iam) {
      frontier.push_back(I);
    }
    else {
      if (I<tree.size()) {
        dfs_traversal(tree,I,frontier);
      }
    }
  }
}

template <typename T> inline void symPACKMatrix<T>::FanBoth() 
{
  scope_timer(a,FACTORIZATION_FB);

  double timeSta, timeEnd;

  SYMPACK_TIMER_START(FB_INIT);
  remDealloc = new upcxx::dist_object<int>(0);

  std::vector<Int> UpdatesToDo = UpdatesToDo_;

  //tmp buffer space
  std::vector<T> src_nzval;
  std::vector<char> src_blocks;


  Int maxheight = 0;
  Int maxwidth = 0;
  for (Int i = 1; i<this->Xsuper_.size(); ++i) {
    Int width =this->Xsuper_[i] - this->Xsuper_[i-1];
    if (width>=maxwidth) {
      maxwidth = width;
    }
    if (UpdateHeight_[i-1]>=maxheight) {
      maxheight = UpdateHeight_[i-1];
    }
  }

  //maxwidth for indefinite matrices
  tmpBufs.Resize(maxwidth + maxheight,maxwidth);
  std::vector< SuperNode<T>* > aggVectors(this->Xsuper_.size()-1,nullptr);

  timeSta =  get_time( );
  //Create a copy of the task graph
  taskGraph graph = taskGraph_;
  {


    std::hash<std::string> hash_fn;

    //This is the important one
#ifdef SP_THREADS
    scheduler_new_->threadInitHandle_ = nullptr;
    scheduler_new_->extraTaskHandle_  = nullptr;
    scheduler_new_->msgHandle = nullptr;

    scheduler_new_->threadInitHandle_ = [this]() {
      std::thread::id tid = std::this_thread::get_id();
      bassert(tid!=this->scheduler_new_->main_tid || Multithreading::NumThread==1);

      std::lock_guard<std::recursive_mutex> lk(scheduler_new_->list_mutex_);
      auto & tmpBuf = tmpBufs_th[tid];
    };

    if (Multithreading::NumThread>1) {
      scheduler_new_->extraTaskHandle_ = [&,this](std::shared_ptr<GenericTask> & pTask)->bool {
        SparseTask & Task = *(SparseTask*)pTask.get();

        Int * meta = reinterpret_cast<Int*>(Task.meta.data());
        Factorization::op_type & type = *reinterpret_cast<Factorization::op_type*>(&meta[3]);

        if (type == Factorization::op_type::FACTOR) {
          return false;
        }
        else {
          int tgt = meta[1];
          SuperNode<T> * tgt_aggreg = nullptr;
          Int iTarget = this->Mapping_->Map(tgt-1,tgt-1);
          if (iTarget == this->iam) {
            tgt_aggreg = snodeLocal(tgt);
          }
          else {
            tgt_aggreg = aggVectors[tgt-1];
            if (tgt_aggreg == nullptr) {
              aggVectors[tgt-1] = CreateSuperNode(this->options_.decomposition);
              tgt_aggreg = aggVectors[tgt-1];
            }
          }

          bool exp = false;
          if (std::atomic_compare_exchange_weak( &tgt_aggreg->in_use, &exp, true )) {
            return false;
          }
          else {
            return true;
          }
        }
      };
    }

#endif

    auto dec_ref = [&graph,this] ( taskGraph::task_iterator taskit, Int loc, Int rem) {
      scope_timer(a,DECREF);
#ifdef SP_THREADS
      std::lock_guard<std::recursive_mutex> lk(scheduler_new_->list_mutex_);
#endif
      taskit->second->local_deps-= loc;
      taskit->second->remote_deps-= rem;

      if (taskit->second->remote_deps==0 && taskit->second->local_deps==0) {

        scheduler_new_->push(taskit->second);
        graph.removeTask(taskit->second->id);
      }

    };

    auto log_task = [&] (taskGraph::task_iterator taskit) {
      SparseTask * tmp = ((SparseTask*)taskit->second.get());
      Int * meta = reinterpret_cast<Int*>(tmp->meta.data());
      Factorization::op_type & type = *reinterpret_cast<Factorization::op_type*>(&meta[3]);
      std::string name;
      switch(type) {
        case Factorization::op_type::FACTOR:
          name="FACTOR";
          break;
        case Factorization::op_type::AGGREGATE:
          name="AGGREGATE";
          break;
        case Factorization::op_type::UPDATE:
          name="UPDATE";
          break;
        default: break; // silence -Wswitch warnings from clang
      }

      logfileptr->OFS()<<" T "<<name<<" "<<meta[0]<<"_"<<meta[1]<<" "<<tmp->local_deps<<" "<<tmp->remote_deps<<std::endl;
    };

    auto log_task_internal = [&] ( const std::shared_ptr<GenericTask> & taskptr) {
      SparseTask * tmp = ((SparseTask*)taskptr.get());
      Int * meta = reinterpret_cast<Int*>(tmp->meta.data());
      Factorization::op_type & type = *reinterpret_cast<Factorization::op_type*>(&meta[3]);
      std::string name;
      switch(type) {
        case Factorization::op_type::FACTOR:
          name="FACTOR";
          break;
        case Factorization::op_type::AGGREGATE:
          name="AGGREGATE";
          break;
        case Factorization::op_type::UPDATE:
          name="UPDATE";
          break;
        default: break; // silence -Wswitch warnings from clang
      }

      std::stringstream sstr;
      sstr<<" Running T "<<name<<" "<<meta[0]<<"_"<<meta[1]<<" "<<tmp->local_deps<<" "<<tmp->remote_deps<<std::endl;
      logfileptr->OFS()<<sstr.str();
    };

#ifdef FANIN_OPTIMIZATION
    if (this->options_.mappingTypeStr ==  "COL2D")
    {
      scope_timer(a,BUILD_SUPETREE);
      ETree SupETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);
      chSupTree_.resize(SupETree.Size()+1);
      for (Int I = 1; I<=SupETree.Size(); I++) {
        Int parent = SupETree.PostParent(I-1);
        chSupTree_[parent].push_back(I);
      }
    }
#endif

    //Finish the initialization of the tasks by creating the lambdas
    {

      scheduler_new_->msgHandle = [hash_fn,&graph,dec_ref,this](std::shared_ptr<IncomingMessage> msg) {
        std::thread::id tid = std::this_thread::get_id();
        bassert(tid==this->scheduler_new_->main_tid);


        scope_timer(a,MSGHANDLE);
        Int iOwner = Mapping_->Map(msg->meta.tgt-1,msg->meta.tgt-1);
        Int iUpdater = Mapping_->Map(msg->meta.tgt-1,msg->meta.src-1);
        //this is an aggregate
        if (iOwner==this->iam && iUpdater!=this->iam) {
#ifdef _DEBUG_PROGRESS_
          logfileptr->OFS()<<"COMM: FETCHED AGGREG MSG("<<msg->meta.src<<","<<msg->meta.tgt<<") from P"<<iUpdater<<std::endl;
#endif

          //insert an aggregation task
          Int I = msg->meta.src;
          Int J = msg->meta.tgt;

          bassert(snodeLocal(J)->Id()==J);

          std::shared_ptr<GenericTask> pTask(new SparseTask);
          SparseTask & Task = *(SparseTask*)pTask.get();
          Task.meta.resize(3*sizeof(Int)+sizeof(Factorization::op_type));
          Int * meta = reinterpret_cast<Int*>(Task.meta.data());
          meta[0] = I;
          meta[1] = J;
          Factorization::op_type & type = *reinterpret_cast<Factorization::op_type*>(&meta[3]);
          type = Factorization::op_type::AGGREGATE;

          Task.remote_deps = 1;
          Task.local_deps = 0;

          Int src = meta[0];
          Int tgt = meta[1];
          Task.execute = [hash_fn,dec_ref,&graph,this,src,tgt,pTask] () {
            scope_timer(a,FB_AGGREGATION_TASK);
            Int iLocalI = snodeLocalIndex(tgt);

#ifdef _DEBUG_PROGRESS_
            logfileptr->OFS()<<"Processing T_AGGREG("<<src<<","<<tgt<<") "<<std::endl;
#endif

            Int src_snode_id = src;
            Int tgt_snode_id = tgt;
            Int I = src_snode_id;
            auto src_snode = LocalSupernodes_[iLocalI -1];
            Int src_first_col = src_snode->FirstCol();
            Int src_last_col = src_snode->LastCol();

            supernode_lock<T> lock(src_snode);
#ifdef _DEBUG_PROGRESS_
            logfileptr->OFS()<<"Processing Supernode "<<I<<std::endl;
#endif


            //Applying aggregates
            //TODO we might have to create a third UPDATE type of task
            //process them one by one
            {
              scope_timer(a,FACT_AGGREGATE);
              assert(pTask->getData().size()==1);
              for (auto msgit = pTask->getData().begin();msgit!=pTask->getData().end();msgit++) {
                auto msgPtr = *msgit;
                assert(msgPtr->IsDone());
                if (!msgPtr->IsLocal()) msgPtr->DeallocRemote(*this->remDealloc);
                char* dataPtr = msgPtr->GetLocalPtr().get();
                auto dist_src_snode = CreateSuperNode(this->options_.decomposition,dataPtr,msgPtr->Size());
                dist_src_snode->InitIdxToBlk();
                src_snode->Aggregate(dist_src_snode);
                delete dist_src_snode;
              }
            }


            char buf[100];
            sprintf(buf,"%d_%d_%d_%d",tgt_snode_id,tgt_snode_id,0,(Int)Factorization::op_type::FACTOR);
            auto id = hash_fn(std::string(buf));

            //this is a particular case : we consider this as a remote dependency
            auto taskit = graph.find_task(id);
            bassert(taskit!=graph.tasks_.end());
            if ( auto ptr = graph.updateTask(id,0,1) ) {
              bassert(ptr!=nullptr);
              scheduler_new_->push(ptr);
            }


          };



          std::stringstream sstr;
          sstr<<meta[0]<<"_"<<meta[1]<<"_"<<0<<"_"<<(Int)type;
          Task.id = hash_fn(sstr.str());
          graph.addTask(pTask);
        }

      };


      for (auto taskit = graph.tasks_.begin(); taskit != graph.tasks_.end(); taskit++) {
        auto & pTask = taskit->second;
        SparseTask & Task = *(SparseTask*)pTask.get();

        Int * meta = reinterpret_cast<Int*>(Task.meta.data());
        Int src = meta[0];
        Int tgt = meta[1];
        Factorization::op_type & type = *reinterpret_cast<Factorization::op_type*>(&meta[3]);

        switch(type) {
          case Factorization::op_type::FACTOR:
            {
              Task.execute = [&,this,src,tgt,pTask] () {
                scope_timer(b,FB_FACTORIZATION_TASK);

                Int iLocalTGT = snodeLocalIndex(tgt);
                Int src_snode_id = src;
                Int tgt_snode_id = tgt;
                Int I = src_snode_id;
                auto src_snode = LocalSupernodes_[iLocalTGT -1];
                Int src_first_col = src_snode->FirstCol();
                Int src_last_col = src_snode->LastCol();

#ifdef _DEBUG_PROGRESS_
                logfileptr->OFS()<<"Factoring Supernode "<<I<<std::endl;
#endif


                SYMPACK_TIMER_START(FACTOR_PANEL);
#ifdef SP_THREADS
                std::thread::id tid = std::this_thread::get_id();
                bassert(tid!=this->scheduler_new_->main_tid || Multithreading::NumThread==1);
                auto & tmpBuf = tmpBufs_th[tid];
                src_snode->Factorize(tmpBuf);
#else
                src_snode->Factorize(tmpBufs);
#endif

                SYMPACK_TIMER_STOP(FACTOR_PANEL);

                //Sending factors and update local tasks
                //Send my factor to my ancestors. 
                std::vector<char> is_factor_sent(this->np);
                SetValue(is_factor_sent,false);

                SnodeUpdate curUpdate;
                SYMPACK_TIMER_START(FIND_UPDATED_ANCESTORS);
                SYMPACK_TIMER_START(FIND_UPDATED_ANCESTORS_FACTORIZATION);
                while (src_snode->FindNextUpdate(curUpdate,this->Xsuper_,this->SupMembership_)) { 
                  Int iTarget = this->Mapping_->Map(curUpdate.tgt_snode_id-1,src_snode->Id()-1);

                  if (iTarget != this->iam) {
                    if (!is_factor_sent[iTarget]) {
                      MsgMetadata meta;

                      //TODO Replace all this by a Serialize function
                      NZBlockDesc & nzblk_desc = src_snode->GetNZBlockDesc(curUpdate.blkidx);
                      Int local_first_row = curUpdate.src_first_row - nzblk_desc.GIndex;
                      Int nzblk_cnt = src_snode->NZBlockCnt() - curUpdate.blkidx;
                      Int nzval_cnt_ = src_snode->Size()*(src_snode->NRowsBelowBlock(curUpdate.blkidx)-local_first_row);
                      T* nzval_ptr = src_snode->GetNZval(nzblk_desc.Offset) + local_first_row*src_snode->Size();

                      upcxx::global_ptr<char> sendPtr = upcxx::to_global_ptr((char*)nzval_ptr);
                      //the size of the message is the number of bytes between sendPtr and the address of nzblk_desc

                      //Send factor 
                      meta.src = curUpdate.src_snode_id;
                      meta.tgt = curUpdate.tgt_snode_id;
                      meta.GIndex = curUpdate.src_first_row;

                      char buf[100];
                      sprintf(buf,"%d_%d_%d_%d",meta.src,meta.tgt,0,(Int)Factorization::op_type::UPDATE);
                      meta.id = hash_fn(std::string(buf));
                      char * last_byte_ptr = (char*)&nzblk_desc + sizeof(NZBlockDesc);
                      size_t msgSize = last_byte_ptr - (char*)nzval_ptr;

                      Int liTarget = this->group_->L2G(iTarget);
                      signal_data(sendPtr, msgSize, liTarget, meta);
                      is_factor_sent[iTarget] = true;
                    }
                  }
                  else {
                    //Update local tasks
                    //find task corresponding to curUpdate

                    char buf[100];
                    sprintf(buf,"%d_%d_%d_%d",curUpdate.src_snode_id,curUpdate.tgt_snode_id,0,(Int)Factorization::op_type::UPDATE);
                    auto id = hash_fn(std::string(buf));
                    auto taskit = graph.find_task(id);
                    bassert(taskit!=graph.tasks_.end());
                    if ( auto ptr = graph.updateTask(id,1,0) ) {
                      bassert(ptr!=nullptr);
                      scheduler_new_->push(ptr);
                    }
                  }
                }
                SYMPACK_TIMER_STOP(FIND_UPDATED_ANCESTORS_FACTORIZATION);
                SYMPACK_TIMER_STOP(FIND_UPDATED_ANCESTORS);
              };




            }
            break;
          case Factorization::op_type::UPDATE:
            {
#ifdef _SEQ_SPECIAL_CASE_
              if (Multithreading::NumThread==1) {
                Task.execute = [&,this,src,tgt,pTask,type] () {
                  scope_timer(a,FB_UPDATE_TASK);

                  Int iLocalTGT = snodeLocalIndex(tgt);
                  Int src_snode_id = src;
                  Int tgt_snode_id = tgt;
                  src_snode_id = abs(src_snode_id);
                  bool is_first_local = src <0;

                  SuperNode<T> * cur_src_snode; 
                  std::shared_ptr<SuperNode<T> > shptr_cur_src_snode = nullptr; 
#ifdef SP_THREADS
                  std::thread::id tid = std::this_thread::get_id();
                  bassert(tid!=this->scheduler_new_->main_tid || Multithreading::NumThread==1);
#endif
                  Int iSrcOwner = this->Mapping_->Map(abs(src_snode_id)-1,abs(src_snode_id)-1);


                  IncomingMessage * structPtr = nullptr;
                  std::shared_ptr<IncomingMessage> msgPtr = nullptr;
                  std::shared_ptr<ChainedMessage<SuperNodeBase<T> > > newMsgPtr = nullptr;
                  //Local or remote factor
                  //we have only one local or one remote incoming aggregate


                  SnodeUpdate curUpdate;
                  bool found = false;

                  if (pTask->getData().size()==0) {
                    cur_src_snode = snodeLocal(src_snode_id);
                  }
                  else {
                    scope_timer(b,FB_UPD_UNPACK_MSG);
                    auto msgit = pTask->getData().begin();
                    msgPtr = *msgit;
                    bassert(msgPtr->IsDone());

                    auto id = pTask->id;

                    bassert(msgPtr->meta.id == id);

                    char* dataPtr = msgPtr->GetLocalPtr().get();
                    {
                      scope_timer(c,FB_UPD_UNPACK_MSG_CREATE);
                      shptr_cur_src_snode.reset(CreateSuperNode(this->options_.decomposition,dataPtr,msgPtr->Size(),msgPtr->meta.GIndex));
                      cur_src_snode = shptr_cur_src_snode.get();
                    }
                    {
                      scope_timer(d,FB_UPD_UNPACK_MSG_INIT_TREE);
                      cur_src_snode->InitIdxToBlk();
                    }

                  }

                  {
                    //Update everything owned locally src_snode_id updates with that factor
                    SYMPACK_TIMER_START(UPDATE_ANCESTORS);
                    SnodeUpdate curUpdate;
                    while (cur_src_snode->FindNextUpdate(curUpdate,this->Xsuper_,this->SupMembership_,this->iam==iSrcOwner)) {

                      //skip if this update is "lower"
                      if (curUpdate.tgt_snode_id<tgt) {
                        continue;
                      }
                      else {
                        Int iUpdater = this->Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.src_snode_id-1);
                        if (iUpdater==this->iam) {
                          SuperNode<T> * tgt_aggreg;
                          Int iTarget = this->Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
                          if (iTarget == this->iam) {
                            //the aggregate std::vector is directly the target snode
                            scope_timer(a,UPD_ANC_Agg_local);
                            tgt_aggreg = snodeLocal(curUpdate.tgt_snode_id);
                            bassert(curUpdate.tgt_snode_id == tgt_aggreg->Id());
                          }
                          else {
                            SYMPACK_TIMER_START(UPD_ANC_Agg_tmp);
                            //Check if src_snode_id already have an aggregate std::vector
                            bool creation_needed = false;
                            if (aggVectors[curUpdate.tgt_snode_id-1]==nullptr) {
                              creation_needed = true;
                            }
                            else if (aggVectors[curUpdate.tgt_snode_id-1]->StorageSize()==0) {
                              creation_needed = true;
                            }

                            if (creation_needed) {
                              SYMPACK_TIMER_START(UPD_ANC_Agg_tmp_creat);
                              //use number of rows below factor as initializer

                              //TODO do a customized version for FANIN as we have all the factors locally
                              // the idea is the following: do a DFS and stop each exploration at the first local descendant of current node
#ifdef FANIN_OPTIMIZION
                              if (this->options_.mappingTypeStr ==  "COL2D")
                              {
                                std::set<Idx> structure;
                                std::list<Int> frontier;
                                {
                                  scope_timer(a,MERGE_STRUCTURE_FANIN);
                                  Idx tgt_fc = this->Xsuper_[curUpdate.tgt_snode_id-1];
                                  dfs_traversal(chSupTree_,curUpdate.tgt_snode_id,frontier);
                                  //process frontier in decreasing order of nodes and merge their structure
                                  for (auto it = frontier.rbegin(); it!=frontier.rend(); it++) {
                                    Int I = *it;
                                    SuperNode<T> * source = snodeLocal(I);
                                    for (Int blkidx=0; blkidx< source->NZBlockCnt();blkidx++) {
                                      NZBlockDesc & nzblk_desc = source->GetNZBlockDesc(blkidx);
                                      Idx fr = nzblk_desc.GIndex;
                                      Idx lr = source->NRows(blkidx) + fr;
                                      for (Idx row = fr; row<lr;row++) {
                                        if (row>=tgt_fc) {
                                          structure.insert(row);
                                        }
                                      }
                                    }
                                  }
                                }

                                if (aggVectors[curUpdate.tgt_snode_id-1]==nullptr) {
                                  aggVectors[curUpdate.tgt_snode_id-1] = CreateSuperNode(this->options_.decomposition);
                                }
                                aggVectors[curUpdate.tgt_snode_id-1]->Init(curUpdate.tgt_snode_id, this->Xsuper_[curUpdate.tgt_snode_id-1], this->Xsuper_[curUpdate.tgt_snode_id-1], this->Xsuper_[curUpdate.tgt_snode_id]-1, this->iSize_,structure,this->options_.panel);
                              } 
                              else
#endif
                              {
                                std::set<Idx> structure;
                                {
                                  scope_timer(a,FETCH_REMOTE_STRUCTURE);
                                  upcxx::global_ptr<char> remoteDesc = std::get<0>(remoteFactors_[curUpdate.tgt_snode_id-1]);
                                  Int block_cnt = std::get<1>(remoteFactors_[curUpdate.tgt_snode_id-1]);


                                  //allocate space to receive block descriptors
                                  char * buffer = (char*)UpcxxAllocator::allocate(sizeof(NZBlockDesc)*block_cnt+ sizeof(SuperNodeDesc));
                                  upcxx::global_ptr<char> remote = remoteDesc;
                                  upcxx::rget(remote, (char*)&buffer[0],block_cnt*sizeof(NZBlockDesc)+sizeof(SuperNodeDesc)).wait();

                                  SuperNodeDesc * pdesc = (SuperNodeDesc*)buffer;
                                  NZBlockDesc * bufferBlocks = (NZBlockDesc*)(pdesc+1);

                                  for (Int i =block_cnt-1;i>=0;i--) {
                                    NZBlockDesc & curdesc = bufferBlocks[i];
                                    size_t end = (i>0)?bufferBlocks[i-1].Offset:pdesc->nzval_cnt_;
                                    Int numRows = (end-curdesc.Offset)/pdesc->iSize_;

                                    for (Idx row = 0; row<numRows;row++) {
                                      structure.insert(curdesc.GIndex+row);
                                    }
                                  }
                                  UpcxxAllocator::deallocate((char*)buffer);
                                }
                                if (aggVectors[curUpdate.tgt_snode_id-1]==nullptr) {
                                  aggVectors[curUpdate.tgt_snode_id-1] = CreateSuperNode(this->options_.decomposition);
                                }
                                aggVectors[curUpdate.tgt_snode_id-1]->Init(curUpdate.tgt_snode_id, this->Xsuper_[curUpdate.tgt_snode_id-1], this->Xsuper_[curUpdate.tgt_snode_id-1], this->Xsuper_[curUpdate.tgt_snode_id]-1, this->iSize_,structure,this->options_.panel);

                              }
                              SYMPACK_TIMER_STOP(UPD_ANC_Agg_tmp_creat);
                            }
                            tgt_aggreg = aggVectors[curUpdate.tgt_snode_id-1];

                            SYMPACK_TIMER_STOP(UPD_ANC_Agg_tmp);
                          }

#ifdef _DEBUG_
                          logfileptr->OFS()<<"RECV Supernode "<<curUpdate.tgt_snode_id<<" is updated by Supernode "<<cur_src_snode->Id()<<" rows "<<curUpdate.src_first_row<<" "<<curUpdate.blkidx<<std::endl;
#endif

                          //Update the aggregate
                          SYMPACK_TIMER_START(UPD_ANC_UPD);
                          tgt_aggreg->UpdateAggregate(cur_src_snode,curUpdate,tmpBufs,iTarget,this->iam);
                          SYMPACK_TIMER_STOP(UPD_ANC_UPD);

                          --UpdatesToDo[curUpdate.tgt_snode_id-1];
#ifdef _DEBUG_
                          logfileptr->OFS()<<UpdatesToDo[curUpdate.tgt_snode_id-1]<<" updates left for Supernode "<<curUpdate.tgt_snode_id<<std::endl;
#endif

                          //TODO if I am the last one updating that target, send it
                          //THIS SHOULD NOT HAVE TO BE PROTECTED OR BE ATOMICAL BECAUSE NO OTHER RUNNING TASK SHOULD UPDATE THE SAME TARGET
                          //Send the aggregate if it's the last
                          //If this is my last update sent it to curUpdate.tgt_snode_id

                          SYMPACK_TIMER_START(UPD_ANC_Agg_Send);
                          if (UpdatesToDo[curUpdate.tgt_snode_id-1]==0) {
                            if (iTarget != this->iam) {
#ifdef _DEBUG_
                              logfileptr->OFS()<<"Remote Supernode "<<curUpdate.tgt_snode_id<<" is updated by Supernode "<<cur_src_snode->Id()<<std::endl;
#endif

                              tgt_aggreg->Shrink();

                              MsgMetadata meta;

                              NZBlockDesc & nzblk_desc = tgt_aggreg->GetNZBlockDesc(0);
                              T* nzval_ptr = tgt_aggreg->GetNZval(0);

                              //this is an aggregate
                              meta.src = curUpdate.src_snode_id;
                              meta.tgt = curUpdate.tgt_snode_id;
                              meta.GIndex = nzblk_desc.GIndex;

                              (*(*this->remDealloc))++;

                              char buf[100];
                              sprintf(buf,"%d_%d_%d_%d",meta.src,meta.tgt,0,(Int)Factorization::op_type::AGGREGATE);
                              meta.id = hash_fn(std::string(buf));


                              upcxx::global_ptr<char> sendPtr = upcxx::to_global_ptr(tgt_aggreg->GetStoragePtr(meta.GIndex));
                              //the size of the message is the number of bytes between sendPtr and the address of nzblk_desc
                              size_t msgSize = tgt_aggreg->StorageSize();
                              Int liTarget = this->group_->L2G(iTarget);
                              signal_data(sendPtr, msgSize, liTarget, meta);
                            }
                          }
                          SYMPACK_TIMER_STOP(UPD_ANC_Agg_Send);

                          SYMPACK_TIMER_START(UPD_ANC_Upd_Deps);
                          if (iTarget == this->iam)
                          {
                            char buf[100];
                            sprintf(buf,"%d_%d_%d_%d",curUpdate.tgt_snode_id,curUpdate.tgt_snode_id,0,(Int)Factorization::op_type::FACTOR);
                            auto id = hash_fn(std::string(buf));
                            if ( auto ptr = graph.updateTask(id,1,0) ) {
                              bassert(ptr!=nullptr);
                              scheduler_new_->push(ptr);
                            }
                          }
                          SYMPACK_TIMER_STOP(UPD_ANC_Upd_Deps);

                          //if local update, push a new task in the queue and stop the while loop
                          if (this->iam==iSrcOwner ) {
                            break;
                          }
                        }
                      }
                    }
                    SYMPACK_TIMER_STOP(UPDATE_ANCESTORS);

                  }
                };
              }
              else
#endif
              {
                Task.execute = [&graph,hash_fn,&UpdatesToDo,&aggVectors,dec_ref,this,src,tgt,pTask,type] () {
                  scope_timer(a,FB_UPDATE_TASK);
                  Int iLocalTGT = snodeLocalIndex(tgt);
                  Int src_snode_id = src;
                  Int tgt_snode_id = tgt;

                  src_snode_id = abs(src_snode_id);
                  bool is_first_local = src <0;

                  SuperNode<T> * cur_src_snode; 
                  std::shared_ptr<SuperNode<T> > shptr_cur_src_snode = nullptr; 
#ifdef SP_THREADS
                  std::thread::id tid = std::this_thread::get_id();
                  bassert(tid!=this->scheduler_new_->main_tid || Multithreading::NumThread==1);
#endif

                  Int iSrcOwner = this->Mapping_->Map(abs(src_snode_id)-1,abs(src_snode_id)-1);

                  {
                    IncomingMessage * structPtr = nullptr;
                    std::shared_ptr<IncomingMessage> msgPtr = nullptr;
                    std::shared_ptr<ChainedMessage<SuperNodeBase<T> > > newMsgPtr = nullptr;
                    //Local or remote factor
                    //we have only one local or one remote incoming aggregate
                    SnodeUpdate curUpdate;
                    bool found = false;

                    if (pTask->getData().size()==0) {
                      cur_src_snode = snodeLocal(src_snode_id);
                    }
                    else {
                      scope_timer(b,FB_UPD_UNPACK_MSG);
                      auto msgit = pTask->getData().begin();
                      msgPtr = *msgit;
                      bassert(msgPtr->IsDone());

                      //GET MY ID
                      auto id = pTask->id;

                      if (msgPtr->meta.id == id) {
                        char* dataPtr = msgPtr->GetLocalPtr().get();
                        {
                          scope_timer(c,FB_UPD_UNPACK_MSG_CREATE);
                          shptr_cur_src_snode.reset(CreateSuperNode(this->options_.decomposition,dataPtr,msgPtr->Size(),msgPtr->meta.GIndex));
                          cur_src_snode = shptr_cur_src_snode.get();
                        }
                        {
                          scope_timer(d,FB_UPD_UNPACK_MSG_INIT_TREE);
                          cur_src_snode->InitIdxToBlk();
                        }

                        //TODO add the message to other local updates and update their remote dependencies

                        {
                          scope_timer(a,ENQUEUING_UPDATE_MSGS);
                          SnodeUpdate localUpdate;

                          {
                            while (cur_src_snode->FindNextUpdate(localUpdate,this->Xsuper_,this->SupMembership_,this->iam==iSrcOwner)) {

                              //skip if this update is "lower"
                              if (localUpdate.tgt_snode_id<tgt) {
                                continue;
                              }
                              else {
                                Int iUpdater = this->Mapping_->Map(localUpdate.tgt_snode_id-1,localUpdate.src_snode_id-1);
                                if (iUpdater==this->iam) {
                                  if (localUpdate.tgt_snode_id==tgt) {
                                    curUpdate = localUpdate;
                                    found = true;
                                  }
                                  else {

                                    char buf[100];
                                    sprintf(buf,"%d_%d_%d_%d",localUpdate.src_snode_id,localUpdate.tgt_snode_id,0,(Int)Factorization::op_type::UPDATE);
                                    auto id = hash_fn(std::string(buf));
                                    auto taskit = graph.find_task(id);
                                    bassert(taskit!=graph.tasks_.end());
                                    if (newMsgPtr==nullptr) {
                                      auto base_ptr = std::static_pointer_cast<SuperNodeBase<T> >(shptr_cur_src_snode);
                                      newMsgPtr = std::make_shared<ChainedMessage<SuperNodeBase<T> > >(  base_ptr  ,msgPtr);
                                    }

                                    //this is where we put the msg in the list
                                    auto base_ptr = std::static_pointer_cast<IncomingMessage>(newMsgPtr);
                                    taskit->second->addData( base_ptr );
                                    //                        dec_ref(taskit,0,1);

                                    if ( auto ptr = graph.updateTask(id,0,1) ) {
                                      bassert(ptr!=nullptr);
                                      scheduler_new_->push(ptr);
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                      else {
                        newMsgPtr = std::dynamic_pointer_cast<ChainedMessage<SuperNodeBase<T> > >(msgPtr);
                        cur_src_snode = dynamic_cast<SuperNode<T> *>(newMsgPtr->data.get());
                      }
                    }

                    //TODO UPDATE do my update here
                    {

                      SYMPACK_TIMER_START(UPDATE_ANCESTORS);
                      if (!found) {
                        while (cur_src_snode->FindNextUpdate(curUpdate,this->Xsuper_,this->SupMembership_,this->iam==iSrcOwner)) {

                          //skip if this update is "lower"
                          if (curUpdate.tgt_snode_id<tgt) {
                            continue;
                          }
                          else {
                            Int iUpdater = this->Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.src_snode_id-1);
                            if (iUpdater==this->iam) {
                              if (curUpdate.tgt_snode_id==tgt) {
                                found = true;
                                break;
                              }
                            }

                            if (curUpdate.tgt_snode_id>tgt) {
                              break;
                            }
                          }
                        }
                      }

                      bassert(found);
                      Int iUpdater = this->Mapping_->Map(curUpdate.tgt_snode_id-1,cur_src_snode->Id()-1);
                      bassert(iUpdater == this->iam);

#ifdef _DEBUG_PROGRESS_
                      logfileptr->OFS()<<"implicit Task: {"<<curUpdate.src_snode_id<<" -> "<<curUpdate.tgt_snode_id<<"}"<<std::endl;
                      logfileptr->OFS()<<"Processing update from Supernode "<<curUpdate.src_snode_id<<" to Supernode "<<curUpdate.tgt_snode_id<<std::endl;
#endif


                      SuperNode<T> * tgt_aggreg;
                      Int iTarget = this->Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
                      if (iTarget == this->iam) {
                        //the aggregate std::vector is directly the target snode
                        SYMPACK_TIMER_START(UPD_ANC_Agg_local);
                        tgt_aggreg = snodeLocal(curUpdate.tgt_snode_id);
                        assert(curUpdate.tgt_snode_id == tgt_aggreg->Id());
                        SYMPACK_TIMER_STOP(UPD_ANC_Agg_local);
                      }
                      else {
                        SYMPACK_TIMER_START(UPD_ANC_Agg_tmp);
                        //Check if src_snode_id already have an aggregate std::vector
                        bool creation_needed = false;
                        if (aggVectors[curUpdate.tgt_snode_id-1]==nullptr) {
                          creation_needed = true;
                        }
                        else if (aggVectors[curUpdate.tgt_snode_id-1]->StorageSize()==0) {
                          creation_needed = true;
                        }

                        if (creation_needed) {
                          SYMPACK_TIMER_START(UPD_ANC_Agg_tmp_creat);
                          //use number of rows below factor as initializer

                          //TODO do a customized version for FANIN as we have all the factors locally
                          // the idea is the following: do a DFS and stop each exploration at the first local descendant of current node
#ifdef FANIN_OPTIMIZATION
                          if (this->options_.mappingTypeStr ==  "COL2D")
                          {
                            std::set<Idx> structure;
                            std::list<Int> frontier;
                            {
                              scope_timer(a,MERGE_STRUCTURE_FANIN);
                              Idx tgt_fc = this->Xsuper_[curUpdate.tgt_snode_id-1];
                              dfs_traversal(chSupTree_,curUpdate.tgt_snode_id,frontier);
                              //process frontier in decreasing order of nodes and merge their structure
                              for (auto it = frontier.rbegin(); it!=frontier.rend(); it++) {
                                Int I = *it;
                                SuperNode<T> * source = snodeLocal(I);
                                for (Int blkidx=0; blkidx< source->NZBlockCnt();blkidx++) {
                                  NZBlockDesc & nzblk_desc = source->GetNZBlockDesc(blkidx);
                                  Idx fr = nzblk_desc.GIndex;
                                  Idx lr = source->NRows(blkidx) + fr;
                                  for (Idx row = fr; row<lr;row++) {
                                    if (row>=tgt_fc) {
                                      structure.insert(row);
                                    }
                                  }
                                }
                              }
                            }

                            if (aggVectors[curUpdate.tgt_snode_id-1]==nullptr) {
                              aggVectors[curUpdate.tgt_snode_id-1] = CreateSuperNode(this->options_.decomposition);
                            }
                            aggVectors[curUpdate.tgt_snode_id-1]->Init(curUpdate.tgt_snode_id,this->Xsuper_[curUpdate.tgt_snode_id-1],
                                this->Xsuper_[curUpdate.tgt_snode_id-1], this->Xsuper_[curUpdate.tgt_snode_id]-1, this->iSize_,structure,this->options_.panel);
                          } 
                          else
#endif
                          {
                            std::set<Idx> structure;
                            {
                              scope_timer(a,FETCH_REMOTE_STRUCTURE);

                              upcxx::global_ptr<char> remoteDesc = std::get<0>(remoteFactors_[curUpdate.tgt_snode_id-1]);

                              Int block_cnt = std::get<1>(remoteFactors_[curUpdate.tgt_snode_id-1]);


                              //allocate space to receive block descriptors
                              char * buffer = (char*)UpcxxAllocator::allocate(sizeof(NZBlockDesc)*block_cnt+ sizeof(SuperNodeDesc));
                              upcxx::global_ptr<char> remote = remoteDesc;
                              upcxx::rget(remote, (char*)&buffer[0],block_cnt*sizeof(NZBlockDesc)+sizeof(SuperNodeDesc)).wait();
                              SuperNodeDesc * pdesc = (SuperNodeDesc*)buffer;
                              NZBlockDesc * bufferBlocks = (NZBlockDesc*)(pdesc+1);

                              for (Int i =block_cnt-1;i>=0;i--) {
                                NZBlockDesc & curdesc = bufferBlocks[i];
                                size_t end = (i>0)?bufferBlocks[i-1].Offset:pdesc->nzval_cnt_;
                                Int numRows = (end-curdesc.Offset)/pdesc->iSize_;

                                for (Idx row = 0; row<numRows;row++) {
                                  structure.insert(curdesc.GIndex+row);
                                }
                              }
                              UpcxxAllocator::deallocate((char*)buffer);
                            }
                            if (aggVectors[curUpdate.tgt_snode_id-1]==nullptr) {
                              aggVectors[curUpdate.tgt_snode_id-1] = CreateSuperNode(this->options_.decomposition);
                            }
                            aggVectors[curUpdate.tgt_snode_id-1]->Init(curUpdate.tgt_snode_id, this->Xsuper_[curUpdate.tgt_snode_id-1],
                                this->Xsuper_[curUpdate.tgt_snode_id-1], this->Xsuper_[curUpdate.tgt_snode_id]-1, this->iSize_,structure,this->options_.panel);

                          }
                          SYMPACK_TIMER_STOP(UPD_ANC_Agg_tmp_creat);
                        }
                        tgt_aggreg = aggVectors[curUpdate.tgt_snode_id-1];

                        SYMPACK_TIMER_STOP(UPD_ANC_Agg_tmp);
                      }

#ifdef _DEBUG_
                      logfileptr->OFS()<<"RECV Supernode "<<curUpdate.tgt_snode_id<<" is updated by Supernode "<<cur_src_snode->Id()<<" rows "<<curUpdate.src_first_row<<" "<<curUpdate.blkidx<<std::endl;
#endif


                      supernode_lock<T> lock(tgt_aggreg);
                      //Update the aggregate
                      SYMPACK_TIMER_START(UPD_ANC_UPD);
#ifdef SP_THREADS
                      if (Multithreading::NumThread>1) {
                        auto & tmpBuf = tmpBufs_th[tid];

                        tgt_aggreg->UpdateAggregate(cur_src_snode,curUpdate,tmpBuf,iTarget,this->iam);
                      }
                      else
#endif
                        tgt_aggreg->UpdateAggregate(cur_src_snode,curUpdate,tmpBufs,iTarget,this->iam);

                      SYMPACK_TIMER_STOP(UPD_ANC_UPD);
                      --UpdatesToDo[curUpdate.tgt_snode_id-1];
#ifdef _DEBUG_
                      logfileptr->OFS()<<UpdatesToDo[curUpdate.tgt_snode_id-1]<<" updates left for Supernode "<<curUpdate.tgt_snode_id<<std::endl;
#endif
                      SYMPACK_TIMER_STOP(UPDATE_ANCESTORS);




                      //TODO if I am the last one updating that target, send it
                      //THIS SHOULD NOT HAVE TO BE PROTECTED OR BE ATOMICAL BECAUSE NO OTHER RUNNING TASK SHOULD UPDATE THE SAME TARGET
                      //Send the aggregate if it's the last
                      //If this is my last update sent it to curUpdate.tgt_snode_id
                      SYMPACK_TIMER_START(UPD_ANC_Agg_Send);
                      if (UpdatesToDo[curUpdate.tgt_snode_id-1]==0) {
                        if (iTarget != this->iam) {
#ifdef _DEBUG_
                          logfileptr->OFS()<<"Remote Supernode "<<curUpdate.tgt_snode_id<<" is updated by Supernode "<<cur_src_snode->Id()<<std::endl;
#endif

                          tgt_aggreg->Shrink();

                          MsgMetadata meta;

                          NZBlockDesc & nzblk_desc = tgt_aggreg->GetNZBlockDesc(0);
                          T* nzval_ptr = tgt_aggreg->GetNZval(0);

                          //this is an aggregate
                          meta.src = curUpdate.src_snode_id;
                          meta.tgt = curUpdate.tgt_snode_id;
                          meta.GIndex = nzblk_desc.GIndex;

                          (*(*this->remDealloc))++;

                          char buf[100];
                          sprintf(buf,"%d_%d_%d_%d",meta.src,meta.tgt,0,(Int)Factorization::op_type::AGGREGATE);
                          meta.id = hash_fn(std::string(buf));

                          upcxx::global_ptr<char> sendPtr = upcxx::to_global_ptr(tgt_aggreg->GetStoragePtr(meta.GIndex));
                          //the size of the message is the number of bytes between sendPtr and the address of nzblk_desc
                          size_t msgSize = tgt_aggreg->StorageSize();
                          Int liTarget = this->group_->L2G(iTarget);
                          signal_data(sendPtr, msgSize, liTarget, meta);
                        }
                      }
                      SYMPACK_TIMER_STOP(UPD_ANC_Agg_Send);

                      SYMPACK_TIMER_START(UPD_ANC_Upd_Deps);
                      if (iTarget == this->iam)
                      {
                        char buf[100];
                        sprintf(buf,"%d_%d_%d_%d",curUpdate.tgt_snode_id,curUpdate.tgt_snode_id,0,(Int)Factorization::op_type::FACTOR);
                        auto id = hash_fn(std::string(buf));


                        auto taskit = graph.find_task(id);
                        bassert(taskit!=graph.tasks_.end());
                        if ( auto ptr = graph.updateTask(id,1,0) ) {
                          bassert(ptr!=nullptr);
                          scheduler_new_->push(ptr);
                        }
                      }
                      SYMPACK_TIMER_STOP(UPD_ANC_Upd_Deps);
                    }

                    if (structPtr!=nullptr) {
                      delete structPtr;
                    }

                  }
                };

              }
            }
            break;
          default: break; // silence -Wswitch warnings from clang
        }
      }
    }
  }
  SYMPACK_TIMER_STOP(FB_INIT);

  if (this->iam==0 && this->options_.verbose) {
    symPACKOS<<"TaskGraph size is: "<<graph.tasks_.size()<<std::endl;
  }
  timeSta = get_time();
  scheduler_new_->run(CommEnv_->MPI_GetComm(),*this->group_,graph,*this->remDealloc,*this->workteam_);
  delete this->remDealloc;
  double timeStop = get_time();
  if (this->iam==0 && this->options_.verbose) {
    symPACKOS<<"Factorization task graph execution time: "<<timeStop - timeSta<<std::endl;
  }

  tmpBufs.Clear();
}

template <typename T> inline void symPACKMatrix<T>::FBGetUpdateCount(std::vector<Int> & UpdatesToDo, std::vector<Int> & AggregatesToRecv,std::vector<Int> & LocalAggregates)
{
  SYMPACK_TIMER_START(FB_GET_UPDATE_COUNT);
  UpdatesToDo.resize(this->Xsuper_.size(),I_ZERO);
  AggregatesToRecv.resize(this->Xsuper_.size(),I_ZERO);
  LocalAggregates.resize(this->Xsuper_.size(),I_ZERO);
  std::vector<Int> marker(this->Xsuper_.size(),I_ZERO);

  //map of map of pairs (supno, count)
  std::map<Idx, std::map<Idx, Idx>  > Updates;
  std::map<Idx, std::map<Idx, Idx>  > sendAfter;

  Int numLocSnode = this->XsuperDist_[this->iam+1]-this->XsuperDist_[this->iam];
  Int firstSnode = this->XsuperDist_[this->iam];
  Int lastSnode = firstSnode + numLocSnode-1;

  for (Int locsupno = 1; locsupno<this->locXlindx_.size(); ++locsupno) {
    Idx s = locsupno + firstSnode-1;

    Int first_col = this->Xsuper_[s-1];
    Int last_col = this->Xsuper_[s]-1;

    Ptr lfi = this->locXlindx_[locsupno-1];
    Ptr lli = this->locXlindx_[locsupno]-1;
    Idx prevSnode = (Idx)-1;
    for (Ptr sidx = lfi; sidx<=lli;sidx++) {
      Idx row = this->locLindx_[sidx-1];
      Int supno = this->SupMembership_[row-1];
      if (supno!=prevSnode) {
        if (marker[supno-1]!=s && supno!=s) {
          marker[supno-1] = s;

          Int iFactorizer = -1;
          Int iUpdater = -1;
          iFactorizer = Mapping_->Map(supno-1,supno-1);
          iUpdater = Mapping_->Map(supno-1,s-1);

          if ( iUpdater==iFactorizer) {
            LocalAggregates[supno-1]++;
          }
          else {
            auto & procSendAfter = sendAfter[iUpdater];
            auto it = procSendAfter.find(supno);
            if (it==procSendAfter.end()) {
              procSendAfter[supno]=s;
            }
          }

          auto & procUpdates = Updates[iUpdater];
          procUpdates[supno]++;
        }
      }
      prevSnode = supno;
    }
  }


  //Build a Alltoallv communication for Updates
  {
    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(std::pair<Idx,Idx>), MPI_BYTE, &type );
    MPI_Type_commit(&type);

    //compute send sizes
    vector<int> ssizes(this->np,0);
    int numPairs = 0;
    for (auto itp = Updates.begin();itp!=Updates.end();itp++) {
      ssizes[itp->first] = itp->second.size();
      numPairs+=itp->second.size();
    }

    //compute send displacements
    vector<int> sdispls(this->np+1,0);
    sdispls[0] = 0;
    std::partial_sum(&ssizes.front(),&ssizes.back(),&sdispls[1]);

    //Build the contiguous array of pairs
    vector<std::pair<Idx, Idx> > sendbuf;
    sendbuf.reserve(numPairs);
    for (auto itp = Updates.begin();itp!=Updates.end();itp++) {
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

    vector<std::pair<Idx, Idx> > recvbuf(rdispls.back());
    MPI_Alltoallv(&sendbuf[0],&ssizes[0],&sdispls[0],type,&recvbuf[0],&rsizes[0],&rdispls[0],type,CommEnv_->MPI_GetComm());


    std::map<Idx, Idx> LocalUpdates;
    for (auto it = recvbuf.begin();it!=recvbuf.end();it++) {
      auto it2 = LocalUpdates.find(it->first);
      if (it2!=LocalUpdates.end()) {
        it2->second += it->second;
      }
      else {
        LocalUpdates[it->first] = it->second;
      }
    }

    MPI_Type_free(&type);


    for (auto it = LocalUpdates.begin();it!=LocalUpdates.end();it++) {
      UpdatesToDo[it->first-1] = it->second;
    }
  }

  //Build a Alltoallv communication for sendAfter
  {
    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(std::pair<Idx,Idx>), MPI_BYTE, &type );
    MPI_Type_commit(&type);
    //compute send sizes
    vector<int> ssizes(this->np,0);
    int numPairs = 0;
    for (auto itp = sendAfter.begin();itp!=sendAfter.end();itp++) {
      ssizes[itp->first] = itp->second.size();
      numPairs+=itp->second.size();
    }

    //compute send displacements
    vector<int> sdispls(this->np+1,0);
    sdispls[0] = 0;
    std::partial_sum(&ssizes.front(),&ssizes.back(),&sdispls[1]);

    //Build the contiguous array of pairs
    vector<std::pair<Idx, Idx> > sendbuf;
    sendbuf.reserve(numPairs);
    for (auto itp = sendAfter.begin();itp!=sendAfter.end();itp++) {
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
    vector<std::pair<Idx, Idx> > recvbuf(rdispls.back());
    MPI_Alltoallv(&sendbuf[0],&ssizes[0],&sdispls[0],type,&recvbuf[0],&rsizes[0],&rdispls[0],type,CommEnv_->MPI_GetComm());



    std::map<Idx, Idx> mLocalSendAfter;
    for (auto it = recvbuf.begin();it!=recvbuf.end();it++) {
      auto it2 = mLocalSendAfter.find(it->first);
      if (it2!=mLocalSendAfter.end()) {
        if (it2->second>it->second) {
          it2->second = it->second;
        }
      }
      else {
        mLocalSendAfter[it->first] = it->second;
      }
    }

    //do an allgatherv ?
    sendbuf.resize(mLocalSendAfter.size());
    std::copy(mLocalSendAfter.begin(),mLocalSendAfter.end(),sendbuf.begin());

    int sendsize = sendbuf.size();
    MPI_Allgather(&sendsize,sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,CommEnv_->MPI_GetComm());
    rdispls[0] = 0;
    std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);

    recvbuf.resize(rdispls.back());
    MPI_Allgatherv(&sendbuf[0],sendsize,type,&recvbuf[0],&rsizes[0],&rdispls[0],type,CommEnv_->MPI_GetComm());

    //rebuild a map structure
    sendAfter.clear();
    for (int p=0;p<this->np;p++) {
      int start = rdispls[p];
      int end = rdispls[p+1];

      for (int idx = start; idx<end;idx++) {
        sendAfter[p][recvbuf[idx].first] = recvbuf[idx].second;
      }
    }

    std::fill(marker.begin(),marker.end(),I_ZERO);
    for (Int locsupno = 1; locsupno<this->locXlindx_.size(); ++locsupno) {
      Idx s = locsupno + firstSnode-1;

      Int first_col = this->Xsuper_[s-1];
      Int last_col = this->Xsuper_[s]-1;

      Ptr lfi = this->locXlindx_[locsupno-1];
      Ptr lli = this->locXlindx_[locsupno]-1;
      Idx prevSnode = (Idx)-1;
      for (Ptr sidx = lfi; sidx<=lli;sidx++) {
        Idx row = this->locLindx_[sidx-1];
        Int supno = this->SupMembership_[row-1];

        if (supno!=prevSnode) {
          if (marker[supno-1]!=s && supno!=s) {
            marker[supno-1] = s;

            Int iFactorizer = -1;
            Int iUpdater = -1;
            iFactorizer = Mapping_->Map(supno-1,supno-1);
            iUpdater = Mapping_->Map(supno-1,s-1);


            if ( iUpdater!=iFactorizer) {
              auto & procSendAfter = sendAfter[iUpdater];
              if (s==procSendAfter[supno]) {
                AggregatesToRecv[supno-1]++;
              }
            }
          }
        }
        prevSnode = supno;
      }
    }
    MPI_Type_free(&type);
  }

  MPI_Allreduce(MPI_IN_PLACE,&AggregatesToRecv[0],AggregatesToRecv.size(),MPI_INT,MPI_SUM,CommEnv_->MPI_GetComm());
  MPI_Allreduce(MPI_IN_PLACE,&LocalAggregates[0],LocalAggregates.size(),MPI_INT,MPI_SUM,CommEnv_->MPI_GetComm());

  SYMPACK_TIMER_STOP(FB_GET_UPDATE_COUNT);
}

#endif //_SYMPACK_MATRIX_IMPL_FB_PULL_HPP_
