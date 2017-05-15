/*
   Copyright (c) 2016 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.  

Author: Mathias Jacquelin

This file is part of symPACK. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to Lawrence Berkeley National
Laboratory, without imposing a separate written license agreement for such
Enhancements, then you hereby grant the following license: a non-exclusive,
royalty-free perpetual license to install, use, modify, prepare derivative
works, incorporate into other computer software, distribute, and sublicense
such enhancements or derivative works thereof, in binary and source code form.
*/
#ifndef _SYMPACK_MATRIX_IMPL_FB_PULL_HPP_
#define _SYMPACK_MATRIX_IMPL_FB_PULL_HPP_

#define FANIN_OPTIMIZATION

template <typename T> inline void symPACKMatrix<T>::dfs_traversal(std::vector<std::list<Int> > & tree,int node,std::list<Int> & frontier){
  for(std::list<Int>::iterator it = tree[node].begin(); it!=tree[node].end(); it++){
    Int I = *it;

    Int iOwner = Mapping_->Map(I-1,I-1);
    if(iOwner == iam){
      frontier.push_back(I);
    }
    else{
      if(I<tree.size()){
        dfs_traversal(tree,I,frontier);
      }
    }
  }
};


template <typename T> inline void symPACKMatrix<T>::FanBoth_New() 
{
  scope_timer(a,FACTORIZATION_FB);

  double timeSta, timeEnd;

  SYMPACK_TIMER_START(FB_INIT);

#ifdef UPCXX_PROGRESS_THREAD
  upcxx::progress_thread_start();
#endif

  std::vector<Int> UpdatesToDo = UpdatesToDo_;

  //tmp buffer space
  std::vector<T> src_nzval;
  std::vector<char> src_blocks;


  Int maxheight = 0;
  Int maxwidth = 0;
  for(Int i = 1; i<Xsuper_.size(); ++i){
    Int width =Xsuper_[i] - Xsuper_[i-1];
    if(width>=maxwidth){
      maxwidth = width;
    }
    if(UpdateHeight_[i-1]>=maxheight){
      maxheight = UpdateHeight_[i-1];
    }
  }

  //maxwidth for indefinite matrices
  tmpBufs.Resize(maxwidth + maxheight,maxwidth);
  std::vector< SuperNode<T>* > aggVectors(Xsuper_.size()-1,NULL);

  timeSta =  get_time( );
  SYMPACK_TIMER_START(BUILD_TASK_LIST);

  //Create a copy of the task graph
  taskGraph graph = taskGraph_New_;
  //std::shared_ptr<Scheduler< std::shared_ptr<GenericTask> > > scheduler(new FIFOScheduler< std::shared_ptr<GenericTask> >( ));

  std::hash<std::string> hash_fn;

//  std::map<Int,Int> factorUser;
//  std::mutex factorinuse_mutex_;

  //This is the important one
#ifdef SP_THREADS
//  std::map<Int,bool> superNodeInUse;
//  std::mutex inuse_mutex_;
  scheduler_new_->threadInitHandle_ = nullptr;
  scheduler_new_->extraTaskHandle_  = nullptr;
  scheduler_new_->msgHandle = nullptr;

  scheduler_new_->threadInitHandle_ = [&,this](){
              std::thread::id tid = std::this_thread::get_id();
              scheduler_new_->list_mutex_.lock();
              auto & tmpBuf = tmpBufs_th[tid];
              scheduler_new_->list_mutex_.unlock();
  };

  if(Multithreading::NumThread>1){
    scheduler_new_->extraTaskHandle_ = [&,this](std::shared_ptr<GenericTask> & pTask)->bool {
      SparseTask & Task = *(SparseTask*)pTask.get();

      Int * meta = reinterpret_cast<Int*>(Task.meta.data());
      Factorization::op_type & type = *reinterpret_cast<Factorization::op_type*>(&meta[3]);

      if(type == Factorization::op_type::FACTOR){
        return false;
      }
      else if(type == Factorization::op_type::TRSM){
        return false;
      }
      else{
        int tgt = meta[1];


        SuperNode<T> * tgt_aggreg = nullptr;
        Int iTarget = this->Mapping_->Map(tgt-1,tgt-1);
        if(iTarget == iam){
          tgt_aggreg = snodeLocal(tgt);
        }
        else{
          tgt_aggreg = aggVectors[tgt-1];
          if(tgt_aggreg == nullptr){
            aggVectors[tgt-1] = CreateSuperNode(options_.decomposition);
            tgt_aggreg = aggVectors[tgt-1];
          }
        }

        //assert(tgt_aggreg!=nullptr);

        if(!tgt_aggreg->in_use){
          tgt_aggreg->in_use = true;
          return false;
        }
        else{
          return true;
        }
      }
    };
  }

#endif

  auto dec_ref = [&] ( taskGraph::task_iterator taskit, Int loc, Int rem) {
            scope_timer(a,DECREF);
#ifdef SP_THREADS
    if(Multithreading::NumThread>1){
      scheduler_new_->list_mutex_.lock();
    }
#endif
    taskit->second->local_deps-= loc;
    taskit->second->remote_deps-= rem;

    if(taskit->second->remote_deps==0 && taskit->second->local_deps==0){
      {
        //SparseTask * tmp = ((SparseTask*)taskit->second.get());
        //Int * meta = reinterpret_cast<Int*>(tmp->meta.data());
        //Factorization::op_type & type = *reinterpret_cast<Factorization::op_type*>(&meta[2]);
        //if(meta[0]==10 && meta[1]==10){
        // gdb_lock();
        //} 
      }
      scheduler_new_->push(taskit->second);
      graph.removeTask(taskit->second->id);
    }

#ifdef SP_THREADS
    if(Multithreading::NumThread>1){
      scheduler_new_->list_mutex_.unlock();
    }
#endif
  };

  auto log_task = [&] (taskGraph::task_iterator taskit){
    SparseTask * tmp = ((SparseTask*)taskit->second.get());
    Int * meta = reinterpret_cast<Int*>(tmp->meta.data());
    Factorization::op_type & type = *reinterpret_cast<Factorization::op_type*>(&meta[3]);
    std::string name;
    switch(type){
      case Factorization::op_type::FACTOR:
        name="FACTOR";
        break;
      case Factorization::op_type::AGGREGATE:
        name="AGGREGATE";
        break;
      case Factorization::op_type::UPDATE:
        name="UPDATE";
        break;
    }

    logfileptr->OFS()<<" T "<<name<<" "<<meta[0]<<"_"<<meta[1]<<" "<<tmp->local_deps<<" "<<tmp->remote_deps<<std::endl;
  };

  auto log_task_internal = [&] ( const std::shared_ptr<GenericTask> & taskptr){
    SparseTask * tmp = ((SparseTask*)taskptr.get());
    Int * meta = reinterpret_cast<Int*>(tmp->meta.data());
    Factorization::op_type & type = *reinterpret_cast<Factorization::op_type*>(&meta[3]);
    std::string name;
    switch(type){
      case Factorization::op_type::FACTOR:
        name="FACTOR";
        break;
      case Factorization::op_type::AGGREGATE:
        name="AGGREGATE";
        break;
      case Factorization::op_type::UPDATE:
        name="UPDATE";
        break;
    }

    std::stringstream sstr;
    sstr<<" Running T "<<name<<" "<<meta[0]<<"_"<<meta[1]<<" "<<tmp->local_deps<<" "<<tmp->remote_deps<<std::endl;
    logfileptr->OFS()<<sstr.str();
  };

#ifdef FANIN_OPTIMIZATION
  if(options_.mappingTypeStr ==  "COL2D")
  {
    scope_timer(a,BUILD_SUPETREE);
    ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
    chSupTree_.resize(SupETree.Size()+1);
    for(Int I = 1; I<=SupETree.Size(); I++){
      Int parent = SupETree.PostParent(I-1);
      chSupTree_[parent].push_back(I);
    }
  }
#endif

  //Finish the initialization of the tasks by creating the lambdas
  {

    scheduler_new_->msgHandle = [&,this](std::shared_ptr<IncomingMessage> msg) {
      Int iOwner = Mapping_->Map(msg->meta.tgt-1,msg->meta.tgt-1);
      Int iUpdater = Mapping_->Map(msg->meta.tgt-1,msg->meta.src-1);
      //this is an aggregate
      if(iOwner==iam && iUpdater!=iam){
#ifdef _DEBUG_PROGRESS_
        logfileptr->OFS()<<"COMM: FETCHED AGGREG MSG("<<msg->meta.src<<","<<msg->meta.tgt<<") from P"<<iUpdater<<std::endl;
#endif

        //insert an aggregation task
        Int I = msg->meta.src;
        Int J = msg->meta.tgt;
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
        Task.execute = [&,this,src,tgt,pTask] () {
          scope_timer(a,FB_AGGREGATION_TASK);
          Int iLocalI = snodeLocalIndex(tgt);

          //log_task_internal(pTask);

#ifdef _DEBUG_PROGRESS_
          logfileptr->OFS()<<"Processing T_AGGREG("<<src<<","<<tgt<<") "<<std::endl;
#endif

          Int src_snode_id = src;
          Int tgt_snode_id = tgt;
          Int I = src_snode_id;
          auto src_snode = LocalSupernodes_[iLocalI -1];
          Int src_first_col = src_snode->FirstCol();
          Int src_last_col = src_snode->LastCol();

#ifdef _DEBUG_PROGRESS_
          logfileptr->OFS()<<"Processing Supernode "<<I<<std::endl;
#endif


          //Applying aggregates
          //TODO we might have to create a third UPDATE type of task
          //process them one by one
          {
            scope_timer(a,FACT_AGGREGATE);
          assert(pTask->getData().size()==1);
          for(auto msgit = pTask->getData().begin();msgit!=pTask->getData().end();msgit++){
            auto msgPtr = *msgit;
            assert(msgPtr->IsDone());

            if(!msgPtr->IsLocal()){
              msgPtr->DeallocRemote();
            }
            char* dataPtr = msgPtr->GetLocalPtr().get();

            auto dist_src_snode = CreateSuperNode(options_.decomposition,dataPtr,msgPtr->Size());

            dist_src_snode->InitIdxToBlk();

            src_snode->Aggregate(dist_src_snode);

            delete dist_src_snode;
            //delete msgPtr;
          }
          }

#ifdef SP_THREADS
                  src_snode->in_use = false;
#endif


          std::stringstream sstr;
          sstr<<tgt_snode_id<<"_"<<tgt_snode_id<<"_"<<0<<"_"<<(Int)Factorization::op_type::FACTOR;
          auto id = hash_fn(sstr.str());
          auto taskit = graph.find_task(id);
          //this is a particular case : we consider this as a remote dependency
          bassert(taskit!=graph.tasks_.end());
          dec_ref(taskit,0,1);

        };



        std::stringstream sstr;
        sstr<<meta[0]<<"_"<<meta[1]<<"_"<<0<<"_"<<(Int)type;
        Task.id = hash_fn(sstr.str());
        graph.addTask(pTask);

        //auto taskit = graph.find_task(Task.id);
        //log_task(taskit);
      }

    };


#if 1
    for(auto taskit = graph.tasks_.begin(); taskit != graph.tasks_.end(); taskit++){
      auto & pTask = taskit->second;
      SparseTask & Task = *(SparseTask*)pTask.get();

      Int * meta = reinterpret_cast<Int*>(Task.meta.data());
      Int src = meta[0];
      Int tgt = meta[1];
      Factorization::op_type & type = *reinterpret_cast<Factorization::op_type*>(&meta[3]);

      switch(type){
        case Factorization::op_type::FACTOR:
          {
            Task.execute = [&,this,src,tgt,pTask] () {

          //log_task_internal(pTask);
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
//              scheduler->list_mutex_.lock();
              auto & tmpBuf = tmpBufs_th[tid];
//              scheduler->list_mutex_.unlock();
              src_snode->Factorize(tmpBuf);
#else
              src_snode->Factorize(tmpBufs);
#endif

              SYMPACK_TIMER_STOP(FACTOR_PANEL);

#if 1
              //unlock the TRSM tasks
              for(Int blkidx=1; blkidx< src_snode->NZBlockCnt();blkidx++){
                NZBlockDesc & nzblk_desc = src_snode->GetNZBlockDesc(blkidx);
                Idx fr = nzblk_desc.GIndex;
                Idx lr = src_snode->NRows(blkidx) + fr;

                std::stringstream sstr;
                sstr<<src_snode_id<<"_"<<fr<<"_"<<0<<"_"<<(Int)Factorization::op_type::TRSM;
                auto id = hash_fn(sstr.str());

      //logfileptr->OFS()<<"UNLOCK "<<sstr.str()<< "-> "<<Task.id<<std::endl;

                auto taskit = graph.find_task(id);
                if(taskit!=graph.tasks_.end()){
                  dec_ref(taskit,1,0);
                }
                else{
                  abort();
                  //TODO send diagonal block to remote rank
                }
              }
#endif


              //Sending factors and update local tasks
              //Send my factor to my ancestors. 
              std::vector<char> is_factor_sent(np);
              SetValue(is_factor_sent,false);

              SnodeUpdate curUpdate;
              SYMPACK_TIMER_START(FIND_UPDATED_ANCESTORS);
              SYMPACK_TIMER_START(FIND_UPDATED_ANCESTORS_FACTORIZATION);
              while(src_snode->FindNextUpdate(curUpdate,Xsuper_,SupMembership_)){ 
                Int iTarget = this->Mapping_->Map(curUpdate.tgt_snode_id-1,src_snode->Id()-1);

                if(iTarget != iam){
                  if(!is_factor_sent[iTarget]){
                    MsgMetadata meta;

                    //TODO Replace all this by a Serialize function
                    NZBlockDesc & nzblk_desc = src_snode->GetNZBlockDesc(curUpdate.blkidx);
                    Int local_first_row = curUpdate.src_first_row - nzblk_desc.GIndex;
                    Int nzblk_cnt = src_snode->NZBlockCnt() - curUpdate.blkidx;
                    Int nzval_cnt_ = src_snode->Size()*(src_snode->NRowsBelowBlock(curUpdate.blkidx)-local_first_row);
                    T* nzval_ptr = src_snode->GetNZval(nzblk_desc.Offset) + local_first_row*src_snode->Size();

                    upcxx::global_ptr<char> sendPtr((char*)nzval_ptr);
                    //the size of the message is the number of bytes between sendPtr and the address of nzblk_desc

                    //Send factor 
                    meta.src = curUpdate.src_snode_id;
                    meta.tgt = curUpdate.tgt_snode_id;
                    meta.GIndex = curUpdate.src_first_row;

                    std::stringstream sstr;
                    sstr<<meta.src<<"_"<<meta.tgt<<"_"<<0<<"_"<<(Int)Factorization::op_type::UPDATE;
                    meta.id = hash_fn(sstr.str());

                    char * last_byte_ptr = (char*)&nzblk_desc + sizeof(NZBlockDesc);
                    size_t msgSize = last_byte_ptr - (char*)nzval_ptr;

                    {
                      signal_data(sendPtr, msgSize, iTarget, meta);
                    }
                    is_factor_sent[iTarget] = true;
                  }
                }
                else{
                  //Update local tasks
                  //find task corresponding to curUpdate

                  std::stringstream sstr;
                  sstr<<curUpdate.src_snode_id<<"_"<<curUpdate.tgt_snode_id<<"_"<<0<<"_"<<(Int)Factorization::op_type::UPDATE;
                  auto id = hash_fn(sstr.str());
                  auto taskit = graph.find_task(id);
                  bassert(taskit!=graph.tasks_.end());
                  dec_ref(taskit,1,0);
                }
              }
              SYMPACK_TIMER_STOP(FIND_UPDATED_ANCESTORS_FACTORIZATION);
              SYMPACK_TIMER_STOP(FIND_UPDATED_ANCESTORS);
            };



          }
          break;
        case Factorization::op_type::TRSM:
          {
            Task.execute = [&,this,src,tgt,pTask] () {
              scope_timer(b,FB_TRSM_TASK);
              Int firstRow = tgt;
              Int tgt_snode_id = SupMembership_[tgt-1];

#if 0
              std::shared_ptr<IncomingMessage> msgPtr = nullptr;
              std::shared_ptr<ChainedMessage<SuperNodeBase<T> > > newMsgPtr = nullptr;
              //check whether this task as a remote incoming dependence
              SuperNode<T> * src_snode = nullptr;
              if (pTask->getData().size()>0){
                //TODO this part is entirely to debug
                src_snode = LocalSupernodes_[iLocalSRC -1];
                scope_timer(b,FB_TRSM_UNPACK_MSG);
                auto msgit = pTask->getData().begin();
                auto & msgPtr = *msgit;
                bassert(msgPtr->IsDone());

                //GET MY ID
                std::stringstream sstr;
                sstr<<src<<"_"<<tgt<<"_"<<0<<"_"<<(Int)type;
                auto id = hash_fn(sstr.str());

                std::shared_ptr<SuperNode<T> > shptr_src_snode = nullptr; 
                if(msgPtr->meta.id == id){
                  char* dataPtr = msgPtr->GetLocalPtr().get();
                  {
                    scope_timer(c,FB_TRSM_UNPACK_MSG_CREATE);
                    shptr_src_snode.reset(CreateSuperNode(options_.decomposition,dataPtr,msgPtr->Size(),msgPtr->meta.GIndex));
                    src_snode = shptr_src_snode.get();
                  }
                  {
                    scope_timer(d,FB_TRSM_UNPACK_MSG_INIT_TREE);
                    src_snode->InitIdxToBlk();
                  }

                  //TODO add the message to other local TRSMs and update their remote dependencies
                  {
                     scope_timer(a,ENQUEUING_TRSM_MSGS);
                     SnodeUpdate localUpdate;

                     while(tgt_snode->FindNextUpdate(localUpdate,Xsuper_,SupMembership_,false)){
                        //skip if this update is "lower"
                        if(localUpdate.tgt_snode_id<tgt){
                           continue;
                        }
                        else{
                           std::stringstream sstr;
                           sstr<<localUpdate.src_snode_id<<"_"<<localUpdate.src_first_row<<"_"<<0<<"_"<<(Int)Factorization::op_type::TRSM;
                           auto id = hash_fn(sstr.str());
                           auto taskit = graph.find_task(id);

                           //If there is a match, I own the task
                           if(taskit!=graph.tasks_.end()){
                                  if(newMsgPtr==nullptr){
                                    auto base_ptr = std::static_pointer_cast<SuperNodeBase<T> >(shptr_src_snode);
                                    newMsgPtr = std::make_shared<ChainedMessage<SuperNodeBase<T> > >(  base_ptr  ,msgPtr);
                                  }

                                  //this is where we put the msg in the list
                                  auto base_ptr = std::static_pointer_cast<IncomingMessage>(newMsgPtr);
                                  taskit->second->addData( base_ptr );
                                  dec_ref(taskit,0,1);
                           }
                           else{
                             gdb_lock();
                             abort();
                           }
                        }
                     }
                  }
                }
                else{
                  newMsgPtr = std::dynamic_pointer_cast<ChainedMessage<SuperNodeBase<T> > >(msgPtr);
                  src_snode = dynamic_cast<SuperNode<T> *>(newMsgPtr->data.get());
                }
              }
              else{
                Int iLocalSRC = snodeLocalIndex(src);
                src_snode = LocalSupernodes_[iLocalSRC -1];
              }

              Int iLocalTGT = snodeLocalIndex(src);
              tgt_snode = LocalSupernodes_[iLocalTGT -1];

              Int blkidx = tgt_snode->FindBlockIdx(firstRow); 
              bassert(blkidx!=-1);

#ifdef _DEBUG_PROGRESS_
              logfileptr->OFS()<<"TRSM Supernode "<<src<<" from row "<<firstRow<<" block "<<<<std::endl;
#endif
              SYMPACK_TIMER_START(FACTOR_TRSM);
              tgt_snode->Factorize_TRSM(src_snode,blkidx);
              SYMPACK_TIMER_STOP(FACTOR_TRSM);

              //TODO decrement count of the asociated local updates (same row) or send tgt_snode if all TRSMs have been done
              //List of updated supernodes is provided in src_snode

              SnodeUpdate curUpdate;
              while(src_snode->FindNextUpdate(curUpdate,Xsuper_,SupMembership_)){ 
                if(curUpdate.tgt_snode_id<tgt){
                   continue;
                }

                //TODO find the rank id on which 
                Int iTarget = this->Mapping_->Map(curUpdate.tgt_snode_id-1,src_snode->Id()-1);

                if(iTarget != iam){
                  if(!is_factor_sent[iTarget]){
                    MsgMetadata meta;

                    //TODO Replace all this by a Serialize function
                    NZBlockDesc & nzblk_desc = src_snode->GetNZBlockDesc(curUpdate.blkidx);
                    Int local_first_row = curUpdate.src_first_row - nzblk_desc.GIndex;
                    Int nzblk_cnt = src_snode->NZBlockCnt() - curUpdate.blkidx;
                    Int nzval_cnt_ = src_snode->Size()*(src_snode->NRowsBelowBlock(curUpdate.blkidx)-local_first_row);
                    T* nzval_ptr = src_snode->GetNZval(nzblk_desc.Offset) + local_first_row*src_snode->Size();

                    upcxx::global_ptr<char> sendPtr((char*)nzval_ptr);
                    //the size of the message is the number of bytes between sendPtr and the address of nzblk_desc

                    //Send factor 
                    meta.src = curUpdate.src_snode_id;
                    meta.tgt = curUpdate.tgt_snode_id;
                    meta.GIndex = curUpdate.src_first_row;

                    std::stringstream sstr;
                    sstr<<meta.src<<"_"<<meta.tgt<<"_"<<0<<"_"<<(Int)Factorization::op_type::UPDATE;
                    meta.id = hash_fn(sstr.str());

                    char * last_byte_ptr = (char*)&nzblk_desc + sizeof(NZBlockDesc);
                    size_t msgSize = last_byte_ptr - (char*)nzval_ptr;

                    {
                      signal_data(sendPtr, msgSize, iTarget, meta);
                    }
                    is_factor_sent[iTarget] = true;
                  }
                }
                else{
                  //Update local tasks
                  //find task corresponding to curUpdate

                  std::stringstream sstr;
                  sstr<<curUpdate.src_snode_id<<"_"<<curUpdate.tgt_snode_id<<"_"<<0<<"_"<<(Int)Factorization::op_type::UPDATE;
                  auto id = hash_fn(sstr.str());
                  auto taskit = graph.find_task(id);
                  bassert(taskit!=graph.tasks_.end());
                  dec_ref(taskit,1,0);
                }
              }
#endif


            };
          }
          break;
        case Factorization::op_type::UPDATE:
          {
#if 0
            if(Multithreading::NumThread==1){
              Task.execute = [&,this,src,tgt,iLocalTGT,pTask,type] () {
                scope_timer(a,FB_UPDATE_TASK);

                Int src_snode_id = src;
                Int tgt_snode_id = tgt;
                src_snode_id = abs(src_snode_id);
                bool is_first_local = src <0;

                SuperNode<T> * cur_src_snode; 
                std::shared_ptr<SuperNode<T> > shptr_cur_src_snode = nullptr; 
#ifdef SP_THREADS
                std::thread::id tid = std::this_thread::get_id();
#endif
                Int iSrcOwner = this->Mapping_->Map(abs(src_snode_id)-1,abs(src_snode_id)-1);


                IncomingMessage * structPtr = NULL;
                std::shared_ptr<IncomingMessage> msgPtr = nullptr;
                std::shared_ptr<ChainedMessage<SuperNodeBase<T> > > newMsgPtr = nullptr;
                //Local or remote factor
                //we have only one local or one remote incoming aggregate


                SnodeUpdate curUpdate;
                bool found = false;

                if(pTask->getData().size()==0){
                  cur_src_snode = snodeLocal(src_snode_id);
                }
                else{
                  scope_timer(b,FB_UPD_UNPACK_MSG);
                  auto msgit = pTask->getData().begin();
                  msgPtr = *msgit;
                  bassert(msgPtr->IsDone());


                  //GET MY ID
                  std::stringstream sstr;
                  sstr<<src<<"_"<<tgt<<"_"<<0<<"_"<<(Int)type;
                  auto id = hash_fn(sstr.str());

                  if(msgPtr->meta.id == id){
                    char* dataPtr = msgPtr->GetLocalPtr();


                    {
                      scope_timer(c,FB_UPD_UNPACK_MSG_CREATE);
                      shptr_cur_src_snode.reset(CreateSuperNode(options_.decomposition,dataPtr,msgPtr->Size(),msgPtr->meta.GIndex));
                      cur_src_snode = shptr_cur_src_snode.get();
                    }
                    {
                      scope_timer(d,FB_UPD_UNPACK_MSG_INIT_TREE);
                      cur_src_snode->InitIdxToBlk();
                    }




                    //TODO add the message to other local updates and update their remote dependencies

                    {
                      SYMPACK_TIMER_START(UPDATE_ANCESTORS);
                      SnodeUpdate curUpdate;

                      while(cur_src_snode->FindNextUpdate(curUpdate,Xsuper_,SupMembership_,iam==iSrcOwner)){

                        //skip if this update is "lower"
                        if(curUpdate.tgt_snode_id<tgt){
                          continue;
                        }
                        else{
                          Int iUpdater = this->Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.src_snode_id-1);
                          if(iUpdater==iam){


#ifdef _DEBUG_PROGRESS_
                            logfileptr->OFS()<<"implicit Task: {"<<curUpdate.src_snode_id<<" -> "<<curUpdate.tgt_snode_id<<"}"<<std::endl;
                            logfileptr->OFS()<<"Processing update from Supernode "<<curUpdate.src_snode_id<<" to Supernode "<<curUpdate.tgt_snode_id<<std::endl;
#endif


                            SuperNode<T> * tgt_aggreg;
                            Int iTarget = this->Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
                            if(iTarget == iam){
                              //the aggregate std::vector is directly the target snode
                              SYMPACK_TIMER_START(UPD_ANC_Agg_local);
                              tgt_aggreg = snodeLocal(curUpdate.tgt_snode_id);
                              assert(curUpdate.tgt_snode_id == tgt_aggreg->Id());
                              SYMPACK_TIMER_STOP(UPD_ANC_Agg_local);
                            }
                            else{
                              SYMPACK_TIMER_START(UPD_ANC_Agg_tmp);
                              //Check if src_snode_id already have an aggregate std::vector
                              bool creation_needed = false;
                              if(aggVectors[curUpdate.tgt_snode_id-1]==nullptr){
                                creation_needed = true;
                              }
                              else if(aggVectors[curUpdate.tgt_snode_id-1]->StorageSize()==0){
                                creation_needed = true;
                              }

                              if(creation_needed){
                                SYMPACK_TIMER_START(UPD_ANC_Agg_tmp_creat);
                                //use number of rows below factor as initializer

                                //TODO do a customized version for FANIN as we have all the factors locally
                                // the idea is the following: do a DFS and stop each exploration at the first local descendant of current node
#ifdef FANIN_OPTIMIZATION
                                if(options_.mappingTypeStr ==  "COL2D")
                                {
                                  std::set<Idx> structure;
                                  std::list<Int> frontier;
                                  {
                                    scope_timer(a,MERGE_STRUCTURE_FANIN);
                                    Idx tgt_fc = Xsuper_[curUpdate.tgt_snode_id-1];
                                    dfs_traversal(chSupTree_,curUpdate.tgt_snode_id,frontier);
                                    //process frontier in decreasing order of nodes and merge their structure
                                    for(auto it = frontier.rbegin(); it!=frontier.rend(); it++){
                                      Int I = *it;
                                      SuperNode<T> * source = snodeLocal(I);
                                      for(Int blkidx=0; blkidx< source->NZBlockCnt();blkidx++){
                                        NZBlockDesc & nzblk_desc = source->GetNZBlockDesc(blkidx);
                                        Idx fr = nzblk_desc.GIndex;
                                        Idx lr = source->NRows(blkidx) + fr;
                                        for(Idx row = fr; row<lr;row++){
                                          if(row>=tgt_fc){
                                            structure.insert(row);
                                          }
                                        }
                                      }
                                    }
                                  }
                                  //                        aggVectors[curUpdate.tgt_snode_id-1] = CreateSuperNode(options_.decomposition,curUpdate.tgt_snode_id, Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id]-1, iSize_,structure);

                                  if(aggVectors[curUpdate.tgt_snode_id-1]==nullptr){
                                    aggVectors[curUpdate.tgt_snode_id-1] = CreateSuperNode(options_.decomposition);
                                  }
                                  aggVectors[curUpdate.tgt_snode_id-1]->Init(curUpdate.tgt_snode_id, Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id]-1, iSize_,structure);
                                } 
                                else
#endif
                                {
                                  std::set<Idx> structure;
                                  {
                                    scope_timer(a,FETCH_REMOTE_STRUCTURE);
                                    upcxx::global_ptr<SuperNodeDesc> remoteDesc = std::get<0>(remoteFactors_[curUpdate.tgt_snode_id-1]);
                                    Int block_cnt = std::get<1>(remoteFactors_[curUpdate.tgt_snode_id-1]);


                                    //allocate space to receive block descriptors
                                    char * buffer = (char*)UpcxxAllocator::allocate(sizeof(NZBlockDesc)*block_cnt+ sizeof(SuperNodeDesc));
                                    upcxx::global_ptr<char> remote = upcxx::global_ptr<char>(remoteDesc);
                                    {
#ifdef SP_THREADS
                                      if(Multithreading::NumThread>1){
                                        std::lock_guard<upcxx_mutex_type> lock(upcxx_mutex);
                                        upcxx::copy(remote, (char*)&buffer[0],block_cnt*sizeof(NZBlockDesc)+sizeof(SuperNodeDesc));
                                      }
                                      else
#endif
                                        upcxx::copy(remote, (char*)&buffer[0],block_cnt*sizeof(NZBlockDesc)+sizeof(SuperNodeDesc));
                                    }
                                    SuperNodeDesc * pdesc = (SuperNodeDesc*)buffer;
                                    NZBlockDesc * bufferBlocks = (NZBlockDesc*)(pdesc+1);

                                    for(Int i =block_cnt-1;i>=0;i--){
                                      NZBlockDesc & curdesc = bufferBlocks[i];
                                      size_t end = (i>0)?bufferBlocks[i-1].Offset:pdesc->nzval_cnt_;
                                      Int numRows = (end-curdesc.Offset)/pdesc->iSize_;

                                      for(Idx row = 0; row<numRows;row++){
                                        structure.insert(curdesc.GIndex+row);
                                      }
                                    }
                                    UpcxxAllocator::deallocate((char*)buffer);
                                  }
                                  //                        aggVectors[curUpdate.tgt_snode_id-1] = CreateSuperNode(options_.decomposition,curUpdate.tgt_snode_id, Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id]-1, iSize_,structure);
                                  if(aggVectors[curUpdate.tgt_snode_id-1]==nullptr){
                                    aggVectors[curUpdate.tgt_snode_id-1] = CreateSuperNode(options_.decomposition);
                                  }
                                  //                        bassert(aggVectors[curUpdate.tgt_snode_id-1]!=nullptr);
                                  aggVectors[curUpdate.tgt_snode_id-1]->Init(curUpdate.tgt_snode_id, Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id]-1, iSize_,structure);

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
#ifdef SP_THREADS
                            if(Multithreading::NumThread>1){
                              //scheduler->list_mutex_.lock();
                              auto & tmpBuf = tmpBufs_th[tid];
                              //scheduler->list_mutex_.unlock();
                              tgt_aggreg->UpdateAggregate(cur_src_snode,curUpdate,tmpBuf,iTarget,iam);
                            }
                            else
#endif
                              tgt_aggreg->UpdateAggregate(cur_src_snode,curUpdate,tmpBufs,iTarget,iam);

                            SYMPACK_TIMER_STOP(UPD_ANC_UPD);

#ifdef SP_THREADS
                            if(Multithreading::NumThread>1){
                              tgt_aggreg->in_use = false;
                            }
#endif

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
                            if(UpdatesToDo[curUpdate.tgt_snode_id-1]==0){
                              if(iTarget != iam){
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

                                std::stringstream sstr;
                                sstr<<meta.src<<"_"<<meta.tgt<<"_"<<0<<"_"<<(Int)Factorization::op_type::AGGREGATE;
                                meta.id = hash_fn(sstr.str());

                                upcxx::global_ptr<char> sendPtr(tgt_aggreg->GetStoragePtr(meta.GIndex));
                                //the size of the message is the number of bytes between sendPtr and the address of nzblk_desc
                                size_t msgSize = tgt_aggreg->StorageSize();
                                {
                                  signal_data(sendPtr, msgSize, iTarget, meta);
                                }
                              }
                            }
                            SYMPACK_TIMER_STOP(UPD_ANC_Agg_Send);

                            SYMPACK_TIMER_START(UPD_ANC_Upd_Deps);
                            if(iTarget == iam)
                            {
                              std::stringstream sstr;
                              sstr<<curUpdate.tgt_snode_id<<"_"<<curUpdate.tgt_snode_id<<"_"<<0<<"_"<<(Int)Factorization::op_type::FACTOR;
                              auto id = hash_fn(sstr.str());
                              auto taskit = graph.find_task(id);
                              bassert(taskit!=graph.tasks_.end());
                              dec_ref(taskit,1,0);

                            }
                            SYMPACK_TIMER_STOP(UPD_ANC_Upd_Deps);


                            if(structPtr!=NULL){
                              delete structPtr;
                            }


                          }
                        }
                      }

                    }

                  }
                  else{
                    abort();
                  }


                }




              };
            }
            else
#endif
            {
              Task.execute = [&,this,src,tgt,pTask,type] () {
                //log_task_internal(pTask);
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
#endif

                Int iSrcOwner = this->Mapping_->Map(abs(src_snode_id)-1,abs(src_snode_id)-1);

                {
                  IncomingMessage * structPtr = NULL;
                  std::shared_ptr<IncomingMessage> msgPtr = nullptr;
                  std::shared_ptr<ChainedMessage<SuperNodeBase<T> > > newMsgPtr = nullptr;
                  //Local or remote factor
                  //we have only one local or one remote incoming aggregate


                  SnodeUpdate curUpdate;
                  bool found = false;

                  if(pTask->getData().size()==0){
                    cur_src_snode = snodeLocal(src_snode_id);
                  }
                  else{
                    scope_timer(b,FB_UPD_UNPACK_MSG);
                    auto msgit = pTask->getData().begin();
                    msgPtr = *msgit;
                    bassert(msgPtr->IsDone());


                    //GET MY ID
                    std::stringstream sstr;
                    sstr<<src<<"_"<<tgt<<"_"<<0<<"_"<<(Int)type;
                    auto id = hash_fn(sstr.str());

                    if(msgPtr->meta.id == id){
                      char* dataPtr = msgPtr->GetLocalPtr().get();

//                      logfileptr->OFS()<<"2 Data for snode id "<<tgt<<" "<<msgPtr->local_ptr.use_count()<<std::endl;
//                      logfileptr->OFS()<<"3 Data for snode id "<<tgt<<" "<<msgPtr->GetLocalPtr().use_count()<<std::endl;

                      {
                        scope_timer(c,FB_UPD_UNPACK_MSG_CREATE);
                        shptr_cur_src_snode.reset(CreateSuperNode(options_.decomposition,dataPtr,msgPtr->Size(),msgPtr->meta.GIndex));
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
                          //std::lock_guard<std::mutex> lock(factorinuse_mutex_);
                          while(cur_src_snode->FindNextUpdate(localUpdate,Xsuper_,SupMembership_,iam==iSrcOwner)){

                            //skip if this update is "lower"
                            if(localUpdate.tgt_snode_id<tgt){
                              continue;
                            }
                            else{
                              Int iUpdater = this->Mapping_->Map(localUpdate.tgt_snode_id-1,localUpdate.src_snode_id-1);
                              if(iUpdater==iam){
                                if(localUpdate.tgt_snode_id==tgt){
                                  curUpdate = localUpdate;
                                  found = true;
                                }
                                else{

                                  std::stringstream sstr;
                                  sstr<<localUpdate.src_snode_id<<"_"<<localUpdate.tgt_snode_id<<"_"<<0<<"_"<<(Int)Factorization::op_type::UPDATE;
                                  auto id = hash_fn(sstr.str());
                                  auto taskit = graph.find_task(id);

                                  bassert(taskit!=graph.tasks_.end());
                                  if(newMsgPtr==nullptr){
                                    auto base_ptr = std::static_pointer_cast<SuperNodeBase<T> >(shptr_cur_src_snode);
                                    newMsgPtr = std::make_shared<ChainedMessage<SuperNodeBase<T> > >(  base_ptr  ,msgPtr);
                                  }

                                  

                                  //                              factorUser[localUpdate.src_snode_id]++;
                                  //this is where we put the msg in the list
                                  auto base_ptr = std::static_pointer_cast<IncomingMessage>(newMsgPtr);
                                  taskit->second->addData( base_ptr );
                                  dec_ref(taskit,0,1);

//                                  auto tmp = new IncomingMessage();
//                                  *tmp = *msgPtr;
//                                  tmp->meta = msgPtr->meta;
//                                  tmp->local_ptr = msgPtr->GetLocalPtr();
//                                  tmp->msg_size = msgPtr->msg_size;
//                      logfileptr->OFS()<<"4 Chained Data for snode id "<<tgt<<" "<<tmp->GetLocalPtr().use_count()<<std::endl;
//                                  delete tmp;
                                }
                              }
                            }
                          }
                        }

                      }

                    }
                    else{
                      newMsgPtr = std::dynamic_pointer_cast<ChainedMessage<SuperNodeBase<T> > >(msgPtr);
                      cur_src_snode = dynamic_cast<SuperNode<T> *>(newMsgPtr->data.get());
//                      logfileptr->OFS()<<"Chained message for snode id "<<cur_src_snode->Id()<<" "<<newMsgPtr->data.use_count()<<std::endl;
//                      logfileptr->OFS()<<"Chained Data for snode id "<<cur_src_snode->Id()<<" "<<newMsgPtr->GetLocalPtr().use_count()<<std::endl;
                    }


                  }

                  //TODO UPDATE do my update here
                  {

                    SYMPACK_TIMER_START(UPDATE_ANCESTORS);
                    if(!found){
                      while(cur_src_snode->FindNextUpdate(curUpdate,Xsuper_,SupMembership_,iam==iSrcOwner)){

                        //skip if this update is "lower"
                        if(curUpdate.tgt_snode_id<tgt){
                          continue;
                        }
                        else{
                          Int iUpdater = this->Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.src_snode_id-1);
                          if(iUpdater==iam){
                            if(curUpdate.tgt_snode_id==tgt){
                              found = true;
                              break;
                            }
                          }

                          if(curUpdate.tgt_snode_id>tgt){
                            break;
                          }
                        }
                      }
                    }

                    bassert(found);
                    Int iUpdater = this->Mapping_->Map(curUpdate.tgt_snode_id-1,cur_src_snode->Id()-1);
                    bassert(iUpdater == iam);

#ifdef _DEBUG_PROGRESS_
                    logfileptr->OFS()<<"implicit Task: {"<<curUpdate.src_snode_id<<" -> "<<curUpdate.tgt_snode_id<<"}"<<std::endl;
                    logfileptr->OFS()<<"Processing update from Supernode "<<curUpdate.src_snode_id<<" to Supernode "<<curUpdate.tgt_snode_id<<std::endl;
#endif


                    SuperNode<T> * tgt_aggreg;
                    Int iTarget = this->Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
                    if(iTarget == iam){
                      //the aggregate std::vector is directly the target snode
                      SYMPACK_TIMER_START(UPD_ANC_Agg_local);
                      tgt_aggreg = snodeLocal(curUpdate.tgt_snode_id);
                      assert(curUpdate.tgt_snode_id == tgt_aggreg->Id());
                      SYMPACK_TIMER_STOP(UPD_ANC_Agg_local);
                    }
                    else{
                      SYMPACK_TIMER_START(UPD_ANC_Agg_tmp);
                      //Check if src_snode_id already have an aggregate std::vector
                      bool creation_needed = false;
                      if(aggVectors[curUpdate.tgt_snode_id-1]==nullptr){
                        creation_needed = true;
                      }
                      else if(aggVectors[curUpdate.tgt_snode_id-1]->StorageSize()==0){
                        creation_needed = true;
                      }

                      if(creation_needed){
                        SYMPACK_TIMER_START(UPD_ANC_Agg_tmp_creat);
                        //use number of rows below factor as initializer

                        //TODO do a customized version for FANIN as we have all the factors locally
                        // the idea is the following: do a DFS and stop each exploration at the first local descendant of current node
#ifdef FANIN_OPTIMIZATION
                        if(options_.mappingTypeStr ==  "COL2D")
                        {
                          std::set<Idx> structure;
                          std::list<Int> frontier;
                          {
                            scope_timer(a,MERGE_STRUCTURE_FANIN);
                            Idx tgt_fc = Xsuper_[curUpdate.tgt_snode_id-1];
                            dfs_traversal(chSupTree_,curUpdate.tgt_snode_id,frontier);
                            //process frontier in decreasing order of nodes and merge their structure
                            for(auto it = frontier.rbegin(); it!=frontier.rend(); it++){
                              Int I = *it;
                              SuperNode<T> * source = snodeLocal(I);
                              for(Int blkidx=0; blkidx< source->NZBlockCnt();blkidx++){
                                NZBlockDesc & nzblk_desc = source->GetNZBlockDesc(blkidx);
                                Idx fr = nzblk_desc.GIndex;
                                Idx lr = source->NRows(blkidx) + fr;
                                for(Idx row = fr; row<lr;row++){
                                  if(row>=tgt_fc){
                                    structure.insert(row);
                                  }
                                }
                              }
                            }
                          }
                          //                        aggVectors[curUpdate.tgt_snode_id-1] = CreateSuperNode(options_.decomposition,curUpdate.tgt_snode_id, Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id]-1, iSize_,structure);

                          if(aggVectors[curUpdate.tgt_snode_id-1]==nullptr){
                            aggVectors[curUpdate.tgt_snode_id-1] = CreateSuperNode(options_.decomposition);
                          }
                          aggVectors[curUpdate.tgt_snode_id-1]->Init(curUpdate.tgt_snode_id,Xsuper_[curUpdate.tgt_snode_id-1],
                                           Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id]-1, iSize_,structure);
                        } 
                        else
#endif
                        {
                          std::set<Idx> structure;
                          {
                            scope_timer(a,FETCH_REMOTE_STRUCTURE);
                            upcxx::global_ptr<SuperNodeDesc> remoteDesc = std::get<0>(remoteFactors_[curUpdate.tgt_snode_id-1]);
                            Int block_cnt = std::get<1>(remoteFactors_[curUpdate.tgt_snode_id-1]);


                            //allocate space to receive block descriptors
                            char * buffer = (char*)UpcxxAllocator::allocate(sizeof(NZBlockDesc)*block_cnt+ sizeof(SuperNodeDesc));
                            upcxx::global_ptr<char> remote = upcxx::global_ptr<char>(remoteDesc);
                            {
#ifdef SP_THREADS
                              if(Multithreading::NumThread>1){
                                std::lock_guard<upcxx_mutex_type> lock(upcxx_mutex);
                                upcxx::copy(remote, (char*)&buffer[0],block_cnt*sizeof(NZBlockDesc)+sizeof(SuperNodeDesc));
                              }
                              else
#endif
                                upcxx::copy(remote, (char*)&buffer[0],block_cnt*sizeof(NZBlockDesc)+sizeof(SuperNodeDesc));
                            }
                            SuperNodeDesc * pdesc = (SuperNodeDesc*)buffer;
                            NZBlockDesc * bufferBlocks = (NZBlockDesc*)(pdesc+1);

                            for(Int i =block_cnt-1;i>=0;i--){
                              NZBlockDesc & curdesc = bufferBlocks[i];
                              size_t end = (i>0)?bufferBlocks[i-1].Offset:pdesc->nzval_cnt_;
                              Int numRows = (end-curdesc.Offset)/pdesc->iSize_;

                              for(Idx row = 0; row<numRows;row++){
                                structure.insert(curdesc.GIndex+row);
                              }
                            }
                            UpcxxAllocator::deallocate((char*)buffer);
                          }
                          //                        aggVectors[curUpdate.tgt_snode_id-1] = CreateSuperNode(options_.decomposition,curUpdate.tgt_snode_id, Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id]-1, iSize_,structure);
                          if(aggVectors[curUpdate.tgt_snode_id-1]==nullptr){
                            aggVectors[curUpdate.tgt_snode_id-1] = CreateSuperNode(options_.decomposition);
                          }
                          //                        bassert(aggVectors[curUpdate.tgt_snode_id-1]!=nullptr);
                          aggVectors[curUpdate.tgt_snode_id-1]->Init(curUpdate.tgt_snode_id, Xsuper_[curUpdate.tgt_snode_id-1],
                                         Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id]-1, iSize_,structure);

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
#ifdef SP_THREADS
                    if(Multithreading::NumThread>1){
                      //scheduler->list_mutex_.lock();
                      auto & tmpBuf = tmpBufs_th[tid];
                      //scheduler->list_mutex_.unlock();
                      tgt_aggreg->UpdateAggregate(cur_src_snode,curUpdate,tmpBuf,iTarget,iam);
                    }
                    else
#endif
                      tgt_aggreg->UpdateAggregate(cur_src_snode,curUpdate,tmpBufs,iTarget,iam);

                    SYMPACK_TIMER_STOP(UPD_ANC_UPD);

#ifdef SP_THREADS
                    if(Multithreading::NumThread>1){
                      tgt_aggreg->in_use = false;
                    }
#endif

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
                    if(UpdatesToDo[curUpdate.tgt_snode_id-1]==0){
                      if(iTarget != iam){
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

                        std::stringstream sstr;
                        sstr<<meta.src<<"_"<<meta.tgt<<"_"<<0<<"_"<<(Int)Factorization::op_type::AGGREGATE;
                        meta.id = hash_fn(sstr.str());

                        upcxx::global_ptr<char> sendPtr(tgt_aggreg->GetStoragePtr(meta.GIndex));
                        //the size of the message is the number of bytes between sendPtr and the address of nzblk_desc
                        size_t msgSize = tgt_aggreg->StorageSize();
                        {
                          signal_data(sendPtr, msgSize, iTarget, meta);
                        }
                      }
                    }
                    SYMPACK_TIMER_STOP(UPD_ANC_Agg_Send);

                    SYMPACK_TIMER_START(UPD_ANC_Upd_Deps);
                    if(iTarget == iam)
                    {
                      std::stringstream sstr;
                      sstr<<curUpdate.tgt_snode_id<<"_"<<curUpdate.tgt_snode_id<<"_"<<0<<"_"<<(Int)Factorization::op_type::FACTOR;
                      auto id = hash_fn(sstr.str());
                      auto taskit = graph.find_task(id);
                      bassert(taskit!=graph.tasks_.end());
                      dec_ref(taskit,1,0);

                    }
                    SYMPACK_TIMER_STOP(UPD_ANC_Upd_Deps);
                  }

                  if(structPtr!=NULL){
                    delete structPtr;
                  }

                }
              };

            }




          }
          break;
      }
    }
#endif
  }

  SYMPACK_TIMER_STOP(FB_INIT);


#ifndef NDEBUG
//  for(auto taskit = graph.tasks_.begin();taskit!=graph.tasks_.end();taskit++){
//    log_task(taskit);
//  }
#endif


  if(iam==0){
    std::cout<<"TaskGraph size is: "<<graph.tasks_.size()<<std::endl;
  }
  timeSta = get_time();
  scheduler_new_->run(CommEnv_->MPI_GetComm(),graph);
  double timeStop = get_time();
  if(iam==0){
    std::cout<<"Factorization task graph execution time: "<<timeStop - timeSta<<std::endl;
  }



  SYMPACK_TIMER_START(BARRIER);
  upcxx::async_wait();
  MPI_Barrier(CommEnv_->MPI_GetComm());

  SYMPACK_TIMER_STOP(BARRIER);
#ifdef UPCXX_PROGRESS_THREAD
  upcxx::progress_thread_stop();
#endif

  tmpBufs.Clear();
}




template <typename T> inline void symPACKMatrix<T>::FanBoth() 
{
  SYMPACK_TIMER_START(FACTORIZATION_FB);

  SYMPACK_TIMER_START(FB_INIT);
  double timeSta, timeEnd;

#ifdef UPCXX_PROGRESS_THREAD
  upcxx::progress_thread_start();
#endif


  //  Int iam = CommEnv_->MPI_Rank();
  //  Int np  = CommEnv_->MPI_Size();

  std::vector<Int> UpdatesToDo = UpdatesToDo_;

  //tmp buffer space
  std::vector<T> src_nzval;
  std::vector<char> src_blocks;


  Int maxheight = 0;
  Int maxwidth = 0;
  for(Int i = 1; i<Xsuper_.size(); ++i){
    Int width =Xsuper_[i] - Xsuper_[i-1];
    if(width>=maxwidth){
      maxwidth = width;
    }
    if(UpdateHeight_[i-1]>=maxheight){
      maxheight = UpdateHeight_[i-1];
    }
  }

  //maxwidth for indefinite matrices
  tmpBufs.Resize(maxwidth + maxheight/*Size()*/,maxwidth);

  std::vector< SuperNode<T> * > aggVectors(Xsuper_.size()-1,NULL);


  timeSta =  get_time( );
  SYMPACK_TIMER_START(BUILD_TASK_LIST);

  //Create a copy of the task graph
  supernodalTaskGraph<FBTask> taskGraph = taskGraph_;
  size_t gTaskCnt = taskGraph.taskLists_.size();

  Int localTaskCount = localTaskCount_;
  //std::vector<std::list<FBTask> * > taskLists;
  //taskLists.resize(origTaskLists_.size(),NULL);
  //size_t taskCnt = taskLists.size();
  //for(int i = 0; i<taskCnt; ++i){
  //  if(origTaskLists_[i] != NULL){
  //    taskLists[i] = new std::list<FBTask>(); 
  //    taskLists[i]->insert(taskLists[i]->end(),origTaskLists_[i]->begin(),origTaskLists_[i]->end());
  //  }
  //}

  for(int i = 0; i<gTaskCnt; ++i){
    if(taskGraph.taskLists_[i] != NULL){
      auto taskit = taskGraph.taskLists_[i]->begin();
      while (taskit != taskGraph.taskLists_[i]->end())
      {


#ifdef PREFETCH_STRUCTURE
        // if target of the update is remote, then get its structure
        Int iFactorizer = this->Mapping_->Map(taskit->tgt_snode_id-1,taskit->tgt_snode_id-1);
        if(iFactorizer!=iam && taskit->type == UPDATE){
          bool needStructure =  taskit->data.size()>0?(taskit->data.front()->meta.GIndex!=-1):true;
          //is the aggregate already allocated ?
          needStructure = needStructure && aggVectors[taskit->tgt_snode_id-1]==NULL;

          if(taskit->remote_deps==0 && taskit->local_deps==0 && needStructure){
            IncomingMessage * msg_ptr = new IncomingMessage();
            MsgMetadata meta;
            meta.src = taskit->src_snode_id;
            meta.tgt = taskit->tgt_snode_id;
            meta.GIndex = -1;

            msg_ptr->meta = meta;

            upcxx::global_ptr<SuperNodeDesc> pRemote_ptr = std::get<0>(remoteFactors_[meta.tgt-1]);
            size_t pMsg_size = std::get<1>(remoteFactors_[meta.tgt-1])*sizeof(NZBlockDesc)+sizeof(SuperNodeDesc);
            msg_ptr->remote_ptr = upcxx::global_ptr<char>(pRemote_ptr);
            msg_ptr->msg_size = pMsg_size;

            //allocate receive buffer
            bool success = msg_ptr->AllocLocal();
            if(success){
              gIncomingRecvAsync.push_back( msg_ptr );
              msg_ptr->AsyncGet();
              taskit->remote_deps++;
            }
            else{
              abort();
              delete msg_ptr;
            }
          }
        }
#endif








        if(taskit->remote_deps==0 && taskit->local_deps==0){
          scheduler2_->push(*taskit);
          auto it = taskit++;
          taskGraph.removeTask(it);
        }
        else
        {
          ++taskit;
        }
      }
    }
  }








#ifdef FANIN_OPTIMIZATION
  if(options_.mappingTypeStr ==  "COL2D")
  {
    scope_timer(a,BUILD_SUPETREE);
    ETree SupETree = ETree_.ToSupernodalETree(Xsuper_,SupMembership_,Order_);
    //logfileptr->OFS()<<"SupETree:"<<SupETree<<std::endl;
    chSupTree_.resize(SupETree.Size()+1);
    for(Int I = 1; I<=SupETree.Size(); I++){
      Int parent = SupETree.PostParent(I-1);
      chSupTree_[parent].push_back(I);
    }
  }
#endif







  SYMPACK_TIMER_STOP(BUILD_TASK_LIST);

  SYMPACK_TIMER_STOP(FB_INIT);


  bool doPrint = true;
  Int prevCnt = -1;
  while(taskGraph.getTaskCount()>0){
    CheckIncomingMessages(taskGraph,aggVectors);

    if(!scheduler2_->done())
    {
      //Pick a ready task
      //auto taskit = scheduler_->top();
      //scheduler_->pop();
      //FBTask & curTask = *taskit;
      FBTask curTask = scheduler2_->top();
      scheduler2_->pop();


      //    assert(find(taskLists_[curTask.tgt_snode_id-1]->begin(),taskLists_[curTask.tgt_snode_id-1]->end(),curTask)!=taskLists_[curTask.tgt_snode_id-1]->end());
#ifdef _DEBUG_PROGRESS_
      logfileptr->OFS()<<"Processing T("<<curTask.src_snode_id<<","<<curTask.tgt_snode_id<<") "<<std::endl;
#endif
      Int iLocalTGT = snodeLocalIndex(curTask.tgt_snode_id);
      switch(curTask.type){
        case FACTOR:
          {
            FBFactorizationTask(taskGraph, curTask, iLocalTGT,aggVectors);
          }
          break;
        case AGGREGATE:
          {
            FBAggregationTask(taskGraph,curTask, iLocalTGT);
          }
          break;
        case UPDATE:
          {
            FBUpdateTask(taskGraph, curTask, UpdatesToDo, aggVectors);
          }
          break;
defaut:
          {
            abort();
          }
          break;
      }

      SYMPACK_TIMER_START(REMOVE_TASK);
      //remove task
      //taskLists[curTask.tgt_snode_id-1]->erase(taskit);
      //localTaskCount--;
      //taskGraph.removeTask(taskit);
      taskGraph.decreaseTaskCount();
      SYMPACK_TIMER_STOP(REMOVE_TASK);


    }

    prevCnt=taskGraph.getTaskCount();
  }

  SYMPACK_TIMER_START(BARRIER);
  upcxx::async_wait();
  MPI_Barrier(CommEnv_->MPI_GetComm());

  SYMPACK_TIMER_STOP(BARRIER);
#ifdef UPCXX_PROGRESS_THREAD
  upcxx::progress_thread_stop();
#endif

  tmpBufs.Clear();

  SYMPACK_TIMER_STOP(FACTORIZATION_FB);
}


template <typename T> inline void symPACKMatrix<T>::FBGetUpdateCount(std::vector<Int> & UpdatesToDo, std::vector<Int> & AggregatesToRecv,std::vector<Int> & LocalAggregates)
{
  SYMPACK_TIMER_START(FB_GET_UPDATE_COUNT);
  UpdatesToDo.resize(Xsuper_.size(),I_ZERO);
  AggregatesToRecv.resize(Xsuper_.size(),I_ZERO);
  LocalAggregates.resize(Xsuper_.size(),I_ZERO);
  std::vector<Int> marker(Xsuper_.size(),I_ZERO);

  //map of map of pairs (supno, count)
  std::map<Idx, std::map<Idx, Idx>  > Updates;
  std::map<Idx, std::map<Idx, Idx>  > sendAfter;

  //  std::vector<bool>isSent(Xsuper_.size()*np,false);
  //Int numLocSnode = ( (Xsuper_.size()-1) / np);
  //Int firstSnode = iam*numLocSnode + 1;

  Int numLocSnode = XsuperDist_[iam+1]-XsuperDist_[iam];
  Int firstSnode = XsuperDist_[iam];
  Int lastSnode = firstSnode + numLocSnode-1;
  bool is2D = dynamic_cast<Mapping2D*>(Mapping_)!=nullptr;

  for(Int locsupno = 1; locsupno<locXlindx_.size(); ++locsupno){
    Idx s = locsupno + firstSnode-1;

    Int first_col = Xsuper_[s-1];
    Int last_col = Xsuper_[s]-1;

    Ptr lfi = locXlindx_[locsupno-1];
    Ptr lli = locXlindx_[locsupno]-1;
    Idx prevSnode = -1;
    for(Ptr sidx = lfi; sidx<=lli;sidx++){
      Idx row = locLindx_[sidx-1];
      Int supno = SupMembership_[row-1];
      if(supno!=prevSnode){
        //Idx supno = locSupLindx_[sidx-1];
        if(marker[supno-1]!=s && supno!=s){
          marker[supno-1] = s;

          Int iFactorizer = -1;
          Int iUpdater = -1;
          if(is2D){
            iFactorizer = Mapping_->Map(Xsuper_[supno-1],supno);
            iUpdater = Mapping_->Map(Xsuper_[supno-1],s);
          }
          else{
            iFactorizer = Mapping_->Map(supno-1,supno-1);
            iUpdater = Mapping_->Map(supno-1,s-1);
          }

          if( iUpdater==iFactorizer){
            LocalAggregates[supno-1]++;
          }
          else{
            auto & procSendAfter = sendAfter[iUpdater];
            auto it = procSendAfter.find(supno);
            if(it==procSendAfter.end()){
              procSendAfter[supno]=s;
            }
            //            if(!isSent[(supno-1)*np+iUpdater]){
            //              AggregatesToRecv[supno-1]++;
            //              isSent[(supno-1)*np+iUpdater]=true;
            //            }
          }

          auto & procUpdates = Updates[iUpdater];
          procUpdates[supno]++;
          //          auto it = procUpdates.find(supno);
          //          if(it==procUpdates.end()){
          //            procUpdates[supno]=1;
          //          }
          //          else{
          //            it->second++;
          //          }

        }
      }
      prevSnode = supno;
    }
  }


  //for(auto itp = Updates.begin();itp!=Updates.end();itp++){
  //  logfileptr->OFS()<<"Updates on P"<<itp->first<<": ";
  //  for(auto it = itp->second.begin();it!=itp->second.end();it++){
  //    logfileptr->OFS()<<it->first<<" = "<<it->second<<" | ";
  //  }
  //  logfileptr->OFS()<<std::endl;
  //}

  //Build a Alltoallv communication for Updates
  {
    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(std::pair<Idx,Idx>), MPI_BYTE, &type );
    MPI_Type_commit(&type);

    //compute send sizes
    vector<int> ssizes(np,0);
    int numPairs = 0;
    for(auto itp = Updates.begin();itp!=Updates.end();itp++){
      ssizes[itp->first] = itp->second.size();
      numPairs+=itp->second.size();
    }

    //compute send displacements
    vector<int> sdispls(np+1,0);
    sdispls[0] = 0;
    std::partial_sum(&ssizes.front(),&ssizes.back(),&sdispls[1]);

    //Build the contiguous array of pairs
    vector<std::pair<Idx, Idx> > sendbuf;
    sendbuf.reserve(numPairs);
    for(auto itp = Updates.begin();itp!=Updates.end();itp++){
      sendbuf.insert(sendbuf.end(),itp->second.begin(),itp->second.end());
    }

    //logfileptr->OFS()<<"Sent sizes: "<<ssizes<<std::endl;
    //logfileptr->OFS()<<"Sent displs: "<<sdispls<<std::endl;
    //logfileptr->OFS()<<"Sent updates: ";
    //for(auto it = sendbuf.begin();it!=sendbuf.end();it++){
    //  logfileptr->OFS()<<it->first<<" = "<<it->second<<" | ";
    //}
    //logfileptr->OFS()<<std::endl;

    //gather receive sizes
    vector<int> rsizes(np,0);
    MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,CommEnv_->MPI_GetComm());
    //for(int p = 0; p<np;p++){
    //  MPI_Gather(&ssizes[p],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,p,CommEnv_->MPI_GetComm());
    //}

    //compute receive displacements
    vector<int> rdispls(np+1,0);
    rdispls[0] = 0;
    std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);

    //logfileptr->OFS()<<"Recv sizes: "<<rsizes<<std::endl;
    //logfileptr->OFS()<<"Recv displs: "<<rdispls<<std::endl;

    //Now do the alltoallv

    vector<std::pair<Idx, Idx> > recvbuf(rdispls.back());
    MPI_Alltoallv(&sendbuf[0],&ssizes[0],&sdispls[0],type,&recvbuf[0],&rsizes[0],&rdispls[0],type,CommEnv_->MPI_GetComm());

    //logfileptr->OFS()<<"Recv updates: ";
    //for(auto it = recvbuf.begin();it!=recvbuf.end();it++){
    //  logfileptr->OFS()<<it->first<<" = "<<it->second<<" | ";
    //}
    //logfileptr->OFS()<<std::endl;

    std::map<Idx, Idx> LocalUpdates;
    for(auto it = recvbuf.begin();it!=recvbuf.end();it++){
      auto it2 = LocalUpdates.find(it->first);
      if(it2!=LocalUpdates.end()){
        it2->second += it->second;
      }
      else{
        LocalUpdates[it->first] = it->second;
      }
    }

    MPI_Type_free(&type);

    //logfileptr->OFS()<<"Local updates: ";
    //for(auto it = LocalUpdates.begin();it!=LocalUpdates.end();it++){
    //  logfileptr->OFS()<<it->first<<" = "<<it->second<<" | ";
    //}
    //logfileptr->OFS()<<std::endl;

    for(auto it = LocalUpdates.begin();it!=LocalUpdates.end();it++){
      UpdatesToDo[it->first-1] = it->second;
    }
  }

  //Build a Alltoallv communication for sendAfter
  {
    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(std::pair<Idx,Idx>), MPI_BYTE, &type );
    MPI_Type_commit(&type);
    //compute send sizes
    vector<int> ssizes(np,0);
    int numPairs = 0;
    for(auto itp = sendAfter.begin();itp!=sendAfter.end();itp++){
      ssizes[itp->first] = itp->second.size();
      numPairs+=itp->second.size();
    }

    //compute send displacements
    vector<int> sdispls(np+1,0);
    sdispls[0] = 0;
    std::partial_sum(&ssizes.front(),&ssizes.back(),&sdispls[1]);

    //Build the contiguous array of pairs
    vector<std::pair<Idx, Idx> > sendbuf;
    sendbuf.reserve(numPairs);
    for(auto itp = sendAfter.begin();itp!=sendAfter.end();itp++){
      sendbuf.insert(sendbuf.end(),itp->second.begin(),itp->second.end());
    }

    //logfileptr->OFS()<<"Sent sizes: "<<ssizes<<std::endl;
    //logfileptr->OFS()<<"Sent displs: "<<sdispls<<std::endl;
    //logfileptr->OFS()<<"Sent updates: ";
    //for(auto it = sendbuf.begin();it!=sendbuf.end();it++){
    //  logfileptr->OFS()<<it->first<<" = "<<it->second<<" | ";
    //}
    //logfileptr->OFS()<<std::endl;

    //gather receive sizes
    vector<int> rsizes(np,0);
    MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,CommEnv_->MPI_GetComm());
    //for(int p = 0; p<np;p++){
    //  MPI_Gather(&ssizes[p],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,p,CommEnv_->MPI_GetComm());
    //}

    //compute receive displacements
    vector<int> rdispls(np+1,0);
    rdispls[0] = 0;
    std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);

    //logfileptr->OFS()<<"Recv sizes: "<<rsizes<<std::endl;
    //logfileptr->OFS()<<"Recv displs: "<<rdispls<<std::endl;

    //Now do the alltoallv
    vector<std::pair<Idx, Idx> > recvbuf(rdispls.back());
    MPI_Alltoallv(&sendbuf[0],&ssizes[0],&sdispls[0],type,&recvbuf[0],&rsizes[0],&rdispls[0],type,CommEnv_->MPI_GetComm());

    //logfileptr->OFS()<<"Recv sendAfter: ";
    //for(auto it = recvbuf.begin();it!=recvbuf.end();it++){
    //  logfileptr->OFS()<<it->first<<" = "<<it->second<<" | ";
    //}
    //logfileptr->OFS()<<std::endl;


    std::map<Idx, Idx> mLocalSendAfter;
    for(auto it = recvbuf.begin();it!=recvbuf.end();it++){
      auto it2 = mLocalSendAfter.find(it->first);
      if(it2!=mLocalSendAfter.end()){
        if(it2->second>it->second){
          it2->second = it->second;
        }
      }
      else{
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
    for(int p=0;p<np;p++){
      int start = rdispls[p];
      int end = rdispls[p+1];

      for(int idx = start; idx<end;idx++){
        sendAfter[p][recvbuf[idx].first] = recvbuf[idx].second;
      }
    }

    std::fill(marker.begin(),marker.end(),I_ZERO);
    for(Int locsupno = 1; locsupno<locXlindx_.size(); ++locsupno){
      Idx s = locsupno + firstSnode-1;

      Int first_col = Xsuper_[s-1];
      Int last_col = Xsuper_[s]-1;

      Ptr lfi = locXlindx_[locsupno-1];
      Ptr lli = locXlindx_[locsupno]-1;
      Idx prevSnode  =-1;
      for(Ptr sidx = lfi; sidx<=lli;sidx++){
        Idx row = locLindx_[sidx-1];
        Int supno = SupMembership_[row-1];
        //Idx supno = locSupLindx_[sidx-1];

        if(supno!=prevSnode){
          if(marker[supno-1]!=s && supno!=s){
            marker[supno-1] = s;

          Int iFactorizer = -1;
          Int iUpdater = -1;
          if(is2D){
            iFactorizer = Mapping_->Map(Xsuper_[supno-1],supno);
            iUpdater = Mapping_->Map(Xsuper_[supno-1],s);
          }
          else{
            iFactorizer = Mapping_->Map(supno-1,supno-1);
            iUpdater = Mapping_->Map(supno-1,s-1);
          }


            if( iUpdater!=iFactorizer){
              auto & procSendAfter = sendAfter[iUpdater];
              if(s==procSendAfter[supno]){
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

  //logfileptr->OFS()<<" REDUCED AggregatesToRecv: "<<AggregatesToRecv<<std::endl;
  //logfileptr->OFS()<<"UpdatesToDo: "<<UpdatesToDo<<std::endl;
  //logfileptr->OFS()<<"LocalAggregates: "<<LocalAggregates<<std::endl;
  SYMPACK_TIMER_STOP(FB_GET_UPDATE_COUNT);
}

template <typename T> inline void symPACKMatrix<T>::FBAggregationTask(supernodalTaskGraph<FBTask> & taskGraph, FBTask & curTask, Int iLocalI, bool is_static)
{
  SYMPACK_TIMER_START(FB_AGGREGATION_TASK);

#ifdef _DEBUG_PROGRESS_
  logfileptr->OFS()<<"Processing T_AGGREG("<<curTask.src_snode_id<<","<<curTask.tgt_snode_id<<") "<<std::endl;
#endif
  Int src_snode_id = curTask.src_snode_id;
  Int tgt_snode_id = curTask.tgt_snode_id;
  Int I = src_snode_id;
  SuperNode<T> * src_snode = LocalSupernodes_[iLocalI -1];
  Int src_first_col = src_snode->FirstCol();
  Int src_last_col = src_snode->LastCol();

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

    if(!msgPtr->IsLocal()){
      msgPtr->DeallocRemote();
    }
    char* dataPtr = msgPtr->GetLocalPtr().get();

    SuperNode<T> * dist_src_snode = CreateSuperNode(options_.decomposition,dataPtr,msgPtr->Size());

    dist_src_snode->InitIdxToBlk();

    src_snode->Aggregate(dist_src_snode);

    delete dist_src_snode;
    delete msgPtr;
  }


  //need to update dependencies of the factorization task
  //auto taskit = find_task(taskLists,tgt_snode_id,tgt_snode_id,FACTOR);
  auto taskit = taskGraph.find_task(tgt_snode_id,tgt_snode_id,FACTOR);
  taskit->remote_deps--;

  if(!is_static){
    if(taskit->remote_deps==0 && taskit->local_deps==0){
      scheduler2_->push(*taskit);    
      taskGraph.removeTask(taskit);
    }
  }
  SYMPACK_TIMER_STOP(FB_AGGREGATION_TASK);
}



template <typename T> inline void symPACKMatrix<T>::FBFactorizationTask(supernodalTaskGraph<FBTask> & taskGraph, FBTask & curTask, Int iLocalI, std::vector< SuperNode<T> * > & aggVectors, bool is_static)
{
  SYMPACK_TIMER_START(FB_FACTORIZATION_TASK);

  Int src_snode_id = curTask.src_snode_id;
  Int tgt_snode_id = curTask.tgt_snode_id;
  Int I = src_snode_id;
  SuperNode<T> * src_snode = LocalSupernodes_[iLocalI -1];
  Int src_first_col = src_snode->FirstCol();
  Int src_last_col = src_snode->LastCol();

#ifdef _DEBUG_PROGRESS_
  logfileptr->OFS()<<"Processing Supernode "<<I<<std::endl;
#endif

#ifdef _DEBUG_PROGRESS_
  logfileptr->OFS()<<"Factoring Supernode "<<I<<std::endl;
#endif

  SYMPACK_TIMER_START(FACTOR_PANEL);
  src_snode->Factorize(tmpBufs);
  SYMPACK_TIMER_STOP(FACTOR_PANEL);

  //Sending factors and update local tasks
  //Send my factor to my ancestors. 
  std::vector<char> is_factor_sent(np);
  SetValue(is_factor_sent,false);

  SnodeUpdate curUpdate;
  SYMPACK_TIMER_START(FIND_UPDATED_ANCESTORS);
  SYMPACK_TIMER_START(FIND_UPDATED_ANCESTORS_FACTORIZATION);
#if 0
#error "this code is not working"
  Int prevSnode = -1;
  //Int src_snode_id = src_snode->Id();
  for(Int blkidx=0; blkidx< src_snode->NZBlockCnt();blkidx++){
    NZBlockDesc & nzblk_desc = src_snode->GetNZBlockDesc(blkidx);
    Idx fr = nzblk_desc.GIndex;
    Idx lr = src_snode->NRows(blkidx) + fr;
    Idx src_first_row = fr;
    while(src_first_row<lr){
      Int tgt_snode_id = SupMembership_[src_first_row];

      if(prevSnode!=tgt_snode_id){
        Int iTarget = this->Mapping_->Map(tgt_snode_id-1,src_snode->Id()-1);
        if(iTarget != iam){
          if(!is_factor_sent[iTarget]){
            MsgMetadata meta;

            Int local_first_row = src_first_row - nzblk_desc.GIndex;
            Int nzblk_cnt = src_snode->NZBlockCnt() - blkidx;
            Int nzval_cnt_ = src_snode->Size()*(src_snode->NRowsBelowBlock(blkidx)-local_first_row);
            T* nzval_ptr = src_snode->GetNZval(nzblk_desc.Offset) + local_first_row*src_snode->Size();

            upcxx::global_ptr<char> sendPtr((char*)nzval_ptr);
            //the size of the message is the number of bytes between sendPtr and the address of nzblk_desc

            //Send factor 
            meta.src = src_snode->Id();
            meta.tgt = tgt_snode_id;
            meta.GIndex = src_first_row;

            char * last_byte_ptr = (char*)&nzblk_desc + sizeof(NZBlockDesc);
            size_t msgSize = last_byte_ptr - (char*)nzval_ptr;
            signal_data(sendPtr, msgSize, iTarget, meta);
            is_factor_sent[iTarget] = true;
          }
        }          
        else{
          //Update local tasks
          //find task corresponding to curUpdate
          auto taskit = taskGraph.find_task(src_snode->Id(),tgt_snode_id,UPDATE);
//#pragma omp atomic
          taskit->local_deps--;
          if(!is_static){
            if(taskit->remote_deps==0 && taskit->local_deps==0){
              scheduler2_->push(*taskit);    
              taskGraph.removeTask(taskit);
            }
          }
        }
        prevSnode = tgt_snode_id;
      }


      src_first_row+=Xsuper_[tgt_snode_id-1];
    }
  }
#else
  while(src_snode->FindNextUpdate(curUpdate,Xsuper_,SupMembership_)){ 
    Int iTarget = this->Mapping_->Map(curUpdate.tgt_snode_id-1,src_snode->Id()-1);

    if(iTarget != iam){
      if(!is_factor_sent[iTarget]){
        MsgMetadata meta;

        //if(curUpdate.tgt_snode_id==29 && curUpdate.src_snode_id==28){gdb_lock(0);}
        //TODO Replace all this by a Serialize function
        NZBlockDesc & nzblk_desc = src_snode->GetNZBlockDesc(curUpdate.blkidx);
        Int local_first_row = curUpdate.src_first_row - nzblk_desc.GIndex;
        Int nzblk_cnt = src_snode->NZBlockCnt() - curUpdate.blkidx;
        Int nzval_cnt_ = src_snode->Size()*(src_snode->NRowsBelowBlock(curUpdate.blkidx)-local_first_row);
        T* nzval_ptr = src_snode->GetNZval(nzblk_desc.Offset) + local_first_row*src_snode->Size();

        upcxx::global_ptr<char> sendPtr((char*)nzval_ptr);
        //the size of the message is the number of bytes between sendPtr and the address of nzblk_desc

        //Send factor 
        meta.src = curUpdate.src_snode_id;
        meta.tgt = curUpdate.tgt_snode_id;
        meta.GIndex = curUpdate.src_first_row;

        char * last_byte_ptr = (char*)&nzblk_desc + sizeof(NZBlockDesc);
        size_t msgSize = last_byte_ptr - (char*)nzval_ptr;
        signal_data(sendPtr, msgSize, iTarget, meta);
        is_factor_sent[iTarget] = true;
      }
    }
    else{
      //Update local tasks
      //find task corresponding to curUpdate

      //auto taskit = find_task(taskLists,curUpdate.src_snode_id,curUpdate.tgt_snode_id,UPDATE);
      auto taskit = taskGraph.find_task(curUpdate.src_snode_id,curUpdate.tgt_snode_id,UPDATE);
      taskit->local_deps--;
      if(!is_static){
#ifdef PREFETCH_STRUCTURE
        // if target of the update is remote, then get its structure
        Int iFactorizer = this->Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
        if(iFactorizer!=iam){
          bool needStructure =  taskit->data.size()>0?(taskit->data.front()->meta.GIndex!=-1):true;
          //is the aggregate already allocated ?
          needStructure = needStructure && aggVectors[curUpdate.tgt_snode_id-1]==NULL;

          if(taskit->remote_deps==0 && taskit->local_deps==0 && needStructure){

            IncomingMessage * msg_ptr = new IncomingMessage();
            MsgMetadata meta;
            meta.src = curUpdate.src_snode_id;
            meta.tgt = curUpdate.tgt_snode_id;
            meta.GIndex = -1;

            msg_ptr->meta = meta;

            upcxx::global_ptr<SuperNodeDesc> pRemote_ptr = std::get<0>(remoteFactors_[curUpdate.tgt_snode_id-1]);
            size_t pMsg_size = std::get<1>(remoteFactors_[curUpdate.tgt_snode_id-1])*sizeof(NZBlockDesc)+sizeof(SuperNodeDesc);
            msg_ptr->remote_ptr = upcxx::global_ptr<char>(pRemote_ptr);
            msg_ptr->msg_size = pMsg_size;

            //allocate receive buffer
            bool success = msg_ptr->AllocLocal();
            if(success){
              gIncomingRecvAsync.push_back( msg_ptr );
              msg_ptr->AsyncGet();
              taskit->remote_deps++;
            }
            else{
              abort();
              delete msg_ptr;
            }
          }
        }
#endif

        if(taskit->remote_deps==0 && taskit->local_deps==0){
          //compute cost 
          //taskit->update_rank();
          //scheduler_->push(taskit);    
          scheduler2_->push(*taskit);    
          taskGraph.removeTask(taskit);
        }
      }
    }
  }
#endif
  SYMPACK_TIMER_STOP(FIND_UPDATED_ANCESTORS_FACTORIZATION);
  SYMPACK_TIMER_STOP(FIND_UPDATED_ANCESTORS);



  SYMPACK_TIMER_STOP(FB_FACTORIZATION_TASK);
}

template <typename T> inline void symPACKMatrix<T>::FBUpdateTask(supernodalTaskGraph<FBTask> & taskGraph, FBTask & curTask, std::vector<Int> & UpdatesToDo, std::vector< SuperNode<T> * > & aggVectors, bool is_static)
{
  SYMPACK_TIMER_START(FB_UPDATE_TASK);
  Int src_snode_id = curTask.src_snode_id;
  Int tgt_snode_id = curTask.tgt_snode_id;

  src_snode_id = abs(src_snode_id);
  bool is_first_local = curTask.src_snode_id <0;
  curTask.src_snode_id = src_snode_id;

  SuperNode<T> * cur_src_snode; 


  Int iSrcOwner = this->Mapping_->Map(abs(curTask.src_snode_id)-1,abs(curTask.src_snode_id)-1);
  //if(iSrcOwner!=iam){gdb_lock();}

  {
    //AsyncComms::iterator it;
    //AsyncComms * cur_incomingRecv = incomingRecvFactArr[tgt_snode_id-1];


    IncomingMessage * structPtr = NULL;
    IncomingMessage * msgPtr = NULL;
    //Local or remote factor
    //we have only one local or one remote incoming aggregate

#ifdef PREFETCH_STRUCTURE
    {
      Int iTarget = this->Mapping_->Map(curTask.tgt_snode_id-1,curTask.tgt_snode_id-1);
      if(curTask.data.size()==0){
        cur_src_snode = snodeLocal(curTask.src_snode_id);
      }
      else if(curTask.data.size()==1){
        if(iTarget!=iam){
          cur_src_snode = snodeLocal(curTask.src_snode_id);
          auto msgit = curTask.data.begin();
          structPtr = *msgit;
          assert(structPtr->meta.GIndex==-1);
        }
        else{
          auto msgit = curTask.data.begin();
          msgPtr = *msgit;
          assert(msgPtr->IsDone());
          char* dataPtr = msgPtr->GetLocalPtr().get();
          cur_src_snode = CreateSuperNode(options_.decomposition,dataPtr,msgPtr->Size(),msgPtr->meta.GIndex);
          cur_src_snode->InitIdxToBlk();
        }
      }
      else if(curTask.data.size()>1){
        assert(iTarget!=iam);

        auto msgit = curTask.data.begin();
        assert(curTask.data.size()>1);
        structPtr = *msgit;
        assert(structPtr->meta.GIndex==-1);
        msgit++;

        msgPtr = *msgit;
        assert(msgPtr->IsDone());
        char* dataPtr = msgPtr->GetLocalPtr().get();
        cur_src_snode = CreateSuperNode(options_.decomposition,dataPtr,msgPtr->Size(),msgPtr->meta.GIndex);
        cur_src_snode->InitIdxToBlk();
      }
    }
#else
    if(curTask.data.size()==0){
      cur_src_snode = snodeLocal(curTask.src_snode_id);
    }
    else{
      auto msgit = curTask.data.begin();
      msgPtr = *msgit;
      assert(msgPtr->IsDone());
      char* dataPtr = msgPtr->GetLocalPtr().get();
      cur_src_snode = CreateSuperNode(options_.decomposition,dataPtr,msgPtr->Size(),msgPtr->meta.GIndex);
      cur_src_snode->InitIdxToBlk();
    }
#endif


    //Update everything src_snode_id own with that factor
    //Update the ancestors
    SnodeUpdate curUpdate;
    SYMPACK_TIMER_START(UPDATE_ANCESTORS);
    while(cur_src_snode->FindNextUpdate(curUpdate,Xsuper_,SupMembership_,iam==iSrcOwner)){

      //skip if this update is "lower"
      if(curUpdate.tgt_snode_id<curTask.tgt_snode_id){continue;}

      Int iUpdater = this->Mapping_->Map(curUpdate.tgt_snode_id-1,cur_src_snode->Id()-1);

      if(iUpdater == iam){
#ifdef _DEBUG_PROGRESS_
        logfileptr->OFS()<<"implicit Task: {"<<curUpdate.src_snode_id<<" -> "<<curUpdate.tgt_snode_id<<"}"<<std::endl;
        logfileptr->OFS()<<"Processing update from Supernode "<<curUpdate.src_snode_id<<" to Supernode "<<curUpdate.tgt_snode_id<<std::endl;
#endif

        SuperNode<T> * tgt_aggreg;

        Int iTarget = this->Mapping_->Map(curUpdate.tgt_snode_id-1,curUpdate.tgt_snode_id-1);
        if(iTarget == iam){
          //the aggregate std::vector is directly the target snode
          SYMPACK_TIMER_START(UPD_ANC_Agg_local);
          tgt_aggreg = snodeLocal(curUpdate.tgt_snode_id);
          assert(curUpdate.tgt_snode_id == tgt_aggreg->Id());
          SYMPACK_TIMER_STOP(UPD_ANC_Agg_local);
        }
        else{
          SYMPACK_TIMER_START(UPD_ANC_Agg_tmp);
          //Check if src_snode_id already have an aggregate std::vector
          if(/*AggregatesDone[curUpdate.tgt_snode_id-1]==0*/aggVectors[curUpdate.tgt_snode_id-1]==NULL){
            SYMPACK_TIMER_START(UPD_ANC_Agg_tmp_creat);
            //use number of rows below factor as initializer

            //TODO do a customized version for FANIN as we have all the factors locally
            // the idea is the following: do a DFS and stop each exploration at the first local descendant of current node
#ifdef FANIN_OPTIMIZATION
            if(options_.mappingTypeStr ==  "COL2D")
            {
              std::set<Idx> structure;
              std::list<Int> frontier;
              {
                scope_timer(a,MERGE_STRUCTURE_FANIN);
                Idx tgt_fc = Xsuper_[curUpdate.tgt_snode_id-1];
                dfs_traversal(chSupTree_,curUpdate.tgt_snode_id,frontier);
                //process frontier in decreasing order of nodes and merge their structure
                for(auto it = frontier.rbegin(); it!=frontier.rend(); it++){
                  Int I = *it;
                  SuperNode<T> * source = snodeLocal(I);
                  for(Int blkidx=0; blkidx< source->NZBlockCnt();blkidx++){
                    NZBlockDesc & nzblk_desc = source->GetNZBlockDesc(blkidx);
                    Idx fr = nzblk_desc.GIndex;
                    Idx lr = source->NRows(blkidx) + fr;
                    for(Idx row = fr; row<lr;row++){
                      if(row>=tgt_fc){
                        structure.insert(row);
                      }
                    }
                  }
                }
              }
              aggVectors[curUpdate.tgt_snode_id-1] = CreateSuperNode(options_.decomposition,curUpdate.tgt_snode_id, Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id]-1, iSize_,structure);
            } 
            else
#endif
            {
              std::set<Idx> structure;
              {
                scope_timer(a,FETCH_REMOTE_STRUCTURE);
#ifdef PREFETCH_STRUCTURE
                char* buffer = structPtr->GetLocalPtr().get();
                SuperNodeDesc * pdesc = (SuperNodeDesc*)buffer;
                NZBlockDesc * bufferBlocks = (NZBlockDesc*)(pdesc+1);
                Int block_cnt = std::get<1>(remoteFactors_[curUpdate.tgt_snode_id-1]);
                logfileptr->OFS()<<*pdesc<<std::endl<<block_cnt<<std::endl;
                assert(pdesc->blocks_cnt_==block_cnt);

                for(Int i =block_cnt-1;i>=0;i--){
                  NZBlockDesc & curdesc = bufferBlocks[i];
                  size_t end = (i>0)?bufferBlocks[i-1].Offset:pdesc->nzval_cnt_;
                  Int numRows = (end-curdesc.Offset)/pdesc->iSize_;

                  for(Idx row = 0; row<numRows;row++){
                    structure.insert(curdesc.GIndex+row);
                  }
                }
#else
                upcxx::global_ptr<SuperNodeDesc> remoteDesc = std::get<0>(remoteFactors_[curUpdate.tgt_snode_id-1]);
                Int block_cnt = std::get<1>(remoteFactors_[curUpdate.tgt_snode_id-1]);


                //allocate space to receive block descriptors
                char * buffer = (char*)UpcxxAllocator::allocate(sizeof(NZBlockDesc)*block_cnt+ sizeof(SuperNodeDesc));
                upcxx::global_ptr<char> remote = upcxx::global_ptr<char>(remoteDesc);
                upcxx::copy(remote, (char*)&buffer[0],block_cnt*sizeof(NZBlockDesc)+sizeof(SuperNodeDesc));
                SuperNodeDesc * pdesc = (SuperNodeDesc*)buffer;
                NZBlockDesc * bufferBlocks = (NZBlockDesc*)(pdesc+1);

                for(Int i =block_cnt-1;i>=0;i--){
                  NZBlockDesc & curdesc = bufferBlocks[i];
                  size_t end = (i>0)?bufferBlocks[i-1].Offset:pdesc->nzval_cnt_;
                  Int numRows = (end-curdesc.Offset)/pdesc->iSize_;

                  for(Idx row = 0; row<numRows;row++){
                    structure.insert(curdesc.GIndex+row);
                  }
                }
                UpcxxAllocator::deallocate((char*)buffer);
#endif
              }
              aggVectors[curUpdate.tgt_snode_id-1] = CreateSuperNode(options_.decomposition,curUpdate.tgt_snode_id, Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id-1], Xsuper_[curUpdate.tgt_snode_id]-1, iSize_,structure);
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
        tgt_aggreg->UpdateAggregate(cur_src_snode,curUpdate,tmpBufs,iTarget,iam);
        SYMPACK_TIMER_STOP(UPD_ANC_UPD);


        --UpdatesToDo[curUpdate.tgt_snode_id-1];
#ifdef _DEBUG_
        logfileptr->OFS()<<UpdatesToDo[curUpdate.tgt_snode_id-1]<<" updates left for Supernode "<<curUpdate.tgt_snode_id<<std::endl;
#endif

        //Send the aggregate if it's the last
        //If this is my last update sent it to curUpdate.tgt_snode_id
        SYMPACK_TIMER_START(UPD_ANC_Agg_Send);
        if(UpdatesToDo[curUpdate.tgt_snode_id-1]==0){
          if(iTarget != iam){
#ifdef _DEBUG_
            logfileptr->OFS()<<"Remote Supernode "<<curUpdate.tgt_snode_id<<" is updated by Supernode "<<cur_src_snode->Id()<<std::endl;
#endif

            tgt_aggreg->Shrink();

            MsgMetadata meta;


            //if(curUpdate.tgt_snode_id==86 && curUpdate.src_snode_id==84){gdb_lock();}

            NZBlockDesc & nzblk_desc = tgt_aggreg->GetNZBlockDesc(0);
            T* nzval_ptr = tgt_aggreg->GetNZval(0);

            //this is an aggregate
            meta.src = curUpdate.src_snode_id;
            meta.tgt = curUpdate.tgt_snode_id;
            meta.GIndex = nzblk_desc.GIndex;

            //            logfileptr->OFS()<<"Remote Supernode "<<curUpdate.tgt_snode_id<<" on P"<<iTarget<<" is updated by Supernode "<<cur_src_snode->Id()<<std::endl;

            upcxx::global_ptr<char> sendPtr(tgt_aggreg->GetStoragePtr(meta.GIndex));
            //the size of the message is the number of bytes between sendPtr and the address of nzblk_desc
            size_t msgSize = tgt_aggreg->StorageSize();
            signal_data(sendPtr, msgSize, iTarget, meta);

          }
        }
        SYMPACK_TIMER_STOP(UPD_ANC_Agg_Send);

        SYMPACK_TIMER_START(UPD_ANC_Upd_Deps);
        if(iTarget == iam)
        {
          //update the dependency of the factorization
          //auto taskit = find_task(taskLists,curUpdate.tgt_snode_id,curUpdate.tgt_snode_id,FACTOR);
          auto taskit = taskGraph.find_task(curUpdate.tgt_snode_id,curUpdate.tgt_snode_id,FACTOR);
          taskit->local_deps--;
          if(!is_static){
            if(taskit->remote_deps==0 && taskit->local_deps==0){
              //compute cost 
              //taskit->update_rank();
              //scheduler_->push(taskit);    
              scheduler2_->push(*taskit);    
              taskGraph.removeTask(taskit);
            }
          }
        }
        SYMPACK_TIMER_STOP(UPD_ANC_Upd_Deps);

        //if local update, push a new task in the queue and stop the while loop
        if(iam==iSrcOwner ){
          break;
        }

      }
    }
    SYMPACK_TIMER_STOP(UPDATE_ANCESTORS);


    if(msgPtr!=NULL){
      delete cur_src_snode;
      delete msgPtr;
    }

    if(structPtr!=NULL){
      delete structPtr;
    }

    //    if(curTask.data.size()>0){
    //      delete cur_src_snode;
    //      auto msgit = curTask.data.begin();
    //      IncomingMessage * msgPtr = *msgit;
    //      delete msgPtr;
    //    }





  }
  SYMPACK_TIMER_STOP(FB_UPDATE_TASK);
}


template <typename T> inline void symPACKMatrix<T>::CheckIncomingMessages(supernodalTaskGraph<FBTask> & taskGraph,std::vector< SuperNode<T> * > & aggVectors,bool is_static)
{
  scope_timer(a,CHECK_MESSAGE);
  //return;

  //call advance

#ifndef UPCXX_PROGRESS_THREAD
  SYMPACK_TIMER_START(UPCXX_ADVANCE);
  upcxx::advance();
  SYMPACK_TIMER_STOP(UPCXX_ADVANCE);
#endif

  bool comm_found = false;
  IncomingMessage * msg = NULL;
  FBTask * curTask = NULL;
  if(is_static){
    curTask = &scheduler2_->top();
#ifdef _DEBUG_PROGRESS_
    logfileptr->OFS()<<"Next Task T("<<curTask->src_snode_id<<","<<curTask->tgt_snode_id<<") "<<std::endl;
#endif
  }


  do{
    msg=NULL;

    SYMPACK_TIMER_START(MV_MSG_SYNC);
    //if we have some room, turn blocking comms into async comms
    if(gIncomingRecvAsync.size() < gMaxIrecv || gMaxIrecv==-1){
      if(is_static){
        while((gIncomingRecvAsync.size() < gMaxIrecv || gMaxIrecv==-1) && !gIncomingRecv.empty()){
          bool success = false;
          for(auto it = gIncomingRecv.begin();it!=gIncomingRecv.end();it++){
            //find a message corresponding to current task
            IncomingMessage * curMsg = *it;
            if(curMsg->meta.tgt==curTask->tgt_snode_id){
              success = (curMsg)->AllocLocal();
              if(success){
                (curMsg)->AsyncGet();
                gIncomingRecvAsync.push_back(curMsg);
                //gIncomingRecv.pop();
                gIncomingRecv.erase(it);
              }
              else{
                abort();
              }
              break;
            }
          }

          if(!success){
            break;
          }
        }
      }
      else{
        while((gIncomingRecvAsync.size() < gMaxIrecv || gMaxIrecv==-1) && !gIncomingRecv.empty()){
          bool success = false;
          //auto it = gIncomingRecv.top();
          auto it = gIncomingRecv.begin();
          //find one which is not done
          while((*it)->IsDone()){it++;}

          success = (*it)->AllocLocal();
          if(success){
            (*it)->AsyncGet();
            gIncomingRecvAsync.push_back(*it);

            gIncomingRecv.erase(it);
            //gIncomingRecv.pop_front();
#ifdef _DEBUG_PROGRESS_
            logfileptr->OFS()<<"TRANSFERRED TO ASYNC COMM"<<std::endl;
#endif
          }
          else{
            //TODO handle out of memory
            abort();
          }
          if(!success){
            break;
          }
        }
      }
    }
    SYMPACK_TIMER_STOP(MV_MSG_SYNC);

#ifdef HANDLE_LOCAL_POINTER
    //find if there is some local (SHMEM) comms
    if(!gIncomingRecvLocal.empty()){
      auto it = gIncomingRecvLocal.begin();
      msg = *it;
      gIncomingRecvLocal.erase(it);

#ifdef _DEBUG_PROGRESS_
      logfileptr->OFS()<<"COMM LOCAL: AVAILABLE MSG("<<msg->meta.src<<","<<msg->meta.tgt<<") from P"<<msg->remote_ptr.where()<<std::endl;
#endif
    }
    else
#endif
    {
      //find if there is some finished async comm
      auto it = TestAsyncIncomingMessage();
      if(it!=gIncomingRecvAsync.end()){

        SYMPACK_TIMER_START(RM_MSG_ASYNC);
        msg = *it;
        gIncomingRecvAsync.erase(it);
        SYMPACK_TIMER_STOP(RM_MSG_ASYNC);
#ifdef _DEBUG_PROGRESS_
        logfileptr->OFS()<<"COMM ASYNC: AVAILABLE MSG("<<msg->meta.src<<","<<msg->meta.tgt<<") from P"<<msg->remote_ptr.where()<<std::endl;
#endif
      }
      else if((is_static || scheduler2_->done()) && !gIncomingRecv.empty()){
        //find a "high priority" task

#if 0
        {
          auto tmp = gIncomingRecv;
          while(!tmp.empty()){
            auto it = tmp.top();
            msg = it;
            tmp.pop();
            logfileptr->OFS()<<"COMM: AVAILABLE MSG("<<msg->meta.src<<","<<msg->meta.tgt<<") from P"<<msg->remote_ptr.where()<<std::endl;
          }
        }
#endif

        SYMPACK_TIMER_START(RM_MSG_SYNC);
        if(is_static){
          msg=NULL;
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
        else{
#if 1
          //find a task we would like to process
          auto it = gIncomingRecv.begin();
          for(auto cur_msg = gIncomingRecv.begin(); 
              cur_msg!= gIncomingRecv.end(); cur_msg++){
            //look at the meta data
            if((*cur_msg)->meta.tgt < (*it)->meta.tgt){
              it = cur_msg;
            }
          }
          msg = *it;
          gIncomingRecv.erase(it);
#else
          auto it = gIncomingRecv.front();
          msg = it;
          gIncomingRecv.pop_front();
#endif
        }
        SYMPACK_TIMER_STOP(RM_MSG_SYNC);
      }

    }

    if(msg!=NULL){
      scope_timer(a,WAIT_AND_UPDATE_DEPS);
      bool success = msg->Wait(); 
      if(!success){
        //          //there has been a problem: only one case: sync get and memory allocation has failed
        //          //try to use the backup buffer
        //          assert(!backupBuffer_.InUse());
        //         
        //          //TODO handle backup mode : messages have to be consumed one by one: can't delay the processing 
        //          msg->SetLocalPtr(backupBuffer_.GetPtr(),false);
        //          success = msg->Wait(); 
        assert(success);
      }

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
        logfileptr->OFS()<<"COMM: FETCHED AGGREG MSG("<<msg->meta.src<<","<<msg->meta.tgt<<") from P"<<iUpdater<<std::endl;
#endif

        if(!is_static){
          //insert an aggregation task
          Int I = msg->meta.src;
          Int J = msg->meta.tgt;
          FBTask curUpdate;
          curUpdate.type=AGGREGATE;
          curUpdate.src_snode_id = I;
          curUpdate.tgt_snode_id = J;

          curUpdate.remote_deps = 1;
          curUpdate.local_deps = 0;

          taskit = taskGraph.addTask(curUpdate);
        }
        else{
          taskit = taskGraph.find_task(msg->meta.tgt,msg->meta.tgt,FACTOR);
        }
      }
      else{
        Int iOwner = Mapping_->Map(msg->meta.src-1,msg->meta.src-1);
#ifdef _DEBUG_PROGRESS_
        logfileptr->OFS()<<"COMM: FETCHED FACTOR MSG("<<msg->meta.src<<","<<msg->meta.tgt<<") from P"<<iOwner<<std::endl;
#endif

#ifdef DYN_TASK_CREATE
        {
          abort();
          FBTask curTask;
          curTask.type=UPDATE;
          curTask.src_snode_id = msg->meta.src;
          curTask.tgt_snode_id = msg->meta.tgt;

          //If I own the factor, it is a local dependency
          curTask.remote_deps = 1;
          curTask.local_deps = 0;

          //create the list if needed

          taskit = taskGraph.addTask(curTask);

          Int J = msg->meta.tgt;
          //            if(taskLists[J-1] == NULL){
          //              taskLists[J-1]=new std::list<FBTask>();
          //            }
          //            taskLists[J-1]->push_back(curTask);
          //            localTaskCount++;
          //            taskit = --taskLists[J-1]->end();

          //compute cost
#if 0
          if(1){
            char* dataPtr = msg->GetLocalPtr();
#ifdef _INDEFINITE_
            SuperNode<T> * dist_src_snode = SuperNodeInd<T>(dataPtr,msg->Size(),msg->meta.GIndex);
#else
            SuperNode<T> * dist_src_snode = SuperNode<T>(dataPtr,msg->Size(),msg->meta.GIndex);
#endif
            dist_src_snode->InitIdxToBlk();
            SnodeUpdate curUpdate;
            taskit->rank = 0.0;
            while(dist_src_snode->FindNextUpdate(curUpdate,Xsuper_,SupMembership_,false)){
              //skip if this update is "lower"
              if(curUpdate.tgt_snode_id<curTask.tgt_snode_id){continue;}
              Int iUpdater = this->Mapping_->Map(curUpdate.tgt_snode_id-1,curTask.src_snode_id-1);
              if(iUpdater == iam){
                Int tgt_snode_width = Xsuper_[curUpdate.tgt_snode_id] - Xsuper_[curUpdate.tgt_snode_id-1];
                taskit->rank += dist_src_snode->NRowsBelowBlock(0)*pow(tgt_snode_width,2.0);
              }
            }
            delete dist_src_snode;
            taskit->rank*=-1;
          }
#endif

        }
#else
        taskit = taskGraph.find_task(msg->meta.src,msg->meta.tgt,UPDATE);
#endif
      } 

      taskit->remote_deps--;
      if(msg->meta.GIndex==-1 && msg->meta.src!=msg->meta.tgt){
        taskit->data.push_front(msg);
      }
      else{
        taskit->data.push_back(msg);
      }

      //if this is a factor task, then we should delete the aggregate
      if(!msg->IsLocal()){
        if( taskit->type==FACTOR || taskit->type==AGGREGATE){
          msg->DeallocRemote();
        }
      }

      if(!is_static){
#ifdef PREFETCH_STRUCTURE
        // if target of the update is remote, then get its structure
        Int iFactorizer = this->Mapping_->Map(msg->meta.tgt-1,msg->meta.tgt-1);
        if(iFactorizer!=iam && taskit->type == UPDATE){
          bool needStructure =  taskit->data.size()>0?(taskit->data.front()->meta.GIndex!=-1):true;
          //is the aggregate already allocated ?
          needStructure = needStructure && aggVectors[msg->meta.tgt-1]==NULL;

          if(taskit->remote_deps==0 && taskit->local_deps==0 && needStructure){
            IncomingMessage * msg_ptr = new IncomingMessage();
            MsgMetadata meta;
            meta.src = msg->meta.src;
            meta.tgt = msg->meta.tgt;
            meta.GIndex = -1;

            msg_ptr->meta = meta;

            upcxx::global_ptr<SuperNodeDesc> pRemote_ptr = std::get<0>(remoteFactors_[meta.tgt-1]);
            size_t pMsg_size = std::get<1>(remoteFactors_[meta.tgt-1])*sizeof(NZBlockDesc)+sizeof(SuperNodeDesc);
            msg_ptr->remote_ptr = upcxx::global_ptr<char>(pRemote_ptr);
            msg_ptr->msg_size = pMsg_size;

            //allocate receive buffer
            bool success = msg_ptr->AllocLocal();
            if(success){
              gIncomingRecvAsync.push_back( msg_ptr );
              msg_ptr->AsyncGet();
              taskit->remote_deps++;
            }
            else{
              abort();
              delete msg_ptr;
            }
          }
        }
#endif
        if(taskit->remote_deps==0 && taskit->local_deps==0){
          //scheduler_->push(taskit);    
          scheduler2_->push(*taskit);    
          taskGraph.removeTask(taskit);
        }
      }
    }
#ifdef UPCXX_PROGRESS_THREAD
    else{
      upcxx::progress_thread_stop();
      upcxx::advance();
      upcxx::progress_thread_start();
    }
#endif

  }while(msg!=NULL);
}


#endif //_SYMPACK_MATRIX_IMPL_FB_PULL_HPP_
