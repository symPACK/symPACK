#ifndef _SYMPACK_MATRIX_IMPL_SOLVE_HPP_
#define _SYMPACK_MATRIX_IMPL_SOLVE_HPP_

#include <sympack/symPACKMatrix.hpp>

template <typename T> inline void symPACKMatrix<T>::solve_(T * RHS, int nrhs,  T * Xptr) {
  scope_timer(a,SPARSE_SOLVE_INTERNAL);
  Int n = this->iSize_;

  if(this->iam<this->np){

    this->remDealloc = new upcxx::dist_object<int>(0);

    Int nsuper = this->TotalSupernodeCnt();
    auto SupETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);

    double timeSta = get_time();
    //compute children size
    std::vector<Int> rem_children(nsuper,0);
    std::vector<Int> loc_children(nsuper,0);
    for(Int I=1;I<nsuper;I++){
      Int parent = SupETree.PostParent(I-1);


      if(parent!=0){
        Int childOwner = this->Mapping_->Map(I-1,I-1);
        Int parentOwner = this->Mapping_->Map(parent-1,parent-1);

        if(parentOwner == childOwner){
          ++loc_children[parent-1];
        }
        else{
          ++rem_children[parent-1];
        }
      }
    }
    double timeStop = get_time();
        if(this->iam==0 && this->options_.verbose){
          symPACKOS<<"Solve getting children count time: "<<timeStop - timeSta<<std::endl;
        }


        {
    timeSta = get_time();
    Contributions2_.resize(LocalSupernodes_.size());

    double timeAlloc = 0.0;
    double timeCopy = 0.0;
    Int nsupLocal = LocalSupernodes_.size();
    for(Int iLocalI = 1; iLocalI <= nsupLocal; iLocalI++){
      auto cur_snode = this->LocalSupernodes_[iLocalI-1];
      Int I = cur_snode->Id();
        //Create the blocks of my contrib with the same nz structure as L
        //and the same width as the final solution
        //MEMORY CONSUMPTION TOO HIGH ?
        timeAlloc -= get_time();
        Contributions2_[iLocalI-1].reset(CreateSuperNode<UpcxxAllocator>(this->options_.decomposition,I,0,1,nrhs, cur_snode->NRowsBelowBlock(0) ,this->iSize_, cur_snode->NZBlockCnt(),this->options_.panel ));
        auto contrib = std::dynamic_pointer_cast< SuperNode<T,UpcxxAllocator> >(Contributions2_[iLocalI-1]);
        timeAlloc += get_time();

        for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){
          NZBlockDesc & cur_desc = cur_snode->GetNZBlockDesc(blkidx);
          contrib->AddNZBlock(cur_snode->NRows(blkidx),cur_desc.GIndex);

          //Copy RHS into contrib if first block
          if(blkidx==0){

        timeCopy -= get_time();
            T * diag_nzval = contrib->GetNZval(0);
            for(Int kk = 0; kk<cur_snode->Size(); ++kk){
              //First, copy the RHS into the contribution
              Int srcRow = this->Order_.perm[cur_desc.GIndex+kk-1];
              for(Int j = 0; j<nrhs;++j){
                diag_nzval[kk*nrhs+j] = RHS[srcRow-1 + j*n];
              }
            }
        timeCopy += get_time();

          }
        }

        if (contrib->NZBlockCnt()>1){
          NZBlockDesc & cur_desc = contrib->GetNZBlockDesc(1);
          Int nRows = contrib->NRowsBelowBlock(1);
          std::fill(contrib->GetNZval(cur_desc.Offset),contrib->GetNZval(cur_desc.Offset)+nRows*nrhs,ZERO<T>());
        }
    }
    timeStop = get_time();
        if(this->iam==0 && this->options_.verbose){
          symPACKOS<<"Solve allocating contributions time / alloc / copy: "<<timeStop - timeSta<<" / "<<timeAlloc<<" / "<<timeCopy<<std::endl;
        }
  }

  scheduler_new_->msgHandle = nullptr;
  scheduler_new_->threadInitHandle_ = nullptr;
  scheduler_new_->extraTaskHandle_  = nullptr;

    taskGraph graph;


    std::hash<std::string> hash_fn;

    auto dec_ref = [&] ( taskGraph::task_iterator taskit, Int loc, Int rem) {
#ifdef SP_THREADS
        std::lock_guard<std::recursive_mutex> lock(scheduler_new_->list_mutex_);
#endif
      taskit->second->local_deps-= loc;
      taskit->second->remote_deps-= rem;
      if(taskit->second->remote_deps==0 && taskit->second->local_deps==0){
        scheduler_new_->push(taskit->second);
        graph.removeTask(taskit->second->id);
      }
    };

    auto log_task = [&] (taskGraph::task_iterator taskit){
        SparseTask * tmp = ((SparseTask*)taskit->second.get());
        Int * meta = reinterpret_cast<Int*>(tmp->meta.data());
        Solve::op_type & type = *reinterpret_cast<Solve::op_type*>(&meta[2]);
        std::string name;
        switch(type){
          case Solve::op_type::FUC:
            name="FUC";
            break;
          case Solve::op_type::BUC:
            name="BUC";
            break;
          case Solve::op_type::BU:
            name="BU";
            break;
          case Solve::op_type::FU:
            name="FU";
            break;
          default: break; // silence -Wswitch warnings from clang
        }
        logfileptr->OFS()<<"  updating "<<name<<" "<<meta[0]<<"_"<<meta[1]<<std::endl;
    };

    auto log_msg = [&] (Int src, Int tgt, Solve::op_type type){
        std::string name;
        switch(type){
          case Solve::op_type::FUC:
            name="FUC";
            break;
          case Solve::op_type::BUC:
            name="BUC";
            break;
          case Solve::op_type::BU:
            name="BU";
            break;
          case Solve::op_type::FU:
            name="FU";
            break;
          default: break; // silence -Wswitch warnings from clang
        }
        logfileptr->OFS()<<"  sending to "<<name<<" "<<src<<"_"<<tgt<<std::endl;
    };

    //TODO Only loop through local supernodes

    timeSta = get_time();
    Int maxNrhsPerTask = std::min(nrhs,nrhs);
    Int numSubTasks = std::ceil(nrhs / maxNrhsPerTask);

    //Do a bottom up traversal

    Int nsupLocal = LocalSupernodes_.size();
    for(Int iLocalI = 1; iLocalI <= nsupLocal; iLocalI++){
      auto cur_snode = this->LocalSupernodes_[iLocalI-1];
      Int I = cur_snode->Id();
      //Add the FUC task
      {
        Int iOwner = Mapping_->Map(I-1,I-1);
        if(this->iam==iOwner){
          for(Int task = 0; task < numSubTasks; task++){
            Int taskNrhs = (task==numSubTasks-1)?nrhs-(numSubTasks-1)*maxNrhsPerTask:maxNrhsPerTask;
            std::shared_ptr<GenericTask> pFUCtask(new SparseTask);
            SparseTask & FUCtask = *(SparseTask*)pFUCtask.get();

            FUCtask.meta.resize(3*sizeof(Int)+sizeof(Solve::op_type));
            Int * meta = reinterpret_cast<Int*>(FUCtask.meta.data());
            meta[0] = I;
            meta[1] = I;
            meta[2] = task;
            Solve::op_type & type = *reinterpret_cast<Solve::op_type*>(&meta[3]);
            type = Solve::op_type::FUC;

            FUCtask.local_deps = loc_children[I-1];
            FUCtask.remote_deps = rem_children[I-1];
            //if this is the last task, add all the other subtasks as dependencies
            Int parent = SupETree.PostParent(I-1);
            if(parent!=0){
              Int parentOwner = this->Mapping_->Map(parent-1,parent-1);
              if(parentOwner!=this->iam){
                if(numSubTasks>1 && task==numSubTasks-1){
                  FUCtask.local_deps += numSubTasks-1;
                }
              }
            }

            //Loop through my children to get the dependencies

            FUCtask.execute = [&,this,I,pFUCtask,task,maxNrhsPerTask,numSubTasks,taskNrhs](){

              scope_timer(a,SOLVE_FUC);

              Int src = I;
              Int tgt = I;
              Int nrhsOffset = task * maxNrhsPerTask;

              SparseTask & FUCtask = *(SparseTask*)pFUCtask.get();

              Int src_local = this->snodeLocalIndex(src); 
              auto cur_snode = this->LocalSupernodes_[src_local-1];
              auto contrib = std::dynamic_pointer_cast< SuperNode<T,UpcxxAllocator> >(this->Contributions2_[src_local-1]);

              //Apply all forward updates first
              //Do the remote ones
              for(auto && msgPtr : FUCtask.getData() ){
                  assert(msgPtr->IsDone());
                  char* dataPtr = msgPtr->GetLocalPtr().get();
                  auto dist_contrib = CreateSuperNode(this->options_.decomposition,dataPtr,msgPtr->Size());
                  Int iOwner = this->Mapping_->Map(dist_contrib->Id()-1,dist_contrib->Id()-1);
                  dist_contrib->InitIdxToBlk();
                  contrib->forward_update(&*dist_contrib,iOwner,this->iam,nrhsOffset,taskNrhs);
                  delete dist_contrib;

                  if(task==0){
                    //link this to the other subtasks
                    //decrement remote dependency count of those tasks
                    for(Int other = 1; other<numSubTasks;other++){
                      std::stringstream sstr;
                      sstr<<tgt<<"_"<<tgt<<"_"<<other<<"_"<<(Int)Solve::op_type::FUC;
                      auto id = hash_fn(sstr.str());
                      auto taskit = graph.find_task(id);

                      taskit->second->addData(msgPtr);
                      dec_ref(taskit,0,1);
                    }
                  }
              }

              //Do the local ones
              Int child_snode_id = src - 1;
              if(child_snode_id>0){
                Int children_found = 0;
                while(child_snode_id>0 && children_found<loc_children[src-1]+rem_children[src-1]){
                  Int parent = SupETree.PostParent(child_snode_id-1);
                  if(parent==src){
                    Int iTarget = this->Mapping_->Map(child_snode_id-1,child_snode_id-1);
                    if(iTarget==this->iam){
                      //do the local thing
                      Int src_local = snodeLocalIndex(child_snode_id); 
                      auto dist_contrib = std::dynamic_pointer_cast< SuperNode<T,UpcxxAllocator> >(Contributions2_[src_local-1]);
                      contrib->forward_update(&*dist_contrib,iTarget,this->iam,nrhsOffset,taskNrhs);
                    }
                    children_found++;
                  }
                  //last column of the prev supernode
                  child_snode_id--;
                }
              }

              contrib->forward_update_contrib(cur_snode,nrhsOffset,taskNrhs);

              //may have to send to parent
              
              Int parent = SupETree.PostParent(src-1);
              if(parent!=0){
                Int parentOwner = this->Mapping_->Map(parent-1,parent-1);
                if(parentOwner!=this->iam){
                  if(task==numSubTasks-1){
                    //Send to the parent if all subtasks are done // TODO
                    {
                      Int src_nzblk_idx = 1;
                      NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);
                      Int src_first_row = pivot_desc.GIndex;
                      T* nzval_ptr = contrib->GetNZval(pivot_desc.Offset);
                      upcxx::global_ptr<char> sendPtr = upcxx::to_global_ptr((char*)nzval_ptr);
                      MsgMetadata meta;
                      meta.src = src;
                      meta.tgt = parent;
                      meta.GIndex = src_first_row;

                      std::stringstream sstr;
                      sstr<<parent<<"_"<<parent<<"_"<<0<<"_"<<(Int)Solve::op_type::FUC;
                      meta.id = hash_fn(sstr.str());

                      char * last_byte_ptr = (char*)&pivot_desc + sizeof(NZBlockDesc);
                      size_t msgSize = last_byte_ptr - (char*)nzval_ptr;

                      {
        Int lparentOwner = this->group_->L2G(parentOwner);
                        signal_data(sendPtr, msgSize, lparentOwner, meta);
                      }
                    }
                  }
                  else{
                    //decrement the last FUC subtask
                    if(task<numSubTasks-1){
                      std::stringstream sstr;
                      sstr<<src<<"_"<<tgt<<"_"<<numSubTasks-1<<"_"<<(Int)Solve::op_type::FUC;
                      auto id = hash_fn(sstr.str());
                      auto taskit = graph.find_task(id);
                      dec_ref(taskit,1,0);
                    }
                  }
                }
                else{
                  //unlock the FUC task of the parent
                  std::stringstream sstr;
                  sstr<<parent<<"_"<<parent<<"_"<<task<<"_"<<(Int)Solve::op_type::FUC;
                  auto id = hash_fn(sstr.str());
                  auto taskit = graph.find_task(id);
                  dec_ref(taskit,1,0);
                }
              }
              else{
                //unlock the BUC task
                std::stringstream sstr;
                sstr<<tgt<<"_"<<tgt<<"_"<<task<<"_"<<(Int)Solve::op_type::BUC;
                auto id = hash_fn(sstr.str());
                auto taskit = graph.find_task(id);
                dec_ref(taskit,1,0);
              }
            };

            std::stringstream sstr;
            sstr<<I<<"_"<<I<<"_"<<task<<"_"<<(Int)type;
            FUCtask.id = hash_fn(sstr.str());

            graph.addTask(pFUCtask);
          }
        }
      }


      //Add the BUC task
      {
        Int iOwner = Mapping_->Map(I-1,I-1);
        if(this->iam==iOwner){
          for(Int task = 0; task < numSubTasks; task++){
            Int taskNrhs = (task==numSubTasks-1)?nrhs-(numSubTasks-1)*maxNrhsPerTask:maxNrhsPerTask;
            std::shared_ptr<GenericTask> pBUCtask(new SparseTask);
            SparseTask & BUCtask = *(SparseTask*)pBUCtask.get();
            BUCtask.meta.resize(3*sizeof(Int)+sizeof(Solve::op_type));
            Int * meta = reinterpret_cast<Int*>(BUCtask.meta.data());
            meta[0] = I;
            meta[1] = I;
            meta[2] = task;
            Solve::op_type & type = *reinterpret_cast<Solve::op_type*>(&meta[3]);
            type = Solve::op_type::BUC;

            BUCtask.local_deps = 1;
            //if this is the last task, add all the other subtasks as dependencies
            if(numSubTasks>1 && task==numSubTasks-1){
              BUCtask.local_deps += numSubTasks-1;
            }

            BUCtask.execute = [&,this,I,pBUCtask,task,maxNrhsPerTask,numSubTasks,taskNrhs] () {
              scope_timer(a,SOLVE_BUC);
              Int src = I;
              Int tgt = I;
              Int nrhsOffset = task * maxNrhsPerTask;
              SparseTask & BUCtask = *(SparseTask*)pBUCtask.get();

              Int * meta = reinterpret_cast<Int*>(BUCtask.meta.data());
              Solve::op_type & type = *reinterpret_cast<Solve::op_type*>(&meta[2]);
              bassert(src==meta[0]);
              bassert(tgt==meta[1]);
 
              Int src_local = snodeLocalIndex(src); 
              auto cur_snode = LocalSupernodes_[src_local-1];
              auto contrib = std::dynamic_pointer_cast< SuperNode<T,UpcxxAllocator> >(Contributions2_[src_local-1]);
              contrib->back_update_contrib(cur_snode,nrhsOffset,taskNrhs);

              //send to my children
              std::vector<char> is_sent(this->np,false);
              Int child_snode_id = src - 1;
              if(child_snode_id>0){
                Int children_found = 0;
                while(child_snode_id>0 && children_found<loc_children[src-1]+rem_children[src-1]){
                  Int parent = SupETree.PostParent(child_snode_id-1);
                  if(parent==src){
                    Int iTarget = this->Mapping_->Map(child_snode_id-1,child_snode_id-1);
                    if(iTarget!=this->iam){
                      if(task == numSubTasks-1){
                        if(!is_sent[iTarget]){

                          Int src_nzblk_idx = 0;
                          NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);
                          Int src_first_row = pivot_desc.GIndex;
                          T* nzval_ptr = contrib->GetNZval(pivot_desc.Offset);
                          upcxx::global_ptr<char> sendPtr = upcxx::to_global_ptr((char*)nzval_ptr);
                          MsgMetadata meta;
                          meta.src = src;
                          meta.tgt = child_snode_id;
                          meta.GIndex = src_first_row;

                          std::stringstream sstr;
                          sstr<<src<<"_"<<child_snode_id<<"_"<<0<<"_"<<(Int)Solve::op_type::BU;
                          meta.id = hash_fn(sstr.str());

                          char * last_byte_ptr = (char*)&pivot_desc + sizeof(NZBlockDesc);
                          size_t msgSize = last_byte_ptr - (char*)nzval_ptr;

                          {
        Int liTarget = this->group_->L2G(iTarget);
                            signal_data(sendPtr, msgSize, liTarget, meta);
                          }
                          is_sent[iTarget] = true;
                        }
                      }
                    }
                    else{
                      std::stringstream sstr;
                      sstr<<src<<"_"<<child_snode_id<<"_"<<task<<"_"<<(Int)Solve::op_type::BU;
                      auto id = hash_fn(sstr.str());
                      auto taskit = graph.find_task(id);
                      dec_ref(taskit,1,0);
                    }
                    children_found++;
                  }
                  //last column of the prev supernode
                  child_snode_id--;
                }
              }

              if(task < numSubTasks-1){
                std::stringstream sstr;
                sstr<<src<<"_"<<tgt<<"_"<<numSubTasks-1<<"_"<<(Int)Solve::op_type::BUC;
                auto id = hash_fn(sstr.str());
                auto taskit = graph.find_task(id);
                dec_ref(taskit,1,0);
              }


            };


            std::stringstream sstr;
            sstr<<I<<"_"<<I<<"_"<<task<<"_"<<(Int)type;
            BUCtask.id = hash_fn(sstr.str());

            graph.addTask(pBUCtask);
          }
        }
      }

      Int parent = SupETree.PostParent(I-1);
      if(parent!=0){
        //Add the BU task
        {
          Int iOwner = Mapping_->Map(I-1,I-1);
          if(this->iam==iOwner){
            for(Int task = 0; task < numSubTasks; task++){
              Int taskNrhs = (task==numSubTasks-1)?nrhs-(numSubTasks-1)*maxNrhsPerTask:maxNrhsPerTask;
              std::shared_ptr<GenericTask> pBUtask(new SparseTask);
              SparseTask & BUtask = *(SparseTask*)pBUtask.get();

              BUtask.meta.resize(3*sizeof(Int)+sizeof(Solve::op_type));
              Int * meta = reinterpret_cast<Int*>(BUtask.meta.data());
              meta[0] = parent;
              meta[1] = I;
              meta[2] = task;



              Solve::op_type & type = *reinterpret_cast<Solve::op_type*>(&meta[3]);
              type = Solve::op_type::BU;

              Int parentOwner = this->Mapping_->Map(parent-1,parent-1);
              if(parentOwner==this->iam){
                BUtask.local_deps = 1;
              }
              else{
                BUtask.remote_deps = 1;

                if(numSubTasks>1){
                  if(task==numSubTasks-1){
                    BUtask.local_deps += numSubTasks - 1;
                  }
                }
              }


              BUtask.execute = [&,this,I,parent,pBUtask,task,maxNrhsPerTask,numSubTasks,taskNrhs] () {
                scope_timer(a,SOLVE_BU);
                SparseTask & BUtask = *(SparseTask*)pBUtask.get();
                Int src = parent;
                Int tgt = I;
                Int nrhsOffset = task * maxNrhsPerTask;
                Int tgt_local = snodeLocalIndex(tgt); 
                auto contrib = std::dynamic_pointer_cast< SuperNode<T,UpcxxAllocator> >(Contributions2_[tgt_local-1]);

                Int iOwner = this->Mapping_->Map(src-1,src-1);
                if(iOwner==this->iam){
                  //local
                  Int src_local = snodeLocalIndex(src); 
                  auto dist_contrib = std::dynamic_pointer_cast< SuperNode<T,UpcxxAllocator> >(Contributions2_[src_local-1]);
                  contrib->back_update(&*dist_contrib,nrhsOffset,taskNrhs);
                }
                else{
                  //remote
                  assert(BUtask.getData().size()==1);
                  auto msgPtr = *BUtask.getData().begin();
                  assert(msgPtr->IsDone());
                  char* dataPtr = msgPtr->GetLocalPtr().get();
                  auto dist_contrib = CreateSuperNode(this->options_.decomposition,dataPtr,msgPtr->Size());
                  dist_contrib->InitIdxToBlk();
                  contrib->back_update(&*dist_contrib,nrhsOffset,taskNrhs);
                  delete dist_contrib;

                  if(task==0){
                    //link this to the other subtasks
                    //decrement remote dependency count of those tasks

                    for(Int other = 1; other<numSubTasks;other++){

                      std::stringstream sstr;
                      sstr<<src<<"_"<<tgt<<"_"<<other<<"_"<<(Int)Solve::op_type::BU;
                      auto id = hash_fn(sstr.str());
                      auto taskit = graph.find_task(id);

abort();
                      taskit->second->addData(msgPtr);
                      dec_ref(taskit,0,1);
                    }
                  }

                  //do I own another children ? if so, put msgPtr in its data and decrement remote_deps
                  bool child_found = false;
                  
                  if(task==numSubTasks-1){
                    Int parent = src;
                    if(parent!=0){

                      Int children_count = rem_children[parent-1];
                      Int child = src -1 ;
                      while(children_count>0 && child>0){
                        if(SupETree.PostParent(child-1)==src){
                          Int iTarget = this->Mapping_->Map(child-1,child-1);
                          if(iTarget==this->iam){
                            children_count--;

                            std::stringstream sstr;
                            sstr<<src<<"_"<<child<<"_"<<0<<"_"<<(Int)Solve::op_type::BU;
                            auto id = hash_fn(sstr.str());
                            auto taskit = graph.find_task(id);

                            if(child==tgt){
                              bassert(taskit==graph.tasks_.end());
                            }

                            if(taskit!=graph.tasks_.end()){
                              child_found = true;
                              logfileptr->OFS()<<"WARNING: this is suboptimal and should use the new ChainedMessage class."<<std::endl;
                              taskit->second->addData(msgPtr);
                              dec_ref(taskit,0,1);
                              break;
                            }
                          }
                        }
                        child--;
                      }
                    }
                  }


                  if(task<numSubTasks-1){
                    //decrement local dependency count of the last task
                    std::stringstream sstr;
                    sstr<<src<<"_"<<tgt<<"_"<<numSubTasks-1<<"_"<<(Int)Solve::op_type::BU;
                    auto id = hash_fn(sstr.str());
                    auto taskit = graph.find_task(id);
                    dec_ref(taskit,1,0);
                  }
                }

                std::stringstream sstr;
                sstr<<tgt<<"_"<<tgt<<"_"<<task<<"_"<<(Int)Solve::op_type::BUC;
                auto id = hash_fn(sstr.str());
                auto taskit = graph.find_task(id);
                dec_ref(taskit,1,0);
              };

              std::stringstream sstr;
              sstr<<parent<<"_"<<I<<"_"<<task<<"_"<<(Int)type;
              BUtask.id = hash_fn(sstr.str());

              graph.addTask(pBUtask);
            }
          }
        }
      }


    }

    timeStop = get_time();
        if(this->iam==0 && this->options_.verbose){
          symPACKOS<<"Solve task graph generation time: "<<timeStop - timeSta<<std::endl;
        }
    

    timeSta = get_time();
    scheduler_new_->run(CommEnv_->MPI_GetComm(),*this->group_,graph,*this->remDealloc,*this->workteam_);
    delete this->remDealloc;
    timeStop = get_time();
        if(this->iam==0 && this->options_.verbose){
          symPACKOS<<"Solve task graph execution time: "<<timeStop - timeSta<<std::endl;
        }

  }

  MPI_Barrier(CommEnv_->MPI_GetComm());




}


#endif // _SYMPACK_MATRIX_IMPL_SOLVE_HPP_
