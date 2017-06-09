#ifndef _SCHEDULER_IMPL_HPP_
#define _SCHEDULER_IMPL_HPP_

//Definitions of the WorkQueue class
namespace symPACK{

  template<typename T, typename Queue >
    inline WorkQueue<T, Queue >::WorkQueue(Int nthreads){
#ifdef THREAD_VERBOSE
      processing_.resize(nthreads,nullptr);
#endif
      for (Int count {0}; count < nthreads; count += 1){
        threads.emplace_back(std::mem_fn<void(Int)>(&WorkQueue::consume ) , this, count);
      }
    }


  template<typename T, typename Queue >
    inline WorkQueue<T, Queue >::WorkQueue(Int nthreads, std::function<void()> & threadInitHandle ){
#ifdef THREAD_VERBOSE
      processing_.resize(nthreads,nullptr);
#endif
      threadInitHandle_ = threadInitHandle;
      for (Int count {0}; count < nthreads; count += 1)
        threads.emplace_back(std::mem_fn<void(Int)>(&WorkQueue::consume ) , this, count);
    }

  template<typename T, typename Queue >
    inline WorkQueue<T, Queue >::~WorkQueue(){
      if(threads.size()>0){
        std::lock_guard<std::mutex> guard(list_mutex_);
        done = true;
        sync.notify_all();
      }
      for (auto &&thread: threads) thread.join();
    }

  template<typename T, typename Queue >
    inline void WorkQueue<T, Queue >::pushTask(T & fut){
      std::unique_lock<std::mutex> lock(list_mutex_);
#ifdef THREAD_VERBOSE
      auto it = std::find(workQueue_.begin(),workQueue_.end(),fut);
      bassert(it==workQueue_.end());
#endif
      workQueue_.push_back(fut);
      sync.notify_one();
    }

  template<typename T, typename Queue >
    inline void WorkQueue<T, Queue >::consume(Int tid)
    {
#ifdef THREAD_VERBOSE
      std::stringstream sstr;
      sstr<<"Thread "<<tid<<std::endl;
      logfileptr->OFS()<<sstr.str();
#endif

      if(threadInitHandle_!=nullptr){
        threadInitHandle_();
      }

      std::unique_lock<std::mutex> lock(list_mutex_);
      while (true) {
        if (not workQueue_.empty()) {
          T func { std::move(workQueue_.front()) };
          workQueue_.pop_front();


#ifdef THREAD_VERBOSE
          processing_[tid] = func;
#endif
          sync.notify_one();
          bool success = false;
          while(!success){
            lock.unlock();
            try{
              func->execute();
              success=true;
              //clear resources
              func->reset();
              lock.lock();
            }
            catch(const MemoryAllocationException & e){
              {
                std::stringstream sstr;
                sstr<<"Task locked on T"<<tid<<std::endl;
                sstr<<e.what();
                logfileptr->OFS()<<sstr.str();
              }

              lock.lock();
              //gdb_lock();
              //wait for a task to be completed before retrying
              //sync.wait(lock);
              sync.wait_for(lock,std::chrono::milliseconds(10));
              {
                std::stringstream sstr;
                sstr<<"Task un locked on T"<<tid<<std::endl;
                logfileptr->OFS()<<sstr.str();
              }

              sync.notify_all();

            }
          }

#ifdef THREAD_VERBOSE
          processing_[tid] = nullptr;
#endif
        } else if (done) {
          break;
        } else {
          sync.wait(lock);
        }
      }
    }
}//end namespace symPACK
//end of definitions of the WorkQueue class


//Definitions of the Scheduler class
namespace symPACK{

  template< typename Task> 
    inline Int Scheduler<Task>::checkIncomingMessages_(taskGraph & graph)
    {
      abort();
      return 0;
    }

  template<>
    inline Int Scheduler<std::shared_ptr<GenericTask> >::checkIncomingMessages_(taskGraph & graph)
    {
      scope_timer(a,CHECK_MESSAGE);

      Int num_recv = 0;

#ifdef SP_THREADS
      if(Multithreading::NumThread>1){
        scope_timer(a,UPCXX_ADVANCE);
        std::lock_guard<upcxx_mutex_type> lock(upcxx_mutex);
        upcxx::advance();
      }
      else
#endif
      {
        scope_timer(a,UPCXX_ADVANCE);
        upcxx::advance();
      }

      bool comm_found = false;
      IncomingMessage * msg = nullptr;

      do{
        msg = nullptr;

        {
#ifdef SP_THREADS
          if(Multithreading::NumThread>1){
            upcxx_mutex.lock();
          }
#endif

#if 1
          {
            //find if there is some finished async comm
            auto it = TestAsyncIncomingMessage();
            if(it!=gIncomingRecvAsync.end()){
              scope_timer(b,RM_MSG_ASYNC);
              msg = (*it);
              gIncomingRecvAsync.erase(it);
            }
            else if(this->done() && !gIncomingRecv.empty()){
              scope_timer(c,RM_MSG_ASYNC);
              //find a "high priority" task
              //find a task we would like to process
              auto it = gIncomingRecv.begin();
              //for(auto cur_msg = gIncomingRecv.begin(); 
              //    cur_msg!= gIncomingRecv.end(); cur_msg++){
              //  //look at the meta data
              //  //TODO check if we can parametrize that
              //  if((*cur_msg)->meta.tgt < (*it)->meta.tgt){
              //    it = cur_msg;
              //  }
              //}
              msg = (*it);
              gIncomingRecv.erase(it);
            }
          }
#else
#endif

#ifdef SP_THREADS
          if(Multithreading::NumThread>1){
            upcxx_mutex.unlock();
          }
#endif
        }

        if(msg!=nullptr){
          scope_timer(a,WAIT_AND_UPDATE_DEPS);
          num_recv++;

          bool success = false;
//          try{
            success = msg->Wait(); 
//          }
//          catch(const MemoryAllocationException & e){
//            //put the message back in the blocking message queue
//#ifdef SP_THREADS
//            if(Multithreading::NumThread>1){
//              upcxx_mutex.lock();
//            }
//#endif
//            //gdb_lock();
//
//            gIncomingRecv.push_back(msg);
//#ifdef SP_THREADS
//            if(Multithreading::NumThread>1){
//              upcxx_mutex.unlock();
//            }
//#endif
//
//
//            msg = nullptr;
//          }

          if(msg!=nullptr){
            std::shared_ptr<IncomingMessage> msgPtr(msg);

            //TODO what are the reasons of failure ?
            bassert(success);

            if(msgHandle!=nullptr){
              scope_timer(a,WAIT_USER_MSG_HANDLE);
              //call user handle
              msgHandle(msgPtr);
            }

            auto taskit = graph.find_task(msg->meta.id);
            bassert(taskit!=graph.tasks_.end());
            {
#ifdef SP_THREADS
              if(Multithreading::NumThread>1){
                list_mutex_.lock();
              }
#endif
              taskit->second->remote_deps--;
              taskit->second->addData(msgPtr);

              if(taskit->second->remote_deps==0 && taskit->second->local_deps==0){
                this->push(taskit->second);    
                graph.removeTask(taskit->second->id);
              }

#ifdef SP_THREADS
              if(Multithreading::NumThread>1){
                list_mutex_.unlock();
              }
#endif
            }
          }
        }
      }while(msg!=nullptr);

#ifdef SP_THREADS
      if(Multithreading::NumThread>1){
        upcxx_mutex.lock();
      }
#endif
      //if we have some room, turn blocking comms into async comms
      if(gIncomingRecvAsync.size() < gMaxIrecv || gMaxIrecv==-1){
        scope_timer(b,MV_MSG_SYNC);
        while((gIncomingRecvAsync.size() < gMaxIrecv || gMaxIrecv==-1) && !gIncomingRecv.empty()){
          bool success = false;
          auto it = gIncomingRecv.begin();
          //find one which is not done
          while((*it)->IsDone()){it++;}

          success = (*it)->AllocLocal();
          if(success){
            (*it)->AsyncGet();
            gIncomingRecvAsync.push_back(*it);
            gIncomingRecv.erase(it);
          }
          else{
            //TODO handle out of memory
            abort();
            break;
          }
        }
      }
#ifdef SP_THREADS
      if(Multithreading::NumThread>1){
        upcxx_mutex.unlock();
      }
#endif

      return num_recv;
    }


  template <class Task > 
    inline void Scheduler<Task>::run(MPI_Comm & workcomm, taskGraph & graph){
    }

  template <> 
    inline void Scheduler<std::shared_ptr<GenericTask> >::run(MPI_Comm & workcomm,taskGraph & graph){

      int np = 1;
      MPI_Comm_size(workcomm,&np);
      int iam = 0;
      MPI_Comm_rank(workcomm,&iam);

      //put rdy tasks in rdy queue
      {
        auto taskit = graph.tasks_.begin();
        while (taskit != graph.tasks_.end()) {
          if(taskit->second->remote_deps==0 && taskit->second->local_deps==0){
            auto it = taskit;
            this->push(it->second);
            taskit++;
            graph.removeTask(it->second->id);
          }
          else{
            taskit++;
          }
        }
      }

      auto log_task = [&] (std::shared_ptr<GenericTask> & task){
        SparseTask * tmp = ((SparseTask*)task.get());
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
        }
        logfileptr->OFS()<<name<<" "<<meta[0]<<"_"<<meta[1]<<std::endl;
      };

#ifdef SP_THREADS
      auto check_handle = [] (std::future< void > const& f) {
        bool retval = f.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready;
        //logfileptr->OFS()<<"thread done ? "<<retval<<std::endl;
        return retval; };

      std::list<std::future<void> > handles;

      if(Multithreading::NumThread>1){
        logfileptr->OFS()<<"Num threads = "<<Multithreading::NumThread<<std::endl;
      }

      //handles.reserve(nthr);

      auto progress_threads = [&]() {
        //logfileptr->OFS()<<"handle size "<<handles.size()<<std::endl;
        scope_timer(a,SOLVE_PROGRESS_THREADS);
        auto it = handles.begin();
        while(it != handles.end()){
          //          it->wait();
          //logfileptr->OFS()<<"checking handle"<<std::endl;
          if(check_handle(*it)){
            //logfileptr->OFS()<<"thread done"<<std::endl;
            it->get();
            it = handles.erase(it);
          }
          else{
            it++;
          }
        }
      };
#endif

      {

        auto log_task = [&] (std::shared_ptr<GenericTask> & taskptr, std::stringstream & sstr){
          SparseTask * tmp = ((SparseTask*)taskptr.get());
          std::string name;
          Int * meta = reinterpret_cast<Int*>(tmp->meta.data());

          if(tmp->meta.size()==(2*sizeof(Int)+sizeof(Factorization::op_type))){
            Factorization::op_type & type = *reinterpret_cast<Factorization::op_type*>(&meta[2]);
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
          }
          else{
          }

          sstr<<"           T "<<name<<" "<<meta[0]<<"_"<<meta[1]<<std::endl;
        };



        Int num_msg =0;
        std::shared_ptr<GenericTask> curTask = nullptr;
        if(Multithreading::NumThread>1){
          WorkQueue<std::shared_ptr<GenericTask> > queue(Multithreading::NumThread,threadInitHandle_);

          while(graph.getTaskCount()>0 || !this->done() || !delayedTasks_.empty() ){

            if(np>1){
              num_msg = checkIncomingMessages_(graph);
            }

            if(extraTaskHandle_!=nullptr)
            {
              auto taskit = delayedTasks_.begin();
              while(taskit!=delayedTasks_.end()){
                bool delay = extraTaskHandle_(*taskit);
                if(!delay){
#ifdef THREAD_VERBOSE
                  std::stringstream sstr;
                  sstr<<"Resuming";
                  log_task(curTask,sstr);
                  logfileptr->OFS()<<sstr.str();
#endif
                  queue.pushTask(*taskit);
                  taskit = delayedTasks_.erase(taskit);
                }
                else{
                  taskit++;
                }
              }
            }

            while(!this->done()){
              curTask = nullptr;
              std::lock_guard<std::mutex> lock(list_mutex_);
              curTask = this->top();
              bassert(curTask!=nullptr);
              this->pop();
              bool delay = false;
              if(extraTaskHandle_!=nullptr){
                delay = extraTaskHandle_(curTask);
                if(delay){
#ifdef THREAD_VERBOSE
                  std::stringstream sstr;
                  sstr<<"Delaying";
                  log_task(curTask,sstr);
                  logfileptr->OFS()<<sstr.str();
#endif
                  delayedTasks_.push_back(curTask);
                }
              }

              if(!delay){
                queue.pushTask(curTask);
              }
            }

#ifdef THREAD_VERBOSE
            {
              std::stringstream sstr;
              queue.list_mutex_.lock();

              sstr<<"======delayed======"<<std::endl;
              for(auto ptr: delayedTasks_){ if(ptr!=nullptr){log_task(ptr,sstr);}}
              sstr<<"======waiting======"<<std::endl;
              for(auto && ptr : queue.workQueue_){ log_task(ptr,sstr); }
              sstr<<"======running======"<<std::endl;
              for(auto ptr: queue.processing_){ if(ptr!=nullptr){log_task(ptr,sstr);}}
              sstr<<"==================="<<std::endl;
              logfileptr->OFS()<<sstr.str();
              queue.list_mutex_.unlock();
            }
#endif
          }
        }
        else{
          std::chrono::time_point<std::chrono::system_clock> start;
          while(graph.getTaskCount()>0 || !this->done()){
            if(np>1){
              num_msg = checkIncomingMessages_(graph);
            }

            //            start = std::chrono::system_clock::now();
            //            while(!this->done() && (std::chrono::system_clock::now() - start < std::chrono::milliseconds(1))   )
            while(!this->done()){// && (std::chrono::system_clock::now() - start < std::chrono::milliseconds(1))   )
              //Pick a ready task
              curTask = this->top();
              this->pop();
              curTask->execute();
              //clear resources
              curTask->reset();
              upcxx::advance();
            }
          }
        }

        upcxx::async_wait();
        MPI_Barrier(workcomm);


      }
    }

}// end namespace symPACK
//end of definitions of the Scheduler class


#endif // _SCHEDULER_IMPL_HPP_

