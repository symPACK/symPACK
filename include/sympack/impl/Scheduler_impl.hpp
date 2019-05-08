#ifndef _SCHEDULER_IMPL_HPP_
#define _SCHEDULER_IMPL_HPP_

#include <sched.h>
#include <errno.h>
#include <unistd.h>
#include <pthread.h>
#include <thread>

//Definitions of the WorkQueue class
namespace symPACK{

  template<typename T, typename Queue >
    inline WorkQueue<T, Queue >::WorkQueue(Int nthreads): done(false){
#ifdef THREAD_VERBOSE
      processing_.resize(nthreads,nullptr);
#endif
      workQueues_.resize(nthreads);
      list_mutexes_.resize(nthreads);
      for (Int count {0}; count < nthreads; count += 1){
        threads.emplace_back(std::mem_fn<void(Int)>(&WorkQueue::consume ) , this, count);
      }

#if 0      
      for(int i = 0; i<threads.size(); i++){
        auto & thread = threads[i];
        // Create a cpu_set_t object representing a set of CPUs. Clear it and mark
        // only CPU i as set.
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(i, &cpuset);
        int rc = pthread_setaffinity_np(thread.native_handle(), sizeof(cpu_set_t), &cpuset);
        if (rc != 0) {
          std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
        }
      }
#endif
    }


  template<typename T, typename Queue >
    inline WorkQueue<T, Queue >::WorkQueue(Int nthreads, std::function<void()> & threadInitHandle ):done(false){
#ifdef THREAD_VERBOSE
      processing_.resize(nthreads,nullptr);
#endif
      workQueues_.resize(nthreads);
      list_mutexes_.resize(nthreads);

      threadInitHandle_ = threadInitHandle;
      for (Int count {0}; count < nthreads; count += 1){
        threads.emplace_back(std::mem_fn<void(Int)>(&WorkQueue::consume ) , this, count);
      }
#if 0
      for(int i = 0; i<threads.size(); i++){
        auto & thread = threads[i];
        // Create a cpu_set_t object representing a set of CPUs. Clear it and mark
        // only CPU i as set.
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(i, &cpuset);
        int rc = pthread_setaffinity_np(thread.native_handle(), sizeof(cpu_set_t), &cpuset);
        if (rc != 0) {
          std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
        }
      }
#endif
    }

  template<typename T, typename Queue >
    inline WorkQueue<T, Queue >::~WorkQueue(){
      if(threads.size()>0){
//        std::lock_guard<std::mutex> guard(list_mutex_);
        done = true;
//        sync.notify_all();
      }
      for (auto &&thread: threads) thread.join();
    }

  template<typename T, typename Queue >
    inline void WorkQueue<T, Queue >::pushTask(T & fut){
//      gdb_lock();
//      std::unique_lock<std::mutex> lock(list_mutex_);
//#ifdef THREAD_VERBOSE
//      auto it = std::find(workQueue_.begin(),workQueue_.end(),fut);
//      bassert(it==workQueue_.end());
//#endif


      auto queue_it = std::min_element(workQueues_.begin(),workQueues_.end(), [](Queue & a, Queue & b){return a.size()<b.size();});
      auto & queue = *queue_it;
      int queue_id = std::distance(workQueues_.begin(),queue_it);
      list_mutexes_[queue_id].lock();
      queue.push_back(fut);
      list_mutexes_[queue_id].unlock();
 //     sync.notify_one();
    }

  template<typename T, typename Queue >
    inline void WorkQueue<T, Queue >::consume(Int tid)
    {
#ifdef THREAD_VERBOSE
      std::stringstream sstr;
      sstr<<"Thread "<<tid<<std::endl;
      logfileptr->OFS()<<sstr.str();
#endif

      
      std::stringstream sstr2;
#if 0
      sstr2 << "Thread #" << tid << ": on CPU " << sched_getcpu() << "\n";
#endif
      logfileptr->OFS()<<sstr2.str();

      if(threadInitHandle_!=nullptr){
        threadInitHandle_();
      }

      auto & queue = workQueues_[tid];
      auto & list_mutex = list_mutexes_[tid];

      //std::unique_lock<std::mutex> lock(list_mutex_);
      while (true) {
        list_mutex.lock();
        bool empty = queue.empty();

        if ( ! empty) {
          T func { std::move(queue.front()) };
          queue.pop_front();
          list_mutex.unlock();


#ifdef THREAD_VERBOSE
          processing_[tid] = func;
#endif
          //sync.notify_one();
          bool success = false;
          while(!success){
            //lock.unlock();
            //try
            {
              func->execute();
              success=true;
              //clear resources
              func->reset();
              //lock.lock();
            }
//            catch(const MemoryAllocationException & e){
//              {
//                std::stringstream sstr;
//                sstr<<"Task locked on T"<<tid<<std::endl;
//                sstr<<e.what();
//                logfileptr->OFS()<<sstr.str();
//              }
//
//              lock.lock();
//              //wait for a task to be completed before retrying
//              //sync.wait(lock);
//              sync.wait_for(lock,std::chrono::milliseconds(10));
//              {
//                std::stringstream sstr;
//                sstr<<"Task un locked on T"<<tid<<std::endl;
//                logfileptr->OFS()<<sstr.str();
//              }
//
//              sync.notify_all();
//
//            }
          }

#ifdef THREAD_VERBOSE
          processing_[tid] = nullptr;
#endif
        } else if (done.load(std::memory_order_acquire) ){
          list_mutex.unlock();
          break;
        }
        else {
          list_mutex.unlock();
          sched_yield();
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
      static symPACK_scope_timer * persistent_timer = nullptr;

      scope_timer(a,CHECK_MESSAGE);

      if(persistent_timer==nullptr){persistent_timer = new symPACK_scope_timer("PERIOD_BETWEEN_MSG");}

      Int num_recv = 0;

#ifdef SP_THREADS
      if(Multithreading::NumThread>1){
        scope_timer(a,UPCXX_ADVANCE);
#ifdef NEW_UPCXX
        //nothing here this is handled by the progress thread
//        upcxx::progress();
#else
        //TODO this is obsolete and has to go`
        //std::lock_guard<upcxx_mutex_type> lock(upcxx_mutex);
        upcxx::advance();
#endif
      }
      else
#endif
      {
        scope_timer(a,UPCXX_ADVANCE);
#ifdef NEW_UPCXX
        upcxx::progress();
#else
        upcxx::advance();
#endif
      }

      bool comm_found = false;
      IncomingMessage * msg = nullptr;

      do{
        msg = nullptr;

        {
#ifdef SP_THREADS
        //  if(Multithreading::NumThread>1){
        //TODO this is obsolete and has to go`
        //    upcxx_mutex.lock();
        //  }
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
        //TODO this is obsolete and has to go`
        //  if(Multithreading::NumThread>1){
        //    upcxx_mutex.unlock();
        //  }
#endif
        }

        if(msg!=nullptr){
          scope_timer(a,WAIT_AND_UPDATE_DEPS);
          num_recv++;

          delete persistent_timer;
          persistent_timer=nullptr;

          bool success = false;
//          try{
        {
          scope_timer(a,MSG_WAIT_ONLY);
            success = msg->Wait(); 
        }
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
        //TODO this is obsolete and has to go`
              if(Multithreading::NumThread>1){ list_mutex_.lock(); }
#endif
              taskit->second->remote_deps--;
              taskit->second->addData(msgPtr);

              if(taskit->second->remote_deps==0 && taskit->second->local_deps==0){
                this->push(taskit->second);    
                graph.removeTask(taskit->second->id);
              }

#ifdef SP_THREADS
        //TODO this is obsolete and has to go`
              if(Multithreading::NumThread>1){ list_mutex_.unlock(); }
#endif
            }
          }
        }
      }while(msg!=nullptr);

#ifdef SP_THREADS
        //TODO this is obsolete and has to go`
      //if(Multithreading::NumThread>1){
      //  upcxx_mutex.lock();
      //}
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
        //TODO this is obsolete and has to go`
      //if(Multithreading::NumThread>1){
      //  upcxx_mutex.unlock();
      //}
#endif

      return num_recv;
    }


  template <class Task > 
#ifdef NEW_UPCXX
    inline void Scheduler<Task>::run(MPI_Comm & workcomm, RankGroup & group , taskGraph & graph, upcxx::dist_object<int> & remDealloc, upcxx::team & workteam)
#else
    inline void Scheduler<Task>::run(MPI_Comm & workcomm, RankGroup & group, taskGraph & graph)
#endif
    {
    }

  template <> 
#ifdef NEW_UPCXX
    inline void Scheduler<std::shared_ptr<GenericTask> >::run(MPI_Comm & workcomm, RankGroup & group , taskGraph & graph, upcxx::dist_object<int> & remDealloc, upcxx::team & workteam)
#else
    inline void Scheduler<std::shared_ptr<GenericTask> >::run(MPI_Comm & workcomm, RankGroup & group , taskGraph & graph)
#endif
    {
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
  //TODO it seems that this is dead code
//      auto check_handle = [] (std::future< void > const& f) {
//        bool retval = f.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready;
//        //logfileptr->OFS()<<"thread done ? "<<retval<<std::endl;
//        return retval; };
//
//      std::list<std::future<void> > handles;
//
//      if(Multithreading::NumThread>1){
//        logfileptr->OFS()<<"Num threads = "<<Multithreading::NumThread<<std::endl;
//      }
//
//      //handles.reserve(nthr);
//
//      auto progress_threads = [&]() {
//        //logfileptr->OFS()<<"handle size "<<handles.size()<<std::endl;
//        scope_timer(a,SOLVE_PROGRESS_THREADS);
//        auto it = handles.begin();
//        while(it != handles.end()){
//          //          it->wait();
//          //logfileptr->OFS()<<"checking handle"<<std::endl;
//          if(check_handle(*it)){
//            //logfileptr->OFS()<<"thread done"<<std::endl;
//            it->get();
//            it = handles.erase(it);
//          }
//          else{
//            it++;
//          }
//        }
//      };
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
        bool done = false;
        std::thread progress_thread;
        if(Multithreading::NumThread>1){
#ifdef NEW_UPCXX
          // declare an agreed upon persona for the progress thread
          //atomic<int> thread_barrier(0);
          // liberate the master persona to allow the progress thread to use it
          symPACK::liberate_master_scope();

          // create the progress thread
          progress_thread = std::thread( [this,&done,&remDealloc,&workteam]() {
              // push the master persona onto this thread's stack
              upcxx::persona_scope scope(upcxx::master_persona());
              // push progress_persona onto this thread's persona stack
              //upcxx::persona_scope scope(this->progress_persona_);
              // progress thread drains progress until work is done
              //gdb_lock();
              while (!done || (*remDealloc) > 0 ){
              sched_yield();
              upcxx::progress();
              }
              //        cout<<"Progress thread on process "<<upcxx::rank_me()<<" is done"<<endl; 

              upcxx::discharge();
//  upcxx::barrier(workteam);
//  while ( (*remDealloc) > 0 ) { sched_yield(); }

              //unlock the other threads
              //thread_barrier += 1;
          });

#endif
        }

        if(Multithreading::NumThread>1){
          WorkQueue<std::shared_ptr<GenericTask> > queue(Multithreading::NumThread,threadInitHandle_);

          while(graph.getTaskCount()>0 || !this->done() || !delayedTasks_.empty() ){

            if(np>1){
              num_msg = checkIncomingMessages_(graph);
            }

            if(extraTaskHandle_!=nullptr)
            {
              auto taskit = delayedTasks_.begin();

#ifdef THREAD_VERBOSE
              for(auto && toto: delayedTasks_){
                std::stringstream sstr;
                sstr<<"delayed";
                log_task(curTask,sstr);
                logfileptr->OFS()<<sstr.str();
              }
#endif

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
  //            std::lock_guard<std::mutex> lock(list_mutex_);
              if ( this->size() > 0 ) {

#ifdef SP_THREADS
      if(Multithreading::NumThread>1){ this->list_mutex_.lock(); }
#endif
                curTask = this->top();
                bassert(curTask==this->top());
                //if(extraTaskHandle_!=nullptr){
                //  if (curTask==nullptr){
                //    gdb_lock();
                //    for(auto && toto: delayedTasks_){
                //      bool delay = extraTaskHandle_(toto);
                //    }
                //  }
                //}


                bassert(curTask!=nullptr);
                this->pop();
#ifdef SP_THREADS
      if(Multithreading::NumThread>1){ this->list_mutex_.unlock(); }
#endif
                bool delay = false;
                if(extraTaskHandle_!=nullptr){
                  delay = extraTaskHandle_(curTask);
                  if(delay){
                    //gdb_lock();
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
              else{
                sched_yield();
              }
            }

#ifdef THREAD_VERBOSE
            {
              std::stringstream sstr;
//              queue.list_mutex_.lock();

              sstr<<"======delayed======"<<std::endl;
              for(auto ptr: delayedTasks_){ if(ptr!=nullptr){log_task(ptr,sstr);}}
              sstr<<"======waiting======"<<std::endl;
//              for(auto && ptr : queue.workQueue_){ log_task(ptr,sstr); }
              sstr<<"======running======"<<std::endl;
              for(auto ptr: queue.processing_){ if(ptr!=nullptr){log_task(ptr,sstr);}}
              sstr<<"==================="<<std::endl;
              logfileptr->OFS()<<sstr.str();
//              queue.list_mutex_.unlock();
            }
#endif
          }

#ifdef NEW_UPCXX
  //double tstop, tstart;
  //upcxx::progress(); //this is handled by the progress thread
  //upcxx::discharge();
  done = true; //this will stop the progress thread
  progress_thread.join();

    // recapture the master persona onto the initial thread's persona stack
    // before calling barrier
          symPACK::capture_master_scope();
  upcxx::barrier(workteam);
#endif
        }
        else{
          std::chrono::time_point<std::chrono::system_clock> start;
          while(graph.getTaskCount()>0 || !this->done()){
            if(np>1){
              num_msg = checkIncomingMessages_(graph);
            }

            while(!this->done()){
              //Pick a ready task
              curTask = this->top();
              this->pop();
              curTask->execute();
              //clear resources
              curTask->reset();
#ifdef NEW_UPCXX
              upcxx::progress();
#else
              upcxx::advance();
#endif
            }
          }
#ifdef NEW_UPCXX
  //double tstop, tstart;
  //upcxx::progress();
  upcxx::discharge();
  while ( (*remDealloc) > 0 ) { upcxx::progress(); }
  //while ( (*remDealloc) > 0 ) { sched_yield(); }
  //done = true; //this will stop the progress thread
  //progress_thread.join();
  //        symPACK::capture_master_scope();
  upcxx::barrier(workteam);
#endif
        }

#ifndef NEW_UPCXX
   double tstart = get_time();
   int barrier_id = get_barrier_id(np);
   signal_exit(barrier_id,np); 
   double tstop = get_time();
   logfileptr->OFS()<<"signal_exit time: "<<tstop-tstart<<std::endl;
   tstart = get_time();
   barrier_wait(barrier_id);
   tstop = get_time();
   logfileptr->OFS()<<"barrier wait: "<<tstop-tstart<<std::endl;
#endif
      }
    }


}// end namespace symPACK
//end of definitions of the Scheduler class


#endif // _SCHEDULER_IMPL_HPP_

