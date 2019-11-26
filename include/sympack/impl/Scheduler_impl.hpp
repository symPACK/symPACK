#ifndef _SCHEDULER_IMPL_HPP_
#define _SCHEDULER_IMPL_HPP_

#include <sched.h>
#include <errno.h>
#include <unistd.h>
#include <pthread.h>
#include <thread>

//Definitions of the WorkQueue class
namespace symPACK{

  template< typename ttask_t >
    class worker_ptrs {
      public: 
        using atomic_ptr = std::atomic<ttask_t *>; 
        atomic_ptr pointers[2];
        int wslot_;
        int rslot_;
        worker_ptrs():pointers{{nullptr},{nullptr}},wslot_(0),rslot_(0){
          bassert( nullptr== pointers[0].load(std::memory_order_acquire));
          bassert( nullptr== pointers[1].load(std::memory_order_acquire));
        }
        virtual ~worker_ptrs() {
          bassert( nullptr== pointers[0].load(std::memory_order_acquire));
          bassert( nullptr== pointers[1].load(std::memory_order_acquire));
        }

        bool set(ttask_t * t){
          ttask_t *  p = pointers[wslot_].load(std::memory_order_relaxed);
          if (nullptr != p) {
            return false;
          }
          pointers[wslot_].store(t, std::memory_order_release);
          wslot_^=1;
          return true;
        }
        
        ttask_t * get(){
          ttask_t * p = pointers[rslot_].load(std::memory_order_acquire);
          if ( p != nullptr ){
            pointers[rslot_].store(nullptr,std::memory_order_relaxed);
            rslot_^=1;
          }
          
          return p;
        }
    };

#ifdef SP_THREADS
  template<typename T, typename Queue >
    inline WorkQueue<T, Queue >::WorkQueue(Int nthreads): done(false){
      workQueues_.resize(nthreads);
      list_mutexes_.resize(nthreads);
      prev_slot_ = 0;
      for (Int count {0}; count < nthreads; count += 1){
        threads.emplace_back(std::mem_fn<void(Int)>(&WorkQueue::consume ) , this, count);
      }
    }

  template<typename T, typename Queue >
    inline WorkQueue<T, Queue >::WorkQueue(Int nthreads, std::function<void()> & threadInitHandle ):done(false){
      workQueues_.resize(nthreads);
      list_mutexes_.resize(nthreads);
      prev_slot_ = 0;

      threadInitHandle_ = threadInitHandle;
      for (Int count {0}; count < nthreads; count += 1){
        threads.emplace_back(std::mem_fn<void(Int)>(&WorkQueue::consume ) , this, count);
      }
    }

  template<typename T, typename Queue >
    inline WorkQueue<T, Queue >::~WorkQueue(){
      if(threads.size()>0){
        done = true;
      }
      for (auto &&thread: threads) thread.join();
    }

  template<typename T, typename Queue >
    inline void WorkQueue<T, Queue >::pushTask(T & fut){
      auto & queue = workQueues_[prev_slot_];
      int queue_id = prev_slot_;
      prev_slot_ = (prev_slot_+1)%workQueues_.size();
      std::lock_guard<std::mutex> lk(list_mutexes_[queue_id]);
      queue.push_back(fut);
    }

  template<typename T, typename Queue >
    inline void WorkQueue<T, Queue >::consume(Int tid) {
      if(threadInitHandle_!=nullptr){
        threadInitHandle_();
      }

      auto & queue = workQueues_[tid];
      auto & list_mutex = list_mutexes_[tid];

      while (true) {
        list_mutex.lock();
        bool empty = queue.empty();

        if ( ! empty) {
          T func { std::move(queue.front()) };
          queue.pop_front();
          list_mutex.unlock();


          bool success = false;
          while(!success){
            //try
            {
              func->execute();
              success=true;
              //clear resources
              func->reset();
            }
          }

        } else if (done.load(std::memory_order_acquire) ){
          list_mutex.unlock();
          break;
        }
        else {
          list_mutex.unlock();
          sched_yield();
        }
      }

      upcxx::discharge();
      upcxx::progress();
    }
#endif
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
      upcxx::progress();
      bool comm_found = false;
      IncomingMessage * msg = nullptr;

      do{
        msg = nullptr;

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
            msg = (*it);
            gIncomingRecv.erase(it);
          }
        }

        if(msg!=nullptr){
          scope_timer(a,WAIT_AND_UPDATE_DEPS);
          num_recv++;

          delete persistent_timer;
          persistent_timer=nullptr;

          bool success = false;
          {
            scope_timer(a,MSG_WAIT_ONLY);
            success = msg->Wait(); 
          }

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
            taskit->second->addData(msgPtr);
            if ( auto ptr = graph.updateTask(taskit,0,1) ) {
              this->push(ptr);
            }
          }
        }
      }while(msg!=nullptr);

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

      return num_recv;
    }


  template <class Task > 
    inline void Scheduler<Task>::run(MPI_Comm & workcomm, RankGroup & group , taskGraph & graph, upcxx::dist_object<int> & remDealloc, upcxx::team & workteam)
    {
    }

  template <> 
    inline void Scheduler<std::shared_ptr<GenericTask> >::run(MPI_Comm & workcomm, RankGroup & group , taskGraph & graph, upcxx::dist_object<int> & remDealloc, upcxx::team & workteam)
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


      {

        main_tid = std::this_thread::get_id();

        Int num_msg =0;
        std::shared_ptr<GenericTask> curTask = nullptr;
#ifdef SP_THREADS
        bool done = false;
        std::thread progress_thread;
        if(Multithreading::NumThread>1){
          // declare an agreed upon persona for the progress thread
          // liberate the master persona to allow the progress thread to use it
          symPACK::liberate_master_scope();

          // create the progress thread
          progress_thread = std::thread( [this,&done,&remDealloc,&workteam]() {
              // push the master persona onto this thread's stack
              upcxx::persona_scope scope(upcxx::master_persona());
              while (!done || (*remDealloc) > 0 ){
              sched_yield();
              upcxx::progress();
              }

              upcxx::discharge();
          });
        }

        if(Multithreading::NumThread>1){
          WorkQueue<std::shared_ptr<GenericTask> > queue(Multithreading::NumThread-2,threadInitHandle_);
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
              if ( this->size() > 0 ) {

#ifdef SP_THREADS
      if(Multithreading::NumThread>1){ this->list_mutex_.lock(); }
#endif
                curTask = this->top();
                bassert(curTask==this->top());
                bassert(curTask!=nullptr);
                this->pop();
#ifdef SP_THREADS
      if(Multithreading::NumThread>1){ this->list_mutex_.unlock(); }
#endif
                bool delay = false;
                if(extraTaskHandle_!=nullptr){
                  delay = extraTaskHandle_(curTask);
                  if(delay){
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

          }

            done = true; //this will stop the progress thread
            progress_thread.join();
            symPACK::capture_master_scope();
            upcxx::discharge();
            while ( (*remDealloc) > 0 ) { upcxx::progress(); }
            upcxx::barrier(workteam);
        }
        else
#endif
        {
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
              upcxx::progress();
            }
          }
            upcxx::discharge();
            while ( (*remDealloc) > 0 ) { upcxx::progress(); }
            upcxx::barrier(workteam);
        }
      }
    }


}// end namespace symPACK
//end of definitions of the Scheduler class


#endif // _SCHEDULER_IMPL_HPP_

