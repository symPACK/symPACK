#ifndef _SCHEDULER_DECL_HPP_
#define _SCHEDULER_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/Types.hpp"
#include "sympack/Task.hpp"
#include "sympack/CommPull.hpp"
#include <list>
#include <queue>



#ifdef SP_THREADS
#include <future>
#include <mutex>
#endif
#include <chrono>

namespace symPACK{

#ifdef SP_THREADS

  template<typename T, typename Queue = std::list<T>  >
    class WorkQueue{
      public:
        std::mutex list_mutex_;
        Queue workQueue_;


        int prev_slot_;


        std::vector<Queue> workQueues_;
        std::deque<std::mutex> list_mutexes_;
        std::vector<std::thread> threads;
        std::condition_variable sync;
        std::function<void()> threadInitHandle_;
        std::atomic<bool> done;
#ifdef THREAD_VERBOSE
        std::vector<T> processing_;  
#endif

        WorkQueue(Int nthreads);
        WorkQueue(Int nthreads, std::function<void()> & threadInitHandle );
        WorkQueue(WorkQueue &&) = default;
        WorkQueue &operator=(WorkQueue &&) = delete;


        ~WorkQueue();

        void pushTask(T & fut);

      protected:
        void consume(Int tid);
    };
#endif

  template <class T >
    class Scheduler{
      public:
        Scheduler(){
        }
        virtual ~Scheduler() = default;

        virtual T& top() =0;
        virtual void pop() =0;

        virtual void push( const T& val) =0;
        virtual unsigned int size() =0;

        virtual bool done() =0;

        std::thread::id main_tid;
        std::list<T> delayedTasks_;
        std::function<void()> threadInitHandle_;
        std::function<bool(T&)> extraTaskHandle_;
        std::function<void( std::shared_ptr<IncomingMessage> )> msgHandle;
        Int checkIncomingMessages_(taskGraph & graph);
        void run(MPI_Comm & workcomm, RankGroup & group , taskGraph & graph, upcxx::dist_object<int> & remDealloc, upcxx::team & workteam);

#ifdef SP_THREADS
        std::recursive_mutex list_mutex_;
#endif
      protected:
    };




  template <class T >
    class DLScheduler: public Scheduler<T>{
      public:
        DLScheduler(){
        }

        virtual T& top(){ 
#ifdef SP_THREADS
        std::lock_guard<std::recursive_mutex> lk(this->list_mutex_);
#endif
          T & ref = const_cast<T&>(readyTasks_.top());
          return ref;
        }
        virtual void pop() {
#ifdef SP_THREADS
        std::lock_guard<std::recursive_mutex> lk(this->list_mutex_);
#endif
          readyTasks_.pop();
        }

        virtual void push( const T& val) { 
#ifdef SP_THREADS
        std::lock_guard<std::recursive_mutex> lk(this->list_mutex_);
#endif
          readyTasks_.push(val);
        }
        virtual unsigned int size() {
          int count = readyTasks_.size();
          return count;
        }

        virtual bool done(){ 
          bool empty = readyTasks_.empty();
          return empty;
        }



        const std::priority_queue<T, std::vector<T>, FBTaskCompare<T> > & GetQueue() {return  readyTasks_;}
      protected:
        std::priority_queue<T, std::vector<T>, FBTaskCompare<T> > readyTasks_;
    };


  template <class T >
    class MCTScheduler: public Scheduler<T>{
      public:
        MCTScheduler(){
        }

        virtual T& top(){
          T & ref = const_cast<T&>(readyTasks_.top());
          return ref;
        }
        virtual void pop() { 
          readyTasks_.pop();
        }

        virtual void push( const T& val) { 
          readyTasks_.push(val);
        }



        virtual unsigned int size() {
          int count = readyTasks_.size();
          return count;
        }

        virtual bool done(){ 
          bool empty = readyTasks_.empty();
          return empty;
        }








        const std::priority_queue<T, std::vector<T>, MCTTaskCompare > & GetQueue(){return  readyTasks_;}

      protected:
        std::priority_queue<T, std::vector<T>, MCTTaskCompare > readyTasks_;
    };


  template <class T >
    class PRScheduler: public Scheduler<T>{
      public:
        PRScheduler(){
        }

        virtual T& top(){
          T & ref = const_cast<T&>(readyTasks_.top());
          return ref;
        }
        virtual void pop() {
          readyTasks_.pop();
        }

        virtual void push( const T& val) { 
          readyTasks_.push(val);
        }


        virtual unsigned int size() {
          int count = readyTasks_.size();
          return count;
        }

        virtual bool done(){ 
          bool empty = readyTasks_.empty();
          return empty;
        }





        const std::priority_queue<T, std::vector<T>, PRTaskCompare > & GetQueue(){return  readyTasks_;}
      protected:
        std::priority_queue<T, std::vector<T>, PRTaskCompare > readyTasks_;
    };


  template <class T >
    class FIFOScheduler: public Scheduler<T>{
      public:
        FIFOScheduler(){
        }

        virtual T& top(){ 
          T & ref = const_cast<T&>(readyTasks_.front());
          return ref;
        }
        virtual void pop() { 
          readyTasks_.pop();
        }

        virtual void push( const T& val) { 
          readyTasks_.push(val);
        }
        virtual unsigned int size() {
          int count = readyTasks_.size();
          return count;
        }

        virtual bool done(){ 
          bool empty = readyTasks_.empty();
          return empty;
        }
        const std::queue<T> & GetQueue(){return  readyTasks_;}
      protected:
        std::queue<T> readyTasks_;
    };

} //end namespace symPACK

#include "sympack/impl/Scheduler_impl.hpp"


#endif // _SCHEDULER_DECL_HPP_
