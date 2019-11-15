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
