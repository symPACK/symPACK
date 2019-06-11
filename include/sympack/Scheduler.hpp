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

template<typename T>
  class mpmc_bounded_queue
  {
    public:
      mpmc_bounded_queue(size_t buffer_size = 1024)
        : buffer_size_(buffer_size), buffer_(new cell_t [buffer_size_])
          , buffer_mask_(buffer_size - 1)
    {
      assert((buffer_size >= 2) &&
          ((buffer_size & (buffer_size - 1)) == 0));
      for (size_t i = 0; i != buffer_size; i += 1)
        buffer_[i].sequence_.store(i, std::memory_order_relaxed);
      enqueue_pos_.store(0, std::memory_order_relaxed);
      dequeue_pos_.store(0, std::memory_order_relaxed);
    }



      mpmc_bounded_queue( const mpmc_bounded_queue & other)
        : buffer_size_(other.buffer_size_), buffer_(new cell_t [buffer_size_])
          , buffer_mask_(buffer_size_ - 1)
    {
      assert((buffer_size_ >= 2) &&
          ((buffer_size_ & (buffer_size_ - 1)) == 0));
      for (size_t i = 0; i != buffer_size_; i += 1)
        buffer_[i].sequence_.store(i, std::memory_order_relaxed);
      enqueue_pos_.store(0, std::memory_order_relaxed);
      dequeue_pos_.store(0, std::memory_order_relaxed);
    }



      ~mpmc_bounded_queue()
      {
        delete [] buffer_;
      }

      bool enqueue(T const& data)
      {
        cell_t* cell;
        size_t pos = enqueue_pos_.load(std::memory_order_relaxed);
        for (;;)
        {
          cell = &buffer_[pos & buffer_mask_];
          size_t seq = 
            cell->sequence_.load(std::memory_order_acquire);
          intptr_t dif = (intptr_t)seq - (intptr_t)pos;
          if (dif == 0)
          {
//            if (enqueue_pos_.compare_exchange_weak
//                (pos, pos + 1, std::memory_order_relaxed))
              enqueue_pos_ = pos + 1;
              break;
          }
          else if (dif < 0)
            return false;
          else
            pos = enqueue_pos_.load(std::memory_order_relaxed);
        }
        cell->data_ = data;
        cell->sequence_.store(pos + 1, std::memory_order_release);
        return true;
      }

      bool dequeue(T& data)
      {
        cell_t* cell;
        size_t pos = dequeue_pos_.load(std::memory_order_relaxed);
        for (;;)
        {
          cell = &buffer_[pos & buffer_mask_];
          size_t seq = 
            cell->sequence_.load(std::memory_order_acquire);
          intptr_t dif = (intptr_t)seq - (intptr_t)(pos + 1);
          if (dif == 0)
          {
            if (dequeue_pos_.compare_exchange_weak
                (pos, pos + 1, std::memory_order_relaxed))
              break;
          }
          else if (dif < 0)
            return false;
          else
            pos = dequeue_pos_.load(std::memory_order_relaxed);
        }
        data = cell->data_;
        cell->sequence_.store
          (pos + buffer_mask_ + 1, std::memory_order_release);
        return true;
      }

    private:
      struct cell_t
      {
        std::atomic<size_t>   sequence_;
        T                     data_;
      };

      static size_t const     cacheline_size = 64;
      typedef char            cacheline_pad_t [cacheline_size];

      cacheline_pad_t         pad4_;
      size_t const            buffer_size_;
      cacheline_pad_t         pad0_;
      cell_t* const           buffer_;
      size_t const            buffer_mask_;
      cacheline_pad_t         pad1_;
      std::atomic<size_t>     enqueue_pos_;
      cacheline_pad_t         pad2_;
      std::atomic<size_t>     dequeue_pos_;
      cacheline_pad_t         pad3_;

      //mpmc_bounded_queue(mpmc_bounded_queue const&);
      void operator = (mpmc_bounded_queue const&);
  };




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
#ifdef NEW_UPCXX
//        void run(MPI_Comm & workcomm, RankGroup & group , taskGraph & graph, std::list< upcxx::future<> > & pFutures);
        void run(MPI_Comm & workcomm, RankGroup & group , taskGraph & graph, upcxx::dist_object<int> & remDealloc, upcxx::team & workteam);
#else
        void run(MPI_Comm & workcomm, RankGroup & group, taskGraph & graph);
#endif

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




        //    virtual void compute_costs(PtrVec & Xlindx, IdxVec & Lindx){
        //
        //
        //        //I is source supernode
        //        Int I = ;
        //
        //            Idx fc = Xsuper_[I-1];
        //            Int width = Xsuper_[I] - Xsuper_[I-1];
        //            Int height = cc_[fc-1];
        //            //cost of factoring curent panel
        //            double local_load = factor_cost(height,width);
        //
        //            //cost of updating ancestors
        //            {
        //              Ptr fi = Xlindx_[I-1];
        //              Ptr li = Xlindx_[I]-1;
        //
        //                  Int m = li-tgt_fr_ptr+1;
        //                  Int n = width;
        //                  Int k = tgt_lr_ptr - tgt_fr_ptr+1; 
        //                  double cost = update_cost(m,n,k);
        //
        //              //parse rows
        //              Int tgt_snode_id = I;
        //              Idx tgt_fr = fc;
        //              Idx tgt_lr = fc;
        //              Ptr tgt_fr_ptr = 0;
        //              Ptr tgt_lr_ptr = 0;
        //              for(Ptr rowptr = fi; rowptr<=li;rowptr++){
        //                Idx row = Lindx_[rowptr-1]; 
        //                if(SupMembership_[row-1]==tgt_snode_id){ continue;}
        //
        //                //we have a new tgt_snode_id
        //                tgt_lr_ptr = rowptr-1;
        //                tgt_lr = Lindx_[rowptr-1 -1];
        //                if(tgt_snode_id !=I){
        //                  Int m = li-tgt_fr_ptr+1;
        //                  Int n = width;
        //                  Int k = tgt_lr_ptr - tgt_fr_ptr+1; 
        //                  double cost = update_cost(m,n,k);
        //                  break;
        //                }
        //                tgt_fr = row;
        //                tgt_fr_ptr = rowptr;
        //                tgt_snode_id = SupMembership_[row-1];
        //              }
        //              //do the last update supernode
        //              //we have a new tgt_snode_id
        //              tgt_lr_ptr = li;
        //              tgt_lr = Lindx_[li -1];
        //              if(tgt_snode_id!=I){
        //                Int m = li-tgt_fr_ptr+1;
        //                Int n = width;
        //                Int k = tgt_lr_ptr - tgt_fr_ptr+1; 
        //                double cost = update_cost(m,n,k);
        //                  if(fan_in_){
        //                    SubTreeLoad[I]+=cost;
        //                    NodeLoad[I]+=cost;
        //                  }
        //                  else{
        //                    SubTreeLoad[tgt_snode_id]+=cost;
        //                    NodeLoad[tgt_snode_id]+=cost;
        //                  }
        //                }
        //            }
        //
        //
        //    }
        //




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
