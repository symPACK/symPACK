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

namespace symPACK{

#ifdef SP_THREADS
#include <atomic>
   
  class SpinLock
  {
    public:
      void lock()
      {
        while(lck.test_and_set(std::memory_order_acquire))
        {}
      }

      void unlock()
      {
        lck.clear(std::memory_order_release);
      }

    private:
      std::atomic_flag lck = ATOMIC_FLAG_INIT;
  };


  template<typename T, typename Queue = std::queue<std::function<T()> > >
    class WorkQueue{
#if 1
      public:
        std::mutex list_mutex_;
        Queue workQueue_;
        std::vector<std::thread> threads;
        std::condition_variable sync;
        bool done = false;

        WorkQueue(Int nthreads){
          for (unsigned int count {0}; count < nthreads; count += 1)
            threads.emplace_back(static_cast<void(WorkQueue<void>::*)()>(&WorkQueue::consume), this);
        }

        WorkQueue(WorkQueue &&) = default;
        WorkQueue &operator=(WorkQueue &&) = delete;


        ~WorkQueue(){
          if(threads.size()>0){
            std::lock_guard<std::mutex> guard(list_mutex_);
            done = true;
            sync.notify_all();
          }
          for (auto &&thread: threads) thread.join();
        }

        void pushTask(std::function<T()> && fut){
          std::unique_lock<std::mutex> lock(list_mutex_);
          workQueue_.push(std::forward<std::function<T()> >(fut));
          sync.notify_one();
        }

      protected:
        void consume()
        {
          std::unique_lock<std::mutex> lock(list_mutex_);
          while (true) {
            if (not workQueue_.empty()) {
              std::function<T()> func { std::move(workQueue_.front()) };
              workQueue_.pop();
              sync.notify_one();
              lock.unlock();
              func();
              lock.lock();
            } else if (done) {
              break;
            } else {
              sync.wait(lock);
            }
          }
        }
#else
      public:
        SpinLock list_lock_;
        Queue workQueue_;
        std::vector<std::thread> threads;
        bool done = false;

        WorkQueue(Int nthreads){
          for (unsigned int count {0}; count < nthreads; count += 1)
            threads.emplace_back(static_cast<void(WorkQueue<void>::*)()>(&WorkQueue::consume), this);
        }

        WorkQueue(WorkQueue &&) = default;
        WorkQueue &operator=(WorkQueue &&) = delete;


        ~WorkQueue(){
          if(threads.size()>0){
            list_lock_.lock();
            done = true;
            list_lock_.unlock();
          }
          for (auto &&thread: threads) thread.join();
        }

        void pushTask(std::function<T()> && fut){
            list_lock_.lock();
          workQueue_.push(std::forward<std::function<T()> >(fut));
            list_lock_.unlock();
        }



      protected:
        void consume()
        {
          while (true) {
            list_lock_.lock();
            if (not workQueue_.empty()) {
              std::function<T()> func { std::move(workQueue_.front()) };
              workQueue_.pop();
              list_lock_.unlock();
              func();
            } else if (done) {
              list_lock_.unlock();
              break;
            } 
            else{
              list_lock_.unlock();
            }
          }
        }

#endif
    };
#endif

  template <class T >
    class Scheduler{
      public:
        Scheduler(){
        }

        virtual T& top() =0;
        virtual void pop() =0;

        virtual void push( const T& val) =0;
        virtual unsigned int size() =0;

        //    virtual const std::priority_queue<T, std::vector<T>, TaskCompare >  GetQueue() =0;

        virtual bool done() =0;


        Int checkIncomingMessages_(Int nthr, taskGraph & graph);
        void run(MPI_Comm & workcomm, Int nthr, taskGraph & graph);


#ifdef SP_THREADS
        std::mutex upcxx_mutex_;
        std::mutex list_mutex_;
#endif
      protected:
    };




  template <class T >
    class DLScheduler: public Scheduler<T>{
      public:
        DLScheduler(){
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


template< typename Task> 
int Scheduler<Task>::checkIncomingMessages_(Int nthr, taskGraph & graph)
{
  abort();
  return 0;
}

template<>
int Scheduler<std::shared_ptr<GenericTask> >::checkIncomingMessages_(Int nthr, taskGraph & graph)
{
  scope_timer(a,CHECK_MESSAGE);

  int num_recv = 0;

  if(nthr>1){
    scope_timer(a,UPCXX_ADVANCE);
#ifdef SP_THREADS
    std::lock_guard<std::mutex> lock(upcxx_mutex_);
#endif
    upcxx::advance();
  }
  else{
    scope_timer(a,UPCXX_ADVANCE);
    upcxx::advance();
  }

  bool comm_found = false;
  IncomingMessage * msg = NULL;

  do{
    msg=NULL;

    {
#ifdef SP_THREADS
      if(nthr>1){
        upcxx_mutex_.lock();
      }
      //std::lock_guard<std::mutex> lock(upcxx_mutex_);
#endif


      {
        //find if there is some finished async comm
        auto it = TestAsyncIncomingMessage();
        if(it!=gIncomingRecvAsync.end()){
          scope_timer(b,RM_MSG_ASYNC);
          msg = *it;
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
          msg = *it;
          gIncomingRecv.erase(it);
        }
      }

#ifdef SP_THREADS
      if(nthr>1){
        upcxx_mutex_.unlock();
      }
#endif
    }

    if(msg!=NULL){
      scope_timer(a,WAIT_AND_UPDATE_DEPS);
      num_recv++;

      bool success = msg->Wait(); 
      //TODO what are the reasons of failure ?
      bassert(success);

      auto taskit = graph.find_task(msg->meta.id);

      bassert(taskit!=graph.tasks_.end());
#if 0
      {
        std::hash<std::string> hash_fn;
        std::stringstream sstr;
        sstr<<50<<"_"<<50<<"_"<<(Int)Solve::op_type::FU;
        auto id = hash_fn(sstr.str());

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
        }
        logfileptr->OFS()<<"  receiving & updating "<<name<<" "<<meta[0]<<"_"<<meta[1]<<std::endl;
      }
#endif
      {
#ifdef SP_THREADS
        //std::lock_guard<std::mutex> lock(list_mutex_);
        if(nthr>1){
          list_mutex_.lock();
        }
#endif
        taskit->second->remote_deps_cnt--;
        taskit->second->data.push_back(msg);

        if(taskit->second->remote_deps_cnt==0 && taskit->second->local_deps_cnt==0){
          this->push(taskit->second);    
          graph.removeTask(taskit->second->id);
        }

#ifdef SP_THREADS
        if(nthr>1){
          list_mutex_.unlock();
        }
#endif
      }
    }
  }while(msg!=NULL);

#ifdef SP_THREADS
      if(nthr>1){
        upcxx_mutex_.lock();
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
      if(nthr>1){
        upcxx_mutex_.unlock();
      }
#endif

  return num_recv;
}


  template <class Task > 
  void Scheduler<Task>::run(MPI_Comm & workcomm, Int nthr,taskGraph & graph){
  }

  template <> 
  void Scheduler<std::shared_ptr<GenericTask> >::run(MPI_Comm & workcomm,Int nthr,taskGraph & graph){

    int np = 1;
    MPI_Comm_size(workcomm,&np);
    int iam = 0;
    MPI_Comm_rank(workcomm,&iam);

    //put rdy tasks in rdy queue
    {
      auto taskit = graph.tasks_.begin();
      while (taskit != graph.tasks_.end()) {
        if(taskit->second->remote_deps_cnt==0 && taskit->second->local_deps_cnt==0){
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
    //Int nthr = (2*std::thread::hardware_concurrency()-1)/np;
    //if(iam==np-1){ nthr= 2*std::thread::hardware_concurrency()-1-(np-1)*nthr; }
    //nthr = 20;

    logfileptr->OFS()<<"Num threads = "<<nthr<<std::endl;
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

#if 1
      {
        WorkQueue<void> queue(nthr);

        Int num_msg =0;
        std::shared_ptr<GenericTask> curTask = nullptr;
        if(nthr>1){
          while(graph.getTaskCount()>0 || !this->done()){
            
            if(np>1){
              num_msg = checkIncomingMessages_(nthr,graph);
            }

            //while(!this->done())
            if(!this->done())
            {
              //Pick a ready task
              {
                std::lock_guard<std::mutex> lock(list_mutex_);
                curTask = this->top();
                this->pop();
              }
              queue.pushTask(std::move(curTask->execute));
            }
          }
        }
        else{
          while(graph.getTaskCount()>0 || !this->done()){
            if(np>1){
              num_msg = checkIncomingMessages_(nthr,graph);
            }
            if(!this->done())
            {
              //Pick a ready task
              curTask = this->top();
              this->pop();
              curTask->execute();
            }
          }
        }

        upcxx::async_wait();
        MPI_Barrier(workcomm);


      }
#else
    while(graph.getTaskCount()>0 || !this->done()){
      int num_msg = checkIncomingMessages_(nthr,graph);

#ifdef SP_THREADS
//      progress_threads();
#endif

#ifdef SP_THREADS
      while(!this->done())
#else
      if(!this->done())
#endif
      {
        //Pick a ready task
        std::shared_ptr<GenericTask> curTask = nullptr;
        {
#ifdef SP_THREADS
        std::lock_guard<std::mutex> lock(list_mutex_);
#endif
        curTask = this->top();
        this->pop();
        }
#if 0
        SparseTask * tmp = ((SparseTask*)curTask.get());
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
        logfileptr->OFS()<<"Running "<<name<<" "<<meta[0]<<"_"<<meta[1]<<std::endl;
#endif

#ifdef SP_THREADS
      progress_threads();
    //logfileptr->OFS()<<"running threads "<<handles.size()<<" vs "<<nthr<<std::endl;
        if( handles.size()<nthr){
            //logfileptr->OFS()<<"thread done"<<std::endl;
          //auto fut = std::async(std::launch::async,curTask->execute);
          //handles.push_back(std::move(fut));
          handles.push_back(std::async(std::launch::async,curTask->execute));
      //progress_threads();
        }
        else{
          curTask->execute();
          break;
        }
#else
        curTask->execute();
#endif
      }
#if 0
      else if (num_msg==0 && graph.getTaskCount()>0){
        logfileptr->OFS()<<"NOT PROGRESSING"<<std::endl;
        for(auto taskit: graph.tasks_){
          log_task(taskit.second);
        }
        logfileptr->OFS()<<"END OF NOT PROGRESSING"<<std::endl;
        upcxx::async_wait();
        //abort();
      }
#endif

    }

#ifdef SP_THREADS
        auto it = handles.begin();
        while(it != handles.end()){
          it->get();
          it = handles.erase(it);
        }
#endif


    upcxx::async_wait();
    MPI_Barrier(workcomm);
#endif

  }





}
#endif // _SCHEDULER_DECL_HPP_
