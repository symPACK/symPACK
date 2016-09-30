#ifndef _SCHEDULER_DECL_HPP_
#define _SCHEDULER_DECL_HPP_

#include "sympack/Types.hpp"
#include "sympack/Task.hpp"
#include <list>
#include <queue>

#ifdef MULTITHREADING
#include <omp.h>
#endif

//#define TASKCOMPARE  MCTTaskCompare
namespace symPACK{

  template <class T >
    class Scheduler{
      public:
        Scheduler(){
#ifdef MULTITHREADING
          omp_init_lock(&scheduler_lock);
#endif
        }

#ifdef MULTITHREADING
        ~Scheduler(){
          omp_destroy_lock(&scheduler_lock);
        }
#endif

        virtual T& top() =0;
        virtual void pop() =0;

        virtual void push( const T& val) =0;
        virtual unsigned int size() =0;

        //    virtual const std::priority_queue<T, std::vector<T>, TaskCompare >  GetQueue() =0;

        virtual bool done() =0;

      protected:
#ifdef MULTITHREADING
        omp_lock_t scheduler_lock;
#endif
    };




  template <class T >
    class DLScheduler: public Scheduler<T>{
      public:
        DLScheduler(){
        }

        virtual T& top(){ 
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          T & ref = const_cast<T&>(readyTasks_.top());
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
          return ref;
        }
        virtual void pop() {
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          readyTasks_.pop();
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
        }

        virtual void push( const T& val) { 
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          readyTasks_.push(val);
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
        }
        virtual unsigned int size() {
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          int count = readyTasks_.size();
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
          return count;
        }

        virtual bool done(){ 
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          bool empty = readyTasks_.empty();
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
          return empty;
        }



        const std::priority_queue<T, std::vector<T>, FBTaskCompare > & GetQueue() {return  readyTasks_;}
      protected:
        std::priority_queue<T, std::vector<T>, FBTaskCompare > readyTasks_;
    };


  template <class T >
    class MCTScheduler: public Scheduler<T>{
      public:
        MCTScheduler(){
        }

        virtual T& top(){
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          T & ref = const_cast<T&>(readyTasks_.top());
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
          return ref;
        }
        virtual void pop() { 
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          readyTasks_.pop();
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
        }

        virtual void push( const T& val) { 
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          readyTasks_.push(val);
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
        }



        virtual unsigned int size() {
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          int count = readyTasks_.size();
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
          return count;
        }

        virtual bool done(){ 
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          bool empty = readyTasks_.empty();
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
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
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          T & ref = const_cast<T&>(readyTasks_.top());
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
          return ref;
        }
        virtual void pop() {
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          readyTasks_.pop();
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
        }

        virtual void push( const T& val) { 
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          readyTasks_.push(val);
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
        }


        virtual unsigned int size() {
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          int count = readyTasks_.size();
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
          return count;
        }

        virtual bool done(){ 
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          bool empty = readyTasks_.empty();
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
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
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          T & ref = const_cast<T&>(readyTasks_.front());
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
          return ref;
        }
        virtual void pop() { 
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          readyTasks_.pop();
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
        }

        virtual void push( const T& val) { 
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          readyTasks_.push(val);
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
        }
        virtual unsigned int size() {
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          int count = readyTasks_.size();
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
          return count;
        }

        virtual bool done(){ 
#ifdef MULTITHREADING
          omp_set_lock(&this->scheduler_lock);
#endif
          bool empty = readyTasks_.empty();
#ifdef MULTITHREADING
          omp_unset_lock(&this->scheduler_lock);
#endif
          return empty;
        }
        const std::queue<T> & GetQueue(){return  readyTasks_;}
      protected:
        std::queue<T> readyTasks_;
    };





}
#endif // _SCHEDULER_DECL_HPP_
