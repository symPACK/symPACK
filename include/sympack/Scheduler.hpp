#ifndef _SCHEDULER_DECL_HPP_
#define _SCHEDULER_DECL_HPP_

#include "sympack/Types.hpp"
#include "sympack/Task.hpp"
#include <list>
#include <queue>

//#define TASKCOMPARE  MCTTaskCompare
namespace SYMPACK{

template <class T >
class Scheduler{
  public:
    Scheduler(){
    }


    virtual const T& top() =0;
    virtual void pop() =0;

    virtual void push( const T& val) =0;
    virtual unsigned int size() =0;

//    virtual const std::priority_queue<T, SYMPACK::vector<T>, TaskCompare >  GetQueue() =0;

    virtual bool done() =0;
};




template <class T >
class DLScheduler: public Scheduler<T>{
  public:
    DLScheduler(){
    }

    virtual const T& top(){ return readyTasks_.top();}
    virtual void pop() { readyTasks_.pop();}

    virtual void push( const T& val) { readyTasks_.push(val);}
    virtual unsigned int size() { return readyTasks_.size();}

    virtual bool done(){ return readyTasks_.empty();}

    const std::priority_queue<T, SYMPACK::vector<T>, FBTaskCompare > & GetQueue() {return  readyTasks_;}
  protected:
      std::priority_queue<T, SYMPACK::vector<T>, FBTaskCompare > readyTasks_;
};


template <class T >
class MCTScheduler: public Scheduler<T>{
  public:
    MCTScheduler(){
    }

    virtual const T& top(){ return readyTasks_.top();}
    virtual void pop() { readyTasks_.pop();}

    virtual void push( const T& val) { readyTasks_.push(val);}
    virtual unsigned int size() { return readyTasks_.size();}

    virtual bool done(){ return readyTasks_.empty();}
    const std::priority_queue<T, SYMPACK::vector<T>, MCTTaskCompare > & GetQueue(){return  readyTasks_;}




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
      std::priority_queue<T, SYMPACK::vector<T>, MCTTaskCompare > readyTasks_;
};


template <class T >
class PRScheduler: public Scheduler<T>{
  public:
    PRScheduler(){
    }

    virtual const T& top(){ return readyTasks_.top();}
    virtual void pop() { readyTasks_.pop();}

    virtual void push( const T& val) { readyTasks_.push(val);}
    virtual unsigned int size() { return readyTasks_.size();}

    virtual bool done(){ return readyTasks_.empty();}
    const std::priority_queue<T, SYMPACK::vector<T>, PRTaskCompare > & GetQueue(){return  readyTasks_;}
  protected:
      std::priority_queue<T, SYMPACK::vector<T>, PRTaskCompare > readyTasks_;
};


template <class T >
class FIFOScheduler: public Scheduler<T>{
  public:
    FIFOScheduler(){
    }

    virtual const T& top(){ return readyTasks_.front();}
    virtual void pop() { readyTasks_.pop();}

    virtual void push( const T& val) { readyTasks_.push(val);}
    virtual unsigned int size() { return readyTasks_.size();}

    virtual bool done(){ return readyTasks_.empty();}
    const std::queue<T> & GetQueue(){return  readyTasks_;}
  protected:
      std::queue<T> readyTasks_;
};





}
#endif // _SCHEDULER_DECL_HPP_
