#ifndef _SCHEDULER_DECL_HPP_
#define _SCHEDULER_DECL_HPP_

#include "ngchol/Types.hpp"
#include "ngchol/Task.hpp"
#include <list>
#include <queue>

//#define TASKCOMPARE  MCTTaskCompare
namespace LIBCHOLESKY{

template <class T >
class Scheduler{
  public:
    Scheduler(){
    }

    virtual const T& top() =0;
    virtual void pop() =0;

    virtual void push( const T& val) =0;
    virtual unsigned int size() =0;

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
    const std::priority_queue<T, vector<T>, FBTaskCompare > & GetQueue(){return  readyTasks_;}
  protected:
      std::priority_queue<T, vector<T>, FBTaskCompare > readyTasks_;
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
    const std::priority_queue<T, vector<T>, MCTTaskCompare > & GetQueue(){return  readyTasks_;}
  protected:
      std::priority_queue<T, vector<T>, MCTTaskCompare > readyTasks_;
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
    const std::priority_queue<T, vector<T>, PRTaskCompare > & GetQueue(){return  readyTasks_;}
  protected:
      std::priority_queue<T, vector<T>, PRTaskCompare > readyTasks_;
};






}
#endif // _SCHEDULER_DECL_HPP_
