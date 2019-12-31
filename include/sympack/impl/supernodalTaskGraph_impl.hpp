#ifndef _SUPERNODAL_TASK_GRAPH_IMPL_HPP_
#define _SUPERNODAL_TASK_GRAPH_IMPL_HPP_

#include "sympack/Environment.hpp"

#include "sympack/Types.hpp"
#include "sympack/Task.hpp"


#include <list>
#include <vector>
#include <limits>
#include <numeric>

//#define SP_THREADS

#ifdef SP_THREADS
#include <future>
#include <mutex>
#endif

//Implementation
namespace symPACK{

  template <typename Task>
    inline supernodalTaskGraph<Task>::supernodalTaskGraph(){
      localTaskCount_=0;
    }
  template <typename Task>
    inline supernodalTaskGraph<Task>::supernodalTaskGraph( const supernodalTaskGraph<Task>& g ){
      (*this) = g;
    }

  template <typename Task>
    inline supernodalTaskGraph<Task>& supernodalTaskGraph<Task>::operator=( const supernodalTaskGraph<Task>& g ){
      localTaskCount_ = g.localTaskCount_;
      taskLists_.resize(g.taskLists_.size(),nullptr);
      for(int i = 0; i<taskLists_.size(); ++i){
        if(g.taskLists_[i] != nullptr){
          taskLists_[i] = new std::list<Task>(); 
          taskLists_[i]->insert(taskLists_[i]->end(),g.taskLists_[i]->begin(),g.taskLists_[i]->end());
        }
      }
      return *this;
    }

  template <typename Task>
    inline supernodalTaskGraph<Task>::~supernodalTaskGraph(){
      for(int i = 0; i<taskLists_.size(); ++i){
        if(taskLists_[i] != nullptr){
          delete taskLists_[i];
        }
      }
    }        

  template <typename Task>
    inline void supernodalTaskGraph<Task>::removeTask(typename std::list<Task>::iterator & taskit){
      throw std::logic_error( "removeTask(std::list<Task>::iterator &taskit) is not implemented for that type.");
    }

  template <>
    inline void supernodalTaskGraph<FBTask>::removeTask(std::list<FBTask>::iterator & taskit){
      taskLists_[taskit->tgt_snode_id-1]->erase(taskit);
    }

  template <>
    inline void supernodalTaskGraph<CompTask>::removeTask(std::list<CompTask>::iterator & taskit){
      Int * meta = reinterpret_cast<Int*>(taskit->meta.data());
      Int tgt = meta[1];
      taskLists_[tgt-1]->erase(taskit);
    }

  template <typename Task>
    inline typename std::list<Task>::iterator supernodalTaskGraph<Task>::addTask(Task & task){
    }

  template <>
    inline std::list<FBTask>::iterator supernodalTaskGraph<FBTask>::addTask(FBTask & task){
      if(taskLists_[task.tgt_snode_id-1] == nullptr){
        taskLists_[task.tgt_snode_id-1]=new std::list<FBTask>();
      }
      taskLists_[task.tgt_snode_id-1]->push_back(task);
      increaseTaskCount();

      std::list<FBTask>::iterator taskit = --taskLists_[task.tgt_snode_id-1]->end();
      return taskit;
    }

  template <>
    inline std::list<CompTask>::iterator supernodalTaskGraph<CompTask>::addTask(CompTask & task){

      Int * meta = reinterpret_cast<Int*>(task.meta.data());
      Int tgt = meta[1];

      if(taskLists_[tgt-1] == nullptr){
        taskLists_[tgt-1]=new std::list<CompTask>();
      }
      taskLists_[tgt-1]->push_back(task);
      increaseTaskCount();

      std::list<CompTask>::iterator taskit = --taskLists_[tgt-1]->end();
      return taskit;
    }



  template <typename Task>
    inline Int supernodalTaskGraph<Task>::getTaskCount()
    {
      Int val = localTaskCount_;
      return val;
    }

  template <typename Task>
    inline Int supernodalTaskGraph<Task>::setTaskCount(Int value)
    {
      localTaskCount_ = value;
      return value;
    }

  template <typename Task>
    inline Int supernodalTaskGraph<Task>::increaseTaskCount()
    {
      Int val;
      val = ++localTaskCount_;
      return val;
    }

  template <typename Task>
    inline Int supernodalTaskGraph<Task>::decreaseTaskCount()
    {
      Int val;
      val = --localTaskCount_;
      return val;
    }

  template <typename Task>
    template <typename Type>
    inline typename std::list<Task>::iterator supernodalTaskGraph<Task>::find_task(Int src, Int tgt, Type type )
    {
      throw std::logic_error( "find_task(Int src, Int tgt, TaskType type) is not implemented for that type.");
      auto taskit = taskLists_[tgt-1]->begin();
      return taskit;
    }

  template <>
    template <>
    inline std::list<CompTask>::iterator supernodalTaskGraph<CompTask>::find_task(Int src, Int tgt, Solve::op_type type )
    {
      scope_timer(a,FB_FIND_TASK);

      //find task corresponding to curUpdate
      auto taskit = taskLists_[tgt-1]->begin();
      for( ;
          taskit!=this->taskLists_[tgt-1]->end();
          taskit++){

        Int * meta = reinterpret_cast<Int*>(taskit->meta.data());
        Int tsrc = meta[0];
        Int ttgt = meta[1];
        const Solve::op_type & ttype = *reinterpret_cast<Solve::op_type*>(&meta[2]);

        if(tsrc==src && ttgt==tgt && ttype==type){
          break;
        }
      }
      return taskit;
    }

  template <>
    template <>
    inline std::list<FBTask>::iterator supernodalTaskGraph<FBTask>::find_task(Int src, Int tgt, TaskType type )
    {
      scope_timer(a,FB_FIND_TASK);
      //find task corresponding to curUpdate
      auto taskit = taskLists_[tgt-1]->begin();
      for(;
          taskit!=this->taskLists_[tgt-1]->end();
          taskit++){
        if(taskit->src_snode_id==src && taskit->tgt_snode_id==tgt && taskit->type==type){
          break;
        }
      }
      return taskit;
    }

}


//Implementation
namespace symPACK{

  inline void taskGraph::removeTask(const GenericTask::id_type & hash){
#ifdef SP_THREADS
    if(Multithreading::NumThread>1){
      std::lock_guard<std::recursive_mutex> lk(list_mutex_);
      this->tasks_.erase(hash);
    }
    else
#endif
      this->tasks_.erase(hash);
  }

  inline void taskGraph::addTask(std::shared_ptr<GenericTask> & task){

#ifdef SP_THREADS
    if(Multithreading::NumThread>1){
      std::lock_guard<std::recursive_mutex> lk(list_mutex_);
      this->tasks_.emplace(task->id,std::move(task));
    }
    else
#endif
      this->tasks_.emplace(task->id,std::move(task));
  }

  inline Int taskGraph::getTaskCount()
  {
    Int val = 0;
#ifdef SP_THREADS
    if(Multithreading::NumThread>1){
      std::lock_guard<std::recursive_mutex> lk(list_mutex_);
      val = tasks_.size();
    }
    else
#endif
      val = tasks_.size();
    return val;
  }

  inline taskGraph::task_iterator taskGraph::find_task( const GenericTask::id_type & hash ){
#ifdef SP_THREADS
    if(Multithreading::NumThread>1){
      std::lock_guard<std::recursive_mutex> lock(list_mutex_);
      return tasks_.find(hash);
    }
    else
#endif
      return tasks_.find(hash);
  }

  inline taskGraph::taskGraph( const taskGraph& g ){
    (*this) = g;
  }

  inline taskGraph& taskGraph::operator=( const taskGraph& g ){
    tasks_ = g.tasks_;
    return *this;
  }




}


#endif // _SUPERNODAL_TASK_GRAPH_DECL_HPP_
