/*
 "symPACK" Copyright (c) 2016, The Regents of the University of California,
 through Lawrence Berkeley National Laboratory (subject to receipt of any
 required approvals from the U.S. Dept. of Energy).  All rights reserved.
  
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
  
 (1) Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
  
 (2) Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
  
 (3) Neither the name of the University of California, Lawrence Berkeley
 National Laboratory, U.S. Dept. of Energy nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.
  
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef _SUPERNODAL_TASK_GRAPH_DECL_HPP_
#define _SUPERNODAL_TASK_GRAPH_DECL_HPP_

#include "sympack/Environment.hpp"

#include "sympack/Types.hpp"
#include "sympack/Task.hpp"

#include <upcxx.h>

#include <list>
#include <vector>
#include <limits>
#include <numeric>

//#define SP_THREADS

#ifdef SP_THREADS
#include <future>
#include <mutex>
#endif

//Declaration
namespace symPACK{


  template <typename Task = FBTask>
      class supernodalTaskGraph{
        template<typename T> friend class symPACKMatrix;
        protected:
        typename std::vector<typename std::list<Task> * > taskLists_;
        Int localTaskCount_;
        
        public:

        supernodalTaskGraph( );
        supernodalTaskGraph( const supernodalTaskGraph& g );
        supernodalTaskGraph& operator=( const supernodalTaskGraph& g );
        ~supernodalTaskGraph();

        void removeTask(typename std::list<Task>::iterator & taskit); 
        typename std::list<Task>::iterator addTask(Task & task);
        Int getTaskCount();
        Int setTaskCount(Int value);
        Int increaseTaskCount();
        Int decreaseTaskCount();

        template<typename Type = TaskType>
        typename std::list<Task>::iterator find_task(Int src, Int tgt, Type type );

      };
}

//Implementation
namespace symPACK{

  template <typename Task>
  supernodalTaskGraph<Task>::supernodalTaskGraph(){
    localTaskCount_=0;
  }
  template <typename Task>
  supernodalTaskGraph<Task>::supernodalTaskGraph( const supernodalTaskGraph<Task>& g ){
    (*this) = g;
  }

  template <typename Task>
  supernodalTaskGraph<Task>& supernodalTaskGraph<Task>::operator=( const supernodalTaskGraph<Task>& g ){
    localTaskCount_ = g.localTaskCount_;
    taskLists_.resize(g.taskLists_.size(),NULL);
    for(int i = 0; i<taskLists_.size(); ++i){
      if(g.taskLists_[i] != NULL){
        taskLists_[i] = new std::list<Task>(); 
        taskLists_[i]->insert(taskLists_[i]->end(),g.taskLists_[i]->begin(),g.taskLists_[i]->end());
      }
    }
  }

  template <typename Task>
  supernodalTaskGraph<Task>::~supernodalTaskGraph(){
    for(int i = 0; i<taskLists_.size(); ++i){
      if(taskLists_[i] != NULL){
        delete taskLists_[i];
      }
    }
  }        

  template <typename Task>
  void supernodalTaskGraph<Task>::removeTask(typename std::list<Task>::iterator & taskit){
    throw std::logic_error( "removeTask(std::list<Task>::iterator &taskit) is not implemented for that type.");
  }

  template <>
  void supernodalTaskGraph<FBTask>::removeTask(std::list<FBTask>::iterator & taskit){
    taskLists_[taskit->tgt_snode_id-1]->erase(taskit);
  }
  
  template <>
  void supernodalTaskGraph<CompTask>::removeTask(std::list<CompTask>::iterator & taskit){
    Int * meta = reinterpret_cast<Int*>(taskit->meta.data());
    Int tgt = meta[1];
    taskLists_[tgt-1]->erase(taskit);
  }

  template <typename Task>
  typename std::list<Task>::iterator supernodalTaskGraph<Task>::addTask(Task & task){
  }

  template <>
  std::list<FBTask>::iterator supernodalTaskGraph<FBTask>::addTask(FBTask & task){
    if(taskLists_[task.tgt_snode_id-1] == NULL){
      taskLists_[task.tgt_snode_id-1]=new std::list<FBTask>();
    }
    taskLists_[task.tgt_snode_id-1]->push_back(task);
    increaseTaskCount();

    std::list<FBTask>::iterator taskit = --taskLists_[task.tgt_snode_id-1]->end();
    return taskit;
  }

  template <>
  std::list<CompTask>::iterator supernodalTaskGraph<CompTask>::addTask(CompTask & task){

    Int * meta = reinterpret_cast<Int*>(task.meta.data());
    Int tgt = meta[1];

    if(taskLists_[tgt-1] == NULL){
      taskLists_[tgt-1]=new std::list<CompTask>();
    }
    taskLists_[tgt-1]->push_back(task);
    increaseTaskCount();

    std::list<CompTask>::iterator taskit = --taskLists_[tgt-1]->end();
    return taskit;
  }



  template <typename Task>
  Int supernodalTaskGraph<Task>::getTaskCount()
  {
    Int val = localTaskCount_;
    return val;
  }

  template <typename Task>
  Int supernodalTaskGraph<Task>::setTaskCount(Int value)
  {
    localTaskCount_ = value;
    return value;
  }

  template <typename Task>
  Int supernodalTaskGraph<Task>::increaseTaskCount()
  {
    Int val;
    val = ++localTaskCount_;
    return val;
  }

  template <typename Task>
  Int supernodalTaskGraph<Task>::decreaseTaskCount()
  {
    Int val;
    val = --localTaskCount_;
    return val;
  }

  template <typename Task>
  template <typename Type>
  typename std::list<Task>::iterator supernodalTaskGraph<Task>::find_task(Int src, Int tgt, Type type )
  {
    throw std::logic_error( "find_task(Int src, Int tgt, TaskType type) is not implemented for that type.");
    auto taskit = taskLists_[tgt-1]->begin();
    return taskit;
  }

  template <>
  template <>
  std::list<CompTask>::iterator supernodalTaskGraph<CompTask>::find_task(Int src, Int tgt, Solve::op_type type )
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
  std::list<FBTask>::iterator supernodalTaskGraph<FBTask>::find_task(Int src, Int tgt, TaskType type )
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


namespace symPACK{

      class taskGraph{
        protected:
        
        public:
          typedef typename std::map<GenericTask::id_type, std::shared_ptr<GenericTask> > task_list; 
          typedef typename task_list::iterator task_iterator; 

          template<typename T> friend class symPACKMatrix;
          task_list tasks_;

        taskGraph(){};
        taskGraph( const taskGraph& g );
        taskGraph& operator=( const taskGraph& g );
        //taskGraph( );
        //supernodalTaskGraph( const supernodalTaskGraph& g );
        //supernodalTaskGraph& operator=( const supernodalTaskGraph& g );
        //~taskGraph();

        void removeTask( const GenericTask::id_type & hash); 
        void addTask(std::shared_ptr<GenericTask> & task);
        Int getTaskCount();
        task_iterator find_task( const GenericTask::id_type & hash );

#ifdef SP_THREADS
        std::mutex list_mutex_;
#endif
      };
}

//Implementation
namespace symPACK{

  void taskGraph::removeTask(const GenericTask::id_type & hash){
#ifdef SP_THREADS
    std::lock_guard<std::mutex> lock(list_mutex_);
#endif
    this->tasks_.erase(hash);
  }

  void taskGraph::addTask(std::shared_ptr<GenericTask> & task){

#ifdef SP_THREADS
    std::lock_guard<std::mutex> lock(list_mutex_);
#endif
//    return this->tasks_[hash] = task;
    this->tasks_.emplace(task->id,task);
  }

  Int taskGraph::getTaskCount()
  {
#ifdef SP_THREADS
    std::lock_guard<std::mutex> lock(list_mutex_);
#endif
    Int val = tasks_.size();
    return val;
  }

  taskGraph::task_iterator taskGraph::find_task( const GenericTask::id_type & hash ){
#ifdef SP_THREADS
    std::lock_guard<std::mutex> lock(list_mutex_);
#endif
    return tasks_.find(hash);
  }

  taskGraph::taskGraph( const taskGraph& g ){
    (*this) = g;
  }

  taskGraph& taskGraph::operator=( const taskGraph& g ){
    tasks_ = g.tasks_;
  }




}


#endif // _SUPERNODAL_TASK_GRAPH_DECL_HPP_
