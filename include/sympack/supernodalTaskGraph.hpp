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

//Declaration
namespace symPACK{

      class supernodalTaskGraph{
        template<typename T> friend class symPACKMatrix;
        protected:
        std::vector<std::list<FBTask> * > taskLists_;
        Int localTaskCount_;
        
        public:

        supernodalTaskGraph( );
        supernodalTaskGraph( const supernodalTaskGraph& g );
        supernodalTaskGraph& operator=( const supernodalTaskGraph& g );
        ~supernodalTaskGraph();

        void removeTask(std::list<FBTask>::iterator & taskit); 
        std::list<FBTask>::iterator addTask(FBTask & task);
        Int getTaskCount();
        Int setTaskCount(Int value);
        Int increaseTaskCount();
        Int decreaseTaskCount();
        std::list<FBTask>::iterator find_task(Int src, Int tgt, TaskType type );
      };
}

//Implementation
namespace symPACK{

  supernodalTaskGraph::supernodalTaskGraph(){
    localTaskCount_=0;
  }
  supernodalTaskGraph::supernodalTaskGraph( const supernodalTaskGraph& g ){
    (*this) = g;
  }

  supernodalTaskGraph& supernodalTaskGraph::operator=( const supernodalTaskGraph& g ){
    localTaskCount_ = g.localTaskCount_;
    taskLists_.resize(g.taskLists_.size(),NULL);
    for(int i = 0; i<taskLists_.size(); ++i){
      if(g.taskLists_[i] != NULL){
        taskLists_[i] = new std::list<FBTask>(); 
        taskLists_[i]->insert(taskLists_[i]->end(),g.taskLists_[i]->begin(),g.taskLists_[i]->end());
      }
    }

  }

  supernodalTaskGraph::~supernodalTaskGraph(){
    for(int i = 0; i<taskLists_.size(); ++i){
      if(taskLists_[i] != NULL){
        delete taskLists_[i];
      }
    }

  }        

  void supernodalTaskGraph::removeTask(std::list<FBTask>::iterator & taskit){
    taskLists_[taskit->tgt_snode_id-1]->erase(taskit);
    //decreaseTaskCount();
  }

  std::list<FBTask>::iterator supernodalTaskGraph::addTask(FBTask & task){
    if(taskLists_[task.tgt_snode_id-1] == NULL){
      taskLists_[task.tgt_snode_id-1]=new std::list<FBTask>();
    }
    taskLists_[task.tgt_snode_id-1]->push_back(task);
    increaseTaskCount();

    std::list<FBTask>::iterator taskit = --taskLists_[task.tgt_snode_id-1]->end();
    return taskit;
  }



  Int supernodalTaskGraph::getTaskCount()
  {
    Int val = localTaskCount_;
    return val;
  }

  Int supernodalTaskGraph::setTaskCount(Int value)
  {
    localTaskCount_ = value;
    return value;
  }

  Int supernodalTaskGraph::increaseTaskCount()
  {
    Int val;
    val = ++localTaskCount_;
    return val;
  }

  Int supernodalTaskGraph::decreaseTaskCount()
  {
    Int val;
    val = --localTaskCount_;
    return val;
  }

  std::list<FBTask>::iterator supernodalTaskGraph::find_task(Int src, Int tgt, TaskType type )
  {
    scope_timer(a,FB_FIND_TASK);
    //find task corresponding to curUpdate
    auto taskit = taskLists_[tgt-1]->begin();
    for(;
        taskit!=taskLists_[tgt-1]->end();
        taskit++){
      if(taskit->src_snode_id==src && taskit->tgt_snode_id==tgt && taskit->type==type){
        break;
      }
    }
    return taskit;
  }

}

#endif // _SUPERNODAL_TASK_GRAPH_DECL_HPP_
