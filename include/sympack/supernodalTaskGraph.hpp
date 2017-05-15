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

        void removeTask( const GenericTask::id_type & hash); 
        void addTask(std::shared_ptr<GenericTask> & task);
        Int getTaskCount();
        task_iterator find_task( const GenericTask::id_type & hash );

#ifdef SP_THREADS
        std::mutex list_mutex_;
#endif
      };
}

#include "sympack/impl/supernodalTaskGraph_impl.hpp"

#endif // _SUPERNODAL_TASK_GRAPH_DECL_HPP_
