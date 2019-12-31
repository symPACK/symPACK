#ifndef _SUPERNODAL_TASK_GRAPH_DECL_HPP_
#define _SUPERNODAL_TASK_GRAPH_DECL_HPP_

#include "sympack/Environment.hpp"

#include "sympack/Types.hpp"
#include "sympack/Task.hpp"

#include <list>
#include <vector>
#include <limits>
#include <numeric>
#include <unordered_map>
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
          typedef typename std::unordered_map<GenericTask::id_type, std::shared_ptr<GenericTask> > task_list; 
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
        std::shared_ptr<GenericTask> updateTask(task_iterator & taskit, int loc, int rem ){
#ifdef SP_THREADS
          std::lock_guard<std::recursive_mutex> lk(list_mutex_);
#endif
          bassert(taskit != tasks_.end());
          taskit->second->local_deps-= loc;
          taskit->second->remote_deps-= rem;

          std::shared_ptr<GenericTask> ptr = nullptr;
          if(taskit->second->remote_deps==0 && taskit->second->local_deps==0){
            ptr = taskit->second;
            removeTask(taskit->first);
          }

          return ptr;
        }

        std::shared_ptr<GenericTask> updateTask(const GenericTask::id_type & hash, int loc, int rem ){
#ifdef SP_THREADS
          std::lock_guard<std::recursive_mutex> lk(list_mutex_);
#endif
          auto taskit = find_task(hash);
          return updateTask(taskit,loc,rem);
        }


#ifdef SP_THREADS
        std::recursive_mutex list_mutex_;
#endif
      };
}

#include "sympack/impl/supernodalTaskGraph_impl.hpp"

#endif // _SUPERNODAL_TASK_GRAPH_DECL_HPP_
