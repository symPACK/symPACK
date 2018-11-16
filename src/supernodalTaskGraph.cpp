#include <sympack/symPACKMatrix.hpp>

namespace symPACK{

   supernodalTaskGraph::supernodalTaskGraph(){
    localTaskCount_=0;
  }

   supernodalTaskGraph::supernodalTaskGraph( const supernodalTaskGraph& g ){
    //(*this) = g;
    localTaskCount_ = g.localTaskCount_;
    taskLists_.resize(g.taskLists_.size(),NULL);
    for(int i = 0; i<taskLists_.size(); ++i){
      if(g.taskLists_[i] != NULL){
        taskLists_[i] = new std::list<FBTask>(); 
        taskLists_[i]->insert(taskLists_[i]->end(),g.taskLists_[i]->begin(),g.taskLists_[i]->end());
      }
    }
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
