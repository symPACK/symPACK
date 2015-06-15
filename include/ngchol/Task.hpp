#ifndef _TASK_DECL_HPP_
#define _TASK_DECL_HPP_

#include "ngchol/Environment.hpp"

namespace LIBCHOLESKY{

  class IncomingMessage;

  enum TaskType {FACTOR, AGGREGATE, UPDATE};


  struct FBTask{
    TaskType type;
    Int src_snode_id;
    Int tgt_snode_id;
    //dependencies
    Int remote_deps;
    Int local_deps;
    //unused but preparing for task scheduling priorities
    Int rank;
    //list of incoming messages
    std::list<IncomingMessage*> data;
    FBTask():rank(-1),remote_deps(0),local_deps(0){}
  };

  struct FBTaskCompare{
    bool operator()(const FBTask & a,const FBTask & b) const
    {

      bool b_factor = b.tgt_snode_id == b.src_snode_id;


      //use the ranks first
      if(a.rank>=0 && b.rank>=0){
        return a.rank<b.rank;
      }
    


      //check whether it is an update or a factorization
      bool a_factor = a.tgt_snode_id == a.src_snode_id;

      //if same type apply the usual priorities
//      if(a_factor == b_factor){

      //use the classic priorities otherwise
      if(a.tgt_snode_id>b.tgt_snode_id){
        return true;
      }
      else if(a.tgt_snode_id==b.tgt_snode_id){
        return a.src_snode_id>b.src_snode_id;
      }
      else{
        return false;
      }

//      }
//      else if (a_factor){
//        if(a.tgt_snode_id
//      }
//      else if (b_factor){
//
//      }
    }
  };


  

  struct SnodeUpdateFB{
    TaskType type;
    Int src_snode_id;
    Int tgt_snode_id;
    //dependencies
    Int remote_deps;
    Int local_deps;
    //unused but preparing for task scheduling priorities
    Int rank;
    SnodeUpdateFB():rank(-1),remote_deps(0),local_deps(0){}
  };


   struct SnodeUpdateFBCompare{
    bool operator()(const SnodeUpdateFB & a,const SnodeUpdateFB & b) const
    {

      bool b_factor = b.tgt_snode_id == b.src_snode_id;


      //use the ranks first
      if(a.rank>=0 && b.rank>=0){
        return a.rank<b.rank;
      }
    


      //check whether it is an update or a factorization
      bool a_factor = a.tgt_snode_id == a.src_snode_id;

      //if same type apply the usual priorities
//      if(a_factor == b_factor){

      //use the classic priorities otherwise
      if(a.tgt_snode_id>b.tgt_snode_id){
        return true;
      }
      else if(a.tgt_snode_id==b.tgt_snode_id){
        return a.src_snode_id>b.src_snode_id;
      }
      else{
        return false;
      }

//      }
//      else if (a_factor){
//        if(a.tgt_snode_id
//      }
//      else if (b_factor){
//
//      }
    }
  };

 



/*
    bool CompareSnodeUpdateFB(const SnodeUpdateFB & a,const SnodeUpdateFB & b)
    {

      bool b_factor = b.tgt_snode_id == b.src_snode_id;


      //use the ranks first
      if(a.rank>=0 && b.rank>=0){
        return a.rank<b.rank;
      }
    


      //check whether it is an update or a factorization
      bool a_factor = a.tgt_snode_id == a.src_snode_id;

      //if same type apply the usual priorities
//      if(a_factor == b_factor){

      //use the classic priorities otherwise
      if(a.tgt_snode_id>b.tgt_snode_id){
        return true;
      }
      else if(a.tgt_snode_id==b.tgt_snode_id){
        return a.src_snode_id>b.src_snode_id;
      }
      else{
        return false;
      }

//      }
//      else if (a_factor){
//        if(a.tgt_snode_id
//      }
//      else if (b_factor){
//
//      }
    }
 
*/

}

#endif
