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
    double rank;
    //list of incoming messages
    std::list<IncomingMessage*> data;
    FBTask():rank(-1.0),remote_deps(0),local_deps(0){}

    void update_rank(){
      if(rank==-1.0 || 1){
            if(type==FACTOR){
              //taskit->rank = 2.0*pow((double)src_snode.Size(),3.0)/3.0 + src_snode.Size()*(src_snode.Size()+1.0)*(src_snode.NRowsBelowBlock(0)-src_snode.Size())/2.0;
              rank = 4.0;
            }
            else if(type == AGGREGATE){
              rank = 2.0;
            }
            else if(type == UPDATE){
              rank = 1.0;
            }
      }
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

  class TaskCompare{
    public:
     virtual bool operator()(const FBTask & a,const FBTask & b) const { return false;}
     virtual bool operator()(const std::list<FBTask>::iterator & a,const std::list<FBTask>::iterator & b) const { return false;}
  };

  class FBTaskCompare: public TaskCompare{
    public:
    virtual bool operator()(const FBTask & a,const FBTask & b) const
    {
      bool retval = false;

      bool b_factor = b.tgt_snode_id == b.src_snode_id;


      //use the ranks first
//      if(a.rank>=0 && b.rank>=0){
//        return a.rank<b.rank;
//      }
    


      //check whether it is an update or a factorization
      bool a_factor = a.tgt_snode_id == a.src_snode_id;

      //if same type apply the usual priorities
//      if(a_factor == b_factor){

      //use the classic priorities otherwise
      if(a.tgt_snode_id>b.tgt_snode_id){
//        return true;
        retval = true;
      }
      else if(a.tgt_snode_id==b.tgt_snode_id){
//        return a.src_snode_id>b.src_snode_id;
        retval = a.src_snode_id>b.src_snode_id;
      }
      else{
//        return false;
        retval = false;
      }

      return retval;
    }
    virtual bool operator()(const std::list<FBTask>::iterator & a,const std::list<FBTask>::iterator & b) const
    {
      return (*this)(*a,*b);
//      return (*this)(*b,*a);
    }

  };




  class MCTTaskCompare: public TaskCompare{
    public:
    double cost(const FBTask & t) const{
      double retVal = 0.0;
      switch(t.type){
        case FACTOR:

        //we have to parse all the aggregates

        //then we need to compute factorization cost

        break;
        case AGGREGATE:
        break;
        case UPDATE:

        //if remote factor, then there are multiple updates to do

        break;
      }

      retVal = (double)t.rank;

      return retVal;
    }

    virtual bool operator()(const FBTask & a,const FBTask & b) const
    {
      return cost(a)<cost(b);
    }
    virtual bool operator()(const std::list<FBTask>::iterator & a,const std::list<FBTask>::iterator & b) const
    {
      return (*this)(*a,*b);
    }

  };

  class PRTaskCompare: public TaskCompare{
    public:
    virtual bool operator()(const FBTask & a,const FBTask & b) const
    {
      bool retval = false;

      //check whether it is an update or a factorization/aggregate
      bool b_factor = b.tgt_snode_id == b.src_snode_id;
      bool a_factor = a.tgt_snode_id == a.src_snode_id;

      if(a.type==FACTOR && b.type==FACTOR){
        retval= a.tgt_snode_id<b.tgt_snode_id;
      }
      else if(a.type==FACTOR){
        retval= a.tgt_snode_id<b.tgt_snode_id;
      }
      else if(b.type==FACTOR){
        retval= a.tgt_snode_id<b.tgt_snode_id;
      }
      else{
        retval= a.tgt_snode_id<b.tgt_snode_id;
      }

      return retval;
    }
    virtual bool operator()(const std::list<FBTask>::iterator & a,const std::list<FBTask>::iterator & b) const
    {
      return (*this)(*a,*b);
    }

  };


}

#endif
