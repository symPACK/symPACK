#ifndef _TASK_DECL_HPP_
#define _TASK_DECL_HPP_

#include "sympack/CommPull.hpp"
#include "sympack/Environment.hpp"
#include "sympack/Types.hpp"

#include <memory>
#include <functional>
#include <vector>

namespace symPACK{


namespace Factorization{
  enum class op_type {UPDATE,AGGREGATE,FACTOR,TRSM,COMM,TRSM_SEND,UPDATE2D_COMP,UPDATE2D_SEND,AGGREGATE2D_SEND,TRSM_RECV,UPDATE2D_RECV,AGGREGATE2D_RECV,UPDATE2D_SEND_OD,UPDATE2D_RECV_OD,DIAG_ENTRIES,UPDATE2D_DIAG_RECV,UPDATE2D_DIAG_SEND,FUC,BUC,FUC_D_SEND,FUC_D_RECV,FUC_DIAG_SEND,FUC_DIAG_RECV,FUC_SEND,FUC_RECV,BUC_D_SEND,BUC_D_RECV,BUC_SEND,BUC_RECV,DLT2D_COMP};
}

namespace Solve{
  enum class op_type {FUC,FU,BUC,BU,FUC_D_SEND,FUC_D_RECV,FUC_DIAG_SEND,FUC_DIAG_RECV,FUC_SEND,FUC_RECV,BUC_D_SEND,BUC_D_RECV,BUC_SEND,BUC_RECV};
}



  class IncomingMessage;

  enum TaskType {FACTOR, AGGREGATE, UPDATE};


  struct FBTask{
    TaskType type;
    Int src_snode_id;
    Int tgt_snode_id;
    //dependencies
    Int remote_deps;
    Int local_deps;
    //list of incoming messages
    std::list<IncomingMessage*> data;
    FBTask():rank(-1.0),remote_deps(0),local_deps(0){}

    //unused but preparing for task scheduling priorities
    double rank;
    void update_rank(){
      if(rank==-1.0 || 0){
        if(type==FACTOR){
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


  class CompTask{
    public:
      //byte storage for tasks 
      std::vector<char> meta;

      //dependencies
      Int remote_deps;
      Int local_deps;
      //list of incoming messages
      std::list<IncomingMessage*> data;

      CompTask( ):remote_deps(0),local_deps(0){}
  };


  class GenericTask{
    protected:
      //list of incoming messages
      std::list< std::shared_ptr<IncomingMessage> > data;
    public:
      typedef size_t id_type;
      id_type id;

      std::function< id_type(char *) > getHash;

      //dependencies
      Int remote_deps;
      Int local_deps;

      void addData(std::shared_ptr<IncomingMessage> & ptr){
        data.push_back(ptr);
      }

      void clearData(){
        data.clear();
      }

      virtual void init(){
      }

      virtual void reset(){
        clearData();
        init();
      }


      std::list<std::shared_ptr<IncomingMessage> > & getData(){
        return std::ref(data);
      }
  
      std::function< void() > execute;

      GenericTask( ):remote_deps(0),local_deps(0){}
  };

  class SparseTask: public GenericTask{
    public:
      //byte storage for tasks 
      std::vector<char> meta;

      SparseTask( ):GenericTask(){}
  };

  inline std::ostream& operator<<( std::ostream& os, const FBTask& t)
  {
    os<<"T("<<t.src_snode_id<<","<<t.tgt_snode_id<<") ";
    return os;
  }



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
    }
  };


  template<typename T = FBTask >
    class TaskCompare{
      public:
        virtual bool operator()(const T & a,const T & b) const { return false;}
        virtual bool operator()(const typename std::list<T>::iterator & a,const typename std::list<T>::iterator & b) const { return false;}
    };



  template<typename T = FBTask >
    class FBTaskCompare: public TaskCompare<T>{
      public:
        virtual bool operator()(const T & a,const T & b) const { return false;};

        virtual bool operator()(const typename std::list<T>::iterator & a,const typename std::list<T>::iterator & b) const
        {
          return (*this)(*a,*b);
        }

    };

  template<>
    inline bool FBTaskCompare<SparseTask>::operator()(const SparseTask & a,const SparseTask & b) const
    {
      bool retval = false;

      const Int * ameta = reinterpret_cast<const Int*>(a.meta.data());
      const Int asrc = ameta[0];
      const Int atgt = ameta[1];

      const Int * bmeta = reinterpret_cast<const Int*>(b.meta.data());
      const Int bsrc = bmeta[0];
      const Int btgt = bmeta[1];

      //use the classic priorities otherwise
      if(atgt>btgt){
        retval = true;
      }
      else if(atgt==btgt){
        retval = asrc>bsrc;
      }
      else{
        retval = false;
      }

      return retval;
    }




  template<>
    inline bool FBTaskCompare<CompTask>::operator()(const CompTask & a,const CompTask & b) const
    {
      bool retval = false;

      const Int * ameta = reinterpret_cast<const Int*>(a.meta.data());
      const Int asrc = ameta[0];
      const Int atgt = ameta[1];

      const Int * bmeta = reinterpret_cast<const Int*>(b.meta.data());
      const Int bsrc = bmeta[0];
      const Int btgt = bmeta[1];

      //use the classic priorities otherwise
      if(atgt>btgt){
        retval = true;
      }
      else if(atgt==btgt){
        retval = asrc>bsrc;
      }
      else{
        retval = false;
      }

      return retval;
    }





  template<>
    inline bool FBTaskCompare<FBTask>::operator()(const FBTask & a,const FBTask & b) const
    {
      bool retval = false;

      bool b_factor = b.tgt_snode_id == b.src_snode_id;
      //check whether it is an update or a factorization
      bool a_factor = a.tgt_snode_id == a.src_snode_id;

      //use the classic priorities otherwise
      if(a.tgt_snode_id>b.tgt_snode_id){
        retval = true;
      }
      else if(a.tgt_snode_id==b.tgt_snode_id){
        retval = a.src_snode_id>b.src_snode_id;
      }
      else{
        retval = false;
      }

      return retval;
    }


  class MCTTaskCompare: public TaskCompare<FBTask>{
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

  class PRTaskCompare: public TaskCompare<FBTask>{
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
#endif //_TASK_DECL_HPP_

