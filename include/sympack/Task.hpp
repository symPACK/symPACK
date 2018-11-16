/*
	 Copyright (c) 2016 The Regents of the University of California,
	 through Lawrence Berkeley National Laboratory.  

   Author: Mathias Jacquelin
	 
   This file is part of symPACK. All rights reserved.

	 Redistribution and use in source and binary forms, with or without
	 modification, are permitted provided that the following conditions are met:

	 (1) Redistributions of source code must retain the above copyright notice, this
	 list of conditions and the following disclaimer.
	 (2) Redistributions in binary form must reproduce the above copyright notice,
	 this list of conditions and the following disclaimer in the documentation
	 and/or other materials provided with the distribution.
	 (3) Neither the name of the University of California, Lawrence Berkeley
	 National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
	 be used to endorse or promote products derived from this software without
	 specific prior written permission.

	 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
	 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
	 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
	 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
	 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

	 You are under no obligation whatsoever to provide any bug fixes, patches, or
	 upgrades to the features, functionality or performance of the source code
	 ("Enhancements") to anyone; however, if you choose to make your Enhancements
	 available either publicly, or directly to Lawrence Berkeley National
	 Laboratory, without imposing a separate written license agreement for such
	 Enhancements, then you hereby grant the following license: a non-exclusive,
	 royalty-free perpetual license to install, use, modify, prepare derivative
	 works, incorporate into other computer software, distribute, and sublicense
	 such enhancements or derivative works thereof, in binary and source code form.
*/
#ifndef _TASK_DECL_HPP_
#define _TASK_DECL_HPP_

#include "sympack/Environment.hpp"

namespace symPACK{

  class IncomingMessage;

  enum TaskType {FACTOR, AGGREGATE, UPDATE/*, ALLOCATE*/};


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
          //taskit->rank = 2.0*pow((double)src_snode.Size(),3.0)/3.0 + src_snode.Size()*(src_snode.Size()+1.0)*(src_snode.NRowsBelowBlock(0)-src_snode.Size())/2.0;
          rank = 4.0;
        }
        else if(type == AGGREGATE){
          rank = 2.0;
        }
        else if(type == UPDATE){
          rank = 1.0;


//          SuperNode2<T> * cur_src_snode; 
//          IncomingMessage * msgPtr = nullptr;
//          //Local or remote factor
//          //we have only one local or one remote incoming aggregate
//          if(curTask.data.size()==0){
//            cur_src_snode = snodeLocal(curTask.src_snode_id);
//          }
//          else{
//            auto msgit = curTask.data.begin();
//            msgPtr = *msgit;
//            assert(msgPtr->IsDone());
//            char* dataPtr = msgPtr->GetLocalPtr();
//
//            cur_src_snode = new SuperNode2<T>(dataPtr,msgPtr->Size(),msgPtr->meta.GIndex);
//            cur_src_snode->InitIdxToBlk();
//          }
//
//
//          //parse cur_src_snode to compute cost
//          while(cur_src_snode->FindNextUpdate(curUpdate,Xsuper_,SupMembership_,iam==iSrcOwner)){
//            //skip if this update is "lower"
//            if(curUpdate.tgt_snode_id<curTask.tgt_snode_id){continue;}
//
//            Int iUpdater = this->Mapping_->Map(curUpdate.tgt_snode_id-1,cur_src_snode->Id()-1);
//            if(iUpdater == iam){
//            }
//          }
        }
      }
    }
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







/************************ UNUSED EXPERIMENTAL CLASS ***********************/

//class Task{
//  public:
//    static const size_t idSize = 100;
//    std::vector<char> id;
//    //dependencies
//    std::list< upcxx::global_ptr<Task> > inputs;
//    std::list< upcxx::global_ptr<Task> > outputs;
//    //data
//    upcxx::global_ptr<char> data;
//
//    //computation
//    std::function<int(std::vector<char>&)> process;
//
//    Task():data(nullptr){
//      id.resize(idSize);
//    }
//};

}
#endif //_TASK_DECL_HPP_

