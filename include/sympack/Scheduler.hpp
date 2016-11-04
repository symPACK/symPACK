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
#ifndef _SCHEDULER_DECL_HPP_
#define _SCHEDULER_DECL_HPP_

#include "sympack/Types.hpp"
#include "sympack/Task.hpp"
#include <list>
#include <queue>


namespace symPACK{

  template <class T >
    class Scheduler{
      public:
        Scheduler(){
        }

        virtual T& top() =0;
        virtual void pop() =0;

        virtual void push( const T& val) =0;
        virtual unsigned int size() =0;

        //    virtual const std::priority_queue<T, std::vector<T>, TaskCompare >  GetQueue() =0;

        virtual bool done() =0;

      protected:
    };




  template <class T >
    class DLScheduler: public Scheduler<T>{
      public:
        DLScheduler(){
        }

        virtual T& top(){ 
          T & ref = const_cast<T&>(readyTasks_.top());
          return ref;
        }
        virtual void pop() {
          readyTasks_.pop();
        }

        virtual void push( const T& val) { 
          readyTasks_.push(val);
        }
        virtual unsigned int size() {
          int count = readyTasks_.size();
          return count;
        }

        virtual bool done(){ 
          bool empty = readyTasks_.empty();
          return empty;
        }



        const std::priority_queue<T, std::vector<T>, FBTaskCompare > & GetQueue() {return  readyTasks_;}
      protected:
        std::priority_queue<T, std::vector<T>, FBTaskCompare > readyTasks_;
    };


  template <class T >
    class MCTScheduler: public Scheduler<T>{
      public:
        MCTScheduler(){
        }

        virtual T& top(){
          T & ref = const_cast<T&>(readyTasks_.top());
          return ref;
        }
        virtual void pop() { 
          readyTasks_.pop();
        }

        virtual void push( const T& val) { 
          readyTasks_.push(val);
        }



        virtual unsigned int size() {
          int count = readyTasks_.size();
          return count;
        }

        virtual bool done(){ 
          bool empty = readyTasks_.empty();
          return empty;
        }








        const std::priority_queue<T, std::vector<T>, MCTTaskCompare > & GetQueue(){return  readyTasks_;}




        //    virtual void compute_costs(PtrVec & Xlindx, IdxVec & Lindx){
        //
        //
        //        //I is source supernode
        //        Int I = ;
        //
        //            Idx fc = Xsuper_[I-1];
        //            Int width = Xsuper_[I] - Xsuper_[I-1];
        //            Int height = cc_[fc-1];
        //            //cost of factoring curent panel
        //            double local_load = factor_cost(height,width);
        //
        //            //cost of updating ancestors
        //            {
        //              Ptr fi = Xlindx_[I-1];
        //              Ptr li = Xlindx_[I]-1;
        //
        //                  Int m = li-tgt_fr_ptr+1;
        //                  Int n = width;
        //                  Int k = tgt_lr_ptr - tgt_fr_ptr+1; 
        //                  double cost = update_cost(m,n,k);
        //
        //              //parse rows
        //              Int tgt_snode_id = I;
        //              Idx tgt_fr = fc;
        //              Idx tgt_lr = fc;
        //              Ptr tgt_fr_ptr = 0;
        //              Ptr tgt_lr_ptr = 0;
        //              for(Ptr rowptr = fi; rowptr<=li;rowptr++){
        //                Idx row = Lindx_[rowptr-1]; 
        //                if(SupMembership_[row-1]==tgt_snode_id){ continue;}
        //
        //                //we have a new tgt_snode_id
        //                tgt_lr_ptr = rowptr-1;
        //                tgt_lr = Lindx_[rowptr-1 -1];
        //                if(tgt_snode_id !=I){
        //                  Int m = li-tgt_fr_ptr+1;
        //                  Int n = width;
        //                  Int k = tgt_lr_ptr - tgt_fr_ptr+1; 
        //                  double cost = update_cost(m,n,k);
        //                  break;
        //                }
        //                tgt_fr = row;
        //                tgt_fr_ptr = rowptr;
        //                tgt_snode_id = SupMembership_[row-1];
        //              }
        //              //do the last update supernode
        //              //we have a new tgt_snode_id
        //              tgt_lr_ptr = li;
        //              tgt_lr = Lindx_[li -1];
        //              if(tgt_snode_id!=I){
        //                Int m = li-tgt_fr_ptr+1;
        //                Int n = width;
        //                Int k = tgt_lr_ptr - tgt_fr_ptr+1; 
        //                double cost = update_cost(m,n,k);
        //                  if(fan_in_){
        //                    SubTreeLoad[I]+=cost;
        //                    NodeLoad[I]+=cost;
        //                  }
        //                  else{
        //                    SubTreeLoad[tgt_snode_id]+=cost;
        //                    NodeLoad[tgt_snode_id]+=cost;
        //                  }
        //                }
        //            }
        //
        //
        //    }
        //




      protected:
        std::priority_queue<T, std::vector<T>, MCTTaskCompare > readyTasks_;
    };


  template <class T >
    class PRScheduler: public Scheduler<T>{
      public:
        PRScheduler(){
        }

        virtual T& top(){
          T & ref = const_cast<T&>(readyTasks_.top());
          return ref;
        }
        virtual void pop() {
          readyTasks_.pop();
        }

        virtual void push( const T& val) { 
          readyTasks_.push(val);
        }


        virtual unsigned int size() {
          int count = readyTasks_.size();
          return count;
        }

        virtual bool done(){ 
          bool empty = readyTasks_.empty();
          return empty;
        }





        const std::priority_queue<T, std::vector<T>, PRTaskCompare > & GetQueue(){return  readyTasks_;}
      protected:
        std::priority_queue<T, std::vector<T>, PRTaskCompare > readyTasks_;
    };


  template <class T >
    class FIFOScheduler: public Scheduler<T>{
      public:
        FIFOScheduler(){
        }

        virtual T& top(){ 
          T & ref = const_cast<T&>(readyTasks_.front());
          return ref;
        }
        virtual void pop() { 
          readyTasks_.pop();
        }

        virtual void push( const T& val) { 
          readyTasks_.push(val);
        }
        virtual unsigned int size() {
          int count = readyTasks_.size();
          return count;
        }

        virtual bool done(){ 
          bool empty = readyTasks_.empty();
          return empty;
        }
        const std::queue<T> & GetQueue(){return  readyTasks_;}
      protected:
        std::queue<T> readyTasks_;
    };





}
#endif // _SCHEDULER_DECL_HPP_
