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
/// @file IntervalTree.hpp
/// @brief Interval Tree.
/// @author Mathias Jacquelin
/// @date 2010-09-27
#ifndef _INTTREE_DECL_HPP_
#define _INTTREE_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/LogFile.hpp"

#ifdef NO_INTRA_PROFILE
#if defined (SPROFILE)
#define SYMPACK_TIMER_START(a) 
#define SYMPACK_TIMER_STOP(a) 
#endif
#endif



//using namespace std;

namespace symPACK{

  class ITree{
    public:
      // Structure to represent an interval
      struct Interval
      {
        Int low, high;
        Int block_idx;
      };


    protected:
      // Structure to represent a node in Interval Search Tree
      class ITNode
      {
        public:
          Interval *i;  // 'i' could also be a normal variable
          Int max;
          Int min;
          Int height;
          ITNode *left, *right;
          ~ITNode(){
            delete i;

            if(left!=NULL){
              delete left;
            }

            if(right!=NULL){
              delete right;
            }
          }
          void Dump(){
            logfileptr->OFS()<< "[" << i->low << ", " << i->high << "] on "<<i->block_idx;
          }
      };


      ITNode *root_;

      // A utility function to create a new Interval Search Tree Node
      inline ITNode * newNode_(Interval & i);

      inline Int height_(ITNode *N);
      inline Int max_(ITNode *N);
      inline Int min_(ITNode *N);


      // A utility function to insert a new Interval Search Tree Node
      // This is similar to BST Insert.  Here the low value of interval
      // is used tomaintain BST property
      virtual inline ITNode *insert_(ITNode *root, Interval & i);
      inline Int recomputeMax_(ITNode * root);
      inline Int recomputeMinMax_(ITNode * root);

      // A utility function to check if given two intervals overlap
      inline bool doOVerlap_(const Interval &i1,const Interval &i2);
      inline bool doOVerlap_(const ITree::Interval &i1, const Int & low, const Int & high);
      inline Interval *intervalSearch_(ITNode *root, const Int & begin, const Int & end);
      inline Interval *intervalSearch_(ITNode *root, const Int & begin, const Int &end, Interval * & closestR, Interval * & closestL);

      inline void inorder_(ITNode *root);

      inline Int getSize_(ITNode *root);


    public:


      ITree():root_(NULL) { }


      virtual ~ITree()
      {
        if(root_!=NULL){
          delete root_;
        }
      }


      virtual inline void Rebalance(){};

      inline void Dump()
      {
        logfileptr->OFS()<<"Height of tree is "<<height_(root_)<<std::endl;
        inorder_(root_);
      }


      inline void Insert(Interval & i)
      {
        SYMPACK_TIMER_START(ITREE_INSERT);
        root_ = insert_(root_,i);
        SYMPACK_TIMER_STOP(ITREE_INSERT);
      }

      inline void Insert(Int i)
      {
        Interval it;  
        it.low = i;
        it.high = i;
        Insert(it);
      }


      inline Interval * IntervalSearch(const Interval & i){
        return intervalSearch_(root_,i.low,i.high);
      }

      inline Interval * IntervalSearch(const Int & low,const Int & high){
        return intervalSearch_(root_,low,high);
      }

      inline Interval * IntervalSearch(const Int & low,const Int & high, Interval * & closestR, Interval * & closestL){
        return intervalSearch_(root_,low,high,closestR,closestL);
      }

      inline Int StorageSize(){
        Int size = 0;
        size = getSize_(root_);
        return size;
      }
  };


  class AVLITree: public ITree{
    public:

    protected:
      // A utility function to right rotate subtree rooted with y
      // See the diagram given above.
      inline ITNode * rightRotate_(ITNode *y);

      inline ITNode * leftRotate_(ITNode *x);
      // Get Balance factor of node N
      inline Int getBalance_(ITNode *N);

      // A utility function to insert a new Interval Search Tree Node
      // This is similar to BST Insert.  Here the low value of interval
      // is used tomaintain BST property
      virtual inline ITNode *insert_(ITNode *root, Interval & i);


    public:


      AVLITree():ITree() {
      }


      virtual ~AVLITree() {
      }
  };

  class DSWITree: public ITree{
    public:

    protected:
      inline Int fullSize_ ( Int size );
      inline void tree_to_vine_ (ITNode * root, Int &size);
      inline void compression_  (ITNode * root, Int count);
      inline void vine_to_tree_ (ITNode * root, Int n);
      inline void correctTree_  (ITNode * root);


    public:


      DSWITree():ITree() { }


      virtual ~DSWITree() {
      }

      virtual inline void Rebalance();
  };


#include "sympack/IntervalTree_impl.hpp"




  namespace UnitTest{
    bool ITree_Test();
  }

}



#ifdef NO_INTRA_PROFILE
#if defined (SPROFILE)
#define SYMPACK_TIMER_START(a) SYMPACK_FSTART(a);
#define SYMPACK_TIMER_STOP(a) SYMPACK_FSTOP(a);
#endif
#endif



#endif //_INTTREE_DECL_HPP_
