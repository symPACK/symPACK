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

  template<typename F = Int>
  class ITree{
    public:
      // Structure to represent an interval
      template<typename T = F>
      struct Interval
      {
        Int low, high;
        T data;
      };


    protected:
      // Structure to represent a node in Interval Search Tree
      template<typename T = F>
      class ITNode
      {
        public:
          Interval<T> *i;  // 'i' could also be a normal variable
          Int max;
          Int min;
          Int height;
          ITNode<T> *left, *right;
          ~ITNode(){
            delete i;

            if(left!=nullptr){
              delete left;
            }

            if(right!=nullptr){
              delete right;
            }
          }
          void Dump(){
            logfileptr->OFS()<< "[" << i->low << ", " << i->high << "] on "<<i->data;//block_idx;
          }
      };


      ITNode<F> *root_;

      // A utility function to create a new Interval Search Tree Node
      inline ITNode<F> * newNode_(Interval<F> & i);

      inline Int height_(ITNode<F> *N);
      inline Int max_(ITNode<F> *N);
      inline Int min_(ITNode<F> *N);


      // A utility function to insert a new Interval Search Tree Node
      // This is similar to BST Insert.  Here the low value of interval
      // is used tomaintain BST property
      virtual inline ITNode<F> *insert_(ITNode<F> *root, Interval<F> & i);
      inline Int recomputeMax_(ITNode<F> * root);
      inline Int recomputeMinMax_(ITNode<F> * root);

      // A utility function to check if given two intervals overlap
      inline bool doOVerlap_(const Interval<F> &i1,const Interval<F> &i2);
      inline bool doOVerlap_(const ITree::Interval<F> &i1, const Int & low, const Int & high);
      inline Interval<F> *intervalSearch_(ITNode<F> *root, const Int & begin, const Int & end);
      inline Interval<F> *intervalSearch_(ITNode<F> *root, const Int & begin, const Int &end, Interval<F> * & closestR, Interval<F> * & closestL);

      inline void inorder_(ITNode<F> *root);

      inline Int getSize_(ITNode<F> *root);


    public:


      ITree():root_(nullptr) { }


      virtual ~ITree()
      {
        if(root_!=nullptr){
          delete root_;
        }
      }


      virtual inline void Rebalance(){};

      inline void Dump()
      {
        logfileptr->OFS()<<"Height of tree is "<<height_(root_)<<std::endl;
        inorder_(root_);
      }


      inline void Insert(Interval<F> & i)
      {
        SYMPACK_TIMER_START(ITREE_INSERT);
        root_ = insert_(root_,i);
        SYMPACK_TIMER_STOP(ITREE_INSERT);
      }

      inline void Insert(Int i)
      {
        Interval<F> it;  
        it.low = i;
        it.high = i;
        Insert(it);
      }


      inline Interval<F> * IntervalSearch(const Interval<F> & i){
        return intervalSearch_(root_,i.low,i.high);
      }

      inline Interval<F> * IntervalSearch(const Int & low,const Int & high){
        return intervalSearch_(root_,low,high);
      }

      inline Interval<F> * IntervalSearch(const Int & low,const Int & high, Interval<F> * & closestR, Interval<F> * & closestL){
        return intervalSearch_(root_,low,high,closestR,closestL);
      }

      inline Int StorageSize(){
        Int size = 0;
        size = getSize_(root_);
        return size;
      }
  };


  template<typename F = Int>
  class AVLITree: public ITree<F>{
    public:

    protected:
      // A utility function to right rotate subtree rooted with y
      // See the diagram given above.
      inline typename ITree<F>::template ITNode<F> * rightRotate_(typename ITree<F>::template ITNode<F> *y);

      inline typename ITree<F>::template ITNode<F> * leftRotate_(typename ITree<F>::template ITNode<F> *x);
      // Get Balance factor of node N
      inline Int getBalance_(typename ITree<F>::template ITNode<F> *N);

      // A utility function to insert a new Interval Search Tree Node
      // This is similar to BST Insert.  Here the low value of interval
      // is used tomaintain BST property
      virtual inline typename ITree<F>::template ITNode<F> *insert_(typename ITree<F>::template ITNode<F> *root, typename ITree<F>::template Interval<F> & i);


    public:


      AVLITree():ITree<F>() {
      }


      virtual ~AVLITree() {
      }
  };

  template<typename F = Int>
  class DSWITree: public ITree<F>{
    public:

    protected:
      inline Int fullSize_ ( Int size );
      inline void tree_to_vine_ (typename ITree<F>::template ITNode<F> * root, Int &size);
      inline void compression_  (typename ITree<F>::template ITNode<F> * root, Int count);
      inline void vine_to_tree_ (typename ITree<F>::template ITNode<F> * root, Int n);
      inline void correctTree_  (typename ITree<F>::template ITNode<F> * root);


    public:


      DSWITree():ITree<F>() { }


      virtual ~DSWITree() {
      }

      virtual inline void Rebalance();
  };


#include "sympack/impl/IntervalTree_impl.hpp"




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
