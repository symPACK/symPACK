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
/// @file SuperNode.hpp
/// @brief TODO
/// @author Mathias Jacquelin
/// @date 2013-09-10
#ifndef _SUPERNODE_FACT_HPP_
#define _SUPERNODE_FACT_HPP_

#include "sympack/Environment.hpp"
#include "sympack/blas.hpp"
#include "sympack/lapack.hpp"
//#include "sympack/NumMat.hpp"
#include "sympack/IntervalTree.hpp"
#include "sympack/CommTypes.hpp"

#include <upcxx.h>
#include <list>
#include <mutex>
#include <atomic>

#ifdef NO_INTRA_PROFILE
#if defined (SPROFILE)
#define SYMPACK_TIMER_START(a) 
#define SYMPACK_TIMER_STOP(a) 
#define scope_timer(b,a)
#endif
#endif


//#define _INDEFINITE_

namespace symPACK{

  struct NZBlockDesc{
    bool Last;
    Int GIndex;
    size_t Offset;
    NZBlockDesc():GIndex(-1),Offset(-1),Last(false){};
    NZBlockDesc(Int aGIndex, size_t aOffset):GIndex(aGIndex),Offset(aOffset),Last(false){};
  };



  struct SuperNodeDesc{
    bool b_own_storage_;
    Int iId_;
    Int iSize_;
    Int iFirstRow_;
    Int iFirstCol_;
    Int iLastCol_;
    Int iN_; 
    Int blocks_cnt_;
    Int nzval_cnt_;
    //Int updrows_cnt_;

    SuperNodeDesc():iId_(-1),iSize_(-1),iFirstRow_(-1),iFirstCol_(-1),
        iLastCol_(-1),iN_(-1),blocks_cnt_(-1),nzval_cnt_(-1)/*,updrows_cnt_(-1)*/{}
  };


  template<typename T>
    class SuperNodeBase{
      public:
        std::mutex lock_;
        virtual ~SuperNodeBase(){};
        virtual inline Int & Id() = 0;
    };

  ////////////////////////////////////////
  /// Class representing a supernode.
  /// Class representing a supernode stored as a collection of 
  /// blocks of contiguous rows in row-major format.
  /////////////////////////////////////////
  template<typename T, class Allocator = UpcxxAllocator>
    class SuperNode: public SuperNodeBase<T>{
      public:
        std::atomic<int> trsm_count;
        //Idx * updrows_;

      protected:


#ifndef ITREE
        Int iLastRow_;
        std::vector<Int> * globalToLocal_;
#else
        ITree<Int> * idxToBlk_;
#endif

        //actual storage
        char * storage_container_;
        char * loc_storage_container_;
        size_t storage_size_;

        //utility pointers
        SuperNodeDesc * meta_;
        T * nzval_;
        NZBlockDesc * blocks_;

      protected:
        inline ITree<Int> * CreateITree();

      public:

        inline void InitIdxToBlk();
        inline Int & Id(){ return meta_->iId_;}
        inline bool OwnDiagonal(){ return FirstRow() == FirstCol();}
        inline Int FirstRow(){ return meta_->iFirstRow_;}
        inline Int FirstCol(){ return meta_->iFirstCol_;}
        inline Int LastCol(){ return meta_->iLastCol_;}
        inline Int Size(){ return meta_->iSize_;}
        inline Int N(){ return meta_->iN_;}
        inline Ptr NNZ(){ return meta_->nzval_cnt_;}
        inline Int NZBlockCnt(){ return meta_->blocks_cnt_;}
        inline NZBlockDesc & GetNZBlockDesc(Int aiLocIndex){ return *(blocks_ -aiLocIndex);}
        inline T* GetNZval(size_t offset){ return &nzval_[offset];}
        inline SuperNodeDesc * GetMeta(){return meta_;}
        inline Int NRows(Int blkidx){
          NZBlockDesc & desc = GetNZBlockDesc(blkidx);
          size_t end = (blkidx<NZBlockCnt()-1)?GetNZBlockDesc(blkidx+1).Offset:meta_->nzval_cnt_;
          Int nRows = (end-desc.Offset)/Size();
          return nRows;
        }

        inline Int NRowsBelowBlock(Int blkidx){
          NZBlockDesc & desc = GetNZBlockDesc(blkidx);
          size_t end = meta_->nzval_cnt_;
          Int nRows = (end-desc.Offset)/Size();
          return nRows;
        }


        inline Int StorageSize(){ return storage_size_;}

        char * GetStoragePtr(Int row){
          //TODO fix this
          return loc_storage_container_;
        }

        SuperNode();
        SuperNode(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices);
        SuperNode(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int ai_num_rows, Int aiN, Int aiNZBlkCnt=-1);
        SuperNode(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int aiN);
        SuperNode(char * storage_ptr,size_t storage_size, Int GIndex = -1);
        ~SuperNode();

        virtual void Init(char * storage_ptr,size_t storage_size, Int GIndex = -1);
        virtual void Init(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices);

        virtual inline void AddNZBlock(Int aiNRows, Int aiNCols, Int aiGIndex);
        inline void Reserve(size_t storage_size);

        inline Int FindBlockIdx(Int aiGIndex);
        inline Int FindBlockIdx(Int aiGIndex,Int & closestR, Int & closestL);
        inline Int FindBlockIdx(Int fr, Int lr, ITree<Int>::Interval<Int> & overlap);

        inline void DumpITree();
        virtual inline Int Shrink();

#ifdef ITREE
        inline bool ITreeInitialized(){return idxToBlk_->StorageSize()!=0;};
#endif

        inline void FindUpdatedFirstCol(SuperNode<T,Allocator> * src_snode, Int pivot_fr, Int & tgt_fc, Int & first_pivot_idx);
        inline void FindUpdatedLastCol(SuperNode<T,Allocator> * src_snode, Int tgt_fc, Int first_pivot_idx, Int & tgt_lc, Int & last_pivot_idx);

        virtual inline Int Aggregate(SuperNode<T,Allocator> * src_snode);
        //this function merge structure of src_snode into the structure of the current supernode
        //right now the destination will have a pretty stupid one line per block structure
        virtual inline Int Merge(SuperNode<T,Allocator> * src_snode, SnodeUpdate &update);

        //Update an Aggregate
        virtual inline Int UpdateAggregate(SuperNode<T,Allocator> * src_snode, SnodeUpdate &update, 
            TempUpdateBuffers<T> & tmpBuffers,Int iTarget, Int iam);

        //Update from a factor
        virtual inline Int Update(SuperNode<T,Allocator> * src_snode, SnodeUpdate &update, 
            TempUpdateBuffers<T> & tmpBuffers);

        //Factorize the supernode
        virtual inline Int Factorize(TempUpdateBuffers<T> & tmpBuffers);
        virtual inline Int Factorize(SuperNode<T,Allocator> * diag_snode, TempUpdateBuffers<T> & tmpBuffers);
        virtual inline Int Factorize_diag(TempUpdateBuffers<T> & tmpBuffers);
        virtual inline Int Factorize_TRSM(SuperNode<T,Allocator> * diag_snode, Int blkidx);

        inline bool FindNextUpdate(SnodeUpdate & nextUpdate, const std::vector<Int> & Xsuper,  const std::vector<Int> & SupMembership,bool isLocal=true); 


        //forward and backward solve phases
        virtual inline void forward_update_contrib( T * RHS, SuperNode<T> * cur_snode, std::vector<Int> & perm);
        virtual inline void forward_update_contrib( SuperNode<T> * cur_snode, Int nrhsOffset = 0,Int pnrhs=-1);
        virtual inline void forward_update(SuperNode<T,Allocator> * src_contrib, Int iOwner,Int iam, Int nrhsOffset = 0,Int pnrhs=-1);

        virtual inline void back_update_contrib(SuperNode<T> * cur_snode, Int nrhsOffset = 0,Int pnrhs=-1);
        virtual inline void back_update(SuperNode<T,Allocator> * src_contrib, Int nrhsOffset = 0,Int pnrhs=-1);


        virtual inline void Serialize(Icomm & buffer, Int first_blkidx=0, Idx first_row=0);
        virtual inline size_t Deserialize(char * buffer, size_t size);


#ifdef SP_THREADS
        std::atomic<bool> in_use;
        
#endif


    };



  template <typename T, class Allocator> inline void Serialize(Icomm & buffer,SuperNode<T,Allocator> & snode, Int first_blkidx=0, Idx first_row=0);
  template <typename T, class Allocator> inline std::ostream& operator<<( std::ostream& os,  SuperNode<T,Allocator>& snode);
  template <typename T, class Allocator> inline size_t Deserialize(char * buffer, size_t size, SuperNode<T,Allocator> & snode);







  inline std::ostream& operator<<( std::ostream& os,  SuperNodeDesc& desc);
  inline std::ostream& operator<<( std::ostream& os,  NZBlockDesc& block);




} // namespace SYMPACK


#include "sympack/impl/SuperNode_impl.hpp"




#ifdef NO_INTRA_PROFILE
#if defined (SPROFILE)
#define SYMPACK_TIMER_START(a) SYMPACK_FSTART(a);
#define SYMPACK_TIMER_STOP(a) SYMPACK_FSTOP(a);
#define scope_timer(b,a) symPACK_scope_timer b(#a)
#endif
#endif



#endif // _SUPERNODE_FACT_HPP_
