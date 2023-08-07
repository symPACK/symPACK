#ifndef _SUPERNODE_FACT_HPP_
#define _SUPERNODE_FACT_HPP_

#include "sympack/Environment.hpp"
#include "sympack/blas.hpp"
#include "sympack/lapack.hpp"
//#include "sympack/NumMat.hpp"
#include "sympack/IntervalTree.hpp"
#include "sympack/CommTypes.hpp"
#ifdef CUDA_MODE
#include "sympack/cuBLAS.hpp"
#endif
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
    NZBlockDesc():GIndex(-1),Offset((size_t)-1),Last(false){};
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
    Int panel_sz_;
    //Int updrows_cnt_;

    SuperNodeDesc():iId_(-1),iSize_(-1),iFirstRow_(-1),iFirstCol_(-1),
        iLastCol_(-1),iN_(-1),panel_sz_(-1),blocks_cnt_(-1),nzval_cnt_(-1)/*,updrows_cnt_(-1)*/{}
  };


  template<typename T>
    class SuperNodeBase{
      public:
#ifdef SP_THREADS
        std::atomic<bool> in_use;
#endif

#ifdef SP_THREADS
        SuperNodeBase():
          in_use(false)
        {}
#else
        SuperNodeBase()
        {}
#endif
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
          size_t nzvalcnt = end - desc.Offset;
          //Int panel = meta_->panel_sz_;
          //if ( panel>0 && desc.GIndex == FirstCol() ) {
          //  gdb_lock();
          //  Int numPanels = std::ceil((double)meta_->iSize_/(double)panel);
          //  nzvalcnt += numPanels*(numPanels-1)*panel*panel/2; 
          //}

          Int nRows = (nzvalcnt)/Size();
          return nRows;
        }

        inline Int NRowsBelowBlock(Int blkidx){
          NZBlockDesc & desc = GetNZBlockDesc(blkidx);
          size_t end = meta_->nzval_cnt_;
          size_t nzvalcnt = end - desc.Offset;
          Int panel = meta_->panel_sz_;
          if ( panel>0 && desc.GIndex == FirstCol() ) {
            Int numPanels = std::ceil((double)meta_->iSize_/(double)panel);
            nzvalcnt += numPanels*(numPanels-1)*panel*panel/2; 
          }
          Int nRows = (nzvalcnt)/Size();
          return nRows;
        }
 
        virtual inline void clear(){ 
          if(NNZ()>0){
            std::fill(nzval_,nzval_+NNZ(),T(0));
          }
        }



        inline Int StorageSize(){ return storage_size_;}

        char * GetStoragePtr(Int row){
          //TODO fix this
          return loc_storage_container_;
        }

        SuperNode();
        SuperNode(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices, Int panel=-1);
        SuperNode(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int ai_num_rows, Int aiN, Int aiNZBlkCnt=-1, Int panel=-1);
        SuperNode(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int aiN, Int panel=-1);
        SuperNode(char * storage_ptr,size_t storage_size, Int GIndex = -1);
        virtual ~SuperNode();

        virtual void Init(char * storage_ptr,size_t storage_size, Int GIndex = -1);
        virtual void Init(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices, Int panel=-1);

        virtual inline void AddNZBlock(Int aiNRows, Int aiGIndex);
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

        inline bool FindNextUpdate(SnodeUpdate & nextUpdate, const std::vector<Int> & Xsuper,  const std::vector<Int> & SupMembership,bool isLocal=true); 


        //forward and backward solve phases
        virtual inline void forward_update_contrib( T * RHS, SuperNode<T> * cur_snode, std::vector<Int> & perm);
        virtual inline void forward_update_contrib( SuperNode<T> * cur_snode, Int nrhsOffset = 0,Int pnrhs=-1);
        virtual inline void forward_update(SuperNode<T,Allocator> * src_contrib, Int iOwner,Int iam, Int nrhsOffset = 0,Int pnrhs=-1);

        virtual inline void back_update_contrib(SuperNode<T> * cur_snode, Int nrhsOffset = 0,Int pnrhs=-1);
        virtual inline void back_update(SuperNode<T,Allocator> * src_contrib, Int nrhsOffset = 0,Int pnrhs=-1);


        virtual inline void Serialize(Icomm & buffer, Int first_blkidx=0, Idx first_row=0);
        virtual inline size_t Deserialize(char * buffer, size_t size);



    };



  template <typename T, class Allocator> inline void Serialize(Icomm & buffer,SuperNode<T,Allocator> & snode, Int first_blkidx=0, Idx first_row=0);
  template <typename T, class Allocator> inline std::ostream& operator<<( std::ostream& os,  SuperNode<T,Allocator>& snode);
  template <typename T, class Allocator> inline size_t Deserialize(char * buffer, size_t size, SuperNode<T,Allocator> & snode);







  inline std::ostream& operator<<( std::ostream& os,  SuperNodeDesc& desc);
  inline std::ostream& operator<<( std::ostream& os,  NZBlockDesc& block);


  template<typename T>
class supernode_lock {
  public:
    supernode_lock(SuperNodeBase<T> * snode) {
      snode_ = snode;
    }

    ~supernode_lock(){
#ifdef SP_THREADS
      if(Multithreading::NumThread>1){
        snode_->in_use = false;
      }
#endif
    }
  protected:
    SuperNodeBase<T> * snode_;
};



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
