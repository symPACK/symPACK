/// @file SuperNode.hpp
/// @brief TODO
/// @author Mathias Jacquelin
/// @date 2013-09-10
#ifndef _SUPERNODE_FACT_HPP_
#define _SUPERNODE_FACT_HPP_

#include "ngchol/Environment.hpp"
#include "ngchol/blas.hpp"
#include "ngchol/lapack.hpp"
#include "ngchol/NumVec.hpp"
//#include "ngchol/NumMat.hpp"
#include "ngchol/IntervalTree.hpp"
#include "ngchol/CommTypes.hpp"

#include <list>

#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) 
#define TIMER_STOP(a) 
#endif
#endif


namespace LIBCHOLESKY{



template <typename T> class NumMat;







struct NZBlockDesc{
    Int GIndex;
    size_t Offset;
    NZBlockDesc():GIndex(-1),Offset(-1){};
    NZBlockDesc(Int aGIndex, size_t aOffset):GIndex(aGIndex),Offset(aOffset){};
};


//TRACE functions
Int trc_Merge(SuperNode<double> & src_snode, SnodeUpdate &update,SuperNode<double> & tgt_snode);
Int trc_Update(SuperNode<double> & src_snode, SnodeUpdate &update, TempUpdateBuffers<double> & tmpBuffers,SuperNode<double> & tgt_snode);
Int trc_UpdateAggregate(SuperNode<double> & src_snode, SnodeUpdate &update, TempUpdateBuffers<double> & tmpBuffers, Int iTarget,SuperNode<double> & tgt_snode);
Int trc_Factorize(SuperNode<double> & tgt_snode);
Int trc_Aggregate(SuperNode<double> & src_snode,SuperNode<double> & tgt_snode);
void trc_SNodeInit(Int aiId, Int aiFc, Int aiLc, NZBlockDesc * a_block_desc, Int a_desc_cnt, double * a_nzval, Int a_nzval_cnt, SuperNode<double> & tgt_snode);

/////////////////////////////////////////
/// Class representing a supernode.
/// Class representing a supernode stored as a collection of 
/// blocks of contiguous rows in row-major format.
/////////////////////////////////////////
template<typename T>
class SuperNode{
  friend Int trc_Merge(SuperNode<double> & src_snode, SnodeUpdate &update,SuperNode<double> & tgt_snode);
  friend Int trc_Update(SuperNode<double> & src_snode, SnodeUpdate &update, TempUpdateBuffers<double> & tmpBuffers,SuperNode<double> & tgt_snode);
  friend Int trc_UpdateAggregate(SuperNode<double> & src_snode, SnodeUpdate &update, TempUpdateBuffers<double> & tmpBuffers, Int iTarget,SuperNode<double> & tgt_snode);
  friend Int trc_Factorize(SuperNode<double> & tgt_snode);
  friend Int trc_Aggregate(SuperNode<double> & src_snode,SuperNode<double> & tgt_snode);
  friend void trc_SNodeInit(Int aiId, Int aiFc, Int aiLc, NZBlockDesc * a_block_desc, Int a_desc_cnt, double * a_nzval, Int a_nzval_cnt, SuperNode<double> & tgt_snode);

  protected:

  ITree * idxToBlk_;

  inline ITree * CreateITree();

  std::vector<T> nzval_container_;
  Int nzval_cnt_;
  T * nzval_;
  std::vector<NZBlockDesc> blocks_container_;
  NZBlockDesc * blocks_;
  Int blocks_cnt_;


  bool b_own_storage_;
  Int iId_;
  Int iSize_;
  Int iFirstCol_;
  Int iLastCol_;
  

  public:

  inline void InitIdxToBlk();

  inline Int & Id(){ return iId_;}
  inline Int FirstCol(){ return iFirstCol_;}
  inline Int LastCol(){ return iLastCol_;}
  inline Int Size(){ return iSize_;}
  inline Int NZBlockCnt(){ return blocks_cnt_;}

  inline NZBlockDesc & GetNZBlockDesc(Int aiLocIndex){ return blocks_[aiLocIndex];}
  inline T* GetNZval(size_t offset){ return &nzval_[offset];}

  inline Int NRows(Int blkidx){
      NZBlockDesc & desc = GetNZBlockDesc(blkidx);
      size_t end = (blkidx<NZBlockCnt()-1)?GetNZBlockDesc(blkidx+1).Offset:nzval_cnt_;
      Int nRows = (end-desc.Offset)/Size();
      return nRows;
  }

  inline Int NRowsBelowBlock(Int blkidx){
      NZBlockDesc & desc = GetNZBlockDesc(blkidx);
      size_t end = nzval_cnt_;
      Int nRows = (end-desc.Offset)/Size();
      return nRows;
  }


  inline Int StorageSize(){ return nzval_cnt_*sizeof(T) + blocks_cnt_*sizeof(NZBlockDesc);}

  
  SuperNode();

  SuperNode(Int aiId, Int aiFc, Int aiLc, Int ai_num_rows);
  SuperNode(Int aiId, Int aiFc, Int aiLc);
  SuperNode(Int aiId, Int aiFc, Int aiLc, IntNumVec & xlindx, IntNumVec & lindx);

  SuperNode(Int aiId, Int aiFc, Int aiLc, NZBlockDesc * a_block_desc, Int a_desc_cnt,
              T * a_nzval, Int a_nzval_cnt);


  ~SuperNode();

    
  void Init(Int aiId, Int aiFc, Int aiLc, NZBlockDesc * a_block_desc, Int a_desc_cnt,
              T * a_nzval, Int a_nzval_cnt);


  inline void AddNZBlock(Int aiNRows, Int aiNCols, Int aiGIndex);


  inline Int FindBlockIdx(Int aiGIndex);
  inline void DumpITree();
  inline Int Shrink();
  
  inline bool ITreeInitialized(){return idxToBlk_->StorageSize()!=0;};

  inline void FindUpdatedFirstCol(SuperNode<T> & src_snode, Int pivot_fr, Int & tgt_fc, Int & first_pivot_idx);
  inline void FindUpdatedLastCol(SuperNode<T> & src_snode, Int tgt_fc, Int first_pivot_idx, Int & tgt_lc, Int & last_pivot_idx);

  inline Int Aggregate(SuperNode<T> & src_snode);
#ifdef COMPACT_AGGREGATES
  //this function merge structure of src_snode into the structure of the current supernode
  //right now the destination will have a pretty stupid one line per block structure
  inline Int Merge(SuperNode<T> & src_snode, SnodeUpdate &update);
#endif

  //Update an Aggregate
  inline Int UpdateAggregate(SuperNode<T> & src_snode, SnodeUpdate &update, 
              TempUpdateBuffers<T> & tmpBuffers,Int iTarget);

  //Update from a factor
  inline Int Update(SuperNode<T> & src_snode, SnodeUpdate &update, 
              TempUpdateBuffers<T> & tmpBuffers);







//  //Update from a factor
//  inline Int Update(SuperNode<T> & src_snode, Int &pivot_idx, 
//              NumMat<T> & tmpBuf,IntNumVec & src_colindx,
//                IntNumVec & src_rowindx,
//                  IntNumVec & src_to_tgt_offset, Int  pivot_fr);

  //Factorize the supernode
  inline Int Factorize();
  inline bool FindNextUpdate(SnodeUpdate & nextUpdate, const IntNumVec & Xsuper,  const IntNumVec & SupMembership,bool isLocal=true); 





};




template <typename T> inline void Serialize(Icomm & buffer,SuperNode<T> & snode, Int first_blkidx, Int first_row);
template <typename T> inline void Serialize(Icomm & buffer,SuperNode<T> & snode, Int first_blkidx, Int first_row, size_t extra_bytespace);
template <typename T> inline void Serialize(Icomm & buffer,SuperNode<T> & snode);
template <typename T> inline std::ostream& operator<<( std::ostream& os,  SuperNode<T>& snode);
template <typename T> inline size_t Deserialize(char * buffer, SuperNode<T> & snode);


#include "ngchol/SuperNode_impl.hpp"

} // namespace LIBCHOLESKY


#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) TAU_FSTART(a);
#define TIMER_STOP(a) TAU_FSTOP(a);
#endif
#endif



#endif // _SUPERNODE_FACT_HPP_
