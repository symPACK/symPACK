/// @file SuperNode.hpp
/// @brief TODO
/// @author Mathias Jacquelin
/// @date 2013-09-10
#ifndef _SUPERNODE_FACT_HPP_
#define _SUPERNODE_FACT_HPP_

#include "Environment.hpp"
#include "NumVec.hpp"
#include "NZBlock.hpp"

namespace LIBCHOLESKY{

struct NZBlockDesc{
    Int GIndex;
    size_t Offset;
    NZBlockDesc():GIndex(-1),Offset(-1){};
    NZBlockDesc(Int aGIndex, size_t aOffset):GIndex(aGIndex),Offset(aOffset){};
};




template<typename T>
class SuperNode{

  protected:

  //store the nnz values of the nz blocks in a contiguous row major buffer
  std::vector<Int> global_to_local_index_;
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


  inline Int Id(){ return iId_;}
  inline Int FirstCol(){ return iFirstCol_;}
  inline Int LastCol(){ return iLastCol_;}
  inline Int Size(){ return iSize_;}
  inline Int NZBlockCnt(){ return blocks_cnt_;}

  inline NZBlockDesc & GetNZBlockDesc(Int aiLocIndex){ return blocks_[aiLocIndex];}
  inline T* GetNZval(size_t offset){ return &nzval_[offset];}

  /*inline*/ Int NRows(Int blkidx){
      NZBlockDesc & desc = GetNZBlockDesc(blkidx);
      size_t end = (blkidx<NZBlockCnt()-1)?GetNZBlockDesc(blkidx+1).Offset:nzval_cnt_;
      Int nRows = (end-desc.Offset)/Size();
      return nRows;
  }

  inline Int NRowsBelowBlock(Int blkidx){
      NZBlockDesc & desc = GetNZBlockDesc(blkidx);
      size_t end = nzval_cnt_;
      Int nRows = (end-desc.Offset)/Size();
  }


  inline Int StorageSize(){ return nzval_cnt_*sizeof(T) + blocks_cnt_*sizeof(NZBlockDesc);}

  SuperNode(Int aiId, Int aiFc, Int aiLc, Int aiN, Int ai_num_rows) :iId_(aiId),iFirstCol_(aiFc),iLastCol_(aiLc) {
     //this is an upper bound
     assert(ai_num_rows>0);

     //compute supernode size / width
     iSize_ = iLastCol_ - iFirstCol_+1;

     nzval_cnt_ = 0;
     nzval_container_.reserve(iSize_*ai_num_rows);
     nzval_ = NULL;

     blocks_container_.reserve(ai_num_rows);
     blocks_ = NULL;
     blocks_cnt_ = 0;

     global_to_local_index_.resize(aiN,-1);


     b_own_storage_ = true;
  }; 










  SuperNode(Int aiId, Int aiFc, Int aiLc, Int aiN, NZBlockDesc * a_block_desc, Int a_desc_cnt,
              T * a_nzval, Int a_nzval_cnt) :
                                  iId_(aiId),iFirstCol_(aiFc), iLastCol_(aiLc)
                                                 {
     //compute supernode size / width
     iSize_ = iLastCol_ - iFirstCol_+1;

    b_own_storage_ = false;
    nzval_= a_nzval;
    nzval_cnt_ = a_nzval_cnt;

    assert(nzval_cnt_ % iSize_ == 0);


    blocks_ = a_block_desc;
    blocks_cnt_ = a_desc_cnt;
    global_to_local_index_.resize(aiN,-1);

    //restore 0-based offsets and compute global_to_local indices
    for(Int blkidx=blocks_cnt_-1; blkidx>=0;--blkidx){
      blocks_[blkidx].Offset -= blocks_[0].Offset;
    }

    for(Int blkidx=0; blkidx<blocks_cnt_;++blkidx){
      for(Int rowidx = 0; rowidx< NRows(blkidx); ++rowidx){
        global_to_local_index_[blocks_[blkidx].GIndex+rowidx-1] = blkidx;
      }
    }


    
  }


  ~SuperNode(){
  }
    

 

 inline void AddNZBlock(Int aiNRows, Int aiNCols, Int aiGIndex){

    //Resize the container if I own the storage
    if(b_own_storage_){
      blocks_container_.push_back(NZBlockDesc(aiGIndex,nzval_cnt_));
      for(Int rowidx = 0; rowidx< aiNRows; ++rowidx){
        global_to_local_index_[aiGIndex+rowidx-1] = blocks_cnt_;
      }
      blocks_cnt_++;
      nzval_cnt_+=aiNRows*iSize_;
      nzval_container_.resize(nzval_cnt_,ZERO<T>());
      nzval_ = &nzval_container_.front();
      blocks_ = &blocks_container_.front();
    }
 }
 


  Int FindBlockIdx(Int aiGIndex){
      return global_to_local_index_[aiGIndex-1];      
  }


  Int Shrink(){
    if(b_own_storage_){
      blocks_container_.shrink_to_fit();
      blocks_ = &blocks_container_.front();

      nzval_container_.resize(nzval_cnt_);
      nzval_ = &nzval_container_.front();
    }
    return StorageSize();
  }



};




template <typename T> inline std::ostream& operator<<( std::ostream& os, const SuperNode<T>& snode){
  os<<"ooooooooooo   Supernode "<<snode.Id()<<" oooooooooooo"<<endl;
  os<<"     size = "<<snode.Size()<<endl;
  os<<"     fc   = "<<snode.FirstCol()<<endl;
  os<<"     lc   = "<<snode.LastCol()<<endl;
  for(Int blkidx =0; blkidx<snode.NZBlockCnt();++blkidx){
     NZBlockDesc & nzblk_desc = snode.GetNZBlockDesc(blkidx);
     T * nzblk_nzval = snode.GetNZval(nzblk_desc.Offset);
     os<<"--- NZBlock "<<nzblk_desc.GIndex<<" ---"<<endl;
     for(Int i = 0; i<snode.NRows(blkidx);++i){
        for(Int j = 0; j<snode.Size();++j){
          os<<nzblk_nzval[i*snode.Size()+j]<<" ";
        }
        os<<endl;
     }
  }
  os<<"oooooooooooooooooooooooooooooooooooooooo"<<endl;
  return os;
}






} // namespace LIBCHOLESKY


#endif // _SUPERNODE_FACT_HPP_
