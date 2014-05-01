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

//  std::vector<char> * storage_lcol_; 
//  Int iCurHeight_;
//  Int iMaxHeight_;
//  NZBlock<T> * head_;
//  char * begin_;
//  char * end_;
//  std::vector<int64_t> Lcol_;
//  std::vector<Int > LocToGlobIndices_;
//
//  inline void advance_head(){head_ = reinterpret_cast<NZBlock<T>*>(reinterpret_cast<char *>(head_) + head_->TotalSize()); };

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
//  inline Int MaxHeight(){ return iMaxHeight_;}
//  inline Int NRows(){ return iCurHeight_;}
//  inline Int LocToGlob(Int aiLocIndex){ return LocToGlobIndices_[aiLocIndex];}
//  inline NZBlock<T> & GetNZBlock(Int aiLocIndex){ return *reinterpret_cast<NZBlock<T>* >(begin_ + Lcol_[aiLocIndex]);}
//  inline char * End(){ return reinterpret_cast<char*>(head_);}


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

//
//     size_t line_size = NZBLOCK_ROW_SIZE<T>(iSize_);
//     size_t storage_size = ai_num_rows*line_size;
//     iMaxHeight_ = ai_num_rows;
//     iCurHeight_ = 0;
//     storage_lcol_ = new std::vector<char>(storage_size);
//     begin_ = &storage_lcol_->front();
//     end_ = &storage_lcol_->back()+1;
//     head_ = reinterpret_cast<NZBlock<T>*>(begin_);
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










///
///
///
///
///  SuperNode(Int aiId, Int aiFc, Int aiLc, 
///              std::vector<char> * a_nzblocks,  bool ab_own_storage = false) :
///                                  iId_(aiId),iFirstCol_(aiFc), iLastCol_(aiLc),
///                                                 b_own_storage_(ab_own_storage){
///
///
///
///
///    storage_lcol_ = a_nzblocks;
///    b_own_storage_ = ab_own_storage;
///
///    begin_ = &storage_lcol_->front();
///    end_ = &storage_lcol_->back()+1;
///    head_ = reinterpret_cast<NZBlock<T>*>(begin_);
///
///
///    //Create a nzblock at the beginning of storage because the first block is incomplete
///    char * storage_ptr = reinterpret_cast<char*>(head_+1);
///    NZBlock<T> * p_first_block = new (head_) NZBlock<T>(storage_ptr);
///
///
///
///    iSize_ = iLastCol_ - iFirstCol_+1;
///    size_t line_size = NZBLOCK_ROW_SIZE<T>(iSize_);
///    Int max_nzcnt = (end_-begin_) / line_size;
///    iMaxHeight_ = max_nzcnt;
///
///    //compute the Lcol array and head_ pointer
///    Lcol_.reserve(max_nzcnt);
///    LocToGlobIndices_.reserve(max_nzcnt+1);
///    iCurHeight_ = 0;
///    while(head_<end_){
///      NZBlock<T> * cur_nzblk = head_;
///
///      //Update the internal pointers
///      char * storage_ptr = reinterpret_cast<char*>(head_+1);
///      cur_nzblk->header_ = reinterpret_cast<NZBlockHeader<T> * >(storage_ptr);
///      cur_nzblk->storage_ = reinterpret_cast<char * >(reinterpret_cast<char*>(cur_nzblk->header_));
///      cur_nzblk->pNzval_ = reinterpret_cast<double * >(cur_nzblk->storage_ + NZBLOCK_HEADER_SIZE<T>() );
///
///     iCurHeight_ += cur_nzblk->NRows();
///    //  logfileptr->OFS()<<*cur_nzblk<<std::endl;
///
///      //Update the offset array
///      Lcol_.push_back(reinterpret_cast<char*>(head_) - begin_);
///      LocToGlobIndices_.push_back(head_->GIndex());
///
///      //advance head_ pointer
///      advance_head();
///    }
///    Int last_nzblk_index = NZBlockCnt()-1;
///    NZBlock<T> & last_nzblk = GetNZBlock(last_nzblk_index);
///    LocToGlobIndices_.push_back(last_nzblk.GIndex()+last_nzblk.NRows());
///
///    Lcol_.shrink_to_fit();
///    LocToGlobIndices_.shrink_to_fit();
///
/////    assert(iCurHeight_<=iMaxHeight_);
///  }
/*
//  SuperNode(SuperNode const& C){
//    iId_ = C.iId_;
//    iSize_ = C.iSize_;
//    iFirstCol_=C.iFirstCol_;
//    iLastCol_=C.iLastCol_;
//    
//    b_own_storage_ = C.b_own_storage;
//    if(b_own_storage_){
//      storage_lcol_ = new std::vector<char>(C.storage_lcol_->size());
//    }
//    else{
//      storage_lcol_ = C.storage_lcol_;
//    }
//
//    Lcol_ = C.Lcol_;
//    LocToGlobIndices_ = C.LocToGlobIndices_;
//
//
//    begin_ = &storage_lcol_->front();
//    end_ = &storage_lcol_->back()+1;
//    //recompute head pointer
//    head_ = reinterpret_cast<NZBlock<T>*>(reinterpret_cast<char*>(C.head_) - &C.begin_ + begin_);
//  }
*/

  ~SuperNode(){
//    if(b_own_storage_){
//      delete storage_lcol_;
//    }
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

    


//
//
//
//    iCurHeight_+=aiNRows;
//    assert(iCurHeight_<=iMaxHeight_);
//
//
//    Lcol_.push_back(reinterpret_cast<char*>(head_) - begin_);
//    
//    //placement new
//    void * storage_ptr = static_cast<void*>(head_+1);
//    void * end_ptr = static_cast<void*>(static_cast<char *>(storage_ptr) + NZBLOCK_HEADER_SIZE<T>() + aiNRows*aiNCols*sizeof(T));
//    assert(end_ptr<= end_);
//
////    logfileptr->OFS()<<"Creating an object at "<<head_<<" and storage address will be "<<storage_ptr<<" with "<<aiNRows*aiNCols<<" nnz"<<std::endl;
//    new (head_) NZBlock<T>(aiNRows, aiNCols, aiGIndex, storage_ptr);
////    logfileptr->OFS()<<"End of object at "<<end_ptr<<std::endl;
//
//    //advance head
//    advance_head();
//
//
//
//
//
//
//
//
//
//    if(LocToGlobIndices_.size()==0){
//      LocToGlobIndices_.push_back(aiGIndex);
//    }
//    else{
//      LocToGlobIndices_.back() = aiGIndex;
//    }
//    //add one extra past last row
//    LocToGlobIndices_.push_back(aiGIndex + aiNRows);

 }
 


  Int FindBlockIdx(Int aiGIndex){
//    Int iNZBlockIdx = 0;
//    for(iNZBlockIdx =0;iNZBlockIdx<LocToGlobIndices_.size();++iNZBlockIdx){ 
//      if(LocToGlob(iNZBlockIdx) == aiGIndex){
//        break;
//      }
//      else if(LocToGlob(iNZBlockIdx) > aiGIndex){
//        --iNZBlockIdx;
//        break;
//      }
//    }
//
//    return iNZBlockIdx;
      return global_to_local_index_[aiGIndex-1];      

  }


  Int Shrink(){
    if(b_own_storage_){
      blocks_container_.shrink_to_fit();
      blocks_ = &blocks_container_.front();

      nzval_container_.resize(nzval_cnt_);
      nzval_ = &nzval_container_.front();
    }
    //   //compute used size
//    size_t used = End() - begin_;
//    storage_lcol_->resize(used);
//    
//    begin_ = &storage_lcol_->front();
//    end_ = &storage_lcol_->back()+1;
//
//    head_ = &GetNZBlock(NZBlockCnt()-1);
//    advance_head();
//
//
//     iMaxHeight_ = iCurHeight_;
//
    return StorageSize();
  }


  //extend the storage size to numBytes
//  Int Extend(Int numBytes){
//    size_t used = End() - begin_;
//    assert(numBytes>=used);
//    storage_lcol_->resize(numBytes);
//
//    begin_ = &storage_lcol_->front();
//    end_ = &storage_lcol_->back()+1;
//    
//  
//    head_ = &GetNZBlock(NZBlockCnt()-1);
//    advance_head();
//
//
//     size_t line_size = NZBLOCK_ROW_SIZE<T>(iSize_);
//     iMaxHeight_ = max(iCurHeight_,(Int)floor(numBytes/line_size));
//
//    return StorageSize();
//
//  }



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
