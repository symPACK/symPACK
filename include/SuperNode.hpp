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

template<typename T>
class SuperNode{
  protected:
  bool b_own_storage_;
  std::vector<char> * storage_lcol_; 
  Int iId_;
  Int iSize_;
  Int iFirstCol_;
  Int iLastCol_;
  Int iMaxHeight_;
  NZBlock<T> * head_;
  char * begin_;
  char * end_;
  std::vector<int64_t> Lcol_;
  std::vector<Int > LocToGlobIndices_;

  inline void advance_head(){head_ = reinterpret_cast<NZBlock<T>*>(reinterpret_cast<char *>(head_) + head_->TotalSize()); };

  public:


  inline Int NZBlockCnt(){ return Lcol_.size();}
  inline Int Size(){ return iSize_;}
  inline Int MaxHeight(){ return iMaxHeight_;}
  inline Int Id(){ return iId_;}
  inline Int FirstCol(){ return iFirstCol_;}
  inline Int LastCol(){ return iLastCol_;}
  inline Int LocToGlob(Int aiLocIndex){ return LocToGlobIndices_[aiLocIndex];}
  inline NZBlock<T> & GetNZBlock(Int aiLocIndex){ return *reinterpret_cast<NZBlock<T>* >(begin_ + Lcol_[aiLocIndex]);}



//  inline size_t BytesToEnd(Int aiLocIndex){ return reinterpret_cast<char*>(head_) - Lcol_[aiLocIndex] - begin_;} 
  inline char * End(){ return reinterpret_cast<char*>(head_);}


  SuperNode(Int aiId, Int aiFc, Int aiLc, Int ai_num_rows) :iId_(aiId),iFirstCol_(aiFc),iLastCol_(aiLc) {
     //this is an upper bound
     assert(ai_num_rows>0);
     iSize_ = iLastCol_ - iFirstCol_+1;
     size_t line_size = NZBLOCK_ROW_SIZE<T>(iSize_);
     size_t storage_size = ai_num_rows*line_size;
     iMaxHeight_ = ai_num_rows;
     storage_lcol_ = new std::vector<char>(storage_size);
     b_own_storage_ = true;
     begin_ = &storage_lcol_->front();
     end_ = &storage_lcol_->back()+1;
     head_ = reinterpret_cast<NZBlock<T>*>(begin_);

     //int * tmp1 = (int*)begin_;
     //int * tmp2 = (int*)end_;
     //logfileptr->OFS()<<"----------------------------------------"<<std::endl<<"Storage will range from address "<<tmp1<<" to "<<tmp2<<std::endl;
  }; 

  SuperNode(Int aiId, Int aiFc, Int aiLc, 
              std::vector<char> * a_nzblocks,  bool ab_own_storage = false) :
                                  iId_(aiId),iFirstCol_(aiFc), iLastCol_(aiLc),
                                                 b_own_storage_(ab_own_storage){




    storage_lcol_ = a_nzblocks;
    b_own_storage_ = ab_own_storage;

    begin_ = &storage_lcol_->front();
    end_ = &storage_lcol_->back()+1;
    head_ = reinterpret_cast<NZBlock<T>*>(begin_);


    //Create a nzblock at the beginning of storage because the first block is incomplete
    char * storage_ptr = reinterpret_cast<char*>(head_+1);
    NZBlock<T> * p_first_block = new (head_) NZBlock<T>(storage_ptr);



    iSize_ = iLastCol_ - iFirstCol_+1;
    size_t line_size = NZBLOCK_ROW_SIZE<T>(iSize_);
    Int max_nzcnt = (end_-begin_) / line_size;
    iMaxHeight_ = max_nzcnt;

    //compute the Lcol array and head_ pointer
    Lcol_.reserve(max_nzcnt);
    LocToGlobIndices_.reserve(max_nzcnt+1);
    while(head_<end_){
      NZBlock<T> * cur_nzblk = head_;

      //Update the internal pointers
      char * storage_ptr = reinterpret_cast<char*>(head_+1);
      cur_nzblk->header_ = reinterpret_cast<NZBlockHeader<T> * >(storage_ptr);
      cur_nzblk->storage_ = reinterpret_cast<char * >(reinterpret_cast<char*>(cur_nzblk->header_));
      cur_nzblk->pNzval_ = reinterpret_cast<double * >(cur_nzblk->storage_ + NZBLOCK_HEADER_SIZE<T>() );

    //  logfileptr->OFS()<<*cur_nzblk<<std::endl;

      //Update the offset array
      Lcol_.push_back(reinterpret_cast<char*>(head_) - begin_);
      LocToGlobIndices_.push_back(head_->GIndex());

      //advance head_ pointer
      advance_head();
    }
    Int last_nzblk_index = NZBlockCnt()-1;
    NZBlock<T> & last_nzblk = GetNZBlock(last_nzblk_index);
    LocToGlobIndices_.push_back(last_nzblk.GIndex()+last_nzblk.NRows());

    Lcol_.shrink_to_fit();
    LocToGlobIndices_.shrink_to_fit();
  }

  SuperNode(SuperNode const& C){
    iId_ = C.iId_;
    iSize_ = C.iSize_;
    iFirstCol_=C.iFirstCol_;
    iLastCol_=C.iLastCol_;
    
    b_own_storage_ = C.b_own_storage;
    if(b_own_storage_){
      storage_lcol_ = new std::vector<char>(C.storage_lcol_->size());
    }
    else{
      storage_lcol_ = C.storage_lcol_;
    }

    Lcol_ = C.Lcol_;
    LocToGlobIndices_ = C.LocToGlobIndices_;


    begin_ = &storage_lcol_->front();
    end_ = &storage_lcol_->back()+1;
    //recompute head pointer
    head_ = reinterpret_cast<NZBlock<T>*>(reinterpret_cast<char*>(C.head_) - &C.begin_ + begin_);
  }

  ~SuperNode(){
    if(b_own_storage_){
      delete storage_lcol_;
    }
  }
    

 inline void AddNZBlock(Int aiNRows, Int aiNCols, Int aiGIndex){
    Lcol_.push_back(reinterpret_cast<char*>(head_) - begin_);
    
    //placement new
    void * storage_ptr = static_cast<void*>(head_+1);
    void * end_ptr = static_cast<void*>(static_cast<char *>(storage_ptr) + NZBLOCK_HEADER_SIZE<T>() + aiNRows*aiNCols*sizeof(T));
    assert(end_ptr<= end_);

//    logfileptr->OFS()<<"Creating an object at "<<head_<<" and storage address will be "<<storage_ptr<<" with "<<aiNRows*aiNCols<<" nnz"<<std::endl;
    new (head_) NZBlock<T>(aiNRows, aiNCols, aiGIndex, storage_ptr);
//    logfileptr->OFS()<<"End of object at "<<end_ptr<<std::endl;

    //advance head
    advance_head();

    if(LocToGlobIndices_.size()==0){
      LocToGlobIndices_.push_back(aiGIndex);
    }
    else{
      LocToGlobIndices_.back() = aiGIndex;
    }
    //add one extra past last row
    LocToGlobIndices_.push_back(aiGIndex + aiNRows);

 }
 


  Int FindBlockIdx(Int aiGIndex){
    Int iNZBlockIdx = 0;
    for(iNZBlockIdx =0;iNZBlockIdx<LocToGlobIndices_.size();++iNZBlockIdx){ 
      if(LocToGlob(iNZBlockIdx) == aiGIndex){
        break;
      }
      else if(LocToGlob(iNZBlockIdx) > aiGIndex){
        --iNZBlockIdx;
        break;
      }
    }

    return iNZBlockIdx;

  }

  Int Shrink(){
 //   //compute used size
 //   char * head = NULL;
 //   if(Lcol_.size()==0){
 //     head = reinterpret_cast<char * > (&storage_lcol_[0]);
 //   }
 //   else{
////      head = reinterpret_cast< char *>(reinterpret_cast<char*>(Lcol_.back()) + Lcol_.back()->ByteSize() + sizeof(NZBlock<T>));
 //     head = reinterpret_cast< char *>(&storage_lcol_[0] + Lcol_.back()  + GetNZBlock(Lcol_.size()-1).ByteSize() + sizeof(NZBlock<T>));
 //   }

 //   Int numBytes = head - &storage_lcol_[0];
 //   logfileptr->OFS()<<"We are using "<<numBytes<<" out of "<<storage_lcol_.size()<<" in the collection"<<std::endl;
 //   storage_lcol_.resize(numBytes);
  }

};



template <typename T> inline std::ostream& operator<<( std::ostream& os, const SuperNode<T>& snode){
  os<<"ooooooooooo   Supernode "<<snode.Id()<<" oooooooooooo"<<endl;
  os<<"     size = "<<snode.Size()<<endl;
  os<<"     fc   = "<<snode.FirstCol()<<endl;
  os<<"     lc   = "<<snode.LastCol()<<endl;
  for(Int i =0; i<snode.NZBlockCnt();++i){
    os<<snode.GetNZBlock(i)<<endl;
  }
  os<<"oooooooooooooooooooooooooooooooooooooooo"<<endl;
  return os;
}







} // namespace LIBCHOLESKY


#endif // _SUPERNODE_FACT_HPP_
