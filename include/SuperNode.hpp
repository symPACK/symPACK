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

#define T double
//template<typename T>
class SuperNode{
  protected:
  Int iId_;
  Int iSize_;
  Int iFirstCol_;
  Int iLastCol_;
  
  std::vector<NZBlock<T> * > Lcol_;
  std::vector<Int > LocToGlobIndices_;
  public:
  SuperNode(){}; 
  SuperNode(Int aiId, Int aiFc, Int aiLc):iId_(aiId),iFirstCol_(aiFc),iLastCol_(aiLc){ iSize_ = iLastCol_ - iFirstCol_+1; }; 

  inline void AddNZBlock(NZBlock<T> * apNewBlock){
    if(Lcol_.size()==0){
      LocToGlobIndices_.push_back(apNewBlock->GIndex());
    }
    else{
      LocToGlobIndices_.back() = apNewBlock->GIndex();
    }
    //add one axtra past last row
    LocToGlobIndices_.push_back(apNewBlock->GIndex()+apNewBlock->NRows());



    Lcol_.push_back(apNewBlock);
  }

  inline void AddNZBlock(Int aiNRows, Int aiNCols, Int aiGIndex){ 
    NZBlock<T> * apNewBlock = new NZBlock<double>(aiNRows , aiNCols,aiGIndex,Lcol_.size()+1);
    AddNZBlock(apNewBlock);
  }
 
  inline Int Size(){ return iSize_;}
  inline Int Id(){ return iId_;}
  inline Int FirstCol(){ return iFirstCol_;}
  inline Int LastCol(){ return iLastCol_;}
  inline Int LocToGlob(Int aiLocIndex){ return LocToGlobIndices_[aiLocIndex];}
  inline NZBlock<T> & GetNZBlock(Int aiLocIndex){ return *Lcol_[aiLocIndex];}
  inline Int NZBlockCnt(){ return Lcol_.size();}
 
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


};

#undef T 

template<typename T>
class SuperNode2{
  protected:
/*
  std::vector<char> storage_lcol2_; 
  Int * piId_;
  Int * piSize_;
  Int * piFirstCol_;
  Int * piLastCol_;
  Int * piArrSize_;
  Int * piNzblkCnt_;
  int64_t * piOffsetNzblkStart_;
  int64_t * piOffsetHead_;

  //byte offset
  int64_t * piOffsets_;
  Int * piLocToGlobIndices_;
*/


  bool b_own_storage_;
  std::vector<char> * storage_lcol_; 
  Int iId_;
  Int iSize_;
  Int iFirstCol_;
  Int iLastCol_;
  NZBlock2<T> * head_;
  char * begin_;
  char * end_;
  std::vector<int64_t> Lcol_;
  std::vector<Int > LocToGlobIndices_;

  public:


  inline Int NZBlockCnt(){ return Lcol_.size();}
  inline Int Size(){ return iSize_;}
  inline Int Id(){ return iId_;}
  inline Int FirstCol(){ return iFirstCol_;}
  inline Int LastCol(){ return iLastCol_;}
  inline Int LocToGlob(Int aiLocIndex){ return LocToGlobIndices_[aiLocIndex];}
  inline NZBlock2<T> & GetNZBlock(Int aiLocIndex){ return *reinterpret_cast<NZBlock2<T>* >(begin_ + Lcol_[aiLocIndex]);}
  inline size_t BytesToEnd(Int aiLocIndex){ return reinterpret_cast<char*>(head_) - Lcol_[aiLocIndex] - begin_;}
 





//  SuperNode2(){}; 
  SuperNode2(Int aiId, Int aiFc, Int aiLc, Int ai_num_rows) :iId_(aiId),iFirstCol_(aiFc),iLastCol_(aiLc) {
     b_own_storage_ = true;
     //this is an upper bound
     assert(ai_num_rows>0);
     iSize_ = iLastCol_ - iFirstCol_+1;
//     size_t storage_size = ai_num_rows*(2*sizeof(Int)+iSize_*sizeof(T) + sizeof(NZBlock2<T>));
     size_t storage_size = ai_num_rows*(iSize_*sizeof(T) + sizeof(NZBlock2<T>));
     storage_lcol_ = new std::vector<char>(storage_size);
     begin_ = &storage_lcol_->front();
     end_ = &storage_lcol_->back()+1;
     head_ = reinterpret_cast<NZBlock2<T>*>(begin_);
     int * tmp1 = (int*)begin_;
     int * tmp2 = (int*)end_;
     logfileptr->OFS()<<"----------------------------------------"<<std::endl<<"Storage will range from address "<<tmp1<<" to "<<tmp2<<std::endl;
     logfileptr->OFS()<<" It corresponds to "<<ai_num_rows<<" NZBlocks of 1 line of size "<<iSize_*sizeof(T) + /*2*sizeof(Int) +*/ sizeof(NZBlock2<T>)<<" bytes ("<<storage_size<<")"<<std::endl;

/*
//if(0)
//{
//     //self contained version
//     assert(ai_num_rows>0);
//     size_t nzblk_size = (2*sizeof(Int)+(aiLc - aiFc +1)*sizeof(T) + sizeof(NZBlock2<T>));
//     storage_lcol2_.resize(6*sizeof(Int) + sizeof(char*) + sizeof(NZBlock2<T>*)
//                             + sizeof(int64_t)*ai_num_rows  + sizeof(Int)*(ai_num_rows+1)
//                              + ai_num_rows*nzblk_size);
//
//     piId_ = reinterpret_cast<Int *>(&storage_lcol2_[0]);
//     piSize_ = piId_ + 1;
//     piFirstCol_ = piId_ + 2;
//     piLastCol_ = piId_ + 3;
//     piArrSize_ = piId_ + 4;
//     piNzblkCnt_ = piId_ + 5;
//     piOffsetNzblkStart_ = reinterpret_cast<int64_t *>(piId_ + 6);
//     piOffsetHead_ = piOffsetNzblkStart_ + 1;
//     piOffsets_ = piOffsetHead_ +1;
//     piLocToGlobIndices_ = reinterpret_cast<Int*>(piOffsets_ + ai_num_rows);
//  
//     logfileptr->OFS() <<"piId_ = "<<piId_<<std::endl
//                       <<"piSize_ = "<<piSize_<<std::endl
//                       <<"piFirstCol_ = "<<piFirstCol_<<std::endl
//                       <<"piLastCol_ = "<<piLastCol_<<std::endl
//                       <<"piArrSize_ = "<<piArrSize_<<std::endl
//                       <<"piNzblkCnt_ = "<<piNzblkCnt_<<std::endl
//                       <<"piOffsetNzblkStart_ = "<<piOffsetNzblkStart_<<std::endl
//                       <<"piOffsetHead_ = "<<piOffsetHead_<<std::endl
//                       <<"piOffsets_ = "<<piOffsets_<<std::endl
//                       <<"piLocToGlobIndices_ = "<<piLocToGlobIndices_<<std::endl
//                       <<"storage ranges from "<< &storage_lcol2_[0] << " to " << &*storage_lcol2_.end() << std::endl;
//
//
//
//
//     *piId_ = aiId;
//     *piFirstCol_ = aiFc;
//     *piLastCol_ = aiLc;
//     *piSize_ = *piLastCol_ - *piFirstCol_+1;
//     *piArrSize_ = ai_num_rows;
//     *piNzblkCnt_ = 0;
//     *piOffsetNzblkStart_ = reinterpret_cast<char *>(piLocToGlobIndices_+ai_num_rows+1)-&storage_lcol2_[0];
//     *piOffsetHead_ = *piOffsetNzblkStart_;
//}
*/




  }; 

  SuperNode2(Int aiId, Int aiFc, Int aiLc, 
              std::vector<char> * a_nzblocks, bool ab_own_storage = false) :
                                  iId_(aiId),iFirstCol_(aiFc), iLastCol_(aiLc),
                                                 b_own_storage_(ab_own_storage){

    //FIXME: ASSUME THAT a_nzblocks is not starting with a full NZBlock but only some given row
    // The Header of the incomplete nzblk is provided as an argument

    storage_lcol_ = a_nzblocks;
    begin_ = &storage_lcol_->front();
    end_ = &storage_lcol_->back()+1;

    iSize_ = iLastCol_ - iFirstCol_+1;
    size_t line_size = (/*4*sizeof(Int)+*/iSize_*sizeof(T) + sizeof(NZBlock2<T>));
    Int max_nzcnt = (end_-begin_) / line_size;

    //compute the Lcol array and head_ pointer
    Lcol_.reserve(max_nzcnt);
    LocToGlobIndices_.reserve(max_nzcnt+1);
    head_ = reinterpret_cast<NZBlock2<T>*>(begin_);
    while(head_<end_){
      NZBlock2<T> & cur_nzblk = *head_;
     
      Lcol_.push_back(reinterpret_cast<char*>(head_) - begin_);
      LocToGlobIndices_.push_back(head_->GIndex());
      //advance head_ pointer
      head_ = reinterpret_cast<NZBlock2<T>*>(reinterpret_cast<char*>(head_)
                                                         + head_->ByteSize())+1;
    }
    Int last_nzblk_index = NZBlockCnt()-1;
    NZBlock2<T> & last_nzblk = GetNZBlock(last_nzblk_index);
    LocToGlobIndices_.push_back(last_nzblk.GIndex()+last_nzblk.NRows());

    Lcol_.shrink_to_fit();
    LocToGlobIndices_.shrink_to_fit();

  }

  SuperNode2(SuperNode2 const& C){
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
    head_ = reinterpret_cast<NZBlock2<T>*>(reinterpret_cast<char*>(C.head_) - &C.begin_ + begin_);
  }

  ~SuperNode2(){
    if(b_own_storage_){
      delete storage_lcol_;
    }
  }
    

//  inline Int NZBlockCnt(){ return *piNzblkCnt_;}
//  inline Int Size(){ return *piSize_;}
//  inline Int Id(){ return *piId_;}
//  inline Int FirstCol(){ return *piFirstCol_;}
//  inline Int LastCol(){ return *piLastCol_;}
//  inline Int LocToGlob(Int aiLocIndex){ return piLocToGlobIndices_[aiLocIndex];}
//  inline NZBlock2<T> & GetNZBlock(Int aiLocIndex){
//       NZBlock2<T> * ptr = reinterpret_cast<NZBlock2<T>* >(&storage_lcol2_[0] + *piOffsetNzblkStart_ + piOffsets_[aiLocIndex]);
//              logfileptr->OFS()<<"ptr is "<<ptr<<" offsetStart is "<<*piOffsetNzblkStart_<<" piOffset is "<<piOffsets_[aiLocIndex]<<" " <<NZBlockCnt()<<std::endl; 
//        return *reinterpret_cast<NZBlock2<T>* >(&storage_lcol2_[0] + *piOffsetNzblkStart_ + piOffsets_[aiLocIndex]);}


 inline void AddNZBlock(Int aiNRows, Int aiNCols, Int aiGIndex){
    Int iLIndex = Lcol_.size()+1;
  
    Lcol_.push_back(reinterpret_cast<char*>(head_) - begin_);
    
    //placement new
    void * storage_ptr = static_cast<void*>(head_+1);
    void * end_ptr = static_cast<void*>(static_cast<char *>(storage_ptr) + /*2*sizeof(Int)+*/ aiNRows*aiNCols*sizeof(T));
    assert(end_ptr<= end_);

    logfileptr->OFS()<<"Creating an object at "<<head_<<" and storage address will be "<<storage_ptr<<" with "<<aiNRows*aiNCols<<" nnz"<<std::endl;
    new (head_) NZBlock2<T>(aiNRows, aiNCols, aiGIndex, iLIndex, storage_ptr);
    logfileptr->OFS()<<"End of object at "<<end_ptr<<std::endl;

    //advance head
    head_ = reinterpret_cast<NZBlock2<T>*>(reinterpret_cast<char *>(head_) + head_->ByteSize())+1; 

    if(LocToGlobIndices_.size()==0){
      LocToGlobIndices_.push_back(aiGIndex);
    }
    else{
      LocToGlobIndices_.back() = aiGIndex;
    }
    //add one extra past last row
    LocToGlobIndices_.push_back(aiGIndex + aiNRows);
/*{
    Int iLIndex = *piNzblkCnt_+1;
    //self contained version
    piOffsets_[iLIndex-1] = *piOffsetHead_ - *piOffsetNzblkStart_;
    
    //placement new
    NZBlock2<T> * head = reinterpret_cast<NZBlock2<T>*>(&storage_lcol2_[0]+ *piOffsetHead_);
    void * storage_ptr = static_cast<void*>(head+1);
    void * end_ptr = static_cast<void*>(static_cast<char *>(storage_ptr) + 2*sizeof(Int)+ aiNRows*aiNCols*sizeof(T));
    assert(end_ptr<= &*storage_lcol2_.end());

    logfileptr->OFS()<<"Creating an object at "<<head<<" and storage address will be "<<storage_ptr<<std::endl;
    new (head) NZBlock2<T>(aiNRows, aiNCols, aiGIndex, iLIndex, storage_ptr);
    logfileptr->OFS()<<"End of object at "<<end_ptr<<std::endl;

    //advance head
    *piOffsetHead_ = *piOffsetHead_ + head->ByteSize() + sizeof(NZBlock2<T>); 

    piLocToGlobIndices_[iLIndex-1] = aiGIndex;
    //add one extra past last row
    piLocToGlobIndices_[iLIndex] = aiGIndex + aiNRows;
    *piNzblkCnt_ = iLIndex;
} 
*/
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
/*{
    for(iNZBlockIdx =0;iNZBlockIdx<*piNzblkCnt_;++iNZBlockIdx){ 
      if(LocToGlob(iNZBlockIdx) == aiGIndex){
        break;
      }
      else if(LocToGlob(iNZBlockIdx) > aiGIndex){
        --iNZBlockIdx;
        break;
      }
    }
}
*/
    return iNZBlockIdx;

  }

  Int Shrink(){
 //   //compute used size
 //   char * head = NULL;
 //   if(Lcol_.size()==0){
 //     head = reinterpret_cast<char * > (&storage_lcol_[0]);
 //   }
 //   else{
////      head = reinterpret_cast< char *>(reinterpret_cast<char*>(Lcol_.back()) + Lcol_.back()->ByteSize() + sizeof(NZBlock2<T>));
 //     head = reinterpret_cast< char *>(&storage_lcol_[0] + Lcol_.back()  + GetNZBlock(Lcol_.size()-1).ByteSize() + sizeof(NZBlock2<T>));
 //   }

 //   Int numBytes = head - &storage_lcol_[0];
 //   logfileptr->OFS()<<"We are using "<<numBytes<<" out of "<<storage_lcol_.size()<<" in the collection"<<std::endl;
 //   storage_lcol_.resize(numBytes);
  }

};










} // namespace LIBCHOLESKY


#endif // _SUPERNODE_FACT_HPP_
