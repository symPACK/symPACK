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


} // namespace LIBCHOLESKY


#endif // _SUPERNODE_FACT_HPP_
