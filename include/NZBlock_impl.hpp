/// @file NZBlock.hpp
/// @brief 
/// @author Mathias Jacquelin
/// @date 2014-03-03
#ifndef _NZBLOCK_IMPL_HPP_ 
#define _NZBLOCK_IMPL_HPP_

namespace LIBCHOLESKY{

  template <typename T> NZBlock<T>::NZBlock(Int aiNzcnt, Int aiGIndex, Int aiLIndex){
    
    this->storage_.resize(sizeof(Int)*2 + sizeof(T)*aiNzcnt);
    this->pGIndex_ = reinterpret_cast<Int*>(&storage_[0]);
    this->pLIndex_ = reinterpret_cast<Int*>(&storage_[0])+1;
    this->pNzval_= reinterpret_cast<T*>(this->pLIndex_+1);
    this->iNzcnt_ = aiNzcnt;
    *this->pGIndex_ = aiGIndex;
    *this->pLIndex_ = aiLIndex;
  }



  template <typename T> NZBlock<T>::NZBlock(T * apNzval, Int aiNzcnt, Int aiGIndex, Int aiLIndex){
    if(apNzval==NULL){
      abort();
    }
    
    this->storage_.resize(sizeof(Int)*2 + sizeof(T)*aiNzcnt);
    this->pGIndex_ = reinterpret_cast<Int*>(&storage_[0]);
    this->pLIndex_ = reinterpret_cast<Int*>(&storage_[0])+1;
    this->pNzval_= reinterpret_cast<T*>(this->pLIndex_+1);
    this->iNzcnt_ = aiNzcnt;
    *this->pGIndex_ = aiGIndex;
    *this->pLIndex_ = aiLIndex;
    std::copy(apNzval, apNzval+aiNzcnt, this->pNzval_);
  }


}

#endif
