/// @file NZBlock.hpp
/// @brief 
/// @author Mathias Jacquelin
/// @date 2014-03-03
#ifndef _NZBLOCK_HPP_ 
#define _NZBLOCK_HPP_

#include "Environment.hpp"

#include <vector>


using namespace std;

namespace LIBCHOLESKY{

  template <typename T>
  class NZBlock{
    protected:
      std::vector<char> storage_;
      Int iNzcnt_;
      Int * pGIndex_;
      Int * pLIndex_;
      //pointer to the nzdata in the storage;
      T * pNzval_;
      
    public:
      NZBlock(T * apNzval, Int aiNzcnt, Int aiGIndex, Int aiLIndex);
      inline Int Nzcnt() { return iNzcnt_; };
      inline size_t ByteSize() { return storage_.size(); };
      inline void * Data() { return static_cast<void*>(&storage_[0]);};
      inline T * Nzval() { return pNzval_;};
      inline T & Nzval(Int i) { if(i>=iNzcnt_){ abort();} return pNzval_[i];};
      inline Int  GIndex() { return *pGIndex_;};
      inline Int  LIndex() { return *pLIndex_;};
  };


}

#include "NZBlock_impl.hpp"

#endif
