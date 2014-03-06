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
      //Int iNzcnt_;
      Int iNRows_;
      Int iNCols_;
      Int * pGIndex_;
      Int * pLIndex_;
      //pointer to the nzdata in the storage;
      T * pNzval_;
      
    public:
      NZBlock(Int aiNRows, Int aiNCols, Int aiGIndex, Int aiLIndex,T * apNzval = NULL);
//      NZBlock(T * apNzval, Int aiNzcnt, Int aiGIndex, Int aiLIndex);
//      NZBlock(Int aiNzcnt, Int aiGIndex, Int aiLIndex);
//      NZBlock(T * apNzval, Int aiNzcnt, Int aiGIndex, Int aiLIndex);
      inline Int Nzcnt() { return iNCols_*iNRows_; };
      inline Int NRows() { return iNRows_; };
      inline Int NCols() { return iNCols_; };

      inline size_t ByteSize() { return storage_.size(); };
      inline void * Data() { return static_cast<void*>(&storage_[0]);};
      inline T * Nzval() { return pNzval_;};
      inline T & Nzval(Int idx);
      inline T & Nzval(Int i, Int j);
      inline Int  GIndex() { return *pGIndex_;};
      inline Int  LIndex() { return *pLIndex_;};

      void Clear();

      NZBlock& Copy(const NZBlock& C);
      NZBlock& operator=(const NZBlock& C);
      inline const T& operator()(Int i, Int j) const;
      inline T& operator()(Int i, Int j);

      inline const T& operator[](Int idx) const;
      inline T& operator[](Int idx);

  };

template <typename T> inline std::ostream& operator<<( std::ostream& os, const NZBlock<T>& vec);


}

#include "NZBlock_impl.hpp"

#endif
