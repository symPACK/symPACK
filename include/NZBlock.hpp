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
      //col major
      inline Int LDA() { return iNRows_; };
      //row major
      //inline Int LDA() { return iNCols_; };

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


  template <typename T>
  class NZBlock2{
    protected:
      //std::vector<char> storage_;
      void * storage_;
      //Int iNzcnt_;
//      Int * piNRows_;
//      Int * piNCols_;
//      Int * pGIndex_;
//      Int * pLIndex_;

      Int iNRows_;
      Int iNCols_;
      Int GIndex_;
      Int LIndex_;

      //pointer to the nzdata in the storage;
      T * pNzval_;
      
    public:
      NZBlock2(Int aiNRows, Int aiNCols, Int aiGIndex, Int aiLIndex, void* apStorage, T * apNzval = NULL);

//      inline Int NRows() { return *piNRows_; };
//      inline Int NCols() { return *piNCols_; };
//      inline Int  GIndex() { return *pGIndex_;};
//      inline Int  LIndex() { return *pLIndex_;};
//      inline Int Nzcnt() { return NRows()*NCols(); };
//      inline size_t ByteSize() { return 4*sizeof(Int)+Nzcnt()*sizeof(T); };


      inline Int NRows() { return iNRows_; };
      inline Int NCols() { return iNCols_; };
      inline Int  GIndex() { return GIndex_;};
      inline Int  LIndex() { return LIndex_;};
      inline Int Nzcnt() { return NRows()*NCols(); };
      inline size_t ByteSize() { return Nzcnt()*sizeof(T); };

      //col major
      inline Int LDA() { return NRows(); };
      //row major
      //inline Int LDA() { return NCols(); };

      inline size_t TotalSize() { return sizeof(*this) + ByteSize();};
      inline void * Data() { return storage_;};
      inline T * Nzval() { return pNzval_;};
      inline T & Nzval(Int idx);
      inline T & Nzval(Int i, Int j);


//      NZBlock2& Copy(const NZBlock2& C);
//      NZBlock2& operator=(const NZBlock2& C);
      inline const T& operator()(Int i, Int j) const;
      inline T& operator()(Int i, Int j);

      inline const T& operator[](Int idx) const;
      inline T& operator[](Int idx);

  };


}

#include "NZBlock_impl.hpp"

#endif
