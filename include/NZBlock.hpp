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

template< typename T> class SuperNode2;

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
  class NZBlockHeader{
    friend class SuperNode2<T>;
    protected:

      Int iNRows_;
      Int iNCols_;
      Int iGIndex_;
      Int iLIndex_;

      
    public:
      NZBlockHeader(Int aiNRows =0, Int aiNCols=0, Int aiGIndex=0, Int aiLIndex=0):
                                       iNRows_(aiNRows),iNCols_(aiNCols),
                                        iGIndex_(aiGIndex),iLIndex_(aiLIndex){};
      inline Int NRows() { return iNRows_; };
      inline Int NCols() { return iNCols_; };
      inline Int  GIndex() { return iGIndex_;};
      inline Int  LIndex() { return iLIndex_;};
  };




  template <typename T>
  class NZBlock2{
    friend class SuperNode2<T>;
    protected:
      //std::vector<char> storage_;
      void * storage_;
      //Int iNRows_;
      //Int iNCols_;
      //Int GIndex_;
      //Int LIndex_;
      bool bOwn_header_ = false;
      bool bHas_header_ = false;
      NZBlockHeader<T> * header_;

      //pointer to the nzdata in the storage;
      T * pNzval_;
      

    public:
      NZBlock2(void* apStorage);
      NZBlock2(Int aiNRows, Int aiNCols, Int aiGIndex, Int aiLIndex, void* apStorage, T * apNzval = NULL);
      ~NZBlock2(){
        if(bHas_header_ && bOwn_header_){
          delete header_;
        }
      };

      void AddHeader(NZBlockHeader<T> * apHeader){
        if(!bHas_header_){
          header_ = apHeader;
          bOwn_header_ = true;
        }
        else{
          abort();
        }
      }; 

//      inline Int NRows() { return iNRows_; };
//      inline Int NCols() { return iNCols_; };
//      inline Int  GIndex() { return GIndex_;};
//      inline Int  LIndex() { return LIndex_;};

      inline Int NRows() { return header_->NRows(); };
      inline Int NCols() { return header_->NCols(); };
      inline Int  GIndex() { return header_->GIndex();};
      inline Int  LIndex() { return header_->LIndex();};
      
      inline Int Nzcnt() { return NRows()*NCols(); };
      inline size_t ByteSize() { return Nzcnt()*sizeof(T); };

      //col major
      inline Int LDA() { return NRows(); };
      //row major
      //inline Int LDA() { return NCols(); };

      inline size_t TotalSize() { return sizeof(*this) + (bHas_header_?sizeof(*header_):0) + ByteSize();};
      inline void * Data() { return storage_;};
      inline T * Nzval() { return pNzval_;};
      inline T & Nzval(Int idx);
      inline T & Nzval(Int i, Int j);


      inline const T& operator()(Int i, Int j) const;
      inline T& operator()(Int i, Int j);

      inline const T& operator[](Int idx) const;
      inline T& operator[](Int idx);

  };



template <typename T> inline std::ostream& operator<<( std::ostream& os, const NZBlockHeader<T>& header);
template <typename T> inline std::ostream& operator<<( std::ostream& os, const NZBlock2<T>& block);



}

#include "NZBlock_impl.hpp"

#endif
