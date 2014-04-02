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

template< typename T> class SuperNode;

  template <typename T>
  class NZBlockHeader{
    friend class SuperNode<T>;
    protected:

      Int iNRows_;
      Int iNCols_;
      Int iGIndex_;
//      Int iLIndex_;

      
    public:
      NZBlockHeader(Int aiNRows =0, Int aiNCols=0, Int aiGIndex=0/*, Int aiLIndex=0*/):
                                       iNRows_(aiNRows),iNCols_(aiNCols),
                                        iGIndex_(aiGIndex)/*,iLIndex_(aiLIndex)*/{};
      inline Int NRows() { return iNRows_; };
      inline Int NCols() { return iNCols_; };
      inline Int  GIndex() { return iGIndex_;};
//      inline Int  LIndex() { return iLIndex_;};
  };




  template <typename T>
  class NZBlock{
    friend class SuperNode<T>;
    protected:
      void * storage_;
      bool bOwn_header_ = false;
      bool bHas_header_ = false;
      NZBlockHeader<T> * header_;

      //pointer to the nzdata in the storage;
      T * pNzval_;
      

    public:


      NZBlock(void* apStorage);
      NZBlock(Int aiNRows, Int aiNCols, Int aiGIndex, void* apStorage, T * apNzval = NULL);
      ~NZBlock(){
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
      inline NZBlockHeader<T> & GetHeader(){ assert(bHas_header_); return *header_;};

      inline Int NRows() { return header_->NRows(); };
      inline Int NCols() { return header_->NCols(); };
      inline Int  GIndex() { return header_->GIndex();};
//      inline Int  LIndex() { return header_->LIndex();};
      
      inline Int Nzcnt() { return NRows()*NCols(); };
      inline size_t ByteSize() { return Nzcnt()*sizeof(T); };

#ifdef ROW_MAJOR
      //row major
      inline Int LDA() { return NCols(); };
#else
      //col major
      inline Int LDA() { return NRows(); };
#endif

      inline size_t TotalSize() { return sizeof(*this) + (bHas_header_?sizeof(*header_):0) + ByteSize();};
      inline void * Data() { return storage_;};
      inline T * Nzval() { return pNzval_;};
      inline T & Nzval(Int idx);
      inline T & Nzval(Int i, Int j, Int LDA = -1 );


//      inline const T& operator()(Int i, Int j) const;
//      inline T& operator()(Int i, Int j);
//
//      inline const T& operator[](Int idx) const;
//      inline T& operator[](Int idx);

  };



template <typename T> inline std::ostream& operator<<( std::ostream& os, const NZBlockHeader<T>& header);
template <typename T> inline std::ostream& operator<<( std::ostream& os, const NZBlock<T>& block);


template <typename T> inline size_t NZBLOCK_HEADER_SIZE(){ return sizeof(NZBlockHeader<T>); };
template <typename T> inline size_t NZBLOCK_OBJ_SIZE() { return sizeof(NZBlock<T>);};
template <typename T> inline size_t NZBLOCK_ROW_SIZE(Int width){ return (width*sizeof(T) + NZBLOCK_HEADER_SIZE<T>() + NZBLOCK_OBJ_SIZE<T>()) ; };

}

#include "NZBlock_impl.hpp"

#endif
