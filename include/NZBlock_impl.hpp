/// @file NZBlock.hpp
/// @brief 
/// @author Mathias Jacquelin
/// @date 2014-03-03
#ifndef _NZBLOCK_IMPL_HPP_ 
#define _NZBLOCK_IMPL_HPP_

namespace LIBCHOLESKY{

  template <typename T> NZBlock<T>::NZBlock(void * apStorage){
    //I dont own the header because it is embedded in the common storage space

    this->header_ = apStorage;
    this->bHas_header_ = true;
    this->bOwn_header_ = false;

    this->storage_ = apStorage;

    this->pNzval_= reinterpret_cast<T*>(apStorage + NZBLOCK_HEADER_SIZE<T>());
  }




  template <typename T> NZBlock<T>::NZBlock(Int aiNRows, Int aiNCols, Int aiGIndex, void * apStorage, T * apNzval){
    Int iNzcnt = aiNRows*aiNCols;

    //create header
    this->header_ = new (apStorage) NZBlockHeader<T>(aiNRows,aiNCols,aiGIndex);
    this->bHas_header_ = true;
    //I dont own the header because it is embedded in the common storage space
    this->bOwn_header_ = false;

    this->storage_ = apStorage;

    this->pNzval_= reinterpret_cast<T*>(apStorage + sizeof(NZBlockHeader<T>));

    if(apNzval!=NULL){
      std::copy(apNzval, apNzval+iNzcnt, this->pNzval_);
    }
  }

  template <typename T> inline T & NZBlock<T>::Nzval(Int i, Int j, Int LDA = -1){
    if( i < 0 || i >= NRows() ||
        j < 0 || j >= NCols() ) {

      std::stringstream ss;
      ss<<"ERROR2: i = "<<i<<" >= "<<NRows()<<" ?  or j = "<<j<<" >= "<<NCols()<<" ?"<<std::endl; 
      logfileptr->OFS()<<ss.str();

#ifdef USE_ABORT
      printf("%s",ss.str().c_str());
      abort();
#endif
      //throw std::logic_error( ss.str().c_str() );
    }


    if(LDA == -1){
      LDA = this->LDA();
    }

#ifdef ROW_MAJOR
    //Row major
    return pNzval_[j+i*LDA];  
#else
    //Col major
    return pNzval_[i+j*LDA];
#endif
  }

//  template <typename T> inline const T& NZBlock<T>::operator()(Int i, Int j) const  { 
//    return const_cast<T&>(this->Nzval(i,j));
//  }

//  template <typename T> inline T& NZBlock<T>::operator()(Int i, Int j) { 
//    return this->Nzval(i,j);
//  }


  template <typename T> inline T & NZBlock<T>::Nzval(Int idx){
    if( idx < 0 || idx >= Nzcnt() ) {

      std::stringstream ss;
      ss<<"ERROR: idx = "<<idx<<" >= "<< Nzcnt()<<" ?"<<std::endl; 
      logfileptr->OFS()<<ss.str();

#ifdef USE_ABORT
      printf("%s",ss.str().c_str());
      abort();
#endif
      //throw std::logic_error( ss.str().c_str() );
    }
    return pNzval_[idx];
  }

//  template <typename T> inline const T& NZBlock<T>::operator[](Int idx) const  { 
//    return const_cast<T&>(this->Nzval(idx));
//  }


//  template <typename T> inline T& NZBlock<T>::operator[](Int idx) { 
//    return this->Nzval(idx);
//  }


template <typename T> inline std::ostream& operator<<( std::ostream& os, const NZBlockHeader<T>& aBlock)
{
  os<<"----Begin of NZBlockHeader----"<<std::endl;
  os<<"Global index is "<<aBlock.GIndex()<<std::endl;
  os<<"NRows is "<<aBlock.NRows()<<std::endl;
  os<<"NCols is "<<aBlock.NCols()<<std::endl;
  os<<"----End of NZBlockHeader----"<<std::endl;
  return os;
}

template <typename T> inline std::ostream& operator<<( std::ostream& os, const NZBlock<T>& aBlock)
{
  os<<"----Begin of NZBlock----"<<std::endl;
  os<<aBlock.GetHeader()<<std::endl;
  os<<"Data pointer is "<<aBlock.Data()<<std::endl;
  os<<"Nzval pointer is "<<aBlock.Nzval()<<std::endl;
  os<<"Nzcnt is "<<aBlock.Nzcnt()<<std::endl;
  os<<"ByteSize is "<<aBlock.ByteSize()<<std::endl;
  os<<"TotalSize is "<<aBlock.TotalSize()<<std::endl;
  os<<"Values are: "<<std::endl;
  for(Int i = 0; i< /*min(2,*/aBlock.NRows()/*)*/; ++i){
    for(Int j = 0; j< /*min(2,*/aBlock.NCols()/*)*/; ++j){
      os<<aBlock.Nzval(i,j)<<" ";
    }
    os<<std::endl;
  }
  os<<"----End of NZBlock----"<<std::endl;
  return os;
}



}

#endif
