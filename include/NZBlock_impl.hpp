/// @file NZBlock.hpp
/// @brief 
/// @author Mathias Jacquelin
/// @date 2014-03-03
#ifndef _NZBLOCK_IMPL_HPP_ 
#define _NZBLOCK_IMPL_HPP_

namespace LIBCHOLESKY{

  //  template <typename T> NZBlock<T>::NZBlock(Int aiNzcnt, Int aiGIndex, Int aiLIndex){
  //    
  //    this->storage_.resize(sizeof(Int)*2 + sizeof(T)*aiNzcnt);
  //    this->pGIndex_ = reinterpret_cast<Int*>(&storage_[0]);
  //    this->pLIndex_ = reinterpret_cast<Int*>(&storage_[0])+1;
  //    this->pNzval_= reinterpret_cast<T*>(this->pLIndex_+1);
  //    this->iNzcnt_ = aiNzcnt;
  //    *this->pGIndex_ = aiGIndex;
  //    *this->pLIndex_ = aiLIndex;
  //  }

  template <typename T> NZBlock<T>::NZBlock(Int aiNRows, Int aiNCols, Int aiGIndex, Int aiLIndex,T * apNzval){
    Int iNzcnt = aiNRows*aiNCols;
    this->storage_.resize(sizeof(Int)*2 + sizeof(T)*iNzcnt);
    this->pGIndex_ = reinterpret_cast<Int*>(&storage_[0]);
    this->pLIndex_ = reinterpret_cast<Int*>(&storage_[0])+1;
    this->pNzval_= reinterpret_cast<T*>(this->pLIndex_+1);
    //this->iNzcnt_ = aiNzcnt;
    this->iNRows_=aiNRows;
    this->iNCols_=aiNCols;
    *this->pGIndex_ = aiGIndex;
    *this->pLIndex_ = aiLIndex;

    if(apNzval!=NULL){
      std::copy(apNzval, apNzval+iNzcnt, this->pNzval_);
    }
  }

  //  template <typename T> NZBlock<T>::NZBlock(T * apNzval, Int aiNzcnt, Int aiGIndex, Int aiLIndex){
  //    if(apNzval==NULL){
  //      abort();
  //    }
  //    
  //    this->storage_.resize(sizeof(Int)*2 + sizeof(T)*aiNzcnt);
  //    this->pGIndex_ = reinterpret_cast<Int*>(&storage_[0]);
  //    this->pLIndex_ = reinterpret_cast<Int*>(&storage_[0])+1;
  //    this->pNzval_= reinterpret_cast<T*>(this->pLIndex_+1);
  //    this->iNzcnt_ = aiNzcnt;
  //    *this->pGIndex_ = aiGIndex;
  //    *this->pLIndex_ = aiLIndex;
  //    std::copy(apNzval, apNzval+aiNzcnt, this->pNzval_);
  //  }

  template <typename T> NZBlock<T>& NZBlock<T>::Copy(const NZBlock& C) {
    storage_ = C.storage_;

    this->pGIndex_ = reinterpret_cast<Int*>(&storage_[0]);
    this->pLIndex_ = reinterpret_cast<Int*>(&storage_[0])+1;
    this->pNzval_= reinterpret_cast<T*>(this->pLIndex_+1);
    this->iNRows_=C.iNRows_;
    this->iNCols_= C.iNCols_;
    *this->pGIndex_ = C.iGIndex_;
    *this->pLIndex_ = C.iLIndex_;

    return *this;
  }


  template <typename T> NZBlock<T>& NZBlock<T>::operator=(const NZBlock& C) {
    return this->Copy(C);
  }


  template <typename T> void NZBlock<T>::Clear()  {
    storage_.clear();
  }

  template <typename T> inline T & NZBlock<T>::Nzval(Int i, Int j){
    if( i < 0 || i >= iNRows_ ||
        j < 0 || j >= iNCols_ ) {

      std::stringstream ss;
      ss<<"ERROR: i = "<<i<<" >= "<<iNRows_<<" ?  or j = "<<j<<" >= "<<iNCols_<<" ?"<<std::endl; 
      logfileptr->OFS()<<ss.str();

#ifdef USE_ABORT
      printf("%s",ss.str().c_str());
      abort();
#endif
      throw std::logic_error( ss.str().c_str() );
    }

    //Col major
    return pNzval_[i+j*iNRows_];
    //Row major
    //return pNzval_[j+i*iNCols_];  
  }

  template <typename T> inline const T& NZBlock<T>::operator()(Int i, Int j) const  { 
    return const_cast<T&>(this->Nzval(i,j));
  }

  template <typename T> inline T& NZBlock<T>::operator()(Int i, Int j) { 
    return this->Nzval(i,j);
  }


  template <typename T> inline T & NZBlock<T>::Nzval(Int idx){
    if( idx < 0 || idx >= iNCols_*iNRows_ ) {

      std::stringstream ss;
      ss<<"ERROR: idx = "<<idx<<" >= "<<iNCols_*iNRows_<<" ?"<<std::endl; 
      logfileptr->OFS()<<ss.str();

#ifdef USE_ABORT
      printf("%s",ss.str().c_str());
      abort();
#endif
      throw std::logic_error( ss.str().c_str() );
    }
    return pNzval_[idx];
  }
  template <typename T> inline const T& NZBlock<T>::operator[](Int idx) const  { 
    return const_cast<T&>(this->Nzval(idx));
  }


  template <typename T> inline T& NZBlock<T>::operator[](Int idx) { 
    return this->Nzval(idx);
  }





template <typename T> inline std::ostream& operator<<( std::ostream& os, const NZBlock<T>& aBlock)
{
  os<<"----Begin of NZBlock----"<<std::endl;
  os<<"Global index is "<<aBlock.GIndex()<<std::endl;
  os<<"Local index is "<<aBlock.LIndex()<<std::endl;
  os<<"Data pointer is "<<aBlock.Data()<<std::endl;
  os<<"Nzval pointer is "<<aBlock.Nzval()<<std::endl;
  os<<"NRows is "<<aBlock.NRows()<<std::endl;
  os<<"NCols is "<<aBlock.NCols()<<std::endl;
  os<<"Nzcnt is "<<aBlock.Nzcnt()<<std::endl;
  //os<<"Values are: "<<std::endl;
  //for(Int i = 0; i< aBlock.NRows(); ++i){
  //  for(Int j = 0; j< aBlock.NCols(); ++j){
  //    os<<aBlock.Nzval(i,j)<<" ";
  //  }
  //  os<<std::endl;
  //}
  os<<"----End of NZBlock----"<<std::endl;
  return os;
}




  template <typename T> NZBlock2<T>::NZBlock2(Int aiNRows, Int aiNCols, Int aiGIndex, Int aiLIndex,void * apStorage, T * apNzval){
    Int iNzcnt = aiNRows*aiNCols;
    this->storage_ = apStorage;
//    this->pGIndex_ = reinterpret_cast<Int*>(storage_);
//    this->pLIndex_ = reinterpret_cast<Int*>(storage_)+1;
//    this->piNRows_ = reinterpret_cast<Int*>(storage_)+2;
//    this->piNCols_ = reinterpret_cast<Int*>(storage_)+3;
//    this->pNzval_= reinterpret_cast<T*>(this->piNCols_+1);
//    *this->piNRows_ = aiNRows;
//    *this->piNCols_ = aiNCols;
//    *this->pGIndex_ = aiGIndex;
//    *this->pLIndex_ = aiLIndex;

    this->iNRows_ = aiNRows;
    this->iNCols_ = aiNCols;
    this->GIndex_ = aiGIndex;
    this->LIndex_ = aiLIndex;
    this->pNzval_= reinterpret_cast<T*>(storage_);

    if(apNzval!=NULL){
      std::copy(apNzval, apNzval+iNzcnt, this->pNzval_);
    }
  }

  //  template <typename T> NZBlock2<T>::NZBlock2(T * apNzval, Int aiNzcnt, Int aiGIndex, Int aiLIndex){
  //    if(apNzval==NULL){
  //      abort();
  //    }
  //    
  //    this->storage_.resize(sizeof(Int)*2 + sizeof(T)*aiNzcnt);
  //    this->pGIndex_ = reinterpret_cast<Int*>(&storage_[0]);
  //    this->pLIndex_ = reinterpret_cast<Int*>(&storage_[0])+1;
  //    this->pNzval_= reinterpret_cast<T*>(this->pLIndex_+1);
  //    this->iNzcnt_ = aiNzcnt;
  //    *this->pGIndex_ = aiGIndex;
  //    *this->pLIndex_ = aiLIndex;
  //    std::copy(apNzval, apNzval+aiNzcnt, this->pNzval_);
  //  }

//  template <typename T> NZBlock2<T>& NZBlock2<T>::Copy(const NZBlock2& C) {
//    storage_ = C.storage_;
//
//    this->pGIndex_ = reinterpret_cast<Int*>(&storage_[0]);
//    this->pLIndex_ = reinterpret_cast<Int*>(&storage_[0])+1;
//    this->pNzval_= reinterpret_cast<T*>(this->pLIndex_+1);
//    this->iNRows_=C.iNRows_;
//    this->iNCols_= C.iNCols_;
//    *this->pGIndex_ = C.iGIndex_;
//    *this->pLIndex_ = C.iLIndex_;
//
//    return *this;
//  }


//  template <typename T> NZBlock2<T>& NZBlock2<T>::operator=(const NZBlock2& C) {
//    return this->Copy(C);
//  }


//  template <typename T> void NZBlock2<T>::Clear()  {
//    storage_.clear();
//  }

  template <typename T> inline T & NZBlock2<T>::Nzval(Int i, Int j){
    if( i < 0 || i >= NRows() ||
        j < 0 || j >= NCols() ) {

      std::stringstream ss;
      ss<<"ERROR: i = "<<i<<" >= "<<NRows()<<" ?  or j = "<<j<<" >= "<<NCols()<<" ?"<<std::endl; 
      logfileptr->OFS()<<ss.str();

#ifdef USE_ABORT
      printf("%s",ss.str().c_str());
      abort();
#endif
      throw std::logic_error( ss.str().c_str() );
    }

    //Col major
    return pNzval_[i+j*iNRows_];
    //Row major
    //return pNzval_[j+i*iNCols_];  
  }

  template <typename T> inline const T& NZBlock2<T>::operator()(Int i, Int j) const  { 
    return const_cast<T&>(this->Nzval(i,j));
  }

  template <typename T> inline T& NZBlock2<T>::operator()(Int i, Int j) { 
    return this->Nzval(i,j);
  }


  template <typename T> inline T & NZBlock2<T>::Nzval(Int idx){
    if( idx < 0 || idx >= Nzcnt() ) {

      std::stringstream ss;
      ss<<"ERROR: idx = "<<idx<<" >= "<< Nzcnt()<<" ?"<<std::endl; 
      logfileptr->OFS()<<ss.str();

#ifdef USE_ABORT
      printf("%s",ss.str().c_str());
      abort();
#endif
      throw std::logic_error( ss.str().c_str() );
    }
    return pNzval_[idx];
  }
  template <typename T> inline const T& NZBlock2<T>::operator[](Int idx) const  { 
    return const_cast<T&>(this->Nzval(idx));
  }


  template <typename T> inline T& NZBlock2<T>::operator[](Int idx) { 
    return this->Nzval(idx);
  }






}

#endif
