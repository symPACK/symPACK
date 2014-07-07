/// @file nummat_decl.hpp
/// @brief Numerical matrix.
/// @author Lexing Ying and Lin Lin
/// @date 2010-09-27
#ifndef _NUMMAT_IMPL_HPP_
#define _NUMMAT_IMPL_HPP_

#include "Environment.hpp"
#include "lapack.hpp"
////#include "NZBlock.hpp"
#include <sstream>

#include <stdlib.h> 

namespace LIBCHOLESKY{




  template <typename F> NumMat<F>::NumMat(Int m, Int n): m_(m), n_(n), owndata_(true) {
#ifdef _ASSERT_
    allocated_=false;
#endif
      alloc_data();
  }



  template <typename F> NumMat<F>::NumMat(Int m, Int n, bool owndata, F* data): m_(m), n_(n), owndata_(owndata) {
#ifdef _ASSERT_
    allocated_=false;
#endif
    if(owndata_) {

      alloc_data();

      if(m_>0 && n_>0) { 
        std::copy(data,data+m_*n_,data_);//for(Int i=0; i<m_*n_; i++) data_[i] = data[i]; 
      }
    } 
    else {
      data_ = data;
    }
  }



  template <typename F> NumMat<F>::NumMat(const NumMat& C){
    m_ = C.m_; n_=C.n_; owndata_=C.owndata_;
#ifdef _ASSERT_
    allocated_=false;
#endif
    if(owndata_) {
      alloc_data();
      if(m_>0 && n_>0) { 
        std::copy(C.data_,C.data_+m_*n_,data_);
      }
    } 
    else {
      data_ = C.data_;
    }


  }


  template <typename F> inline void NumMat<F>::alloc_data(){
    //logfileptr->OFS()<<"ALLOCATING DATA"<<std::endl;
#ifdef _ASSERT_
    if(!allocated_) {
      if(owndata_) {
        if(m_>0 && n_>0) { 
          data_ = new F[m_*n_]; 
          if( data_ == NULL ) 
            throw std::runtime_error("Cannot allocate memory."); 
        } 
        else 
          data_=NULL;

        allocated_=true;
      }

    }
    else{

#ifdef USE_ABORT
      printf("%s","Data already allocated.");
      abort();
#endif
          throw std::runtime_error("Data already allocated."); 
    }
#else
      if(owndata_) {
        if(m_>0 && n_>0) { 
          data_ = new F[m_*n_]; 
          if( data_ == NULL ) 
            throw std::runtime_error("Cannot allocate memory."); 
        } 
        else 
          data_=NULL;
      }
#endif
  }

  template <typename F> inline void NumMat<F>::delete_data(){
    if(owndata_) {
#ifdef _ASSERT_
      if(allocated_) {
        if(m_>0 && n_>0) {
          delete[] data_;
          data_ = NULL; 
        }
        allocated_=false;
      }
#else
        if(m_>0 && n_>0) {
          delete[] data_;
          data_ = NULL; 
        }
#endif
      m_=0;
      n_=0;
    }
  }

  template <typename F> NumMat<F>::~NumMat() {
    delete_data();
  }

template <typename F> NumMat<F>& NumMat<F>::Copy(const NumMat& C) {
    delete_data();

    m_ = C.m_; n_=C.n_; owndata_=C.owndata_;
    
    if(owndata_) {
      alloc_data();
      if(m_>0 && n_>0) { 
        std::copy(C.data_,C.data_+m_*n_,data_);
      }
    } 
    else {
      data_ = C.data_;
    }

    return *this;
  }


  template <typename F> NumMat<F>& NumMat<F>::operator=(const NumMat& C) {
    return this->Copy(C);
  }


  template <typename F> NumMat<F>::NumMat(const SuperNode<F>& S){
#ifdef _ASSERT_
    allocated_=false;
#endif
    this->SnodeToDense(S);
  }





template <typename F> NumMat<F>& NumMat<F>::SnodeToDense(const SuperNode<F>& S){
//    delete_data();
//
//    owndata_ = true;
//
//    //compute m_ and n_
//    n_ = S.Size();
//    m_=0;
//    for(Int blkidx=0; blkidx < S.NZBlockCnt(); ++blkidx){
//      m_+=S.GetNZBlock(blkidx).NRows();
//    }
//    
//    alloc_data();
//
//    if(m_>0 && n_>0) {
//      Int head =0;
//      for(Int blkidx=0; blkidx < S.NZBlockCnt(); ++blkidx){
//        NZBlock<F> & cur_block = S.GetNZBlock(blkidx);
//        lapack::Lacpy('N',cur_block.NRows(),cur_block.NCols(),
//                            cur_block.Nzval(),cur_block.NRows(), 
//                                                        &at(head,0), m_);
//        head+=cur_block.NRows();
//      }
//    }
//
    return *this; 
}

  template <typename F> NumMat<F>& NumMat<F>::operator=(const SuperNode<F>& S) {
    return this->SnodeToDense(S);
  }



  template <typename F> void NumMat<F>::Clear()  {
		if( owndata_ == false ){
			throw std::logic_error("Matrix being cleared must own data.");
		}
      delete_data();
  }




  template <typename F> void NumMat<F>::Resize(Int m, Int n)  {
		if( owndata_ == false ){
			throw std::logic_error("Matrix being resized must own data.");
		}
    if(m_!=m || n_!=n) {
      delete_data();
      m_ = m; n_ = n;
      if(m*n>0){
        alloc_data();
      }
    }
  }


  template<typename F> void NumMat<F>::error_message(std::stringstream & ss,Int i, Int j){
      ss<<"Index ("<<i<<","<<j<<") is out of bound. Dimensions are ("<<m_<<","<<n_<<")."<<std::endl;
  }
  template<typename F> void NumMat<F>::error_message(std::stringstream & ss, Int j){
      ss<<"Index (*,"<<j<<") is out of bound. Dimensions are ("<<m_<<","<<n_<<")."<<std::endl;
  }


template <typename F> inline const F& NumMat<F>::at(Int i, Int j) const  { 
		if( i < 0 || i >= m_ ||
				j < 0 || j >= n_ ) {
      
      std::stringstream ss;
      //error_message(ss,i,j);
      logfileptr->OFS()<<ss.str();

#ifdef USE_ABORT
      printf("%s",ss.str().c_str());
      abort();
#endif
			throw std::logic_error( ss.str().c_str() );
		}
    return data_[i+j*m_];
  }


template <typename F> inline F& NumMat<F>::at(Int i, Int j)  { 
		if( i < 0 || i >= m_ ||
				j < 0 || j >= n_ ) {

      std::stringstream ss;
      //error_message(ss,i,j);

      logfileptr->OFS()<<ss.str();

#ifdef USE_ABORT
      printf("%s",ss.str().c_str());
      abort();
#endif
			throw std::logic_error( ss.str().c_str() );
		}
    return data_[i+j*m_];
  }
  




template <typename F> inline const F& NumMat<F>::operator()(Int i, Int j) const  { 
		if( i < 0 || i >= m_ ||
				j < 0 || j >= n_ ) {
      
      std::stringstream ss;
      //error_message(ss,i,j);
      logfileptr->OFS()<<ss.str();

#ifdef USE_ABORT
      printf("%s",ss.str().c_str());
      abort();
#endif
			throw std::logic_error( ss.str().c_str() );
		}
    return data_[i+j*m_];
  }

template <typename F> inline F& NumMat<F>::operator()(Int i, Int j)  { 
		if( i < 0 || i >= m_ ||
				j < 0 || j >= n_ ) {

      std::stringstream ss;
      //error_message(ss,i,j);

      logfileptr->OFS()<<ss.str();

#ifdef USE_ABORT
      printf("%s",ss.str().c_str());
      abort();
#endif
			throw std::logic_error( ss.str().c_str() );
		}
    return data_[i+j*m_];
  }
  
  template <typename F> F* NumMat<F>::Data() const { return data_; }
  template <typename F> F* NumMat<F>::VecData(Int j)  const 
	{ 
		if( j < 0 || j >= n_ ) {

      std::stringstream ss;
      //error_message(ss,j);

#ifdef USE_ABORT
      printf("%s",ss.str().c_str());
      abort();
#endif
			throw std::logic_error( ss.str().c_str() );
		}
		return &(data_[j*m_]); 
	}


template <typename F> inline void SetValue(NumMat<F>& M, F val)
{
	F *ptr = M.data_;
	for (Int i=0; i < M.m()*M.n(); i++) *(ptr++) = val;
}

template <typename F> void
Transpose ( const NumMat<F>& A, NumMat<F>& B )
{
	if( A.m() != B.n() || A.n() != B.m() ){
		B.Resize( A.n(), A.m() );
	}

	F* Adata = A.Data();
	F* Bdata = B.Data();
	Int m = A.m(), n = A.n();

	for( Int i = 0; i < m; i++ ){
		for( Int j = 0; j < n; j++ ){
			Bdata[ j + n*i ] = Adata[ i + j*m ];
		}
	}


	return ;
}		// -----  end of function Transpose  ----- 

template <typename F> void
Symmetrize( NumMat<F>& A )
{
	if( A.m() != A.n() ){
		throw std::logic_error( "The matrix to be symmetrized should be a square matrix." );
	}

	NumMat<F> B;
	Transpose( A, B );

	F* Adata = A.Data();
	F* Bdata = B.Data();

	F  half = (F) 0.5;

	for( Int i = 0; i < A.m() * A.n(); i++ ){
		*Adata = half * (*Adata + *Bdata);
		Adata++; Bdata++;
	}


	return ;
}		// -----  end of function Symmetrize ----- 
}


#endif // _NUMMAT_IMPL_HPP_
