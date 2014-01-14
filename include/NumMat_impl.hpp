/// @file nummat_decl.hpp
/// @brief Numerical matrix.
/// @author Lexing Ying and Lin Lin
/// @date 2010-09-27
#ifndef _NUMMAT_IMPL_HPP_
#define _NUMMAT_IMPL_HPP_

#include "Environment.hpp"
#include <sstream>

#ifdef UPCXX
using namespace upcxx;
#endif

namespace LIBCHOLESKY{


  template <class F> NumMat<F>::NumMat(Int m, Int n): m_(m), n_(n), owndata_(true) {
      alloc_data();
  }



  template <class F> NumMat<F>::NumMat(Int m, Int n, bool owndata, F* data): m_(m), n_(n), owndata_(owndata) {
    if(owndata_) {

      alloc_data();

      if(m_>0 && n_>0) { 
        std::copy(data,data+m_*n_,data_);//for(Int i=0; i<m_*n_; i++) data_[i] = data[i]; 
      }
    } 
    else {
#ifdef UPCXX
      gdata_= global_ptr<F>(data);
#endif
      data_ = data;
    }
  }



  template <class F> NumMat<F>::NumMat(const NumMat& C): m_(C.m_), n_(C.n_), owndata_(C.owndata_) {
    if(owndata_) {
      alloc_data();

      if(m_>0 && n_>0) { 
        std::copy(C.data_,C.data_+m_*n_,data_);//for(Int i=0; i<m_*n_; i++) data_[i] = data[i]; 
      }
    } 
    else {
#ifdef UPCXX
      gdata_= global_ptr<F>(C.data_);
#endif

      data_ = C.data_;
    }
  }


  template <class F> inline void NumMat<F>::alloc_data(){
    if(owndata_) {
      if(m_>0 && n_>0) { 
#ifdef UPCXX
      gdata_ = allocate<F>(MYTHREAD,m_*n_);     
      data_ = (F*)gdata_; 
#else
        data_ = new F[m_*n_]; 
#endif
        if( data_ == NULL ) 
          throw std::runtime_error("Cannot allocate memory."); 
      } 
      else 
        data_=NULL;
    }
  }

  template <class F> inline void NumMat<F>::delete_data(){
    if(owndata_) {
      if(m_>0 && n_>0) {
#ifdef UPCXX
  deallocate<F>(gdata_);
#else
        delete[] data_;
#endif
        data_ = NULL; 
      }
      m_=0;
      n_=0;
    }
  }


  template <class F> NumMat<F>::~NumMat() {
    delete_data();
  }

template <class F> NumMat<F>& NumMat<F>::Copy(const NumMat& C) {
    delete_data();

    m_ = C.m_; n_=C.n_; owndata_=C.owndata_;

    if(owndata_) {
      alloc_data();
      if(m_>0 && n_>0) { 
        std::copy(C.data_,C.data_+m_*n_,data_);//for(Int i=0; i<m_*n_; i++) data_[i] = data[i]; 
      }
    } 
    else {
#ifdef UPCXX
      gdata_= global_ptr<F>(C.data_);
#endif
      data_ = C.data_;
    }

    return *this;
  }


  template <class F> NumMat<F>& NumMat<F>::operator=(const NumMat& C) {

    delete_data();
    m_ = C.m_; n_=C.n_; owndata_=C.owndata_;

    if(owndata_) {
      alloc_data();
      if(m_>0 && n_>0) { 
        std::copy(C.data_,C.data_+m_*n_,data_);//for(Int i=0; i<m_*n_; i++) data_[i] = data[i]; 
      }
    } 
    else {
#ifdef UPCXX
      gdata_= global_ptr<F>(C.data_);
#endif
      data_ = C.data_;
    }

    return *this;
  }


  template <class F> void NumMat<F>::Clear()  {
		if( owndata_ == false ){
			throw std::logic_error("Matrix being cleared must own data.");
		}
      delete_data();
  }




  template <class F> void NumMat<F>::Resize(Int m, Int n)  {
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



template <class F> inline const F& NumMat<F>::operator()(Int i, Int j) const  { 
		if( i < 0 || i >= m_ ||
				j < 0 || j >= n_ ) {
      
      std::stringstream ss;

#ifdef UPCXX
      ss<<"P"<<MYTHREAD<<" ";
#endif
      ss<<"Index ("<<i<<","<<j<<") is out of bound. Dimensions are ("<<m_<<","<<n_<<")."<<std::endl;
			throw std::logic_error( ss.str().c_str() );
		}
    return data_[i+j*m_];
  }

template <class F> inline F& NumMat<F>::operator()(Int i, Int j)  { 
		if( i < 0 || i >= m_ ||
				j < 0 || j >= n_ ) {

      std::stringstream ss;

#ifdef UPCXX
      ss<<"P"<<MYTHREAD<<" ";
#endif
      ss<<"Index ("<<i<<","<<j<<") is out of bound. Dimensions are ("<<m_<<","<<n_<<")."<<std::endl;
			throw std::logic_error( ss.str().c_str() );
		}
    return data_[i+j*m_];
  }
  
  template <class F> F* NumMat<F>::Data() const { return data_; }
  template <class F> F* NumMat<F>::VecData(Int j)  const 
	{ 
		if( j < 0 || j >= n_ ) {

      std::stringstream ss;

#ifdef UPCXX
      ss<<"P"<<MYTHREAD<<" ";
#endif
      ss<<"Index (*,"<<j<<") is out of bound. Dimensions are ("<<m_<<","<<n_<<")."<<std::endl;
			throw std::logic_error( ss.str().c_str() );
		}
		return &(data_[j*m_]); 
	}

#ifdef UPCXX
  template <class F> global_ptr<F> NumMat<F>::GData() const { return gdata_; }
  template <class F> global_ptr<F> NumMat<F>::GVecData(Int j)  const 
	{ 
		if( j < 0 || j >= n_ ) {

      std::stringstream ss;

#ifdef UPCXX
      ss<<"P"<<MYTHREAD<<" ";
#endif
      ss<<"Index (*,"<<j<<") is out of bound. Dimensions are ("<<m_<<","<<n_<<")."<<std::endl;
			throw std::logic_error( ss.str().c_str() );
		}
		return (global_ptr<F>(&data_[j*m_])); 
	}
#endif


template <class F> inline void SetValue(NumMat<F>& M, F val)
{
	F *ptr = M.data_;
	for (Int i=0; i < M.m()*M.n(); i++) *(ptr++) = val;
}

template <class F> void
Transpose ( const NumMat<F>& A, NumMat<F>& B )
{
#ifndef _RELEASE_
	PushCallStack("Transpose");
#endif
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

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function Transpose  ----- 

template <class F> void
Symmetrize( NumMat<F>& A )
{
#ifndef _RELEASE_
	PushCallStack("Symmetrize");
#endif
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

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function Symmetrize ----- 
}

#endif // _NUMMAT_IMPL_HPP_
