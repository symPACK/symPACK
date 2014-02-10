/// @file NumMat_upcxx_impl.hpp
/// @brief Numerical matrix.
/// @author Mathias Jacquelin
/// @date 2014-02-05
#ifndef _NUMMAT_UPCXX_IMPL_HPP_
#define _NUMMAT_UPCXX_IMPL_HPP_

#include "Environment.hpp"
#include <sstream>

using namespace upcxx;

namespace LIBCHOLESKY{



  template <typename F> NumMat_upcxx<F>::NumMat_upcxx(Int m, Int n, bool owndata, F* data) {
#ifdef _ASSERT_
    this->allocated_=false;
#endif

      this->m_=m;
      this->n_=n;
      this->owndata_=true;

    if(this->owndata_) {

      alloc_data();

      if(this->m_>0 && this->n_>0) { 
        std::copy(data,data+this->m_*this->n_,this->data_);
      }
    } 
    else {
      gdata_= global_ptr<F>(data);
      this->data_ = data;
    }
  }


  template <typename F> NumMat_upcxx<F>::NumMat_upcxx(Int m, Int n) {
#ifdef _ASSERT_
    this->allocated_=false;
#endif
      this->m_=m;
      this->n_=n;
      this->owndata_=true;

      alloc_data();
  }

  template <typename F> NumMat_upcxx<F>::NumMat_upcxx(const NumMat_upcxx<F>& C) {

#ifdef _ASSERT_
    this->allocated_=false;
#endif

    this->m_ = C.m_;
    this->n_=C.n_;
    this->owndata_=C.owndata_;

    if(this->owndata_) {
      alloc_data();
      if(this->m_>0 && this->n_>0) { 
        std::copy(C.data_,C.data_+this->m_*this->n_,this->data_);
      }
    } 
    else {
      gdata_= global_ptr<F>(C.data_);
      this->data_ = C.data_;
    }
  }


  template <typename F> NumMat_upcxx<F>::~NumMat_upcxx() {
    delete_data();
  }


  template <typename F> inline void NumMat_upcxx<F>::alloc_data(){
    //logfileptr->OFS()<<"ALLOCATING UPCXX DATA"<<std::endl;

#ifdef _ASSERT_
    if(!this->allocated_) {
      if(this->owndata_) {
        if(this->m_>0 && this->n_>0) { 
          gdata_ = allocate<F>(MYTHREAD,this->m_*this->n_);     
          this->data_ = (F*)gdata_; 

          if( this->data_ == NULL ) 
            throw std::runtime_error("Cannot allocate memory."); 
        } 
        else 
          this->data_=NULL;

        this->allocated_=true;
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
      if(this->owndata_) {
        if(this->m_>0 && this->n_>0) { 
          gdata_ = allocate<F>(MYTHREAD,this->m_*this->n_);     
          this->data_ = (F*)gdata_; 

          if( this->data_ == NULL ) 
            throw std::runtime_error("Cannot allocate memory."); 
        } 
        else 
          this->data_=NULL;

      }
#endif
  }

  template <typename F> inline void NumMat_upcxx<F>::delete_data(){

#ifdef _ASSERT_
    if(this->owndata_) {
      if(this->allocated_) {
        if(this->m_>0 && this->n_>0) {
          deallocate<F>(gdata_);
          this->data_ = NULL; 
        }
        this->allocated_=false;
      }
      this->m_=0;
      this->n_=0;
    }
#else
    if(this->owndata_) {
      if(this->m_>0 && this->n_>0) {
        deallocate<F>(gdata_);
        this->data_ = NULL; 
      }
      this->m_=0;
      this->n_=0;
    }
#endif
  }



template <typename F> NumMat_upcxx<F>& NumMat_upcxx<F>::Copy(const NumMat_upcxx<F>& C) {
    delete_data();

    this->m_ = C.m_; this->n_=C.n_; this->owndata_=C.owndata_;

    if(this->owndata_) {
      alloc_data();
      if(this->m_>0 && this->n_>0) { 
        std::copy(C.data_,C.data_+this->m_*this->n_,this->data_);
      }
    } 
    else {
      gdata_= global_ptr<F>(C.data_);
      this->data_ = C.data_;
    }

    return *this;
  }






  template<typename F> void NumMat_upcxx<F>::error_message(std::stringstream & ss,Int i, Int j){
      ss<<"P"<<MYTHREAD<<" ";
      NumMat<F>::error_message(ss,i,j);
  }
  template<typename F> void NumMat_upcxx<F>::error_message(std::stringstream & ss, Int j){
      ss<<"P"<<MYTHREAD<<" ";
      NumMat<F>::error_message(ss,j);
  }






  template <typename F> global_ptr<F> NumMat_upcxx<F>::GData() const { return gdata_; }

  template <typename F> global_ptr<F> NumMat_upcxx<F>::GVecData(Int j)  const 
	{ 
		if( j < 0 || j >= this->n_ ) {

      std::stringstream ss;

      error_message(ss,j);
      
#ifdef USE_ABORT
      printf("%s",ss.str().c_str());
      abort();
#endif

			throw std::logic_error( ss.str().c_str() );
		}
		return (global_ptr<F>(&this->data_[j*this->m_])); 
	}

template <typename F> void NumMat_upcxx<F>::Clear()  {
		if(this->owndata_ == false ){
			throw std::logic_error("Matrix being cleared must own data.");
		}
      delete_data();
  }




  template <typename F> void NumMat_upcxx<F>::Resize(Int m, Int n)  {
		if( this->owndata_ == false ){
			throw std::logic_error("Matrix being resized must own data.");
		}
    if(this->m_!=m || this->n_!=n) {
      delete_data();
      this->m_ = m; this->n_ = n;
      if(m*n>0){
        alloc_data();
      }
    }
  }



}

#endif // _NUMMAT_UPCXX_IMPL_HPP_
