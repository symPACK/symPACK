/// @file numvec_impl.hpp
/// @brief Implementation of Numerical Vector.
/// @author Lexing Ying and Lin Lin
/// @date 2010-09-27
#ifndef _NUMVEC_IMPL_HPP_
#define _NUMVEC_IMPL_HPP_

namespace LIBCHOLESKY{

// Templated form of numerical vectors
//
// The main advantage of this portable NumVec structure is that it can
// either own (owndata == true) or view (owndata == false) a piece of
// data.

  template <class F, class TIdx> inline void NumVec<F,TIdx>::alloc_data()
{
#ifdef _ASSERT_
    if(!allocated_) {
      if(owndata_) {
        if(m_>0) { 
		      data_ = (F*) malloc(m_*sizeof(F));
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
        if(m_>0 ) { 
		      data_ = (F*) malloc(m_*sizeof(F));
          if( data_ == NULL ) 
            throw std::runtime_error("Cannot allocate memory."); 
        } 
        else 
          data_=NULL;
      }
#endif
  }



template <class F, class TIdx> NumVec<F,TIdx>::NumVec	( TIdx m ) : m_(m), owndata_(true)
{
  alloc_data();
} 		// -----  end of method NumVec<F,TIdx>::NumVec  ----- 

template <class F, class TIdx> NumVec<F,TIdx>::NumVec	( TIdx m, bool owndata, F* data ) : m_(m), owndata_(owndata)
{
	if( owndata_ ){
    alloc_data();
		if( m_ > 0 ) {
      std::copy(data,data+m_,data_);
		}
	}
	else{
		data_ = data;
	}
} 		// -----  end of method NumVec<F,TIdx>::NumVec  ----- 

template <class F, class TIdx> NumVec<F,TIdx>::NumVec	( const NumVec<F,TIdx>& C ) : m_(C.m_), owndata_(C.owndata_)
{
	if( owndata_ ){
    alloc_data();
		if( m_ > 0 ) {
      std::copy(C.data_,C.data_+m_,data_);
		}
	}
	else{
		data_ = C.data_;
	}
} 		// -----  end of method NumVec<F,TIdx>::NumVec  ----- 


template < class F, class TIdx > NumVec<F,TIdx>::~NumVec	(  )
{
	if( owndata_ ){
		if( m_ > 0 ){
      free(data_);
			data_ = NULL;
		}
	}

} 		// -----  end of method NumVec<F,TIdx>::~NumVec  ----- 


template < class F, class TIdx > NumVec<F,TIdx>& NumVec<F,TIdx>::operator =	( const NumVec& C  )
{
	if( owndata_ ){
		if( m_ > 0 ){
      free(data_);
			data_ = NULL;
		}
	}
	m_ = C.m_;
	owndata_ = C.owndata_;

	if( owndata_ ) {
    alloc_data();
		if( m_ > 0 ){
      std::copy(C.data_,C.data_+m_,data_);
		}
	}
	else{
		data_ = C.data_;
	}


	return *this;
} 		// -----  end of method NumVec<F,TIdx>::operator=  ----- 


template < class F, class TIdx > void NumVec<F,TIdx>::Resize	( const TIdx m )
{
	if( owndata_ == false ){

#ifdef USE_ABORT
      printf("%s","Vector being resized must own data.");
      abort();
#endif
		throw std::logic_error("Vector being resized must own data.");
	}
	if( m != m_ ){
      F* newdata = (F*) realloc((void *)data_, m*sizeof(F));
      if (newdata==NULL){
        free(data_);

#ifdef USE_ABORT
      printf("%s","Cannot reallocate memory.");
      abort();
#endif
        throw std::runtime_error("Cannot reallocate memory.");
      }
      else{
        data_ = newdata;
        m_ = m;
      }
	}

	return ;
} 		// -----  end of method NumVec<F,TIdx>::Resize  ----- 


template <class F, class TIdx> F& NumVec<F,TIdx>::operator()	( TIdx i )
{
#ifndef _RELEASE_
	if( i < 0 || i >= m_ ){

#ifdef USE_ABORT
      printf("Index %d is out of bound.",i);
      abort();
#endif
		throw std::logic_error( "Index is out of bound." );
	}
#endif  // ifndef _RELEASE_
	return data_[i];

} 		// -----  end of method NumVec<F,TIdx>::operator()  ----- 


template <class F, class TIdx> const F&
NumVec<F,TIdx>::operator()	( TIdx i ) const
{
#ifndef _RELEASE_
	if( i < 0 || i >= m_ ){

#ifdef USE_ABORT
      printf("Index %d is out of bound.",i);
      abort();
#endif
		throw std::logic_error( "Index is out of bound." );
	}
#endif  // ifndef _RELEASE_
	return data_[i];

} 		// -----  end of method NumVec<F,TIdx>::operator()  ----- 


template <class F, class TIdx> F& NumVec<F,TIdx>::operator[]	( TIdx i )
{
#ifndef _RELEASE_
	if( i < 0 || i >= m_ ){

#ifdef USE_ABORT
      printf("Index %d is out of bound.",i);
      abort();
#endif
		throw std::logic_error( "Index is out of bound." );
	}
#endif  // ifndef _RELEASE_
	return data_[i];
} 		// -----  end of method NumVec<F,TIdx>::operator[]  ----- 


template <class F, class TIdx> const F& NumVec<F,TIdx>::operator[]	( TIdx i ) const
{
#ifndef _RELEASE_
	if( i < 0 || i >= m_ ){

#ifdef USE_ABORT
      printf("Index %d is out of bound.",i);
      abort();
#endif
		throw std::logic_error( "Index is out of bound." );
	}
#endif  // ifndef _RELEASE_
	return data_[i];

} 		// -----  end of method NumVec<F,TIdx>::operator[]  ----- 


 template <class F, class TIdx> void NumVec<F,TIdx>::Clear()  
{
		if( owndata_ == false ){
			throw std::logic_error("Vector being cleared must own data.");
		}

      if (data_==NULL){
        free(data_);
      }
  }



template <class F, class TIdx> void SetValue( NumVec<F,TIdx>& vec, F val )
{
  std::fill(&vec[0],&vec[0]+vec.m(),val);
}

}

#endif
