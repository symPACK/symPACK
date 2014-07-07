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


template <class F> NumVec<F>::NumVec	( Int m ) : m_(m), owndata_(true)
{
	if(m_>0) { 
//		data_ = new F[m_]; 
		data_ = (F*) malloc(m_*sizeof(F));
		if( data_ == NULL ){

#ifdef USE_ABORT
      printf("%s","Cannot allocate memory.");
      abort();
#endif
			throw std::runtime_error("Cannot allocate memory.");
		}
	} 
	else 
		data_=NULL;
} 		// -----  end of method NumVec<F>::NumVec  ----- 

template <class F> NumVec<F>::NumVec	( Int m, bool owndata, F* data ) : m_(m), owndata_(owndata)
{
	if( owndata_ ){
		if( m_ > 0 ) { 
//			data_ = new F[m_]; 
  		data_ = (F*) malloc(m_*sizeof(F));
			if( data_ == NULL ){

#ifdef USE_ABORT
      printf("%s","Cannot allocate memory.");
      abort();
#endif
				throw std::runtime_error("Cannot allocate memory.");
			}
		}
		else
			data_ = NULL;

		if( m_ > 0 ) {
			for( Int i = 0; i < m_; i++ ){
				data_[i] = data[i];
			}
		}
	}
	else{
		data_ = data;
	}
} 		// -----  end of method NumVec<F>::NumVec  ----- 

template <class F> NumVec<F>::NumVec	( const NumVec<F>& C ) : m_(C.m_), owndata_(C.owndata_)
{
	if( owndata_ ){
		if( m_ > 0 ) { 
			//data_ = new F[m_]; 
		  data_ = (F*) malloc(m_*sizeof(F));
			if( data_ == NULL ){

#ifdef USE_ABORT
      printf("%s","Cannot allocate memory.");
      abort();
#endif
				throw std::runtime_error("Cannot allocate memory.");
			}
		}
		else
			data_ = NULL;

		if( m_ > 0 ) {
			for( Int i = 0; i < m_; i++ ){
				data_[i] = C.data_[i];
			}
		}
	}
	else{
		data_ = C.data_;
	}
} 		// -----  end of method NumVec<F>::NumVec  ----- 


template < class F > NumVec<F>::~NumVec	(  )
{
	if( owndata_ ){
		if( m_ > 0 ){
//			delete[] data_;  
      free(data_);
			data_ = NULL;
		}
	}

} 		// -----  end of method NumVec<F>::~NumVec  ----- 


template < class F > NumVec<F>& NumVec<F>::operator =	( const NumVec& C  )
{
	if( owndata_ ){
		if( m_ > 0 ){
//			delete[]  data_;
      free(data_);
			data_ = NULL;
		}
	}
	m_ = C.m_;
	owndata_ = C.owndata_;

	if( owndata_ ) {
		if( m_ > 0 ){
//			data_ = new F[m_];
		  data_ = (F*) malloc(m_*sizeof(F));
			if( data_ == NULL ){

#ifdef USE_ABORT
      printf("%s","Cannot allocate memory.");
      abort();
#endif
				throw std::runtime_error("Cannot allocate memory.");
			}
		}
		else{
			data_ = NULL;
		}

		if( m_ > 0 ){
			for( Int i = 0; i < m_; i++ ){
				data_[i] = C.data_[i];
			}
		}
	}
	else{
		data_ = C.data_;
	}


	return *this;
} 		// -----  end of method NumVec<F>::operator=  ----- 


template < class F > void NumVec<F>::Resize	( const Int m )
{
	if( owndata_ == false ){

#ifdef USE_ABORT
      printf("%s","Vector being resized must own data.");
      abort();
#endif
		throw std::logic_error("Vector being resized must own data.");
	}
	if( m != m_ ){
//    if( m_ > 0 ){
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
//    }

//		if( m_ > 0 ){
//			delete[] data_;
//			data_ = NULL;
//		}
//		m_ = m;
//		if( m_ > 0 ){
//			data_ = new F[m_];
//			if( data_ == NULL ){
//				throw std::runtime_error("Cannot allocate memory.");
//			}
//		}
	}

	return ;
} 		// -----  end of method NumVec<F>::Resize  ----- 


template <class F> F& NumVec<F>::operator()	( Int i )
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

} 		// -----  end of method NumVec<F>::operator()  ----- 


template <class F> const F&
NumVec<F>::operator()	( Int i ) const
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

} 		// -----  end of method NumVec<F>::operator()  ----- 


template <class F> F& NumVec<F>::operator[]	( Int i )
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
} 		// -----  end of method NumVec<F>::operator[]  ----- 


template <class F> const F& NumVec<F>::operator[]	( Int i ) const
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

} 		// -----  end of method NumVec<F>::operator[]  ----- 



template <class F> void SetValue( NumVec<F>& vec, F val )
{
	for(Int i=0; i<vec.m(); i++)
		vec(i) = val;
}

template <class F> Real Energy( const NumVec<F>& vec )
{
	Real sum = 0;
	for(Int i=0; i<vec.m(); i++)
		sum += abs(vec(i)*vec(i));
	return sum;
}  
}

#endif
