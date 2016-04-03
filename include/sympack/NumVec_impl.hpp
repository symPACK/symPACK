/// @file numvec_impl.hpp
/// @brief Implementation of Numerical Vector.
/// @author Lexing Ying and Lin Lin
/// @date 2010-09-27
#ifndef _NUMVEC_IMPL_HPP_
#define _NUMVEC_IMPL_HPP_

namespace SYMPACK{

// Templated form of numerical SYMPACK::vectors
//
// The main advantage of this portable NumVec structure is that it can
// either own (owndata == true) or view (owndata == false) a piece of
// data.

  template <class F, class TIdx> inline void NumVec<F,TIdx>::alloc_data()
  {
    if(owndata_) {
      if(m_>0 ) { 
//        data_ = (F*) malloc(m_*sizeof(F));
        data_ = (F*) new F[m_];
        if( data_ == NULL ) 
          throw std::runtime_error("Cannot allocate memory."); 
      } 
      else 
        data_=NULL;
    }
    begin_.p = data_;
    end_.p = data_+head_;
  }


template <class F, class TIdx> NumVec<F,TIdx>::NumVec	() : m_(0), owndata_(true),head_(0),data_(NULL){
      begin_.p = data_;
      end_.p = data_+head_;
}

template <class F, class TIdx> NumVec<F,TIdx>::NumVec	( TIdx m ) : m_(m), owndata_(true),head_(m_)
{
  alloc_data();
} 		// -----  end of method NumVec<F,TIdx>::NumVec  ----- 

template <class F, class TIdx> NumVec<F,TIdx>::NumVec	( TIdx m, bool owndata, F* data ) : m_(m), owndata_(owndata),head_(m_)
{
	if( owndata_ ){
    alloc_data();
		if( m_ > 0 ) {
      std::copy(data,data+m_,data_);
		}
	}
	else{
		data_ = data;
      begin_.p = data_;
      end_.p = data_+head_;
	}
} 		// -----  end of method NumVec<F,TIdx>::NumVec  ----- 

template <class F, class TIdx> NumVec<F,TIdx>::NumVec	( const NumVec<F,TIdx>& C ) : m_(C.m_), owndata_(C.owndata_),head_(C.head_)
{
	if( owndata_ ){
    alloc_data();
		if( m_ > 0 ) {
      std::copy(C.data_,C.data_+m_,data_);
		}
	}
	else{
		data_ = C.data_;
      begin_.p = data_;
      end_.p = data_+head_;
	}
} 		// -----  end of method NumVec<F,TIdx>::NumVec  ----- 


template < class F, class TIdx > NumVec<F,TIdx>::~NumVec	(  )
{
	if( owndata_ ){
		if( m_ > 0 ){
      //free(data_);
      delete [] data_;
		}
	}

} 		// -----  end of method NumVec<F,TIdx>::~NumVec  ----- 


template < class F, class TIdx > NumVec<F,TIdx>& NumVec<F,TIdx>::operator =	( const NumVec& C  )
{
	if( owndata_ ){
		if( m_ > 0 ){
      //free(data_);
      delete [] data_;
			data_ = NULL;
		}
	}
	m_ = C.m_;
	owndata_ = C.owndata_;
  head_ = C.head_;
	if( owndata_ ) {
    alloc_data();
		if( m_ > 0 ){
      std::copy(C.data_,C.data_+m_,data_);
		}
	}
	else{
		data_ = C.data_;
      begin_.p = data_;
      end_.p = data_+head_;
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
      //F* newdata = (F*) realloc((void *)data_, m*sizeof(F));
      F* newdata = (F*) new F[m];
      if (newdata==NULL){
        //free(data_);
        delete [] data_;
			  data_ = NULL;

#ifdef USE_ABORT
      printf("%s","Cannot reallocate memory.");
      abort();
#endif
        throw std::runtime_error("Cannot reallocate memory.");
      }
      else{
        std::copy(data_,data_+min(m,m_),newdata);
        delete [] data_;
        data_ = newdata;
        m_ = m;
        head_=m_;
      begin_.p = data_;
      end_.p = data_+head_;
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

      if (data_!=NULL){
        //free(data_);
        delete [] data_;
        m_ = 0;
        head_=0;
        data_=NULL;
        begin_.p = data_;
        end_.p = data_+head_;
      }
  }



template <class F, class TIdx> void SetValue( NumVec<F,TIdx>& vec, F val )
{
  std::fill(&vec[0],&vec[0]+vec.m(),val);
}

}

#endif
