/// @file numvec_decl.hpp
/// @brief  Numerical SYMPACK::vector.
/// @author Lexing Ying and Lin Lin
/// @date 2010-09-27
#ifndef _NUMVEC_DECL_HPP_
#define _NUMVEC_DECL_HPP_

#include <iterator>
#include "sympack/Environment.hpp"

namespace SYMPACK{

/// @class NNumVec
///
/// @brief Numerical SYMPACK::vector.
/// 
/// NumVec is a portable encapsulation of a pointer to represent a 1D
/// SYMPACK::vector. The main difference between NumVec<F> and SYMPACK::vector<F> is
/// that NumVec<F> allows the SYMPACK::vector to not owning the data, by
/// specifying (owndata_ == false).

template <class F, class TIdx  > class NumVec;
// *********************************************************************
// Utility functions
// *********************************************************************
/// @brief SetValue sets a numerical SYMPACK::vector to a constant val.
template <class F, class TIdx> void SetValue( NumVec<F,TIdx>& vec, F val );




template <class F, class TIdx = Int > class NumVec
{
public:
	/// @brief The size of the SYMPACK::vector.
	TIdx  m_;                                
	///
	/// @brief Whether it owns the data.
	bool owndata_;                          

	/// @brief The pointer for the actual data.
	F* data_;        
                        
  void alloc_data();
public:
	NumVec();
	NumVec(TIdx m);
	NumVec(TIdx m, bool owndata, F* data);
	NumVec(const NumVec& C);
	~NumVec();

	NumVec& operator=(const NumVec& C);

	void Resize ( TIdx m );
	void Clear();


	const F& operator()(TIdx i) const;  
	F& operator()(TIdx i);  
	const F& operator[](TIdx i) const;
	F& operator[](TIdx i);

	bool IsOwnData() const { return owndata_; }
	F*   Data() const { return data_; }
	TIdx  m () const { return m_; }

  //compatibility with std vector 
//	NumVec(TIdx m , int & value):NumVec(m){SetValue(*this,(F)value);};
	template <typename V> NumVec(TIdx m , V value):NumVec((TIdx)m){SetValue(*this,(F)value);};


  const F& at(TIdx i) const {return (*this)[i];}
  F& at(TIdx i) {return (*this)[i];}

  F& front()  {return data_[0];} 
  F& back()  {return  data_[head_-1];}
  const F& front() const {return const_cast<const F&>(data_[0]);} 
  const F& back() const  {return  const_cast<const F&>(data_[head_-1]);}
  TIdx size() const { return head_;}
  TIdx capacity() const { return m();}
  bool empty() const { return head_==0;}
	inline void resize ( TIdx m ){Resize(m);};
	inline void resize ( TIdx m, F value ){ TIdx prevHead = head_; Resize(m); std::fill(Data()+prevHead,Data()+head_,value); };
	inline void clear(){Clear();};
  inline void swap(NumVec &C){F* tmpData = data_; TIdx tmpM = m_;  TIdx tmpHead = head_; head_=C.head_; m_=C.m_; data_ = C.data_; C.head_=tmpHead; C.m_=tmpM; C.data_=tmpData;}
  inline void reserve(TIdx newSize){ TIdx prevHead = head_; Resize(newSize); head_ = prevHead; end_.p = data_+head_; };
  //should be initialized in the constructor
  TIdx head_;
  inline void push_back(F value){
    if(head_==capacity()){
      reserve(capacity()+1);
    }
    data_[head_++] = value;
  }; 
  inline void pop_back(){
      assert(size()>0);
      resize(size()-1);
      head_--;
  }

  using value_type = F;
  using size_type = TIdx;
  using difference_type = ptrdiff_t;
  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;
//  using iterator = std::iterator<std::random_access_iterator_tag,value_type>;
//  using const_iterator = std::iterator<std::random_access_iterator_tag,const value_type>;

  class iterator : public std::iterator<std::random_access_iterator_tag,value_type>{
    //template <value_type, class TIdx> 
    friend class NumVec;

    pointer p;
    public:
      iterator() :p(NULL) {}
      iterator(value_type* x) :p(x) {}
      //iterator
      iterator(const iterator& mit) : p(mit.p) {}
      iterator& operator=(const iterator &rhs) {p = rhs.p; return *this;}
      iterator& operator=(iterator &rhs) {p = rhs.p; return *this;}
      iterator& operator=(const pointer rhs) {p = rhs; return *this;}
      iterator& operator++() {++p;return *this;}
      reference operator*() const {return *p;}

      //input iterator
      iterator operator++(int) {iterator tmp(*this); operator++(); return tmp;}
      pointer operator->() const {return p;}
      friend bool operator==(const iterator& lhs, const iterator& rhs) {return lhs.p==rhs.p;}
      friend bool operator!=(const iterator& lhs, const iterator& rhs) {return lhs.p!=rhs.p;}

      //output iterator
      //nothing extra

      //bidirectional_iterator
      iterator& operator--() {--p;return *this;}
      iterator operator--(int) {iterator tmp(*this); operator--(); return tmp;}

      //random_access_iterator
      friend bool  operator>(const iterator &lhs, const iterator& rhs)  {return lhs.p > rhs.p;}
      friend bool  operator<(const iterator &lhs, const iterator& rhs)  {return lhs.p < rhs.p;}
      friend bool operator>=(const iterator &lhs, const iterator& rhs) {return lhs.p >= rhs.p;}
      friend bool operator<=(const iterator &lhs, const iterator& rhs) {return lhs.p <= rhs.p;}

      iterator& operator+=(size_type rhs){ p+=rhs; return *this;}
      friend iterator operator+(size_type lhs, const iterator& rhs) {return iterator(lhs+rhs.p);}
      friend iterator operator+(const iterator& lhs, size_type rhs) {return iterator(rhs+lhs.p);}


      iterator& operator-=(size_type rhs){ p-=rhs; return *this;}
      friend iterator operator-(const iterator& lhs, size_type rhs) {return iterator(lhs.p-rhs);}
      friend difference_type operator-(iterator lhs, iterator rhs) {return lhs.p-rhs.p;}

      reference operator[](size_type rhs) const {return p[rhs];}

  };
//
//   class const_iterator : public std::iterator<std::random_access_iterator_tag,const value_type>{
//    const value_type* p;
//    public:
//      const_iterator(const value_type* x) :p(x) {}
//      const_iterator(const const_iterator& mit) : p(mit.p) {}
//      const_iterator& operator++() {++p;return *this;}
//      const_iterator operator++(int) {const_iterator tmp(*this); operator++(); return tmp;}
//      const_iterator& operator--() {--p;return *this;}
//      const_iterator operator--(int) {const_iterator tmp(*this); operator--(); return tmp;}
//
//        inline const_iterator operator+(const const_iterator& rhs) {return const_iterator(p+rhs.p);}
//        inline const_iterator operator-(const const_iterator& rhs) {return const_iterator(p-rhs.p);}
//        inline const_iterator operator+(const int& rhs) {return const_iterator(p+rhs);}
//        inline const_iterator operator-(const difference_type& rhs) {return const_iterator(p-rhs);}
//        friend inline const_iterator operator+(const int& lhs, const const_iterator& rhs) {return const_iterator(lhs+rhs.p);}
//        friend inline const_iterator operator-(const int& lhs, const const_iterator& rhs) {return const_iterator(lhs-rhs.p);}
//
//
//
//        inline bool operator==(const const_iterator& rhs) {return p == rhs.p;}
//        inline bool operator!=(const const_iterator& rhs) {return p != rhs.p;}
//        inline bool operator>(const const_iterator& rhs)  {return p > rhs.p;}
//        inline bool operator<(const const_iterator& rhs)  {return p < rhs.p;}
//        inline bool operator>=(const const_iterator& rhs) {return p >= rhs.p;}
//        inline bool operator<=(const const_iterator& rhs) {return p <= rhs.p;}
//
//
//        inline const_iterator& operator=(const value_type* rhs) {p = rhs; return *this;}
//        inline const_iterator& operator=(const const_iterator &rhs) {p = rhs.p; return *this;}
//        inline const_iterator& operator+=(const difference_type& rhs) {p += rhs; return *this;}
//        inline const_iterator& operator-=(const difference_type& rhs) {p -= rhs; return *this;}
//        inline const value_type& operator*() {return *p;}
//        inline const value_type* operator->() {return p;}
//        inline const value_type& operator[](const int& rhs) {return p[rhs];}
//  };

  iterator begin_;
  iterator end_;
 
//  iterator& begin() {return (begin_);} 
//  iterator& end() {return (end_);} 
  iterator begin() {return iterator(Data());} 
  iterator end() {return iterator(Data()+head_);} 
  template <class InputIterator>
  void insert (iterator position, InputIterator first, InputIterator last){
    size_t numElem = std::distance(first,last);
    size_t pos = std::distance(begin(),position);
    iterator ibeg = begin();
    iterator iend = end();
    resize(m()+numElem);
    std::copy(position,iend,position+numElem);
    std::copy(first,last,position);
  }


  void assign (size_type n, const value_type& val){
    resize(n);
    SetValue(*this,val); 
  }

};

// Commonly used
//template <class F> using NumVec = NumVec<F,Int>;

typedef NumVec<bool>       BolNumVec;
typedef NumVec<Int>        IntNumVec;
typedef NumVec<Real>       DblNumVec;
typedef NumVec<Complex>    CpxNumVec;

//  typedef NumVec<Idx64,Idx64> PtrVec;
//  typedef NumVec<Idx64,Idx64> IdxVec;

} // namespace SYMPACK

#include "sympack/NumVec_impl.hpp"

#endif // _NUMVEC_DECL_HPP_
