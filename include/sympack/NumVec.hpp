/// @file numvec_decl.hpp
/// @brief  Numerical vector.
/// @author Lexing Ying and Lin Lin
/// @date 2010-09-27
#ifndef _NUMVEC_DECL_HPP_
#define _NUMVEC_DECL_HPP_

#include "sympack/Environment.hpp"

namespace SYMPACK{

/// @class NNumVec
///
/// @brief Numerical vector.
/// 
/// NumVec is a portable encapsulation of a pointer to represent a 1D
/// vector. The main difference between NumVec<F> and std::vector<F> is
/// that NumVec<F> allows the vector to not owning the data, by
/// specifying (owndata_ == false).

template <class F, class TIdx = Int > class NumVec
{
public:
	/// @brief The size of the vector.
	TIdx  m_;                                
	///
	/// @brief Whether it owns the data.
	bool owndata_;                          

	/// @brief The pointer for the actual data.
	F* data_;        
                        
  void alloc_data();
public:
	NumVec(TIdx m = 0);
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
};

// Commonly used
//template <class F> using NumVec = NumVec<F,Int>;

typedef NumVec<bool>       BolNumVec;
typedef NumVec<Int>        IntNumVec;
typedef NumVec<Real>       DblNumVec;
typedef NumVec<Complex>    CpxNumVec;

//  typedef NumVec<Idx64,Idx64> PtrVec;
//  typedef NumVec<Idx64,Idx64> IdxVec;

// *********************************************************************
// Utility functions
// *********************************************************************
/// @brief SetValue sets a numerical vector to a constant val.
template <class F, class TIdx> void SetValue( NumVec<F,TIdx>& vec, F val );


} // namespace SYMPACK

#include "sympack/NumVec_impl.hpp"

#endif // _NUMVEC_DECL_HPP_
