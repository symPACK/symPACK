/// @file numvec_decl.hpp
/// @brief  Numerical vector.
/// @author Lexing Ying and Lin Lin
/// @date 2010-09-27
#ifndef _NUMVEC_DECL_HPP_
#define _NUMVEC_DECL_HPP_

#include "Environment.hpp"

namespace LIBCHOLESKY{

/// @class NumVec
///
/// @brief Numerical vector.
/// 
/// NumVec is a portable encapsulation of a pointer to represent a 1D
/// vector. The main difference between NumVec<F> and std::vector<F> is
/// that NumVec<F> allows the vector to not owning the data, by
/// specifying (owndata_ == false).
template <class F> class NumVec
{
public:
	/// @brief The size of the vector.
	Int  m_;                                
	///
	/// @brief Whether it owns the data.
	bool owndata_;                          

	/// @brief The pointer for the actual data.
	F* data_;                                
public:
	NumVec(Int m = 0);
	NumVec(Int m, bool owndata, F* data);
	NumVec(const NumVec& C);
	~NumVec();

	NumVec& operator=(const NumVec& C);

	void Resize ( Int m );

	const F& operator()(Int i) const;  
	F& operator()(Int i);  
	const F& operator[](Int i) const;
	F& operator[](Int i);

	bool IsOwnData() const { return owndata_; }
	F*   Data() const { return data_; }
	Int  m () const { return m_; }
};

// Commonly used
typedef NumVec<bool>       BolNumVec;
typedef NumVec<Int>        IntNumVec;
typedef NumVec<Real>       DblNumVec;
typedef NumVec<Complex>    CpxNumVec;


// *********************************************************************
// Utility functions
// *********************************************************************
/// @brief SetValue sets a numerical vector to a constant val.
template <class F> void SetValue( NumVec<F>& vec, F val );

/// @brief Energy computes the L2 norm of a vector.
template <class F> Real Energy( const NumVec<F>& vec );


} // namespace LIBCHOLESKY

#include "NumVec_impl.hpp"

#endif // _NUMVEC_DECL_HPP_
