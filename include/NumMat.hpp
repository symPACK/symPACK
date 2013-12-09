/// @file nummat_decl.hpp
/// @brief Numerical matrix.
/// @author Lexing Ying and Lin Lin
/// @date 2010-09-27
#ifndef _NUMMAT_DECL_HPP_
#define _NUMMAT_DECL_HPP_

#include "Environment.hpp"

#ifdef UPCXX
#include <upcxx.h>
#endif

namespace LIBCHOLESKY{

/// @class NumMat
///
/// @brief Numerical matrix.
///
/// NumMat is a portable encapsulation of a pointer to represent a 2D
/// matrix, which can either own (owndata == true) or view (owndata ==
/// false) a piece of data.  
template <class F>
class NumMat
{
public:
	/// @brief The size of the first dimension.
  Int m_; 
	
	/// @brief The size of second dimension.
	Int n_;
	
	/// @brief Whether it owns the data.
  bool owndata_;
	
	/// @brief The pointer for the actual data.
  F* data_;

#ifdef UPCXX
  upcxx::global_ptr<F> gdata_;
#endif

  inline void alloc_data();
  inline void delete_data();

public:
  NumMat(Int m=0, Int n=0);
  NumMat(Int m, Int n, bool owndata, F* data);
  NumMat(const NumMat& C);
  ~NumMat();
  NumMat& Copy(const NumMat& C);
  void Resize(Int m, Int n);
  void Clear();

  NumMat& operator=(const NumMat& C);
  inline const F& operator()(Int i, Int j) const;
  inline F& operator()(Int i, Int j);

 
#ifdef UPCXX
  upcxx::global_ptr<F> GData() const;
  upcxx::global_ptr<F> GVecData(Int j) const;
#endif
  F* Data() const;
  F* VecData(Int j) const;
  Int m() const { return m_; }
  Int n() const { return n_; }
};

// Commonly used
typedef NumMat<bool>     BolNumMat;
typedef NumMat<Int>      IntNumMat;
typedef NumMat<Real>     DblNumMat;
typedef NumMat<Complex>  CpxNumMat;

// *********************************************************************
// Utility functions
// *********************************************************************
/// @brief SetValue sets a numerical matrix to a constant val.
template <class F> inline void SetValue(NumMat<F>& M, F val);

/// @brief Energy computes the L2 norm of a matrix (treated as a vector).
template <class F> Real Energy(const NumMat<F>& M);
 
template <class F> inline void Transpose ( const NumMat<F>& A, NumMat<F>& B );
template <class F> inline void Symmetrize( NumMat<F>& A );
} // namespace LIBCHOLESKY

#include "NumMat_impl.hpp"

#endif // _NUMMAT_DECL_HPP_
