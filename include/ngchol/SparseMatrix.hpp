/// @file SparseMatrix.hpp
/// @brief Sparse matrix and Distributed sparse matrix in compressed
/// column format.
/// @author Lin Lin
/// @date 2012-11-10
#ifndef _SPARSE_MATRIX_DECL_HPP_
#define _SPARSE_MATRIX_DECL_HPP_

#include "ngchol/Environment.hpp"
#include "ngchol/NumVec.hpp"

namespace LIBCHOLESKY{

/// @struct SparseMatrix
/// 
/// @brief SparseMatrix describes a sequential sparse matrix saved in
/// compressed sparse column format.
///
/// Note
/// ----
///
/// Since in PEXSI and PPEXSI only symmetric matrix is considered, the
/// compressed sparse row format will also be represented by the
/// compressed sparse column format.
template <class F> struct SparseMatrix{
	Int          size;                            // Matrix dimension
	Int          nnz;                             // Number of nonzeros
	IntNumVec    colptr;                          // Column index pointer
	IntNumVec    rowind;                          // Starting row index pointer
	NumVec<F>    nzval;                           // Nonzero values for the sparse matrix
};




// Commonly used
typedef SparseMatrix<Real>       DblSparseMatrix;
typedef SparseMatrix<Complex>    CpxSparseMatrix;



} // namespace LIBCHOLESKY

#include "ngchol/SparseMatrix_impl.hpp"

#endif // _SPARSE_MATRIX_DECL_HPP_
