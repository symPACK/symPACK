/// @file sparse_matrix_decl.hpp
/// @brief Sparse matrix and Distributed sparse matrix in compressed
/// column format.
/// @author Lin Lin
/// @date 2012-11-10
#ifndef _DIST_SPARSE_MATRIX_DECL_HPP_
#define _DIST_SPARSE_MATRIX_DECL_HPP_

#include "Environment.hpp"
#include "NumVec.hpp"
#include "ETree.hpp"
//#include "utility.hpp"

namespace LIBCHOLESKY{



struct SparseMatrixStructure{
	Int          size;                            // Matrix dimension
	Int          nnz;                             // Number of nonzeros
	IntNumVec    colptr;                          // Column index pointer
	IntNumVec    rowind;                          // Starting row index pointer
};




/// @class DistSparseMatrix
///
/// @brief DistSparseMatrix describes a Sparse matrix in the compressed
/// sparse column format (CSC) and distributed with column major partition. 
///
/// Note
/// ----
/// 
/// Since in PEXSI and PPEXSI only symmetric matrix is considered, the
/// compressed sparse row format will also be represented by the
/// compressed sparse column format.
///
/// TODO Add the parameter of numColLocal
template <class F> class DistSparseMatrix{
  protected:
  bool globalAllocated = false;


  public:
	/// @brief Matrix dimension.
	Int          size;         

	/// @brief Total number of nonzeros elements.
	Int          nnz;                             

//	/// @brief Local number of local nonzeros elements on this processor.
//	Int          nnzLocal;                        
//
//	/// @brief Dimension numColLocal + 1, storing the pointers to the
//	/// nonzero row indices and nonzero values in rowptrLocal and
//	/// nzvalLocal, respectively.  numColLocal is the number
//	/// of local columns saved on this processor. The indices are 1-based
//	/// (FORTRAN-convention), i.e.  colptrLocal[0] = 1. 
//	IntNumVec    colptrLocal;                     
//
//	/// @brief Dimension nnzLocal, storing the nonzero indices.
//	/// The indices are 1-based (FORTRAN-convention), i.e. the first row
//	/// index is 1. 
//	IntNumVec    rowindLocal;                    
	




	/// @brief Dimension nnzLocal, storing the nonzero values.
	NumVec<F>    nzvalLocal;                      

	/// @brief MPI communicator
	MPI_Comm     comm;        

  SparseMatrixStructure Local_;
  SparseMatrixStructure Global_;



  void ToGlobalStruct();

  void ConstructETree(ETree & tree);
  void GetLColRowCount(ETree & tree, IntNumVec & cc, IntNumVec & rc);
  void FindSupernodes(ETree& tree, IntNumVec & cc, IntNumVec & xsuper);
};

// Commonly used
typedef DistSparseMatrix<Real>       DblDistSparseMatrix;
typedef DistSparseMatrix<Complex>    CpxDistSparseMatrix;



} // namespace LIBCHOLESKY

#include "DistSparseMatrix_impl.hpp"

#endif // _DIST_SPARSE_MATRIX_DECL_HPP_
