/// @file sparse_matrix_decl.hpp
/// @brief Sparse matrix and Distributed sparse matrix in compressed
/// column format.
/// @author Lin Lin
/// @date 2012-11-10
#ifndef _DIST_SPARSE_MATRIX_DECL_HPP_
#define _DIST_SPARSE_MATRIX_DECL_HPP_

#include <vector>
#include "ngchol/ETree.hpp"
#include "ngchol/SparseMatrixStructure.hpp"

#include <mpi.h>



using namespace std;
namespace LIBCHOLESKY{



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
template <typename F> class DistSparseMatrix{
  protected:
  bool globalAllocated;
  SparseMatrixStructure Local_;
  SparseMatrixStructure Global_;


  public:
	/// @brief Matrix dimension.
	int          size;

	/// @brief Total number of nonzeros elements.
	int          nnz;                             

	/// @brief Dimension nnzLocal, storing the nonzero values.
	vector<F>    nzvalLocal;                      

	/// @brief MPI communicator
	MPI_Comm     comm;        


  DistSparseMatrix(MPI_Comm aComm){globalAllocated=false; comm = aComm;};
  void CopyData(const int n, const int nnz, const int * colptr, const int * rowidx, const F * nzval,bool onebased=false );
  DistSparseMatrix(const int n, const int nnz, const int * colptr, const int * rowidx, const F * nzval , MPI_Comm oComm);
  SparseMatrixStructure  GetGlobalStructure();
  SparseMatrixStructure  GetLocalStructure() const;
  //const SparseMatrixStructure & GetLocalStructure() const;

};




} // namespace LIBCHOLESKY

#include "ngchol/DistSparseMatrix_impl.hpp"

#endif // _DIST_SPARSE_MATRIX_DECL_HPP_
