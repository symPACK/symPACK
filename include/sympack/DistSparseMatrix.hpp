/// @file sparse_matrix_decl.hpp
/// @brief Sparse matrix and Distributed sparse matrix in compressed
/// column format.
/// @author Lin Lin
/// @date 2012-11-10
#ifndef _DIST_SPARSE_MATRIX_DECL_HPP_
#define _DIST_SPARSE_MATRIX_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/ETree.hpp"
#include "sympack/SparseMatrixStructure.hpp"

#include <mpi.h>


namespace symPACK{


class DistSparseMatrixGraph;
template <typename T> class symPACKMatrix;



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
  friend class symPACKMatrix<F>;
  friend class symPACKMatrix<F>;
  //friend functions
  friend void ReadDistSparseMatrixFormatted ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm );
  friend void ReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm );
  friend void ParaWriteDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm );
  friend void ParaReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm );

  friend void ReadDistSparseMatrixFormatted ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm );
  friend void ReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm );
  friend void ParaWriteDistSparseMatrix ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm );
  friend void ParaReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm );

  template <typename SCALAR, typename INSCALAR >
  friend int ReadHB_PARA(std::string & filename, DistSparseMatrix<SCALAR> & HMat);

  template <typename SCALAR, typename INSCALAR >
  friend int ReadHB_PARA_MPIIO(std::string & filename, DistSparseMatrix<SCALAR> & HMat);

  protected:
  bool globalAllocated;
  SparseMatrixStructure Local_;
  SparseMatrixStructure Global_;

  //SparseMatrixGraph Globalg_;

  public:
  DistSparseMatrixGraph Localg_;
	/// @brief Matrix dimension.
	Int          size;

	/// @brief Total number of nonzeros elements.
	//Idx64          nnz;                             
	Ptr          nnz;                             

	/// @brief Dimension nnzLocal, storing the nonzero values.
	std::vector<F>    nzvalLocal;                      

	/// @brief MPI communicator
	MPI_Comm     comm;        

  /// @brief Permutation currently used in the matrix. Empty vector means unpermuted matrix
  std::vector<Int> cinvp;

  inline bool isPermuted(Int * invp, Idx * vertexDist = NULL);

  DistSparseMatrix();
  DistSparseMatrix(MPI_Comm aComm);
  void CopyData(const int n, const int nnz, const int * colptr, const int * rowidx, const F * nzval,bool onebased=false);
  DistSparseMatrix(const int n, const int nnz, const int * colptr, const int * rowidx, const F * nzval , MPI_Comm oComm);
  
  template <typename T> void ConvertData(const int n, const int nnz, const int * colptr, const int * rowidx, const T * nzval,bool onebased=false);

  SparseMatrixStructure  GetGlobalStructure();
  SparseMatrixStructure  GetLocalStructure() const;
  const DistSparseMatrixGraph & GetLocalGraph() const;
  void SetLocalGraph(const DistSparseMatrixGraph & pgraph);
  //const SparseMatrixStructure & GetLocalStructure() const;

  void Dump() const;
  void DumpMatlab() const;

  void Permute(Int * invp, Idx * newVertexDist = NULL);
  void ExpandSymmetric();
  void ToLowerTriangular();
  void SortGraph();
};

// Commonly used
typedef DistSparseMatrix<Real>       DblDistSparseMatrix;
typedef DistSparseMatrix<Complex>    CpxDistSparseMatrix;



} // namespace SYMPACK

#include "sympack/DistSparseMatrix_impl.hpp"

#endif // _DIST_SPARSE_MATRIX_DECL_HPP_
