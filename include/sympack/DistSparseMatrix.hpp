#ifndef _DIST_SPARSE_MATRIX_DECL_HPP_
#define _DIST_SPARSE_MATRIX_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/ETree.hpp"

#include <mpi.h>


namespace symPACK{



class DistSparseMatrixGraph;
class symPACKMatrixBase;
template <typename T> class symPACKMatrix;

class DistSparseMatrixBase{
  public:
    virtual ~DistSparseMatrixBase(){};
};


template <typename F> class DistSparseMatrix: public DistSparseMatrixBase{
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

  inline bool isPermuted(Int * invp, Idx * vertexDist = nullptr);

  DistSparseMatrix();
  DistSparseMatrix(MPI_Comm aComm);
  void CopyData(const int n, const int nnz, const int * colptr, const int * rowidx, const F * nzval,bool onebased=false);
  DistSparseMatrix(const int n, const int nnz, const int * colptr, const int * rowidx, const F * nzval , MPI_Comm oComm);
  
  template <typename T> void ConvertData(const int n, const int nnz, const int * colptr, const int * rowidx, const T * nzval,bool onebased=false);

  const DistSparseMatrixGraph & GetLocalGraph() const;
  DistSparseMatrixGraph & GetLocalGraph();

  void SetLocalGraph(const DistSparseMatrixGraph & pgraph);

  void Dump() const;
  void DumpMatlab() const;

  void Permute(Int * invp, Idx * newVertexDist = nullptr);
  void ExpandSymmetric();
  void ToLowerTriangular();
  void SortGraph();
};

// Commonly used
typedef DistSparseMatrix<Real>       DblDistSparseMatrix;
typedef DistSparseMatrix<Complex>    CpxDistSparseMatrix;

template <typename T>
  void symPACK_LoadGraph( T & pHMat, int n, int * colptr , int * rowind);






} // namespace SYMPACK

#include "sympack/impl/DistSparseMatrix_impl.hpp"

#endif // _DIST_SPARSE_MATRIX_DECL_HPP_
