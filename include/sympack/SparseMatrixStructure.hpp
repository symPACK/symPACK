#ifndef _DIST_SPARSE_MATRIX_STRUCTURE_HPP_
#define _DIST_SPARSE_MATRIX_STRUCTURE_HPP_

#include "sympack/Environment.hpp"
#include "sympack/Ordering.hpp"
#include "sympack/Types.hpp"

namespace SYMPACK{
class ETree;
class Ordering;
class DistSparseMatrixGraph;
template <typename F> class DistSparseMatrix;

class SparseMatrixStructure{
  friend class Ordering;
  friend class DistSparseMatrixGraph;
  template <typename F> friend class DistSparseMatrix;


  protected:
  bool bIsGlobal;
  bool bIsExpanded;
  

  public:
  int baseval;
  int keepDiag;
  int sorted;
	Int          size;                            // Matrix dimension
	Int          nnz;                             // Number of nonzeros
	SYMPACK::vector<Ptr>  colptr;                          // Column index pointer
	SYMPACK::vector<Idx>  rowind;                          // Starting row index pointer

	SYMPACK::vector<Ptr>  expColptr;                          // Column index pointer expanded
	SYMPACK::vector<Idx>  expRowind;                          // Starting row index pointer expanded

  SparseMatrixStructure();

  //Distributed graph to global SparseMatrixStructure
  SparseMatrixStructure( const DistSparseMatrixGraph & G );


  bool IsExpanded() const { return bIsExpanded;}


  void ClearExpandedSymmetric();

  void ExpandSymmetric(MPI_Comm * workcomm = NULL);

  void ToGlobal(SparseMatrixStructure & pGlobal,MPI_Comm & comm);


};



}

#endif
