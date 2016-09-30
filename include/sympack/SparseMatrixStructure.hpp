#ifndef _DIST_SPARSE_MATRIX_STRUCTURE_HPP_
#define _DIST_SPARSE_MATRIX_STRUCTURE_HPP_

#include "sympack/Environment.hpp"
#include "sympack/Ordering.hpp"
#include "sympack/Types.hpp"

namespace symPACK{
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
	std::vector<Ptr>  colptr;                          // Column index pointer
	std::vector<Idx>  rowind;                          // Starting row index pointer

	std::vector<Ptr>  expColptr;                          // Column index pointer expanded
	std::vector<Idx>  expRowind;                          // Starting row index pointer expanded

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
