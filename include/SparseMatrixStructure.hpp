#ifndef _DIST_SPARSE_MATRIX_STRUCTURE_HPP_
#define _DIST_SPARSE_MATRIX_STRUCTURE_HPP_

#include "Environment.hpp"
#include "NumVec.hpp"

namespace LIBCHOLESKY{
class SparseMatrixStructure{
  protected:
  bool bIsGlobal=false;
  bool bIsExpanded = false;
  public:
	Int          size;                            // Matrix dimension
	Int          nnz;                             // Number of nonzeros
	IntNumVec    colptr;                          // Column index pointer
	IntNumVec    rowind;                          // Starting row index pointer

	IntNumVec    expColptr;                          // Column index pointer expanded
	IntNumVec    expRowind;                          // Starting row index pointer expanded

  void ClearExpandedSymmetric();

  void ExpandSymmetric();

  void ToGlobal(SparseMatrixStructure & pGlobal);


};




}

#endif
