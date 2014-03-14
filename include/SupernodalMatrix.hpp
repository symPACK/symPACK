#ifndef _SUPERNODAL_MATRIX_DECL_HPP_
#define _SUPERNODAL_MATRIX_DECL_HPP_

#include "Environment.hpp"
#include "SuperNode.hpp"

#include "NumVec.hpp"
#include "DistSparseMatrix.hpp"
#include "ETree.hpp"
#include "FBMatrix.hpp"



#include <vector>
//#include "utility.hpp"



namespace LIBCHOLESKY{

template <typename T> class SupernodalMatrix{
  protected:
  bool globalAllocated = false;
  IntNumVec Xsuper_;
 
  SparseMatrixStructure Local_;
  SparseMatrixStructure Global_;
  ETree ETree_;
  Int iSize_;

  public:
  std::vector<SuperNode > LocalSupernodes_;

	/// @brief MPI communicator
	//MPI_Comm     comm = MPI_COMM_NULL;        
  //SparseMatrixStructure Local_;
  //SparseMatrixStructure Global_;


  SupernodalMatrix();
//  SupernodalMatrix(const DistSparseMatrix<T> & pMat);
  SupernodalMatrix(const DistSparseMatrix<T> & pMat, FBMatrix* Afactptr );
};



} // namespace LIBCHOLESKY

#include "SupernodalMatrix_impl.hpp"

#endif 
