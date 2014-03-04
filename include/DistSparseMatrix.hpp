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

extern "C" {
#include <bebop/util/config.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/csc_matrix.h>
#include <bebop/smc/sparse_matrix_ops.h>

#include <bebop/util/get_options.h>
#include <bebop/util/init.h>
#include <bebop/util/malloc.h>
#include <bebop/util/timer.h>
#include <bebop/util/util.h>
}





namespace LIBCHOLESKY{



class SparseMatrixStructure{
  public:
	Int          size;                            // Matrix dimension
	Int          nnz;                             // Number of nonzeros
	IntNumVec    colptr;                          // Column index pointer
	IntNumVec    rowind;                          // Starting row index pointer

	IntNumVec    expColptr;                          // Column index pointer expanded
	IntNumVec    expRowind;                          // Starting row index pointer expanded
  bool expanded = false;

  inline void ClearExpandedSymmetric(){
    expColptr.Resize(0);
    expRowind.Resize(0);
    expanded=false;
  }

  inline void ExpandSymmetric(){
    if(!expanded){
      //code from sparsematrixconverter
      /* set-up */

      IntNumVec cur_col_nnz(size);
      IntNumVec new_col_nnz(size);

      /*
       * Scan A and count how many new non-zeros we'll need to create.
       *
       * Post:
       *   cur_col_nnz[i] == # of non-zeros in col i of the original symmetric 
       *                     matrix.
       *   new_col_nnz[i] == # of non-zeros to be stored in col i of the final 
       *                     expanded matrix.
       *   new_nnz == total # of non-zeros to be stored in the final expanded 
       *              matrix.
       */
      Int new_nnz = 0; 
      for (Int i = 0; i < size; i++) 
      {    
        cur_col_nnz[i] = colptr[i+1] - colptr[i];
        new_col_nnz[i] = cur_col_nnz[i];
        new_nnz += new_col_nnz[i];

      }    





      for (Int i = 0; i < size; i++) 
      {    
        Int k;
        for (k = colptr[i]; k < colptr[i+1]; k++) 
        {
          Int j = rowind[k-1]-1;
          if (j != i)
          {
            new_col_nnz[j]++;
            new_nnz++;
          }
        }
      }

      /*
       *  Initialize row pointers in expanded matrix.
       *
       *  Post:
       *    new_colptr initialized to the correct, final values.
       *    new_col_nnz[i] reset to be equal to cur_col_nnz[i].
       */
      expColptr.Resize(size+1);
      expColptr[0] = 1;
      for (Int i = 1; i <= size; i++)
      {
        expColptr[i] = expColptr[i-1] + new_col_nnz[i-1];
        new_col_nnz[i-1] = cur_col_nnz[i-1];
      }
      expColptr[size] = new_nnz+1;

      expRowind.Resize(new_nnz);

      /*
       *  Complete expansion of A to full storage.
       *
       *  Post:
       *    (new_colptr, new_rowind, new_values) is the full-storage equivalent of A.
       *    new_col_nnz[i] == # of non-zeros in col i of A.
       */



      for (Int i = 0; i < size; i++)
      {
        Int cur_nnz = cur_col_nnz[i];
        Int k_cur   = colptr[i] -1;
        Int k_new   = expColptr[i] -1;

        /* copy current non-zeros from old matrix to new matrix */
        std::copy(rowind.Data() + k_cur, rowind.Data() + k_cur + cur_nnz , expRowind.Data() + k_new);

        /* fill in the symmetric "missing" values */
        while (k_cur < colptr[i+1]-1)
        {
          /* non-zero of original matrix */
          Int j = rowind[k_cur]-1;

          if (j != i)  /* if not a non-diagonal element ... */
          {
            /* position of this transposed element in new matrix */
            k_new = expColptr[j]-1 + new_col_nnz[j];

            /* store */
            expRowind[k_new] = i+1;
            /*  update so next element stored at row j will appear
             *  at the right place.  */
            new_col_nnz[j]++;
          }

          k_cur++;
        }
      }
      expanded =true;
    }

  }



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
	MPI_Comm     comm = MPI_COMM_NULL;        

  SparseMatrixStructure Local_;
  SparseMatrixStructure Global_;

  void CopyData(const csc_matrix_t * cscptr);
  DistSparseMatrix(const csc_matrix_t * cscptr);
  DistSparseMatrix(MPI_Comm oComm, const csc_matrix_t * cscptr);

  void ToGlobalStruct();

  void ConstructETree(ETree & tree);
  void ConstructETreeBis(ETree & tree);
  void GetLColRowCount(ETree & tree, IntNumVec & cc, IntNumVec & rc);
  void FindSupernodes(ETree& tree, IntNumVec & cc, IntNumVec & xsuper);
  void SymbolicFactorization(ETree& tree,const IntNumVec & cc,const IntNumVec & xsuper, IntNumVec & xlindx, IntNumVec & xlnz,  IntNumVec & lindx, DblNumVec & lnz);
};

// Commonly used
typedef DistSparseMatrix<Real>       DblDistSparseMatrix;
typedef DistSparseMatrix<Complex>    CpxDistSparseMatrix;



} // namespace LIBCHOLESKY

#include "DistSparseMatrix_impl.hpp"

#endif // _DIST_SPARSE_MATRIX_DECL_HPP_
