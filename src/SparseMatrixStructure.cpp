#include "SparseMatrixStructure.hpp"

namespace LIBCHOLESKY{
//  void SparseMatrixStructure::ClearExpandedSymmetric(){
//    expColptr.Resize(0);
//    expRowind.Resize(0);
//    expanded=false;
//  }
//
//  void SparseMatrixStructure::ExpandSymmetric(){
//    if(!expanded){
//      //code from sparsematrixconverter
//      /* set-up */
//
//      IntNumVec cur_col_nnz(size);
//      IntNumVec new_col_nnz(size);
//
//      /*
//       * Scan A and count how many new non-zeros we'll need to create.
//       *
//       * Post:
//       *   cur_col_nnz[i] == # of non-zeros in col i of the original symmetric 
//       *                     matrix.
//       *   new_col_nnz[i] == # of non-zeros to be stored in col i of the final 
//       *                     expanded matrix.
//       *   new_nnz == total # of non-zeros to be stored in the final expanded 
//       *              matrix.
//       */
//      Int new_nnz = 0; 
//      for (Int i = 0; i < size; i++) 
//      {    
//        cur_col_nnz[i] = colptr[i+1] - colptr[i];
//        new_col_nnz[i] = cur_col_nnz[i];
//        new_nnz += new_col_nnz[i];
//
//      }    
//
//
//
//
//
//      for (Int i = 0; i < size; i++) 
//      {    
//        Int k;
//        for (k = colptr[i]; k < colptr[i+1]; k++) 
//        {
//          Int j = rowind[k-1]-1;
//          if (j != i)
//          {
//            new_col_nnz[j]++;
//            new_nnz++;
//          }
//        }
//      }
//
//      /*
//       *  Initialize row pointers in expanded matrix.
//       *
//       *  Post:
//       *    new_colptr initialized to the correct, final values.
//       *    new_col_nnz[i] reset to be equal to cur_col_nnz[i].
//       */
//      expColptr.Resize(size+1);
//      expColptr[0] = 1;
//      for (Int i = 1; i <= size; i++)
//      {
//        expColptr[i] = expColptr[i-1] + new_col_nnz[i-1];
//        new_col_nnz[i-1] = cur_col_nnz[i-1];
//      }
//      expColptr[size] = new_nnz+1;
//
//      expRowind.Resize(new_nnz);
//
//      /*
//       *  Complete expansion of A to full storage.
//       *
//       *  Post:
//       *    (new_colptr, new_rowind, new_values) is the full-storage equivalent of A.
//       *    new_col_nnz[i] == # of non-zeros in col i of A.
//       */
//
//
//
//      for (Int i = 0; i < size; i++)
//      {
//        Int cur_nnz = cur_col_nnz[i];
//        Int k_cur   = colptr[i] -1;
//        Int k_new   = expColptr[i] -1;
//
//        /* copy current non-zeros from old matrix to new matrix */
//        std::copy(rowind.Data() + k_cur, rowind.Data() + k_cur + cur_nnz , expRowind.Data() + k_new);
//
//        /* fill in the symmetric "missing" values */
//        while (k_cur < colptr[i+1]-1)
//        {
//          /* non-zero of original matrix */
//          Int j = rowind[k_cur]-1;
//
//          if (j != i)  /* if not a non-diagonal element ... */
//          {
//            /* position of this transposed element in new matrix */
//            k_new = expColptr[j]-1 + new_col_nnz[j];
//
//            /* store */
//            expRowind[k_new] = i+1;
//            /*  update so next element stored at row j will appear
//             *  at the right place.  */
//            new_col_nnz[j]++;
//          }
//
//          k_cur++;
//        }
//      }
//      expanded =true;
//    }
//
//  }



}

