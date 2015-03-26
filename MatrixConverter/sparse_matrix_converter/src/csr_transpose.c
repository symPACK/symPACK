/**
 * @file csr_transpose.c
 * @author Mark Hoemmen
 * @since 20 Jun 2006
 * @date Time-stamp: <2008-07-16 11:16:46 mhoemmen>
 * 
 * Implementation of CSR sparse matrix transpose kernels.
 *
 * @note A result of refactoring code and needing a separate transpose routine.
 *
 * Copyright (c) 2008, Regents of the University of California 
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the
 * following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright 
 *   notice, this list of conditions and the following disclaimer in 
 *   the documentation and/or other materials provided with the 
 *   distribution.
 *
 * * Neither the name of the University of California, Berkeley, nor
 *   the names of its contributors may be used to endorse or promote
 *   products derived from this software without specific prior
 *   written permission.  
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include <bebop/smc/csr_transpose.h>

#include <bebop/util/complex.h>
#include <bebop/util/malloc.h>
#include <bebop/util/util.h>

#include <assert.h>

/**
 * Returns A^T, in which A is an m x n CSR format sparse matrix with 
 * double precision real values.  Assumes that A is not stored to 
 * exploit symmetry.
 *
 * @param val [OUT] Pointer to array storing the nonzero values of the
 *                  output matrix
 * @param ind [OUT] Pointer to array storing the column indices of the
 *                  output matrix
 * @param ptr [OUT] Pointer to array storing the row pointers of the
 *                  output matrix
 * @param m [IN]     Number of rows in the input matrix A
 * @param n [IN]     Number of columns in the input matrix A
 * @param Aval [IN]  Array storing the nonzero values of the input matrix
 * @param Aind [IN]  Array storing the column indices of the input matrix
 * @param Aptr [IN]  Array storing the column pointers of the input matrix
 *
 * @return 0 if no errors, else nonzero
 */
int
csr_matrix_transpose_kernel_double_real (double** val, int** ind, int** ptr,
					 const int m, const int n,
					 const double Aval[], const int Aind[],
					 const int Aptr[])
{
  const int nnz = Aptr[m];
  int* col_nnz = NULL; /* # of nonzeros in each column */
  int i;

  double* _val;
  int* _ind;
  int* _ptr;

  /* Allocate space for the transposed matrix */
  _ptr = bebop_malloc ((n+1) * sizeof (int));
  _ind = bebop_malloc (nnz * sizeof (int));
  _val = bebop_malloc (nnz * sizeof (double));


  /* Count # of nonzeros per column */
  col_nnz = bebop_calloc (n, sizeof (int));
  for (i = 0; i < nnz; i++)
    col_nnz[Aind[i]]++;

  /* 
   * Initialize row pointers of transposed matrix.
   * Reset col_nnz to zero, so that we can use it 
   * later to keep track of the number of nonzeros
   * added to each row.
   */
  _ptr[0] = 0;
  for (i = 1; i <= n; i++)
    {
      _ptr[i] = _ptr[i-1] + col_nnz[i-1];
      col_nnz[i-1] = 0;
    }

  for (i = 0; i < m; i++)
    {
      int k;
      int nnz_row; /* # of nonzeros in row i of A */
      
      nnz_row = Aptr[i+1] - Aptr[i];
      for (k = Aptr[i]; k < Aptr[i] + nnz_row; k++)
	{
	  int j = Aind[k];              /* col index */
	  double a = Aval[k];           /* nonzero value */
	  int h = _ptr[j] + col_nnz[j]; /* nonzero position */

	  /* Add the nonzero A(i,j) to the transposed matrix */
	  _ind[h] = i;
	  _val[h] = a;

	  /* Keep track of how many entries we've added thus far */
	  col_nnz[j]++;
	}
    }

  /* Free the temporary storage */
  bebop_free (col_nnz);

  /* Return the matrix and indicate success */
  *val = _val;
  *ind = _ind;
  *ptr = _ptr;
  return 0;
}
				    


/**
 * Returns A^T, in which A is an m x n CSR format sparse matrix with 
 * double precision complex values.  Assumes that A is not stored to 
 * exploit symmetry.
 *
 * @param val [OUT] Pointer to array storing the nonzero values of the
 *                  output matrix
 * @param ind [OUT] Pointer to array storing the column indices of the
 *                  output matrix
 * @param ptr [OUT] Pointer to array storing the row pointers of the
 *                  output matrix
 * @param m [IN]     Number of rows in the input matrix A
 * @param n [IN]     Number of columns in the input matrix A
 * @param Aval [IN]  Array storing the nonzero values of the input matrix
 * @param Aind [IN]  Array storing the column indices of the input matrix
 * @param Aptr [IN]  Array storing the column pointers of the input matrix
 *
 * @return 0 if no errors, else nonzero
 */
int
csr_matrix_transpose_kernel_double_complex (double** val, int** ind, int** ptr,
					    const int m, const int n,
					    const double Aval[], const int Aind[],
					    const int Aptr[])
{
  /* Not yet implemented */
  assert (0); 
  return -1;
}


/**
 * Returns A^T, in which A is an m x n CSR format sparse matrix with 
 * no values, just a nonzero pattern.  Assumes that A is not stored to 
 * exploit symmetry.
 *
 * @param ind [OUT] Pointer to array storing the column indices of the
 *                  output matrix
 * @param ptr [OUT] Pointer to array storing the row pointers of the
 *                  output matrix
 * @param m [IN]     Number of rows in the input matrix A
 * @param n [IN]     Number of columns in the input matrix A
 * @param Aind [IN]  Array storing the column indices of the input matrix
 * @param Aptr [IN]  Array storing the column pointers of the input matrix
 *
 * @return 0 if no errors, else nonzero
 */
int
csr_matrix_transpose_kernel_pattern (int** ind, int** ptr, 
				     const int m, const int n,
				     const int Aind[], const int Aptr[])
{
  /* Not yet implemented */
  assert (0); 
  return -1;
}

