#ifndef _csr_transpose_h
#define _csr_transpose_h
/**
 * @file csr_transpose.h
 * @author Mark Hoemmen
 * @since 20 Jun 2006
 * @date Time-stamp: <2008-07-16 10:52:58 mhoemmen>
 * 
 * Implements matrix transpose kernels.
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
#include <bebop/util/complex.h>

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
					 const int Aptr[]);
			    
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
					    const int Aptr[]);

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
				     const int Aind[], const int Aptr[]);

#endif /* _csr_transpose_h */
