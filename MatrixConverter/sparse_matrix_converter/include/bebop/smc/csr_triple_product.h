#ifndef _csr_triple_product_h
#define _csr_triple_product_h
/**
 * @file csr_triple_product.h
 * @author Mark Hoemmen
 * @since 20 Jun 2006
 * @date Time-stamp: <2008-07-16 10:53:18 mhoemmen>
 * 
 * Kernels used for the CSR sparse matrix triple product operation
 * which is useful in algebraic multigrid.
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
 * Returns in val, ind and ptr the CSR sparse matrix triple product
 * R^T * A * P, with double_Complex values.
 */
int
csr_matrix_triple_product_kernel_double_complex (const int m, const int n,
						 const double_Complex RT_val[], 
						 const int RT_ind[],
						 const int RT_ptr[],
						 const double_Complex A_val[],
						 const int A_ind[],
						 const int A_ptr[],
						 const double_Complex P_val[],
						 const int P_ind[],
						 const int P_ptr[],
						 double_Complex** val,
						 int** ind,
						 int** ptr);

/**
 * Returns in val, ind and ptr the CSR sparse matrix triple product
 * R^T * A * P, for pattern matrices (no values).
 */
int
csr_matrix_triple_product_kernel_double_pattern (const int m, const int n,
						 const int RT_ind[],
						 const int RT_ptr[],
						 const int A_ind[],
						 const int A_ptr[],
						 const int P_ind[],
						 const int P_ptr[],
						 int** ind,
						 int** ptr);

/**
 * Performs the "triple product" kernel used for algebraic multigrid:
 * Returns R^T   *   A    *   P, in which R^T, A and P are CSR format 
 *        m x n    n x n    n x m
 * sparse matrices.  Does not exploit possible symmetry in A.  R^T and
 * P may be the same matrix.  A temporary copy of R = R^T^T is 
 * constructed and used internally, and freed after use.
 */
int
csr_matrix_triple_product_kernel_double_real (const int m, const int n,
					      const double RT_val[], 
					      const int RT_ind[],
					      const int RT_ptr[],
					      const double A_val[],
					      const int A_ind[],
					      const int A_ptr[],
					      const double P_val[],
					      const int P_ind[],
					      const int P_ptr[],
					      double** val,
					      int** ind,
					      int** ptr);


/**
 * Computes P' * A * P.
 * Extra memory usage is at worst O(n_coarse * nnz-per-row(P)^2).
 */
int
ptap_csr_dr (const int n_fine, 
	     const int n_coarse,
	     const double A_val[],
	     const int A_ind[],
	     const int A_ptr[],
	     const double P_val[],
	     const int P_ind[],
	     const int P_ptr[],
	     double** val,
	     int** ind,
	     int** ptr,
	     int* nnz);


#endif /* _csr_triple_product_h */
