#ifndef _csr_matmatmult_h
#define _csr_matmatmult_h
/**
 * @file csr_matmatmult.h
 * @author Mark Hoemmen <mhoemmen@cs.berkeley.edu>
 * @date Time-stamp: <2008-07-16 10:51:02 mhoemmen>
 * 
 * Copyright (c) 2008, Regents of the University of California All
 * rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the
 * distribution.
 * 
 * * Neither the name of the University of California, Berkeley, nor
 * the names of its contributors may be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
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
 * Computes C := alpha * A * B, in which A, B and C are CSR format
 * sparse matrices storing double precision real values.  A is
 * dimension m x p and B is dimension p x n.
 *
 * @param pCptr [OUT] pointer to the ptr array in the CSR structure of C
 * @param pCind [OUT] pointer to the ind array in the CSR structure of C
 * @param pCval [OUT] pointer to the val array in the CSR structure of C
 * @param pCnnz [OUT] pointer to the number of (structural) nonzeros in C
 * @param alpha [IN]  scalar multiplier
 * @param Aptr [IN]   ptr array in the CSR structure of A
 * @param Aind [IN]   ind array in the CSR structure of A
 * @param Aval [IN]   val array in the CSR structure of A
 * @param Bptr [IN]   ptr array in the CSR structure of B
 * @param Bind [IN]   ind array in the CSR structure of B
 * @param Bval [IN]   val array in the CSR structure of B
 * @param m [IN]      number of rows in A and C 
 * @param p [IN]      number of columns in A, and number of rows in B
 * @param n [IN]      number of columns in B, and number of columns in C
 *
 * @return Error code (zero if no error)
 */
int
csr_matmatmult_double_real (int** pCptr, int** pCind, double** pCval, int* pCnnz,
			    double alpha, int* Aptr, int* Aind, double* Aval, 
			    int* Bptr, int* Bind, double* Bval, 
			    const int m, const int p, const int n);

/**
 * Computes C := alpha * A * B, in which A, B and C are CSR format
 * sparse matrices storing double precision complex values.  A is
 * dimension m x p and B is dimension p x n.
 *
 * @param pCptr [OUT] pointer to the ptr array in the CSR structure of C
 * @param pCind [OUT] pointer to the ind array in the CSR structure of C
 * @param pCval [OUT] pointer to the val array in the CSR structure of C
 * @param pCnnz [OUT] pointer to the number of (structural) nonzeros in C
 * @param alpha [IN]  scalar multiplier
 * @param Aptr [IN]   ptr array in the CSR structure of A
 * @param Aind [IN]   ind array in the CSR structure of A
 * @param Aval [IN]   val array in the CSR structure of A
 * @param Bptr [IN]   ptr array in the CSR structure of B
 * @param Bind [IN]   ind array in the CSR structure of B
 * @param Bval [IN]   val array in the CSR structure of B
 * @param m [IN]      number of rows in A and C 
 * @param p [IN]      number of columns in A, and number of rows in B
 * @param n [IN]      number of columns in B, and number of columns in C
 *
 * @return Error code (zero if no error)
 */
int
csr_matmatmult_double_complex (int** pCptr, int** pCind, 
			       double_Complex** pCval, 
			       int* pCnnz, double_Complex alpha, int* Aptr, 
			       int* Aind, double_Complex* Aval, int* Bptr, 
			       int* Bind, double_Complex* Bval, 
			       const int m, const int p, const int n);

/**
 * Computes C := A * B, in which A, B and C are CSR format
 * sparse matrices storing patterns (no nonzero values, just indices).  
 * A is dimension m x p and B is dimension p x n.
 *
 * @param pCptr [OUT] pointer to the ptr array in the CSR structure of C
 * @param pCind [OUT] pointer to the ind array in the CSR structure of C
 * @param pCnnz [OUT] pointer to the number of (structural) nonzeros in C
 * @param Aptr [IN]   ptr array in the CSR structure of A
 * @param Aind [IN]   ind array in the CSR structure of A
 * @param Bptr [IN]   ptr array in the CSR structure of B
 * @param Bind [IN]   ind array in the CSR structure of B
 * @param m [IN]      number of rows in A and C 
 * @param p [IN]      number of columns in A, and number of rows in B
 * @param n [IN]      number of columns in B, and number of columns in C
 *
 * @return Error code (zero if no error)
 */
int
csr_matmatmult_pattern (int** pCptr, int** pCind, int* pCnnz, int* Aptr, 
			int* Aind, int* Bptr, int* Bind, 
			const int m, const int p, const int n);

#endif /* _csr_matmatmult_h */
