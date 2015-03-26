#ifndef _csr_weighted_jacobi_h
#define _csr_weighted_jacobi_h
/** 
 * @file csr_weighted_jacobi.h
 * @author Mark Hoemmen
 * @since 23 Jun 2006
 * @date Time-stamp: <2008-07-16 10:55:13 mhoemmen>
 *
 * Implementation of weighted Jacobi smoothing with a CSR format sparse matrix.
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
#include <bebop/util/config.h>
#include <bebop/util/complex.h>

/**
 * val,ind,ptr is a CSR format sparse matrix with n rows.
 */
void
csr_weighted_jacobi_kernel_real_double (double y[], 
					double w, 
					const double val[], 
					const int ind[], 
					const int ptr[], 
					const int n, 
					const double x[], 
					const double b[]);


/**
 * val,ind,ptr is a CSR format sparse matrix with n rows.
 */
void
csr_weighted_jacobi_kernel_complex_double (double_Complex y[], 
					   const double_Complex w, 
					   const double_Complex val[], 
					   const int  ind[], 
					   const int  ptr[], 
					   const int n, 
					   const double_Complex x[],
					   const double_Complex b[]);

#endif /* _csr_weighted_jacobi_h */
