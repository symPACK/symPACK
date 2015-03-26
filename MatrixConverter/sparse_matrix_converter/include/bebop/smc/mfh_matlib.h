#ifndef _mfh_matlib_h
#define _mfh_matlib_h
/**
 * @file mfh_matlib.h
 * @author Mark Hoemmmen
 * @date Time-stamp: <2008-07-16 10:56:38 mhoemmen>
 * @since 06/08/04 11:05:44 PDT
 * @version 1.2
 * 
 * Some basic sparse and dense matrix manipulation functions.
 *
 * Version 1.1 (mfh 1 July 2004): Moved csc_matrix_t struct definition and 
 * member functions into their own files.
 * Version 1.2 (mfh 26 May 2005): Moved coo_matrix_t and coord_elem_t into
 * their own files.
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

#include <stdio.h> /* FILE */


/**
 * Prints the given dense vector x[0:length-1] to the given output stream.
 */
void
print_dense_vector (FILE* out, const double* x, const int length);


/**
 * Prints the given dense column-oriented m x n matrix to the given output 
 * stream.
 *
 * @param out [OUT]   Output stream
 * @param A   [IN]    Dense column-oriented m x n matrix
 * @param m   [IN]    Number of rows in A
 * @param n   [IN]    Number of columns in A
 * @param trans [IN]  If "T", prints the transpose of the input matrix A.
 *                    The order of the m and n parameters should not be 
 *                    changed.  This option can be used to print a row-
 *                    oriented matrix.
 */
void
print_dense_matrix (FILE* out, const double* A, const int m, const int n, 
					const char* trans);


/** 
 * Equivalent to Matlab's ``A = ones (m,n)''.  Treats A as a flat (1-D) 
 * array; it may be either column-oriented or row-oriented.
 */
void
ones (double* A, const int m, const int n);


/** 
 * Equivalent to Matlab's ``A = zeros (m,n)''.  Treats A as a flat (1-D) 
 * array; it may be either column-oriented or row-oriented.
 */
void
zeros (double* A, const int m, const int n);


/**
 * Fills v[0:n-1] with the which-th unit vector.
 */
void
unit_vector (double* v, const int n, const int which);


/**
 * Returns the transpose of the given column-oriented dense matrix.
 */
double* 
transpose_dense_col_oriented_matrix (const double* A, 
				     const int m, const int n);


#endif /* _mfh_matlib_h */
