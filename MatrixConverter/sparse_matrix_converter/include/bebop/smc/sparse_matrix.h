#ifndef _sparse_matrix_h
#define _sparse_matrix_h
/**
 * @file sparse_matrix.h
 * @author Mark Hoemmen
 * @since 26 May 2005
 * @date Time-stamp: <2008-07-16 10:57:06 mhoemmen>
 *
 * A wrapper representation of a general sparse matrix.
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

/**
 * Enum indicating the internal storage format of the sparse matrix.
 */
enum
sparse_matrix_storage_format_t
{
  COO, CSC, CSR, BCOO, BCSR, JAD
};

/**
 * The supported file formats for sparse matrices.
 */
enum
sparse_matrix_file_format_t
{
  HARWELL_BOEING, MATRIX_MARKET, MATLAB
};

/**
 * Wrapper struct for a sparse matrix.  
 */
struct
sparse_matrix_t
{
  /** The internal storage format of the sparse matrix. */
  enum sparse_matrix_storage_format_t format;

  /** 
   * Pointer to the internal representation of the sparse matrix.  For example,
   * If format==CSC, this is a "struct csc_matrix_t*"; if format==COO, this is 
   * a "struct coo_matrix_t*". 
   */
  void* repr;
};

/**
 * Allocates a wrapper object with shallow copies of the given format and 
 * representation.
 */
struct sparse_matrix_t*
create_sparse_matrix (enum sparse_matrix_storage_format_t format, 
		      void* repr);

/** 
 * If no error, deallocates the internal representation of the given sparse
 * matrix, as well as the sparse_matrix_t struct itself, and returns zero.  If
 * error, returns nonzero without trying to deallocate anything.
 *
 * @param A [IN/OUT]   Pointer to sparse matrix object
 *
 * @return Zero if no error, else nonzero.
 */
int 
destroy_sparse_matrix (struct sparse_matrix_t* A);

/**
 * Returns a string version of the internal storage format of the given
 * sparse matrix.
 */
const char* 
sparse_matrix_format_string (struct sparse_matrix_t* A);

/**
 * Returns C := B * A, in which C is sparse with the same formats as B and A.
 * Precondition:  B and A are in the same format.
 * 
 * @note If sparse matrix-matrix multiplication is not supported for
 * the given format, returns NULL.
 */
struct sparse_matrix_t*
sparse_matrix_matmatmult (struct sparse_matrix_t* B, struct sparse_matrix_t* A);


#endif /* _sparse_matrix_h */
