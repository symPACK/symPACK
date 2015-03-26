#ifndef _sparse_matrix_ops_h
#define _sparse_matrix_ops_h
/**
 * @file sparse_matrix_ops.h
 * @author Mark Hoemmen
 * @since 26 May 2005
 * @date Time-stamp: <2008-07-16 10:57:15 mhoemmen>
 *
 * Operations on general sparse matrices.
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
#include <bebop/smc/sparse_matrix.h>
#include <stdio.h>

/**
 * Prints the given sparse matrix to out in Matrix Market format.
 * 
 * FIXME: should allow any valid output format!
 */
int
print_sparse_matrix (FILE* out, struct sparse_matrix_t* A);

/**
 * Saves the sparse matrix A to the given file, using the given file format.
 *
 * @param filename [IN]  Name of the file to which to save the matrix
 * @param A [IN]         The matrix to save
 * @param file_format [IN]  Output file format
 */
int
save_sparse_matrix (const char* const filename, 
		    struct sparse_matrix_t* A, 
		    enum sparse_matrix_file_format_t file_format);

/**
 * Loads a sparse matrix from the given file, with the assumption that 
 * the file is in the format file_format.
 *
 * @param file_format [IN]      Input file format 
 * @param matrix_filename [IN]  Name of the file from which to load the matrix
 *
 * @return If no error, the sparse matrix; else NULL.
 */
struct sparse_matrix_t*
load_sparse_matrix (enum sparse_matrix_file_format_t file_format, 
		    const char *const matrix_filename);

/**
 * If A is stored in a space-saving symmetric format, expands A in place so 
 * that all the nonzero elements are stored explicitly.
 *
 * @param A [IN/OUT]  Pointer to sparse matrix object
 *
 * @return Zero if successful, else nonzero.
 */
int
sparse_matrix_expand_symmetric_storage (struct sparse_matrix_t* A);

/**
 * Converts the sparse matrix A from its current format into the new format
 * output_format.
 *
 * @param A [IN/OUT]   Sparse matrix to convert
 * @param output_format [IN]  New storage format
 *
 * @return Zero if successful, else nonzero.
 */
int
sparse_matrix_convert (struct sparse_matrix_t* A, 
		       enum sparse_matrix_storage_format_t output_format);

/**
 * Returns 0 if the given sparse matrix has a valid representation, else 
 * returns nonzero.
 */
int
valid_sparse_matrix (struct sparse_matrix_t* A);

/**
 * Given a string, returns the corresponding sparse matrix 
 * (internal) storage format enum. 
 */
enum sparse_matrix_storage_format_t 
sparse_matrix_storage_format_string_to_enum (const char* const s);

/**
 * Given a string, returns the corresponding sparse matrix 
 * file format enum. 
 */
enum sparse_matrix_file_format_t
sparse_matrix_file_format_string_to_enum (const char* const s);

#endif /* _sparse_matrix_ops_h */
