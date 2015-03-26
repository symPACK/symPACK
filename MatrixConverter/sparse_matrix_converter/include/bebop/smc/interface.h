#ifndef _interface_h
#define _interface_h
/**
 * @file interface.h
 * @author Mark Hoemmen
 * @since 09 Jun 2006
 * @date Time-stamp: <2008-07-16 10:55:32 mhoemmen>
 * 
 * Basic sparse matrix functionality.  Intended for use by a 
 * Foreign Function Interface, SWIG, or some other such method 
 * for calling C functions from other languages.  This is really
 * only intended for basic "load a file, convert it and save it"
 * functionality.  For more extensive functionality, the files
 * sparse_matrix.h and sparse_matrix_ops.h are recommended.
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

/* 
 * Forward declaration of what is here just an opaque type. 
 * See sparse_matrix.h for a definition.
 */
struct sparse_matrix_t;

/**
 * Loads a sparse matrix of the given type from the given file.
 *
 * @param path [IN]  Path to the sparse matrix file
 * @param fmt [IN]   File format.  Currently accepted values
 *                   are "MATRIX_MARKET", "HARWELL_BOEING" and
 *                   "MATLAB".
 * @return Pointer to the sparse matrix loaded from the file,
 * or NULL if the load failed (e.g. incorrect path or file format).
 */
struct sparse_matrix_t*
sp_load (const char* path, const char* fmt);

/**
 * Saves the given sparse matrix to the given path in the given format.
 *
 * @param A [IN]  Pointer to a sparse matrix
 * @param path [IN]  Path to the desired output file
 * @param fmt [IN]  File format.  Currently accepted values
 *                  are "MATRIX_MARKET", "HARWELL_BOEING" and "MATLAB".
 * @return 0 if successful, else nonzero.
 */
int
sp_save (struct sparse_matrix_t* A, const char* path, 
      const char* fmt);

/**
 * Prints the internal storage format of the given sparse matrix
 * to stdout.
 *
 * @param A [IN]  Pointer to a sparse matrix
 */
void
sp_format (struct sparse_matrix_t* A);

/**
 * Converts the internal storage format of the given sparse matrix
 * to the given type.
 *
 * @param A [IN/OUT]  Pointer to a sparse matrix
 * @param type [IN]   The target type, e.g. "COO", "CSR", "JAD"
 *
 * @return 0 if conversion was successful, else nonzero.
 *
 * @warn Not all conversions are supported.
 */
int
sp_convert (struct sparse_matrix_t* A, const char* type);

/**
 * Returns the sparse matrix which is the product B * A.
 *
 * @param B [IN]  Pointer to a sparse matrix
 * @param A [IN]  Pointer to a sparse matrix
 *
 * @return The product matrix B * A
 *
 * @warn Experimental and not supported for all matrix types!
 */
struct sparse_matrix_t* 
sp_mult (struct sparse_matrix_t* B, struct sparse_matrix_t* A);

/**
 * Returns the sparse matrix R*A*P, in which R is given as R^T.
 *
 * @warn Experimental and not supported for all matrix types!
 */
struct sparse_matrix_t* 
sp_triprod (struct sparse_matrix_t* RT, 
	    struct sparse_matrix_t* A, 
	    struct sparse_matrix_t* P);
/**
 * Prints A to stdout in Matrix Market format.
 */
void
sp_print (struct sparse_matrix_t* A);

/**
 * Deallocates the given sparse matrix handle and its internal storage 
 * of the sparse matrix.
 *
 * @param A [IN/OUT] Pointer to a sparse matrix
 */
void
sp_destroy (struct sparse_matrix_t* A);

#endif /* _interface_h */
