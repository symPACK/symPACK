/**
 * @file csc_matrix.h
 * @author Mark Hoemmen
 * @since 07/01/04 11:08:15 PDT
 * @date Time-stamp: <2008-07-16 10:49:58 mhoemmen>
 * 
 * CSC (compressed sparse column, a.k.a. Harwell-Boeing) format sparse matrix 
 * member functions.
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
#ifndef _csc_matrix_h
#define _csc_matrix_h

#include <bebop/util/enumerations.h>
#include <stdio.h> /* FILE */


/** Evaluates to true iff A is a square matrix */
#define SQUARE_P( A )   (A->m == A->n)


/**
 * @struct csc_matrix_t
 * @author Mark Hoemmen
 * @since 06/08/04 11:11:48 PDT
 *
 * Sparse matrix in Compressed Sparse Column (also known as
 * Harwell-Boeing) format.
 */
struct
csc_matrix_t
{
  /** Number of rows in the matrix */
  int m;
  
  /** Number of columns in the matrix */
  int n;
  
  /** 
   * Number of stored (nonzero) entries.  If the matrix is stored in a 
   * symmetric (or skew, etc.) format, nnz only refers to the number of 
   * stored entries, not the actual number of nonzeros. 
   */
  int nnz;

  /** Array of stored (nonzero) entries of the matrix */
  void* values;

  /** Array of row indices of the stored (nonzero) entries of the matrix */
  int* rowidx;

  /** Array of indices into the rowidx and values arrays, for each column */
  int* colptr;

  /**
   * Symmetry type of the matrix.
   */
  enum symmetry_type_t symmetry_type;

  /**
   * If the matrix has a kind of symmetry (or skew-symmetry): Where the actual 
   * elements of the matrix are stored: in the lower triangle or the upper 
   * triangle.
   */
  enum symmetric_storage_location_t symmetric_storage_location;

  /** 
   * Indicates the type of the entries in val:  REAL means "double", COMPLEX
   * means "double _Complex", PATTERN means the "val" array is NULL and contains
   * no entries.
   */ 
  enum value_type_t value_type;


  /**
   * The ownership mode: tells whether this library or the user is
   * responsible for deallocating input arrays.
   */
  enum ownership_mode_t ownership;

  /**
   * The deallocation function to be called on the values, colidx and
   * rowptr arrays, if ownership == LIBRARY_DEALLOCATES.
   */
  void (*deallocator) (void*);
};


struct sparse_vector_t; /* forward declaration */


/**
 * Fills in the data structure with shallow copies of the given arguments.
 *
 * @param A [OUT]
 * @param m [IN]
 * @param n [IN]
 * @param nnz [IN]
 * @param values [IN]
 * @param rowidx [IN]
 * @param colptr [IN]
 * @param symmetry_type [IN]
 * @param symmetric_storage_location [IN]
 * @param value_type [IN]
 */
void
pack_csc_matrix (struct csc_matrix_t* A, 
		 const int m, const int n, const int nnz, 
		 void* values, int* rowidx, int* colptr,
		 enum symmetry_type_t symmetry_type,
		 enum symmetric_storage_location_t symmetric_storage_location,
		 enum value_type_t value_type,
		 enum ownership_mode_t ownership,
		 void (*deallocator) (void*),
		 enum copy_mode_t copy_mode);

/**
 * Fills in the data structure with the given arguments.
 *
 * @param A [OUT]
 * @param m [IN]
 * @param n [IN]
 * @param nnz [IN]
 * @param values [IN]
 * @param rowidx [IN]
 * @param colptr [IN]
 * @param symmetry_type [IN]
 * @param symmetric_storage_location [IN]
 * @param value_type [IN]
 */
void
init_csc_matrix (struct csc_matrix_t* A, 
		 const int m, const int n, const int nnz, 
		 void* values, int* rowidx, int* colptr,
		 enum symmetry_type_t symmetry_type,
		 enum symmetric_storage_location_t symmetric_storage_location,
		 enum value_type_t value_type,
		 enum ownership_mode_t ownership,
		 void (*deallocator) (void*),
		 enum copy_mode_t copy_mode);

/**
 * Dynamically allocates a new csc_matrix_t object, and fills it in with
 * the given arguments.
 *
 * @param m [IN]
 * @param n [IN]
 * @param nnz [IN]
 * @param values [IN]
 * @param rowidx [IN]
 * @param colptr [IN]
 * @param symmetry_type [IN]
 * @param symmetric_storage_location [IN]
 * @param value_type [IN]
 *
 * @return Pointer to freshly allocated csc_matrix_t object
 */
struct csc_matrix_t*
create_csc_matrix (const int m, const int n, const int nnz, 
		   void* values, int* rowidx, int* colptr,
		   enum symmetry_type_t symmetry_type,
		   enum symmetric_storage_location_t symmetric_storage_location,
		   enum value_type_t value_type,
		   enum ownership_mode_t ownership,
		   void (*deallocator) (void*),
		   enum copy_mode_t copy_mode);

/**
 * Unpacks the CSC format sparse matrix into the given data structures.
 *
 * @param A       [IN]
 * @param m      [OUT]
 * @param n      [OUT]
 * @param nnz    [OUT]
 * @param values [OUT]
 * @param rowidx [OUT]
 * @param colptr [OUT]
 * @param symmetry_type [OUT]
 * @param symmetric_storage_location [OUT]
 * @param value_type [OUT]
 *
 * @warn FIXME: should include the deallocator and the ownership mode!
 */
void
unpack_csc_matrix (const struct csc_matrix_t* A,
		   int* m, int* n, int* nnz,
		   void** values, int** rowidx, int** colptr,
		   enum symmetry_type_t* symmetry_type,
		   enum symmetric_storage_location_t* symmetric_storage_location,
		   enum value_type_t* value_type);

/**
 * Deep-copies src into dest.
 */
void
copy_csc_matrix (struct csc_matrix_t* dest, const struct csc_matrix_t* src);


/**
 * Deallocates the data structures used by A.
 */
void
dealloc_csc_matrix (struct csc_matrix_t* A);


/**
 * Deallocates the data structures used by A, and the csc_matrix_t data
 * structure itself.
 */
void
destroy_csc_matrix (struct csc_matrix_t* A);


/**
 * Returns 0 if the CSC-format matrix is invalid, else returns nonzero.
 *
 * @param A [IN]
 */
int
valid_csc_matrix_p (const struct csc_matrix_t* A);


/**
 * Returns nonzero if the matrices have the same dimensions and nonzero 
 * pattern, else returns zero.
 *
 * @param A [IN]
 * @param B [IN]
 */
int
same_structure_csc_matrix (const struct csc_matrix_t* A, 
			   const struct csc_matrix_t* B);



/**
 * Prints the given CSC-format matrix to the given output stream in MatrixMarket format.
 *
 * @param out [OUT]  Valid (open) output stream
 * @param A [IN]     CSC-format sparse matrix
 *
 * @return Zero if no error, else nonzero.
 */
int
print_csc_matrix_in_matrix_market_format (FILE* out, const struct csc_matrix_t* A);


/**
 * An ASCII-based version of ``spy'' for CSC-format sparse matrices; prints
 * a human-readable visual representation of the nonzero structure of A to 
 * the given output stream.
 *
 * @param out [OUT]  Valid (open) output stream
 * @param A [IN]     CSC-format sparse matrix
 */
void
spy_csc_matrix (FILE* out, const struct csc_matrix_t* A);


/**
 * Writes the given matrix in Harwell-Boeing format to the given file.
 *
 * @param filename [IN] The file to write
 * @param A        [IN] The matrix (in Harwell-Boeing format)
 *
 * @return Nonzero if something went wrong, else zero.
 */
int
save_csc_matrix_in_harwell_boeing_format (const char* const filename, 
					  const struct csc_matrix_t* A);

/**
 * Reads the given matrix in Harwell-Boeing format from the given file.
 *
 * @param filename [IN]  The file to write
 * @param A        [OUT] The matrix (in Harwell-Boeing format)
 *
 * @return Nonzero if something went wrong, else zero.
 */
int
read_harwell_boeing_mat_double (const char* filename, struct csc_matrix_t* A);


/**
 * Generates an m x n CSC-format sparse matrix with diagonal `diag'.
 *
 * @param A [OUT]    CSC-format sparse matrix, uninitialized
 * @param diag [IN]  Diagonal to copy into the matrix: length min(m,n)
 * @param m [IN]
 * @param n [IN]
 *
 * @return Zero if successful, nonzero if an error occurred 
 *         (in which case A is invalid).
 */
int
diag_csc_matrix (struct csc_matrix_t* A, 
		 const double* diag, const int m, const int n);

/**
 * Saves the given CSC-format sparse matrix A to the given filename, 
 * in MatrixMarket format.
 */
int 
save_csc_matrix_in_matrix_market_format (const char* const filename, 
					 const struct csc_matrix_t* A);

/**
 * Prints the given CSC-format matrix to the given output stream in Matlab format.
 *
 * @param out [OUT]  Valid (open) output stream
 * @param A [IN]     CSC-format sparse matrix
 *
 * @return Zero if no error, else nonzero.
 */
int
print_csc_matrix_in_matlab_format (FILE* out, struct csc_matrix_t* A);


/**
 * Saves the given CSC-format sparse matrix A to the given filename, 
 * in Matlab format.
 */
int 
save_csc_matrix_in_matlab_format (const char* const filename, 
				  struct csc_matrix_t* A);



/**
 * If the sparse matrix A is stored in a symmetric space-saving format,
 * expands out A so that all the nonzeros ((i,j) and (j,i)) are explicitly
 * stored.
 *
 * @warn Changes ownership to LIBRARY_DEALLOCATES and the deallocator
 * to the default.
 *
 * @param A [IN/OUT]  Sparse matrix in CSC format
 *
 * @return Nonzero if error, else zero.
 */
int
csc_matrix_expand_symmetric_storage (struct csc_matrix_t* A);

#endif /* _csc_matrix_h */
