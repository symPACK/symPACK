#ifndef _coo_matrix_h
#define _coo_matrix_h
/**
 * @file coo_matrix.h
 * @author Mark Hoemmen
 * @since 26 May 2005
 * @date Time-stamp: <2008-07-16 10:50:16 mhoemmen>
 *
 * Wrapper for a coordinate (ijv) format sparse matrix.
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
#include <bebop/util/enumerations.h>
#include <stdio.h> /* FILE */


/**
 * A coordinate (COO) format sparse matrix.  This uses the ``struct of
 * arrays'' storage paradigm.  An example of the opposite storage
 * paradigm, ``array of structs,'' would be an array of coord_elem_t
 * (which see).
 */
struct
coo_matrix_t
{
  /** Number of rows */
  int m;
  
  /** Number of columns */
  int n;
  
  /** 
   * Number of stored (nonzero) entries.  If the matrix is stored in a 
   * symmetric (or skew, etc.) format, nnz only refers to the number of 
   * stored entries, not the actual number of nonzeros. 
   */
  int nnz;
 
  /** For entry e, II[e] is its row (we can't call it I because that is reserved as a macro in C99 for sqrt{-1}) */
  int* II;

  /** For entry e, JJ[e] is its column */
  int* JJ;

  /** For entry e, val[e] is its value.  Type of the entries depends 
   * on the value of value_type. */
  void* val;

  /** Base of indexing (either one or zero) */
  enum index_base_t index_base;

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

  enum ownership_mode_t ownership;
  void (*deallocator) (void*);
};


/**
 * Initializes the given struct A to be the given COO format sparse matrix.
 *
 * @param A [OUT]  Valid pointer to uninitialized struct.
 * @param m [IN]   Number of rows
 * @param n [IN]   Number of columns
 * @param nnz [IN]  Number of nonzero entries
 * @param II [IN]    Array of row indices
 * @param JJ [IN]    Array of column indices
 * @param val [IN]  Array of nonzero values
 * @param index_base [IN]   Base of the indices (ZERO or ONE).
 * @param symmetry_type [IN]  Symmetry type of the matrix
 * @param symmetric_storage_location [IN]  If the matrix has symmetry, where 
 * 	the nonzeros are stored.  If the matrix is to be stored in unsymmetric 
 * 	format, the value of this enum is irrelevant.
 * @param value_type [IN]  
 */
void
init_coo_matrix (struct coo_matrix_t *A, int m, int n, int nnz, int *II, 
		 int *JJ, void *val, enum index_base_t index_base, 
		 enum symmetry_type_t symmetry_type, 
		 enum symmetric_storage_location_t symmetric_storage_location,
                 enum value_type_t value_type,
		 enum ownership_mode_t ownership,
		 void (*deallocator) (void*),
		 enum copy_mode_t copy_mode);

struct coo_matrix_t*
create_coo_matrix (int m, int n, int nnz, int *II, 
		   int *JJ, void *val, enum index_base_t index_base, 
		   enum symmetry_type_t symmetry_type, 
		   enum symmetric_storage_location_t symmetric_storage_location,
		   enum value_type_t value_type,
		   enum ownership_mode_t ownership,
		   void (*deallocator) (void*),
		   enum copy_mode_t copy_mode);


/**
 * Creates and returns a coo_matrix_t object, uninitialized, with nnz space reserved.
 *
 * @note automatically sets ownership to LIBRARY_DEALLOCATES and the
 * deallocator to free().  If ownership was previously
 * LIBRARY_DEALLOCATES but the deallocator not the default, then the
 * non-default deallocator is used to deallocate the original arrays,
 * but the new arrays are allocated with the default deallocator
 * (malloc).
 */
struct coo_matrix_t*
reserve_coo_matrix (const int m, const int n, const int nnz, 
		    enum index_base_t index_base,
		    enum symmetry_type_t symmetry_type, 
		    enum symmetric_storage_location_t symmetric_storage_location,
                    enum value_type_t value_type);

/**
 * Returns true (nonzero) if A is a valid COO format sparse matrix, else returns zero.
 *
 * @param A [IN]  Coordinate (COO) format sparse matrix
 */
int
valid_coo_matrix_p (struct coo_matrix_t *A);


/**
 * Deallocates the struct A and the storage it contains.
 *
 * @param A [IN/OUT]  Valid pointer to a dynamically allocated
 *                    coo_matrix_t object.
 */
void
destroy_coo_matrix (struct coo_matrix_t* A);


/**
 * Returns a deep copy of the given COO format sparse matrix.
 *
 * @param A [IN]  Valid pointer to a coo_matrix_t object.
 */
struct coo_matrix_t*
copy_coo_matrix (const struct coo_matrix_t* A);


/**
 * Returns nonzero if the given COO format sparse matrices are identical
 * (i.e. deep copies of each other), else returns zero.  As a side effect,
 * sorts the entries of A and B.
 *
 * @param A [IN]  COO format sparse matrix
 * @param B [IN]  COO format sparse matrix
 */
int
coo_matrix_equal_p (struct coo_matrix_t *A, struct coo_matrix_t *B);


struct csc_matrix_t; /* forward declaration */


/**
 * Given a sparse matrix of CSC (compressed sparse column, aka Harwell-Boeing)
 * format, creates a copy of the matrix in COO (coordinate) format and returns
 * a pointer to the freshly allocated copy.
 *
 * @param A [IN]   Sparse matrix in compressed sparse column format
 * @param index_base [IN]  ZERO or ONE; the index base (C-style or
 *                         Fortran-style) of the output COO format 
 *                         sparse matrix.
 *
 * @return Freshly allocated copy of the input matrix, in COO format
 */
struct coo_matrix_t*
csc_to_coo_matrix (const struct csc_matrix_t* A, int index_base);
 

/**
 * Converts the indices in A from C indexing (zero-based) to Fortran
 * indexing (one-based).
 *
 * @param A [IN/OUT]  Sparse matrix in coordinate (COO) format
 */
void
coo_c_to_fortran (struct coo_matrix_t* A);


/**
 * Converts the indices in A from Fortran indexing (one-based) to C indexing 
 * (zero-based).
 *
 * @param A [IN/OUT]  Sparse matrix in coordinate (COO) format
 */
void
coo_fortran_to_c (struct coo_matrix_t* A);


/**
 * Allocates storage in the given struct for the matrix.
 * 
 * @param A [IN]  Coordinate format sparse matrix
 * @param m [IN]  Number of rows in the matrix A
 * @param n [IN]  Number of columns in the matrix A
 * @param nnz [IN]  Number of nonzero entries in the matrix A
 *
 * @note This automatically sets ownership to LIBRARY_DEALLOCATES 
 *       and the deallocator to the default.
 */
void
alloc_coo_matrix (struct coo_matrix_t* A, const int m, const int n, const int nnz,
		  enum index_base_t index_base, 
		  enum symmetry_type_t symmetry_type, 
		  enum symmetric_storage_location_t symmetric_storage_location,
                  enum value_type_t value_type);

/**
 * Deallocates the dynamic memory used to store the matrix.
 */
void 
dealloc_coo_matrix (struct coo_matrix_t* A);


/** 
* Checks whether the indices in the given COO-format sparse matrix are 
* within their valid ranges, for Fortran-style (1-based) indices.
*
* @param m [IN]    Number of rows in the matrix
* @param n [IN]    Number of columns in the matrix
* @param nnz [IN]  Number of structural nonzeros (stored elements)
* @param II [IN]    Array of row indices
* @param JJ [IN]    Array of column indices
*
* @return Nonzero if indices in range, 0 if an index is out of range.
*/
int
coo_matrix_in_fortran_format_p (const int m, const int n, const int nnz, 
				const int II[], const int JJ[]);


/**
 * Prints the given COO matrix to the given file stream, in MatrixMarket format.
 *
 * @param out [OUT]   Valid file stream
 * @param A [IN]      Sparse matrix in COO format
 *
 * @return Nonzero if error, else zero.
 */
int
print_coo_matrix_in_matrix_market_format (FILE* out, const struct coo_matrix_t* A);


/**
 * Save the given COO matrix to the given file, in Harwell-Boeing format.
 *
 * @param filename [IN]   Filename to which to save the given sparse matrix
 * @param A [IN]          Sparse matrix in COO format to save
 *
 * @return Nonzero if error, else zero.
 */
int
save_coo_matrix_in_harwell_boeing_format (const char* const filename, 
					  const struct coo_matrix_t* A);


/**
 * Save the given COO matrix to the given file, in MatrixMarket format.
 *
 * @param filename [IN]   Filename to which to save the given sparse matrix
 * @param A [IN]          Sparse matrix in COO format to save
 *
 * @return Nonzero if error, else zero.
 */
int
save_coo_matrix_in_matrix_market_format (const char* const filename, 
					 const struct coo_matrix_t* A);


/**
 * Prints the given COO matrix to the given file stream, in Matlab format.
 *
 * @param out [OUT]   Valid file stream
 * @param A [IN]      Sparse matrix in COO format
 *
 * @return Nonzero if error, else zero.
 */
int
print_coo_matrix_in_matlab_format (FILE* out, struct coo_matrix_t* A);

/**
 * Save the given COO matrix to the given file, in Matlab ASCII format.
 *
 * @param filename [IN]   Filename to which to save the given sparse matrix
 * @param A [IN]          Sparse matrix in COO format to save
 *
 * @return Nonzero if error, else zero.
 */
int
save_coo_matrix_in_matlab_format (const char* const filename,
				  struct coo_matrix_t* A);

/**
 * Loads a COO-format sparse matrix from the given file, in Matlab ASCII format.
 *
 * @param filename [IN]  Name of file from which to load the sparse matrix
 * 
 * @return NULL if error, else pointer to COO-format sparse matrix
 */
struct coo_matrix_t* 
load_coo_matrix_in_matlab_format (const char* const filename);


/**
 * Takes a sparse matrix in COO format that uses a symmetric form of storage, 
 * and expands it "in place" into an unsymmetric storage format (so that all 
 * the nonzeros are stored explicitly).
 *
 * @param A [IN/OUT]   Sparse matrix in COO format
 * 
 * @return Nonzero if error, else zero.
 */
int 
coo_matrix_expand_symmetric_storage (struct coo_matrix_t* A);

/**
 * Takes a sparse matrix in COO format that uses a symmetric form of storage, 
 * and expands it into an unsymmetric storage format (so that all the nonzeros 
 * are stored explicitly).
 *
 * @param A [IN/OUT]   Sparse matrix in COO format
 * 
 * @return NULL if error, else the expanded sparse matrix in COO format
 */
struct coo_matrix_t*  
coo_matrix_expand_symmetric_storage_copy (struct coo_matrix_t* A);

 
/* 
 * Matlab 7 actually won't save the imaginary parts of complex variables to
 * ASCII format data files.  This is obviously a MASSIVE BUG on Matlab's part. 
 */

#endif /* _coo_matrix_h */
