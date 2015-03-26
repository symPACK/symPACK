#ifndef _bcsr_matrix_h
#define _bcsr_matrix_h
/**
 * @file bcsr_matrix.h
 * @author Mark Hoemmen
 * @since 22 Feb 2005
 * @date Time-stamp: <2008-07-16 10:50:11 mhoemmen>
 * 
 * BCSR (register blocked CSR) format sparse matrix data structure and 
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

#include <bebop/util/enumerations.h>
#include <stdio.h> /* FILE */


/**
 * Sparse matrix in Block Compressed Sparse Row (BCSR) format.
 */
struct
bcsr_matrix_t
{
  /** 
   * Number of block rows in the matrix (actual number of rows is r*bm) 
   */
  int bm;
  
  /** 
   * Number of block columns in the matrix (actual number of columns is
   * c*bn)
   */
  int bn;

  /**
   * Number of rows in each (dense) subblock
   */
  int r;

  /**
   * Number of columns in each (dense) subblock
   */
  int c;
  
  /** 
   * Number of stored (nonzero) blocks in the matrix.  If the matrix is 
   * stored in a symmetric (or skew, etc.) format, nnz only refers to the 
   * number of stored entries, not the actual number of nonzeros. 
   */
  int nnzb;

  /** Array of stored (nonzero) entries of the matrix */
  void* values;

  /** 
   * Array of column indices of the upper left corners of the nonzero blocks 
   * of the matrix.  These are the actual column indices and not the "block 
   * indices". 
   */
  int* colind;

  /** Array of indices into the colind and values arrays, for each row */
  int* rowptr;

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

  /** 0 if the nonzero blocks are row-oriented, 1 if column-oriented */
  int col_oriented_p;

  enum ownership_mode_t ownership;
  void (*deallocator) (void*);
};



/**
 * Fills in the data structure with shallow copies of the given arguments.
 *
 * @param A [OUT]     The struct to fill in
 * @param bm [IN]     Number of block rows
 * @param bn [IN]     Number of block columns
 * @param r [IN]      Number of rows in each dense subblock
 * @param c [IN]      Number of columns in each dense subblock
 * @param nnzb [IN]   Number of nonzero dense subblocks
 * @param values [IN] Array of values (length nnzb * r * c, or NULL if this is a PATTERN matrix)
 * @param colind [IN] Array of column indices (length nnzb)
 * @param rowptr [IN] Array of indices for each row (length bm)
 * @param symmetry_type [IN]
 * @param symmetric_storage_location [IN]
 * @param value_type [IN]
 * @param col_oriented_p [IN]  0 if the nonzero blocks are row-oriented,
 *                             1 if column-oriented
 */
void
init_bcsr_matrix (struct bcsr_matrix_t* A, const int bm, const int bn, 
		  const int r, const int c, const int nnzb, void* values, 
		  int* colind, int* rowptr, 
		  const enum symmetry_type_t symmetry_type,
		  const enum symmetric_storage_location_t symmetric_storage_location,
		  const enum value_type_t value_type, 
		  const int col_oriented_p,
		  enum ownership_mode_t ownership,
		  void (*deallocator) (void*),
		  enum copy_mode_t copy_mode);


/**
 * Unpacks the CSR format sparse matrix into the given data structures.
 *
 * @param A [IN]     The struct from which to unpack
 * @param bm [OUT]     Number of block rows
 * @param bn [OUT]     Number of block columns
 * @param r [OUT]      Number of rows in each dense subblock
 * @param c [OUT]      Number of columns in each dense subblock
 * @param nnzb [OUT]   Number of nonzero dense subblocks
 * @param values [OUT] Array of values (length nnzb * r * c, unless a PATTERN matrix)
 * @param colind [OUT] Array of column indices (length nnzb)
 * @param rowptr [OUT] Array of indices for each row (length bm)
 * @param symmetry_type [OUT]
 * @param symmetric_storage_location [OUT]
 * @param value_type [OUT]
 * @param col_oriented_p [OUT]  0 if the nonzero blocks are row-oriented,
 *                              1 if column-oriented
 * 
 * @warn FIXME: does not yet incorporate things like ownership and the
 *       deallocator!
 */
void
unpack_bcsr_matrix (struct bcsr_matrix_t* A, int* bm, int* bn, int* r, int* c,
		    int* nnzb, void** values, int** colind, int** rowptr,
		    enum symmetry_type_t* symmetry_type,
		    enum symmetric_storage_location_t* symmetric_storage_location,
		    enum value_type_t* value_type,
		    int* col_oriented_p);



/**
 * Allocates a handle of a CSR-format sparse matrix.
 */
struct bcsr_matrix_t* 
create_bcsr_matrix_handle ();


/**
 * Destroys a handle of a CSR-format sparse matrix, without touching the
 * internal storage of the matrix.
 */
void
destroy_bcsr_matrix_handle (struct bcsr_matrix_t* A);

/**
 * Creates a handle of a BCSR-format sparse matrix with the given attributes.
 *
 * @param bm [IN]     Number of block rows
 * @param bn [IN]     Number of block columns
 * @param r [IN]      Number of rows in each dense subblock
 * @param c [IN]      Number of columns in each dense subblock
 * @param nnzb [IN]   Number of nonzero dense subblocks
 * @param values [IN] Array of values (length nnzb * r * c)
 * @param colind [IN] Array of column indices (length nnzb)
 * @param rowptr [IN] Array of indices for each row (length bm)
 * @param symmetry_type [IN]
 * @param symmetric_storage_location [IN]
 * @param value_type [IN]
 * @param col_oriented_p [IN]  0 if the nonzero blocks are 
 *                      row-oriented, 1 if column-oriented 
 *
 * @return Pointer to dynamically allocated bcsr_matrix_t object 
 */
struct bcsr_matrix_t*
create_bcsr_matrix (const int bm, const int bn, const int r, const int c, 
		    const int nnzb, void* values, int* colind, int* rowptr,
		    const enum symmetry_type_t symmetry_type,
		    const enum symmetric_storage_location_t symmetric_storage_location, 
		    const enum value_type_t value_type,
		    const int col_oriented_p,
		    enum ownership_mode_t ownership,
		    void (*deallocator) (void*),
		    enum copy_mode_t copy_mode);

/**
 * Deallocates the internal storage of A, and then the handle A itself.
 */
void
destroy_bcsr_matrix (struct bcsr_matrix_t* A);

/**
 * Returns a deep copy of the sparse matrix A.
 */
struct bcsr_matrix_t*
clone_bcsr_matrix (struct bcsr_matrix_t* A);

/**
 * Prints the given BCSR-format sparse matrix A to the given output stream, in
 * MatrixMarket format.
 *
 * @return Zero if no error, else nonzero.
 */
int
print_bcsr_matrix_in_matrix_market_format (FILE* out, const struct bcsr_matrix_t* A);

/**
 * Returns 1 if A is a valid BCSR format sparse matrix, else returns 0.
 */
int
valid_bcsr_matrix_p (const struct bcsr_matrix_t* A);


/**
 * Dumps the contents of A to the given valid output stream.
 * Mainly useful for debugging.
 *
 * @param out [OUT]   Valid output stream
 * @param A   [IN]    Pointer to BCSR format sparse matrix object
 */
void 
dump_bcsr_matrix (FILE* out, const struct bcsr_matrix_t* A);


/**
 * Returns the n x n identity matrix, with 1 x 1 blocks.
 */
struct bcsr_matrix_t*
identity_bcsr_matrix (const int n);

/**
 * Returns a generalized bm*r x bn*c "identity" matrix, with r x c blocks.
 * Each block consists of all ones.
 */
struct bcsr_matrix_t*
extended_identity_bcsr_matrix (const int bm, const int bn, const int r, const int c);

/**
 * Saves the given BCSR-format sparse matrix A to the given filename, 
 * in MatrixMarket format.
 */
int 
save_bcsr_matrix_in_matrix_market_format (const char* const filename, 
					  struct bcsr_matrix_t* A);





#endif /* _bcsr_matrix_h */
