#ifndef _bcoo_matrix_h
#define _bcoo_matrix_h
/**
 * @file bcoo_matrix.h
 * @author Mark Hoemmen
 * @since 21 Feb 2005
 * @date Time-stamp: <2008-07-16 10:50:06 mhoemmen>
 *
 * Wrapper for a block coordinate format sparse matrix.
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
 * A block coordinate (COO) format sparse matrix.
 *
 * \note This data structure is suited for low amortized cost of dynamic addition 
 * of entries, but isn't necessarily as space-efficient as it could be.
 */
struct
bcoo_matrix_t
{
  /** Number of block rows */
  int bm;

  /** Number of block columns */
  int bn;

  /** Number of rows in each dense subblock */
  int r;

  /** Number of columns in each dense subblock */
  int c;

  /** 
   * Number of stored (nonzero) blocks.  If the matrix is stored in a 
   * symmetric (or skew, etc.) format, nnz only refers to the number of 
   * stored blocks, not the actual number of nonzero blocks. 
   */
  int nnzb;

  /** 
   * Maximum number of stored (nonzero) blocks that the data structure can 
   * currently hold.  Invariant:  nnzb <= nnzb_upper_bound.
   */
  int nnzb_upper_bound;

  /** 
   * For block e, II[e] is the actual (not block) row index of the upper left 
   * corner of the block (we can't call the array "I" because that is reserved 
   * as a macro in C99 for sqrt{-1}).
   */
  int *II;

  /** 
   * For block e, JJ[e] is the actual (not block) column index of the upper 
   * left corner of the block.
   */
  int *JJ;

  /** 
   * For entry e, val[r*c*e : r*c*e + r*c - 1] holds its values.  Type of the 
   * entries depends on the value of value_type.  If value_type == PATTERN, 
   * then val is NULL (and holds nothing).
   */
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

  /** 
   * 0 if the values in each nonzero block are stored row-oriented, 1 if 
   * column-oriented.  If value_type==PATTERN, then this field is ignored. 
   */
  int col_oriented_p;

  /**
   * Indicates who is responsible for deallocating the input arrays.
   */
  enum ownership_mode_t ownership;

  /**
   * Function for deallocating the input arrays, if LIBRARY_DEALLOCATES.
   */
  void (*deallocator) (void*);
};


/**
 * Prints the given BCOO matrix in MatrixMarket coordinate format to the 
 * given output file.
 *
 * @param out [OUT]  Valid output file
 * @param A   [IN]   BCOO format matrix 
 *
 * @return Zero if no error, else nonzero (nothing is printed if there is an 
 * error).
 */
int
print_bcoo_matrix_in_matrix_market_format (FILE *out, struct bcoo_matrix_t* A);


/**
 * Creates an initially empty bcoo_matrix object with the given 
 * characteristics.
 *
 * @param bm [IN]  Number of block rows
 * @param bn [IN]  Number of block columns
 * @param r [IN]   Number of rows in each dense subblock
 * @param c [IN]   Number of columns in each dense subblock
 *
 * @return Dynamically allocated bcoo_matrix object
 */
struct bcoo_matrix_t*
create_bcoo_matrix (const int bm, 
		    const int bn, 
		    const int r, 
		    const int c, 
		    const int nnzb, 
		    const int nnzb_upper_bound, 
		    int *II, 
		    int *JJ, 
		    void* val, 
		    enum index_base_t index_base, 
		    enum symmetry_type_t symmetry_type, 
		    enum symmetric_storage_location_t symmetric_storage_location, 
		    enum value_type_t value_type,
		    enum ownership_mode_t ownership,
		    void (*deallocator) (void*),
		    enum copy_mode_t copy_mode);

/**
 * Ownership is of course LIBRARY_DEALLOCATES and copy_mode is NO_COPY.
 */
struct bcoo_matrix_t*
create_empty_bcoo_matrix (const int bm, 
			  const int bn, 
			  const int r, 
			  const int c, 
			  enum index_base_t index_base, 
			  enum symmetry_type_t symmetry_type, 
			  enum symmetric_storage_location_t symmetric_storage_location, 
			  enum value_type_t value_type);

/**
 * Destroys the storage used by the given bcoo_matrix_t, 
 * and destroys the handle as well.
 *
 * @param A [IN/OUT]
 */
void
destroy_bcoo_matrix (struct bcoo_matrix_t* A);


/**
 * Reallocates the II, JJ and val arrays in the bcoo_matrix struct to be 
 * of length newmaxlength.  Adjusts nnzb only if necessary (i.e. only if the 
 * size of the arrays is to be reduced).  This is a helper function for 
 * bcoo_matrix_resize.  
 *
 * @warn If ownership is currently USER_DEALLOCATES, this function changes 
 * the ownership to LIBRARY_DEALLOCATES.
 *
 * @see bcoo_matrix_resize
 *
 * @param A [IN/OUT] 
 * @param newmaxlength [IN]   New length of the II and JJ arrays
 */
void
bcoo_matrix_realloc (struct bcoo_matrix_t* A, const int newmaxlength);

/**
 * Resizes the II, JJ and val arrays in the bcoo_matrix_t struct to have 
 * the new given length.  Allocates more memory (thus adjusting 
 * nnzb_upper_bound) only if necessary.
 *
 * @warn If ownership is currently USER_DEALLOCATES, this function changes 
 * the ownership to LIBRARY_DEALLOCATES.
 *
 * @param A [IN/OUT]
 * @param newlength [IN]  New lengths of (i.e. number of entries stored in) 
 *                        the II, JJ and val arrays.  Has nothing to do with the 
 *                        actual amount of allocated space for these arrays.  
 *                        nnzb_upper_bound space is allocated, but nnzb entries
 *                        are stored.
 */
void
bcoo_matrix_resize (struct bcoo_matrix_t* A, const int newlength);

/**
 * Adds an entry to the given sparse matrix structure representation.
 *
 * @param A [IN/OUT]  Sparse matrix data structure in block coordinate 
 *                    (BCOO) format
 * @param bi [IN]   Block row index at which to add the block (block row 
 *                  index of the upper left corner of the block)
 * @param bj [IN]   Block column index at which to add the block (block 
 *                  column index of the upper left corner of the block)
 * @param value [IN]  Pointer to the value to add (type depends on 
 *                    A->value_type).  If value_type is PATTERN, you can
 *                    set value==NULL.
 *
 * @warn If ownership is currently USER_DEALLOCATES, this function changes 
 * the ownership to LIBRARY_DEALLOCATES, as increasing the length of the
 * internal arrays may require reallocation.
 */
void
bcoo_matrix_add_entry (struct bcoo_matrix_t* A, const int bi, const int bj, 
		       const void* value);

/**
 * If A is valid, returns 1, else returns 0.  By "valid," we mean that 
 * all the dimension parameters are within reasonable ranges and all
 * the indices are within their appropriate ranges.
 *
 * @param A [IN]   Pointer to a BCOO matrix structure object
 *
 * @return One if A is valid, else zero.
 */
int
valid_bcoo_matrix_p (struct bcoo_matrix_t* A);


struct bcsr_matrix_t; /* forward declaration */

/**
 * Given a sparse matrix in block coordinate (BCOO) format, generates a block 
 * compressed sparse row (BCSR) matrix with the same structure and real-valued
 * random entries from a uniform [-1,1] distribution.  As a side effect, 
 * coalesces duplicate entries from the BCOO matrix.
 *
 * @param A [IN]  Sparse matrix structure in BCOO format
 *
 * @return BCSR sparse matrix with the same structure as the input and
 *         real-valued random entries; if error, returns NULL.
 */
struct bcsr_matrix_t*
bcoo_matrix_to_random_bcsr_matrix (struct bcoo_matrix_t* A);


/**
 * Given a sparse matrix in block coordinate (BCOO) format, generates a copy 
 * in BCSR (block compressed sparse row) format.  As a side effect, coalesces 
 * duplicate entries from the BCOO matrix.
 *
 * @param A [IN]  Sparse matrix in BCOO format
 *
 * @return BCSR format sparse matrix, which is a copy of the input matrix
 */
struct bcsr_matrix_t*
bcoo_to_bcsr (struct bcoo_matrix_t* A);

struct bcoo_matrix_t*
bcsr_to_bcoo (struct bcsr_matrix_t* A);

/**
 * Saves the given BCOO format sparse matrix to the given file in 
 * MatrixMarket format.
 *
 * @param filename [IN]  Name of file to which to save the matrix
 * @param A [IN]         The matrix to save
 *
 * @return Nonzero if error, else zero.
 */
int
save_bcoo_matrix_in_matrix_market_format (const char* const filename, 
					  struct bcoo_matrix_t* A);

/**
 * Coalesces duplicate entries in the given BCOO format sparse matrix.
 * For example, if the original matrix has two entries:  (30, 45, 1.23) and 
 * (30, 45, 4.56), the resulting matrix will have one entry:  (30, 45, 5.79).
 *
 * @param A [IN/OUT]          Sparse matrix in BCOO format
 * @param num_removed [OUT]   The number of duplicate entries removed from A
 */
void
bcoo_matrix_coalesce_duplicate_entries (struct bcoo_matrix_t* A, int* num_removed);


#endif /* _bcoo_matrix_h */
