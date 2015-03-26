#ifndef _spvec_h
#define _spvec_h
/**
 * @file spvec.h
 * @author Mark Hoemmen
 * @since 07 Jun 2006
 * @date Time-stamp: <2008-07-16 10:57:40 mhoemmen>
 * 
 * A sparse vector object.
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
 * A sparse vector.  It can be thought of as a degenerate case of a CSC or CSR
 * matrix, where colptr[] resp. rowptr[] is no longer necessary, since there is
 * only one column resp. row.  idx[j] is the index of the value val[j] in the
 * vector.  Entries may not necessarily be sorted by index; this depends on how
 * one constructs the vector.  
 *
 * The vector is initially allocated with maxlen space; maxlen >= len, which 
 * facilitates appending to the end of the vector (append_to_spvec()).
 * This makes the append operation amortized O(1).  This is a desirable 
 * property since we use the sparse vector to extract a row from a sparse 
 * matrix; the matrix may be very large, so we want to allocate the array 
 * space dynamically, rather than trying to reserve some huge space at the 
 * beginning and return it once we have figured out exactly how much we need.
 *
 * @note Even though the struct itself is lightweight, it needs to be passed
 * around by pointer, so that the internal length, etc. fields have the correct
 * pass-by-reference semantics.
 */
struct 
spvec_t
{
  /** Number of entries in the vector */
  int     len;

  /** Number of entries reserved for potential use */
  int     maxlen;

  /** Array of values */
  void* val;

  /** Array of indices */
  int* idx;
  
  /** 
   * Type of the values: REAL means doubles, COMPLEX means
   * double_Complex, and PATTERN means no values are stored (just
   * indices). 
   */
  enum value_type_t value_type;
};


/**
 * Scatters the given sparse vector into a dense vector, performing
 * dest[src->idx] = dest[src->idx] + src->val.  Assumes that the 
 * entries of the sparse vector are sorted in increasing order of index.
 * We assume that the values in dest are of the same type as the values
 * in src.  If src->value_type is PATTERN, this is a noop.
 *
 * @param dest [OUT]   Dense 1-D vector
 * @param src [IN]     Sparse vector to scatter into dest
 * @param n [IN]       Length of the dense vector (including structural zeros)
 */
void
scatter_spvec (void* dest, const struct spvec_t* src, 
	       const int n);

/**
 * Equivalent to ``v->idx[i] = index_value'', with bounds checking.
 */
void
set_spvec_index (struct spvec_t* v, const int i, 
		 const int index_value);


/**
 * Equivalent to the rvalue ``v->idx[i]'', with bounds checking.
 */
int
get_spvec_index (struct spvec_t* v, const int i);


/**
 * Allocates the struct itself, then allocates initial_length entries in it.
 * Returns a pointer to the struct.
 *
 * @return Pointer to the allocated sparse vector object
 */
struct spvec_t* 
create_spvec (const int initial_length);

/**
 * Returns a pointer to a deep copy of the given sparse vector.
 */
struct spvec_t*
clone_spvec (const struct spvec_t* src);

/**
 * Frees the contents of v, and then frees the struct v itself.
 *
 * @param v [OUT]  
 *
 * @warning Pointer invalid after calling this function. 
 */
void
destroy_spvec (struct spvec_t* v);



/** 
 * Allocates space in the given sparse vector to hold initial_length entries.
 * Does not allocate the struct itself.
 *
 * @param v [IN/OUT]  Pointer to the sparse vector in which space is to be
 *                    allocated.
 * @param initial_length [IN]  Initial number of entries in the sparse vector.
 */
void
init_spvec (struct spvec_t* v, const int initial_length);

/**
 * Frees the contents of v, but not v itself.
 *
 * @param v [OUT]
 */
void
deinit_spvec (struct spvec_t* v);

/**
 * Returns the number of entries in the sparse vector.
 *
 * @param v [IN]
 *
 * @return Number of entries in v
 */
int
length_spvec (const struct spvec_t* v);

/**
 * Prints the contents of v to the given output stream.
 *
 * @param out [OUT]  Valid output stream
 * @param v [IN]     Sparse vector
 */
void
print_spvec (FILE* out, const struct spvec_t* v);


/**
 * Prints the contents of v to the given output stream, 
 * prepending each line with line_starter.
 *
 * @param out         [OUT]  Valid output stream
 * @param v            [IN]  Sparse vector
 * @param line_starter [IN]  zero-terminated string by which to prepend 
 *                           each line of output
 */
void
print_spvec_with_line_starter (FILE* out, 
			       const struct spvec_t* v, 
			       const char* line_starter);


/**
 * Resizes the given sparse vector to have a given number of entries.
 *
 * @param v [IN/OUT]      Previously initialized sparse vector
 * @param newlength [IN]  New number of entries in the vector
 *
 * @warning Sparse vector must be initialized first!
 */
void
resize_spvec (struct spvec_t* v, const int newlength);


/**
 * Appends the given (value,index) pairs to the end of the given sparse vector.
 *
 * @param v [OUT]
 * @param values [IN]
 * @param indices [IN]
 * @param num_to_append [IN]
 */
void
append_to_spvec (struct spvec_t* v, const void* values, 
		 const int* indices, const int num_to_append);

/**
 * Given a sparse vector x and a tolerance ``tol'', returns a new sparse 
 * vector which is a copy of the original x, but with all elements smaller 
 * than ``tol'' in absolute value removed.
 *
 * \warn UNTESTED!!!
 */
struct spvec_t* 
filter_small_elements (const struct spvec_t* x, const double tol);

/**
 * Swaps the insides of the two data structures (a shallow swap).
 * Assumes that both a and b are valid (initialized) sparse vector objects.
 *
 * @param a [IN/OUT]  Sparse vector which will get the contents of b
 * @param b [IN/OUT]  Sparse vector which will get the contents of a
 */
void
swap_spvec (struct spvec_t *a, struct spvec_t *b);


/**
 * Coalesces entries of the sparse vector which have the same index.
 * For example, if (3.0,10) and (-5.0,10) are two (val,idx) pairs in
 * the sparse vector, then the two entries are combined into a single
 * entry (-2.0,10).  As a side effect, sorts the sparse vector's entries
 * by index.
 *
 * @param v [IN/OUT]  Sparse vector
 */
void
coalesce_spvec_entries_with_common_indices (struct spvec_t *v);

/**
 * Sorts the entries of v in increasing order of their indices.
 *
 * @param v [IN/OUT]  Sparse vector
 */
void
sort_spvec_by_index (struct spvec_t *v);


/**
 * z := a*x + b*y 
 */
void
add_spvec (struct spvec_t* z,
	   const void* a,
	   struct spvec_t* x,
	   const void* b,
	   struct spvec_t* y);

/**
 * Copies the contents of the given sparse vector x into the arrays
 * ind (for indices) and val (for values), and returns the number of
 * (index, value) pairs written.
 */
int
export_spvec (int* ind, void* val, struct spvec_t* x);


void
append_spvec_to_csr_matrix (int** ptr, int** ind, void** val,
			    int* current_nnz,
			    int* nnz_upper_bound,
			    struct spvec_t* row,
			    const int which_row);



#endif /* _spvec_h */
