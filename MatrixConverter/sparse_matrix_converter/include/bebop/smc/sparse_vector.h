#ifndef _sparse_vector_h
#define _sparse_vector_h
/**
 * @file sparse_vector.h
 * @author Mark Hoemmen
 * @since 06/16/04 16:28:39 PDT
 * @date Time-stamp: <2008-07-16 10:57:31 mhoemmen>
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
 * facilitates appending to the end of the vector (append_to_sparse_vector()).
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
sparse_vector_t
{
  /** Number of entries in the vector */
  int     len;

  /** Number of entries reserved for potential use */
  int     maxlen;

  /** Array of values */
  double* val;

  /** Array of indices */
  int*    idx;
};


/**
 * Scatters the given sparse vector into a dense vector.  Assumes that the 
 * entries of the sparse vector are sorted in increasing order of index.
 *
 * @param dest [OUT]   Dense 1-D vector
 * @param src [IN]     Sparse vector to scatter into dest
 * @param n [IN]       Length of the dense vector (including structural zeros)
 */
void
scatter_sparse_vector (double* dest, const struct sparse_vector_t* src, 
					   const int n);

/**
 * Gathers the given dense vector into a sparse vector, using the given
 * tolerance to decide if an entry of the dense vector counts as zero.
 *
 *
 * @param sv [OUT]  Sparse vector
 * @param dv [IN]   Dense vector
 * @param n  [IN]   Length of dense vector
 * @param tol [IN]  Tolerance for dropping nearly-zero entries of dv.
 *                  Values less than or equal to this are considered zero 
 *                  entries and are not inserted into the sparse vector.
 *                  Must be non-negative.
 */
void
gather_sparse_vector (struct sparse_vector_t* sv, const double* dv, 
					  const int n, const double tol);
  

/**
 * Equivalent to ``v->idx[i] = index_value'', with bounds checking.
 */
void
set_sparse_vector_index (struct sparse_vector_t* v, const int i, 
						 const int index_value);


/**
 * Equivalent to the rvalue ``v->idx[i]'', with bounds checking.
 */
int
get_sparse_vector_index (struct sparse_vector_t* v, const int i);



/**
 * Allocates the struct itself, then allocates initial_length entries in it.
 * Returns a pointer to the struct.
 *
 * @return Pointer to the allocated sparse vector object
 */
struct sparse_vector_t* 
create_sparse_vector (const int initial_length);

/**
 * Returns a pointer to a deep copy of the given sparse vector.
 */
struct sparse_vector_t*
clone_sparse_vector (const struct sparse_vector_t* src);

/**
 * Frees the contents of v, and then frees the struct v itself.
 *
 * @param v [OUT]  
 *
 * @warning Pointer invalid after calling this function. 
 */
void
destroy_sparse_vector (struct sparse_vector_t* v);



/** 
 * Allocates space in the given sparse vector to hold initial_length entries.
 * Does not allocate the struct itself.
 *
 * @param v [IN/OUT]  Pointer to the sparse vector in which space is to be
 *                    allocated.
 * @param initial_length [IN]  Initial number of entries in the sparse vector.
 */
void
init_sparse_vector (struct sparse_vector_t* v, const int initial_length);

/**
 * Frees the contents of v, but not v itself.
 *
 * @param v [OUT]
 */
void
deinit_sparse_vector (struct sparse_vector_t* v);

/**
 * Returns the number of entries in the sparse vector.
 *
 * @param v [IN]
 *
 * @return Number of entries in v
 */
int
length_sparse_vector (const struct sparse_vector_t* v);

/**
 * Prints the contents of v to the given output stream.
 *
 * @param out [OUT]  Valid output stream
 * @param v [IN]     Sparse vector
 */
void
print_sparse_vector (FILE* out, const struct sparse_vector_t* v);


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
print_sparse_vector_with_line_starter (FILE* out, 
									   const struct sparse_vector_t* v, 
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
resize_sparse_vector (struct sparse_vector_t* v, const int newlength);


/**
 * Appends the given (value,index) pair to the end of the given sparse vector.
 *
 * @param v [OUT]
 * @param value [IN]
 * @param index [IN]
 */
void
append_to_sparse_vector (struct sparse_vector_t* v, 
						 const double value, const int index);

/**
 * Equivalent to ``v->val[i] = value''.
 */
void
set_sparse_vector_value (struct sparse_vector_t* v, const int i, 
						 const double value);

/**
 * Equivalent to the r-value ``v->val[i]''.
 */
double
get_sparse_vector_value (const struct sparse_vector_t* v, const int i);


/**
 * Returns the dot product of a sparse vector (represented by (val,idx) with
 * len elements) and a dense vector x, in the index range [start, end).
 *
 * @param val [IN]  Values of sparse vector
 * @param idx [IN]  Indices of sparse vector
 * @param x   [IN]  Dense vector
 *
 * @return The dot product of the sparse vector and the dense vector,
 *         unless the index range is empty, in which case it returns zero.
 */
double
ddot_svdv (const double* val, const int* idx, const double* x,
           const int start, const int end);

/**
 * Returns the infinity norm of the given sparse vector.
 *
 * \warn UNTESTED!!!
 */
double
infinity_norm_sparse_vector (const struct sparse_vector_t* v);


/**
 * Given a sparse vector x and a tolerance ``tol'', returns a new sparse 
 * vector which is a copy of the original x, but with all elements smaller 
 * than ``tol'' removed.
 *
 * \warn UNTESTED!!!
 */
struct sparse_vector_t* 
filter_small_elements (const struct sparse_vector_t* x, const double tol);

/**
 * Swaps the insides of the two data structures (a shallow swap).
 * Assumes that both a and b are valid (initialized) sparse vector objects.
 *
 * @param a [IN/OUT]  Sparse vector which will get the contents of b
 * @param b [IN/OUT]  Sparse vector which will get the contents of a
 */
void
swap_sparse_vector (struct sparse_vector_t *a, struct sparse_vector_t *b);


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
coalesce_sparse_vector_entries_with_common_indices (struct sparse_vector_t *v);

/**
 * Sorts the entries of v in increasing order of their indices.
 *
 * @param v [IN/OUT]  Sparse vector
 */
void
sort_sparse_vector_by_index (struct sparse_vector_t *v);



#endif /* _sparse_vector_h */
