#ifndef _jad_matrix_h
#define _jad_matrix_h
/**
 * @file jad_matrix.h
 * @author Mark Hoemmen
 * @since 05 Jul 2005
 * @date Time-stamp: <2008-07-16 10:56:16 mhoemmen>
 *
 * Jagged diagonal (JAD) format sparse matrix: struct and member functions.
 *
 * @note Original author: Ankit Jain (ankit@berkeley.edu), later revised and 
 * expanded by Mark Hoemmen (mhoemmen@cs.berkeley.edu).
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

struct jad_matrix_t 
{
  /**
   * Number of "jagged diagonals."
   */
  int n_jagged_diagonals;

  /** 
   * Stores the nonzero values of the matrix.  It is NULL if value_type == 
   * PATTERN, it holds doubles if value_type == REAL, and it holds 
   * double_Complex if value_type == COMPLEX. 
   */
  void* jad_value;

  int * jad_col_index;
  int* jad_diag_start;
  int* jad_prm_nto;
  /**
   * The permuted y vector.
   */
  double* Py;

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
   * Indicates the type of the entries in jad_value:  REAL means "double", COMPLEX
   * means "double _Complex", PATTERN means the jad_value array is NULL and contains
   * no entries.
   */ 
  enum value_type_t value_type;

  /**
   * The ownership mode: tells whether this library or the user is
   * responsible for deallocating input arrays.
   */
  enum ownership_mode_t ownership;

  /**
   * The deallocation function to be called on the storage of this 
   * matrix, if ownership == LIBRARY_DEALLOCATES.
   */
  void (*deallocator) (void*);
};

struct csr_matrix_t; /* forward declaration */

/**
 * Returns a copy of the given sparse matrix A in jagged diagonal format.
 */
struct jad_matrix_t* 
csr_to_jad (struct csr_matrix_t* A);

/**
 * Returns a copy of the given sparse matrix A in CSR format.
 */
struct csr_matrix_t*
jad_to_csr (struct jad_matrix_t* A);

/**
 * Deallocates the storage used by the given jad_matrix_t object, and frees
 * the jS pointer.
 */
void 
destroy_jad_matrix (struct jad_matrix_t* jS);

/**
 * Deallocates the storage used by the given jad_matrix_t object, without 
 * actually freeing the jS pointer.
 */
void
dealloc_jad_matrix (struct jad_matrix_t* jS);


#endif /* _jad_matrix_h */
