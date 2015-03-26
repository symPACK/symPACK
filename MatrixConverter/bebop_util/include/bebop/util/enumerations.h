#ifndef _enumerations_h
#define _enumerations_h
/**
 * @file enumerations.h
 * @author Mark Hoemmen
 * @since 26 May 2005
 * @date Time-stamp: <2008-07-16 10:16:25 mhoemmen>
 *
 * An include file for all the enums and #defines that are globally
 * useful.  Many of these pertain to the storage of sparse matrices.
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
 * Enum for base of array indices:  either zero-based (C-style) or one-based 
 * (Fortran-style).
 */
enum
index_base_t {ZERO = 0, ONE = 1};

/**
 * Symmetry type of the matrix.  Unfortunately MatrixMarket doesn't have a 
 * skew-Hermitian type, otherwise we would include one.
 */
enum
symmetry_type_t 
{ 
  UNSYMMETRIC = 0, SYMMETRIC, SKEW_SYMMETRIC, HERMITIAN 
};

/**
 * For matrices with a type of symmetry:  Where the elements of the matrix 
 * are actually stored: in the lower triangle or the upper triangle.  If 
 * the matrix is UNSYMMETRIC, this value is ignored.
 */
enum
symmetric_storage_location_t
{
  UPPER_TRIANGLE, LOWER_TRIANGLE
};

/**
 * REAL means that the entries of the matrix are real-valued (IEEE 754 double
 * precision).
 * COMPLEX means that they are complex-valued (using C99 "double _Complex"
 * or an equivalent defined type in __complex.h).
 * PATTERN means that no values are stored, only the nonzero pattern.
 */
enum
value_type_t
{
  REAL, COMPLEX, PATTERN
};

/**
 * \brief Input matrix ownership modes.
 *
 * On matrix handle creation, users can specify how they want the
 * library to deal with the input arrays defining the sparse matrix.
 * Note that any conversion of the sparse matrix (e.g. from CSR to
 * CSC) means that the original input arrays are no longer a
 * representation of the sparse matrix pointed to by the handle.  The
 * library therefore takes responsibility for deallocating this new
 * representation upon destruction of the matrix handle.
 *
 * - LIBRARY_DEALLOCATES:  the Sparse Matrix Converter library is 
 *   responsible for deallocating the input arrays
 * - USER_DEALLOCATES:  the user is responsible for deallocating the
 *   input arrays
 */
enum
ownership_mode_t
{
  LIBRARY_DEALLOCATES,
  USER_DEALLOCATES
};

/**
 * \brief Input matrix copy modes.
 *
 * On matrix handle creation using input arrays (e.g. calling
 * create_csr_matrix by giving it val, ind and ptr arrays), NO_COPY
 * indicates that the input arrays should be used directly, and COPY
 * indicates that the input arrays should be copied immediately (thus
 * allowing the user to free() them as soon as e.g. create_csr_matrix
 * returns).
 */
enum
copy_mode_t
  {
    NO_COPY,
    COPY
  };




#endif /* _enumerations_h */
