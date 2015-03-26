/**
 * @file csr_matrix.c
 * @author Mark Hoemmen
 * @since 31 May 2005
 * @date Time-stamp: <2008-07-16 11:15:52 mhoemmen>
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
#include <bebop/util/config.h>
#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/csr_expand_to_dense.h>
#include <bebop/smc/csr_matmatmult.h>
#include <bebop/smc/csr_spmv.h>
#include <bebop/smc/csr_triple_product.h>
#include <bebop/smc/csr_trisolve.h>
#include <bebop/smc/csr_weighted_jacobi.h>
#include <bebop/smc/coo_matrix.h>
#include <bebop/smc/csc_matrix.h>
#include <bebop/smc/iohb.h>
#include <bebop/smc/read_mm.h>
#include <bebop/smc/sparse_vector.h>

#include <bebop/util/complex.h>
#include <bebop/util/malloc.h>
#include <bebop/util/util.h>
#include <bebop/util/log.h>

#include <assert.h>
#include <math.h>   /* fabs */
#include <string.h> /* memcpy */

#ifdef USING_VALGRIND_CLIENT_REQUESTS
#  include <valgrind/memcheck.h>
#endif /* USING_VALGRIND_CLIENT_REQUESTS */



/**
 * An (index,value) pair.  In our context of sparse matrices in CSR format, 
 * it represents a single row index and its associated nonzero value.
 */
struct index_real_value_pair_t
{
  int index;
  double value;
};

struct index_complex_value_pair_t
{
  int index;
  double_Complex value;
};


/**
 * Compares two (index,value) pairs by their indices.
 */
static int
compare_index_real_value_pairs (const void *pa, const void *pb)
{
  const struct index_real_value_pair_t* a = (const struct index_real_value_pair_t*) pa;
  const struct index_real_value_pair_t* b = (const struct index_real_value_pair_t*) pb;

  if (a->index > b->index)
    return 1;
  else if (a->index < b->index)
    return -1;
  else
    return 0;
}

/**
 * Compares two (index,value) pairs by their indices.
 */
static int
compare_index_complex_value_pairs (const void *pa, const void *pb)
{
  const struct index_complex_value_pair_t* a = (const struct index_complex_value_pair_t*) pa;
  const struct index_complex_value_pair_t* b = (const struct index_complex_value_pair_t*) pb;

  if (a->index > b->index)
    return 1;
  else if (a->index < b->index)
    return -1;
  else
    return 0;
}

static int
compare_ints (const void *pa, const void *pb)
{
  int a = *((int*) pa);
  int b = *((int*) pb);

  if (a > b)
    return 1;
  else if (a < b)
    return -1;
  else 
    return 0;
}




/** 
 * Takes row number i of the given sparse matrix A, and copies all the 
 * indices (and values, if there are any values) of that column into "col".
 *
 * @param A [IN]    Sparse matrix in CSR format
 * @param i [IN]    Index of the row to copy (zero-based indices)
 * @param col [OUT] An array of either index_real_value_pair_t, 
 *                  index_complex_value_pair_t or int, depending on whether 
 *                  the value type of the matrix is REAL, COMPLEX, or PATTERN,
 *                  respectively.
 * @param max_nnz [IN]  The max number of nonzeros in any column of A
 */
static void
copy_row2pairs (const struct csr_matrix_t* A, int i, void* row, int max_nnz)
{
  int a = A->rowptr[i];
  int b = A->rowptr[i+1];
  int nnz = b - a;
  int k;

  assert (nnz <= max_nnz);

  if (A->value_type == REAL)
    {
      struct index_real_value_pair_t* _row = (struct index_real_value_pair_t*) row;
      const double* const values = (const double* const) (A->values);

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  _row[k].index = A->colidx[a+k];
	  _row[k].value = values[a+k];
	}
    }
  else if (A->value_type == COMPLEX)
    {
      struct index_complex_value_pair_t* _row = (struct index_complex_value_pair_t*) row;
      const double_Complex* const values = (const double_Complex* const) (A->values);

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  _row[k].index = A->colidx[a+k];
	  _row[k].value = values[a+k];
	}
    }
  else if (A->value_type == PATTERN)
    {
      int* _row = (int*) row;

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  _row[k] = A->colidx[a+k];
	}
    }
}


/** 
 * Given a sparse matrix A in CSR format, a row number i, and a list of
 * (index,value) pairs, copies the list of (index,value) pairs back into that
 * row of the matrix.
 */
static void
copy_pairs2row (const void* row, int max_nnz, struct csr_matrix_t* A, int i)
{
  int a = A->rowptr[i];
  int b = A->rowptr[i+1];
  int nnz = b - a;

  int k;

  assert (nnz <= max_nnz);

  if (A->value_type == REAL)
    {
      struct index_real_value_pair_t* _row = (struct index_real_value_pair_t*) row;
      double* values = (double*) (A->values);

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  A->colidx[a+k] = _row[k].index;
	  values[a+k] = _row[k].value;
	}
    }
  else if (A->value_type == COMPLEX)
    {
      struct index_complex_value_pair_t* _row = (struct index_complex_value_pair_t*) row;
      double_Complex* values = (double_Complex*) (A->values);

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  A->colidx[a+k] = _row[k].index;
	  values[a+k] = _row[k].value;
	}
    }
  else if (A->value_type == PATTERN)
    {
      int* _row = (int*) row;

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  A->colidx[a+k] = _row[k];
	}
    }
}


/**
 * Sorts the column indices within each column of the sparse matrix A.
 */
static void
csr_matrix_sort_colidx (struct csr_matrix_t* A)
{
  int i;
  int max_nnz;  /* Will be the max # of nonzeros in each column */

  bebop_log (2, "=== csr_matrix_sort_colidx ===\n");

#ifdef USING_VALGRIND_CLIENT_REQUESTS
  {
    int retval = 0;
    retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
    if (retval != 0)
      {
	bebop_log (0, "*** csr_matrix_sort_colidx: Pointer A is invalid! ***\n");
	bebop_exit (EXIT_FAILURE);
      }
    bebop_log (2, "Pointer A is valid.\n");
  }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

  if (A->m <= 0)
    {
      bebop_log (2, "m = %d <= 0\n", A->m);
      bebop_log (2, "=== Done with csr_matrix_sort_colidx ===\n");
      return;
    }

  /* Find the max # of nonzeros in each row.  This lets us allocate 
   * workspace large enough to accomodate any row (otherwise we would
   * have to allocate it separately for _each_ row, which would waste
   * a lot of malloc calls). */
  max_nnz = A->rowptr[1] - A->rowptr[0];
  for (i = 1; i < A->m; i++)
    {
      int nnz = A->rowptr[i+1] - A->rowptr[i];
      max_nnz = MAX(nnz, max_nnz);
    }

  /* Sort the column indices in each column, reordering the corresponding values accordingly. */
  if (A->value_type == REAL)
    {
      /* row: Workspace for sorting (index,value) pairs */
      struct index_real_value_pair_t* row = bebop_malloc (max_nnz * sizeof (struct index_real_value_pair_t));

      for (i = 0; i < A->m; i++)
	{
	  int nnz = A->rowptr[i+1] - A->rowptr[i];
	  /* Copy the (index,value) pairs into the temp location "row". */
	  copy_row2pairs (A, i, row, max_nnz);
	  /* Do the sorting in "row". */
	  qsort (row, nnz, sizeof (struct index_real_value_pair_t), compare_index_real_value_pairs);
	  /* Copy the (index,value) pairs back from "row" into the matrix. */
	  copy_pairs2row (row, max_nnz, A, i);
	}

      /* Free up the temp storage. */
      bebop_free (row);
    }
  else if (A->value_type == COMPLEX)
    {
      struct index_complex_value_pair_t* row = bebop_malloc (max_nnz * sizeof (struct index_complex_value_pair_t));
      for (i = 0; i < A->m; i++)
	{
	  int nnz = A->rowptr[i+1] - A->rowptr[i];
	  copy_row2pairs (A, i, row, max_nnz);
	  qsort (row, nnz, sizeof (struct index_real_value_pair_t), compare_index_complex_value_pairs);
	  copy_pairs2row (row, max_nnz, A, i);
	}

      bebop_free (row);
    }
  else if (A->value_type == PATTERN)
    {
      int* row = bebop_malloc (max_nnz * sizeof (int));

      for (i = 0; i < A->m; i++)
	{
	  int nnz = A->rowptr[i+1] - A->rowptr[i];
	  copy_row2pairs (A, i, row, max_nnz);
	  qsort (row, nnz, sizeof (int), compare_ints);
	  copy_pairs2row (row, max_nnz, A, i);
	}

      bebop_free (row);
    }
  bebop_log (2, "=== Done with csr_matrix_sort_colidx ===\n");
}




/*=====================================================================*/
void
pack_csr_matrix (struct csr_matrix_t* A, 
		 const int m, const int n, const int nnz, 
		 void* values, int* colidx, int* rowptr,
		 enum symmetry_type_t symmetry_type,
		 enum symmetric_storage_location_t symmetric_storage_location,
		 enum value_type_t value_type,
		 enum ownership_mode_t ownership,
		 void (*deallocator) (void*),
		 enum copy_mode_t copy_mode)
{
  bebop_log (3, "=== pack_csr_matrix ===\n");

#ifdef USING_VALGRIND_CLIENT_REQUESTS
  {
    int retval = 0;
    retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
    if (retval != 0)
      {
	bebop_log (0, "*** pack_csr_matrix: Pointer A is invalid! ***\n");
	bebop_exit (EXIT_FAILURE);
      }
  }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

  A->m = m;
  A->n = n;
  A->nnz = nnz;
  if (copy_mode == NO_COPY)
    {
      A->values = values;
      A->colidx = colidx;
      A->rowptr = rowptr;
    }
  else /* copy mode: this matrix gets (deep) copies of the user's input arrays */
    {
      /* In copy mode, the library is responsible for deallocating the arrays */
      if (ownership != LIBRARY_DEALLOCATES)
	{
	  bebop_log (0,  "*** WARNING: pack_csr_matrix: in copy mode, the library is responsible for deallocating the arrays: the ownership parameter has been changed to LIBRARY_DEALLOCATES ***\n");
	  A->ownership = LIBRARY_DEALLOCATES;
	}
      if (value_type == REAL) 
	{
	  A->values = bebop_malloc (nnz * sizeof (double));
	  memcpy (A->values, values, nnz * sizeof (double));
	}
      else if (value_type == COMPLEX)
	{
	  A->values = bebop_malloc (nnz * sizeof (double_Complex));
	  memcpy (A->values, values, nnz * sizeof (double_Complex));
	}
      else
	A->values = NULL;

      A->colidx = bebop_malloc (nnz * sizeof (int));
      memcpy (A->colidx, colidx, nnz * sizeof (int));

      A->rowptr = bebop_malloc ((n+1) * sizeof (int));
      memcpy (A->rowptr, rowptr, (n+1) * sizeof (int));
    }
  A->symmetry_type = symmetry_type;
  A->symmetric_storage_location = symmetric_storage_location;
  A->value_type = value_type;
  A->ownership = ownership;
  if (deallocator == NULL)
    A->deallocator = &free;
  else
    A->deallocator = deallocator;

  bebop_log (3, "=== Done with pack_csr_matrix ===\n");
}


/*=====================================================================*/
void
init_csr_matrix (struct csr_matrix_t* A, 
		 const int m, const int n, const int nnz, 
		 void* values, int* colidx, int* rowptr,
		 enum symmetry_type_t symmetry_type,
		 enum symmetric_storage_location_t symmetric_storage_location,
		 enum value_type_t value_type,
		 enum ownership_mode_t ownership,
		 void (*deallocator) (void*),
		 enum copy_mode_t copy_mode)
{
  pack_csr_matrix (A, m, n, nnz, values, colidx, rowptr, symmetry_type, 
		   symmetric_storage_location, value_type, ownership,
		   deallocator, copy_mode);
}



/*======================================================================*/
void
dealloc_csr_matrix (struct csr_matrix_t* A)
{
#ifdef USING_VALGRIND_CLIENT_REQUESTS
  int retval = 0;
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

  bebop_log (2, "=== dealloc_csr_matrix ===\n");

#ifdef USING_VALGRIND_CLIENT_REQUESTS
  {
    retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
    if (retval != 0)
      {
	bebop_log (0, "*** dealloc_csr_matrix: Pointer A is invalid! ***\n");
	bebop_exit (EXIT_FAILURE);
      }
  }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

  if (A->ownership == USER_DEALLOCATES)
    {
      if (bebop_debug_level () > 1)
	{
	  fprintf (stderr, "User responsible for deallocation; we do nothing\n");
	  fprintf (stderr, "=== Done with dealloc_csr_matrix ===\n");
	}
      return;
    }

  if (A->deallocator == NULL)
    A->deallocator = &free;

  if (A->values != NULL && A->nnz > 0) 
    {
#ifdef USING_VALGRIND_CLIENT_REQUESTS
      if (A->value_type != PATTERN && A->nnz > 0)
	{
	  bebop_log (2, "\tChecking if values pointer is valid\n");
	  if (A->value_type == REAL)
	    retval = VALGRIND_CHECK_READABLE( A->values, A->nnz * sizeof (double) );
	  else if (A->value_type == COMPLEX)
	    retval = VALGRIND_CHECK_READABLE( A->values, A->nnz * sizeof (double_Complex) );
	  else 
	    {
	      bebop_log (0, "*** WARNING: dealloc_csr_matrix: Invalid value type "
			"%d ***\n", A->value_type);
	      retval = 0; /* sensible default */
	    }

	  if (retval != 0)
	    {
	      bebop_log (0, "*** dealloc_csr_matrix: Valgrind says values pointer is invalid! ***\n");
	      bebop_exit (EXIT_FAILURE);
	    }
	}
#endif /* USING_VALGRIND_CLIENT_REQUESTS */
      bebop_log (2, "\tFreeing values pointer...\n");
      (A->deallocator) (A->values);
      bebop_log (2, "...done.\n");
      A->values = NULL;
    }

  if (A->colidx != NULL && A->nnz > 0)
    {
#ifdef USING_VALGRIND_CLIENT_REQUESTS
      bebop_log (2, "\tChecking if colidx pointer is valid\n");
      retval = VALGRIND_CHECK_READABLE( A->colidx, A->nnz * sizeof (int) );
      if (retval != 0)
	{
	  bebop_log (0, "*** dealloc_csr_matrix: Valgrind says colidx pointer is invalid! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

      bebop_log (2, "\tFreeing colidx pointer...\n");
      (A->deallocator) (A->colidx);
      bebop_log (2, "...done.\n");
      A->colidx = NULL;
    }

  if (A->rowptr != NULL && A->n > 0)
    {
#ifdef USING_VALGRIND_CLIENT_REQUESTS
      bebop_log (2, "\tChecking if rowptr pointer is valid\n");
      retval = VALGRIND_CHECK_READABLE( A->rowptr, (A->m + 1) * sizeof (int) );
      if (retval != 0)
	{
	  bebop_log (0, "*** dealloc_csr_matrix: Valgrind says rowptr pointer is invalid! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

      bebop_log (2, "\tFreeing rowptr pointer...\n");
      (A->deallocator) (A->rowptr);
      bebop_log (2, "...done.\n");
      A->rowptr = NULL;
    }

  A->m = 0;
  A->n = 0;
  A->nnz = 0;

  bebop_log (2, "=== Done with dealloc_csr_matrix ===\n");
}


struct csr_matrix_t*
clone_csr_matrix (const struct csr_matrix_t* src)
{
  struct csr_matrix_t* dest;
  bebop_log (3, "=== clone_csr_matrix ===\n");

#ifdef USING_VALGRIND_CLIENT_REQUESTS
  {
    int retval = 0;
    retval = VALGRIND_CHECK_READABLE( src, sizeof (struct csr_matrix_t) );
    if (retval != 0)
      {
	bebop_log (0, "*** clone_csr_matrix: pointer to input matrix is invalid! ***\n");
	bebop_exit (EXIT_FAILURE);
      }
  }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

  dest = bebop_calloc (1, sizeof (struct csr_matrix_t));
  copy_csr_matrix (dest, src);
  bebop_log (3, "=== Done with clone_csr_matrix ===\n");
  return dest;
}

/*======================================================================*/
void
copy_csr_matrix (struct csr_matrix_t* dest, const struct csr_matrix_t* src)
{
  const int m = src->m;
  const int n = src->n;
  const int nnz = src->nnz;

  bebop_log (2, "=== copy_csr_matrix ===\n");

#ifdef USING_VALGRIND_CLIENT_REQUESTS
  {
    int retval = 0;
    retval = VALGRIND_CHECK_READABLE( src, sizeof (struct csr_matrix_t) );
    if (retval != 0)
      {
	bebop_log (0, "*** copy_csr_matrix: pointer to src is invalid! ***\n");
	bebop_exit (EXIT_FAILURE);
      }
    retval = VALGRIND_CHECK_READABLE( dest, sizeof (struct csr_matrix_t) );
    if (retval != 0)
      {
	bebop_log (0, "*** copy_csr_matrix: pointer to dest is invalid! ***\n");
	bebop_exit (EXIT_FAILURE);
      }
  }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

  if (bebop_debug_level () > 1)
    {
      bebop_log (2, "Checking if input matrix is valid\n");
      if (! valid_csr_matrix_p (src))
	{
	  bebop_log (0, "*** copy_csr_matrix: input matrix is invalid! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
    }

  dest->m   = m;
  dest->n   = n;
  dest->nnz = nnz;

  if (src->value_type == REAL)
    {
      dest->values = bebop_malloc (nnz * sizeof (double));
      memcpy ((double*) (dest->values), (double*) (src->values), 
	      nnz * sizeof (double));
    }
  else if (src->value_type == COMPLEX)
    {
      dest->values = bebop_malloc (nnz * sizeof (double_Complex));
      memcpy ((double_Complex*) (dest->values), 
	      (double_Complex*) (src->values), 
	      nnz * sizeof (double_Complex));
    }
  else if (src->value_type == PATTERN)
    dest->values = NULL;

  dest->colidx = bebop_malloc (nnz * sizeof (int));
  dest->rowptr = bebop_malloc ((m+1) * sizeof (int));

  memcpy (dest->colidx, src->colidx, nnz * sizeof (int));
  memcpy (dest->rowptr, src->rowptr, (m+1) * sizeof (int));

  dest->symmetry_type = src->symmetry_type;
  dest->symmetric_storage_location = src->symmetric_storage_location;
  dest->value_type = src->value_type;
  dest->ownership = LIBRARY_DEALLOCATES; /* not the user's copy anymore */
  dest->deallocator = &free;

  if (bebop_debug_level () > 1)
    {
      if (! valid_csr_matrix_p (dest))
	{
	  bebop_log (0, "*** copy_csr_matrix: BUG: output matrix is invalid! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
    }

  bebop_log (2, "=== Done with copy_csr_matrix ===\n");
}


/*======================================================================*/
void
unpack_csr_matrix (const struct csr_matrix_t* A,
		   int* m, int* n, int* nnz,
		   void** values, int** colidx, int** rowptr,
		   enum symmetry_type_t* symmetry_type,
		   enum symmetric_storage_location_t* symmetric_storage_location,
		   enum value_type_t* value_type)
{
  bebop_log (3, "=== unpack_csr_matrix ===\n");

#ifdef USING_VALGRIND_CLIENT_REQUESTS
  {
    int retval = 0;
    retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
    if (retval != 0)
      {
	bebop_log (0, "*** unpack_csr_matrix: pointer to A is invalid! ***\n");
	bebop_exit (EXIT_FAILURE);
      }
  }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

  *m = A->m;
  *n = A->n;
  *nnz = A->nnz;

  *values = A->values;
  *colidx = A->colidx;
  *rowptr = A->rowptr;

  *symmetry_type = A->symmetry_type;
  *symmetric_storage_location = A->symmetric_storage_location;
  *value_type = A->value_type;

  bebop_log (3, "=== Done with unpack_csr_matrix ===\n");
}



/*======================================================================*/
struct csr_matrix_t*
create_csr_matrix (const int m, const int n, const int nnz, 
		   void* values, int* colidx, int* rowptr,
		   enum symmetry_type_t symmetry_type,
		   enum symmetric_storage_location_t symmetric_storage_location,
		   enum value_type_t value_type,
		   enum ownership_mode_t ownership,
		   void (*deallocator) (void*),
		   enum copy_mode_t copy_mode)
{
  struct csr_matrix_t *A;

  bebop_log (3, "=== create_csr_matrix ===\n");
  A = bebop_malloc (sizeof (struct csr_matrix_t));
  init_csr_matrix (A, m, n, nnz, values, colidx, rowptr, symmetry_type, 
		   symmetric_storage_location, value_type, ownership, 
		   deallocator, copy_mode);
  bebop_log (3, "=== Done with create_csr_matrix ===\n");
  return A;
}


/*======================================================================*/
void
destroy_csr_matrix (struct csr_matrix_t* A)
{
  bebop_log (2, "=== destroy_csr_matrix ===\n"); 
  if (A != NULL)
    {
#ifdef USING_VALGRIND_CLIENT_REQUESTS
      int retval = 0;
      retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
      if (retval != 0)
	{
	  bebop_log (0, "*** destroy_csr_matrix: Valgrind says that the "
		    "struct A itself is an invalid pointer! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
#endif /* USING_VALGRIND_CLIENT_REQUESTS */
      dealloc_csr_matrix (A);
      /* 
       * NOTE: We don't use the deallocator in A to free the struct
       * itself; that struct was allocated using bebop_malloc. 
       */
      bebop_free (A);
    }
  bebop_log (2, "=== Done with destroy_csr_matrix ===\n"); 
}

struct csc_matrix_t*
csr_to_csc (const struct csr_matrix_t* A)
{
  int *col_nnz = NULL; /* # of nonzeros in each column */
  int i;
  struct csc_matrix_t* B = bebop_calloc (1, sizeof (struct csc_matrix_t));

  bebop_log (2, "=== csr_to_csc ===\n");

  if (bebop_debug_level () > 1)
    {
#ifdef USING_VALGRIND_CLIENT_REQUESTS
      int retval = 0;
      retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
      if (retval != 0)
	{
	  bebop_log (0, "*** csr_to_csc: Valgrind says that the "
		    "struct A itself is an invalid pointer! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
#endif /* USING_VALGRIND_CLIENT_REQUESTS */
      if (! valid_csr_matrix_p (A))
	{
	  bebop_log (0, "*** csr_to_csc: the input matrix A is invalid! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
    }

  B->m = A->m;
  B->n = A->n;
  B->nnz = A->nnz;
  B->symmetry_type = A->symmetry_type;
  B->symmetric_storage_location = A->symmetric_storage_location;
  B->value_type = A->value_type;

  B->colptr = bebop_malloc ((A->n + 1) * sizeof (int));
  B->rowidx = bebop_malloc (A->nnz * sizeof (int));
  if (A->value_type == REAL)
    B->values = bebop_malloc (A->nnz * sizeof (double));
  else if (A->value_type == COMPLEX)
    B->values = bebop_malloc (A->nnz * sizeof (double_Complex));
  else if (A->value_type == PATTERN)
    B->values = NULL;

  B->ownership = LIBRARY_DEALLOCATES;
  B->deallocator = &free;

  if (A->symmetry_type == SYMMETRIC)
    {
      /* Just give it the transpose, since $A^T = A$ for a symmetric matrix. 
       * TODO:  Think about a clever way to do this for a skew-symmetric or
       * a Hermitian matrix!  */

      assert (A->m == A->n);

      memcpy (B->colptr, A->rowptr, (A->n + 1) * sizeof (int));
      memcpy (B->rowidx, A->colidx, A->nnz * sizeof (int));

      if (A->value_type == REAL)
	memcpy ((double*) (B->values), (double*) (A->values), 
		A->nnz * sizeof (double));
      else if (A->value_type == COMPLEX)
	memcpy ((double_Complex*) (B->values), 
		(double_Complex*) (A->values), 
		A->nnz * sizeof (double));

      bebop_log (2, "=== Done with csr_to_csc ===\n");
      return B;
    }

  /* count # of non-zeros per column */

  col_nnz = bebop_calloc (A->n, sizeof(int));

  for (i = 0; i < A->nnz; i++)
    {
      int k = A->colidx[i];
      col_nnz[k]++;
    }

  /*
   *  initialize CSC's column pointers.
   *  reset col_nnz to zero.
   */
  B->colptr[0] = 0;
  for (i = 1; i <= A->n; i++)
    {
      B->colptr[i] = B->colptr[i-1] + col_nnz[i-1];
      col_nnz[i-1] = 0;
    }

  /*
   *  convert from CSR to CSC.
   *  use col_nnz to keep track of the number of non-zeros
   *  added to each column.
   */
  if (A->value_type == REAL)
    {
      double* in_values  = (double*) (A->values);
      double* out_values = (double*) (B->values);

      for (i = 0; i < A->m; i++)
	{
	  int k;
	  int nnz_row;    /* # of non-zeros in row i of A */

	  nnz_row = A->rowptr[i+1] - A->rowptr[i];

	  for (k = A->rowptr[i]; k < (A->rowptr[i]+nnz_row); k++)
	    {
	      int j = A->colidx[ k ];             /* col index */
	      double a = in_values[ k ];          /* non-zero value */
	      int h = B->colptr[j] + col_nnz[j];  /* non-zero position */

	      /* add the non-zero A(i,j) to B */

	      B->rowidx[ h ] = i;
	      out_values[ h ] = a;

	      col_nnz[j]++;
	    }
	}
    }
  else if (A->value_type == COMPLEX)
    {
      double_Complex* in_values  = (double_Complex*) (A->values);
      double_Complex* out_values = (double_Complex*) (B->values);

      for (i = 0; i < A->m; i++)
	{
	  int k;
	  int nnz_row;    /* # of non-zeros in row i of A */

	  nnz_row = A->rowptr[i+1] - A->rowptr[i];

	  for (k = A->rowptr[i]; k < (A->rowptr[i]+nnz_row); k++)
	    {
	      int j = A->colidx[ k ];               /* col index */
	      double_Complex a = in_values[ k ];   /* non-zero value */
	      int h = B->colptr[j] + col_nnz[j];    /* non-zero position */

	      /* add the non-zero A(i,j) to B */

	      B->rowidx[ h ] = i;
	      out_values[ h ] = a;

	      col_nnz[j]++;
	    }
	}
    }
  else if (A->value_type == PATTERN)
    {
      for (i = 0; i < A->m; i++)
	{
	  int k;
	  int nnz_row;    /* # of non-zeros in row i of A */

	  nnz_row = A->rowptr[i+1] - A->rowptr[i];

	  for (k = A->rowptr[i]; k < (A->rowptr[i]+nnz_row); k++)
	    {
	      int j = A->colidx[ k ];             /* col index */
	      int h = B->colptr[j] + col_nnz[j];  /* non-zero position */

	      /* add the non-zero A(i,j) to B */
	      B->rowidx[ h ] = i;
	      col_nnz[j]++;
	    }
	}
    }

  /* clean-up */
  bebop_free (col_nnz);

  bebop_log (2, "=== Done with csr_to_csc ===\n");
  return B;
}


struct csr_matrix_t*
csc_to_csr (struct csc_matrix_t* A)
{
  int *row_nnz = NULL; /* # of nonzeros in each row */
  int i, j;
  struct csr_matrix_t* B = bebop_calloc (1, sizeof (struct csr_matrix_t));

  bebop_log (2, "=== csc_to_csr ===\n");
  B->m = A->m;
  B->n = A->n;
  B->nnz = A->nnz;
  B->symmetry_type = A->symmetry_type;
  B->symmetric_storage_location = A->symmetric_storage_location;
  B->value_type = A->value_type;

  B->rowptr = bebop_malloc ((A->m + 1) * sizeof (int));
  B->colidx = bebop_malloc (A->nnz * sizeof (int));
  if (A->value_type == REAL)
    B->values = bebop_malloc (A->nnz * sizeof (double));
  else if (A->value_type == COMPLEX)
    B->values = bebop_malloc (A->nnz * sizeof (double_Complex));
  else if (A->value_type == PATTERN)
    B->values = NULL;

  B->ownership = LIBRARY_DEALLOCATES;
  B->deallocator = &free;

  if (A->symmetry_type == SYMMETRIC)
    {
      /* Just give it the transpose, since $A^T = A$ for a symmetric matrix. 
       * TODO:  Think about a clever way to do this for a skew-symmetric or
       * a Hermitian matrix!  */

      assert (A->m == A->n);

      memcpy (B->rowptr, A->colptr, (A->m + 1) * sizeof (int));
      memcpy (B->colidx, A->rowidx, A->nnz * sizeof (int));

      if (A->value_type == REAL)
	memcpy ((double*) (B->values), (double*) (A->values), 
		A->nnz * sizeof (double));
      else if (A->value_type == COMPLEX)
	memcpy ((double_Complex*) (B->values), 
		(double_Complex*) (A->values), 
		A->nnz * sizeof (double));

      return B;
    }

  /* count # of non-zeros per row */

  row_nnz = bebop_calloc (A->m, sizeof(int));

  for (i = 0; i < A->nnz; i++)
    {
      int k = A->rowidx[i];
      row_nnz[k]++;
    }

  /*
   *  initialize CSR's row pointers.
   *  reset row_nnz to zero.
   */
  B->rowptr[0] = 0;
  for (i = 1; i <= A->m; i++)
    {
      B->rowptr[i] = B->rowptr[i-1] + row_nnz[i-1];
      row_nnz[i-1] = 0;
    }

  /*
   *  convert from CSC to CSR.
   *  use row_nnz to keep track of the number of non-zeros
   *  added to each row.
   */
  if (A->value_type == REAL)
    {
      double* in_values  = (double*) (A->values);
      double* out_values = (double*) (B->values);

      for (j = 0; j < A->n; j++)
	{
	  int k;
	  int nnz_col;    /* # of non-zeros in column j of A */

	  nnz_col = A->colptr[j+1] - A->colptr[j];

	  for (k = A->colptr[j]; k < (A->colptr[j]+nnz_col); k++)
	    {
	      int i = A->rowidx[ k ];             /* row index */
	      double a = in_values[ k ];          /* non-zero value */
	      int h = B->rowptr[i] + row_nnz[i];  /* non-zero position */

	      /* add the non-zero A(i,j) to B */

	      B->colidx[ h ] = j;
	      out_values[ h ] = a;

	      row_nnz[i]++;
	    }
	}
    }
  else if (A->value_type == COMPLEX)
    {
      double_Complex* in_values  = (double_Complex*) (A->values);
      double_Complex* out_values = (double_Complex*) (B->values);

      for (j = 0; j < A->n; j++)
	{
	  int k;
	  int nnz_col;    /* # of non-zeros in column j of A */

	  nnz_col = A->colptr[j+1] - A->colptr[j];

	  for (k = A->colptr[j]; k < (A->colptr[j]+nnz_col); k++)
	    {
	      int i = A->rowidx[ k ];               /* row index */
	      double_Complex a = in_values[ k ];   /* non-zero value */
	      int h = B->rowptr[i] + row_nnz[i];    /* non-zero position */

	      /* add the non-zero A(i,j) to B */

	      B->colidx[ h ] = j;
	      out_values[ h ] = a;

	      row_nnz[i]++;
	    }
	}
    }
  else if (A->value_type == PATTERN)
    {
      for (j = 0; j < A->n; j++)
	{
	  int k;
	  int nnz_col;    /* # of non-zeros in column j of A */

	  nnz_col = A->colptr[j+1] - A->colptr[j];

	  for (k = A->colptr[j]; k < (A->colptr[j]+nnz_col); k++)
	    {
	      int i = A->rowidx[ k ];  /* row index */
	      int h = B->rowptr[i] + row_nnz[i];  /* non-zero position */

	      /* add the non-zero A(i,j) to B */
	      B->colidx[ h ] = j;
	      row_nnz[i]++;
	    }
	}
    }

  /* clean-up */
  bebop_free (row_nnz);

  bebop_log (2, "=== Done with csc_to_csr ===\n");
  return B;
}


int
save_csr_matrix_in_harwell_boeing_format (const char* const filename, 
					  const struct csr_matrix_t* A)
{
  struct csc_matrix_t *B = NULL;

  bebop_log (2, "=== save_csr_matrix_in_harwell_boeing_format ===\n");
  if (bebop_debug_level () > 1)
    {
#ifdef USING_VALGRIND_CLIENT_REQUESTS
      {
	int retval = 0;
	retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
	if (retval != 0)
	  {
	    bebop_log (0, "*** save_csr_matrix_in_harwell_boeing_format: "
		      "pointer to A is invalid! ***\n");
	    bebop_exit (EXIT_FAILURE);
	  }
      }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

      if (! valid_csr_matrix_p (A))
	{
	  bebop_log (0, "*** save_csr_matrix_in_harwell_boeing_format: A "
		    "is an invalid CSR matrix! ***\n");
	  return -1;
	}

    }
  
  B = csr_to_csc (A);
  if (B == NULL)
    {
      bebop_log (2, "*** save_csr_matrix_in_harwell_boeing_format: "
		"Conversion from CSC to CSR failed! ***\n");
      return -1;
    }
  else
    {
      const int errcode = save_csc_matrix_in_harwell_boeing_format (filename, B);
      destroy_csc_matrix (B);
      bebop_log (2, "=== Done with save_csr_matrix_in_harwell_boeing_format ===\n");
      return errcode;
    }
}


int
print_csr_matrix_in_matrix_market_format (FILE* out, const struct csr_matrix_t* A)
{
  int start, end;
  int i, k;
  char value_type_label[20];
  char symmetry_type_label[20];

  bebop_log (2, "=== print_csr_matrix_in_matrix_market_format ===\n"); 
  if (bebop_debug_level () > 1)
    {
      if (! valid_csr_matrix_p (A))
	{
	  bebop_log (0, "*** print_csr_matrix_in_matrix_market_format: A is an "
		    "invalid CSR matrix! ***\n");
	  return -1;
	}
    }

  if (A->value_type == REAL)
    strncpy (value_type_label, "real", 19);
  else if (A->value_type == COMPLEX)
    strncpy (value_type_label, "complex", 19);
  else if (A->value_type == PATTERN)
    strncpy (value_type_label, "pattern", 19);
  else 
    {
      bebop_log (0, "*** print_csr_matrix_in_matrix_market_format: Unsupported "
		"value type %d; supported value types are %d, %d and %d ***\n", 
		A->value_type, REAL, COMPLEX, PATTERN);
      bebop_log (2, "=== Done with print_csr_matrix_in_matrix_market_format ===\n"); 
      return -1;
    }

  if (A->symmetry_type == UNSYMMETRIC)
    strncpy (symmetry_type_label, "general", 19);
  else if (A->symmetry_type == SYMMETRIC)
    strncpy (symmetry_type_label, "symmetric", 19);
  else if (A->symmetry_type == SKEW_SYMMETRIC)
    strncpy (symmetry_type_label, "skew-symmetric", 19);
  else if (A->symmetry_type == HERMITIAN)
    strncpy (symmetry_type_label, "hermitian", 19);
  else
    {
      bebop_log (0, "*** print_csr_matrix_in_matrix_market_format: Unsuppo"
		"rted symmetry type %d; supported symmetry types are %d"
		", %d, %d and %d ***\n", A->symmetry_type, UNSYMMETRIC, 
		SYMMETRIC, SKEW_SYMMETRIC, HERMITIAN);
      bebop_log (2, "=== Done with print_csr_matrix_in_matrix_market_format ===\n"); 
      return -1;
    }

  fprintf (out, "%%%%MatrixMarket matrix coordinate %s %s\n", 
	   value_type_label, symmetry_type_label);

  if (bebop_debug_level() > 1)
    {
      fprintf (out, "%% rowptr[%d]: ", A->m + 1);
      for (i = 0; i <= A->m; i++)
	fprintf (out, " %d", A->rowptr[i]);

      fprintf (out, "\n%% colidx[%d]: ", A->nnz);
      for (i = 0; i < A->nnz; i++)
	fprintf (out, " %d", A->colidx[i]);

      if (A->value_type != PATTERN)
	{
	  fprintf (out, "\n%% values[%d]: ", A->nnz);

	  if (A->value_type == REAL)
	    {
	      const double* const values = (const double* const) (A->values);

	      for (i = 0; i < A->nnz; i++)
		fprintf (out, " %g", values[i]);
	    }
	  else if (A->value_type == COMPLEX)
	    {
	      const double_Complex* const values = (const double_Complex* const) (A->values);

	      for (i = 0; i < A->nnz; i++)
		fprintf (out, " %g+I*%g", double_Complex_real_part(values[i]), double_Complex_imag_part(values[i]));
	    }
	}

      fprintf (out, "\n");
    }

  fprintf (out, "%d %d %d\n", A->m, A->n, A->nnz);

  if (A->value_type == REAL)
    {
      const double* const values = (const double* const) (A->values);

      for (i = 0; i < A->m; i++)
	{
	  start = A->rowptr[i];
	  end   = A->rowptr[i+1];

	  /* MatrixMarket files use 1-based indices. */
	  for (k = start; k < end; k++)
	    fprintf (out, "%d %d %.13e\n", i+1, A->colidx[k] + 1, values[k]);
	}
    }
  else if (A->value_type == COMPLEX)
    {
      const double_Complex* const values = (const double_Complex* const) (A->values);

      for (i = 0; i < A->m; i++)
	{
	  start = A->rowptr[i];
	  end   = A->rowptr[i+1];

	  /* MatrixMarket files use 1-based indices. */
	  for (k = start; k < end; k++)
	    fprintf (out, "%d %d %.13e %.13e\n", i+1, A->colidx[k] + 1, 
		     double_Complex_real_part(values[k]), 
		     double_Complex_imag_part(values[k]));
	}
    }
  else if (A->value_type == PATTERN)
    {
      for (i = 0; i < A->m; i++)
	{
	  start = A->rowptr[i];
	  end   = A->rowptr[i+1];

	  /* MatrixMarket files use 1-based indices. */
	  for (k = start; k < end; k++)
	    fprintf (out, "%d %d\n", i+1, A->colidx[k] + 1);
	}
    }

  bebop_log (2, "=== Done with print_csc_matrix_in_matrix_market_format ===\n");
  return 0;
}



int 
save_csr_matrix_in_matrix_market_format (const char* const filename, 
					 const struct csr_matrix_t* A)
{
  FILE* out = NULL;
  int errcode = 0;

  bebop_log (2, "=== save_csr_matrix_in_matrix_market_format ===\n");

#ifdef USING_VALGRIND_CLIENT_REQUESTS
  {
    int retval = 0;
    retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
    if (retval != 0)
      {
	bebop_log (0, "*** save_csr_matrix_in_matrix_market_format: "
		  "Pointer A is invalid! ***\n");
	bebop_exit (EXIT_FAILURE);
      }
  }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

  out = fopen (filename, "w");
  if (out == NULL)
    {
      fprintf (stderr, "*** save_csr_matrix_in_matrix_market_format: failed "
	       "to open output file \"%s\" ***\n", filename);
      return -1;
    }

  errcode = print_csr_matrix_in_matrix_market_format (out, A);
  if (0 != fclose (out))
    {
      fprintf (stderr, "*** save_csr_matrix_in_matrix_market_format: failed "
	       "to close output file \"%s\" ***\n", filename);
      return -1;
    }
  bebop_log (2, "=== Done with save_csr_matrix_in_matrix_market_format ===\n");
  return errcode;
}


int
csr_matrix_expand_symmetric_storage (struct csr_matrix_t* A)
{
  /* # of non-zeros in each row of original symmetric matrix */
  int* cur_row_nnz;   
  /* # of non-zeros in each row of final expanded matrix */
  int* new_row_nnz;   
  /* total # of non-zeros in final expanded matrix */
  int new_nnz;        
  int i;
  int* new_rowptr = NULL;
  int* new_colidx = NULL;
  void* new_values = NULL;

  bebop_log (2, "=== csr_matrix_expand_symmetric_storage ===\n");
  if (bebop_debug_level () > 1)
    {
#ifdef USING_VALGRIND_CLIENT_REQUESTS
      {
	int retval = 0;
	retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
	if (retval != 0)
	  {
	    bebop_log (0, "*** csr_matrix_expand_symmetric_storage: "
		      "Pointer A is invalid! ***\n");
	    bebop_exit (EXIT_FAILURE);
	  }
      }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

      if (! valid_csr_matrix_p (A))
	{
	  bebop_log (0, "*** csr_matrix_expand_symmetric_storage: "
		    "input matrix is invalid! ***\n");
	  return -1;
	}
    }

  if (A->m != A->n)
    {
      fprintf (stderr, "*** csr_matrix_expand_symmetric_storage: A is not "
	       "square! ***\n");
      return -1;
    }
  else if (A->symmetry_type == UNSYMMETRIC)
    {
      bebop_log (1, "*** csr_matrix_expand_symmetric_storage: A is "
		"already stored in an unsymmetric format, so we don\'t "
		"need to do anything ***\n");
      /* Not an error */
      return 0; 
    }

  /* set-up */
  cur_row_nnz = bebop_calloc (A->m, sizeof (int));
  new_row_nnz = bebop_calloc (A->m, sizeof(int));

  /*
   * Scan A and count how many new non-zeros we'll need to create.
   *
   * Post:
   *   cur_row_nnz[i] == # of non-zeros in row i of the original symmetric 
   *                     matrix.
   *   new_row_nnz[i] == # of non-zeros to be stored in row i of the final 
   *                     expanded matrix.
   *   new_nnz == total # of non-zeros to be stored in the final expanded 
   *              matrix.
   */
  new_nnz = 0;
  for (i = 0; i < A->m; i++)
    {
      cur_row_nnz[i] = A->rowptr[i+1] - A->rowptr[i];
      new_row_nnz[i] = cur_row_nnz[i];
      new_nnz += new_row_nnz[i];
    }
  for (i = 0; i < A->m; i++)
    {
      int k;
      for (k = A->rowptr[i]; k < A->rowptr[i+1]; k++)
	{
	  int j = A->colidx[k];
	  if (j != i)
	    {
	      new_row_nnz[j]++;
	      new_nnz++;
	    }
	}
    }

  /*
   *  Initialize row pointers in expanded matrix.
   *
   *  Post:
   *    new_rowptr initialized to the correct, final values.
   *    new_row_nnz[i] reset to be equal to cur_row_nnz[i].
   */
  new_rowptr = bebop_calloc (A->m + 1, sizeof(int));
  for (i = 1; i <= A->m; i++)
    {
      new_rowptr[i] = new_rowptr[i-1] + new_row_nnz[i-1];
      new_row_nnz[i-1] = cur_row_nnz[i-1];
    }
  new_rowptr[A->m] = new_nnz;

  new_colidx = bebop_calloc (new_nnz, sizeof (int));
  if (A->value_type == REAL)
    new_values = bebop_malloc (new_nnz * sizeof (double));
  else if (A->value_type == COMPLEX)
    new_values = bebop_malloc (new_nnz * sizeof (double_Complex));
  else if (A->value_type == PATTERN)
    new_values = NULL;
 
  /*
   *  Complete expansion of A to full storage.
   *
   *  Post:
   *    (new_rowptr, new_colidx, new_values) is the full-storage equivalent of A.
   *    new_row_nnz[i] == # of non-zeros in row i of A.
   */
  for (i = 0; i < A->m; i++)
    {
      int cur_nnz = cur_row_nnz[i];
      int k_cur   = A->rowptr[i];
      int k_new   = new_rowptr[i];

      /* copy current non-zeros from old matrix to new matrix */
      memcpy (new_colidx + k_new, A->colidx + k_cur, cur_nnz * sizeof(int));
      /* mfh 24 Feb 2006:  fixed bug (had forgotten to check for value type) */
      if (A->value_type == REAL)
	memcpy ((double*)new_values + k_new, 
		(double*)new_values + k_cur, 
		cur_nnz * sizeof(double));
      else if (A->value_type == COMPLEX)
	memcpy ((double_Complex*)new_values + k_new, 
		(double_Complex*)new_values + k_cur, 
		cur_nnz * sizeof(double_Complex));

      /* fill in the symmetric "missing" values */
      while (k_cur < A->rowptr[i+1])
	{
	  /* non-zero of original matrix */
	  int j = A->colidx[k_cur];

	  if (j != i)  /* if not a non-diagonal element ... */
	    {
	      /* position of this transposed element in new matrix */
	      k_new = new_rowptr[j] + new_row_nnz[j];

	      /* store */
	      new_colidx[k_new] = i;
	      if (A->value_type == REAL)
		{
		  double* old_values = (double*) (A->values);
		  double* __new_values = (double*) new_values;
		  double a = old_values[k_cur];

		  if (A->symmetry_type == SYMMETRIC)
		    __new_values[k_new] = a;
		  else if (A->symmetry_type == SKEW_SYMMETRIC)
		    __new_values[k_new] = -a;
		  else if (A->symmetry_type == HERMITIAN)
		    __new_values[k_new] = a;
		}
	      else if (A->value_type == COMPLEX)
		{
		  double_Complex* old_values = (double_Complex*) (A->values);
		  double_Complex* __new_values = (double_Complex*) new_values;
		  double_Complex a = old_values[k_cur];

		  if (A->symmetry_type == SYMMETRIC)
		    __new_values[k_new] = a;
		  else if (A->symmetry_type == SKEW_SYMMETRIC)
		    __new_values[k_new] = double_Complex_negate (a);
		  else if (A->symmetry_type == HERMITIAN)
		    __new_values[k_new] = double_Complex_conj (a);
		}
		    
	      /*  update so next element stored at row j will appear
	       *  at the right place.  */
	      new_row_nnz[j]++;
	    }

	  k_cur++;
	}
    }

  /* clean-up */
  free (cur_row_nnz);
  free (new_row_nnz);

  A->nnz = new_nnz;
  bebop_free (A->values);
  A->values = new_values;
  bebop_free (A->colidx);
  A->colidx = new_colidx;
  bebop_free (A->rowptr);
  A->rowptr = new_rowptr;

  csr_matrix_sort_colidx (A);
  /* 
   * Since the library reallocated the arrays, it's now responsible for 
   * deallocating them.
   */
  A->ownership = LIBRARY_DEALLOCATES;
  A->deallocator = &free;
  return 0;
}


int
print_csr_matrix_in_matlab_format (FILE* out, struct csr_matrix_t* A)
{
  int errcode = 0;
  struct coo_matrix_t* B = csr_to_coo (A);
  if (B == NULL)
    {
      bebop_log (1, "*** print_csr_matrix_in_matlab_format: Failed to "
		"convert CSR-format sparse matrix into COO format for "
		"printing! ***\n");
      return -1;
    }
  errcode = print_coo_matrix_in_matlab_format (out, B);
  destroy_coo_matrix (B);
  return errcode;
}


int
save_csr_matrix_in_matlab_format (const char* const filename, struct csr_matrix_t* A)
{
  FILE* out = NULL;
  int errcode = 0;

  bebop_log (2, "=== save_csr_matrix_in_matlab_format ===\n");
  out = fopen (filename, "w");
  if (out == NULL)
    {
      fprintf (stderr, "*** save_csr_matrix_in_matlab_format: failed "
	       "to open output file \"%s\" ***\n", filename);
      bebop_log (2, "=== Done with save_csr_matrix_in_matlab_format ===\n");
      return -1;
    }

  errcode = print_csr_matrix_in_matlab_format (out, A);
  if (0 != fclose (out))
    {
      fprintf (stderr, "*** save_csr_matrix_in_matlab_format: failed "
	       "to close output file \"%s\" ***\n", filename);
      bebop_log (2, "=== Done with save_csr_matrix_in_matlab_format ===\n");
      return -1;
    }
  bebop_log (2, "=== Done with save_csr_matrix_in_matlab_format ===\n");
  return errcode;
}


int
valid_csr_matrix_p (const struct csr_matrix_t* A)
{
  int m;
  int n;
  int nnz;
  int* rowptr;
  int* colidx;
  int reached_nnz_in_rowptr_p = 0;
  int i, j;

  bebop_log (2, "=== valid_csr_matrix_p ===\n"); 

#ifdef USING_VALGRIND_CLIENT_REQUESTS
      {
	int retval = 0;
	retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
	if (retval != 0)
	  {
	    bebop_log (0, "*** valid_csr_matrix_p: Pointer A "
		      "is invalid! ***\n");
	    bebop_exit (EXIT_FAILURE);
	  }
      }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

  m = A->m;
  n = A->n;
  nnz = A->nnz;
  colidx = (int*) (A->colidx);
  rowptr = (int*) (A->rowptr);

  if (m < 1)
    {
      bebop_log (1, "*** valid_csr_matrix_p: Matrix has m = %d < 1 ***\n", m);
      bebop_log (2, "=== Done with valid_csr_matrix_p ===\n"); 
      return 0;
    }
  else if (n < 1)
    {
      bebop_log (1, "*** valid_csr_matrix_p: Matrix has n = %d < 1 ***\n", nnz);
      bebop_log (2, "=== Done with valid_csr_matrix_p ===\n"); 
      return 0;
    }
  else if (nnz < 0)
    {
      bebop_log (1, "*** valid_csr_matrix_p: Matrix has nnz = %d < 0 ***\n", nnz);
      bebop_log (2, "=== Done with valid_csr_matrix_p ===\n"); 
      return 0;
    }

#ifdef USING_VALGRIND_CLIENT_REQUESTS
 {
   int retval = 0;
   /* See if rowptr has m+1 entries */
   retval = VALGRIND_CHECK_READABLE( A->rowptr, (A->m + 1) * sizeof (int) );
   if (retval != 0)
     {
       bebop_log (1, "*** valid_csr_matrix_p: Valgrind says A->rowptr does not "
		 "have m+1 = %d entries! ***\n", A->m + 1);
       return 0;
     }
   /* See if colidx has nnz entries */
   retval = VALGRIND_CHECK_READABLE( A->colidx, (A->nnz) * sizeof (int) );
   if (retval != 0)
     {
       bebop_log (1, "*** valid_csr_matrix_p: Valgrind says A->colidx does not "
		 "have nnz = %d entries! ***\n", A->nnz);
       return 0;
     }
   /* See if values has nnz entries */
   if (A->value_type == REAL)
     {
       retval = VALGRIND_CHECK_READABLE( A->values, (A->nnz) * sizeof (double) );
       if (retval != 0)
	 {
	   bebop_log (1, "*** valid_csr_matrix_p: Valgrind says A->values does not "
		     "have nnz = %d entries! ***\n", A->nnz);
	   return 0;
	 }
     }
   else if (A->value_type == COMPLEX)
     {
       retval = VALGRIND_CHECK_READABLE( A->values, (A->nnz) * sizeof (double_Complex) );
       if (retval != 0)
	 {
	   bebop_log (1, "*** valid_csr_matrix_p: Valgrind says A->values does not "
		     "have nnz = %d entries! ***\n", A->nnz);
	   return 0;
	 }
     } 
 }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

  for (i = 0; i < m; i++)
    {
      const int cur = rowptr[i];
      const int next = rowptr[i+1];

      if (cur == nnz)
	{
	  if (! reached_nnz_in_rowptr_p)
	    reached_nnz_in_rowptr_p = 1;
	  if (next != cur)
	    {
	      bebop_log (1, "*** valid_csr_matrix_p: Although rowptr[%d]==nnz==%d, "
			"rowptr[%d+1]==%d != rowptr[%d] ***\n", 
			i, nnz, i, next, i);
	      bebop_log (2, "=== Done with valid_csr_matrix_p ===\n"); 
	      return 0;
	    }
	}

      /* Verify that rowptr entries are in order */
      if (cur > next)
	{
	  bebop_log (1, "*** valid_csr_matrix_p: Matrix: rowptr[%d] = %d > "
		    "rowptr[%d] = %d ***\n", i, cur, i+1, next);
	  bebop_log (2, "=== Done with valid_csr_matrix_p ===\n"); 
	  return 0;
	}

      /* Verify that current rowptr entry doesn't exceed nnz */
      if (cur >= nnz) /** \bug Should be cur >= nnz, right???? */
	{
	  bebop_log (1, "*** valid_csr_matrix_p: Matrix: At col %d, rowptr[i] = "
		    "%d >= nnz = %d ***\n", i, rowptr[i], nnz);
	  bebop_log (2, "=== Done with valid_csr_matrix_p ===\n"); 
	  return 0;
	}

      /* Verify that colidx entries are in range */
      for (j = cur; j < next; j++)
	{
	  if (colidx[j] < 0 || colidx[j] >= n)
	    {
	      bebop_log (1, "*** valid_csr_matrix_p: Matrix: at col %d, "
		       "colidx[%d]=%d out of range [%d,%d) ***\n", 
		       i, j, colidx[j], 0, n);
	      bebop_log (2, "=== Done with valid_csr_matrix_p ===\n"); 
	      return 0;
	    }
	}
    }

  /* Verify that the last rowptr entry is nnz */
  if (rowptr[m] != nnz)
    {
      bebop_log (1, "*** valid_csr_matrix_p: rowptr[m=%d] = %d != nnz=%d ***\n", 
		m, rowptr[m], nnz);
      bebop_log (2, "=== Done with valid_csr_matrix_p ===\n"); 
      return 0;
    }

  bebop_log (1, "valid_csr_matrix_p: matrix is valid\n");
  bebop_log (2, "=== Done with valid_csr_matrix_p ===\n");
  return 1;
}


struct csr_matrix_t*
csr_matrix_matmatmult (struct csr_matrix_t* B, struct csr_matrix_t* A)
{
  /* B is m x p and A is p x n. */
  const int m = B->m;
  const int p = B->n;
  const int n = A->n;
  const enum value_type_t value_type = B->value_type;
  int errcode = 0;

  bebop_log (2, "=== csr_matrix_matmatmult ===\n");
  if (bebop_debug_level () > 1)
    {
      int failure = 0;

#ifdef USING_VALGRIND_CLIENT_REQUESTS
      {
	int retval = 0;
	retval = VALGRIND_CHECK_READABLE( B, sizeof (struct csr_matrix_t) );
	if (retval != 0)
	  {
	    bebop_log (0, "*** valid_csr_matrix_p: Pointer B (first input "
		      "matrix) is invalid! ***\n");
	    bebop_exit (EXIT_FAILURE);
	  }
	retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
	if (retval != 0)
	  {
	    bebop_log (0, "*** valid_csr_matrix_p: Pointer A (second "
		      "input matrix) is invalid! ***\n");
	    bebop_exit (EXIT_FAILURE);
	  }
      }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

      if (! valid_csr_matrix_p (B))
	{
	  bebop_log (0, "*** csr_matrix_matmatmult: first input matrix is invalid! ***\n");
	  failure = 1;
	}
      if (! valid_csr_matrix_p (A))
	{
	  bebop_log (0, "*** csr_matrix_matmatmult: second input matrix is invalid! ***\n");
	  failure = 1;
	}
      if (failure)
	return NULL;
    }

  if (p != A->m)
    {
      bebop_log (0, "*** csr_matrix_matmatmult: incompatible dimensions for "
		"B*A: B is %d by %d and A is %d by %d ***\n", 
		B->m, B->n, A->m, A->n);
      return NULL;
    }
  else if (B->symmetry_type != UNSYMMETRIC || A->symmetry_type != UNSYMMETRIC)
    {
      bebop_log (0, "*** csr_matrix_matmatmult: currently sparse matrix-matrix "
		"multiplication is only supported for matrices using unsymmetr"
		"ic storage ***\n");
      return NULL;
    }
  else if (B->value_type != A->value_type)
    {
      bebop_log (0, "*** csr_matrix_matmatmult: currently sparse matrix-matrix "
		"multiplication is only supported for matrices of the same val"
		"ue type ***\n");
      return NULL;
    }

  if (value_type == REAL)
    {
      double* Cval = NULL;
      int* Cptr = NULL; 
      int* Cind = NULL;
      int nnz = 0;
      double* Bvalues = (double*) B->values;
      double* Avalues = (double*) A->values;

      errcode = csr_matmatmult_double_real (&Cptr, &Cind, &Cval, &nnz, 1.0, 
					    B->rowptr, B->colidx, Bvalues,
					    A->rowptr, A->colidx, Avalues,
					    m, p, n);
      if (errcode != 0)
	{
	  bebop_log (0, "*** csr_matrix_matmatmult: matrix-matrix multiplication "
		    "kernel failed with error code %d ***\n", errcode);
	  return NULL;
	}

      if (nnz > 0) 
	{
	  assert (Cval != NULL);
	  assert (Cind != NULL);
	}
      if (m > 0)
	assert (Cptr != NULL);

      {
        struct csr_matrix_t* C = create_csr_matrix (m, n, nnz, Cval, Cind, Cptr, 
						    UNSYMMETRIC, 0, REAL, 
						    LIBRARY_DEALLOCATES, NULL, 
						    NO_COPY);
	assert (C != NULL);
	bebop_log (2, "=== Done with csr_matrix_matmatmult ===\n");
	return C;
      }
    }
  else if (value_type == COMPLEX)
    {
      struct csr_matrix_t* C = NULL;
      double_Complex* Cval = NULL;
      int* Cptr = NULL; 
      int* Cind = NULL;
      int nnz = 0;
      double_Complex* Bvalues = (double_Complex*) B->values;
      double_Complex* Avalues = (double_Complex*) A->values;

      errcode = csr_matmatmult_double_complex (&Cptr, &Cind, &Cval, &nnz, 
					       new_double_Complex (1.0, 0.0), 
					       B->rowptr, B->colidx, Bvalues,
					       A->rowptr, A->colidx, Avalues,
					       m, p, n);
      if (errcode != 0)
	{
	  bebop_log (0, "*** csr_matrix_matmatmult: error code %d ***\n", 
		    errcode);
	  return NULL;
	}
      C = create_csr_matrix (m, n, nnz, Cval, Cind, Cptr, UNSYMMETRIC, 0, 
			     COMPLEX, LIBRARY_DEALLOCATES, NULL, NO_COPY);
      bebop_log (2, "=== Done with csr_matrix_matmatmult ===\n");
      return C;
    }
  else /* value_type == PATTERN */
    {
      struct csr_matrix_t* C = bebop_malloc (sizeof (struct csr_matrix_t));

      errcode = csr_matmatmult_pattern (&(C->rowptr), &(C->colidx), &(C->nnz), 
					B->rowptr, B->colidx, 
					A->rowptr, A->colidx, 
					m, p, n);
      if (errcode != 0)
	{
	  bebop_log (0, "*** csr_matrix_matmatmult: error code %d ***\n", errcode);
	  return NULL;
	}
      C->m = m;
      C->n = n;
      C->symmetry_type = UNSYMMETRIC;
      C->symmetric_storage_location = 0;
      C->value_type = PATTERN;
      C->ownership = LIBRARY_DEALLOCATES;
      C->deallocator = &free;
      return C;
    }
}

static void
csr_matrix_scale_by_rowsums (struct csr_matrix_t* A)
{
  const int m = A->m;
  const int* ptr = A->rowptr;
  const enum value_type_t value_type = A->value_type;
  int i, k;

  bebop_log (2, "=== csr_matrix_scale_by_rowsums ===\n");

#ifdef USING_VALGRIND_CLIENT_REQUESTS
  {
    int retval = 0;
    retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
    if (retval != 0)
      {
	bebop_log (0, "*** csr_matrix_scale_by_rowsums: Pointer A "
		  "is invalid! ***\n");
	bebop_exit (EXIT_FAILURE);
      }
  }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

  /* FIXME: I have to finish the symmetric, etc. cases yet */
  assert (A->symmetry_type == UNSYMMETRIC);

  if (value_type == REAL)
    { 
      double* val = (double*) (A->values);
      for (i = 0; i < m; i++) 
	{
	  const int start = ptr[i];
	  const int end = ptr[i+1];
	  double s = 0.0;
	  
	  /* Figure out the row sum */
	  for (k = start; k < end; k++)
	    {
	      s += val[k];
	    }
	  
	  /* Scale the row by the row sum */
	  for (k = start; k < end; k++)
	    val[k] /= s;
	}
    }
  else if (value_type == COMPLEX)
    {
      double_Complex* val = (double_Complex*) (A->values);
      for (i = 0; i < m; i++) 
	{
	  const int start = ptr[i];
	  const int end = ptr[i+1];
	  double_Complex s = double_Complex_ZERO;
	  
	  /* Figure out the row sum */
	  for (k = start; k < end; k++)
	    {
	      s = double_Complex_add (s, val[k]);
	    }
	  
	  /* Scale the row by the row sum */
	  for (k = start; k < end; k++)
	    val[k] = double_Complex_divide (val[k], s);
	}
    }

  bebop_log (2, "=== Done with csr_matrix_scale_by_rowsums ===\n");
}


struct csr_matrix_t*
csr_matrix_triple_product (struct csr_matrix_t* RT,
			   struct csr_matrix_t* A,
			   struct csr_matrix_t* P)
{
  const enum value_type_t value_type = A->value_type;
  const int m = P->n; 
  const int n = A->n; 
  void* val = NULL; 
  int* ind = NULL; 
  int* ptr = NULL;
  int errcode = 0; 
  
  bebop_log (2, "=== csr_matrix_triple_product ===\n");

#ifdef USING_VALGRIND_CLIENT_REQUESTS
  if (bebop_debug_level () > 1)
    {
      int retval = 0;
      bebop_log (2, "Checking if RT pointer is valid\n");
      retval = VALGRIND_CHECK_READABLE( RT, sizeof (struct csr_matrix_t) );
      if (retval != 0)
	{
	  bebop_log (2, "*** csr_matrix_triple_product: RT pointer "
		    "is invalid! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
      bebop_log (2, "Checking if A pointer is valid\n");
      retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
      if (retval != 0)
	{
	  bebop_log (2, "*** csr_matrix_triple_product: A pointer "
		    "is invalid! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
      bebop_log (2, "Checking if P pointer is valid\n");
      retval = VALGRIND_CHECK_READABLE( P, sizeof (struct csr_matrix_t) );
      if (retval != 0)
	{
	  bebop_log (2, "*** csr_matrix_triple_product: P pointer "
		    "is invalid! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
      bebop_log (2, "Pointer validity checks succeeded.\n");
    }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */


   if (A->n != A->m) 
     return NULL; 
   else if (A->n != P->m) 
     return NULL; 
   else if (RT->m != P->m) 
     return NULL; 
   else if (RT->n != P->n) 
     return NULL; 
   else if (value_type != RT->value_type || value_type != P->value_type) 
     return NULL; 
   else if (A->symmetry_type != UNSYMMETRIC)
     return NULL;

   if (RT == P)
     {
       double* __val = NULL;
       int nnz = 0;

       if (value_type != REAL)
	 return NULL;

       errcode = ptap_csr_dr (A->m, P->n, 
			      A->values, A->colidx, A->rowptr,
			      P->values, P->colidx, P->rowptr, 
			      &__val, &ind, &ptr, &nnz);
       val = (void*) __val;
     }
   else
     {
       if (value_type == REAL) 
	 { 
	   double* _val = (double*) val; 
	   errcode = csr_matrix_triple_product_kernel_double_real (m, n,  
								   RT->values, RT->colidx,  
								   RT->rowptr, A->values,  
								   A->colidx, A->rowptr, 
								   P->values, P->colidx,  
								   P->rowptr, 
								   &_val, &ind, &ptr); 
	   val = (void*) _val; 
	 } 
       else if (value_type == COMPLEX)
	 {
	   double_Complex* _val = (double_Complex*) val;
	   errcode = csr_matrix_triple_product_kernel_double_complex (m, n, RT->values, RT->colidx,
								      RT->rowptr, A->values, A->colidx,
								      A->rowptr, P->values, P->colidx,
								      P->rowptr, &_val, &ind, &ptr);
	   val = (void*) _val;
	 }
       else if (value_type == PATTERN) 
	 errcode = csr_matrix_triple_product_kernel_double_pattern (m, n,  
								    RT->colidx, RT->rowptr,  
								    A->colidx, A->rowptr,
								    P->colidx, P->rowptr,  
								    &ind, &ptr); 
       else  
	 return NULL; 
     }

   if (errcode != 0) 
     return NULL; 
   else 
     return create_csr_matrix (m, m, ptr[m], val, ind, ptr, 
 			      UNSYMMETRIC, 0, value_type,
			      LIBRARY_DEALLOCATES, &free, 
 			      NO_COPY); 
}


int
csr_matrix_expand_to_dense (void* B, 
			    const int col_oriented, 
			    const int ldb, 
			    const struct csr_matrix_t* A)
{
  const enum value_type_t value_type = A->value_type;
  const enum symmetry_type_t symmetry_type = A->symmetry_type;

  bebop_log (2, "=== csr_matrix_expand_to_dense ===\n");

  if (bebop_debug_level () > 1)
    {
#ifdef USING_VALGRIND_CLIENT_REQUESTS
      {
	int retval = 0;
	retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
	if (retval != 0)
	  {
	    bebop_log (0, "*** csr_matrix_expand_to_dense: Pointer A "
		      "is invalid! ***\n");
	    bebop_exit (EXIT_FAILURE);
	  }
      }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

      if (! valid_csr_matrix_p (A))
	{
	  bebop_log (0, "*** csr_matrix_expand_to_dense: input CSR "
		    "matrix is invalid! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
    }

  if (symmetry_type != UNSYMMETRIC)
    {
      /* This is inefficient but it works */
      struct csr_matrix_t* expandedA = clone_csr_matrix (A);
      if (expandedA == NULL)
	{
	  bebop_log (0, "*** csr_matrix_expand_to_dense: failed "
		    "to expand A into unsymmetric form! ***\n");
	  return -1;
	}
      
      if (value_type == REAL)
	csr_expand_to_dense_kernel_unsymmetric_real_double ((double*) B, 
							    col_oriented, 
							    expandedA->m, 
							    expandedA->n, 
							    ldb, 
							    expandedA->values, 
							    expandedA->colidx, 
							    expandedA->rowptr);
      else if (value_type == COMPLEX)
	csr_expand_to_dense_kernel_unsymmetric_complex_double ((double_Complex*) B, 
							       col_oriented, 
							       expandedA->m, 
							       expandedA->n, 
							       ldb, 
							       expandedA->values, 
							       expandedA->colidx, 
							       expandedA->rowptr);
      destroy_csr_matrix (expandedA);
    }
  else /* unsymmetric-stored matrix */
    {
      if (value_type == REAL)
	csr_expand_to_dense_kernel_unsymmetric_real_double ((double*) B, col_oriented, 
							    A->m, A->n, ldb, A->values, 
							    A->colidx, A->rowptr);
      else if (value_type == COMPLEX)
	csr_expand_to_dense_kernel_unsymmetric_complex_double ((double_Complex*) B, 
							       col_oriented, 
							       A->m, A->n, ldb, 
							       A->values, A->colidx, 
							       A->rowptr);
    }

  bebop_log (2, "=== Done with csr_matrix_expand_to_dense ===\n");
  return 0;
}

void
csr_matrix_weighted_jacobi (void* y, 
			    const void* weight, 
			    const struct csr_matrix_t* A, 
			    const void* x,
			    const void* r)
{
  const enum value_type_t value_type = A->value_type;

  bebop_log (2, "=== csr_matrix_weighted_jacobi ===\n");
  if (bebop_debug_level () > 1)
    {
#ifdef USING_VALGRIND_CLIENT_REQUESTS
      int retval = 0;
      retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
      if (retval != 0)
	{
	  bebop_log (0, "*** csr_matrix_weighted_jacobi: Pointer A "
		    "is invalid! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

      if (! valid_csr_matrix_p (A))
	{
	  bebop_log (0, "*** csr_matrix_weighted_jacobi: invalid input "
		    "matrix! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}

#ifdef USING_VALGRIND_CLIENT_REQUESTS
      if (A->value_type == REAL)
	{ 
	  retval = VALGRIND_CHECK_READABLE( weight, sizeof (double) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** csr_matrix_weighted_jacobi: Valgrind says "
			"weight is an invalid pointer! ***\n");
	      bebop_exit (EXIT_FAILURE);
	    } 
	  retval = VALGRIND_CHECK_READABLE( r, (A->m) * sizeof (double) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** csr_matrix_weighted_jacobi: Valgrind says "
			"r is an invalid pointer! ***\n");
	      bebop_exit (EXIT_FAILURE);
	    } 
	  retval = VALGRIND_CHECK_READABLE( y, (A->n) * sizeof (double) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** csr_matrix_weighted_jacobi: Valgrind says "
			"y is an invalid pointer! ***\n");
	      bebop_exit (EXIT_FAILURE);
	    } 
	}
      else if (A->value_type == COMPLEX)
	{
	  retval = VALGRIND_CHECK_READABLE( weight, sizeof (double_Complex) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** csr_matrix_weighted_jacobi: Valgrind says "
			"weight is an invalid pointer! ***\n");
	      bebop_exit (EXIT_FAILURE);
	    } 
	  retval = VALGRIND_CHECK_READABLE( r, (A->m) * sizeof (double_Complex) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** csr_matrix_weighted_jacobi: Valgrind says "
			"r is an invalid pointer! ***\n");
	      bebop_exit (EXIT_FAILURE);
	    } 
	  retval = VALGRIND_CHECK_READABLE( y, (A->n) * sizeof (double_Complex) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** csr_matrix_weighted_jacobi: Valgrind says "
			"y is an invalid pointer! ***\n");
	      bebop_exit (EXIT_FAILURE);
	    } 
	}
#endif /* USING_VALGRIND_CLIENT_REQUESTS */
    }

  if (value_type == REAL)
    {
      double w = *((double*) weight);
      csr_weighted_jacobi_kernel_real_double ((double*) y, w, A->values, 
					      A->colidx, A->rowptr, A->n, 
					      (const double*) x,
					      (const double*) r);
    }
  else if (value_type == COMPLEX)
    {
      double_Complex w = *((double_Complex*) weight);
      csr_weighted_jacobi_kernel_complex_double ((double_Complex*) y, w, A->values, 
					         A->colidx, A->rowptr, A->n, 
						 (const double_Complex*) x,
					         (const double_Complex*) r);
    }
  bebop_log (2, "=== Done with csr_matrix_weighted_jacobi ===\n");
}


int 
csr_matrix_num_rows (const struct csr_matrix_t* A)
{
#ifdef USING_VALGRIND_CLIENT_REQUESTS
  {
    int retval = 0;
    retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
    if (retval != 0)
      {
	bebop_log (0, "*** csr_matrix_num_rows: input pointer A is "
		  "invalid! ***\n");
	bebop_exit (EXIT_FAILURE);
      }
  }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */
  return A->m;
}

int 
csr_matrix_num_cols (const struct csr_matrix_t* A)
{
#ifdef USING_VALGRIND_CLIENT_REQUESTS
  {
    int retval = 0;
    retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
    if (retval != 0)
      {
	bebop_log (0, "*** csr_matrix_num_cols: input pointer A is "
		  "invalid! ***\n");
	bebop_exit (EXIT_FAILURE);
      }
  }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */
  return A->n;
}


void
csr_spmv (void* y, const void* beta, const void* alpha, 
	  const struct csr_matrix_t* A, const void* x)
{
  enum value_type_t value_type;

  bebop_log (2, "=== csr_spmv ===\n");

#ifdef USING_VALGRIND_CLIENT_REQUESTS
  {
    int retval = 0;
    retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
    if (retval != 0)
      {
	bebop_log (0, "*** csr_spmv: input pointer A is "
		  "invalid! ***\n");
	bebop_exit (EXIT_FAILURE);
      }
  }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

  value_type = A->value_type;
  if (value_type == REAL)
    {
      const double _alpha = *( (const double*) alpha );
      const double _beta = *( (const double*) beta );

      csr_spmv_double_real ((double*) y, A->m, _beta, _alpha, 
			    (const double*) A->values, A->colidx, A->rowptr, 
			    (double*) x, A->n);
    }
  else if (value_type == COMPLEX)
    {
      const double_Complex _alpha = *( (const double_Complex*) alpha );
      const double_Complex _beta = *( (const double_Complex*) beta );

      csr_spmv_double_complex ((double_Complex*) y, A->m, _beta, _alpha,
			       (const double_Complex*) A->values,
			       A->colidx, A->rowptr, (double_Complex*) x,
			       A->n);
    }
  bebop_log (2, "=== Done with csr_spmv ===\n");
}

void 
csr_upper_trisolve (void* x, const struct csr_matrix_t* A)
{
  const enum value_type_t value_type = A->value_type;

  if (A->m != A->n)
    {
      bebop_log (0, "*** csr_upper_trisolve: A is not a square matrix! ***\n");
      bebop_exit (EXIT_FAILURE);
    }
  if (value_type == REAL)
    csr_upper_trisolve_double_real ((double*) x, A->values, A->colidx, A->rowptr, A->m);
  else if (value_type == COMPLEX)
    csr_upper_trisolve_double_complex ((double_Complex*) x, A->values, A->colidx, A->rowptr, A->m);
}

void
csr_lower_trisolve (void* x, const struct csr_matrix_t* A) 
{
  const enum value_type_t value_type = A->value_type;
  if (A->m != A->n)
    {
      bebop_log (0, "*** csr_lower_trisolve: A is not a square matrix! ***\n");
      bebop_exit (EXIT_FAILURE);
    }
  if (value_type == REAL)
    csr_lower_trisolve_double_real ((double*) x, A->values, A->colidx, A->rowptr, A->m);
  else if (value_type == COMPLEX)
    csr_lower_trisolve_double_complex ((double_Complex*) x, A->values, A->colidx, A->rowptr, A->m);
}

void
csr_spmv_transpose (void* y, const void* beta, const void* alpha, 
		    const struct csr_matrix_t* A, const void* x)
{
  enum value_type_t value_type;

  bebop_log (2, "=== csr_spmv_transpose ===\n");

  if (bebop_debug_level () > 1)
    {
#ifdef USING_VALGRIND_CLIENT_REQUESTS
      int retval = 0;
      retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
      if (retval != 0)
	{
	  bebop_log (0, "*** csr_spmv_transpose: pointer to matrix "
		    "is invalid! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
#endif /* USING_VALGRIND_CLIENT_REQUESTS */
      if (! valid_csr_matrix_p (A))
	{
	  bebop_log (0, "*** csr_spmv_transpose: input matrix is "
		    "invalid! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
    }

  value_type = A->value_type;
  if (value_type == REAL)
    {
      const double _alpha = *( (const double*) alpha );
      const double _beta = *( (const double*) beta );

      csr_spmv_transpose_double_real ((double*) y, 
				      A->m, 
				      _beta, 
				      _alpha, 
				      (const double*) A->values, 
				      A->colidx, 
				      A->rowptr, 
				      (double*) x, 
				      A->n);
    }
  else if (value_type == COMPLEX)
    {
      const double_Complex _alpha = *( (const double_Complex*) alpha );
      const double_Complex _beta = *( (const double_Complex*) beta );

      csr_spmv_transpose_double_complex ((double_Complex*) y, 
					 A->m, 
					 _beta, 
					 _alpha, 
					 (const double_Complex*) A->values,
					 A->colidx, 
					 A->rowptr, 
					 (double_Complex*) x,
					 A->n);
    }
  bebop_log (2, "=== Done with csr_spmv_transpose ===\n");
}


int
display_csr_matrix_in_matrix_market_format (const struct csr_matrix_t* A)
{
  return print_csr_matrix_in_matrix_market_format (stdout, A);
}



void
csr_matrix_restrict (void* rc, const struct csr_matrix_t* RT, const void* r)
{
  int i, j, k;

  bebop_log (2, "=== csr_matrix_restrict ===\n");
  if (bebop_debug_level () > 1)
    {
#ifdef USING_VALGRIND_CLIENT_REQUESTS
      int retval = 0;
      retval = VALGRIND_CHECK_READABLE( RT, sizeof (struct csr_matrix_t) );
      if (retval != 0)
	{
	  bebop_log (0, "*** csr_matrix_restrict: pointer to matrix RT "
		    "is invalid! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

      if (! valid_csr_matrix_p (RT))
	{
	  bebop_log (0, "*** csr_matrix_restrict: input matrix is invalid! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
    }

  enum value_type_t value_type = RT->value_type;
  if (value_type == REAL)
    {
      const double* val = RT->values;
      const int* ind = RT->colidx;
      const int* ptr = RT->rowptr;
      double* _rc = (double*) rc;
      double* _r  = (double*) r;

      /* Fill rc with zeros */
      for (i = 0; i < RT->n; i++)
	_rc[i] = 0.0;
      
      /* For each column of R, do a daxpy / scatter */
      for (j = 0; j < RT->m; j++)
	{
	  const int start = ptr[j];
	  const int end = ptr[j+1];
	  const double r_j = _r[j];
	  for (k = start; k < end; k++)
	    {
	      const int row = ind[k];
	      _rc[row] = _rc[row] + val[k] * r_j;
	    }
	}
    }
  else if (value_type == COMPLEX)
    {
      const double_Complex* val = RT->values;
      const int* ind = RT->colidx;
      const int* ptr = RT->rowptr;
      double_Complex* _rc = (double_Complex*) rc;
      double_Complex* _r  = (double_Complex*) r;

      /* Fill rc with zeros */
      for (i = 0; i < RT->n; i++)
	_rc[i] = double_Complex_ZERO; 
      
      /* For each column of R, do a daxpy / scatter */
      for (j = 0; j < RT->m; j++)
	{
	  const int start = ptr[j];
	  const int end = ptr[j+1];
	  const double_Complex r_j = _r[j];
	  for (k = start; k < end; k++)
	    {
	      const int row = ind[k];
	      _rc[row] = double_Complex_add (_rc[row], double_Complex_multiply (val[k], r_j));
	    }
	}
    }

  bebop_log (2, "=== Done with csr_matrix_restrict ===\n");
}



struct csr_matrix_t* 
csr_matrix_transpose (const struct csr_matrix_t* A)
{
  /* NOTE: Code borrowed from Hypre (seq_mv/csr_matop.c, version
     1.11.1b). */
  int i, j;
  enum value_type_t value_type;

  bebop_log (2, "=== csr_matrix_transpose ===\n");
  if (bebop_debug_level () > 1)
    {
#ifdef USING_VALGRIND_CLIENT_REQUESTS
      {
	int retval = 0;
	bebop_log (2, "Checking if pointer A is valid\n");
	retval = VALGRIND_CHECK_READABLE( A, sizeof (struct csr_matrix_t) );
	if (retval != 0)
	  {
	    bebop_log (0, "*** csr_matrix_transpose: Pointer A is invalid! ***\n");
	    bebop_exit (EXIT_FAILURE);
	  }
	bebop_log (2, "Pointer A is valid.\n");
      }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */
      if (! valid_csr_matrix_p (A))
	{
	  bebop_log (0, "*** csr_matrix_transpose: input matrix is invalid! ***\n");
	  return NULL;
	}
    }

  value_type = A->value_type;
  if (A->symmetry_type == SYMMETRIC)
    {
      bebop_log (2, "Matrix is symmetric, returning copy of A\n");
      return clone_csr_matrix (A);
    }
  else if (A->symmetry_type != UNSYMMETRIC) 
    {
      fprintf (stderr, "*** csr_matrix_transpose: operation not yet "
	       "implemented for matrices of symmetry type %d ***\n",
	       A->symmetry_type);
      return NULL;
    }

  if (value_type == REAL)
    {
      int* ptr = (int*) bebop_calloc (A->n + 1, sizeof (int));
      int* ind = (int*) bebop_calloc (A->nnz, sizeof (int));
      double* val = (double*) bebop_calloc (A->nnz, sizeof (double));

      /* 
       * Count the number of entries in each column if A (row of A^T)
       * and fill the ptr array (of A^T).
       */

      /* 
       * Count the number of entries in each column of A (row of A^T).  
       * Store the counts in the new ptr array of A^T.  Later we will
       * fix up ptr so that it contains row pointers instead of counts.
       */
      bebop_log (2, "\tCounting # of entries in each col of A\n");
      for (i = 0; i < A->nnz; i++)
	{
	  int j = A->colidx[i];
	  if (j < 0 || j >= A->n)
	    {
	      bebop_log (0, "*** csr_matrix_transpose: A->colidx[%d] = %d "
			"is out of range [%d,%d] ***\n", i, j, 
			0, A->n - 1);
	      bebop_free (ptr);
	      bebop_free (ind);
	      bebop_free (val);
	      return NULL;
	    }
	  assert (j+1 >= 0 && j+1 < A->n + 1);
	  ptr[j+1] = ptr[j+1] + 1;
	}
      /*
       * Fix up ptr so that ptr[i+1] - ptr[i] is the number of nonzeros
       * in row i of A^T.
       */
      bebop_log (2, "\tFixing up ptr\n");
#ifdef USING_VALGRIND_CLIENT_REQUESTS
      {
	int e = 0;
	bebop_log (2, "Checking if ptr is valid\n");
	e = VALGRIND_CHECK_READABLE(ptr, (A->n + 1) * sizeof (int));
	if (e != 0)
	  {
	    bebop_log (0, "*** csr_matrix_transpose: ptr is an invalid pointer -- it does not point to n+1 = %d integers! ***\n", A->n + 1);
	    bebop_exit (EXIT_FAILURE);
	  }
	bebop_log (2, "ptr is valid.\n");
      }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */
      for (i = 2; i < A->n; i++)
	{
	  ptr[i] += ptr[i-1];
	}
      bebop_log (2, "\tDone fixing up ptr\n");

      /* 
       * Load the data and column numbers of A^T. 
       */
      bebop_log (2, "\tLoading data and column indices of A transpose\n");
      for (i = 0; i < A->m; i++)
	{
	  for (j = A->rowptr[i]; j < A->rowptr[i+1]; j++)
	    {
	      if (ptr[A->colidx[j]] < 0 || ptr[A->colidx[j]] >= A->nnz)
		{
		  bebop_log (0, "*** csr_matrix_transpose: ptr[%d] = %d "
			    "for j = %d is out of range [%d,%d] ***\n", 
			    A->colidx[j], ptr[A->colidx[j]], j, 0, A->nnz - 1);
	          bebop_free (ptr);
	          bebop_free (ind);
	          bebop_free (val);
		}
	      ind[ptr[A->colidx[j]]] = i;
	      val[ptr[A->colidx[j]]] = ((double*) (A->values))[j];
	      ptr[A->colidx[j]]++;
	    }
	}

      /*
       * ptr[j] now points to the *end* of the j-th row of entries
       * instead of the beginning.  Restore ptr to front of row.
       */
      bebop_log (2, "\tFixing up ptr for the second time\n");
      for (i = A->n; i > 0; i--)
	{
	  ptr[i] = ptr[i-1];
	}
      ptr[0] = 0;

      {
	struct csr_matrix_t* Retval = create_csr_matrix (A->n, A->m, A->nnz, (void*) val, ind, ptr, A->symmetry_type, 1 - A->symmetric_storage_location, A->value_type, LIBRARY_DEALLOCATES, NULL, NO_COPY);
	bebop_log (2, "=== Done with csr_matrix_transpose ===\n");
	return Retval;
      }
    }
  else if (value_type == COMPLEX)
    {
      int* ptr = bebop_calloc (A->n + 1, sizeof (int));
      int* ind = bebop_calloc (A->nnz, sizeof (int));
      double_Complex* val = bebop_calloc (A->nnz, sizeof (double_Complex));

      /* 
       * Count the number of entries in each column if A (row of A^T)
       * and fill the ptr array (of A^T).
       */

      /* 
       * Count the number of entries in each column of A (row of A^T).  
       * Store the counts in the new ptr array of A^T.  Later we will
       * fix up ptr so that it contains row pointers instead of counts.
       */
      for (i = 0; i < A->nnz; i++)
	{
	  assert (A->colidx[i] >= 0 && A->colidx[i] < A->n);
	  ++ptr[A->colidx[i]+1];
	}
      /*
       * Fix up ptr so that ptr[i+1] - ptr[i] is the number of nonzeros
       * in row i of A^T.
       */
      for (i = 2; i < A->n; i++)
	{
	  ptr[i] += ptr[i-1];
	}

      /* 
       * Load the data and column numbers of A^T. 
       */
      for (i = 0; i < A->m; i++)
	{
	  for (j = A->rowptr[i]; j < A->rowptr[i+1]; j++)
	    {
	      assert (ptr[A->colidx[j]] >= 0);
	      assert (ptr[A->colidx[j]] < A->nnz);
	      ind[ptr[A->colidx[j]]] = i;
	      val[ptr[A->colidx[j]]] = ((double_Complex*) (A->values))[j];
	      ptr[A->colidx[j]]++;
	    }
	}

      /*
       * ptr[j] now points to the *end* of the j-th row of entries
       * instead of the beginning.  Restore ptr to front of row.
       */
      for (i = A->n; i > 0; i--)
	{
	  ptr[i] = ptr[i-1];
	}
      ptr[0] = 0;

      bebop_log (2, "=== Done with csr_matrix_transpose ===\n");
      return create_csr_matrix (A->n, A->m, A->nnz, (void*) val, ind, ptr, 
				A->symmetry_type, 1 - A->symmetric_storage_location, 
				A->value_type, LIBRARY_DEALLOCATES, NULL, NO_COPY);
    }
  else if (value_type == PATTERN)
    {
      int* ptr = bebop_calloc (A->n + 1, sizeof (int));
      int* ind = bebop_calloc (A->nnz, sizeof (int));

      /* 
       * Count the number of entries in each column if A (row of A^T)
       * and fill the ptr array (of A^T).
       */

      /* 
       * Count the number of entries in each column of A (row of A^T).  
       * Store the counts in the new ptr array of A^T.  Later we will
       * fix up ptr so that it contains row pointers instead of counts.
       */
      for (i = 0; i < A->nnz; i++)
	{
	  assert (A->colidx[i] >= 0 && A->colidx[i] < A->n);
	  ++ptr[A->colidx[i]+1];
	}
      /*
       * Fix up ptr so that ptr[i+1] - ptr[i] is the number of nonzeros
       * in row i of A^T.
       */
      for (i = 2; i < A->n; i++)
	{
	  ptr[i] += ptr[i-1];
	}

      /* 
       * Load the column numbers of A^T. 
       */
      for (i = 0; i < A->m; i++)
	{
	  for (j = A->rowptr[i]; j < A->rowptr[i+1]; j++)
	    {
	      assert (ptr[A->colidx[j]] >= 0);
	      assert (ptr[A->colidx[j]] < A->nnz);
	      ind[ptr[A->colidx[j]]] = i;
	      ptr[A->colidx[j]]++;
	    }
	}

      /*
       * ptr[j] now points to the *end* of the j-th row of entries
       * instead of the beginning.  Restore ptr to front of row.
       */
      for (i = A->n; i > 0; i--)
	{
	  ptr[i] = ptr[i-1];
	}
      ptr[0] = 0;

      bebop_log (2, "=== Done with csr_matrix_transpose ===\n");
      return create_csr_matrix (A->n, A->m, A->nnz, (void*) NULL, ind, ptr, 
				A->symmetry_type, 1 - A->symmetric_storage_location, 
				A->value_type, LIBRARY_DEALLOCATES, NULL, NO_COPY);
    }

 return NULL; /* sentinel */
}



int
csr_matrix_equal_p (const struct csr_matrix_t* A, const struct csr_matrix_t* B)
{
  struct csr_matrix_t Acopy;
  struct csr_matrix_t Bcopy;
  int retval = 0;
  int i, j;
  const enum value_type_t value_type = A->value_type;

  bebop_log (2, "=== csr_matrix_equal_p ===\n");
  if (bebop_debug_level () > 1)
    {
      int failure = 0;
      if (! valid_csr_matrix_p (A))
	{
	  bebop_log (0, "*** csr_matrix_equal_p: First input matrix is invalid! ***\n");
	  failure = 1;
	}
      if (! valid_csr_matrix_p (B))
	{
	  bebop_log (0, "*** csr_matrix_equal_p: Second input matrix is invalid! ***\n");
	  failure = 1;
	}
      if (failure)
	bebop_exit (EXIT_FAILURE);
    }

  copy_csr_matrix (&Acopy, A);
  copy_csr_matrix (&Bcopy, B);

  assert (0 == csr_matrix_expand_symmetric_storage (&Acopy));
  assert (0 == csr_matrix_expand_symmetric_storage (&Bcopy));

  /* Sort the rows by column indices so we can do comparisons easily */
  csr_matrix_sort_colidx (&Acopy);
  csr_matrix_sort_colidx (&Bcopy);

  /* If Not Equal Return Zero */
#define INERZ( thing ) do { \
    if (Acopy.thing != Bcopy.thing) { \
      bebop_log (2, "### csr_matrix_equal_p: %s not equal ###\n", #thing ); \
      goto cleanup; \
    } \
  } while(0)

  INERZ(m);
  INERZ(n);
  INERZ(nnz);
  INERZ(value_type);

  /* If Not Equal Array element Return Zero */
#define INEARZ( array, idx ) do { \
    if (Acopy.array[idx] != Bcopy.array[idx]) { \
      bebop_log (2, "### csr_matrix_equal_p: array %s element %d not equal ###\n", #array , idx ); \
      goto cleanup; \
    } \
  } while(0)

  if (value_type == REAL)
    {
      for (i = 0; i < Acopy.m; i++)
	{
	  int start, end;

	  INEARZ(rowptr, i);
	  INEARZ(rowptr, i+1);
	  start = Acopy.rowptr[i];
	  end = Acopy.rowptr[i+1];
	  for (j = start; j < end; j++)
	    {
	      INEARZ(colidx, j);

	      if (((double*) (Acopy.values))[j] != ((double*) (Bcopy.values))[j])
		{
		  bebop_log (2, "### csr_matrix_equal_p: values[%d] not equal ###\n", j);
		  goto cleanup;
		}
	    }
	}
    }
  else if (value_type == COMPLEX)
    {
      for (i = 0; i < Acopy.m; i++)
	{
	  int start, end;

	  INEARZ(rowptr, i);
	  INEARZ(rowptr, i+1);
	  start = Acopy.rowptr[i];
	  end = Acopy.rowptr[i+1];
	  for (j = start; j < end; j++)
	    {
	      INEARZ(colidx, j);
	      if (double_Complex_not_equal (((double_Complex*) (Acopy.values))[j], 
					    ((double_Complex*) (Bcopy.values))[j]))
		{
		  bebop_log (2, "### csr_matrix_equal_p: values[%d] not equal ###\n", j);
		  goto cleanup;
		}
	    }
	}
    }
  else if (value_type == PATTERN)
    {
      for (i = 0; i < Acopy.m; i++)
	{
	  int start, end;

	  INEARZ(rowptr, i);
	  INEARZ(rowptr, i+1);
	  start = Acopy.rowptr[i];
	  end = Acopy.rowptr[i+1];
	  for (j = start; j < end; j++)
	    {
	      INEARZ(colidx, j);
	    }
	}
    }
  
  retval = 1; /* passed all the equality tests */
  bebop_log (2, "### csr_matrix_equal_p: matrices are equal ###\n");
 cleanup:
  dealloc_csr_matrix (&Acopy);
  dealloc_csr_matrix (&Bcopy);
  bebop_log (2, "=== Done with csr_matrix_equal_p ===\n");
  return retval;
}


static void
csr_gauss_seidel_double_real_unsymmetric (double* x, 
					  const int m,
					  const int* ptr,
					  const int* ind,
					  const double* val,
					  const double* b, 
					  const char* direction)
{
  int i, k;

  if (direction[0] == 'F' || direction[0] == 'f')
    {
      for (i = 0; i < m; i++)
	{
	  const int start = ptr[i];
	  const int end = ptr[i+1];
	  double s = 0.0;
	  double d = 0.0;
	  
	  for (k = start; k < end; k++)
	    {
	      const int j = ind[k];
	      if (i != j)
		s += val[k] * x[j]; 
	      else
		d += val[k];
	    }
	  if (d == 0.0)
	    {
	      bebop_log (0, "*** csr_gauss_seidel_double_real"
			"_unsymmetric: at i = %d, diagonal "
			"element is zero for a row with %d "
			"elements, when doing forward Gauss-"
			"Seidel! ***\n", i, end - start);
	      bebop_exit (EXIT_FAILURE);
	    }
	  x[i] = (b[i] - s) / d;
	}
    }
  else if (direction[0] == 'B' || direction[0] == 'b')
    {
      for (i = m - 1; i >= 0; i--)
	{
	  const int start = ptr[i];
	  const int end = ptr[i+1];
	  double s = 0.0;
	  double d = 0.0;
	  
	  for (k = start; k < end; k++)
	    {
	      const int j = ind[k];
	      if (i != j)
		s += val[k] * x[j]; 
	      else
		d += val[k];
	    }
	  if (d == 0.0)
	    {
	      bebop_log (0, "*** csr_gauss_seidel_double_real"
			"_unsymmetric: at i = %d, diagonal "
			"element is zero for a row with %d "
			"elements, when doing backward Gauss-"
			"Seidel! ***\n", i, end - start);
	      bebop_exit (EXIT_FAILURE);
	    }
	  x[i] = (b[i] - s) / d;
	}
    }
  else if (direction[0] == 's' || direction[0] == 'S')
    {
      csr_gauss_seidel_double_real_unsymmetric (x, m, ptr, ind, val, b, "forward");
      csr_gauss_seidel_double_real_unsymmetric (x, m, ptr, ind, val, b, "backward");
    }
}

static void
csr_gauss_seidel_double_complex_unsymmetric (double_Complex* x, const int m, 
					     const int* ptr, const int* ind, 
					     const double_Complex* val,
					     const double_Complex* b, 
					     const char* direction)
{
  int i, k;

  if (direction[0] == 'F' || direction[0] == 'f')
    {
      for (i = 0; i < m; i++)
	{
	  const int start = ptr[i];
	  const int end = ptr[i+1];
	  double_Complex s = double_Complex_ZERO;
	  double_Complex d = double_Complex_ZERO;
	  
	  for (k = start; k < end; k++)
	    {
	      const int j = ind[k];
	      if (i != j)
		s = double_Complex_add (s, double_Complex_multiply (val[k], x[j]));
	      else
		d = double_Complex_add (d, val[k]); /* could be more than one diagonal entry; sum them */
	    }
	  x[i] = double_Complex_divide (double_Complex_subtract (b[i], s), d);
	}
    }
  else if (direction[0] == 'B' || direction[0] == 'b')
    {
      for (i = m - 1; i >= 0; i--)
	{
	  const int start = ptr[i];
	  const int end = ptr[i+1];
	  double_Complex s = double_Complex_ZERO;
	  double_Complex d = double_Complex_ZERO;
	  
	  for (k = start; k < end; k++)
	    {
	      const int j = ind[k];
	      if (i != j)
		s = double_Complex_add (s, double_Complex_multiply (val[k], x[j]));
	      else
		d = double_Complex_add (d, val[k]); /* could be more than one diagonal entry; sum them */
	    }
	  x[i] = double_Complex_divide (double_Complex_subtract (b[i], s), d);
	}
    }
  else if (direction[0] == 's' || direction[0] == 'S')
    {
      csr_gauss_seidel_double_complex_unsymmetric (x, m, ptr, ind, val, b, "forward");
      csr_gauss_seidel_double_complex_unsymmetric (x, m, ptr, ind, val, b, "backward");
    }
}

void
csr_gauss_seidel (void* x, const struct csr_matrix_t* A, const void* b, const char* direction)
{
  enum value_type_t value_type = A->value_type;
  bebop_log (2, "=== csr_gauss_seidel ===\n");
  assert (direction != NULL);

  /* FIXME: we haven't finished the nonsymmetric case yet */
  assert (A->symmetry_type == UNSYMMETRIC);

  if (A->symmetry_type == UNSYMMETRIC)
    {
      if (value_type == REAL)
	csr_gauss_seidel_double_real_unsymmetric ((double*) x, A->m, 
						  A->rowptr, A->colidx, 
						  (const double*) (A->values), 
						  (const double*) b, direction);
      else if (value_type == COMPLEX)
	csr_gauss_seidel_double_complex_unsymmetric ((double_Complex*) x, 
						     A->m, A->rowptr, A->colidx,
						     (const double_Complex*) (A->values),
						     (const double_Complex*) b, direction);
    }
 
  bebop_log (2, "=== Done with csr_gauss_seidel ===\n");
}
