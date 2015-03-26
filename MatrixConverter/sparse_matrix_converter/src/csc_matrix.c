/**
 * @file csc_matrix.c
 * @author Mark Hoemmen
 * @since 07/01/04 11:25:31 PDT
 * @date Time-stamp: <2008-07-16 11:15:01 mhoemmen>
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
#include <bebop/smc/csc_matrix.h>
#include <bebop/smc/coo_matrix.h>
#include <bebop/smc/iohb.h>
#include <bebop/smc/read_mm.h>   /* csc_to_coo */
#include <bebop/smc/sparse_vector.h>

#include <bebop/util/complex.h>
#include <bebop/util/malloc.h>
#include <bebop/util/util.h>
#include <bebop/util/log.h>

#include <assert.h>
#include <math.h>   /* fabs */
#include <string.h> /* memcpy */

/**
 * Works like (void*) (&array[idx]), in which array is treated as an array
 * of objects, each of which satisfies sizeof(object) == size.  "size"
 * must be a variable or constant in scope ("deliberate variable
 * capture").
 */
#ifdef VOIDAREF
#  undef VOIDAREF
#endif /* VOIDAREF */
#define VOIDAREF( array, idx )  ((void*) ((char*) (array) + size*(idx)))  

/**
 * An (index,value) pair.  In our context of sparse matrices in CSC format, 
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
 * Takes column number j of the given sparse matrix A, and copies all the 
 * indices (and values, if there are any values) of that column into "col".
 *
 * @param A [IN]    Sparse matrix in CSC format
 * @param j [IN]    Index of the column to copy (zero-based indices)
 * @param col [OUT] An array of either index_real_value_pair_t, 
 *                  index_complex_value_pair_t or int, depending on whether 
 *                  the value type of the matrix is REAL, COMPLEX, or PATTERN,
 *                  respectively.
 * @param max_nnz [IN]  The max number of nonzeros in any column of A
 */
static void
copy_col2pairs (const struct csc_matrix_t* A, int j, void* col, int max_nnz)
{
  int a = A->colptr[j];
  int b = A->colptr[j+1];
  int nnz = b - a;
  int k;

  assert (nnz <= max_nnz);

  if (A->value_type == REAL)
    {
      struct index_real_value_pair_t* _col = (struct index_real_value_pair_t*) col;
      const double* const values = (const double* const) (A->values);

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  _col[k].index = A->rowidx[a+k];
	  _col[k].value = values[a+k];
	}
    }
  else if (A->value_type == COMPLEX)
    {
      struct index_complex_value_pair_t* _col = (struct index_complex_value_pair_t*) col;
      const double_Complex* const values = (const double_Complex* const) (A->values);

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  _col[k].index = A->rowidx[a+k];
	  _col[k].value = values[a+k];
	}
    }
  else if (A->value_type == PATTERN)
    {
      int* _col = (int*) col;

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  _col[k] = A->rowidx[a+k];
	}
    }
}


/** 
 * Given a sparse matrix A in CSC format, a column number j, and a list of
 * (index,value) pairs, copies the list of (index,value) pairs back into that
 * column of the matrix.
 */
static void
copy_pairs2col (const void* col, int max_nnz, struct csc_matrix_t* A, int j)
{
  int a = A->colptr[j];
  int b = A->colptr[j+1];
  int nnz = b - a;

  int k;

  assert (nnz <= max_nnz);

  if (A->value_type == REAL)
    {
      struct index_real_value_pair_t* _col = (struct index_real_value_pair_t*) col;
      double* values = (double*) (A->values);

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  A->rowidx[a+k] = _col[k].index;
	  values[a+k] = _col[k].value;
	}
    }
  else if (A->value_type == COMPLEX)
    {
      struct index_complex_value_pair_t* _col = (struct index_complex_value_pair_t*) col;
      double_Complex* values = (double_Complex*) (A->values);

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  A->rowidx[a+k] = _col[k].index;
	  values[a+k] = _col[k].value;
	}
    }
  else if (A->value_type == PATTERN)
    {
      int* _col = (int*) col;

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  A->rowidx[a+k] = _col[k];
	}
    }
}


/**
 * Sorts the row indices within each column of the sparse matrix A.
 */
static void
csc_matrix_sort_rowidx (struct csc_matrix_t* A)
{
  int j;
  int max_nnz = 0;  /* Will be the max # of nonzeros in each column */

  bebop_log (2, "=== csc_matrix_sort_rowidx ===\n");
  bebop_log (2, "Matrix is %d x %d with %d nonzeros\n", A->m, A->n, A->nnz);

  if (A->n <= 0)
    {
      bebop_log (2, "=== Done with csc_matrix_sort_rowidx ===\n");
      return;
    }

  /* Find the max # of nonzeros in each column.  This lets us allocate 
   * workspace large enough to accomodate any column (otherwise we would
   * have to allocate it separately for _each_ column, which would waste
   * a lot of malloc calls). */
  bebop_log (2, "Finding max # of nonzeros in each column...");
  max_nnz = A->colptr[1] - A->colptr[0];
  assert (A->colptr[1] - A->colptr[0] >= 0);
  for (j = 1; j < A->n; j++)
    {
      int nnz = A->colptr[j+1] - A->colptr[j];
      if (nnz < 0)
	{
	  bebop_log (0, "*** csc_matrix_sort_rowidx: at column %d, A->colptr[%d]="
		    "%d < A->colptr[%d]=%d ***\n", 
		    j, j+1, A->colptr[j+1], j, A->colptr[j]);
	  bebop_exit (EXIT_FAILURE);
	}
      max_nnz = MAX(nnz, max_nnz);
    }
  bebop_log (2, "it is %d\n", max_nnz);

  /* Sort the row indices in each column, reordering the corresponding values accordingly. */
  if (A->value_type == REAL)
    {
      /* col: Workspace for sorting (index,value) pairs */
      struct index_real_value_pair_t* col = bebop_calloc (max_nnz, sizeof (struct index_real_value_pair_t));

      bebop_log (2, "value_type==REAL\n");

      for (j = 0; j < A->n; j++)
	{
	  int nnz = A->colptr[j+1] - A->colptr[j];

	  if (nnz < 0)
	    {
	      bebop_log (0, "*** csc_matrix_sort_rowidx: At column %d, "
			"nnz = %d < 0 ***\n", j, nnz);
	      bebop_exit (EXIT_FAILURE);
	    }

	  /* Copy the (index,value) pairs into the temp location "col". */
	  copy_col2pairs (A, j, col, max_nnz);
	  /* Do the sorting in "col". */
	  qsort (col, nnz, sizeof (struct index_real_value_pair_t), compare_index_real_value_pairs);
	  /* Copy the (index,value) pairs back from "col" into the matrix. */
	  copy_pairs2col (col, max_nnz, A, j);
	}

      /* Free up the temp storage. */
      bebop_free (col);
    }
  else if (A->value_type == COMPLEX)
    {
      struct index_complex_value_pair_t* col = bebop_calloc (max_nnz, sizeof (struct index_complex_value_pair_t));
      bebop_log (2, "value_type==COMPLEX\n");

      for (j = 0; j < A->n; j++)
	{
	  int nnz = A->colptr[j+1] - A->colptr[j];

	  if (nnz < 0)
	    {
	      bebop_log (0, "*** csc_matrix_sort_rowidx: At column %d, "
			"nnz = %d < 0 ***\n", j, nnz);
	      bebop_exit (EXIT_FAILURE);
	    }

	  copy_col2pairs (A, j, col, max_nnz);
	  qsort (col, nnz, sizeof (struct index_real_value_pair_t), compare_index_complex_value_pairs);
	  copy_pairs2col (col, max_nnz, A, j);
	}

      bebop_free (col);
    }
  else if (A->value_type == PATTERN)
    {
      int* col = bebop_calloc (max_nnz, sizeof (int));

      bebop_log (2, "value_type==PATTERN\n");
      for (j = 0; j < A->n; j++)
	{
	  int nnz = A->colptr[j+1] - A->colptr[j];

	  if (nnz < 0)
	    {
	      bebop_log (0, "*** csc_matrix_sort_rowidx: At column %d, "
			"nnz = %d < 0 ***\n", j, nnz);
	      bebop_exit (EXIT_FAILURE);
	    }

	  copy_col2pairs (A, j, col, max_nnz);
	  qsort (col, nnz, sizeof (int), compare_ints);
	  copy_pairs2col (col, max_nnz, A, j);
	}

      bebop_free (col);
    }

  bebop_log (2, "=== Done with csc_matrix_sort_rowidx ===\n");
}



/**
 * Writes the given matrix in Harwell-Boeing format to the given file.
 *
 * @param filename [IN] The file to write
 * @param title    [IN] Title of the matrix
 * @param key      [IN] Key for the matrix
 * @param A        [IN] The matrix (in Harwell-Boeing format)
 *
 * @return Nonzero if something went wrong, else zero.
 */
int
write_harwell_boeing_mat_double (const char* filename, 
				 const char* title, const char* key,
				 const struct csc_matrix_t* A);

int
save_csc_matrix_in_harwell_boeing_format (const char* const filename, 
					  const struct csc_matrix_t* A)
{
  return write_harwell_boeing_mat_double (filename, "No label", "No key", A);
}

static void* 
malloc_and_copy (void* origarray, size_t nelt, size_t eltsize)
{
  void* newarray = bebop_malloc (nelt * eltsize);
  memcpy (newarray, origarray, nelt * eltsize);
  return newarray;
}


/*=====================================================================*/
void
pack_csc_matrix (struct csc_matrix_t* A, 
		 const int m, const int n, const int nnz, 
		 void* values, int* rowidx, int* colptr,
		 enum symmetry_type_t symmetry_type,
		 enum symmetric_storage_location_t symmetric_storage_location,
		 enum value_type_t value_type,
		 enum ownership_mode_t ownership,
		 void (*deallocator) (void*),
		 enum copy_mode_t copy_mode)
{
  bebop_log (2, "=== pack_csc_matrix ===\n");

  A->m = m;
  A->n = n;
  A->nnz = nnz;
  A->symmetry_type = symmetry_type;
  A->symmetric_storage_location = symmetric_storage_location;
  A->value_type = value_type;
  A->ownership = ownership;
  A->deallocator = deallocator == NULL ? &free : deallocator;
  if (copy_mode == NO_COPY)
    {
      A->values = values;
      A->rowidx = rowidx;
      A->colptr = colptr;
    }
  else
    {
      A->ownership = LIBRARY_DEALLOCATES;
      A->deallocator = &free;

      if (value_type == REAL)
	A->values = malloc_and_copy (values, nnz, sizeof (double));
      else if (value_type == COMPLEX)
	A->values = malloc_and_copy (values, nnz, sizeof (double_Complex));
      else
	A->values = NULL;

      A->rowidx = malloc_and_copy (rowidx, (m+1), sizeof (int));
      A->colptr = malloc_and_copy (colptr, nnz, sizeof (int));
    }

  bebop_log (2, "=== Done with pack_csc_matrix ===\n");
}


/*=====================================================================*/
void
init_csc_matrix (struct csc_matrix_t* A, 
		 const int m, const int n, const int nnz, 
		 void* values, int* rowidx, int* colptr,
		 enum symmetry_type_t symmetry_type,
		 enum symmetric_storage_location_t symmetric_storage_location,
		 enum value_type_t value_type,
		 enum ownership_mode_t ownership,
		 void (*deallocator) (void*),
		 enum copy_mode_t copy_mode)
{
  pack_csc_matrix (A, m, n, nnz, values, rowidx, colptr, symmetry_type, 
		   symmetric_storage_location, value_type, ownership,
		   deallocator, copy_mode);
}


/*========================================================================*/
int
valid_csc_matrix_p (const struct csc_matrix_t* A)
{
  const int m = A->m;
  const int n = A->n;
  const int nnz = A->nnz;
  const int* colptr = A->colptr;
  const int* rowidx = A->rowidx;
  int reached_nnz_in_colptr_p = 0;
  int i, j;

  bebop_log (2, "=== valid_csc_matrix_p ===\n");

  if (m < 1)
    {
      bebop_log (1, "*** valid_csc_matrix_p: Matrix has m = %d < 1 ***\n", m);
      return 0;
    }
  else if (n < 1)
    {
      bebop_log (1, "*** valid_csc_matrix_p: Matrix has n = %d < 1 ***\n", nnz);
      return 0;
    }
  else if (nnz < 0)
    {
      bebop_log (1, "*** valid_csc_matrix_p: Matrix has nnz = %d < 0 ***\n", nnz);
      return 0;
    }

  for (i = 0; i < n; i++)
    {
      const int cur = colptr[i];
      const int next = colptr[i+1];

      if (cur == nnz)
	{
	  if (! reached_nnz_in_colptr_p)
	    {
	      reached_nnz_in_colptr_p = 1;
	    }
	  if (next != cur)
	    {
	      bebop_log (1, "*** valid_csc_matrix_p: Although colptr[%d]==nnz==%d, "
		       "colptr[%d+1]==%d != colptr[%d] ***\n", 
		       i, nnz, i, next, i);
	      return 0;
	    }
	}

      /* Verify that colptr entries are in order */
      if (cur > next)
	{
	  bebop_log (1, "*** valid_csc_matrix_p: Matrix: colptr[%d] = %d > "
		   "colptr[%d] = %d ***\n", i, cur, i+1, next);
	  return 0;
	}

      /* Verify that current colptr entry doesn't exceed nnz */
      if (cur >= nnz) 
	{
	  bebop_log (1, "*** valid_csc_matrix_p: Matrix: At col %d, colptr[i] = "
		   "%d >= nnz = %d ***\n", i, colptr[i], nnz);
	  return 0;
	}

      /* Verify that rowidx entries are in range */
      for (j = cur; j < next; j++)
	{
	  if (rowidx[j] < 0 || rowidx[j] >= m)
	    {
	      bebop_log (1, "*** valid_csc_matrix_p: Matrix: at col %d, "
		       "rowidx[%d]=%d out of range [%d,%d) ***\n", 
		       i, j, rowidx[j], 0, m);
	      return 0;
	    }
	}
    }

  /* Verify that the last colptr entry is nnz */
  if (colptr[n] != nnz)
    {
      bebop_log (1, "*** valid_csc_matrix_p: colptr[n=%d] = %d != nnz=%d ***\n", 
	       n, colptr[n], nnz);
      return 0;
    }

  bebop_log (2, "=== Done with valid_csc_matrix_p ===\n");
  return 1;
}


/*=======================================================================*/
int
diag_csc_matrix (struct csc_matrix_t* A, 
		 const double* diag, const int m, const int n)
{
  const int nnz = m < n ? m : n;
  double *values;
  int i;

  bebop_log (2, "=== diag_csc_matrix ===\n");
  bebop_log (2, "\t%d x %d diagonal matrix, nnz = %d\n", m, n, nnz);

  A->nnz = nnz;
  A->m   = m;
  A->n   = n;

  A->colptr = (int*) bebop_calloc ((n+1), sizeof (int));
  A->rowidx = (int*) bebop_calloc (nnz, sizeof (int));
  values = (double*) bebop_calloc (nnz, sizeof (double));

  for (i = 0; i < nnz; i++)
    {
      A->colptr[i] = i;
      A->rowidx[i] = i;
      values[i] = diag[i];

      bebop_log (2, "\t\tA(%d,%d) = %g\n", i, i, diag[i]);
    }
  for (i = nnz; i <= n; i++)  /* in case m < n */
    {
      A->colptr[i] = nnz;
      bebop_log (2, "\t\tcolptr[%d] = nnz = %d\n", i, nnz);
    }

  A->symmetry_type = UNSYMMETRIC;
  A->value_type = REAL;
  A->values = (void*) values;
  A->ownership = LIBRARY_DEALLOCATES;
  A->deallocator = &free;
  bebop_log (2, "=== Done with diag_csc_matrix ===\n");
  return 0;
}


/*=====================================================================*/
int
print_csc_matrix_in_matrix_market_format (FILE* out, const struct csc_matrix_t* A)
{
  int start, end;
  int i, j;
  char value_type_label[20];
  char symmetry_type_label[20];

  bebop_log (2, "=== print_csc_matrix_in_matrix_market_format ===\n");

  if (A->value_type == REAL)
    strncpy (value_type_label, "real", 19);
  else if (A->value_type == COMPLEX)
    strncpy (value_type_label, "complex", 19);
  else if (A->value_type == PATTERN)
    strncpy (value_type_label, "pattern", 19);
  else 
    {
      bebop_log (0, "*** print_csc_matrix_in_matrix_market_format: Unsupported value type! ***\n");
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
      bebop_log (0, "*** print_csc_matrix_in_matrix_market_format: Unsupported symmetry type! ***\n");
      return -1;
    }

  fprintf (out, "%%%%MatrixMarket matrix coordinate %s %s\n", 
	   value_type_label, symmetry_type_label);

  if (bebop_debug_level() > 1)
    {
      fprintf (out, "%% colptr[%d]: ", A->n + 1);
      for (i = 0; i <= A->n; i++)
	fprintf (out, " %d", A->colptr[i]);

      fprintf (out, "\n%% rowidx[%d]: ", A->nnz);
      for (i = 0; i < A->nnz; i++)
	fprintf (out, " %d", A->rowidx[i]);

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
		fprintf (out, " %g+I*%g", 
			 double_Complex_real_part(values[i]), 
			 double_Complex_imag_part(values[i]));
	    }
	}

      fprintf (out, "\n");
    }

  fprintf (out, "%d %d %d\n", A->m, A->n, A->nnz);

  if (A->value_type == REAL)
    {
      const double* const values = (const double* const) (A->values);

      for (j = 0; j < A->n; j++)
	{
	  start = A->colptr[j];
	  end   = A->colptr[j+1];

	  /* MatrixMarket files use 1-based indices. */
	  for (i = start; i < end; i++)
	    fprintf (out, "%d %d %.13e\n", A->rowidx[i] + 1, j + 1, values[i]);
	}
    }
  else if (A->value_type == COMPLEX)
    {
      const double_Complex* const values = (const double_Complex* const) (A->values);

      for (j = 0; j < A->n; j++)
	{
	  start = A->colptr[j];
	  end   = A->colptr[j+1];

	  /* MatrixMarket files use 1-based indices. */
	  for (i = start; i < end; i++)
	    fprintf (out, "%d %d %.13e %.13e\n", A->rowidx[i] + 1, j + 1, double_Complex_real_part(values[i]), double_Complex_imag_part(values[i]));
	}
    }
  else if (A->value_type == PATTERN)
    {
      for (j = 0; j < A->n; j++)
	{
	  start = A->colptr[j];
	  end   = A->colptr[j+1];

	  /* MatrixMarket files use 1-based indices. */
	  for (i = start; i < end; i++)
	    fprintf (out, "%d %d\n", A->rowidx[i] + 1, j + 1);
	}
    }

  bebop_log (2, "=== Done with print_csc_matrix_in_matrix_market_format ===\n");
  return 0;
}


/*========================================================================*/
int
same_structure_csc_matrix (const struct csc_matrix_t* A, 
			   const struct csc_matrix_t* B)
{
  int i, j, m, n, nnz;

  bebop_log (2, "=== same_structure_csc_matrix ===\n");

  if (A->m != B->m)
    {
      bebop_log (1, "*** A->m = %d != B->m = %d\n ***\n", A->m, B->m);
      return 0;
    }
  m = A->m;

  if (A->n != B->n)
    {
      bebop_log (1, "*** A->n = %d != B->n = %d\n ***\n", A->n, B->n);
      return 0;
    }
  n = A->n;

  if (A->nnz != B->nnz)
    {
      bebop_log (1, "*** A->nnz = %d != B->nnz = %d\n ***\n", A->nnz, B->nnz);
      return 0;
    }
  nnz = A->nnz;

  /*
   * Make sure that the matrices have their last colptr element set to nnz.
   */
  if (A->colptr[n] != nnz)
    {
      bebop_log (0, "*** same_structure_csc_matrix: CSC sparse matrix data "
		"structures invalid:\n" 
		"\t{A,B}->colptr[n = %d] = %d != nnz = %d\n",
		n, A->colptr[n], nnz);
      bebop_exit (EXIT_FAILURE);
    }

  for (i = 0; i < n; i++)
    {
      int a_start = A->colptr[i];
      int a_end   = A->colptr[i+1];
      int b_start = B->colptr[i];
      int b_end   = B->colptr[i+1];

      if (a_start != b_start)
	{
	  bebop_log (1, "*** At column %d: A->colptr[%d] = %d "
		    "!= B->colptr[%d] = %d ***\n", 
		    i, i, A->colptr[i], i, B->colptr[i]);
	  return 0;
	}
      if (a_end != b_end)
	{
	  bebop_log (1, "*** At column %d: A->colptr[%d] = %d "
		    "!= B->colptr[%d] = %d ***\n", 
		    i, i+1, A->colptr[i+1], i+1, B->colptr[i+1]);
	  return 0;
	}
      for (j = a_start; j < a_end; j++)
	{
	  if (A->rowidx[j] != B->rowidx[j])
	    {
	      bebop_log (1, "*** At column %d:  A->rowidx[%d] "
			"= %d != B->rowidx[%d] = %d ***\n",
			i, j, A->rowidx[j], j, B->rowidx[j]);
	      return 0;
	    }
	}
    }
  if (A->colptr[n] != B->colptr[n])
    {
      bebop_log (1, "*** A->colptr[n = %d] = %d != B->colptr[n"
		" = %d] = %d ***\n", 
		n, A->colptr[n], n, B->colptr[n]);
      return 0;
    }


  bebop_log (2, "@@@ Matrix structures of A and B match @@@\n");
  bebop_log (2, "=== Done with same_structure_csc_matrix ===\n");
  return 1;
}


/*======================================================================*/
void
dealloc_csc_matrix (struct csc_matrix_t* A)
{
  bebop_log (2, "=== dealloc_csc_matrix ===\n");

  if (A == NULL)
    return;

  if (A->ownership == USER_DEALLOCATES)
    return;

  if (A->deallocator == NULL)
    {
      if (A->values != NULL && A->nnz > 0) 
	{
	  bebop_free (A->values);
	  A->values = NULL;
	}

      if (A->rowidx != NULL && A->nnz > 0)
	{
	  bebop_free (A->rowidx);
	  A->rowidx = NULL;
	}

      if (A->colptr != NULL && A->n > 0)
	{
	  bebop_free (A->colptr);
	  A->colptr = NULL;
	}
    }
  else
    {
      if (A->values != NULL && A->nnz > 0) 
	{
	  (A->deallocator) (A->values);
	  A->values = NULL;
	}

      if (A->rowidx != NULL && A->nnz > 0)
	{
	  (A->deallocator) (A->rowidx);
	  A->rowidx = NULL;
	}

      if (A->colptr != NULL && A->n > 0)
	{
	  (A->deallocator) (A->colptr);
	  A->colptr = NULL;
	}
    }

  A->m = 0;
  A->n = 0;
  A->nnz = 0;

  bebop_log (2, "=== Done with dealloc_csc_matrix ===\n");
}


/*======================================================================*/
void
copy_csc_matrix (struct csc_matrix_t* dest, const struct csc_matrix_t* src)
{
  const int m = src->m;
  const int n = src->n;
  const int nnz = src->nnz;

  bebop_log (2, "=== copy_csc_matrix ===\n");
  dest->m   = m;
  dest->n   = n;
  dest->nnz = nnz;

  if (src->value_type == REAL)
    {
      dest->values = bebop_malloc (nnz * sizeof (double));
      memcpy (dest->values, src->values, nnz * sizeof (double));
    }
  else if (src->value_type == COMPLEX)
    {
      dest->values = bebop_malloc (nnz * sizeof (double_Complex));
      memcpy (dest->values, src->values, nnz * sizeof (double_Complex));
    }
  else if (src->value_type == PATTERN)
    dest->values = NULL;

  dest->rowidx = bebop_malloc (nnz * sizeof (int));
  dest->colptr = bebop_malloc ((n+1) * sizeof (int));

  memcpy (dest->rowidx, src->rowidx, nnz * sizeof (int));
  memcpy (dest->colptr, src->colptr, (n+1) * sizeof (int));

  dest->symmetry_type = src->symmetry_type;
  dest->symmetric_storage_location = src->symmetric_storage_location;
  dest->value_type = src->value_type;
  dest->ownership = LIBRARY_DEALLOCATES;
  dest->deallocator = &free;

  bebop_log (2, "=== Done with copy_csc_matrix ===\n");
}


/*======================================================================*/
void
unpack_csc_matrix (const struct csc_matrix_t* A,
		   int* m, int* n, int* nnz,
		   void** values, int** rowidx, int** colptr,
		   enum symmetry_type_t* symmetry_type,
		   enum symmetric_storage_location_t* symmetric_storage_location,
		   enum value_type_t* value_type)
{
  bebop_log (2, "=== unpack_csc_matrix ===\n");

  *m = A->m;
  *n = A->n;
  *nnz = A->nnz;

  *values = A->values;
  *rowidx = A->rowidx;
  *colptr = A->colptr;

  *symmetry_type = A->symmetry_type;
  *symmetric_storage_location = A->symmetric_storage_location;
  *value_type = A->value_type;

  bebop_log (2, "=== Done with unpack_csc_matrix ===\n");
}


/*======================================================================*/
void
spy_csc_matrix (FILE* out, const struct csc_matrix_t* A)
{
  struct sparse_vector_t* rowvec;
  int m           = A->m; 
  int n           = A->n;
  int     *rowind = A->rowidx;
  int     *colptr = A->colptr;
  int which_row, j, k;


  bebop_log (2, "=== spy_csc_matrix ===\n");

  rowvec = create_sparse_vector (0);

  if (A->colptr[n] == 0)
    {
      bebop_log (0, "*** spy_csc_matrix: colptr[n=%d]=0, which means matrix is empty! ***\n", n);
      return;
    }

  /* 
   * For each row, extract row from sparse matrix and print its nonzero
   * structure.  We can assume that the entries in the row are sorted by 
   * column index.  
   */

  for (which_row = 0; which_row < m; which_row++) /* for each row */
    {
      resize_sparse_vector (rowvec, 0);

      for (j = 0; j < n; j++)
	{
	  const int start = colptr[j];
	  const int end   = colptr[j+1];

	  for (k = start; k < end; k++)
	    {
	      int row_index = rowind[k];

	      if (row_index == which_row) 
		/* We just insert a dummy value instead of the actual nonzero value */
		append_to_sparse_vector (rowvec, 1.0, j);
	    } /* for loop over entries in column j */
	} /* for loop over columns */

      /*
       * Print nonzero structure of the current row.
       */

      if (n > 0)
	{
	  int curidx = 0;
	  j = 0;

	  while (j < n)
	    {
	      int index;

	      if (curidx < length_sparse_vector (rowvec))
		{
		  index = get_sparse_vector_index (rowvec, curidx);

		  if (index == j)
		    {
		      curidx++;
		      if (j < n-1)
			fprintf (out, "X ");
		      else
			fprintf (out, "X");
		    }
		  else
		    {
		      if (j < n-1)
			fprintf (out, "  ");
		      else
			fprintf (out, " ");
		    }
		  j++;
		}
	      else
		{
		  if (j < n-1)
		    fprintf (out, "  ");
		  else
		    fprintf (out, " ");

		  j++;
		}
	    }
	}
      fprintf (out, "]\n");
    } /* end of for loop over rows */

  destroy_sparse_vector (rowvec);
  bebop_log (2, "=== Done with spy_csc_matrix ===\n");
}



/*========================================================================*/
int
read_harwell_boeing_mat_double (const char* filename, struct csc_matrix_t* A)
{
  /* 
   * values: The H-B interface functions always put read-in values into arrays 
   *         of doubles.  In the complex case, the real and imaginary parts 
   *         are interlaced (in that order: a1_real, a1_imag, a2_real, a2_imag, 
   *         ...). 
   * type:   The three-letter "type" field of the Harwell-Boeing format.
   */
  int retcode;
  int m, n, nnz;
  int* colptr;
  int* rowidx;
  double* values; 
  char* type = NULL;
  int nrhs = 0;
  enum symmetry_type_t symmetry_type = -1; /* flag */
  enum symmetric_storage_location_t symmetric_storage_location = -1; /* flag */
  enum value_type_t value_type = REAL; /* Default */

  bebop_log (2, "=== read_harwell_boeing_mat_double ===\n");

  /* Figure out the type of the matrix.  "type" is malloc'd, so we need to 
   * free it afterwards.  We won't be using nrhs, even if right-hand sides
   * are supplied in the file. */
  retcode = readHB_info (filename, &m, &n, &nnz, &type, &nrhs);

  /* HB functions return 0 if something went wrong, which is non-standard
   * in the C world.  */
  if (retcode == 0)
    {
      bebop_log (0, "*** read_harwell_boeing_mat_double: Failed to read "
	       "matrix info from file \"%s\" ***\n", filename);
      return -1;
    }

  /* Make sure that the matrix is a type that we can handle. */ 

  if (type[2] != 'A' && type[2] != 'a')
    {
      bebop_log (0, "*** read_harwell_boeing_mat_double: we don\'t know "
	       "how to handle unassembled matrices! ***\n");
      free (type);
      return -1;
    }
  
  /* Use "type" to figure out what kind of symmetry the matrix exhibits.  */
  if (type[1] == 'U' || type[1] == 'u' || type[1] == 'R' || type[1] == 'r')
    symmetry_type = UNSYMMETRIC; 
  else if (type[1] == 'S' || type[1] == 's')
    symmetry_type = SYMMETRIC; 
  else if (type[1] == 'Z' || type[1] == 'z')
    symmetry_type = SKEW_SYMMETRIC; 
  else
    {
      bebop_log (0, "*** read_harwell_boeing_mat_double: we don\'t know ho"
	       "w to handle matrices of symmetry type \'%c\' ***\n", type[1]);
      free (type);
      return -1;
    }

  /* Use "type" to figure out the type of values stored in the matrix. */
  if (type[0] == 'R' || type[0] == 'r')
    value_type = REAL;
  else if (type[0] == 'C' || type[0] == 'c')
    value_type = COMPLEX;
  else if (type[0] == 'P' || type[0] == 'p')
    value_type = PATTERN;
  else 
    {
      bebop_log (0, "*** read_harwell_boeing_mat_double: file specifies unknown value type %c ***\n", type[0]);
      free (type);
      return -1;
    }

  /* If the matrix is symmetric (skew-symmetric, etc.), we have to read the 
   * matrix first before we can tell whether the nonzeros are stored in the 
   * upper or the lower triangle. */ 

  /* Read in the matrix. */
  retcode = readHB_newmat_double (filename, &m, &n, &nnz, 
				  &colptr, &rowidx, &values);
  /* HB functions return 0 if something went wrong, which is non-standard
   * in the C world.  */
  if (! retcode)
    {
      bebop_log (0, "*** read_harwell_boeing_mat_double: Failed to read "
	       "matrix from file \"%s\" ***\n", filename);
      free (type);
      return -1;
    }

  /* If symmetric (skew-symmetric, etc.), figure out if the nonzeros are 
   * stored in the upper or lower triangle.  For robustness, we check for 
   * correctness, i.e. if there are entries in both triangles (which means 
   * there is something wrong with the matrix file). */
  if (symmetry_type != UNSYMMETRIC)
    {
      int j, k;

      symmetric_storage_location = -1; /* flag for "not yet set" */

      for (j = 0; j < n; j++)
	{
	  const int start = colptr[j];
	  const int end   = colptr[j+1];

	  for (k = start; k < end; k++)
	    {
	      if (rowidx[k] < j)
		{
		  if (symmetric_storage_location == -1)
		    {
		      if (bebop_debug_level() > 1)
			bebop_log (0, "first off-diagonal entry is (%d,%d), which is "
				 "in the upper triangle\n", rowidx[k], j);
		      symmetric_storage_location = UPPER_TRIANGLE;
		    }
		  else
		    {
		      if (symmetric_storage_location != UPPER_TRIANGLE)
			{
			  bebop_log (0, "\n*** WARNING: read_harwell_boeing"
				    "_mat_double: matrix has a symmetric or "
				    "skew-symmetric symmetry type, but conta"
				    "ins entries in both the upper and lower"
				    " triangles! ***\n");
			  bebop_log (0, "Offending entry: (%d,%d) at k = %d\n", 
				    rowidx[k], j, k);
			}
		    }
		}
	      else if (rowidx[k] > j)
		{
		  if (symmetric_storage_location == -1)
		    {
		      if (bebop_debug_level() > 1)
			bebop_log (0, "first off-diagonal entry is (%d,%d),"
				  " which is in the lower triangle\n", 
				  rowidx[k], j);
		      symmetric_storage_location = LOWER_TRIANGLE;
		    }
		  else
		    {
		      if (symmetric_storage_location != LOWER_TRIANGLE)
			{
			  bebop_log (0, "*** WARNING: read_harwell_boeing"
				    "_mat_double: matrix has a symmetric"
				    " or skew-symmetric symmetry type, b"
				    "ut contains entries in both the"
				   " upper and lower triangles! ***\n");
			  bebop_log (0, "Offending entry: (%d,%d) at k = "
				    "%d\n", rowidx[k], j, k);
			}
		    }
		}
	    }
	}

      /* If the matrix is diagonal, then we arbitrarily set the storage 
       * location to be the upper triangle. */
      if (symmetric_storage_location == -1)
	symmetric_storage_location = UPPER_TRIANGLE; 
    }

  if (value_type == REAL)
    pack_csc_matrix (A, m, n, nnz, values, rowidx, colptr, symmetry_type, 
		     symmetric_storage_location, value_type, 
		     LIBRARY_DEALLOCATES, &free, NO_COPY);
  else if (value_type == COMPLEX)
    {
      /* We have to copy the (real,imag) interlaced values in the "values" 
       * array into an array of "double_Complex".  It could be that "double 
       * _Complex" is represented in the same way, but we can't count on that. */
      double_Complex* __values = bebop_malloc (nnz * sizeof (double_Complex));
      int i;

      for (i = 0; i < nnz; i++)
	__values[i] = new_double_Complex(values[2*i], values[2*i+1]);

      free (values);
      pack_csc_matrix (A, m, n, nnz, __values, rowidx, colptr, symmetry_type, 
		       symmetric_storage_location, value_type,
		       LIBRARY_DEALLOCATES, &free, NO_COPY);
    }
  else if (value_type == PATTERN)
    {
      pack_csc_matrix (A, m, n, nnz, NULL, rowidx, colptr, symmetry_type, 
		       symmetric_storage_location, value_type,
		       LIBRARY_DEALLOCATES, &free, NO_COPY);
    }

  free (type);
  bebop_log (2, "=== Done with read_harwell_boeing_mat_double ===\n");
  return 0;
}


/*========================================================================*/
int
write_harwell_boeing_mat_double (const char* filename, 
				 const char* title, const char* key,
				 const struct csc_matrix_t* A)
{
  const char* rhstype = "F";  /* Full storage (we aren't using RHS anyway */
  char type[5]; /* An extra byte "just in case" */
  int retcode = 0;

  bebop_log (2, "=== write_harwell_boeing_mat_double ===\n");

  /* Default type: Real, Rectangular, Assembled.  "Rectangular" means m != n in general. */
  strncpy (type, "RRA", 4);

  if (A->value_type == REAL)
    type[0] = 'R';
  else if (A->value_type == COMPLEX)
    type[0] = 'C';
  else if (A->value_type == PATTERN)
    type[0] = 'P';
  else 
    {
      bebop_log (0, "*** write_harwell_boeing_mat_double: unsupported value type! ***\n");
      return -1;
    }

  if (A->m == A->n)
    {
      if (A->symmetry_type == SYMMETRIC)
	type[1] = 'S'; 
      else if (A->symmetry_type == SKEW_SYMMETRIC)
	type[1] = 'Z';
      else if (A->symmetry_type == HERMITIAN)
	type[1] = 'H';
      else
	type[1] = 'U'; /* Unsymmetric */
    }

  if (A->value_type == REAL)
    {
      const double* const values = (const double* const) (A->values);
      retcode = writeHB_mat_double (filename, 
				    A->m, A->n, A->nnz, 
				    A->colptr, A->rowidx, values,
				    0, NULL, NULL, NULL, 
				    title, key, type, 
				    NULL, NULL, NULL, NULL,
				    rhstype);
    }
  else if (A->value_type == COMPLEX)
    {
      /* The Harwell-Boeing routines want interlaced values */
      int k;
      const double_Complex* const values = (const double_Complex* const) (A->values);
      double *interlaced_values = bebop_malloc (A->nnz * 2 * sizeof(double));
      for (k = 0; k < A->nnz; k++)
	{
	  interlaced_values[2*k] = double_Complex_real_part (values[k]);
	  interlaced_values[2*k + 1] = double_Complex_imag_part (values[k]);
	}
      retcode = writeHB_mat_double (filename, 
				    A->m, A->n, A->nnz, 
				    A->colptr, A->rowidx, 
				    interlaced_values,
				    0, NULL, NULL, NULL, 
				    title, key, type, 
				    NULL, NULL, NULL, NULL,
				    rhstype);
      bebop_free (interlaced_values);
    }
  else if (A->value_type == PATTERN)
    {
      retcode = writeHB_mat_double (filename, 
				    A->m, A->n, A->nnz, 
				    A->colptr, A->rowidx, NULL,
				    0, NULL, NULL, NULL, 
				    title, key, type, 
				    NULL, NULL, NULL, NULL,
				    rhstype);
    }

  /* 
   * HB functions return 0 if something went wrong, which is non-standard
   * in the C world.
   */
  if (! retcode)
    {
      bebop_log (0, "*** ERROR: write_harwell_boeing_mat_double: Failed "
		"to write matrix to file %s ***\n", filename);
      return -1;
    }

  /* else */
  bebop_log (2, "=== Done with write_harwell_boeing_mat_double ===\n");
  return 0;
}


/*======================================================================*/
struct csc_matrix_t*
create_csc_matrix (const int m, const int n, const int nnz, 
		   void* values, int* rowidx, int* colptr,
		   enum symmetry_type_t symmetry_type,
		   enum symmetric_storage_location_t symmetric_storage_location,
		   enum value_type_t value_type,
		   enum ownership_mode_t ownership,
		   void (*deallocator) (void*),
		   enum copy_mode_t copy_mode)
{
  struct csc_matrix_t *A = NULL;

  bebop_log (2, "=== create_csc_matrix ===\n");
  A = bebop_malloc (sizeof (struct csc_matrix_t));
  init_csc_matrix (A, m, n, nnz, values, rowidx, colptr, symmetry_type, 
		   symmetric_storage_location, value_type, ownership,
		   deallocator, copy_mode);
  bebop_log (2, "=== Done with create_csc_matrix ===\n");
  return A;
}


/*======================================================================*/
void
destroy_csc_matrix (struct csc_matrix_t* A)
{
  bebop_log (2, "=== destroy_csc_matrix ===\n");
  if (A != NULL)
    {
      dealloc_csc_matrix (A);
      bebop_free (A);
    }
  else
    {
      bebop_log (0, "*** destroy_csc_matrix: A is NULL! ***\n");
    }
  bebop_log (2, "=== Done with destroy_csc_matrix ===\n");
}



int 
save_csc_matrix_in_matrix_market_format (const char* const filename, 
					 const struct csc_matrix_t* A)
{
  FILE* out = NULL;
  int errcode = 0;

  bebop_log (2, "=== save_csc_matrix_in_matrix_market_format ===\n");
  out = fopen (filename, "w");
  if (out == NULL)
    {
      bebop_log (0, "*** save_csc_matrix_in_matrix_market_format: failed "
	       "to open output file \"%s\" ***\n", filename);
      return -1;
    }

  errcode = print_csc_matrix_in_matrix_market_format (out, A);
  if (0 != fclose (out))
    {
      bebop_log (0, "*** save_csc_matrix_in_matrix_market_format: failed "
	       "to close output file \"%s\" ***\n", filename);
      return -1;
    }
  bebop_log (2, "=== Done with save_csc_matrix_in_matrix_market_format ===\n");
  return errcode;
}




int 
csc_matrix_expand_symmetric_storage (struct csc_matrix_t* A) 
{
  /* Code borrowed from Rich Vuduc's Sparsity-0.1 suite (spmv_util.c,
   * spmv_expand_symmcsr) and adapted for CSC-format matrices (the original was
   * for CSR-format matrices). */

  int *cur_col_nnz = NULL; /* # of nonzeros currently in each col of the matrix A */
  int *new_col_nnz = NULL; /* # of nonzeros in each col of the expanded matrix */
  int  new_nnz;            /* Total number of nonzeros in the expanded matrix */
  int j;                   /* Current column index */

  int *colptr = NULL;      /* Will be the new colptr array of A */
  int *rowidx = NULL;      /* Will be the new rowidx array of A */

  bebop_log (2, "=== csc_matrix_expand_symmetric_storage ===\n");
  if (A == NULL)
    {
      bebop_log (0, "*** csc_matrix_expand_symmetric_storage: A is NULL! ***\n");
      return -1;
    }
  /* The valid_csc_matrix_p check takes a while, so only do it if we are in debug mode */
  else if (bebop_debug_level() > 0 && !valid_csc_matrix_p (A))
    {
      bebop_log (0, "*** csc_matrix_expand_symmetric_storage: A is invalid! ***\n");
      return -1;
    }
  else if (A->m != A->n)
    {
      bebop_log (0, "*** csc_matrix_expand_symmetric_storage: A is not square! ***\n");
      return -2;
    }
  else if (A->symmetry_type == UNSYMMETRIC)
    {
      bebop_log (0, "*** csc_matrix_expand_symmetric_storage: "
		"A is already stored in an unsymmetric format! ***\n");
      /* Not an error -- it just means that we don't have any work to do */
      return 0;
    }

  bebop_log (2, "%d x %d matrix with %d nonzeros in symmetric storage\n", 
	    A->m, A->n, A->nnz);

  cur_col_nnz = (int*) bebop_calloc (A->n, sizeof (int));
  new_col_nnz = (int*) bebop_calloc (A->n, sizeof (int));

  /* 
   * Scan A and count how many new nonzeros we will need to create. 
   * When we are done scanning:
   *
   * cur_col_nnz[j] == # of nonzeros in column j of original matrix
   * new_col_nnz[j] == # of nonzeros to be stored in col j of matrix
   * new_nnz        == Total # of nonzeros to be stored in matrix
   */
  new_nnz = 0;
  /* First, count the number of nonzeros currently in each column. */
  bebop_log (2, "Counting # nnz currently in each column...");
  for (j = 0; j < A->n; j++)
    {
      /* j: index of the current column */
      cur_col_nnz[j] = A->colptr[j+1] - A->colptr[j];
      new_col_nnz[j] = cur_col_nnz[j];
      new_nnz += new_col_nnz[j];
    }
  bebop_log (2, "done.\n");

  if (bebop_debug_level() > 1)
    {
      for (j = 0; j < A->n; j++)
	{
	  if (cur_col_nnz[j] < 0)
	    {
	      bebop_log (0, "*** csc_matrix_expand_symmetric_storage: cur_"
		       "col_nnz[%d] = %d < 0 ***", j, cur_col_nnz[j]);
	      bebop_exit (EXIT_FAILURE);
	    }
	  else if (cur_col_nnz[j] > A->nnz)
	    {
	      bebop_log (0, "*** csc_matrix_expand_symmetric_storage: cur_co"
		       "l_nnz[%d] = %d > nnz = %d ***", j, cur_col_nnz[j], A->nnz);
	      bebop_exit (EXIT_FAILURE);
	    }
	  else if (A->colptr[j] > A->colptr[j+1])
	    {
	      bebop_log (0, "*** csc_matrix_expand_symmetric_storage: A->c"
		       "olptr[%d] = %d > A->colptr[%d] = %d ***\n", 
		       j, A->colptr[j], j+1, A->colptr[j+1]);
	      bebop_exit (EXIT_FAILURE);
	    }
	}
    }
  bebop_log (2, "Figuring out how many elts to reflect across diagonal...");
  /* Now figure out how many elements we need to reflect across the diagonal. */
  for (j = 0; j < A->n; j++)
    {
      int k;
      for (k = A->colptr[j]; k < A->colptr[j+1]; k++)
	{
	  int ii = A->rowidx[k];

	  /* Reflect off-diagonal elements across the diagonal; don't count
	   * diagonal elements twice.  Element (ii,j) goes into column ii (it
	   * is reflected into element (j,ii)). */
	  if (ii != j)
	    {
	      new_col_nnz[ii]++;
	      new_nnz++;
	    }
	}
    }
  bebop_log (2, "done.\n");

  /* Initialize new arrays.  These will hold the expanded version of the matrix. */
  bebop_log (2, "Initializing new arrays and cloning colptr...");
  colptr = (int*) bebop_malloc ((A->n + 1) * sizeof (int));
  rowidx = (int*) bebop_malloc (new_nnz * sizeof (int));

  /* Copy old colptr into new colptr */
  memcpy (colptr, A->colptr, (A->n + 1) * sizeof (int));
  bebop_log (2, "done.\n");

  /* 
   * Initialize column pointers in A.  After we are done:
   *
   * colptr initialized to the correct, final values.
   * new_col_nnz[j] reset to be equal to cur_col_nnz[j].
   */
  bebop_log (2, "Initializing new colptr...");
  for (j = 1; j <= A->n; j++)
    {
      colptr[j] = colptr[j-1] + new_col_nnz[j-1];
      new_col_nnz[j-1] = cur_col_nnz[j-1];
    }
  colptr[A->n] = new_nnz;
  bebop_log (2, "done.\n");

  if (bebop_debug_level() > 1)
    {
      for (j = 0; j < A->n; j++)
	{
	  if (colptr[j] < 0)
	    {
	      bebop_log (0, "*** csc_matrix_expand_symmetric_storage: upda"
			"ted colptr[%d] = %d < 0 ***\n", j, colptr[j]);
	      bebop_exit (EXIT_FAILURE);
	    }
	  else if (colptr[j] < A->colptr[j])
	    {
	      bebop_log (0, "*** csc_matrix_expand_symmetric_storage: colp"
		       "tr[%d] = %d < A->colptr[%d] = %d ***\n", j, colptr[j], 
		       j, A->colptr[j]);
	      bebop_exit (EXIT_FAILURE);
	    }
	  else if (colptr[j] > new_nnz)
	    {
	      bebop_log (0, "*** csc_matrix_expand_symmetric_storage: colp"
		       "tr[%d] = %d > new_nnz = %d ***\n", j, colptr[j], new_nnz);
	      bebop_exit (EXIT_FAILURE);
	    }
	  else if (colptr[j+1] - colptr[j] < 0)
	    {
	      bebop_log (0, "*** csc_matrix_expand_symmetric_storage: colp"
		       "tr[%d] - colptr[%d] = %d < 0 ***\n", j+1, j, colptr[j+1] - colptr[j]);
	      bebop_exit (EXIT_FAILURE);
	    }
	}
    }

  if (bebop_debug_level() > 1)
    {
      for (j = 0; j < A->n; j++)
	{
	  if (colptr[j+1] - colptr[j] < 0)
	    {
	      bebop_log (0, "*** csc_matrix_expand_symmetric_storage: colp"
		       "tr[%d] - colptr[%d] = %d < 0, in which colptr[%d]=%d "
		       "and colptr[%d]=%d ***\n", j+1, j, colptr[j+1] - colptr[j], 
		       j, colptr[j], j+1, colptr[j+1]);
	      bebop_exit (EXIT_FAILURE);
	    }
	}
    }

  /* 
   * Complete expansion of the matrix to full storage.  After we are done:
   *
   * (colptr, rowidx, values) is the full-storage equivalent of A.
   * new_col_nnz[j] == # of nonzeros in col j of the (expanded) matrix.
   */
  bebop_log (2, "Completing expansion into full storage...");
  if (A->value_type == REAL)
    {
      double* values = (double*) bebop_malloc (new_nnz * sizeof (double));
      for (j = 0; j < A->n; j++)
	{
	  int cur_nnz = cur_col_nnz[j]; /* number of nonzeros in current row of old matrix */
	  int k_cur = A->colptr[j];    /* current position in old matrix */
	  int k_new = colptr[j];       /* current position in expanded matrix */

	  if (bebop_debug_level() > 1)
	    {
	      int failure = 0;
	      assert (k_cur >= 0);
	      assert (k_new >= 0);
	      if (k_cur >= A->nnz)
		{
		  bebop_log (0, "*** csc_matrix_expand_symmetric_storage: "
			   "for column %d, k_cur = %d >= A->nnz = %d ***\n", 
			   j, k_cur, A->nnz);
		  failure = 1;
		}
	      else if (k_new >= new_nnz)
		{
		  bebop_log (0, "*** csc_matrix_expand_symmetric_storage: "
			   "for column %d, k_new = %d >= new_nnz = %d ***\n", 
			   j, k_new, new_nnz);
		  failure = 1;
		}
	      if (failure)
		bebop_exit (EXIT_FAILURE);
	    }

	  /* Copy current nonzeros from old matrix to new matrix */
	  memcpy (&rowidx[k_new], &(A->rowidx[k_cur]), cur_nnz * sizeof (int));
	  memcpy (values + k_new, ((double*) (A->values)) + k_cur, cur_nnz * sizeof (double));

	  /* Fill in the symmetric "missing" values */
	  while (k_cur < A->colptr[j+1])
	    {
	      /* Nonzero of A */
	      int ii = A->rowidx[k_cur];

	      if (ii != j) /* If not a diagonal element */
		{
		  /* Get the current element from the old matrix */
		  double a = ((double*) A->values)[k_cur];

		  /* Position of this transposed element in A */
		  k_new = colptr[ii] + new_col_nnz[ii];

		  /* Store the nonzero */
		  rowidx[k_new] = j;
		  if (A->symmetry_type == SYMMETRIC)
		    values[k_new] = a;
		  else if (A->symmetry_type == SKEW_SYMMETRIC)
		    values[k_new] = -a;
		  else if (A->symmetry_type == HERMITIAN)
		    /* Strictly speaking, a real symmetric matrix is also Hermitian. */
		    values[k_new] = a; 

		  /* Update so that the next element stored at column ii will
		   * appear at the right place */
		  new_col_nnz[ii]++;
		}

	      k_cur++;
	    }
	}

      /* Cleanup */
      bebop_free (cur_col_nnz);
      bebop_free (new_col_nnz);

      /* Free the old arrays and assign the new ones */
      bebop_free (A->colptr);
      bebop_free (A->rowidx);
      bebop_free (A->values);

      A->colptr = colptr;
      A->rowidx = rowidx;
      A->values = values;
      A->nnz = new_nnz;
    }
  else if (A->value_type == COMPLEX)
    {
      double_Complex* values = (double_Complex*) bebop_malloc (new_nnz * sizeof (double_Complex));
      for (j = 0; j < A->n; j++)
	{
	  int cur_nnz = cur_col_nnz[j]; /* number of nonzeros in current row of old matrix */
	  int k_cur = A->colptr[j];    /* current position in old matrix */
	  int k_new = colptr[j];       /* current position in expanded matrix */

	  /* Copy current nonzeros from old matrix to new matrix */
	  memcpy (&rowidx[k_new], &(A->rowidx[k_cur]), cur_nnz * sizeof (int));
	  memcpy (values + k_new, ((double_Complex*) A->values) + k_cur, cur_nnz * sizeof (double_Complex));

	  /* Fill in the symmetric "missing" values */
	  while (k_cur < A->colptr[j+1])
	    {
	      /* Nonzero of A */
	      int ii = A->rowidx[k_cur];

	      if (ii != j) /* If not a diagonal element */
		{
		  /* Get current element from old matrix */
		  double_Complex a = ((double_Complex*) A->values)[k_cur];

		  /* Position of this transposed element in A */
		  k_new = colptr[ii] + new_col_nnz[ii];

		  /* Store the nonzero */
		  rowidx[k_new] = j;
		  if (A->symmetry_type == SYMMETRIC)
		    values[k_new] = a;
		  else if (A->symmetry_type == SKEW_SYMMETRIC)
		    values[k_new] = double_Complex_negate(a);
		  else if (A->symmetry_type == HERMITIAN)
		    values[k_new] = double_Complex_conj(a);

		  /* Update so that the next element stored at column ii will
		   * appear at the right place */
		  new_col_nnz[ii]++;
		}

	      k_cur++;
	    }
	}

      /* Cleanup */
      bebop_free (cur_col_nnz);
      bebop_free (new_col_nnz);

      /* Free the old arrays and assign the new ones */
      bebop_free (A->colptr);
      bebop_free (A->rowidx);
      bebop_free (A->values);

      A->colptr = colptr;
      A->rowidx = rowidx;
      A->values = values;
      A->nnz = new_nnz;
    }
  else if (A->value_type == PATTERN)
    {
      for (j = 0; j < A->n; j++)
	{
	  int cur_nnz = cur_col_nnz[j];
	  int k_cur = A->colptr[j];    /* current position in old matrix */
	  int k_new = colptr[j];       /* current position in expanded matrix */

	  /* Copy current nonzeros from old matrix to new matrix */
	  /* mfh 1 Mar 2006:  BUG!!!  cur_nnz is too many values to copy 
	     for the current j, right? */
	  memcpy (&rowidx[k_new], &(A->rowidx[k_cur]), cur_nnz * sizeof (int));

	  /* Fill in the symmetric "missing" values */
	  while (k_cur < A->colptr[j+1])
	    {
	      /* Nonzero of A */
	      int ii = A->rowidx[k_cur];

	      if (ii != j) /* If not a diagonal element */
		{
		  /* Position of this transposed element in A */
		  k_new = colptr[ii] + new_col_nnz[ii];

		  /* Store the nonzero */
		  rowidx[k_new] = j;

		  /* Update so that the next element stored at column ii will
		   * appear at the right place */
		  new_col_nnz[ii]++;
		}

	      k_cur++;
	    }
	}

      /* Cleanup */
      bebop_free (cur_col_nnz);
      bebop_free (new_col_nnz);

      /* Free the old arrays and assign the new ones */
      bebop_free (A->colptr);
      bebop_free (A->rowidx);

      A->colptr = colptr;
      A->rowidx = rowidx;
      A->values = NULL;
      A->nnz = new_nnz;
    }
  bebop_log (2, "done.\n");

  if (bebop_debug_level() > 1)
    {
      for (j = 0; j < A->n; j++)
	{
	  if (colptr[j+1] - colptr[j] < 0)
	    {
	      bebop_log (0, "*** csc_matrix_expand_symmetric_storage: colp"
		       "tr[%d] - colptr[%d] = %d < 0, in which colptr[%d]=%d "
		       "and colptr[%d]=%d ***\n", 
		       j+1, j, colptr[j+1] - colptr[j], j, colptr[j], j+1, 
		       colptr[j+1]);
	      bebop_exit (EXIT_FAILURE);
	    }
	}
    }
  /* Sort the row indices in the matrix */
  bebop_log (2, "\tsorting row indices...");
  csc_matrix_sort_rowidx (A);
  bebop_log (2, "done.\n");

  /* Now A uses unsymmetric storage */
  A->symmetry_type = UNSYMMETRIC;

  bebop_log (2, "=== Done with csc_matrix_expand_symmetric_storage ===\n");
  return 0;
}


int
print_csc_matrix_in_matlab_format (FILE* out, struct csc_matrix_t* A)
{
  int errcode = 0;
  struct coo_matrix_t* B = NULL;

  bebop_log (2, "=== print_csc_matrix_in_matlab_format ===\n");
  
  bebop_log (2, "Converting matrix to COO format\n");
  B = csc_to_coo (A);
  if (B == NULL)
    {
      bebop_log (0, "*** print_csc_matrix_in_matlab_format: Failed "
		"to convert CSC-format sparse matrix into COO form"
		"at for printing! ***\n");
      return -1;
    }
  bebop_log (2, "Printing COO matrix in Matlab format\n");
  errcode = print_coo_matrix_in_matlab_format (out, B);
  destroy_coo_matrix (B);
  bebop_log (2, "=== Done with print_csc_matrix_in_matlab_format ===\n");
  return errcode;
}


int
save_csc_matrix_in_matlab_format (const char* const filename, struct csc_matrix_t* A)
{
  FILE* out = NULL;
  int errcode = 0;

  bebop_log (2, "=== save_csc_matrix_in_matlab_format ===\n");
  out = fopen (filename, "w");
  if (out == NULL)
    {
      bebop_log (0, "*** save_csc_matrix_in_matlab_format: failed "
	       "to open output file \"%s\" ***\n", filename);
      bebop_log (2, "=== Done with save_csc_matrix_in_matlab_format ===\n");
      return -1;
    }

  errcode = print_csc_matrix_in_matlab_format (out, A);
  if (0 != fclose (out))
    {
      bebop_log (0, "*** save_csc_matrix_in_matlab_format: failed "
	       "to close output file \"%s\" ***\n", filename);
      bebop_log (2, "=== Done with save_csc_matrix_in_matlab_format ===\n");
      return -1;
    }
  bebop_log (2, "=== Done with save_csc_matrix_in_matlab_format ===\n");
  return errcode;
}


