/**
 * @file read_mm.c
 * @author Mark Hoemmen
 * @since 10 Jun 2004
 * @date Time-stamp: <2009-05-16 15:42:45 mhoemmen>
 * 
 * Functions for loading and converting MatrixMarket - format sparse matrix 
 * files.
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
#include <bebop/smc/coo_matrix.h>
#include <bebop/smc/coord_elem.h>
#include <bebop/smc/csc_matrix.h>
#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/mmio.h>
#include <bebop/smc/read_mm.h>

#include <bebop/util/complex.h>
#include <bebop/util/merge_sort.h>
#include <bebop/util/malloc.h>
#include <bebop/util/util.h>
#include <bebop/util/log.h>

#include <assert.h>
#include <stdio.h>
#include <string.h>


/**
 * Compares two coordinate-array elements by their row indices, and returns 
 * values like strcmp does.  Comparison function for sorting.
 */
int
compare_coord_elem_by_row_real (const void* a, const void* b)
{
  struct coord_elem_t
  {
    int r;
    int c;
    double val;
  };

  const struct coord_elem_t* x = ((struct coord_elem_t*) a);
  const struct coord_elem_t* y = ((struct coord_elem_t*) b);

  if (x->r < y->r) 
    return -1;
  else if (x->r > y->r)
    return +1;

  /* else */
  return 0;
}

int
compare_coord_elem_by_row_complex (const void* a, const void* b)
{
  struct coord_elem_t
  {
    int r;
    int c;
    double_Complex val;
  };

  const struct coord_elem_t* x = ((struct coord_elem_t*) a);
  const struct coord_elem_t* y = ((struct coord_elem_t*) b);

  if (x->r < y->r) 
    return -1;
  else if (x->r > y->r)
    return +1;

  /* else */
  return 0;
}

int
compare_coord_elem_by_row_pattern (const void* a, const void* b)
{
  struct coord_elem_t
  {
    int r;
    int c;
  };

  const struct coord_elem_t* x = ((struct coord_elem_t*) a);
  const struct coord_elem_t* y = ((struct coord_elem_t*) b);

  if (x->r < y->r) 
    return -1;
  else if (x->r > y->r)
    return +1;

  /* else */
  return 0;
}


/**
 * Compares two coordinate-array elements by their column indices, and 
 * returns values like strcmp does.  Comparison function for sorting.
 */
int
compare_coord_elem_by_col_real (const void* a, const void* b)
{
  struct coord_elem_t
  {
    int r;
    int c;
    double val;
  };

  const struct coord_elem_t* x = ((struct coord_elem_t*) a);
  const struct coord_elem_t* y = ((struct coord_elem_t*) b);

  if (x->c < y->c) 
    return -1;
  else if (x->c > y->c)
    return +1;

  /* else */
  return 0;
}

int
compare_coord_elem_by_col_complex (const void* a, const void* b)
{
  struct coord_elem_t
  {
    int r;
    int c;
    double_Complex val;
  };

  const struct coord_elem_t* x = ((struct coord_elem_t*) a);
  const struct coord_elem_t* y = ((struct coord_elem_t*) b);

  if (x->c < y->c) 
    return -1;
  else if (x->c > y->c)
    return +1;

  /* else */
  return 0;
}

int
compare_coord_elem_by_col_pattern (const void* a, const void* b)
{
  struct coord_elem_t
  {
    int r;
    int c;
  };

  const struct coord_elem_t* x = ((struct coord_elem_t*) a);
  const struct coord_elem_t* y = ((struct coord_elem_t*) b);

  if (x->c < y->c) 
    return -1;
  else if (x->c > y->c)
    return +1;

  /* else */
  return 0;
}

/*========================================================================*/
void
sort_coord_elem_array_for_csr_conversion (void* coord_array, 
					  const int length,
					  enum value_type_t value_type)
{ 
  bebop_log (2, "=== sort_coord_elem_array_for_csr_conversion ===\n");

  if (value_type == REAL)
    {
      struct coord_elem_t
      {
	int r;
	int c;
	double val;
      };
      merge_sort (coord_array, length, sizeof (struct coord_elem_t), 
		  compare_coord_elem_by_col_real);
      merge_sort (coord_array, length, sizeof (struct coord_elem_t), 
		  compare_coord_elem_by_row_real);
    }
  else if (value_type == COMPLEX)
    {
      struct coord_elem_t
      {
	int r;
	int c;
	double_Complex val;
      };
      merge_sort (coord_array, length, sizeof (struct coord_elem_t), 
		  compare_coord_elem_by_col_complex);
      merge_sort (coord_array, length, sizeof (struct coord_elem_t), 
		  compare_coord_elem_by_row_complex);
    }
  else if (value_type == PATTERN)
    {
      struct coord_elem_t
      {
	int r;
	int c;
      };
      merge_sort (coord_array, length, sizeof (struct coord_elem_t), 
		  compare_coord_elem_by_col_pattern);
      merge_sort (coord_array, length, sizeof (struct coord_elem_t), 
		  compare_coord_elem_by_row_pattern);
    }
  bebop_log (2, "=== Done with sort_coord_elem_array_for_csr_conversion ===\n");
}



/*========================================================================*/
void
sort_coord_elem_array_for_csc_conversion (void* coord_array, 
					  const int length,
					  enum value_type_t value_type)
{ 
  bebop_log (2, "=== sort_coord_elem_array_for_csc_conversion ===\n");

  if (value_type == REAL)
    {
      struct coord_elem_t
      {
	int r;
	int c;
	double val;
      };
      merge_sort (coord_array, length, sizeof (struct coord_elem_t), 
		  compare_coord_elem_by_row_real);
      merge_sort (coord_array, length, sizeof (struct coord_elem_t), 
		  compare_coord_elem_by_col_real);
    }
  else if (value_type == COMPLEX)
    {
      struct coord_elem_t
      {
	int r;
	int c;
	double_Complex val;
      };
      merge_sort (coord_array, length, sizeof (struct coord_elem_t), 
		  compare_coord_elem_by_row_complex);
      merge_sort (coord_array, length, sizeof (struct coord_elem_t), 
		  compare_coord_elem_by_col_complex);
    }
  else if (value_type == PATTERN)
    {
      struct coord_elem_t
      {
	int r;
	int c;
      };
      merge_sort (coord_array, length, sizeof (struct coord_elem_t), 
		  compare_coord_elem_by_row_pattern);
      merge_sort (coord_array, length, sizeof (struct coord_elem_t), 
		  compare_coord_elem_by_col_pattern);
    }
  bebop_log (2, "=== Done with sort_coord_elem_array_for_csc_conversion ===\n");
}


/*========================================================================*/
void
coo_matrix_to_coord_elem_array (void** p_coord_array, 
				int *p_length, 
				const struct coo_matrix_t* A)
{
  const int nnz = A->nnz;
  int k;

  bebop_log (2, "=== coo_matrix_to_coord_elem_array ===\n");
  *p_length = nnz;

  if (A->value_type == REAL)
    {
      struct coord_elem_t
      {
	int r;
	int c;
	double val;
      };
      struct coord_elem_t* coord_array = bebop_malloc (nnz * sizeof (struct coord_elem_t));
      double* values = (double*) (A->val);

      for (k = 0; k < nnz; k++)
	{
	  coord_array[k].r = A->II[k];
	  coord_array[k].c = A->JJ[k];
	  coord_array[k].val = values[k];
	}
      *p_coord_array = coord_array;
    }
  else if (A->value_type == COMPLEX)
    {
      struct coord_elem_t
      {
	int r;
	int c;
	double_Complex val;
      };
      struct coord_elem_t* coord_array = bebop_malloc (nnz * sizeof (struct coord_elem_t));
      double_Complex* values = (double_Complex*) (A->val);

      for (k = 0; k < nnz; k++)
	{
	  coord_array[k].r = A->II[k];
	  coord_array[k].c = A->JJ[k];
	  coord_array[k].val = values[k];
	}
      *p_coord_array = coord_array;
    }
  else if (A->value_type == PATTERN)
    {
      struct coord_elem_t
      {
	int r;
	int c;
      };
      struct coord_elem_t* coord_array = bebop_malloc (nnz * sizeof (struct coord_elem_t));

      for (k = 0; k < nnz; k++)
	{
	  coord_array[k].r = A->II[k];
	  coord_array[k].c = A->JJ[k];
	}
      *p_coord_array = coord_array;
    }

  bebop_log (2,"=== Done with coo_matrix_to_coord_elem_array ===\n");
}



/*========================================================================*/
void
coo_to_csc_matrix (struct csc_matrix_t* A, const struct coo_matrix_t* B)
{
  int m   = B->m;
  int n   = B->n;
  int nnz = B->nnz;
  enum symmetry_type_t symmetry_type = B->symmetry_type;
  enum symmetric_storage_location_t symmetric_storage_location = B->symmetric_storage_location;
  int    *colptr;
  int    *rowidx;
  void   *__values = NULL;
  void   *__coord_array;
  int i, j, curcol = 0;
  int index_base = B->index_base;

  bebop_log (2, "=== coo_to_csc_matrix ===\n");
  bebop_log (2, "\tm = %d, n = %d, nnz = %d\n", m, n, nnz);


  if (B->value_type == REAL)
    __values = bebop_calloc (nnz, sizeof (double));
  else if (B->value_type == COMPLEX)
    __values = bebop_calloc (nnz, sizeof (double_Complex));
  else if (B->value_type == PATTERN)
    __values = NULL;
  else 
    {
      bebop_log (0, "*** coo_to_csc_matrix: input matrix in COO format has"
		" invalid value type %d ***\n", B->value_type);
      bebop_exit (EXIT_FAILURE);
    }

  colptr = bebop_calloc ((n+1), sizeof (int));
  rowidx = bebop_calloc (nnz,   sizeof (int));

  if (nnz == 0) 
    {
      /* calloc fills in colptr with zeros, which is all 
         we need if the matrix is empty */
      bebop_log (2, "\tMatrix is empty, done\n");
      init_csc_matrix (A, m, n, nnz, __values, rowidx, colptr, UNSYMMETRIC, 
		       UPPER_TRIANGLE, B->value_type, LIBRARY_DEALLOCATES,
		       &free, NO_COPY);
      bebop_log (2, "=== Done with coo_to_csc_matrix ===\n");
      return; 
    }

  /* Intermediate conversion to coordinate array, so we can sort the entries */
  bebop_log (2, "\tIntermediate conversion to coordinate array format\n");
  coo_matrix_to_coord_elem_array (&__coord_array, &nnz, B);

  bebop_log (2, "\tSorting elements of coordinate array\n");
  sort_coord_elem_array_for_csc_conversion (__coord_array, nnz, B->value_type);

  if (B->value_type == REAL)
    {
      struct coord_array_t
      {
	int r;
	int c;
	double val;
      };
      struct coord_array_t* coord_array = (struct coord_array_t*) __coord_array;
      double* values = (double*) __values;

      /*
       * Having sorted the elements of the coordinate array, first by initial 
       * column, then by initial row, now the first column with nonzeros is
       * coord_array[0].c.
       */
      curcol = coord_array[0].c - index_base;
      for (i = 0; i <= curcol; i++)
	{
	  /* 
	   * Until we get to first column with an entry in it, all the colptrs 
	   * before then are zero.  The colptr for that first column is also zero. 
	   */
	  colptr[i] = 0;  
	}
      /* For each coord_array entry i, add it to the CSC matrix */
      for (i = 0; i < nnz; i++)
	{
	  if (coord_array[i].c - index_base> curcol)
	    {
	      /* 
	       * We may jump more than one column at a time, so set the colptr 
	       * entries for the empty columns in between.
	       */
	      for (j = curcol+1; j <= coord_array[i].c - index_base; j++) 
		colptr[j] = i;

	      curcol = coord_array[i].c - index_base;
	    }

	  values[i] = coord_array[i].val;   
	  rowidx[i] = coord_array[i].r - index_base;
	}

      /* Set the last entries in colptr appropriately */
      for (j = curcol+1; j <= n; j++)
	colptr[j] = nnz;

      init_csc_matrix (A, m, n, nnz, __values, rowidx, colptr, symmetry_type, 
		       symmetric_storage_location, B->value_type, 
		       LIBRARY_DEALLOCATES, &free, NO_COPY);
    }
  else if (B->value_type == COMPLEX)
    {
      struct coord_array_t
      {
	int r;
	int c;
	double_Complex val;
      };
      struct coord_array_t* coord_array = (struct coord_array_t*) __coord_array;
      double_Complex* values = (double_Complex*) __values;

      /*
       * Having sorted the elements of the coordinate array, first by initial 
       * column, then by initial row, now the first column with nonzeros is
       * coord_array[0].c.
       */
      curcol = coord_array[0].c - index_base;
      for (i = 0; i <= curcol; i++)
	{
	  /* 
	   * Until we get to first column with an entry in it, all the colptrs 
	   * before then are zero.  The colptr for that first column is also zero. 
	   */
	  colptr[i] = 0;  
	}
      /* For each coord_array entry i, add it to the CSC matrix */
      for (i = 0; i < nnz; i++)
	{
	  if (coord_array[i].c - index_base > curcol)
	    {
	      /* 
	       * We may jump more than one column at a time, so set the colptr 
	       * entries for the empty columns in between.
	       */
	      for (j = curcol+1; j <= coord_array[i].c - index_base; j++) 
		colptr[j] = i;

	      curcol = coord_array[i].c - index_base;
	    }

	  values[i] = coord_array[i].val;   
	  rowidx[i] = coord_array[i].r - index_base;
	}

      /* Set the last entries in colptr appropriately */
      for (j = curcol+1; j <= n; j++)
	colptr[j] = nnz;

      init_csc_matrix (A, m, n, nnz, __values, rowidx, colptr, symmetry_type, 
		       symmetric_storage_location, B->value_type, 
		       LIBRARY_DEALLOCATES, &free, NO_COPY);
    }
  else if (B->value_type == PATTERN)
    {
      struct coord_array_t
      {
	int r;
	int c;
      };
      struct coord_array_t* coord_array = (struct coord_array_t*) __coord_array;

      /*
       * Having sorted the elements of the coordinate array, first by initial 
       * column, then by initial row, now the first column with nonzeros is
       * coord_array[0].c.
       */
      curcol = coord_array[0].c - index_base;
      for (i = 0; i <= curcol; i++)
	{
	  /* 
	   * Until we get to first column with an entry in it, all the colptrs 
	   * before then are zero.  The colptr for that first column is also zero. 
	   */
	  colptr[i] = 0;  
	}
      /* For each coord_array entry i, add it to the CSC matrix */
      for (i = 0; i < nnz; i++)
	{
	  if (coord_array[i].c - index_base > curcol)
	    {
	      /* 
	       * We may jump more than one column at a time, so set the colptr 
	       * entries for the empty columns in between.
	       */
	      for (j = curcol+1; j <= coord_array[i].c - index_base; j++) 
		colptr[j] = i;

	      curcol = coord_array[i].c - index_base;
	    }

	  rowidx[i] = coord_array[i].r - index_base;
	}

      /* Set the last entries in colptr appropriately */
      for (j = curcol+1; j <= n; j++)
	colptr[j] = nnz;

      init_csc_matrix (A, m, n, nnz, NULL, rowidx, colptr, symmetry_type, 
		       symmetric_storage_location, B->value_type, 
		       LIBRARY_DEALLOCATES, &free, NO_COPY);
    }

  bebop_free (__coord_array);
  bebop_log (2, "=== Done with coo_to_csc_matrix ===\n");
}


int
read_matrix_market_sparse (const char* filename, struct coo_matrix_t* A)
{
  int ret_code;
  char _matcode[4];
  FILE *f;
  int M, N, nnz;   
  int *II; 
  int *JJ;
  void *val;
  int i;
  enum value_type_t value_type;
  /* Set to nonzero if the MatrixMarket file specifies that symmetric storage is used. */
  int symmetric_p = 0;
  enum symmetry_type_t symmetry_type = UNSYMMETRIC;
  enum symmetric_storage_location_t symmetric_storage_location = -1;

  char* matcode = (char*) _matcode;

  bebop_log (2, "=== read_matrix_market_sparse ===\n");

  bebop_log (2, "\tOpening file\n");
  f = fopen (filename, "r");
  if (f == NULL)
    {
      bebop_log (0, "*** read_matrix_market_sparse: Failed to open Mat"
		"rixMarket file %s ***\n", filename);
      return MM_COULD_NOT_READ_FILE;
    }

  bebop_log (2, "\tReading MatrixMarket banner\n");
  if (mm_read_banner (f, &matcode) != 0)
    {
      bebop_log (0, "*** read_matrix_market_sparse: Could not process "
		"MatrixMarket banner ***\n");
      return MM_NO_HEADER;
    }

  if (! mm_is_matrix (matcode))
    {
      bebop_log (0, "*** read_matrix_market_sparse: Matrix isn\'t a "
		"matrix!!! ***\n" );
      return MM_UNSUPPORTED_TYPE;
    }
  else if (! mm_is_sparse (matcode))
    {
      bebop_log (0, "*** read_matrix_market_sparse: Matrix is dense, "
		"but we only want sparse matrices ***\n" );
      return MM_UNSUPPORTED_TYPE;
    }

  if (mm_is_real (matcode))
    value_type = REAL;
  else if (mm_is_complex (matcode))
    value_type = COMPLEX;
  else if (mm_is_pattern (matcode))
    value_type = PATTERN;
  else
    {
      bebop_log (0, "*** read_matrix_market_sparse: Unsupported value type ***\n");
      return MM_UNSUPPORTED_TYPE;
    }

  if (mm_is_general (matcode))
    symmetry_type = UNSYMMETRIC;
  else 
    {
      if (mm_is_symmetric (matcode))
	{
	  symmetry_type = SYMMETRIC;
	  symmetric_p = 1;
	}
      else if (mm_is_skew (matcode))
	{
	  symmetry_type = SKEW_SYMMETRIC;
	  symmetric_p = 1;
	}
      else 
	{
	  bebop_log (0, "*** read_matrix_market_real_sparse: We can only handle "
		   "general, symmetric or skew-symmetric matrices; some other kind "
		   "of symmetry is present here ***\n");
	  return MM_UNSUPPORTED_TYPE;
	}
    }

  /* 
   * Find out size of sparse matrix.  If the matrix is symmetric, then only 
   * the nonzeros on one side of the matrix are shown.  If we wanted to know
   * how many nonzeros there really are, we would need to examine all of them
   * -- in particular, we would need to know how many nonzeros are along the 
   * diagonal (the diagonal may not be full).  But we don't care, since we
   * are only storing the symmetric part.  What we do need to do is figure 
   * out on what side the nonzeros are stored, whether in the upper triangle
   * or in the lower triangle.
   */
  bebop_log (2, "\tReading sparse matrix size...");
  ret_code = mm_read_mtx_crd_size (f, &M, &N, &nnz);
  if (ret_code != 0)
    {
      bebop_log (0, "\n*** read_matrix_market_sparse: Failed to read "
	       "matrix size ***\n" );
      return -1;
    }
  bebop_log (2, "M = %d, N = %d, nnz = %d\n", M, N, nnz);

  if (M < 0 || N < 0 || nnz < 0)
    {
      bebop_log (0, "\n*** read_matrix_market_real_sparse: Matrix dimensions or "
		"nnz out of range:  m=%d, n=%d, nnz=%d ***\n", M, N, nnz);
      return -2;
    }

  /* 
   * Reserve memory for matrix.  
   */
  bebop_log (2, "\tAllocating memory for matrix\n");
  II = (int *) bebop_malloc (nnz * sizeof(int));
  JJ = (int *) bebop_malloc (nnz * sizeof(int));
  if (value_type == REAL)
    val = (double *) bebop_malloc (nnz * sizeof(double));
  else if (value_type == COMPLEX)
    val = (double_Complex *) bebop_malloc (nnz * sizeof(double_Complex));
  else
    val = NULL;

  /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
  /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
  /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

  bebop_log (2, "\tReading matrix entries from file\n");
  if (value_type == REAL)
    {
      for (i = 0; i < nnz; i++)
	{
	  double x;
	  double* __val = (double*) (val);

	  fscanf (f, "%d %d %lg\n", &II[i], &JJ[i], &x);
	  __val[i] = x;
	  II[i]--;  /* adjust from 1-based to 0-based */
	  JJ[i]--;
	}
    }
  else if (value_type == COMPLEX)
    {
      for (i = 0; i < nnz; i++)
	{
	  double real, imag;
	  double_Complex* __val = (double_Complex*) (val);

	  fscanf (f, "%d %d %lg %lg\n", &II[i], &JJ[i], &real, &imag);
	  __val[i] = new_double_Complex(real, imag);
	  II[i]--;  /* adjust from 1-based to 0-based */
	  JJ[i]--;
	}
    }
  else if (value_type == PATTERN)
    {
      for (i = 0; i < nnz; i++)
	{
	  fscanf (f, "%d %d", &II[i], &JJ[i]);
	  II[i]--;  /* adjust from 1-based to 0-based */
	  JJ[i]--;
	}
    }

  if (symmetric_p)
    {
      /* Figure out whether the nonzeros are stored in the upper triangle or 
       * the lower triangle.  We assume correctness here.  TODO:  For 
       * robustness, we should have a flag that, if set, requires that all 
       * the stored nonzeros be checked to make sure that they are in the 
       * same triangle. */
      for (i = 0; i < nnz; i++)
	{
	  if (II[i] < JJ[i])
	    {
	      symmetric_storage_location = UPPER_TRIANGLE;
	      break;
	    }
	  else if (II[i] > JJ[i])
	    {
	      symmetric_storage_location = LOWER_TRIANGLE;
	      break;
	    }
	}

      /* If the matrix is diagonal (i.e. if i reaches nnz), then we arbitrarily
       * decide that the nonzeros are in the upper triangle (which they are,
       * though they are also in the lower triangle). */
      if (i == nnz)
	symmetric_storage_location = UPPER_TRIANGLE;
    }

  if (f != stdin) 
    {
      bebop_log (2, "\tClosing file\n");
      if (fclose (f))
	{
	  bebop_log (0, "*** read_matrix_market_sparse: Failed to close file ***\n");
	  bebop_free (II);
	  bebop_free (JJ);
	  bebop_free (val);
	  A->II = NULL;
	  A->JJ = NULL;
	  A->val = NULL;
	  A->m = 0;
	  A->n = 0;
	  A->nnz = 0;
	  return -3;
	}
    }

  A->II   = II;
  A->JJ   = JJ;
  A->val = val;
  A->m   = M;
  A->n   = N;
  A->nnz = nnz;
  A->symmetry_type = symmetry_type;
  A->symmetric_storage_location = symmetric_storage_location;
  A->value_type = value_type;
  A->index_base = 0;
  A->ownership = LIBRARY_DEALLOCATES;
  A->deallocator = &free;

  bebop_log (2, "=== Done with read_matrix_market_sparse ===\n");
  return 0;
}



/*========================================================================*/
int
read_matrix_market_real_sparse (const char* filename, struct coo_matrix_t* A)
{
  int ret_code;
  char _matcode[4];
  FILE *f;
  int M, N, nnz;   
  int *II; 
  int *JJ;
  double *val;
  int i;
  /* Set to nonzero if the MatrixMarket file specifies that symmetric or 
   * skew-symmetric storage is used. */
  int symmetric_p = 0;
  enum symmetry_type_t symmetry_type = UNSYMMETRIC;
  enum symmetric_storage_location_t symmetric_storage_location = -1;
  enum value_type_t value_type = REAL; 

  char* matcode = (char*) _matcode;


  bebop_log (2, "=== read_matrix_market_real_sparse ===\n");

  bebop_log (2, "\tOpening file\n");
  f = fopen (filename, "r");
  if (f == NULL)
    {
      bebop_log (0, "*** read_matrix_market_real_sparse: Failed to open Mat"
	      "rixMarket file %s ***\n", filename);
      return MM_COULD_NOT_READ_FILE;
    }

  bebop_log (2,"\tReading MatrixMarket banner\n");
  if (mm_read_banner (f, &matcode) != 0)
    {
      bebop_log (0, "*** read_matrix_market_real_sparse: Could not process "
	      "MatrixMarket banner ***\n");
      return MM_NO_HEADER;
    }

  if (! mm_is_matrix (matcode))
    {
      bebop_log (0, "*** read_matrix_market_real_sparse: Matrix isn\'t a "
	   "matrix!!! ***\n" );
	  return MM_UNSUPPORTED_TYPE;
    }
  else if (! mm_is_sparse (matcode))
    {
      bebop_log (0, "*** read_matrix_market_real_sparse: Matrix is dense, "
	   "but we only want sparse matrices ***\n" );
      return MM_UNSUPPORTED_TYPE;
    }
  else if (! mm_is_real (matcode))
    {
      bebop_log (0, "*** read_matrix_market_real_sparse: We can only handle "
			   "real (i.e. non-complex) matrices ***\n");
	  return MM_UNSUPPORTED_TYPE;
    }

  if (mm_is_general (matcode))
    symmetry_type = UNSYMMETRIC;
  else 
    {
      if (mm_is_symmetric (matcode))
	{
	  symmetry_type = SYMMETRIC;
	  symmetric_p = 1;
	}
      else if (mm_is_skew (matcode))
	{
	  symmetry_type = SKEW_SYMMETRIC;
	  symmetric_p = 1;
	}
      else
	{
	  bebop_log (0, "*** read_matrix_market_real_sparse: We can only handle "
		   "general, symmetric or skew-symmetric matrices; some other kind "
		   "of symmetry is present here ***\n");
	  return MM_UNSUPPORTED_TYPE;
	}
    }

  /* 
   * Find out size of sparse matrix.  If the matrix is symmetric, then only 
   * the nonzeros on one side of the matrix are shown.  If we wanted to know
   * how many nonzeros there really are, we would need to examine all of them
   * -- in particular, we would need to know how many nonzeros are along the 
   * diagonal (the diagonal may not be full).  But we don't care, since we
   * are only storing the symmetric part.  What we do need to do is figure 
   * out on what side the nonzeros are stored, whether in the upper triangle
   * or in the lower triangle.
   */
  bebop_log (2, "\tReading sparse matrix size...");
  ret_code = mm_read_mtx_crd_size (f, &M, &N, &nnz);
  if (ret_code != 0)
    {
      bebop_log (0, "\n*** read_matrix_market_real_sparse: Failed to read "
	       "matrix size ***\n" );
      return -1;
    }
  bebop_log (2, "M = %d, N = %d, nnz = %d\n", M, N, nnz);

  if (M < 0 || N < 0 || nnz < 0)
    {
      bebop_log (0, "\n*** read_matrix_market_real_sparse: Matrix dimensions or "
	       "nnz out of range:  m=%d, n=%d, nnz=%d ***\n", M, N, nnz);
      return -2;
    }

  /* 
   * Reserve memory for matrix.  
   */
  II = (int *) bebop_malloc (nnz * sizeof(int));
  JJ = (int *) bebop_malloc (nnz * sizeof(int));
  val = (double *) bebop_malloc (nnz * sizeof(double));

  /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
  /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
  /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

  bebop_log (2, "\tReading matrix entries from file\n");
  for (i = 0; i < nnz; i++)
    {
      double x;

      fscanf (f, "%d %d %lg\n", &II[i], &JJ[i], &x);
      /* This means val[i], but we have to go through contortions because val is a void*. */
      *((double*) val + i*sizeof(double)) = x;
      II[i]--;  /* adjust from 1-based to 0-based */
      JJ[i]--;
    }

  if (symmetric_p)
    {
      /* Figure out whether the nonzeros are stored in the upper triangle or 
       * the lower triangle.  We assume correctness here.  TODO:  For 
       * robustness, we should have a flag that, if set, requires that all 
       * the stored nonzeros be checked to make sure that they are in the 
       * same triangle. */
      for (i = 0; i < nnz; i++)
	{
	  if (II[i] < JJ[i])
	    {
	      symmetric_storage_location = UPPER_TRIANGLE;
	      break;
	    }
	  else if (II[i] > JJ[i])
	    {
	      symmetric_storage_location = LOWER_TRIANGLE;
	      break;
	    }
	}

      /* If the matrix is diagonal (i.e. if i reaches nnz), then we arbitrarily
       * decide that the nonzeros are in the upper triangle (which they are,
       * though they are also in the lower triangle). */
      if (i == nnz)
	symmetric_storage_location = UPPER_TRIANGLE;
    }

  if (f != stdin) 
    {
      bebop_log (2, "\tClosing file\n");
      if (fclose (f))
	{
	  bebop_log (0, "\n*** read_matrix_market_real_sparse: Failed to close file ***\n");
	  bebop_free (II);
	  bebop_free (JJ);
	  bebop_free (val);
	  A->II = NULL;
	  A->JJ = NULL;
	  A->val = NULL;
	  A->m = 0;
	  A->n = 0;
	  A->nnz = 0;
	  return -3;
	}
    }

  A->II   = II;
  A->JJ   = JJ;
  A->val = val;
  A->m   = M;
  A->n   = N;
  A->nnz = nnz;
  A->symmetry_type = symmetry_type;
  A->symmetric_storage_location = symmetric_storage_location;
  A->index_base = 0;
  A->value_type = value_type;
  A->ownership = LIBRARY_DEALLOCATES;
  A->deallocator = &free;

  bebop_log (2, "=== Done with read_matrix_market_real_sparse ===\n");
  return 0;
}


/*======================================================================*/
int
read_matrix_market_real_general_dense (const char* filename, 
				       int *m, int *n, double** A)
{
  int ret_code;
  char _matcode[4];
  FILE *f;
  int M, N, i;
  double* B;
  char* matcode = (char*) _matcode;

  bebop_log (2, "=== read_matrix_market_real_general_dense ===\n");

  bebop_log (2, "\tOpening file\n");
  f = fopen (filename, "r");
  if (f == NULL)
    {
      bebop_log (0, "*** read_matrix_market_real_general_dense: Failed to open "
	       "MatrixMarket file %s ***\n", filename);
      return MM_COULD_NOT_READ_FILE;
    }

  bebop_log (2, "\tReading MatrixMarket banner\n");
  if (mm_read_banner (f, &matcode) != 0)
    {
      bebop_log (0, "*** read_matrix_market_real_general_dense: Could not process "
	       "MatrixMarket banner ***\n");
      fclose (f);
      return MM_NO_HEADER;
    }

  if (! mm_is_matrix (matcode))
    {
      bebop_log (0, "*** read_matrix_market_real_general_dense: Object type is not "
		   "\'matrix\' ***\n");
      return MM_UNSUPPORTED_TYPE;
    }
  else if (! mm_is_dense (matcode))
    {
      bebop_log (0, "*** read_matrix_market_real_general_dense: Matrix is not dense ***\n");
      return MM_UNSUPPORTED_TYPE;
    }
  else if (! mm_is_real (matcode))
    {
      bebop_log (0, "*** read_matrix_market_real_general_dense: Matrix is not real ***\n");
      return MM_UNSUPPORTED_TYPE;
    }
  else if (! mm_is_general (matcode))
    {
      bebop_log (0, "*** read_matrix_market_real_general_dense: Matrix is not general ***\n");
      return MM_UNSUPPORTED_TYPE;
    }

  /* 
   * Find out size of dense matrix 
   */
  bebop_log (2, "\tReading dense matrix dimensions...");
  ret_code = mm_read_mtx_array_size (f, &M, &N);
  if (ret_code != 0)
    {
      bebop_log (0, "\n*** read_matrix_market_real_general_dense: Failed to read matrix dimensions ***\n");
      return -1;
    }
  bebop_log (2, "M = %d, N = %d\n", M, N);

  if (M < 0 || N < 0)
    {
      bebop_log (0, "\n*** read_matrix_market_real_general_dense: Matrix dimensions "
	       "out of valid range:  m=%d, n=%d ***\n", M, N);
      return -2;
    }

  /* 
   * Reserve memory for matrix
   */
  bebop_log (2, "\tAllocating memory for matrix\n");
  B = (double *) bebop_malloc (M * N * sizeof(double));

  /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
  /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
  /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

  bebop_log (2, "\tReading matrix entries from file\n");
  for (i = 0; i < M * N; i++)
    fscanf (f, "%lg\n", &B[i]);

  if (f != stdin) 
    {
      bebop_log (2, "\tClosing file\n");
      if (fclose (f))
	{
	  bebop_log (0, "*** read_matrix_market_real_general_dense: "
		    "Failed to close input file %s ***\n", filename);
	  /* Make sure that memory is deallocated and that the matrix's 
	   * dimensions are set to reasonable error-case values */
	  bebop_free (B);
	  *A = NULL;
	  *m = 0;
	  *n = 0;
	  return -3;
	}
    }

  *m = M;
  *n = N;
  *A = B;
  bebop_log (2, "=== Done with read_matrix_market_real_general_dense ===\n");
  return 0;
}



struct csc_matrix_t*
coo_to_csc (struct coo_matrix_t* A)
{
  struct csc_matrix_t* B = bebop_malloc (sizeof (struct csc_matrix_t));
  coo_to_csc_matrix (B, A);
  return B;
}



struct csr_matrix_t*
coo_to_csr (struct coo_matrix_t* A)
{
  int m   = A->m;
  int n   = A->n;
  int nnz = A->nnz;
  enum symmetry_type_t symmetry_type = A->symmetry_type;
  enum symmetric_storage_location_t symmetric_storage_location = A->symmetric_storage_location;
  enum value_type_t value_type = A->value_type;
  int    *rowptr = NULL;
  int    *colidx = NULL;
  void   *__values = NULL;
  void   *__coord_array = NULL;
  int i, j, currow = 0;
  struct csr_matrix_t* B = bebop_calloc (1, sizeof (struct csr_matrix_t));
  int index_base = A->index_base;

  /* bebop_set_debug_level (2); */

  bebop_log (2, "=== coo_to_csr ===\n");
  bebop_log (2, "\tm = %d, n = %d, nnz = %d, value_type = %d\n", m, n, nnz, value_type);

  if (A->value_type == REAL)
    __values = bebop_calloc (nnz, sizeof (double));
  else if (A->value_type == COMPLEX)
    __values = bebop_calloc (nnz, sizeof (double_Complex));
  else if (A->value_type == PATTERN)
    __values = NULL;
  else 
    {
      bebop_log (0, "*** coo_to_csr: input matrix in COO format has"
	       " invalid value type %d ***\n", A->value_type);
      /* 
       * FIXME: should do something smarter here, like marking A 
       * as "invalid" so that we don't try to do any other operations
       * with it.
       */
      bebop_exit (EXIT_FAILURE);
    }

  rowptr = bebop_calloc ((m+1), sizeof (int));
  colidx = bebop_calloc (nnz,   sizeof (int));

  if (nnz == 0) 
    {
      /* calloc fills in rowptr with zeros, which is all 
         we need if the matrix is empty */
      bebop_log (2, "\tMatrix is empty, done\n");
      init_csr_matrix (B, m, n, nnz, __values, colidx, rowptr, UNSYMMETRIC, 
		       UPPER_TRIANGLE, value_type, LIBRARY_DEALLOCATES,
		       &free, NO_COPY);
      bebop_log (2, "=== Done with coo_to_csr_matrix ===\n");
      return B; 
    }

  /* Intermediate conversion to coordinate array, so we can sort the entries */
  bebop_log (2, "\tIntermediate conversion to coordinate array format\n");
  coo_matrix_to_coord_elem_array (&__coord_array, &nnz, A);

  bebop_log (2, "\tSorting elements of coordinate array\n");
  sort_coord_elem_array_for_csr_conversion (__coord_array, nnz, value_type);

  if (value_type == REAL)
    {
      struct coord_array_t
      {
	int r;
	int c;
	double val;
      };
      struct coord_array_t* coord_array = (struct coord_array_t*) __coord_array;
      double* values = (double*) __values;

      bebop_log (2, "\tReal value specialization\n");

      /*
       * Having sorted the elements of the coordinate array, first by initial 
       * row, then within each row by initial column, now the first row with 
       * nonzeros is coord_array[0].r.
       */
      bebop_log (2, "\tInitializing start of rowptr\n");
      currow = coord_array[0].r - index_base;
      for (i = 0; i <= currow; i++)
	{
	  /* 
	   * Until we get to first row with an entry in it, all the rowptrs 
	   * before then are zero.  The rowptr for that first column is also 
	   * zero. 
	   */
	  rowptr[i] = 0;  
	}
      /* For each coord_array entry i, add it to the CSR matrix */
      bebop_log (2, "\tAdding entries to CSR matrix\n");
      for (i = 0; i < nnz; i++)
	{
	  bebop_log (3, "\t\ti = %d of %d\n", i, nnz);
	  if (coord_array[i].r - index_base > currow)
	    {
	      /* 
	       * We may jump more than one row at a time, so set the rowptr 
	       * entries for the empty rows in between.
	       */
	      for (j = currow+1; j <= coord_array[i].r - index_base; j++) 
		{
		  if (j - index_base < 0 || j - index_base > m)
		    {
		      bebop_log (0, "*** At entry %d of %d, "
				"j = %d is out of the valid range "
				"[%d,%d] ***\n", 
				i, nnz, j - index_base, 0, m);
		      bebop_free (values);
		      bebop_free (rowptr);
		      bebop_free (colidx);
		      bebop_free (B);
		      return NULL;
		    }
		  rowptr[j] = i;
		}

	      currow = coord_array[i].r - index_base;
	    }

	  values[i] = coord_array[i].val;   
	  colidx[i] = coord_array[i].c - index_base;
	}

      /* Set the last entries in rowptr appropriately */
      for (j = currow+1; j <= m; j++)
	rowptr[j] = nnz;

      init_csr_matrix (B, m, n, nnz, __values, colidx, rowptr, symmetry_type, 
		       symmetric_storage_location, value_type, 
		       LIBRARY_DEALLOCATES, &free, NO_COPY);
    }
  else if (value_type == COMPLEX)
    {
      struct coord_array_t
      {
	int r;
	int c;
	double_Complex val;
      };
      struct coord_array_t* coord_array = (struct coord_array_t*) __coord_array;
      double_Complex* values = (double_Complex*) __values;

      /*
       * Having sorted the elements of the coordinate array, first by initial 
       * row, then within each row by initial column, now the first row with nonzeros is
       * coord_array[0].r.
       */
      currow = coord_array[0].r - index_base;
      for (i = 0; i <= currow; i++)
	{
	  /* 
	   * Until we get to first row with an entry in it, all the rowptrs 
	   * before then are zero.  The rowptr for that first column is also zero. 
	   */
	  rowptr[i] = 0;  
	}
      /* For each coord_array entry i, add it to the CSR matrix */
      for (i = 0; i < nnz; i++)
	{
	  if (coord_array[i].r - index_base > currow)
	    {
	      /* 
	       * We may jump more than one row at a time, so set the rowptr
	       * entries for the empty rows in between.
	       */
	      for (j = currow+1; j <= coord_array[i].r - index_base; j++) 
		rowptr[j] = i;

	      currow = coord_array[i].r - index_base;
	    }

	  values[i] = coord_array[i].val;   
	  colidx[i] = coord_array[i].c - index_base;
	}

      /* Set the last entries in colptr appropriately */
      for (j = currow+1; j <= m; j++)
	rowptr[j] = nnz;

      init_csr_matrix (B, m, n, nnz, __values, colidx, rowptr, symmetry_type, 
		       symmetric_storage_location, value_type,
		       LIBRARY_DEALLOCATES, &free, NO_COPY);
    }
  else if (value_type == PATTERN)
    {
      struct coord_array_t
      {
	int r;
	int c;
      };
      struct coord_array_t* coord_array = (struct coord_array_t*) __coord_array;

      /*
       * Having sorted the elements of the coordinate array, first by initial 
       * row, then within each row by initial column, now the first row with nonzeros is
       * coord_array[0].r.
       */
      currow = coord_array[0].r - index_base;
      for (i = 0; i <= currow; i++)
	{
	  /* 
	   * Until we get to first row with an entry in it, all the rowptrs
	   * before then are zero.  The rowptr for that first row is also zero. 
	   */
	  rowptr[i] = 0;  
	}
      /* For each coord_array entry i, add it to the CSR matrix */
      for (i = 0; i < nnz; i++)
	{
	  if (coord_array[i].c - index_base > currow)
	    {
	      /* 
	       * We may jump more than one row at a time, so set the rowptr 
	       * entries for the empty rows in between.
	       */
	      for (j = currow+1; j <= coord_array[i].r - index_base; j++) 
		rowptr[j] = i;

	      currow = coord_array[i].r - index_base;
	    }

	  colidx[i] = coord_array[i].c - index_base;
	}

      /* Set the last entries in rowptr appropriately */
      for (j = currow+1; j <= m; j++)
	rowptr[j] = nnz;

      init_csr_matrix (B, m, n, nnz, NULL, colidx, rowptr, symmetry_type, 
		       symmetric_storage_location, value_type,
		       LIBRARY_DEALLOCATES, &free, NO_COPY);
    }

  bebop_free (__coord_array);
  bebop_log (2, "=== Done with coo_to_csr ===\n");
  return B;
}



struct coo_matrix_t*
csr_to_coo (struct csr_matrix_t* A)
{
  struct coo_matrix_t* B = bebop_calloc (1, sizeof (struct coo_matrix_t));
  int *II, *JJ;
  int i, k;

  bebop_log (2, "=== csr_to_coo ===\n");

  B->m = A->m;
  B->n = A->n;
  B->nnz = A->nnz;
  B->symmetry_type = A->symmetry_type;
  B->symmetric_storage_location = A->symmetric_storage_location;
  B->value_type = A->value_type;
  B->index_base = 0;

  II = bebop_malloc (A->nnz * sizeof (int));
  JJ = bebop_malloc (A->nnz * sizeof (int));

  if (A->value_type == REAL)
    {
      B->val = bebop_malloc (A->nnz * sizeof (double));
      memcpy (B->val, A->values, A->nnz * sizeof (double));
    }
  else if (A->value_type == COMPLEX)
    {
      B->val = bebop_malloc (A->nnz * sizeof (double_Complex));
      memcpy (B->val, A->values, A->nnz * sizeof (double_Complex));
    }
  else if (A->value_type == PATTERN)
    {
      B->val = NULL;
    }

  k = 0;
  for (i = 0; i < A->m; i++)
    {
      const int start = A->rowptr[i];
      const int end   = A->rowptr[i+1];

      for (k = start; k < end; k++)
	{
	  II[k] = i;
	  JJ[k] = A->colidx[k];
	}
    }

  B->II = II;
  B->JJ = JJ;
  B->ownership = LIBRARY_DEALLOCATES;
  B->deallocator = &free;

  bebop_log (2, "=== Done with csr_to_coo ===\n");
  return B;
}



struct coo_matrix_t*
csc_to_coo (struct csc_matrix_t* A)
{
  struct coo_matrix_t* B = bebop_calloc (1, sizeof (struct coo_matrix_t));
  int *II, *JJ;
  int j, k;

  bebop_log (2, "=== csc_to_coo ===\n");

  B->m = A->m;
  B->n = A->n;
  B->nnz = A->nnz;
  B->symmetry_type = A->symmetry_type;
  B->symmetric_storage_location = A->symmetric_storage_location;
  B->value_type = A->value_type;
  B->index_base = 0;

  II = bebop_malloc (A->nnz * sizeof (int));
  JJ = bebop_malloc (A->nnz * sizeof (int));

  bebop_log (2, "Transferring values...");
  if (A->value_type == REAL)
    {
      B->val = bebop_malloc (A->nnz * sizeof (double));
      memcpy (B->val, A->values, A->nnz * sizeof (double));
    }
  else if (A->value_type == COMPLEX)
    {
      B->val = bebop_malloc (A->nnz * sizeof (double_Complex));
      memcpy (B->val, A->values, A->nnz * sizeof (double_Complex));
    }
  else if (A->value_type == PATTERN)
    {
      B->val = NULL;
    }
  bebop_log (2, "done.\n");

  k = 0;
  for (j = 0; j < A->n; j++)
    {
      const int start = A->colptr[j];
      const int end   = A->colptr[j+1];

      for (k = start; k < end; k++)
	{
	  II[k] = A->rowidx[k];
	  JJ[k] = j;
	}
    }

  B->II = II;
  B->JJ = JJ;
  B->ownership = LIBRARY_DEALLOCATES;
  B->deallocator = &free;

  bebop_log (2, "=== Done with csc_to_coo ===\n");
  return B;
}

