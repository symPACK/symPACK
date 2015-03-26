/**
 * @file csr_triple_product.c
 * @author Mark Hoemmen
 * @since 19 Jun 2006
 * @date Time-stamp: <2008-07-16 11:17:06 mhoemmen>
 *
 * Implementation of the CSR sparse matrix matrix triple product
 * kernel used for algebraic multigrid.  Algorithm taken from the
 * Hypre source code (rap.c).
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
#include <bebop/smc/csr_triple_product.h>
#include <bebop/smc/csr_transpose.h>

#include <bebop/util/malloc.h>
#include <bebop/util/util.h>
#include <bebop/util/list.h>

#include <assert.h>


int
csr_matrix_triple_product_kernel_double_complex (const int m, const int n,
						 const double_Complex RT_val[], 
						 const int RT_ind[],
						 const int RT_ptr[],
						 const double_Complex A_val[],
						 const int A_ind[],
						 const int A_ptr[],
						 const double_Complex P_val[],
						 const int P_ind[],
						 const int P_ptr[],
						 double_Complex** val,
						 int** ind,
						 int** ptr)
{
  /* Not yet implemented */
  assert (0);
  return -1;
}

int
csr_matrix_triple_product_kernel_double_pattern (const int m, const int n,
						 const int RT_ind[],
						 const int RT_ptr[],
						 const int A_ind[],
						 const int A_ptr[],
						 const int P_ind[],
						 const int P_ptr[],
						 int** ind,
						 int** ptr)
{
  /* Not yet implemented */
  assert (0);
  return -1;
}



typedef struct {
  int j; /* column index */
  double v; /* value */
} dr_entry_t;

dr_entry_t*
new_dr_entry (int j, double v)
{
  dr_entry_t* e = bebop_malloc (sizeof (dr_entry_t));
  e->j = j;
  e->v = v;
  return e;
}

int
dr_entry_comp (void* __e1, void* __e2)
{
  const dr_entry_t* e1 = (dr_entry_t*) __e1;
  const dr_entry_t* e2 = (dr_entry_t*) __e2;
  const int j1 = e1->j;
  const int j2 = e2->j;

  if (j1 < j2)
    return -1;
  else if (j1 > j2)
    return +1;
  else 
    return 0;
}

int
dr_entry_equal (void* e1, void* e2)
{
  return (dr_entry_comp (e1, e2) == 0);
}


int
ptap_csr_dr (const int n_fine, 
	     const int n_coarse,
	     const double A_val[],
	     const int A_ind[],
	     const int A_ptr[],
	     const double P_val[],
	     const int P_ind[],
	     const int P_ptr[],
	     double** val,
	     int** ind,
	     int** ptr,
	     int* nnz)
{
  alist_t* entries = bebop_malloc (n_coarse * sizeof (alist_t));
  int i, j, l, k, ii, jj, ll;
  int* __ptr = bebop_malloc ((n_coarse + 1) * sizeof (int));
  int* __ind = NULL;
  double* __val = NULL;
  int count = 0;

  for (i = 0; i < n_coarse; i++)
    entries[i] = alist_create ();

  for (l = 0; l < n_fine; l++)
    {
      const int start = A_ptr[l];
      const int end = A_ptr[l+1];

      for (ll = start; ll < end; ll++)
	{
	  const double a_value = A_val[ll];
	  k = A_ind[ll];
	  
	  /* Look in rows l,k of P */
	  for (ii = P_ptr[l]; ii < P_ptr[l+1]; ii++)
	    {
	      const double AP_val_ii = a_value * P_val[ii];
	      alist_t L;

	      i = P_ind[ii];
	      L = entries[i];
	      for (jj = P_ptr[k]; jj < P_ptr[k+1]; jj++)
		{
		  list_node_t* node;
		  /* We don't allocate e on the heap yet because
		   * we may not need to make a new entry; we want
		   * to save heap allocation (consing) for when
		   * we really need it. */
		  dr_entry_t e;

		  j = P_ind[jj];
		  e.j = j;
		  e.v = AP_val_ii * P_val[jj];
		  node = alist_find (L, &e, &dr_entry_equal);
		  if (node == NULL)
		    L = alist_push_item (L, new_dr_entry (e.j, e.v));
		  else
		    {
		      dr_entry_t* existing = (dr_entry_t*) (node->item);
		      assert (existing->j == j);
		      existing->v += e.v;
		    }
		}
	      entries[i] = L;
	    }
	}
    }

  /* Count nnz_coarse and set up the ptr array */
  count = 0;
  __ptr[0] = 0;
  for (i = 0; i < n_coarse; i++)
    {
      alist_t L = entries[i];
      const int len = alist_length (L);
      
      L = alist_sort (L, &dr_entry_comp);
      __ptr[i+1] = count + len;
      count += len;

      entries[i] = L; /* make sure we keep the sorted list */
    }

  /* Allocate space for the result matrix */
  __val = bebop_malloc (count * sizeof (double));
  __ind = bebop_malloc (count * sizeof (int));

  /* Fill in the values and column indices of the sparse matrix */
  count = 0;
  for (i = 0; i < n_coarse; i++)
    {
      alist_t L = entries[i];
      while (! alist_empty_p (L))
	{
	  const dr_entry_t* e = (dr_entry_t*) alist_car (L);
	  __ind[count] = e->j;
	  __val[count] = e->v;
	  count++;
	}
      entries[i] = alist_destroy (L);
    }
  bebop_free (entries);

  *nnz = count;
  *val = __val;
  *ind = __ind;
  *ptr = __ptr;
  return 0;
}




int
csr_matrix_triple_product_kernel_double_real (const int m, const int n,
					      const double RT_val[], 
					      const int RT_ind[],
					      const int RT_ptr[],
					      const double A_val[],
					      const int A_ind[],
					      const int A_ptr[],
					      const double P_val[],
					      const int P_ind[],
					      const int P_ptr[],
					      double** val,
					      int** ind,
					      int** ptr)
{
  int errcode = 0;
  double* _val;
  int* _ptr;
  int* _ind;
  int nnz = 0;

  int* P_marker;
  int* A_marker;

  int ic, i, i1, i2, i3, jj1, jj2, jj3;

  int jj_counter, jj_row_beginning;
  int start_indexing = 0;  /* base of the indexing (zero-based or one-based) */
  
  int n_fine = n;
  int n_coarse = m;

  double r_entry, r_a_product, r_a_p_product;

   /*
    *  Copy RT into R so that we have row-wise access to restriction.
    */
  double* R_val = NULL;
  int* R_ind = NULL;
  int* R_ptr = NULL;
  errcode = csr_matrix_transpose_kernel_double_real (&R_val, &R_ind, &R_ptr, 
						     n_fine, n_coarse, 
						     RT_val, RT_ind, RT_ptr);
  if (errcode != 0)
    return errcode;
  
  /* 
   * Allocate space for R^T * A * P and marker arrays.
   */
  _ptr = bebop_malloc ((n_coarse + 1) * sizeof (int));
  P_marker = bebop_malloc (n_coarse * sizeof (int)); 
  A_marker = bebop_malloc (n_fine * sizeof (int)); 


  /* 
   * First pass: determine number of nonzeros in R^T * A * P (nnz) and set
   * up the ptr array for that matrix.
   */

  jj_counter = start_indexing;
  for (ic = 0; ic < n_coarse; ic++)
    P_marker[ic] = -1;
  for (i = 0; i < n_fine; i++)
    A_marker[i] = -1;

  /* Loop over the coarse points */
  for (ic = 0; ic < n_coarse; ic++)
    {
      /* Set marker for diagonal entry (ic,ic) of result matrix */
      
      P_marker[ic] = jj_counter;
      jj_row_beginning = jj_counter;
      jj_counter++;

      /* Loop over entries in row ic of R */
      for (jj1 = R_ptr[ic]; jj1 < R_ptr[ic+1]; jj1++)
	{
	  i1 = R_ind[jj1];
	  
	  /* Loop over entries in row i1 of A */
	  for (jj2 = A_ptr[i1]; jj2 < A_ptr[i1+1]; jj2++)
	    {
	      i2 = A_ind[jj2];

	      /* 
	       * Check A_marker to see if point i2 has been previously
	       * visited. New entries in R^T * A * P only occur from
	       * unmarked points.
	       */
	      if (A_marker[i2] != ic)
		{
		  /* Mark i2 as visited */
		  A_marker[i2] = ic;

		  /* Loop over entries in row i2 of P */
		  for (jj3 = P_ptr[i2]; jj3 < P_ptr[i2+1]; jj3++)
		    {
		      i3 = P_ind[jj3];

		      /* 
		       * Check P_marker to see that entry (ic,i3) 
		       * of the result matrix has not already been
		       * accounted for.  If not, mark it and 
		       * increment counter.
		       */
		      if (P_marker[i3] < jj_row_beginning)
			{
			  P_marker[i3] = jj_counter;
			  jj_counter++;
			}
		    }
		}
	    }
	}
      /* Set R^T*A*P 's ptr for this row */
      _ptr[ic] = jj_row_beginning;
    }

  _ptr[n_coarse] = jj_counter; /* == nnz in R^T * A * P */

  /* 
   * Now we know how many entries will be in the result matrix, so
   * allocate the arrays 
   */
  nnz = jj_counter;
  _val = bebop_malloc (nnz * sizeof (double));
  _ind = bebop_malloc (nnz * sizeof (int));

  /* 
   * Second pass: fill in _val and _ind arrays 
   */

  jj_counter = start_indexing;
  for (ic = 0; ic < n_coarse; ic++)
    P_marker[ic] = -1;
  for (i = 0; i < n_fine; i++)
    A_marker[i] = -1;

  /* Loop over coarse points */
  for (ic = 0; ic < n_coarse; ic++)
    {
      /* Create diagonal entry of result matrix */
      P_marker[ic] = jj_counter;
      jj_row_beginning = jj_counter;
      _val[jj_counter] = 0.0;
      _ind[jj_counter] = ic;
      jj_counter++;

      /* Loop over entries in row ic of R */
      for (jj1 = R_ptr[ic]; jj1 < R_ptr[ic+1]; jj1++)
	{
	  i1 = R_ind[jj1];
	  r_entry = R_val[jj1];

	  /* Loop over entries in row i1 of A */
	  for (jj2 = A_ptr[i1]; jj2 < A_ptr[i1+1]; jj2++)
	    {
	      i2 = A_ind[jj2];
	      r_a_product = r_entry * A_val[jj2];

	      /* 
	       * Check A_marker to see if point i2 has been previously
	       * visited.  New entries in R^T * A * P only occur from
	       * unmarked points. 
	       */
	      if (A_marker[i2] != ic)
		{
		  /* Mark i2 as visited */
		  A_marker[i2] = ic;

		  /* Loop over entries in row i2 of P */
		  for (jj3 = P_ptr[i2]; jj3 < P_ptr[i2+1]; jj3++)
		    {
		      i3 = P_ind[jj3];
		      r_a_p_product = r_a_product * P_val[jj3];

		      /* 
		       * Check P_marker to see that entry (ic,i3) in
		       * the result matrix has not already been
		       * accounted for.  If it has not, create a new
		       * entry.  If it has, add new contribution. 
		       */
		      if (P_marker[i3] < jj_row_beginning)
			{
			  P_marker[i3] = jj_counter;
			  _val[jj_counter] = r_a_p_product;
			  _ind[jj_counter] = i3;
			  jj_counter++;
			}
		      else
			{
			  _val[P_marker[i3]] += r_a_p_product;
			}
		    }
		}
	      /* 
	       * If i2 is previously visited (A_marker[i2] == ic) it
	       * yields no new entries in R^T*A*P, and we can just add
	       * the new contributions. 
	       */
	      else
		{
		  for (jj3 = P_ptr[i2]; jj3 < P_ptr[i2+1]; jj3++)
		    {
		      i3 = P_ind[jj3];
		      r_a_p_product = r_a_product * P_val[jj3];
		      _val[P_marker[i3]] += r_a_p_product;
		    }
		}
	    }
	}
    }

  /* Free R and marker arrays */
  bebop_free (R_val);
  bebop_free (R_ind);
  bebop_free (R_ptr);
  bebop_free (P_marker);
  bebop_free (A_marker);

  /* Return the new matrix */
  *val = _val;
  *ind = _ind;
  *ptr = _ptr;
  return 0;
}



