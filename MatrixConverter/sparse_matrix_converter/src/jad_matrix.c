/*
 * @file jad_matrix.c
 * @author Mark Hoemmen
 * @since 5 Jul 2005
 * @date Time-stamp: <2008-07-16 11:23:46 mhoemmen>
 *
 * Jagged diagonal (JAD) format sparse matrix: struct and member functions.
 * 
 * @note First author: Ankit Jain (ankit@berkeley.edu), later revised
 * by Mark Hoemmen (mhoemmen@cs.berkeley.edu).
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
#include <bebop/smc/jad_matrix.h>
#include <bebop/smc/csr_matrix.h>

#include <bebop/util/complex.h>
#include <bebop/util/malloc.h>
#include <bebop/util/util.h>
#include <bebop/util/log.h>


typedef struct row_len_str 
{
  int len;
  int orig_pos;
} RowStr;


/** 
 * In the given array of RowStr, swap the elements at positions iOne and iTwo.
 */
static void 
swap (RowStr *rs, int iOne, int iTwo)
{
  RowStr rsTemp;

  rsTemp.len = rs[iOne].len;
  rsTemp.orig_pos = rs[iOne].orig_pos;
  rs[iOne].len = rs[iTwo].len;
  rs[iOne].orig_pos = rs[iTwo].orig_pos;
  rs[iTwo].len = rsTemp.len;
  rs[iTwo].orig_pos = rsTemp.orig_pos;
}



/** 
 * Divide rs[iLeft...iRight] into two partitions so the elements in the first 
 * partition are greater than or equal to the elements in the second partition.
 * Return the index of the element that marks the partition boundary. 
 */
static int 
partition (RowStr *rs, int iLeft, int iRight)
{
   int iFirst = iLeft-1;
   int iLast = iRight;

   while (1)
   {
      while (rs[++iFirst].len - rs[iRight].len > 0)
         ;
      while (rs[iRight].len-rs[--iLast].len > 0)
         if (iLast == iLeft)
            break;
      if (iFirst >= iLast)
         break;
      swap (rs, iFirst, iLast);
   }
   swap(rs, iFirst, iRight);
   return iFirst;
}


/**
 * Run Quicksort on the given interval in the given array of RowStr.
 */
static void 
myQsort (RowStr* rs, int iLeft, int iRight)
{
  if (iRight > iLeft)
    {
      int iMid = partition (rs, iLeft, iRight);
      myQsort (rs, iLeft, iMid - 1);
      myQsort (rs, iMid + 1, iRight);
    }
}

/*
struct csr_matrix_t*
jad_to_csr (struct jad_matrix_t* A)
{
  int i;
  void* val;
  int* ind;
  int* ptr;

  for (i = 0; i < A->n_jagged_diagonals; i++)
    {

    }

  return create_csr_matrix (A->m, A->n, A->nnz, val, ind, ptr, UNSYMMETRIC, 0, A->value_type, LIBRARY_DEALLOCATES, NULL, NO_COPY);  
}
*/

/** 
 * this will return a struct jad_matrix_t with all the structures needed
 * for the Jagged Diagonal Algorithm to perform SpMV
 */
struct jad_matrix_t* 
csr_to_jad (struct csr_matrix_t* A)
{
  int num_rows = 0;
  int num_cols = 0;
  int n_nz = 0;
  int* row_start = NULL;
  int* col_ind = NULL;
  
  RowStr* row_lengths = NULL; 
  int* perm = NULL;  
  int* perm_amounts = NULL; 
  void* jdiag = NULL;
  int* jcol_ind = NULL;  
  int* jd_ptr = NULL;
  void* permedy = NULL; 
  
  int i, j, k, maxval;
  struct jad_matrix_t* jS = NULL;

  bebop_log (2, "=== csr_to_jad ===\n");

  if (A == NULL)
    {
      bebop_log (0, "*** csr_to_jad: A is NULL! ***\n");
      return NULL;
    }
  else if (A->m <= 0)
    {
      bebop_log (0, "*** csr_to_jad: A->m <= 0 ***\n");
      return NULL;
    }
  else if (A->n <= 0)
    {
      bebop_log (0, "*** csr_to_jad: A->n <= 0 ***\n");
      return NULL;
    }
  else if (A->nnz < 0)
    {
      bebop_log (0, "*** csr_to_jad: A->nnz < 0 ***\n");
      return NULL;
    }

  num_rows = A->m;
  num_cols = A->n;
  n_nz = A->nnz;
  row_start = A->rowptr;
  col_ind = A->colidx;

  row_lengths  = (RowStr*) bebop_malloc (sizeof(RowStr)*num_rows);
  perm         = (int*)    bebop_malloc (sizeof(int)*num_rows);
  perm_amounts = (int*)    bebop_malloc (sizeof(int)*num_rows);
  /* mfh 30 Mar 2006: shouldn't we use "sizeof(int)" here? */
  jcol_ind     = (int*)    bebop_malloc (sizeof(double)*n_nz);


  if (A->value_type == REAL)
    {
      jdiag   = bebop_malloc (sizeof (double) * n_nz);
      permedy = bebop_malloc (sizeof(double) * num_cols);
    }
  else if (A->value_type == COMPLEX)
    {
      jdiag = bebop_malloc (sizeof (double_Complex) * n_nz);
      permedy = bebop_malloc (sizeof(double_Complex) * num_cols);
    }
  else if (A->value_type == PATTERN)
    {
      jdiag = NULL;
      permedy = NULL;
    }
  else 
    {
      bebop_log (0, "*** csr_to_jad: A has invalid value type! ***\n");
      return NULL;
    }

  
  /* Set up array of row lengths.  The row index corresponding to each row
   * length is also stored, so that in the next step we can set up a 
   * permutation array. */
  for(i = 0; i < num_rows; i++)
    {
      row_lengths[i].len = row_start[i+1]-row_start[i];
      row_lengths[i].orig_pos = i;
    }

  
  /* Sort the rows of the matrix by nnz per row (using row_lengths), and thus 
   * set up permutation array. */
  myQsort (row_lengths, 0, num_rows-1);
  for(i = 0; i < num_rows; i++)
    {
      perm[i] = row_lengths[i].orig_pos;
      perm_amounts[i] = row_lengths[i].len;
    }

  /* 
   * Set up jdiag and jcol_ind
   */
  k = 0;
  maxval = perm_amounts[0];
  jd_ptr = (int*) bebop_malloc (sizeof(int)*(maxval+1));

  if (A->value_type == REAL)
    {
      double *inval = (double*) (A->values);
      double *outval = (double*) jdiag;

      for(i = 0; i < maxval; i++)
	{   
	  jd_ptr[i] = k;
	  for(j = 0; j < num_rows; j++)
	    {
	      if(i < perm_amounts[j])
		{
		  outval[k] = inval[row_start[perm[j]]+i];
		  jcol_ind[k] = col_ind[row_start[perm[j]]+i];
		  k++;
		}
	    }
	}
    }
  else if (A->value_type == COMPLEX)
    {
      double_Complex *inval = (double_Complex*) (A->values);
      double_Complex *outval = (double_Complex*) jdiag;

      for (i = 0; i < maxval; i++)
	{   
	  jd_ptr[i] = k;
	  for(j = 0; j < num_rows; j++)
	    {
	      if(i < perm_amounts[j])
		{
		  outval[k] = inval[row_start[perm[j]]+i];
		  jcol_ind[k] = col_ind[row_start[perm[j]]+i];
		  k++;
		}
	    }
	}
    }
  else if (A->value_type == PATTERN)
    {
      for (i = 0; i < maxval; i++)
	{   
	  jd_ptr[i] = k;
	  for(j = 0; j < num_rows; j++)
	    {
	      if(i < perm_amounts[j])
		{
		  /* outval[k] = inval[row_start[perm[j]]+i]; */
		  jcol_ind[k] = col_ind[row_start[perm[j]]+i];
		  k++;
		}
	    }
	}
    }
  else 
    {
      bebop_log (0, "*** csr_to_jad: A has invalid value type! ***\n");
      return NULL;
    }

  jd_ptr[maxval] = n_nz+1;

  jS = bebop_calloc (sizeof(struct jad_matrix_t), 1);
  
  jS->n_jagged_diagonals = maxval;
  jS->jad_value = jdiag;
  jS->jad_col_index = jcol_ind;
  jS->jad_diag_start = jd_ptr;
  jS->jad_prm_nto = perm;
  jS->Py = permedy;

  jS->value_type = A->value_type;
  jS->symmetry_type = A->symmetry_type;
  jS->symmetric_storage_location = A->symmetric_storage_location;
  jS->ownership = LIBRARY_DEALLOCATES;
  jS->deallocator = &free;
  
  bebop_free (row_lengths);
  bebop_free (perm_amounts);
  
  return jS;

  /*remember to free the following at the end of jagged diagonal
   * perm(:)
   * jdiag
   * jcol_ind
   */
}



void
dealloc_jad_matrix (struct jad_matrix_t* jS)
{
  if (jS == NULL)
    return;

  if (jS->ownership == USER_DEALLOCATES)
    return;

  if (jS->deallocator == NULL)
    {
      if (jS->jad_value != NULL && jS->value_type != PATTERN)
	bebop_free (jS->jad_value);

      bebop_free (jS->jad_col_index);
      bebop_free (jS->jad_diag_start);
      bebop_free (jS->jad_prm_nto);
      bebop_free (jS->Py);
    }
  else
    {
      if (jS->jad_value != NULL && jS->value_type != PATTERN)
	(jS->deallocator) (jS->jad_value);

      (jS->deallocator) (jS->jad_col_index);
      (jS->deallocator) (jS->jad_diag_start);
      (jS->deallocator) (jS->jad_prm_nto);
      (jS->deallocator) (jS->Py);
    }
}

void 
destroy_jad_matrix (struct jad_matrix_t* jS)
{
  if (jS != NULL)
    {
      dealloc_jad_matrix (jS);
      bebop_free (jS);
    }
}


/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/


  /*eof*/
  
  


