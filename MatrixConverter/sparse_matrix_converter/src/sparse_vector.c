/**
 * @file sparse_vector.c
 * @author Mark Hoemmen
 * @since 16 Jun 2004
 * @date Time-stamp: <2008-07-16 11:25:18 mhoemmen>
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
#include <bebop/smc/sparse_vector.h>

#include <bebop/util/random_number.h>
#include <bebop/util/malloc.h>
#include <bebop/util/util.h>
#include <bebop/util/log.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


const int SPARSE_VECTOR_INIT_LENGTH = 16;


/**
 * Expands the reserve space in the sparse vector, without changing the number
 * of elements (unless newmaxlength < current number of elements, in which case
 * we do change the number of elements to newmaxlength).
 */
void
realloc_sparse_vector (struct sparse_vector_t* v, const int newmaxlength);



/*======================================================================*/
void
init_sparse_vector (struct sparse_vector_t* v, const int initial_length)
{
  int maxlen;

  /* WITH_DEBUG2(fprintf(stderr,"=== init_sparse_vector ===\n")); */
  
  if (initial_length == 0)
    maxlen = SPARSE_VECTOR_INIT_LENGTH;
  else
    maxlen = 2 * initial_length;

  v->len = initial_length;
  v->maxlen = maxlen;
  v->val = (double*) bebop_calloc (maxlen, sizeof (double));
  v->idx = (int*)    bebop_calloc (maxlen, sizeof (int));
}


/*======================================================================*/
void
deinit_sparse_vector (struct sparse_vector_t* v)
{
  /* WITH_DEBUG2(fprintf(stderr,"=== deinit_sparse_vector ===\n")); */

  bebop_free (v->val);
  bebop_free (v->idx);

  v->len = 0;
  v->maxlen = 0;
}


/*======================================================================*/
void
destroy_sparse_vector (struct sparse_vector_t* v)
{
  if (v != NULL)
	deinit_sparse_vector (v);
  bebop_free (v);
}


/*======================================================================*/
int
length_sparse_vector (const struct sparse_vector_t* v)
{
  return v->len;
}


/*======================================================================*/
void
print_sparse_vector (FILE* out, const struct sparse_vector_t* v)
{
  int i;
  const int     length = length_sparse_vector (v);
  const int*    idx    = v->idx;
  const double* val    = v->val;

  for (i = 0; i < length; i++)
	fprintf (out, "%d %E\n", idx[i], val[i]);
}



/*======================================================================*/
void
print_sparse_vector_with_line_starter (FILE* out, 
									   const struct sparse_vector_t* v, 
									   const char* line_starter)
{
  int i;
  const int     length = length_sparse_vector (v);
  const int*    idx    = v->idx;
  const double* val    = v->val;

  for (i = 0; i < length; i++)
	fprintf (out, "%s %d %E\n", line_starter, idx[i], val[i]);
}


/*======================================================================*/
void
resize_sparse_vector (struct sparse_vector_t* v, const int newlength)
{
  /* WITH_DEBUG2(fprintf(stderr, "=== resize_sparse_vector ===")); */
  /* WITH_DEBUG2(fprintf(stderr,"\tResizing from %d to %d ===\n", v->len, newlength)); */

  if (v->len == newlength)
	{
	  /* WITH_DEBUG2(fprintf(stderr,"\t\t--- No need to resize ---\n")); */
      return;
	}

  v->len = newlength;

  /* 
   * It's worthwhile to reallocate when we either extend the vector (so that 
   * we must reallocate), or we shrink it so the new size max(1/2 * original 
   * max length, minimum max length of a sparse vector) (so it's worthwhile to 
   * reallocate).
   */
  /** \bug (mfh 22 July 2004) This code is iffy!!! */
  /*
  if (newlength > v->maxlen || 
	  newlength < MAX(v->maxlen / 2, SPARSE_VECTOR_INIT_LENGTH)) 
	  */
  if (newlength > v->maxlen)
    {
      const int new_maxlen = 2 * newlength;

	  realloc_sparse_vector (v, new_maxlen);
    }

  /* WITH_DEBUG2(fprintf(stderr, "=== Done with resize_sparse_vector ===")); */
}


/*======================================================================*/
void
realloc_sparse_vector (struct sparse_vector_t* v, const int newmaxlength)
{
  /* WITH_DEBUG2(fprintf(stderr, "=== realloc_sparse_vector ===")); */
  /* WITH_DEBUG2(fprintf(stderr,"\t--- maxlen from %d to %d ---\n", v->maxlen, newmaxlength)); */

  v->maxlen = newmaxlength;
  v->val = (double*) bebop_realloc (v->val, newmaxlength * sizeof (double));
  v->idx = (int*)    bebop_realloc (v->idx, newmaxlength * sizeof (int));

  if (v->len > newmaxlength)
    {
      bebop_log (0, "*** WARNING:  realloc_sparse_vector: truncating vector from "
	       "length %d to length %d ***\n", v->len, newmaxlength);
      v->len = newmaxlength;
    }
  /* WITH_DEBUG2(fprintf(stderr, "=== Done with realloc_sparse_vector ===")); */
}


/*======================================================================*/
void
append_to_sparse_vector (struct sparse_vector_t* v, 
			 const double value, const int index)
{
  int oldlen, newlen;

  /* WITH_DEBUG2(fprintf(stderr, "=== append_to_sparse_vector ===")); */
  oldlen = v->len;
  newlen = oldlen + 1;

  /* WITH_DEBUG2(fprintf(stderr,"\t--- (%E, %d) at %d ---\n", value, index, oldlen)); */

  /* This makes more space if needed */
  resize_sparse_vector (v, newlen);

  v->val[oldlen] = value;
  v->idx[oldlen] = index;

  /* WITH_DEBUG2(fprintf(stderr, "=== Done with append_to_sparse_vector ===")); */
  /* WITH_DEBUG2(fprintf(stderr,"\n")); */
}


/*======================================================================*/
struct sparse_vector_t* 
create_sparse_vector (const int initial_length)
{
  struct sparse_vector_t* s;
  /* WITH_DEBUG2(fprintf(stderr,"=== create_sparse_vector ===\n")); */

  s = bebop_malloc (sizeof (struct sparse_vector_t));
  init_sparse_vector (s, initial_length);

  /* WITH_DEBUG2(fprintf(stderr, "=== Done with create_sparse_vector ===")); */
  return s;
}


/*======================================================================*/
struct sparse_vector_t*
clone_sparse_vector (const struct sparse_vector_t* src)
{
  struct sparse_vector_t* dest;
  /* WITH_DEBUG2(fprintf(stderr,"=== clone_sparse_vector ===\n")); */

  dest = bebop_malloc (sizeof (struct sparse_vector_t));

  if (src->len > src->maxlen)
	{
	  bebop_log (0, "*** ERROR: clone_sparse_vector: src->len = %d > "
		   "src->maxlen = %d ***", src->len, src->maxlen);
	  die (-1);
	}
 
  dest->len = src->len;
  dest->maxlen = src->maxlen;

  dest->val = (double*) bebop_calloc (src->maxlen, sizeof (double));
  memcpy (dest->val, src->val, src->len * sizeof (double)); /* Only need to copy first len elements */

  dest->idx = (int*)    bebop_calloc (src->maxlen, sizeof (int));
  memcpy (dest->idx, src->idx, src->len * sizeof (int)); /* Only need to copy first len elements */

  /* WITH_DEBUG2(fprintf(stderr,"=== Done with clone_sparse_vector ===\n")); */
  /* WITH_DEBUG2(fprintf(stderr,"\n")); */
  return dest;
}


/*======================================================================*/
void
scatter_sparse_vector (double* dest, const struct sparse_vector_t* src, const int n)
{
  const int len = length_sparse_vector (src);
  int i;


  /* WITH_DEBUG2(fprintf(stderr, "=== scatter_sparse_vector ===")); */
  /* print_sparse_vector (stderr, src); */

  /* First fill dest with zeros */
  for (i = 0; i < n; i++)
    dest[i] = 0.0; 

  /* Now go back and fill in all the nonzero elements */
  for (i = 0; i < len; i++)
    {
      const int index = src->idx[i];
      const double value = src->val[i];

      dest[index] = value;
    }

  /* WITH_DEBUG2(fprintf(stderr, "=== Done with scatter_sparse_vector ===")); */
}



/*======================================================================*/
void
gather_sparse_vector (struct sparse_vector_t* sv, const double* dv, 
					  const int n, const double tol)
{
  int i;

  /* WITH_DEBUG2(fprintf(stderr, "=== gather_sparse_vector ===")); */
  /* WITH_DEBUG2(fprintf(stderr,"\tn = %d, tol = %g\n", n, tol)); */
  if (tol < 0.0)
    {
      fprintf(stderr, "*** gather_sparse_vector: Called with negative drop "
	     "tolerance tol = %g ***\n", tol);
      die (-1);
    }
  resize_sparse_vector (sv, 0);

  for (i = 0; i < n; i++)
    {
      const double value = dv[i];

      if (fabs (value) > tol) 
	append_to_sparse_vector (sv, value, i);
    }

  /* WITH_DEBUG2(fprintf(stderr, "=== Done with gather_sparse_vector ===")); */
  /* WITH_DEBUG2(fprintf(stderr,"\n")); */
}


/*======================================================================*/
void
set_sparse_vector_index (struct sparse_vector_t* v, const int i, 
			 const int index_value)
{
  if (i < 0 || i >= length_sparse_vector (v))
    {
      fprintf(stderr,"*** set_sparse_vector_index:  i=%d out of range [0,%d) "
	     "***\n", i, length_sparse_vector (v));
      die (-1);
    }
  v->idx[i] = index_value;
}


/*======================================================================*/
int
get_sparse_vector_index (struct sparse_vector_t* v, const int i)
{
  if (i < 0 || i >= length_sparse_vector (v))
    {
      fprintf(stderr,"*** get_sparse_vector_index:  i=%d out of range [0,%d) "
	  "***\n", i, length_sparse_vector (v));
      die(-1);
    }
  return v->idx[i];
}


/*======================================================================*/
void
set_sparse_vector_value (struct sparse_vector_t* v, const int i, 
						 const double value)
{
  if (i < 0 || i >= length_sparse_vector (v))
    {
      fprintf(stderr,"*** set_sparse_vector_value:  i=%d out of range [0,%d) "
	  "***\n", i, length_sparse_vector (v));
      die(-1);
    }
  v->val[i] = value;
}


/*======================================================================*/
double
get_sparse_vector_value (const struct sparse_vector_t* v, const int i)
{
  if (i < 0 || i >= length_sparse_vector (v))
    {
      fprintf(stderr,"*** get_sparse_vector_value:  i=%d out of range [0,%d) "
	  "***\n", i, length_sparse_vector (v));
      die(-1);
    }
  return v->val[i];
}


/*======================================================================*/
double
ddot_svdv (const double* val, const int* idx, const double* x,
           const int start, const int end)
{
  double dot = 0.0;
  int j;

  for (j = start; j < end; j++)
    dot += val[j] * x[idx[j]];

  return dot;
}



/*======================================================================*/
double
infinity_norm_sparse_vector (const struct sparse_vector_t* v)
{
  const int     len = length_sparse_vector (v);
  const double *val = v->val;
  double max = 0.0;
  int i;

  for (i = 0; i < len; i++)
    {
      double curval = fabs (val[i]);

      if (curval > max)
	max = curval;
    }

  return max;
}


/*======================================================================*/
struct sparse_vector_t* 
filter_small_elements (const struct sparse_vector_t* x, const double tol)
{
  const int     len = length_sparse_vector (x);
  const double *val = x->val;
  const int    *idx = x->idx;
  int i;
  struct sparse_vector_t* y = create_sparse_vector (0);

  for (i = 0; i < len; i++)
    {
      const double newval = val[i];

      if (fabs (newval) >= tol)
	append_to_sparse_vector (y, newval, idx[i]);
    }

  return y;
}


/*======================================================================*/
void
swap_sparse_vector (struct sparse_vector_t *a, struct sparse_vector_t *b)
{
  int len = a->len;
  int maxlen = a->maxlen;
  double *val = a->val;
  int *idx = a->idx;

  a->len = b->len;
  b->len = len;

  a->maxlen = b->maxlen;
  b->maxlen = maxlen;

  a->val = b->val;
  b->val = val;

  a->idx = b->idx;
  b->idx = idx;
}


/**
 * Does the partitioning for quicksorting the entries of the given sparse 
 * vector in increasing order of their indices.
 *
 * @param v [IN/OUT]   Sparse vector to partition by indices
 * @param low [IN]     Low index of the array
 * @param high [IN]    High index of the array 
 *
 * @return Pivot index
 */
static int 
quicksort_partition_sparse_vector (struct sparse_vector_t *v, int low, int high)
{
  int    *idx = v->idx;
  double *val = v->val;
  int     pivot_idx, left, right;
  double  pivot_val;
  int tmp_idx;
  double tmp_val;

  /* WITH_DEBUG2(fprintf(stderr,"=== quicksort_partition_sparse_vector: [%d,%d] ===\n", low, high)); */
  /* mfh 21 Oct 2005: For randomized quicksort, pick a random integer i 
   * between low and high inclusive, and swap v[low] and v[i]. */
  {
    int i = bebop_random_integer (low, high);
   
    tmp_idx = idx[low];
    idx[low] = idx[i];
    idx[i] = tmp_idx;

    tmp_val = val[low];
    val[low] = val[i];
    val[i] = tmp_val;
  }      

  pivot_idx = idx[low];
  pivot_val = val[low];
  left = low;
  right = high;

  while (left < right)
    {
      /* Move left while item < pivot */
      while (idx[left] <= pivot_idx) left++;

      /* Move right while item > pivot */
      while (idx[right] > pivot_idx) right--;

      if (left < right && left >= 0)
	{
          tmp_idx = idx[left];
          tmp_val = val[left];
	  idx[left] = idx[right];
	  val[left] = val[right];
	  idx[right] = tmp_idx;
	  val[right] = tmp_val;
	}
    }

  /* right is final position for the pivot */
  if (right <= high)
    {
      idx[low] = idx[right];
      val[low] = val[right];
      idx[right] = pivot_idx;
      val[right] = pivot_val;
    }
  return right;
}


/**
 * Run quicksort on the given sparse vector to sort the entries by index.
 * Only does entries p:r of the sparse vector.  The initial call to this
 * function should be 
 * quicksort_sparse_vector (v,0,length_sparse_vector(v) - 1).
 */
static void
quicksort_sparse_vector (struct sparse_vector_t *v, int low, int high)
{
  /* WITH_DEBUG2(fprintf(stderr,"=== quicksort_sparse_vector: [%d,%d] ===\n", low, high)); */
  if (low < high)
    {
      int pivot = quicksort_partition_sparse_vector (v, low, high);
      quicksort_sparse_vector (v, low, pivot - 1);
      quicksort_sparse_vector (v, pivot + 1, high);
    }
  /* WITH_DEBUG2(fprintf(stderr, "=== Done with quicksort_sparse_vector ===")); */
}

/**
 * Returns true (nonzero) if the entries of v are sorted by index, false 
 * (nonzero) otherwise.  Duplicate indices are allowed.
 */
static int
sparse_vector_sorted_by_index_p (struct sparse_vector_t *v)
{
  int i;
  int n = length_sparse_vector (v);
  int prev_idx;

  if (n == 0) 
    return 1;

  prev_idx = v->idx[0];
  for (i = 1; i < n; i++)
    {
      int cur_idx = v->idx[i];
      if (cur_idx < prev_idx)
	return 0;
      prev_idx = cur_idx;
    }
  return 1;
}

void
sort_sparse_vector_by_index (struct sparse_vector_t *v)
{
  /* WITH_DEBUG2(fprintf(stderr, "=== sort_sparse_vector_by_index ===")); */
  quicksort_sparse_vector (v, 0, length_sparse_vector (v) - 1);
  assert (sparse_vector_sorted_by_index_p (v));
  /* WITH_DEBUG2(fprintf(stderr, "=== Done with sort_sparse_vector_by_index ===")); */
}



void
coalesce_sparse_vector_entries_with_common_indices (struct sparse_vector_t *v)
{
  int i, prev_idx;
  int *idx = v->idx;
  double *val = v->val;
  int n = length_sparse_vector (v);
  int spot = 0;

  /* WITH_DEBUG2(fprintf(stderr, "=== coalesce_sparse_vector_entries_with_common_indices ===")); */
  sort_sparse_vector_by_index (v);

  if (n < 2)
    return;  /* We're done in that case */

  prev_idx = idx[0];
  for (i = 1; i < n; i++)
    {
      if (idx[i] > prev_idx)
	{
	  spot++;
	  idx[spot] = idx[i];
	  val[spot] = val[i];
	}
      else
	val[spot] += val[i];
    }

  v->len = spot + 1;
  /* WITH_DEBUG2(fprintf(stderr, "=== Done with coalesce_sparse_vector_entries_with_common_indices ===")); */
}




