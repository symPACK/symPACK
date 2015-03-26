/**
 * @file spvec.c
 * @author Mark Hoemmen
 * @since 07 Jun 2006
 * @date Time-stamp: <2008-07-16 11:25:26 mhoemmen>
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
#include <bebop/smc/spvec.h>

#include <bebop/util/random_number.h>
#include <bebop/util/malloc.h>
#include <bebop/util/util.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


const int SPVEC_INIT_LENGTH = 4;


/**
 * Expands the reserve space in the sparse vector, without changing the number
 * of elements (unless newmaxlength < current number of elements, in which case
 * we do change the number of elements to newmaxlength).
 */
void
realloc_spvec (struct spvec_t* v, const int newmaxlength);



/*======================================================================*/
void
init_spvec (struct spvec_t* v, const int initial_length, 
	    const enum value_type_t value_type)
{
  int maxlen;

  /* WITH_DEBUG2(fprintf(stderr,"=== init_spvec ===\n")); */
  
  if (initial_length == 0)
    maxlen = SPVEC_INIT_LENGTH;
  else
    maxlen = 2 * initial_length;

  v->value_type = value_type;
  v->len = initial_length;
  v->maxlen = maxlen;
  if (value_type == REAL)
    v->val = (void*) bebop_malloc (maxlen * sizeof (double));
  else if (value_type == COMPLEX)
    v->val = (void*) bebop_malloc (maxlen * sizeof (double_Complex));
  else
    v->val = NULL;

  v->idx = (int*) bebop_calloc (maxlen, sizeof (int));
}


/*======================================================================*/
void
deinit_spvec (struct spvec_t* v)
{
  /* WITH_DEBUG2(fprintf(stderr,"=== deinit_spvec ===\n")); */

  bebop_free (v->val);
  bebop_free (v->idx);

  v->len = 0;
  v->maxlen = 0;
}


/*======================================================================*/
void
destroy_spvec (struct spvec_t* v)
{
  if (v != NULL)
    deinit_spvec (v);
  bebop_free (v);
}


/*======================================================================*/
int
length_spvec (const struct spvec_t* v)
{
  return v->len;
}


/*======================================================================*/
void
print_spvec (FILE* out, const struct spvec_t* v)
{
  int i;
  const int     length = length_spvec (v);
  const int*    idx    = v->idx;
  const enum value_type_t value_type = v->value_type;

  if (value_type == REAL)
    {
      const double* val = (double*) (v->val);
      for (i = 0; i < length; i++)
	fprintf (out, "%d %E\n", idx[i], val[i]);
    }
  else if (value_type == COMPLEX)
    {
      const double_Complex* val = (double_Complex*) (v->val);
      for (i = 0; i < length; i++)
	fprintf (out, "%d %E %E\n", idx[i], 
		 double_Complex_real_part(val[i]), 
		 double_Complex_imag_part(val[i]));
    }
  else 
    {
      for (i = 0; i < length; i++)
	fprintf (out, "%d\n", idx[i]);
    }
}



/*======================================================================*/
void
print_spvec_with_line_starter (FILE* out, 
			       const struct spvec_t* v, 
			       const char* line_starter)
{
  int i;
  const int     length = length_spvec (v);
  const int*    idx    = v->idx;
  const enum value_type_t value_type = v->value_type;

  if (value_type == REAL)
    {
      const double* val = (double*) (v->val);
      for (i = 0; i < length; i++)
	fprintf (out, "%s %d %E\n", line_starter, idx[i], val[i]);
    }
  else if (value_type == COMPLEX)
    {
      const double_Complex* val = (double_Complex*) (v->val);
      for (i = 0; i < length; i++)
	fprintf (out, "%s %d %E %E\n", line_starter, idx[i], 
		 double_Complex_real_part(val[i]), 
		 double_Complex_imag_part(val[i]));
    }
  else 
    {
      for (i = 0; i < length; i++)
	fprintf (out, "%s %d\n", line_starter, idx[i]);
    }
}


/*======================================================================*/
void
resize_spvec (struct spvec_t* v, const int newlength)
{
  /* WITH_DEBUG2(fprintf(stderr, "=== resize_spvec ===")); */
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
  if (newlength > v->maxlen || newlength < MAX(v->maxlen / 2, SPVEC_INIT_LENGTH)) 
  */
  if (newlength > v->maxlen)
    {
      const int new_maxlen = 2 * newlength;

      realloc_spvec (v, new_maxlen);
    }
  else if (newlength == 0)
    {
      bebop_free (v->val);
      bebop_free (v->idx);
      v->val = NULL;
      v->idx = NULL;
    }

  /* WITH_DEBUG2(fprintf(stderr, "=== Done with resize_spvec ===")); */
}


/*======================================================================*/
void
realloc_spvec (struct spvec_t* v, const int newmaxlength)
{
  /* WITH_DEBUG2(fprintf(stderr, "=== realloc_spvec ===")); */
  /* WITH_DEBUG2(fprintf(stderr,"\t--- maxlen from %d to %d ---\n", v->maxlen, newmaxlength)); */

  v->maxlen = newmaxlength;
  if (v->value_type == REAL)
    v->val = (void*) bebop_realloc (v->val, newmaxlength * sizeof (double));
  else if (v->value_type == COMPLEX)
    v->val = (void*) bebop_realloc (v->val, 
				   newmaxlength * sizeof (double_Complex));
   
  v->idx = (int*) bebop_realloc (v->idx, newmaxlength * sizeof (int));

  if (v->len > newmaxlength)
    {
      fprintf (stderr, "*** WARNING:  realloc_spvec: truncating vector from "
	       "length %d to length %d ***\n", v->len, newmaxlength);
      v->len = newmaxlength;
    }
  /* WITH_DEBUG2(fprintf(stderr, "=== Done with realloc_spvec ===")); */
}

/**
 * Append a single real value and index to the given sparse vector.
 *
 * Precondition:  v->value_type == REAL
 */
static void
append_to_spvec_double (struct spvec_t* v, const double value, const int index)
{
  int oldlen, newlen;
  double* dest = (double*) (v->val);

  oldlen = v->len;
  newlen = oldlen + 1;
  resize_spvec (v, newlen);

  dest[oldlen] = value;
  v->idx[oldlen] = index;
}

/**
 * Append an array of real values and corresponding indices to the given 
 * sparse vector.
 *
 * Precondition:  v->value_type == REAL
 */
static void 
append_to_spvec_double_array (struct spvec_t* v, const double* values,
			      const int* indices, const int num_to_append)
{
  int oldlen, newlen, i;
  double* dest;

  /* Bail out now if we're not appending anything, so that we don't
     dereference v->val at an invalid place. */
  if (num_to_append < 1)
    return;

  oldlen = v->len;
  newlen = oldlen + num_to_append;
  resize_spvec (v, newlen);

  dest = (double*) &(v->val[oldlen]);

  for (i = 0; i < num_to_append; i++)
    {
      dest[i] = values[i];
      v->idx[oldlen + i] = indices[i];
    }
}



/**
 * Append a single complex value and index to the given sparse vector.
 */
void
append_to_spvec_complex (struct spvec_t* v, const double_Complex value, 
			 const int index)
{
  int oldlen, newlen;
  double_Complex* dest = (double_Complex*) (v->val);

  oldlen = v->len;
  newlen = oldlen + 1;
  resize_spvec (v, newlen);

  dest[oldlen] = value;
  v->idx[oldlen] = index;
}

/**
 * Append an array of complex values and corresponding indices to the given 
 * sparse vector.
 *
 * Precondition:  v->value_type == COMPLEX
 */
static void 
append_to_spvec_complex_array (struct spvec_t* v, const double_Complex* values,
			       const int* indices, const int num_to_append)
{
  int oldlen, newlen, i;
  double_Complex* dest;

  /* Bail out now if we're not appending anything, so that we don't
     dereference v->val at an invalid place. */
  if (num_to_append < 1)
    return;

  oldlen = v->len;
  newlen = oldlen + num_to_append;
  resize_spvec (v, newlen);

  dest = (double_Complex*) &(v->val[oldlen]);

  for (i = 0; i < num_to_append; i++)
    {
      dest[i] = values[i];
      v->idx[oldlen + i] = indices[i];
    }
}

void
append_to_spvec_pattern (struct spvec_t* v, const int index)
{
  int oldlen, newlen;

  oldlen = v->len;
  newlen = oldlen + 1;
  resize_spvec (v, newlen);

  v->idx[oldlen] = index;
}

/**
 * Append an array of indices to the given pattern sparse vector.
 *
 * Precondition:  v->value_type == PATTERN
 */
static void 
append_to_spvec_pattern_array (struct spvec_t* v, const int* indices, 
			       const int num_to_append)
{
  int oldlen, newlen, i;

  /* Bail out now if we're not appending anything, so that we don't
     dereference v->val at an invalid place. */
  if (num_to_append < 1)
    return;

  oldlen = v->len;
  newlen = oldlen + num_to_append;
  resize_spvec (v, newlen);

  for (i = 0; i < num_to_append; i++)
    v->idx[oldlen + i] = indices[i];
}

/*======================================================================*/
void
append_to_spvec (struct spvec_t* v, const void* values, 
		 const int* indices, const int num_to_append)
{
  int oldlen, newlen, i;

  /* Bail out now if we're not appending anything, so that we don't
     dereference v->val at an invalid place. */
  if (num_to_append < 1)
    return;

  oldlen = v->len;
  newlen = oldlen + num_to_append;

  /* This makes more space if needed */
  resize_spvec (v, newlen);

  if (v->value_type == REAL)
    {
      double* dest = (double*) &(v->val[oldlen]);
      double* src = (double*) values;

      for (i = 0; i < num_to_append; i++)
	dest[i] = src[i];
    }
  else if (v->value_type == COMPLEX)
    {
      double_Complex* val = (double_Complex*) &(v->val[oldlen]);

      for (i = 0; i < num_to_append; i++)
	dest[i] = src[i];
    }
  for (i = 0; i < num_to_append; i++)
    v->idx[oldlen + i] = indices[i];
}


/*======================================================================*/
struct spvec_t* 
create_spvec (const int initial_length)
{
  struct spvec_t* s;
  /* WITH_DEBUG2(fprintf(stderr,"=== create_spvec ===\n")); */

  s = bebop_malloc (sizeof (struct spvec_t));
  init_spvec (s, initial_length);

  /* WITH_DEBUG2(fprintf(stderr, "=== Done with create_spvec ===")); */
  return s;
}


/*======================================================================*/
struct spvec_t*
clone_spvec (const struct spvec_t* src)
{
  struct spvec_t* dest;
  /* WITH_DEBUG2(fprintf(stderr,"=== clone_spvec ===\n")); */

  dest = bebop_malloc (sizeof (struct spvec_t));

  if (src->len > src->maxlen)
    {
      fprintf (stderr, "*** ERROR: clone_spvec: src->len = %d > "
	       "src->maxlen = %d ***", src->len, src->maxlen);
      die (-1);
    }
 
  dest->len = src->len;
  dest->maxlen = src->maxlen;
  dest->value_type = src->value_type;

  if (src->value_type == REAL)
    { 
      dest->val = (void*) bebop_malloc (src->maxlen * sizeof (double));
      /* Only need to copy first len elements */
      memcpy (dest->val, src->val, src->len * sizeof (double)); 
    }
  else if (src->value_type == COMPLEX)
    {
      dest->val = (void*) bebop_malloc (src->maxlen * sizeof (double_Complex));
      /* Only need to copy first len elements */
      memcpy (dest->val, src->val, src->len * sizeof (double_Complex)); 
    }
  else
    dest->val = NULL;

  dest->idx = (int*) bebop_malloc (src->maxlen * sizeof (int));
  memcpy (dest->idx, src->idx, src->len * sizeof (int)); /* Only need to copy first len elements */

  /* WITH_DEBUG2(fprintf(stderr,"=== Done with clone_spvec ===\n")); */
  /* WITH_DEBUG2(fprintf(stderr,"\n")); */
  return dest;
}


/*======================================================================*/
void
scatter_spvec (void* dest, const struct spvec_t* src, const int n)
{
  const int len = length_spvec (src);
  const int* src_idx = src->idx;
  int i;

  if (src->value_type == REAL)
    {
      const double* src_val = (double*) (src->val);
      double* destination = (double*) dest;
      
      /* DON'T FILL dest WITH ZEROS HERE!  That's the caller's job. */
      
      /* Fill in all the nonzero elements */
      for (i = 0; i < len; i++)
	{
	  const int index = src_idx[i];
	  const double value = src_val[i];

	  destination[index] = destination[index] + value;
	}
    }
  else if (src->value_type == COMPLEX)
    {
      const double_Complex* src_val = (double_Complex*) (src->val);
      double_Complex* destination = (double_Complex*) dest;
      
      /* DON'T FILL dest WITH ZEROS HERE!  That's the caller's job. */
      
      /* Fill in all the nonzero elements */
      for (i = 0; i < len; i++)
	{
	  const int index = src_idx[i];
	  const double_Complex value = src_val[i];

	  destination[index] = double_Complex_add (destination[index], value);
	}
    }
}


/*======================================================================*/
void
set_spvec_index (struct spvec_t* v, const int i, 
		 const int index_value)
{
  if (i < 0 || i >= length_spvec (v))
    {
      fprintf(stderr,"*** set_spvec_index:  i=%d out of range [0,%d) "
	     "***\n", i, length_spvec (v));
      die (-1);
    }
  v->idx[i] = index_value;
}


/*======================================================================*/
int
get_spvec_index (struct spvec_t* v, const int i)
{
  if (i < 0 || i >= length_spvec (v))
    {
      fprintf(stderr,"*** get_spvec_index:  i=%d out of range [0,%d) "
	  "***\n", i, length_spvec (v));
      die(-1);
    }
  return v->idx[i];
}



/*======================================================================*/
struct spvec_t* 
filter_small_elements (const struct spvec_t* x, const double tol)
{
  const int len = length_spvec (x);
  const int *idx = x->idx;
  struct spvec_t* y = create_spvec (0);
  int i;

  if (x->value_type == REAL)
    {
      const double *val = (double*) (x->val);

      for (i = 0; i < len; i++)
	{
	  const double newval = val[i];
	  
	  if (fabs (newval) >= tol)
	    append_to_spvec (y, &newval, idx[i]);
	}
    }
  else if (x->value_type == COMPLEX)
    {
      const double_Complex* val = (double_Complex*) (x->val);
      
      for (i = 0; i < len; i++)
	{
	  const double_Complex newval = val[i];

	  if (double_Complex_cabs (newval) >= tol)
	    append_to_spvec (y, &newval, idx[i]);
	}
    }

  return y;
}


/*======================================================================*/
void
swap_spvec (struct spvec_t *a, struct spvec_t *b)
{
  int len = a->len;
  int maxlen = a->maxlen;
  void *val = a->val;
  int *idx = a->idx;
  enum value_type_t value_type = a->value_type;

  a->len = b->len;
  b->len = len;

  a->value_type = b->value_type;
  b->value_type = value_type;

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
quicksort_partition_spvec (struct spvec_t *v, int low, int high)
{
  if (v->value_type == REAL)
    {
      int    *idx = v->idx;
      double *val = (double*) (v->val);
      int     pivot_idx, left, right;
      double  pivot_val;
      int tmp_idx;
      double tmp_val;

      /* WITH_DEBUG2(fprintf(stderr,"=== quicksort_partition_spvec: [%d,%d] ===\n", low, high)); */
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
  else if (v->value_type == COMPLEX)
    {
      int    *idx = v->idx;
      double_Complex *val = (double_Complex*) (v->val);
      int     pivot_idx, left, right;
      double_Complex pivot_val;
      int tmp_idx;
      double_Complex tmp_val;

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
  else
    {
      int    *idx = v->idx;
      int     pivot_idx, left, right;
      int tmp_idx;
      
      {
	int i = bebop_random_integer (low, high);
   
	tmp_idx = idx[low];
	idx[low] = idx[i];
	idx[i] = tmp_idx;
      }      

      pivot_idx = idx[low];
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
	      idx[left] = idx[right];
	      idx[right] = tmp_idx;
	    }
	}

      /* right is final position for the pivot */
      if (right <= high)
	{
	  idx[low] = idx[right];
	  idx[right] = pivot_idx;
	}
      return right;
    }
}


/**
 * Run quicksort on the given sparse vector to sort the entries by index.
 * Only does entries p:r of the sparse vector.  The initial call to this
 * function should be 
 * quicksort_spvec (v,0,length_spvec(v) - 1).
 */
static void
quicksort_spvec (struct spvec_t *v, int low, int high)
{
  /* WITH_DEBUG2(fprintf(stderr,"=== quicksort_spvec: [%d,%d] ===\n", low, high)); */
  if (low < high)
    {
      int pivot = quicksort_partition_spvec (v, low, high);
      quicksort_spvec (v, low, pivot - 1);
      quicksort_spvec (v, pivot + 1, high);
    }
  /* WITH_DEBUG2(fprintf(stderr, "=== Done with quicksort_spvec ===")); */
}

/**
 * Returns true (nonzero) if the entries of v are sorted by index, false 
 * (nonzero) otherwise.  Duplicate indices are allowed.
 */
static int
spvec_sorted_by_index_p (struct spvec_t *v)
{
  int i;
  int n = length_spvec (v);
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
sort_spvec_by_index (struct spvec_t *v)
{
  quicksort_spvec (v, 0, length_spvec (v) - 1);
  assert (spvec_sorted_by_index_p (v));
}




void
coalesce_spvec_entries_with_common_indices (struct spvec_t *v)
{
  int i, prev_idx;
  int *idx = v->idx;
  int n = length_spvec (v);
  int spot = 0;

  sort_spvec_by_index (v);

  if (n < 2)
    return;  /* We're done in that case */

  if (v->value_type == REAL)
    {
      double *val = (double*) (v->val);
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
    }
  else if (v->value_type == COMPLEX)
    {
      double_Complex *val = (double_Complex*) (v->val);
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
    } 
  else
    {
      prev_idx = idx[0];
      for (i = 1; i < n; i++)
	{
	  if (idx[i] > prev_idx)
	    {
	      spot++;
	      idx[spot] = idx[i];
	    }
	}
    }
  v->len = spot + 1;
}


static void
add_spvec_double (struct spvec_t* z,
		  const double alpha,
		  struct spvec_t* x,
		  const double beta,
		  struct spvec_t* y)
{
  const double* xval = (double*) (x->val);
  const double* yval = (double*) (y->val);
  const int xL = length_spvec (x);
  const int yL = length_spvec (y);
  int xi = 0;
  int yi = 0;

  resize_spvec (z, MAX( xL, yL )); /* here's a place to start */
  if (z->len == 0)
    return;

  sort_spvec_by_index (x);
  sort_spvec_by_index (y);

  while (xi < xL && yi < yL)
    {
      int xind = x->ind[xi];
      int yind = y->ind[yi];
      if (xind < yind)
	{
	  append_to_spvec_double (z, alpha * xval[xi], xind);
	  xi++;
	}
      else if (xind > yind)
	{
	  append_to_spvec_double (z, beta * yval[yi], yind);
	  yi++;
	}
      else 
	{
	  append_to_spvec_double (z, alpha * xval[xi] + beta * yval[yi], xind);
	  xi++;
	  yi++;
	}
    }
  if (xi < xL)
    append_to_spvec_double_array (z, &xval[xi], &(x->ind[xi]), xL - xi);
  else if (yi < yL)
    append_to_spvec_double_array (z, &yval[yi], &(y->ind[yi]), yL - yi);
}

static void
add_spvec_double_Complex (struct spvec_t* z,
			  const double_Complex alpha,
			  struct spvec_t* x,
			  const double_Complex beta,
			  struct spvec_t* y)
{
  const double_Complex* xval = (double_Complex*) (x->val);
  const double_Complex* yval = (double_Complex*) (y->val);
  const int xL = length_spvec (x);
  const int yL = length_spvec (y);
  int xi = 0;
  int yi = 0;

  resize_spvec (z, MAX( xL, yL )); /* here's a place to start */
  if (z->len == 0)
    return;

  sort_spvec_by_index (x);
  sort_spvec_by_index (y);

  while (xi < xL && yi < yL)
    {
      int xind = x->ind[xi];
      int yind = y->ind[yi];
      if (xind < yind)
	{
	  append_to_spvec_complex (z, double_Complex_multiply (alpha, xval[xi]),
				   xind);
	  xi++;
	}
      else if (xind > yind)
	{
	  append_to_spvec_complex (z, double_Complex_multiply (beta, yval[yi]),
				   yind);
	  yi++;
	}
      else 
	{
	  append_to_spvec_complex (z, double_Complex_add (double_Complex_multiply (alpha, xval[xi]), double_Complex_multiply (beta, yval[yi])), xind);
	  xi++;
	  yi++;
	}
    }
  if (xi < xL)
    append_to_spvec_complex_array (z, &xval[xi], &(x->ind[xi]), xL - xi);
  else if (yi < yL)
    append_to_spvec_complex_array (z, &yval[yi], &(y->ind[yi]), yL - yi);
}

static void
add_spvec_pattern (struct spvec_t* z, struct spvec_t* x, struct spvec_t* y)
{
  const int xL = length_spvec (x);
  const int yL = length_spvec (y);
  int xi = 0;
  int yi = 0;

  resize_spvec (z, MAX( xL, yL ));
  if (z->len == 0)
    return;

  sort_spvec_by_index (x);
  sort_spvec_by_index (y);

  while (xi < xL && yi < yL)
    {
      int xind = x->ind[xi];
      int yind = y->ind[yi];
      if (xind < yind)
	{
	  append_to_spvec_pattern (z, xind);
	  xi++;
	}
      else if (xind > yind)
	{
	  append_to_spvec_pattern (z, yind);
	  yi++;
	}
      else 
	{
	  append_to_spvec_pattern (z, xind);
	  xi++;
	  yi++;
	}
    }
  if (xi < xL)
    append_to_spvec_pattern_array (z, &(x->ind[xi]), xL - xi);
  else if (yi < yL)
    append_to_spvec_pattern_array (z, &(y->ind[yi]), yL - yi);
}


void
add_spvec (struct spvec_t* z,
	   const void* a,
	   struct spvec_t* x,
	   const void* b,
	   struct spvec_t* y)
{
  const enum value_type_t value_type = x->value_type;
  assert (value_type == y->value_type);
  z->value_type = value_type;

  if (value_type == REAL)
    add_spvec_double (z, *((double*) a), x, *((double*) b), y);
  else if (value_type == COMPLEX)
    add_spvec_complex (z, *((double_Complex*) a), x, *((double_Complex*) b), y);
  else if (value_type == PATTERN)
    add_spvec_pattern (z, x, y);
}


static int
export_spvec_double (int* ind, double* val, struct spvec_t* x)
{
  int* srcind = x->idx;
  double* srcval = (double*) (x->val);
  const int L = length_spvec (x);
  int i;

  for (i = 0; i < L; i++)
    {
      ind[i] = srcind[i];
      val[i] = srcval[i];
    }
  
  return L;
}


static int
export_spvec_complex (int* ind, double_Complex* val, struct spvec_t* x)
{
  int* srcind = x->idx;
  double_Complex* srcval = (double_Complex*) (x->val);
  const int L = length_spvec (x);
  int i;

  for (i = 0; i < L; i++)
    {
      ind[i] = srcind[i];
      val[i] = srcval[i];
    }

  return L;
}


static int
export_spvec_pattern (int* ind, struct spvec_t* x)
{
  int* srcind = x->idx;
  const int L = length_spvec (x);
  int i;

  for (i = 0; i < L; i++)
    ind[i] = srcind[i];

  return L;
}


int
export_spvec (int* ind, void* val, struct spvec_t* x)
{
  enum value_type_t value_type = x->value_type;

  if (value_type == REAL)
    return export_spvec_double (ind, (double*) val, x);
  else if (value_type == COMPLEX)
    return export_spvec_complex (ind, (double_Complex*) val, x);
  else if (value_type == PATTERN)
    return export_spvec_pattern (ind, x);
  else
    return 0; /* dummy to pacify compiler */
}


static void
append_spvec_to_csr_matrix_double (int** ptr, int** ind, double** val,
				   int* current_nnz,
				   int* nnz_upper_bound,
				   struct spvec_t* row,
				   const int which_row,
				   const int total_num_cols)
{
  const int L = length_spvec (row);
  int new_nnz, new_nnz_upper_bound;

  /* Superfluous */
  ptr[which_row] = *current_nnz;

  new_nnz = L + (*current_nnz);
  /* See if we need to make more space in the CSR matrix */
  if (new_nnz > *nnz_upper_bound)
    {
      /* heuristic to avoid too frequent reallocation:  if there are 
       * k columns left to append after this one, then make space for 
       * L * k additional nonzeros */
      const int num_cols_left = total_num_cols - 1 - which_row;
      new_nnz_upper_bound = new_nnz + L * num_cols_left;

      *ind = bebop_realloc (*ind, new_nnz_upper_bound * sizeof (int));
      *val = bebop_realloc (*val, new_nnz_upper_bound * sizeof (double));
    }

  export_spvec_double (*ind, *val, row);
  ptr[which_row+1] = new_nnz;
  *current_nnz = new_nnz;
  *nnz_upper_bound = new_nnz_upper_bound;
}

static void
append_spvec_to_csr_matrix_complex (int** ptr, int** ind, double_Complex** val,
				    int* current_nnz,
				    int* nnz_upper_bound,
				    struct spvec_t* row,
				    const int which_row,
				    const int total_num_cols)
{
  const int L = length_spvec (row);
  int new_nnz, new_nnz_upper_bound;

  /* Superfluous */
  ptr[which_row] = *current_nnz;

  new_nnz = L + (*current_nnz);
  /* See if we need to make more space in the CSR matrix */
  if (new_nnz > *nnz_upper_bound)
    {
      /* heuristic to avoid too frequent reallocation:  if there are 
       * k columns left to append after this one, then make space for 
       * L * k additional nonzeros */
      const int num_cols_left = total_num_cols - 1 - which_row;
      new_nnz_upper_bound = new_nnz + L * num_cols_left;

      *ind = bebop_realloc (*ind, new_nnz_upper_bound * sizeof (int));
      *val = bebop_realloc (*val, new_nnz_upper_bound * sizeof (double_Complex));
    }

  export_spvec_complex (*ind, *val, row);
  ptr[which_row+1] = new_nnz;
  *current_nnz = new_nnz;
  *nnz_upper_bound = new_nnz_upper_bound;
}

static void
append_spvec_to_csr_matrix_pattern (int** ptr, int** ind,
				    int* current_nnz,
				    int* nnz_upper_bound,
				    struct spvec_t* row,
				    const int which_row,
				    const int total_num_cols)
{
  const int L = length_spvec (row);
  int new_nnz, new_nnz_upper_bound;

  /* Superfluous */
  ptr[which_row] = *current_nnz;

  new_nnz = L + (*current_nnz);
  /* See if we need to make more space in the CSR matrix */
  if (new_nnz > *nnz_upper_bound)
    {
      /* heuristic to avoid too frequent reallocation:  if there are 
       * k columns left to append after this one, then make space for 
       * L * k additional nonzeros */
      const int num_cols_left = total_num_cols - 1 - which_row;
      new_nnz_upper_bound = new_nnz + L * num_cols_left;

      *ind = bebop_realloc (*ind, new_nnz_upper_bound * sizeof (int));
    }

  export_spvec_pattern (*ind, *val, row);
  ptr[which_row+1] = new_nnz;
  *current_nnz = new_nnz;
  *nnz_upper_bound = new_nnz_upper_bound;
}


void
append_spvec_to_csr_matrix (int** ptr, int** ind, void** val,
			    int* current_nnz,
			    int* nnz_upper_bound,
			    struct spvec_t* row,
			    const int which_row,
			    const int total_num_cols)
{
  const enum value_type_t value_type = x->value_type;
  
  assert (which_row >= 0 && which_row < total_num_cols);

  if (value_type == REAL)
    append_spvec_to_csr_matrix_double (ptr, ind, (double**) val, current_nnz, 
				       nnz_upper_bound, row, which_row);
  else if (value_type == COMPLEX)
    append_spvec_to_csr_matrix_complex (ptr, ind, (double_Complex**) val, 
					current_nnz, nnz_upper_bound, row, 
					which_row);
  else if (value_type == PATTERN)
    append_spvec_to_csr_matrix_pattern (ptr, ind, current_nnz, nnz_upper_bound, 
					row, which_row);
}




