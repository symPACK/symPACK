/**
 * @file merge_sort.c
 * @author Mark Hoemmen
 * @since 10 June 2004
 * @date Time-stamp: <2008-07-16 10:10:48 mhoemmen>
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
#include <bebop/util/merge_sort.h>
#include <bebop/util/malloc.h>
#include <bebop/util/util.h>

#include <stdio.h>
#include <string.h> 
#include <stdlib.h>

/**
 * Computes the offset memory location &array[idx], without
 * dereferencing "array". Assumes that `array' is a void*, an array of
 * elements, each of which has size `size'.  `idx' is the index into
 * the array.  This macro is needed because `array' is a void pointer,
 * so array[idx] is an invalid expression (since the compiler doesn't
 * know how big the elements of `array' are).  So we do the pointer
 * math ourselves.
 *
 * @note The usual `AREF' macro dereferences the pointer; this macro
 * does NOT.
 */
#define VOIDAREF( array, idx )  ((void*) ((unsigned long)(array) + size*(idx)))


/**
 * Array for use in merges.
 */
void* v1 = NULL;


/**
 * Array for use in merges.
 */
void* v2 = NULL;

 
/**
 * @function merge
 *
 * For merge_sort():  Merge the subarray v[start..middle] with v[middle..end],
 * placing the result back into v.
 */
void
merge (void* v, 
       size_t nmemb,
       size_t size,
       int (*compar) (const void*, const void*),
       int start, int middle, int end);

/** 
 * @function merge_sort_helper
 *
 * Sort the subarray from start to end.  A recursive-enabled version of 
 * merge_sort(), and is used by merge_sort() internally.
 *
 * @note Code borrowed directly from some online computer science
 * lectures, in particular Lectures 2 (merge sort) and 10 (stable
 * sorting): See 
 * <a href="http://camino.rutgers.edu/ut/utsa/cs3343/lecture10.html">
 * this web site. </a>
 */
void 
merge_sort_helper (void* v, 
		   size_t nmemb,
		   size_t size,
		   int (*compar) (const void*, const void*),
		   int start, int end);


void 
merge_sort (void* v, 
	    size_t nmemb,
	    size_t size,
	    int (*compar) (const void*, const void*))
{
  /* 
   * Allocate space for doing merges 
   *
   * mfh 3 Aug 2004:  We allocate before use and deallocate after use, so 
   * that we don't have to remember to deallocate them at the end of our 
   * program.
   */

  v1 = bebop_malloc (nmemb * size);
  v2 = bebop_malloc (nmemb * size);

  merge_sort_helper (v, nmemb, size, compar, 0, nmemb);

  bebop_free (v1);
  v1 = NULL;
  bebop_free (v2);
  v2 = NULL;
}


void 
merge_sort_helper (void* v, 
		   size_t nmemb,
		   size_t size,
		   int (*compar) (const void*, const void*),
		   int start, int end) 
{
  int middle;  /* the middle of the subarray */

  if (start == end) return;  /* no elements to sort */

  if (start == end - 1) return;	/* one element; already sorted! */

  /* find the middle of the array, splitting it into two subarrays */
  middle = (start + end) / 2;

  /* sort the subarray from start..middle */
  merge_sort_helper (v, nmemb, size, compar, start, middle);

  /* sort the subarray from middle..end */
  merge_sort_helper (v, nmemb, size, compar, middle, end);

  /* merge the two sorted halves */
  merge (v, nmemb, size, compar, start, middle, end);
}


void 
merge (void* v, 
       size_t nmemb,
       size_t size,
       int (*compar) (const void*, const void*),
       int start, int middle, int end) 
{
  int	v1_n, v2_n, v1_index, v2_index, i;

  /* number of elements in first subarray */
  v1_n = middle - start;

  /* number of elements in second subarray */
  v2_n = end - middle;

  /* fill v1 and v2 with the elements of the first and second
   * subarrays, respectively
   */
  for (i=0; i<v1_n; i++) 
    memcpy (VOIDAREF(v1, i), VOIDAREF(v, start+i), size);
  for (i=0; i<v2_n; i++) 
    memcpy (VOIDAREF(v2, i), VOIDAREF(v, middle+i), size);

  /* v1_index and v2_index will index into v1 and v2, respectively... */
  v1_index = 0;
  v2_index = 0;

  /* ... as we pick elements from one or the other to place back
   * into v
   */
  for (i=0; (v1_index < v1_n) && (v2_index < v2_n); i++) 
    {
      /* current v1 element less than current v2 element? */
      /* MFH: Note modification (from < to <=) to ensure sort is stable */

      /* if (v1[v1_index] <= v2[v2_index])  */
      if (compar (VOIDAREF(v1, v1_index), VOIDAREF(v2, v2_index)) <= 0) 
	{
	  /* v[start + i] = v1[v1_index++]; */
	  memcpy (VOIDAREF(v, start+i), VOIDAREF(v1, v1_index), size);
	  v1_index++;
	}
      else
	{
	  /* otherwise, the element from v2 belongs there */
	  /* v[start + i] = v2[v2_index++]; */
	  memcpy (VOIDAREF(v, start+i), VOIDAREF(v2, v2_index), size);
	  v2_index++;
	}
    }

  /* clean up; either v1 or v2 may have stuff left in it */
  for (; v1_index < v1_n; i++) 
    {
      /* v[start + i] = v1[v1_index++]; */
      memcpy (VOIDAREF(v, start+i), VOIDAREF(v1, v1_index), size);
      v1_index++;
    }
  for (; v2_index < v2_n; i++) 
    {
      /* v[start + i] = v2[v2_index++]; */
      memcpy (VOIDAREF(v, start+i), VOIDAREF(v2, v2_index), size);
      v2_index++;
    }
}

