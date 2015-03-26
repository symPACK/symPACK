/**
 * @file sort_joint_arrays.c
 * @author Mark Hoemmen
 * @since 22 Feb 2005
 * @date Time-stamp: <2008-07-16 10:11:43 mhoemmen>
 *
 * Function for sorting two arrays "jointly," using a comparison function 
 * that takes both arrays' values into account.
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
#include <bebop/util/sort_joint_arrays.h>
#include <bebop/util/random_number.h>
#include <bebop/util/malloc.h>
#include <bebop/util/util.h>
#include <bebop/util/log.h>

#include <stdio.h>
#include <string.h> 
#include <stdlib.h>

/**
 * Assumes that `array' is a void*, an array of elements, each of
 * which has size `size'.  `idx' is the index into the array.  This
 * macro is needed because `array' is a void pointer, so array[idx] is
 * an invalid expression (since the compiler doesn't know how big the
 * elements of `array' are).  So we do the pointer math ourselves.
 *
 * @note The usual `AREF' macro dereferences the pointer; this macro
 * does NOT.
 *
 * @warn This won't work if sizeof(char) != 1.  The reason we cast to
 * char* is to get a T* such that sizeof(T) == 1.
 */
#ifdef VOIDAREF
#  undef VOIDAREF
#endif /* VOIDAREF */
#define VOIDAREF( array, idx )  ((void*) ((char*) (array) + size*(idx)))

#ifdef SWAP
#  undef SWAP
#endif /* SWAP */
#define SWAP( array, i, j )  do { \
  memcpy (swapspace, VOIDAREF(array, i), size); \
  memcpy (VOIDAREF(array, i), VOIDAREF(array, j), size); \
  memcpy (VOIDAREF(array, j), swapspace, size); \
} while(0)

#ifdef ASSIGN
#  undef ASSIGN
#endif /* ASSIGN */
#define ASSIGN( voidp1, voidp2 )  memcpy(voidp1, voidp2, size)



/* If all compilers would support C99, we could do something more friendly 
 * like this.  Unfortunately we have to deal with lame ancient compilers.
 * uintptr_t is defined in <stdint.h>, but is required (according to the 
 * POSIX standard) only on XSI-conformant systems.
 */
/*
#define VOIDAREF( array, idx )  ((void*) ((uintptr_t) (array) + (uintptr_t) size*(idx)))
*/

/* Here is something we could do if we have to deal with Microsoft compilers.
#ifdef _MSC_VER
  typedef __int64 int64_t
#else
  #include <stdint.h>
#endif 
*/


/*========================================================================*/
/*========================================================================*/

static int
partition (void* a1, 
	   void* a2, 
	   size_t size, 
	   int (*compar) (const void*, const void*, const void*, const void*), 
	   int low, 
	   int high, 
	   void *workspace)
{
  int left, right, pivot;
  void *pivot_item1 = workspace;
  void *pivot_item2 = VOIDAREF( workspace, 1 );
  void *swapspace   = VOIDAREF( workspace, 2 );

  /*
  WITH_DEBUG2(fprintf(stderr, "=== partition ===\n"));
  WITH_DEBUG2(fprintf(stderr, "[low,high] = [%d,%d]\n", low, high));
  */


  /* mfh 6 July 2005: For randomized quicksort, pick a random integer i 
   * between low and high inclusive, and swap A[low] and A[i]. */
  {
    int i = bebop_random_integer (low, high);

    /* SWAP( a1, low, i ); */
    memcpy (swapspace, VOIDAREF(a1, low), size); 
    memcpy (VOIDAREF(a1, low), VOIDAREF(a1, i), size); 
    memcpy (VOIDAREF(a1, i), swapspace, size); 

    SWAP( a2, low, i );
  }

  pivot = left = low; 
  right = high;

  ASSIGN( pivot_item1, VOIDAREF(a1, pivot) );
  ASSIGN( pivot_item2, VOIDAREF(a2, pivot) );


  while (left < right)
    {
      /* Move left while item < pivot */
      while (compar (VOIDAREF(a1, left), pivot_item1, VOIDAREF(a2, left), pivot_item2) <= 0 && left <= high) 
	left++;

      /* Move right while item > pivot */
      while (compar (VOIDAREF(a1, right), pivot_item1, VOIDAREF(a2, right), pivot_item2) > 0 && right >= low) 
	right--;

      if (left < right)
	{
	  SWAP( a1, left, right );
	  SWAP( a2, left, right );
	}
    }

  /* right is final position for the pivot */
  ASSIGN( VOIDAREF(a1, low), VOIDAREF(a1, right) );
  ASSIGN( VOIDAREF(a2, low), VOIDAREF(a2, right) );

  ASSIGN( VOIDAREF(a1, right), pivot_item1 );
  ASSIGN( VOIDAREF(a2, right), pivot_item2 );

  return right;
}


static void
quicksort (void* a1, 
	   void* a2, 
	   size_t size, 
	   int (*compar) (const void*, const void*, const void*, const void*), 
	   int low, 
	   int high, 
	   void* workspace)
{
  int pivot;

  /*
  WITH_DEBUG2(fprintf(stderr, "=== quicksort ===\n"));
  WITH_DEBUG2(fprintf(stderr, "[low,high] = [%d,%d]\n", low, high));
  */

  if (high > low)
    {
      pivot = partition (a1, a2, size, compar, low, high, workspace);
      quicksort (a1, a2, size, compar, low, pivot - 1, workspace);
      quicksort (a1, a2, size, compar, pivot + 1, high, workspace);
    }
}
	   


int
compare_by_first_array_increasing (const void* a1x, const void* a1y, 
				   const void* a2x, const void* a2y)
{
  int x1 = *((int*) a1x);
  int y1 = *((int*) a1y);

  if (x1 < y1)
    return -1;
  else if (x1 == y1)
    return 0;

  return +1;
}

int
compare_by_second_array_increasing (const void* a1x, const void* a1y, 
				    const void* a2x, const void* a2y)
{
  int x2 = *((int*) a2x);
  int y2 = *((int*) a2y);

  if (x2 < y2)
    return -1;
  else if (x2 == y2)
    return 0;

  return +1;
}

int 
compare_lexicographically (const void* a1x, const void* a1y, 
			   const void* a2x, const void* a2y)
{
  int x1 = *((int*) a1x);
  int y1 = *((int*) a1y);
  int x2 = *((int*) a2x);
  int y2 = *((int*) a2y);


  if (x1 < y1)
    return -1;
  else if (x1 > y1)
    return +1;
  else
    {
      if (x2 < y2)
	return -1;
      else if (x2 > y2)
	return +1;
      else 
	return 0;
    }

  return 0; /* to pacify the compiler */
}


int 
compare_reverse_lexicographically (const void* a1x, const void* a1y, 
				   const void* a2x, const void* a2y)
{
  int x1 = *((int*) a1x);
  int y1 = *((int*) a1y);
  int x2 = *((int*) a2x);
  int y2 = *((int*) a2y);


  if (x2 < y2)
    return -1;
  else if (x2 > y2)
    return +1;
  else
    {
      if (x1 < y1)
	return -1;
      else if (x1 > y1)
	return +1;
      else 
	return 0;
    }

  return 0; /* to pacify the compiler */
}

static void
__sort_joint_arrays_quicksort (void* a1, void* a2, size_t nmemb, size_t size, 
          int (*compar) (const void*, const void*, const void*, const void*))
{
  void* workspace;

  workspace = (void*) bebop_calloc (3, size); 
  quicksort (a1, a2, size, compar, 0, nmemb - 1, workspace);
  bebop_free (workspace);
}

#if 0
static void
__sort_joint_arrays_insertion_sort (void* a1, void* a2, size_t nmemb, size_t size, 
          int (*compar) (const void*, const void*, const void*, const void*))
{
  int i, j;
  void* swapspace;

  swapspace = (void*) bebop_calloc (1, size);

  for (i = 0; i < nmemb; i++)
    for (j = i; j > 0; j--)
      {
        if (compar (VOIDAREF(a1, j), VOIDAREF(a1, j-1), 
		    VOIDAREF(a2, j), VOIDAREF(a2, j-1)) < 0)
	  {
	    SWAP( a1, j-1, j );
	    SWAP( a2, j-1, j );
	  }
      }

  bebop_free (swapspace);
}
#endif

void
sort_joint_arrays (void* a1, void* a2, size_t nmemb, size_t size, 
		   int (*compar) (const void*, const void*, const void*, const void*))
{
  /* bebop_log (2, "=== sort_joint_arrays ===\n"); */
  __sort_joint_arrays_quicksort (a1, a2, nmemb, size, compar);
  /*__sort_joint_arrays_insertion_sort (a1, a2, nmemb, size, compar);*/
  /* bebop_log (2, "=== Done with sort_joint_arrays ===\n"); */
}



