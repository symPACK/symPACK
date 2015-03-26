/**
 * @file bcoo_matrix.c
 * @author Mark Hoemmen
 * @since 21 Feb 2005
 * @date Time-stamp: <2008-07-16 11:14:23 mhoemmen>
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
#include <bebop/smc/bcoo_matrix.h>
#include <bebop/smc/bcsr_matrix.h>

#include <bebop/util/complex.h>
#include <bebop/util/random_number.h>
#include <bebop/util/malloc.h>
#include <bebop/util/util.h>
#include <bebop/util/sort_joint_arrays.h>
#include <bebop/util/log.h>

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memcpy */

#include <stddef.h> /* ptrdiff_t */

#ifdef USING_VALGRIND_CLIENT_REQUESTS
#  include <valgrind/memcheck.h>
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

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
 * Macro returning the absolute value of x.  
 *
 * Assumes the following operations make sense:
 * - Negation:  -(x)
 * - less-than comparison with the integer zero:  (x) < 0
 *
 * Unfortunately, the only way in C to avoid evaluating x in this macro more
 * than once is to rely on gcc-specific typing tricks.  There's no way (as far
 * as I know) to get away with it in ANSI C99.  If we were coding in Common 
 * Lisp it would be easy ;) 
 */
#ifdef ABS
#  undef ABS
#endif /* #ifdef ABS */
#define ABS( x )  ((x) < 0 ? -(x) : (x))


/** 
 * Initial length allocated for the II and JJ arrays in bcoo_matrix_t.
 *
 * @see bcoo_matrix_t
 */
const int BCOO_ARRAY_INIT_LENGTH = 16;



static int 
my_compare_lexicographically (const int x1, const int y1, const int x2, const int y2)
{
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



/**
 * For pattern values, workspace should be 2 * sizeof (int) * r * c.
 */
static int
partition_pattern (int II[], int JJ[], 
		   const int r, const int c, 
		   int low, int high, void* workspace)
{
  int left, right, pivot;
  int int_pivot_item1 = 0;
  int int_pivot_item2 = 0;
  int int_swapspace = 0;
  void *pivot_item = NULL; 
  void *swapspace = NULL;   
  const size_t size = 0;
  const int block_size = r * c;

  pivot_item = workspace;
  swapspace  = VOIDAREF( workspace,  block_size );

#ifdef SWAP_INTS
#  undef SWAP_INTS
#endif /* SWAP_INTS */
#define SWAP_INTS( array, i, j )  do { \
  int_swapspace = array[i]; \
  array[i] = array[j]; \
  array[j] = int_swapspace; \
} while(0)

  /* mfh 21 Oct 2005: For randomized quicksort, pick a random integer i 
   * between low and high inclusive, and swap A[low] and A[i]. */
  {
    int i = bebop_random_integer (low, high);

    SWAP_INTS( II, low, i );
    SWAP_INTS( JJ, low, i );
  }         

  int_pivot_item1 = II[low];
  int_pivot_item2 = JJ[low];

  pivot = left = low;
  right = high;

  while (left < right)
    {
      /* Move left while item < pivot */
      while (my_compare_lexicographically (II[left], int_pivot_item1, 
					   JJ[left], int_pivot_item2) <= 0 && 
	     left <= high)
	left++;

      /* Move right while item > pivot */
      while (my_compare_lexicographically (II[right], int_pivot_item1, 
					   JJ[right], int_pivot_item2) > 0 && 
	     right >= low)
	right--;

      if (left < right)
	{
	  SWAP_INTS( II, left, right );
	  SWAP_INTS( JJ, left, right );
	}
    }

  /* right is final position for the pivot */

  II[low] = II[right];
  JJ[low] = JJ[right];

  II[right] = int_pivot_item1;
  JJ[right] = int_pivot_item2;

  return right;
}

/**
 * For complex values, workspace should be 2*max (sizeof (double_Complex), sizeof (int)) * r * c.
 */
static int
partition_complex_double (int II[], int JJ[], double_Complex val[], 
			  const int r, const int c, 
		          int low, int high, void* workspace)
{
  int left, right, pivot;
  int int_pivot_item1 = 0;
  int int_pivot_item2 = 0;
  int int_swapspace = 0;
  void *pivot_item = NULL; 
  void *swapspace = NULL;   
  const size_t size = sizeof (double_Complex);
  const int block_size = r * c;

  pivot_item = workspace;
  swapspace  = VOIDAREF( workspace,  block_size );

#ifdef SWAP_INTS
#  undef SWAP_INTS
#endif /* SWAP_INTS */
#define SWAP_INTS( array, i, j )  do { \
  int_swapspace = array[i]; \
  array[i] = array[j]; \
  array[j] = int_swapspace; \
} while(0)

#ifdef SWAP
#  undef SWAP
#endif /* SWAP */
#define SWAP( array, i, j )  do { \
  double_Complex* __v = (double_Complex*) (array);                     \
  memcpy (swapspace, &__v[block_size * i], block_size * sizeof(double_Complex));            \
  memcpy (&__v[block_size * i], &__v[block_size * j], block_size * sizeof(double_Complex)); \
  memcpy (&__v[block_size * j], swapspace, block_size * sizeof(double_Complex));         \
} while(0)

#ifdef ASSIGN
#  undef ASSIGN
#endif /* ASSIGN */
#define ASSIGN( voidp1, voidp2 )  do { \
  memcpy ((voidp1), (voidp2), block_size * sizeof(double_Complex)); \
} while(0)

  /* mfh 21 Oct 2005: For randomized quicksort, pick a random integer i 
   * between low and high inclusive, and swap A[low] and A[i]. */
  {
    int i = bebop_random_integer (low, high);

    SWAP_INTS( II, low, i );
    SWAP_INTS( JJ, low, i );
    SWAP( val, low, i );
  }         

  int_pivot_item1 = II[low];
  int_pivot_item2 = JJ[low];
  ASSIGN (pivot_item, VOIDAREF( val, block_size*low ));

  pivot = left = low;
  right = high;

  while (left < right)
    {
      /* Move left while item < pivot */
      while (my_compare_lexicographically (II[left], int_pivot_item1, 
					   JJ[left], int_pivot_item2) <= 0 && 
	     left <= high)
	left++;

      /* Move right while item > pivot */
      while (my_compare_lexicographically (II[right], int_pivot_item1, 
					   JJ[right], int_pivot_item2) > 0 && 
	     right >= low)
	right--;

      if (left < right)
	{
	  SWAP_INTS( II, left, right );
	  SWAP_INTS( JJ, left, right );
	  SWAP( val, left, right );
	}
    }

  /* right is final position for the pivot */

  II[low] = II[right];
  JJ[low] = JJ[right];
  ASSIGN( VOIDAREF(val, block_size*low), VOIDAREF(val, block_size*right) );

  II[right] = int_pivot_item1;
  JJ[right] = int_pivot_item2;
  ASSIGN( VOIDAREF(val, block_size*right), pivot_item );

  return right;
}




/**
 * For real values, workspace should be 2*max (sizeof (double), sizeof (int)) * r * c.
 */
static int
partition_real_double (int II[], int JJ[], double val[],
		       const int r, const int c, 
		       int low, int high, void* workspace)
{
  int left, right, pivot;
  int int_pivot_item1 = 0;
  int int_pivot_item2 = 0;
  int int_swapspace = 0;
  void *pivot_item = NULL; 
  void *swapspace = NULL;   
  const size_t size = sizeof (double);
  const int block_size = r * c;

  pivot_item = workspace;
  swapspace  = VOIDAREF( workspace,  block_size );

#ifdef SWAP_INTS
#  undef SWAP_INTS
#endif /* SWAP_INTS */
#define SWAP_INTS( array, i, j )  do { \
  int_swapspace = array[i]; \
  array[i] = array[j]; \
  array[j] = int_swapspace; \
} while(0)

#ifdef SWAP
#  undef SWAP
#endif /* SWAP */
#define SWAP( array, i, j )  do { \
  double* __v = (double*) (array);                     \
  memcpy (swapspace, &__v[block_size * i], block_size * sizeof(double));            \
  memcpy (&__v[block_size * i], &__v[block_size * j], block_size * sizeof(double)); \
  memcpy (&__v[block_size * j], swapspace, block_size * sizeof(double));         \
} while(0)

#ifdef ASSIGN
#  undef ASSIGN
#endif /* ASSIGN */
#define ASSIGN( voidp1, voidp2 )  do { \
  memcpy ((voidp1), (voidp2), block_size * sizeof(double)); \
} while(0)

  /* mfh 21 Oct 2005: For randomized quicksort, pick a random integer i 
   * between low and high inclusive, and swap A[low] and A[i]. */
  {
    int i = bebop_random_integer (low, high);

    SWAP_INTS( II, low, i );
    SWAP_INTS( JJ, low, i );
    SWAP( val, low, i );
  }         

  int_pivot_item1 = II[low];
  int_pivot_item2 = JJ[low];
  ASSIGN (pivot_item, VOIDAREF( val, block_size*low ));

  pivot = left = low;
  right = high;

  while (left < right)
    {
      /* Move left while item < pivot */
      while (my_compare_lexicographically (II[left], int_pivot_item1, 
					   JJ[left], int_pivot_item2) <= 0 && 
	     left <= high)
	left++;

      /* Move right while item > pivot */
      while (my_compare_lexicographically (II[right], int_pivot_item1, 
					   JJ[right], int_pivot_item2) > 0 && 
	     right >= low)
	right--;

      if (left < right)
	{
	  SWAP_INTS( II, left, right );
	  SWAP_INTS( JJ, left, right );
	  SWAP( val, left, right );
	}
    }

  /* right is final position for the pivot */

  II[low] = II[right];
  JJ[low] = JJ[right];
  ASSIGN( VOIDAREF(val, block_size*low), VOIDAREF(val, block_size*right) );

  II[right] = int_pivot_item1;
  JJ[right] = int_pivot_item2;
  ASSIGN( VOIDAREF(val, block_size*right), pivot_item );

  return right;
}


/**
 * If A has value_type REAL, workspace should be 
 * 2*max(sizeof(double), sizeof(int)) * r * c;
 * if A has value_type COMPLEX, workspace should be
 * 2*max(sizeof(double_Complex), sizeof(int)) * r * c.
 */
static int
partition (struct bcoo_matrix_t* A,
	   int low, 
	   int high, 
	   void *workspace)
{
  int left, right, pivot;
  int int_pivot_item1 = 0;
  int int_pivot_item2 = 0;
  int int_swapspace = 0;
  void *pivot_item = NULL; 
  void *swapspace = NULL;   
  size_t size = 0;
  const int r = A->r;
  const int c = A->c;
  const int block_size = r * c;

  bebop_log (4, "=== partition ===\n");
  bebop_log (4, "[low,high] = [%d,%d]\n", low, high);

  if (A->value_type == REAL)
    size = sizeof (double);
  else if (A->value_type == COMPLEX)
    size = sizeof (double_Complex);
  else 
    size = 0;

  pivot_item = workspace;
  swapspace  = VOIDAREF( workspace,  block_size );

#ifdef SWAP_INTS
#  undef SWAP_INTS
#endif /* SWAP_INTS */
#define SWAP_INTS( array, i, j )  do { \
  int_swapspace = array[i]; \
  array[i] = array[j]; \
  array[j] = int_swapspace; \
} while(0)

#ifdef SWAP
#  undef SWAP
#endif /* SWAP */
#define SWAP( array, i, j )  do { \
  if (A->value_type == REAL)  \
    { \
      double* __v = (double*) (array);                     \
      memcpy (swapspace, &__v[block_size * i], block_size * sizeof(double));            \
      memcpy (&__v[block_size * i], &__v[block_size * j], block_size * sizeof(double)); \
      memcpy (&__v[block_size * j], swapspace, block_size * sizeof(double));         \
    } \
  else if (A->value_type == COMPLEX) \
    { \
      double_Complex* __v = (double_Complex*) (array);     \
      memcpy (swapspace, &__v[block_size * i], block_size * sizeof(double_Complex)); \
      memcpy (&__v[block_size * i], &__v[block_size * j], block_size * sizeof(double_Complex));   \
      memcpy (&__v[block_size * j], swapspace, block_size * sizeof(double_Complex)); \
    } \
} while(0)

#ifdef ASSIGN
#  undef ASSIGN
#endif /* ASSIGN */
#define ASSIGN( voidp1, voidp2 )  do { \
  if (A->value_type == REAL)        \
    memcpy ((voidp1), (voidp2), block_size * sizeof(double)); \
  else if (A->value_type == COMPLEX) \
    memcpy ((voidp1), (voidp2), block_size * sizeof(double_Complex)); \
} while(0)

#ifdef USING_VALGRIND_CLIENT_REQUESTS
  if (bebop_debug_level() > 1)
    {
      int retval = 0;

      retval = VALGRIND_CHECK_READABLE( A->II + low, sizeof (int) );
      if (retval != 0)
	{
	  bebop_log (0, "*** partition: low=%d, high=%d:  Valgrind says A-"
		   ">II[%d] is not readable, before doing partitioning (readi"
		   "ng from A->II[low])! ***\n", low, high, low);
	  bebop_exit (EXIT_FAILURE);
	}

      retval = VALGRIND_CHECK_READABLE( A->JJ + low, sizeof (int) );
      if (retval != 0)
	{
	  bebop_log (0, "*** partition: low=%d, high=%d:  Valgrind says A->"
		   "JJ[%d] is not readable, before doing partitioning (reading"
		   " from A->JJ[low])! ***\n", low, high, low);
	  bebop_exit (EXIT_FAILURE);
	}
    }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */


  /* mfh 21 Oct 2005: For randomized quicksort, pick a random integer i 
   * between low and high inclusive, and swap A[low] and A[i]. */
  {
    int i = bebop_random_integer (low, high);

    SWAP_INTS( (A->II), low, i );
    SWAP_INTS( (A->JJ), low, i );
    SWAP( (A->val), low, i );
  }         

  int_pivot_item1 = A->II[low];
  int_pivot_item2 = A->JJ[low];
  ASSIGN (pivot_item, VOIDAREF( A->val, block_size*low ));

  pivot = left = low;
  right = high;

  while (left < right)
    {
#ifdef USING_VALGRIND_CLIENT_REQUESTS
      if (bebop_debug_level() > 1)
	{
	  int retval = 0;

	  retval = VALGRIND_CHECK_READABLE( A->II + left, sizeof (int) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** partition: low=%d, high=%d:  Valgrind says A-"
		       ">II[%d] is not readable, in the middle of partitioning (l"
		       "eft comparison)! ***\n", 
		       low, high, left);
	      bebop_exit (EXIT_FAILURE);
	    }

	  retval = VALGRIND_CHECK_READABLE( A->JJ + left, sizeof (int) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** partition: low=%d, high=%d:  Valgrind says A->"
		       "JJ[%d] is not readable, in the middle of partitioning (lef"
		       "t comparison)! ***\n", 
		       low, high, left);
	      bebop_exit (EXIT_FAILURE);
	    }
	}
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

      /* Move left while item < pivot */
      while (my_compare_lexicographically (A->II[left], int_pivot_item1, 
					   A->JJ[left], int_pivot_item2) <= 0 && 
	     left <= high)
	left++;

#ifdef USING_VALGRIND_CLIENT_REQUESTS
      if (bebop_debug_level() > 1)
	{
	  int retval = 0;

	  retval = VALGRIND_CHECK_READABLE( A->II + right, sizeof (int) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** partition: low=%d, high=%d:  Valgrind says A-"
			">II[%d] is not readable, in the middle of partitioning (r"
			"ight comparison)! ***\n", 
			low, high, right);
	      bebop_exit (EXIT_FAILURE);
	    }

	  retval = VALGRIND_CHECK_READABLE( A->JJ + right, sizeof (int) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** partition: low=%d, high=%d:  Valgrind says A->"
			"JJ[%d] is not readable, in the middle of partitioning (rig"
			"ht comparison)! ***\n", 
			low, high, right);
	      bebop_exit (EXIT_FAILURE);
	    }
	}
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

      /* Move right while item > pivot */
      while (my_compare_lexicographically (A->II[right], int_pivot_item1, 
					   A->JJ[right], int_pivot_item2) > 0 && 
	     right >= low)
	right--;

      if (left < right)
	{
#ifdef USING_VALGRIND_CLIENT_REQUESTS
	  if (bebop_debug_level() > 1)
	    {
	      int retval = 0;

	      retval = VALGRIND_CHECK_WRITABLE( A->II + left, sizeof (int) );
	      if (retval != 0)
		{
		  bebop_log (0, "*** partition: low=%d, high=%d:  Valgrind says A->"
			    "II[%d] is not writable, in the middle of partitioning "
			    "(left position, at left-right swap)! ***\n", 
			    low, high, left);
		  bebop_exit (EXIT_FAILURE);
		}

	      retval = VALGRIND_CHECK_WRITABLE( A->JJ + left, sizeof (int) );
	      if (retval != 0)
		{
		  bebop_log (0, "*** partition: low=%d, high=%d:  Valgrind says A->"
			   "JJ[%d] is not writable, in the middle of partitioning (lef"
			   "t position, at left-right swap)! ***\n", low, high, left);
		  bebop_exit (EXIT_FAILURE);
		}

	      retval = VALGRIND_CHECK_WRITABLE( A->II + right, sizeof (int) );
	      if (retval != 0)
		{
		  bebop_log (0, "*** partition: low=%d, high=%d:  Valgrind says A->"
			   "II[%d] is not writable, in the middle of partitioning "
			    "(right position, at left-right swap)! ***\n", 
			    low, high, right);
		  bebop_exit (EXIT_FAILURE);
		}

	      retval = VALGRIND_CHECK_WRITABLE( A->JJ + right, sizeof (int) );
	      if (retval != 0)
		{
		  bebop_log (0, "*** partition: low=%d, high=%d:  Valgrind says A->"
			    "JJ[%d] is not writable, in the middle of partitioning "
			    "(right position, at left-right swap)! ***\n", 
			    low, high, right);
		  bebop_exit (EXIT_FAILURE);
		}
	    }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

	  SWAP_INTS( (A->II), left, right );
	  SWAP_INTS( (A->JJ), left, right );
	  SWAP( (A->val), left, right );
	}
    }

  /* right is final position for the pivot */

#ifdef USING_VALGRIND_CLIENT_REQUESTS
  if (bebop_debug_level() > 1)
    {
      int retval = 0;

      retval = VALGRIND_CHECK_WRITABLE( A->II + low, sizeof (int) );
      if (retval != 0)
	{
	  bebop_log (0, "*** partition: low=%d, high=%d:  Valgrind says A-"
		    ">II[%d] is not writable, at the end of partitioning (assi"
		    "gnment to A->II[low]) ***\n", low, high, low);
	  bebop_exit (EXIT_FAILURE);
	}

      retval = VALGRIND_CHECK_WRITABLE( A->JJ + low, sizeof (int) );
      if (retval != 0)
	{
	  bebop_log (0, "*** partition: low=%d, high=%d:  Valgrind says A->"
		    "JJ[%d] is not writable, at the end of partitioning (assign"
		    "ment to A->II[low]) ***\n", low, high, low);
	  bebop_exit (EXIT_FAILURE);
	}

      retval = VALGRIND_CHECK_READABLE( A->II + right, sizeof (int) );
      if (retval != 0)
	{
	  bebop_log (0, "*** partition: low=%d, high=%d:  Valgrind says A->"
		   "II[%d] is not readable, at the end of partitioning (readin"
		   "g from A->II[right] in order to assign to A->II[low])! ***\n",
		   low, high, right);
	  bebop_exit (EXIT_FAILURE);
	}

      retval = VALGRIND_CHECK_READABLE( A->JJ + right, sizeof (int) );
      if (retval != 0)
	{
	  bebop_log (0, "*** partition: low=%d, high=%d:  Valgrind says A->"
		   "JJ[%d] is not readable, at the end of partitioning (readin"
		   "g from A->JJ[right] in order to assign to A->JJ[low])! ***\n",
		   low, high, right);
	  bebop_exit (EXIT_FAILURE);
	}
    }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

  A->II[low] = A->II[right];
  A->JJ[low] = A->JJ[right];
  ASSIGN( VOIDAREF(A->val, block_size*low), VOIDAREF(A->val, block_size*right) );


#ifdef USING_VALGRIND_CLIENT_REQUESTS
  if (bebop_debug_level() > 1)
    {
      int retval = 0;

      retval = VALGRIND_CHECK_WRITABLE( A->II + right, sizeof (int) );
      if (retval != 0)
	{
	  bebop_log (0, "*** partition: low=%d, high=%d:  Valgrind says A->"
		   "II[%d] is not readable, at the end of partitioning (assign"
		   "ing int_pivot_item1 to A->II[right])! ***\n",
		   low, high, right);
	  bebop_exit (EXIT_FAILURE);
	}

      retval = VALGRIND_CHECK_WRITABLE( A->JJ + right, sizeof (int) );
      if (retval != 0)
	{
	  bebop_log (0, "*** partition: low=%d, high=%d:  Valgrind says A->"
		   "JJ[%d] is not readable, at the end of partitioning (assign"
		   "ing int_pivot_item2 to A->JJ[right])! ***\n",
		   low, high, right);
	  bebop_exit (EXIT_FAILURE);
	}
    }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */


  A->II[right] = int_pivot_item1;
  A->JJ[right] = int_pivot_item2;
  ASSIGN( VOIDAREF(A->val, block_size*right), pivot_item );

  return right;
}

/** 
 * Sort the entries of A in the range [low, high] (note: inclusive range). 
 */
static void
quicksort (struct bcoo_matrix_t* A,
	   int low, 
	   int high, 
	   void* workspace)
{
  int pivot;

  bebop_log (4, "=== quicksort ===\n");
  bebop_log (4, "[low,high] = [%d,%d]\n", low, high);

  if (high > low)
    {
      /*
       * Verify A->II, A->JJ and A->val (the latter only if A is a matrix 
       * that actually contains values, as opposed to a pattern matrix) 
       * before we start doing memcpy operations with them to sort.
       */
#ifdef USING_VALGRIND_CLIENT_REQUESTS
      if (bebop_debug_level() > 1)
	{
	  int retval = 0;
	  int found_error_p = 0;
	  int k;


	  if (! valid_bcoo_matrix_p (A))
	    {
	      bebop_log (0, "*** quicksort(A,%d,%d,workspace): A is invalid! "
			"***\n", low, high);
	      bebop_exit (EXIT_FAILURE);
	    }

	  for (k = low; k <= high; k++)
	    {
	      retval = VALGRIND_CHECK_READABLE( A->II + k, sizeof(int) );
	      if (retval != 0)
		{
		  found_error_p = 1;
		  retval = 0;
		  bebop_log (0, "*** quicksort: Valgrind says A->II[%d] is"
			    " not readable! nnzb=%d, low=%d, high=%d ***\n", k, 
			    A->nnzb, low, high);
		}

	      retval = VALGRIND_CHECK_READABLE( A->JJ + k, sizeof(int) );
	      if (retval != 0)
		{
		  found_error_p = 1;
		  retval = 0;
		  bebop_log (0, "*** quicksort: Valgrind says A->JJ[%d] is"
			    " not readable! nnzb=%d, low=%d, high=%d ***\n", k, 
			    A->nnzb, low, high);
		}
	    }
	  if (found_error_p)
	    {
	      bebop_exit (EXIT_FAILURE);
	    }

	  retval = VALGRIND_CHECK_READABLE( A->II + low, (high - low + 1) * sizeof(int) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** quicksort: Valgrind says A->II[%d:%d] is"
			" not readable! nnzb=%d ***\n", low, high, A->nnzb);
	      bebop_exit (EXIT_FAILURE);
	    }
	  retval = VALGRIND_CHECK_READABLE( A->JJ + low, (high - low + 1) * sizeof(int) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** quicksort: Valgrind says A->JJ[%d:%d] is"
		       " not readable! nnzb=%d ***\n", low, high, A->nnzb);
	      bebop_exit (EXIT_FAILURE);
	    }

	  retval = VALGRIND_CHECK_WRITABLE( A->II + low, (high - low + 1) * sizeof(int) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** quicksort: Valgrind says A->II[%d:%d] is"
		       " not writable! nnzb=%d ***\n", low, high, A->nnzb);
	      bebop_exit (EXIT_FAILURE);
	    }
	  retval = VALGRIND_CHECK_WRITABLE( A->JJ + low, (high - low + 1) * sizeof(int) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** quicksort: Valgrind says A->JJ[%d:%d] is"
		       " not writable! nnzb=%d ***\n", low, high, A->nnzb);
	      bebop_exit (EXIT_FAILURE);
	    }

	  if (A->value_type == REAL)
	    {
	      retval = VALGRIND_CHECK_READABLE( A->val + low, (high - low + 1) * sizeof(double) );
	      if (retval != 0)
		{
		  bebop_log (0, "*** quicksort: Valgrind says A->val in th"
			   "e interval [%d,%d] is not readable! nnzb=%d, r=%d"
			   ", c=%d ***\n", low, high, A->nnzb, A->r, A->c);
		  bebop_exit (EXIT_FAILURE);
		}

	      retval = VALGRIND_CHECK_WRITABLE( A->val + low, (high - low + 1) * sizeof(double) );
	      if (retval != 0)
		{
		  bebop_log (0, "*** quicksort: Valgrind says A->val in th"
			   "e interval [%d,%d] is not writable! nnzb=%d, r=%d"
			   ", c=%d ***\n", low, high, A->nnzb, A->r, A->c);
		  bebop_exit (EXIT_FAILURE);
		}
	    }
	  else if (A->value_type == COMPLEX)
	    {
	      VALGRIND_CHECK_READABLE( A->val + low, (high - low + 1) * sizeof(double_Complex) );
	      if (retval != 0)
		{
		  bebop_log (0, "*** quicksort: Valgrind says A->val in th"
			   "e interval [%d,%d] is not readable! nnzb=%d, r=%d"
			   ", c=%d ***\n", low, high, A->nnzb, A->r, A->c);
		  bebop_exit (EXIT_FAILURE);
		}

	      VALGRIND_CHECK_WRITABLE( A->val + low, (high - low + 1) * sizeof(double_Complex) );
	      if (retval != 0)
		{
		  bebop_log (0, "*** quicksort: Valgrind says A->val in th"
			    "e interval [%d,%d] is not writable! nnzb=%d, r=%d"
			    ", c=%d ***\n", low, high, A->nnzb, A->r, A->c);
		  bebop_exit (EXIT_FAILURE);
		}
	    }
	}			   
#endif /* USING_VALGRIND_CLIENT_REQUESTS */


      pivot = partition (A, low, high, workspace);
      quicksort (A, low, pivot - 1, workspace);
      quicksort (A, pivot + 1, high, workspace);
    }
}
	
static int
sort_bcoo_matrix_by_rows_then_columns (struct bcoo_matrix_t* A)
{
  void* workspace = NULL;

  bebop_log (2, "=== sort_bcoo_matrix_by_rows_then_columns ===\n");

  if (bebop_debug_level() > 0)
    {
      if (! valid_bcoo_matrix_p (A))
	{
	  bebop_log (0, "*** sort_bcoo_matrix_by_rows_then_columns: A is "
		    "invalid on entry to routine! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
    }

  if (A->value_type == REAL)
    workspace = bebop_calloc (3 * A->r * A->c, MAX(sizeof(double), sizeof(int)));
  else if (A->value_type == COMPLEX)
    workspace = bebop_calloc (3 * A->r * A->c, MAX(sizeof(double_Complex), sizeof(int)));
  else if (A->value_type == PATTERN)
    workspace = bebop_calloc (3 * A->r * A->c, sizeof(int));
  else 
    {
      bebop_log (0, "*** sort_bcoo_matrix_by_rows_then_columns: A has "
		"an invalid value_type! ***\n");
      return -1;
    }
  /*
   * Verify A->II, A->JJ and A->val (the latter only if A is a matrix 
   * that actually contains values, as opposed to a pattern matrix) 
   * before we sort them.
   */
#ifdef USING_VALGRIND_CLIENT_REQUESTS
  if (bebop_debug_level() > 1)
    {
      int retval = 0;

      retval = VALGRIND_CHECK_READABLE( A->II, A->nnzb * sizeof(int) );
      if (retval != 0)
	{
	  bebop_log (0, "*** sort_bcoo_matrix_by_rows_then_columns: Valgrind says A->II[0:A->nnzb-1] is invalid! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
      retval = VALGRIND_CHECK_READABLE( A->JJ, A->nnzb * sizeof(int) );
      if (retval != 0)
	{
	  bebop_log (0, "*** sort_bcoo_matrix_by_rows_then_columns: Valgrind says A->JJ[0:A->nnzb-1] is invalid! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
      if (A->value_type == REAL)
	{
	  retval = VALGRIND_CHECK_READABLE( A->val, A->nnzb * A->r * A->c * sizeof(double) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** sort_bcoo_matrix_by_rows_then_columns: Valgrind says A->val is invalid! ***\n");
	      bebop_exit (EXIT_FAILURE);
	    }
	}
      else if (A->value_type == COMPLEX)
	{
	  retval = VALGRIND_CHECK_READABLE( A->val, A->nnzb * A->r * A->c * sizeof(double_Complex) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** sort_bcoo_matrix_by_rows_then_columns: Valgrind says A->val is invalid! ***\n");
	      bebop_exit (EXIT_FAILURE);
	    }
	}
    }			   
#endif /* USING_VALGRIND_CLIENT_REQUESTS */



  quicksort (A, 0, A->nnzb - 1, workspace);
  if (workspace != NULL)
    {
      bebop_log (3, "\tsort_bcoo_matrix_by_rows_then_columns: freeing workspace\n");
      bebop_free (workspace);
    }
  bebop_log (2, "=== Done with sort_bcoo_matrix_by_rows_then_columns ===\n");

  return 0;
}


int
print_bcoo_matrix_in_matrix_market_format (FILE *out, struct bcoo_matrix_t* A)
{
  const int nnzb = A->nnzb;
  const int r = A->r;
  const int c = A->c;
  const int block_size = r * c;
  const int bm = A->bm;
  const int bn = A->bn;
  int i;
  char symmetry_type_label[20];
  char value_type_label[20];
  
  bebop_log (2, "=== print_bcoo_matrix_in_matrix_market_format ===\n");

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
      bebop_log (0, "*** print_bcoo_matrix_in_matrix_market_format: "
	       "Invalid symmetry type %d of the given bcoo_matrix_t sparse "
	       "matrix! ***\n", A->symmetry_type);
      bebop_log (2, "=== Done with print_bcoo_matrix_in_matrix_market_format ===\n");
      return -1;
    }

  if (A->value_type == REAL)
    strncpy (value_type_label, "real", 19);
  else if (A->value_type == COMPLEX)
    strncpy (value_type_label, "complex", 19);
  else if (A->value_type == PATTERN)
    strncpy (value_type_label, "pattern", 19);
  else
    {
      bebop_log (0, "*** print_bcoo_matrix_in_matrix_market_format: Unsupported value type! ***\n");
      bebop_log (2, "=== Done with print_bcoo_matrix_in_matrix_market_format ===\n");
      return -1;
    }

  fprintf (out, "%%%%MatrixMarket matrix coordinate %s %s\n", value_type_label, symmetry_type_label);
  fprintf (out, "%d %d %d\n", bm * r, bn * c, nnzb * block_size);

  if (A->value_type == REAL)
    {
      double* val = (double*) (A->val);

      /* MatrixMarket format uses 1-based indices, so we have to convert
       * to 1-based when we print out the matrix.  */
      for (i = 0; i < nnzb; i++, val += block_size)
	{
	  int ii, jj;
	  const int rowindexbase = A->II[i] + (1 - A->index_base);
	  const int colindexbase = A->JJ[i] + (1 - A->index_base);

	  if (A->col_oriented_p)
	    {
	      for (ii = 0; ii < A->r; ii++)
		for (jj = 0; jj < A->c; jj++)
		  fprintf (out, "%d %d %.13e\n", rowindexbase + ii, 
			   colindexbase + jj, val[jj + ii*c]);
	    }
	  else 
	    {
	      for (ii = 0; ii < A->r; ii++)
		for (jj = 0; jj < A->c; jj++)
		  fprintf (out, "%d %d %.13e\n", rowindexbase + ii, 
			   colindexbase + jj, val[ii + jj*r]);
	    }
	}
    }
  else if (A->value_type == COMPLEX)
    {
      double_Complex* val = (double_Complex*) (A->val);

      /* MatrixMarket format uses 1-based indices, so we have to convert
       * to 1-based when we print out the matrix.  */
      for (i = 0; i < nnzb; i++, val += block_size)
	{
	  int ii, jj;
	  const int rowindexbase = A->II[i] + (1 - A->index_base);
	  const int colindexbase = A->JJ[i] + (1 - A->index_base);

	  if (A->col_oriented_p)
	    {
	      for (ii = 0; ii < r; ii++)
		for (jj = 0; jj < c; jj++)
		  fprintf (out, "%d %d %.13e %.13e\n", rowindexbase + ii, 
			   colindexbase + jj, 
			   double_Complex_real_part (val[jj + ii*c]),
			   double_Complex_imag_part (val[jj + ii*c]));
	    }
	  else 
	    {
	      for (ii = 0; ii < r; ii++)
		for (jj = 0; jj < c; jj++)
		  fprintf (out, "%d %d %.13e %.13e\n", rowindexbase + ii, 
			   colindexbase + jj, 
			   double_Complex_real_part (val[block_size*i + ii + jj*r]),
			   double_Complex_imag_part (val[block_size*i + ii + jj*r]));
	    }
	}
    }
  else if (A->value_type == PATTERN)
    {
      /* MatrixMarket format uses 1-based indices, so we have to convert
       * to 1-based when we print out the matrix.  */
      for (i = 0; i < nnzb; i++)
	{
	  int ii, jj;
	  const int rowindexbase = A->II[i] + (1 - A->index_base);
	  const int colindexbase = A->JJ[i] + (1 - A->index_base);

	  /* It doesn't matter if the blocks are row-oriented or 
	   * column-oriented, since there are no values. */

	  for (ii = 0; ii < r; ii++)
	    for (jj = 0; jj < c; jj++)
	      fprintf (out, "%d %d\n", rowindexbase + ii, colindexbase + jj);
	}
    }

  bebop_log (2, "=== Done with print_bcoo_matrix_in_matrix_market_format ===\n");
  return 0;
}


struct bcoo_matrix_t*
create_bcoo_matrix (const int bm, 
		    const int bn, 
		    const int r, 
		    const int c, 
		    const int nnzb, 
		    const int nnzb_upper_bound, 
		    int *II, 
		    int *JJ, 
		    void* val, 
		    enum index_base_t index_base, 
		    enum symmetry_type_t symmetry_type, 
		    enum symmetric_storage_location_t symmetric_storage_location, 
		    enum value_type_t value_type,
		    enum ownership_mode_t ownership,
		    void (*deallocator) (void*),
		    enum copy_mode_t copy_mode)
{
  struct bcoo_matrix_t* A = NULL;

  bebop_log (2, "=== create_bcoo_matrix ===\n");
  
  A = bebop_malloc (sizeof (struct bcoo_matrix_t));

  assert (bm > 0);
  assert (bn > 0);
  A->bm = bm;
  A->bn = bn;
  assert (r > 0);
  assert (c > 0);
  A->r = r;
  A->c = c;
  assert (nnzb >= 0);
  assert (nnzb <= nnzb_upper_bound);
  A->nnzb = nnzb;
  A->nnzb_upper_bound = nnzb_upper_bound;
  A->index_base = index_base;
  A->symmetry_type = symmetry_type;
  A->symmetric_storage_location = symmetric_storage_location;
  A->value_type = value_type;
  A->ownership = ownership;
  if (deallocator != NULL)
    A->deallocator = deallocator;
  else 
    A->deallocator = &free;

  if (copy_mode == NO_COPY)
    {
      A->II = II;
      A->JJ = JJ;
      A->val = val;
    }
  else
    {
      const int nnz = nnzb * r * c;
      /* In copy mode, the library always deallocates */
      A->ownership = LIBRARY_DEALLOCATES;
      A->deallocator = &free;

      A->II = bebop_malloc (nnz * sizeof (int));
      memcpy (A->II, II, nnz * sizeof (int));
      A->JJ = bebop_malloc (nnz * sizeof (int));
      memcpy (A->JJ, JJ, nnz * sizeof (int));

      if (value_type == REAL) 
	{
	  A->val = bebop_malloc (nnz * sizeof (double));
	  memcpy (A->val, val, nnz * sizeof (double));
	}
      else if (value_type == COMPLEX)
	{
	  A->val = bebop_malloc (nnz * sizeof (double));
	  memcpy (A->val, val, nnz * sizeof (double));
	}
      else
	A->val = NULL;
    }

  return A;
}

struct bcoo_matrix_t*
create_empty_bcoo_matrix (const int bm, 
			  const int bn, 
			  const int r, 
			  const int c, 
			  enum index_base_t index_base, 
			  enum symmetry_type_t symmetry_type, 
			  enum symmetric_storage_location_t symmetric_storage_location,
			  enum value_type_t value_type)
{
  int *II = NULL;
  int *JJ = NULL;
  void *val = NULL;

  /*
  II = bebop_calloc (BCOO_ARRAY_INIT_LENGTH, sizeof (int));
  JJ = bebop_calloc (BCOO_ARRAY_INIT_LENGTH, sizeof (int));
  */
  II = bebop_malloc (BCOO_ARRAY_INIT_LENGTH * sizeof (int));
  JJ = bebop_malloc (BCOO_ARRAY_INIT_LENGTH * sizeof (int));

#ifdef USING_VALGRIND_CLIENT_REQUESTS
  if (bebop_debug_level() > 1)
    {
      /* Make sure that II and JJ are writable up to the upper bound */
      int k;
      int retval = 0;
      for (k = 0; k < BCOO_ARRAY_INIT_LENGTH; k++)
	{
	  /* Recall that if we say II + k, C automatically uses the type of 
	   * the pointer to figure out how much to add.  We don't have to say 
	   * II + k*sizeof(int); in fact if we do, we're adding sizeof(int) 
	   * times more than we should. */
	  retval = VALGRIND_CHECK_WRITABLE( II + k, sizeof(int) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** create_empty_bcoo_matrix: II[0:%d] sho"
		       "uld be writable, but II[%d] is not! ***\n", 
		       BCOO_ARRAY_INIT_LENGTH - 1, k);
	      bebop_exit (EXIT_FAILURE);
	    }
	  retval = VALGRIND_CHECK_WRITABLE( II + k, sizeof(int) );
	  if (retval != 0)
	    {
	      bebop_log (0, "*** create_empty_bcoo_matrix: JJ[0:%d] sho"
		       "uld be writable, but JJ[%d] is not! ***\n", 
		       BCOO_ARRAY_INIT_LENGTH - 1, k);
	      bebop_exit (EXIT_FAILURE);
	    }
	}
    }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

  if (bebop_debug_level() > 1)
    {
      long diff = (unsigned long) JJ - (unsigned long) II;
      long diff2 = BCOO_ARRAY_INIT_LENGTH * sizeof (int);
      if (diff2 > diff)
	{
	  bebop_log (0, "*** create_empty_bcoo_matrix:  Not enough space "
		    "between II and JJ!  need space for BCOO_ARRAY_LENGTH=%d "
		    "entries, but only have space for %lu ***\n", 
		    BCOO_ARRAY_INIT_LENGTH, diff / sizeof(int));
	  bebop_log (0, "*** diff = %lu ***\n", diff);
	  bebop_log (0, "*** II = %p, JJ = %p ***\n", II, JJ);
	  bebop_exit (EXIT_FAILURE);
	}
      
    }
  if (value_type == REAL)
    val = bebop_calloc (BCOO_ARRAY_INIT_LENGTH, sizeof (double));
  else if (value_type == COMPLEX)
    val = bebop_calloc (BCOO_ARRAY_INIT_LENGTH, sizeof (double_Complex));
  else if (value_type == PATTERN)
    val = NULL;
  else
    {
      bebop_log (0, "*** create_empty_bcoo_matrix: value_type invalid! ***\n");
      return NULL;
    }

  return create_bcoo_matrix (bm, bn, r, c, 0, BCOO_ARRAY_INIT_LENGTH, II, JJ, 
			     val, index_base, symmetry_type, 
			     symmetric_storage_location, value_type,
			     LIBRARY_DEALLOCATES, &free, NO_COPY);
}



void
destroy_bcoo_matrix (struct bcoo_matrix_t* A)
{
  bebop_log (2, "=== destroy_bcoo_matrix ===\n");
 
  if (A != NULL && A->ownership == LIBRARY_DEALLOCATES)
    {
      if (A->deallocator == NULL)
	A->deallocator = &free;

      if (A->II != NULL)
	(A->deallocator) (A->II);
      if (A->JJ != NULL)
	(A->deallocator) (A->JJ);
      if (A->val != NULL && A->value_type != PATTERN) /* just in case */
	(A->deallocator) (A->val);

      bebop_free (A);
    }

  bebop_log (2, "=== Done with destroy_bcoo_matrix\n");
}

void
bcoo_matrix_realloc (struct bcoo_matrix_t* A, const int newmaxlength)
{
  bebop_log (2, "=== bcoo_matrix_realloc ===\n");
  bebop_log (2, "\told maxlength = %d, new maxlength = %d\n",
	    A->nnzb_upper_bound, newmaxlength);

  A->ownership = LIBRARY_DEALLOCATES;
  A->deallocator = &free;

  A->nnzb_upper_bound = newmaxlength;
  A->II = (int*) bebop_realloc (A->II, newmaxlength * sizeof (int));
  A->JJ = (int*) bebop_realloc (A->JJ, newmaxlength * sizeof (int));

  if (A->value_type == REAL)
    A->val = bebop_realloc (A->val, newmaxlength * sizeof (double));
  else if (A->value_type == COMPLEX)
    A->val = bebop_realloc (A->val, newmaxlength * sizeof (double_Complex));
  
  if (A->nnzb > newmaxlength)
    {
      bebop_log (0, "*** bcoo_matrix_realloc: WARNING: "
	       "truncating sparse matrix from nnzb %d to nnzb %d ***\n", 
	       A->nnzb, newmaxlength);
      A->nnzb = newmaxlength;
    }
}

void
bcoo_matrix_resize (struct bcoo_matrix_t* A, const int newlength)
{
  bebop_log (2, "=== bcoo_matrix_resize ===\n");
  bebop_log (2, "\told upper bound:  %d\n", A->nnzb_upper_bound);

  if (A->nnzb == newlength)
    {
      bebop_log (2, "\tNo need to resize\n");
      bebop_log (2, "=== Done with bcoo_matrix_resize ===\n");
      return;
    }

  A->nnzb = newlength;

  if (newlength > A->nnzb_upper_bound)
    {
      const int new_maxlen = 2 * newlength;

      bebop_log (2, "\tResizing to new length %d, new upper bound %d\n", 
		newlength, 2*newlength);
      bcoo_matrix_realloc (A, new_maxlen);
    }
  else
    {
      bebop_log (2, "\tNo need to resize: newlength=%d, upper bound = %d\n", 
		newlength, A->nnzb_upper_bound);
    }

  bebop_log (2, "=== Done with bcoo_matrix_resize ===\n");
}


void
bcoo_matrix_add_entry (struct bcoo_matrix_t* A, 
		       const int bi, const int bj, const void* value)
{
  int oldlen, newlen;

  A->ownership = LIBRARY_DEALLOCATES;
  A->deallocator =  &free;

  bebop_log (2, "=== bcoo_matrix_add_entry ===\n");
  bebop_log (2, "\tcurrent nnzb = %d; bi,bj = %d,%d\n", A->nnzb, bi, bj);

  oldlen = A->nnzb;
  newlen = oldlen + 1;

  if (bebop_debug_level() > 1)
    {
      if (! valid_bcoo_matrix_p (A))
	{
	  bebop_log (0, "*** bcoo_matrix_add_entry: before adding entry at"
		    " index %d, data structure A is invalid! ***\n", oldlen);
	  bebop_exit (EXIT_FAILURE);
	}
    }

  /* This makes more space if needed */
  bcoo_matrix_resize (A, newlen);

  /* At this point, A is temporarily invalid, because it contains an 
   * uninitialized entry.  Valgrind will report that entry as not readable.
   * Therefore, we don't run valid_bcoo_matrix_p(A) until we have assigned 
   * to that entry. */

#ifdef USING_VALGRIND_CLIENT_REQUESTS
  if (bebop_debug_level() > 1)
    {
      int retval = 0;

      retval = VALGRIND_CHECK_WRITABLE( A->II + oldlen, sizeof (int) );
      if (retval != 0)
	{
	  bebop_log (0, "*** bcoo_matrix_add_entry: Valgrind says A->II[oldlen=%d] is not writable! ***\n", oldlen);
	  bebop_exit (EXIT_FAILURE);
	}

      retval = VALGRIND_CHECK_WRITABLE( A->JJ + oldlen, sizeof (int) );
      if (retval != 0)
	{
	  bebop_log (0, "*** bcoo_matrix_add_entry: Valgrind says A->JJ[oldlen=%d] is not writable! ***\n", oldlen);
	  bebop_exit (EXIT_FAILURE);
	}
    }
#endif /* USING_VALGRIND_CLIENT_REQUESTS */


  A->II[oldlen] = bi * A->r;
  A->JJ[oldlen] = bj * A->c;

  if (A->value_type == REAL)
    {
      double* val = (double*) (A->val);
      assert (value != NULL);
      val[oldlen] = *((double*) value);
    }
  else if (A->value_type == COMPLEX)
    {
      double_Complex* val = (double_Complex*) (A->val);
      assert (value != NULL);
      val[oldlen] = *((double_Complex*) value);
    }

  if (bebop_debug_level() > 1)
    {
      if (! valid_bcoo_matrix_p (A))
	{
	  bebop_log (0, "*** bcoo_matrix_add_entry: data structure A has b"
		   "een made invalid by assignment after array expansion! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
    }
}


int
valid_bcoo_matrix_p (struct bcoo_matrix_t* A)
{
  int k;
#ifdef USING_VALGRIND_CLIENT_REQUESTS
  int retval = 0;
  int found_error_p = 0;
#endif /* USING_VALGRIND_CLIENT_REQUESTS */

  if (A->bm < 0)
    {
      bebop_log (1, "*** valid_bcoo_matrix_p: A->bm=%d < 0 ***\n", A->bm);
      return 0;
    }
  else if (A->bn < 0)
    {
      bebop_log (1, "*** valid_bcoo_matrix_p: A->bn=%d < 0 ***\n", A->bn);
      return 0;
    }
  else if (A->r < 1)
    { 
      bebop_log (1, "*** valid_bcoo_matrix_p: A->r=%d < 0 ***\n", A->r);
      return 0;
    }
  else if (A->c < 1)
    {
      bebop_log (1, "*** valid_bcoo_matrix_p: A->c=%d < 0 ***\n", A->c);
      return 0;
    }
  else if (A->bm < A->r)
    {
      bebop_log (1, "*** valid_bcoo_matrix_p: A->bm=%d < A->r=%d ***\n", A->bm, A->r);
      return 0;
    }
  else if (A->bn < A->c)
    {
      bebop_log (1, "*** valid_bcoo_matrix_p: A->bn=%d < A->c=%d ***\n", A->bn, A->c);
      return 0;
    }
  else if (A->nnzb < 0)
    {
      bebop_log (1, "*** valid_bcoo_matrix_p: A->nnzb=%d < 0 ***\n", A->nnzb);
      return 0;
    }
  else if (A->nnzb > A->nnzb_upper_bound)
    {
      bebop_log (1, "*** valid_bcoo_matrix_p: A->nnzb=%d > A->nnzb_uppe"
		"r_bound=%d ***\n", A->nnzb, A->nnzb_upper_bound);
      return 0;
    }
  else if (A->nnzb > 0)
    {
      /* int retval = 0; */
      long diff = 0;
      /* int found_error_p = 0; */

      if (A->II == NULL)
	{
          bebop_log (1, "*** valid_bcoo_matrix_p: A->II == NULL, but A"
		    "->nnzb=%d > 0 ***\n", A->nnzb);
	  return 0;
	}
      else if (A->JJ == NULL)
	{
          bebop_log (1, "*** valid_bcoo_matrix_p: A->JJ == NULL, but A->"
		    "nnzb=%d > 0 ***\n", A->nnzb);
	  return 0;
	}
      else if ((A->value_type == REAL || A->value_type == COMPLEX) && A->val == NULL)
	{
	  bebop_log (1, "*** valid_bcoo_matrix_p: A->value_type is REAL or COMPLEX, but A->val == NULL ***\n");
	  return 0;
	}

      /* Make sure that II and JJ have enough room for all the values they 
       * contain.  In particular, make sure that II and JJ don't overlap.  */
      if ((unsigned long) (A->JJ) > (unsigned long) (A->II))
        diff = (unsigned long) (A->JJ) - (unsigned long) (A->II);
      else 
        diff = (unsigned long) (A->II) - (unsigned long) (A->JJ);

      if (diff < A->nnzb * sizeof (int))
	{
	  if (bebop_debug_level () > 0)
	    {
	      bebop_log (0, "*** valid_bcoo_matrix_p: A->II and A->JJ "
			"are too close together:  diff=%lu but "
			"A->nnzb=%d and sizeof(int)=%lu ***\n", 
			(long unsigned int) diff, 
			A->nnzb, 
			(long unsigned int) (sizeof (int)));
            }	
	  return 0;
	}

#ifdef USING_VALGRIND_CLIENT_REQUESTS
      for (k = 0; k < A->nnzb; k++)
	{
	  retval = VALGRIND_CHECK_READABLE( A->II + k, sizeof(int) );
	  if (retval != 0)
	    {
	      found_error_p = 1;
	      retval = 0;
	      bebop_log (1, "*** valid_bcoo_matrix_p: Valgrind says A->II[%d]"
			" is not readable! nnzb=%d ***\n", k, A->nnzb);
	    }

	  retval = VALGRIND_CHECK_READABLE( A->JJ + k, sizeof(int) );
	  if (retval != 0)
	    {
	      found_error_p = 1;
	      retval = 0;
	      bebop_log (1, "*** valid_bcoo_matrix_p: Valgrind "
			"says A->JJ[%d] is not readable! nnzb=%d ***\n",
			k, A->nnzb);
	    }
	}
      if (found_error_p)
	{
	  return 0;
	}

      /* Check that II, JJ and val are long enough */
      VALGRIND_CHECK_READABLE( A->II, A->nnzb * sizeof (int) );
      VALGRIND_CHECK_READABLE( A->JJ, A->nnzb * sizeof (int) );
      VALGRIND_CHECK_WRITABLE( A->II, A->nnzb * sizeof (int) );
      VALGRIND_CHECK_WRITABLE( A->JJ, A->nnzb * sizeof (int) );
      if (A->value_type == REAL)
	{
          VALGRIND_CHECK_READABLE( A->val, A->nnzb * A->r * A->c * sizeof (double) );
          VALGRIND_CHECK_WRITABLE( A->val, A->nnzb * A->r * A->c * sizeof (double) );
	}
      else if (A->value_type == COMPLEX)
	{
          VALGRIND_CHECK_READABLE( A->val, A->nnzb * A->r * A->c * sizeof (double_Complex) );
          VALGRIND_CHECK_WRITABLE( A->val, A->nnzb * A->r * A->c * sizeof (double_Complex) );
	}
#endif /* USING_VALGRIND_CLIENT_REQUESTS */
    }

  for (k = 0; k < A->nnzb; k++)
    {
      const int i = A->II[k];
      const int j = A->JJ[k];

      if (i < A->index_base || i >= A->bm * A->r + A->index_base)
	{
	  bebop_log (1, "*** valid_bcoo_matrix_p: entry %d = (%d"
		    ",%d) is out of valid range (%d,%d) ***\n", 
		    k, i, j, A->bm * A->r + A->index_base, 
		    A->bn * A->c + A->index_base);
	  return 0;
	}
      else if (j < A->index_base || j >= A->bn * A->c + A->index_base)
	{
	  bebop_log (1, "*** valid_bcoo_matrix_p: entry %d = (%d"
		    ",%d) is out of valid range (%d,%d) ***\n", 
		    k, i, j, A->bm * A->r + A->index_base, 
		    A->bn * A->c + A->index_base);
	  return 0;
	}
    }

  return 1;
}


struct bcsr_matrix_t*
bcoo_matrix_to_random_bcsr_matrix (struct bcoo_matrix_t* A)
{
  int i = 0;
  int j = 0;
  int currow = 0;
  int curentry = 0;
  int *rowptr, *colidx;
  double *values = NULL;
  int num_removed = 0;

  bebop_log (2, "=== bcoo_matrix_to_random_bcsr_matrix ===\n");

  /* Verify correctness of the BCOO matrix structure */
  if (bebop_debug_level() > 0)
    {
      if (! valid_bcoo_matrix_p (A))
	{
	  bebop_log (0, "*** bcoo_matrix_to_random_bcsr_matrix: "
		    "input BCOO matrix is invalid! ***\n");
	  return NULL;
	}
    }

  /* Don't allocate storage for the BCSR matrix until after sort_joint_arrays
   * is done.  The sorting routine needs its own storage, and we want to avoid
   * out-of-memory issues. */

  /* 
   * First sort II and JJ by the II indices, then sort them by the JJ indices.  
   * That will put them in a suitable order so we can convert the matrix 
   * structure into BCSR format. 
   */
  bebop_log (2, "\tSorting II and JJ by (i,j) lexicographically\n");

  if (bebop_debug_level() > 0)
    {
      /* Check to make sure that the II and JJ arrays do not overlap.  
       * We should always have abs(diff) >= sizeof(int) * A->nnzb, 
       * otherwise the arrays overlap.  */
      long diff = (unsigned long) (A->JJ) - (unsigned long) (A->II);
      if (ABS(diff) < sizeof(int) * A->nnzb)
	{
	  /* If sizeof (sizeof(int)) > sizeof (int), then we need to use 
	   * %lu in the format string; otherwise we use %u.  Unfortunately
	   * the C preprocessor doesn't allow us to use "sizeof" operators 
	   * in an #if directive. */
	  bebop_log (0, "*** bcoo_matrix_to_random_bcsr_matrix: A->II and "
		    "A->JJ overlap by %ld elements! ***\n", 
		    (long) (A->nnzb) - ABS(diff/sizeof(int)));
	  bebop_log (0, "*** nnzb = %d, diff/sizeof(int) = %ld ***\n", 
		    A->nnzb, diff/sizeof(int));
	  bebop_exit (EXIT_FAILURE);
	}
    }
  sort_joint_arrays (A->II, A->JJ, A->nnzb, sizeof(int), 
		     compare_lexicographically);

  bebop_log (2, "\n\n");
  WITH_DEBUG2( print_bcoo_matrix_in_matrix_market_format(stderr, A) );
  bebop_log (2, "\n\n");

  /* 
   * Now coalesce duplicate entries in the BCOO matrix. 
   */
  bebop_log (2, "\tCoalescing duplicate entries in BCOO matrix\n");
  bcoo_matrix_coalesce_duplicate_entries (A, &num_removed);


  bebop_log (2, "\tRemoved %d duplicate entrie(s)\n", num_removed);

  bebop_log (2, "\n\n"); 
  WITH_DEBUG2( print_bcoo_matrix_in_matrix_market_format(stderr, A) );
  bebop_log (2, "\n\n"); 

  /* Allocate storage for the BCSR matrix */
  bebop_log (2, "\tAllocating storage for the BCSR matrix\n" );
  rowptr = bebop_malloc ((A->bm + 1) * sizeof (int));
  colidx = bebop_malloc (A->nnzb * sizeof (int));
  values = bebop_malloc (A->nnzb * A->r * A->c * sizeof (double));

  /* Fill up the values array with random numbers (uniform [-1,1]) */
  bebop_log (2, "\tFilling the values array with random numbers\n" );
  for (i = 0; i < A->nnzb * A->r * A->c; i++)
    values[i] = bebop_random_double (-1.0, 1.0);

  /* Copy over the BCOO matrix structure into the BCSR matrix */
  bebop_log (2, "\tCopying over the BCOO matrix structure into the BCSR matrix\n");

  currow = 0;    /* _block_ index of the current row */
  curentry = 0;  /* Current index in the colidx array */
  rowptr[0] = 0;

  /* Add each block entry i in the matrix structure to the BCSR matrix */
  for (i = 0; i < A->nnzb; i++)
    {
      if ((A->II[i] - A->index_base) / (A->r) > currow)
	{
	  /*
	   * We may jump more than one block row at a time, so set the rowptr
	   * entries for the empty block rows in between.
	   */
	  for (j = currow+1; j <= (A->II[i] - A->index_base) / (A->r); j++)
	    rowptr[j] = curentry;

	  currow = (A->II[i] - A->index_base) / (A->r);
	}

      colidx[i] = A->JJ[i] - A->index_base;
      curentry++;
    }

#ifdef USING_VALGRIND_CLIENT_REQUESTS
  /* Check if rowptr is long enough */
  VALGRIND_CHECK_WRITABLE( rowptr, (A->bm + 1) * sizeof (int) );
#endif /* USING_VALGRIND_CLIENT_REQUESTS */


  /* Set the last entries in rowptr appropriately */
  for (j = currow+1; j <= A->bm; j++)
    rowptr[j] = A->nnzb;

  /* Verify correctness of the conversion */
  if (bebop_debug_level() > 0)
    {
      /* The nonzero blocks in A are row-oriented, hence the zero at the end */
      struct bcsr_matrix_t* B = 
	create_bcsr_matrix (A->bm, A->bn, A->r, A->c, A->nnzb, values, colidx, 
			    rowptr, A->symmetry_type, A->symmetric_storage_location,
			    REAL, 0, LIBRARY_DEALLOCATES, &free, NO_COPY);

      if (! valid_bcsr_matrix_p (B))
	{
	  bebop_log (0, "*** bcoo_matrix_to_random_bcsr_matrix: "
		    "conversion from BCOO matrix structure to BCSR matrix fail"
		    "ed! ***\n");
	  dump_bcsr_matrix (stderr, B);
	  destroy_bcsr_matrix (B);
	  return NULL;
	}
      else
	{
	  dump_bcsr_matrix (stderr, B);
	  return B;
	}
    }

  /* The nonzero blocks in A are row-oriented, hence the zero at the end */
  return create_bcsr_matrix (A->bm, A->bn, A->r, A->c, A->nnzb, values, colidx, 
			     rowptr, A->symmetry_type, A->symmetric_storage_location,
			     REAL, 0, LIBRARY_DEALLOCATES, &free, NO_COPY);
}



struct bcsr_matrix_t*
bcoo_to_bcsr (struct bcoo_matrix_t* A)
{
  int i = 0;
  int j = 0;
  int currow = 0;
  int curentry = 0;
  int *rowptr, *colidx;
  void *values = NULL;
  int num_removed = 0;
  int errcode = 0;

  bebop_log (2, "=== bcoo_to_bcsr ===\n");

  /* Verify correctness of the BCOO matrix structure */
  if (bebop_debug_level() > 0)
    {
      if (! valid_bcoo_matrix_p (A))
	{
	  bebop_log (0, "*** bcoo_matrix_to_random_bcsr_matrix: "
		   "input BCOO matrix is invalid! ***\n");
	  return NULL;
	}
    }

  /* Don't allocate storage for the BCSR matrix until after sort_joint_arrays
   * is done.  The sorting routine needs its own storage, and we want to avoid
   * out-of-memory issues. */

  /* 
   * First sort II and JJ by the II indices, then sort them by the JJ indices.  
   * Sort the val array along with them.  That will put the entries in a 
   * suitable order so we can convert the matrix into BCSR format. 
   */
  errcode = sort_bcoo_matrix_by_rows_then_columns (A);

  if (errcode != 0)
    {
      bebop_log (0, "*** bcoo_to_bcsr: Failed to sort BCOO matrix! "
		"error code %d ***\n", errcode);
      return NULL;
    }

  bebop_log (2, "\n\n");
  WITH_DEBUG2( print_bcoo_matrix_in_matrix_market_format(stderr, A) );
  bebop_log (2, "\n\n");

  /* 
   * Now coalesce duplicate entries in the BCOO matrix. 
   */
  bebop_log (2, "\tCoalescing duplicate entries in BCOO matrix\n");
  bcoo_matrix_coalesce_duplicate_entries (A, &num_removed);
  bebop_log (2, "\tRemoved %d duplicate entrie(s)\n", num_removed);

  bebop_log (2, "\n\n"); 
  WITH_DEBUG2( print_bcoo_matrix_in_matrix_market_format(stderr, A) ); 
  bebop_log (2, "\n\n");

  /* Allocate storage for the BCSR matrix */
  bebop_log (2, "\tAllocating storage for the BCSR matrix\n");
  rowptr = bebop_malloc ((A->bm + 1) * sizeof (int));
  colidx = bebop_malloc (A->nnzb * sizeof (int));
  values = bebop_malloc (A->nnzb * A->r * A->c * sizeof (double));

  /* Allocate storage for the BCSR matrix */
  bebop_log (2, "\tAllocating storage for the BCSR matrix\n");
  rowptr = bebop_malloc (A->bm * sizeof (int));
  colidx = bebop_malloc (A->nnzb * sizeof (int));
  if (A->value_type == REAL)
    values = bebop_malloc (A->nnzb * A->r * A->c * sizeof (double));
  else if (A->value_type == COMPLEX)
    values = bebop_malloc (A->nnzb * A->r * A->c * sizeof (double_Complex));
  else if (A->value_type == PATTERN)
    values = NULL;

  /* Copy over the values */
  if (A->value_type == REAL) 
    {
      double *v_out = (double*) values;
      double *v_in = (double*) (A->val);

      for (i = 0; i < A->nnzb * A->r * A->c; i++) 
	v_out[i] = v_in[i];
    }
  else if (A->value_type == COMPLEX)
    {
      double_Complex *v_out = (double_Complex*) values;
      double_Complex *v_in = (double_Complex*) (A->val);

      for (i = 0; i < A->nnzb * A->r * A->c; i++) 
	v_out[i] = v_in[i];
    }

  /* Copy over the BCOO matrix structure into the BCSR matrix */
  bebop_log (2, "\tCopying over the BCOO matrix structure into the BCSR matrix\n");

  currow = 0;    /* _block_ index of the current row */
  curentry = 0;  /* Current index in the colidx array */
  rowptr[0] = 0;

  /* Add each block entry i in the matrix structure to the BCSR matrix */
  for (i = 0; i < A->nnzb; i++)
    {
      if ((A->II[i] - A->index_base) / (A->r) > currow)
	{
	  /*
	   * We may jump more than one block row at a time, so set the rowptr
	   * entries for the empty block rows in between.
	   */
	  for (j = currow+1; j <= (A->II[i] - A->index_base) / (A->r); j++)
	    rowptr[j] = curentry;

	  currow = (A->II[i] - A->index_base) / (A->r);
	}

      colidx[i] = A->JJ[i] - A->index_base;
      curentry++;
    }

  /* Set the last entries in rowptr appropriately */
  for (j = currow+1; j <= A->bm; j++)
    rowptr[j] = A->nnzb;

  /* Verify correctness of the conversion */
  if (bebop_debug_level() > 0)
    {
      /* The nonzero blocks in A are row-oriented, hence the zero at the end */
      struct bcsr_matrix_t* B = 
	create_bcsr_matrix (A->bm, A->bn, A->r, A->c, A->nnzb, values, colidx, 
			    rowptr, A->symmetry_type, A->symmetric_storage_location,
			    A->value_type, A->col_oriented_p, LIBRARY_DEALLOCATES,
			    &free, NO_COPY);

      if (! valid_bcsr_matrix_p (B))
	{
	  bebop_log (0, "*** bcoo_to_bcsr: conversion from BCOO matrix to "
		   "BCSR matrix failed! ***\n");
	  dump_bcsr_matrix (stderr, B);
	  destroy_bcsr_matrix (B);
	  return NULL;
	}
      else
	{
	  dump_bcsr_matrix (stderr, B);
	  return B;
	}
    }

  return create_bcsr_matrix (A->bm, A->bn, A->r, A->c, A->nnzb, values, colidx, 
			     rowptr, A->symmetry_type, A->symmetric_storage_location,
			     A->value_type, A->col_oriented_p, LIBRARY_DEALLOCATES,
			     &free, NO_COPY);
}


int
save_bcoo_matrix_in_matrix_market_format (const char* const filename, 
					  struct bcoo_matrix_t* A)
{
  FILE* out = NULL;
  int errcode = 0;

  bebop_log (2, "=== save_bcoo_matrix_in_matrix_market_format ===\n");
  out = fopen (filename, "w");
  if (out == NULL)
    {
      bebop_log (0, "*** save_bcoo_matrix_in_matrix_market_format: failed "
	       "to open output file \"%s\" ***\n", filename);
      bebop_log (2, "=== Done with save_bcoo_matrix_in_matrix_market_format ===\n");
      return -1;
    }

  errcode = print_bcoo_matrix_in_matrix_market_format (out, A);
  if (0 != fclose (out))
    {
      bebop_log (0, "*** save_bcoo_matrix_in_matrix_market_format: failed "
	       "to close output file \"%s\" ***\n", filename);
      bebop_log (2, "=== Done with save_bcoo_matrix_in_matrix_market_format ===\n");
      return -1;
    }
  bebop_log (2, "=== Done with save_bcoo_matrix_in_matrix_market_format ===\n");
  return errcode;
}

void
bcoo_matrix_coalesce_duplicate_entries (struct bcoo_matrix_t* A, int* num_removed)
{
  int __num_removed = 0;
  int prev = 0;
  int cur = 1;

  bebop_log (2, "=== bcoo_matrix_coalesce_duplicate_entries ===\n");


  if (A == NULL)
    {
      *num_removed = 0;
      bebop_log (2, "=== Done with bcoo_matrix_coalesce_duplicate_entries ===\n");
      return;
    }
  else if (A->nnzb <= 0)
    {
      *num_removed = 0;
      bebop_log (2, "=== Done with bcoo_matrix_coalesce_duplicate_entries ===\n");
      return;
    }

  if (bebop_debug_level() > 1)
    {
      /* Make sure that the sorting process didn't damage A */
      if (! valid_bcoo_matrix_p (A))
	{
	  bebop_log (0, "*** bcoo_matrix_coalesce_duplicate_entries: A "
		   "is invalid on entry to routine! ***\n");
	  bebop_exit (EXIT_FAILURE);
	}

      bebop_log (2, "\tOriginal nnzb:  %d\n", A->nnzb);
      bebop_log (2, "\tSorting matrix first by rows, then by columns\n");
    }
  sort_bcoo_matrix_by_rows_then_columns (A);

  if (bebop_debug_level() > 1)
    {
      /* Make sure that the sorting process didn't damage A */
      if (! valid_bcoo_matrix_p (A))
        {
	  bebop_log (0, "*** bcoo_matrix_coalesce_duplicate_entries: sorti"
		   "ng A by rows then columns seems to have damaged its BCOO "
		   "data structure integrity ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
    }

  for (cur = 1; cur < A->nnzb; cur++)
    {
      if (A->II[cur] != A->II[prev] || A->JJ[cur] != A->JJ[prev])
	{
	  A->II[prev+1] = A->II[cur];
	  A->JJ[prev+1] = A->JJ[cur];

	  if (A->value_type == REAL)
	    {
	      double* val = A->val;
	      val[prev+1] = val[prev+1] + val[cur];
	    }
	  else if (A->value_type == COMPLEX)
	    {
	      double_Complex* val = A->val;
	      val[prev+1] = double_Complex_add (val[prev+1], val[cur]);
	    }

	  prev++;
	}
      else
	__num_removed++;
    }

  bebop_log (2, "\tRemoved %d entries\n", __num_removed);

  /* We've reduced the number of nonzero entries by removing duplicates,
   * so update nnzb. */
  if (__num_removed + prev + 1 != A->nnzb)
    {
      bebop_log (0, "*** bcoo_matrix_coalesce_duplicate_entries: __num_"
	       "removed=%d + prev=%d + 1 != A->nnzb=%d ***\n", 
	       __num_removed, prev, A->nnzb);
      bebop_exit (EXIT_FAILURE);
    }
  A->nnzb = prev + 1;

  bebop_log (2, "\tNew nnzb:  %d\n", A->nnzb);

  *num_removed = __num_removed;
}


