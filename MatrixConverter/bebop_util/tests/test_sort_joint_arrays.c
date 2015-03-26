/**
 * @file test_sort_joint_arrays.c
 * @author Mark Hoemmen
 * @since 7 Jul 2005
 * @date Time-stamp: <2008-07-16 10:22:03 mhoemmen>
 *
 * Driver program for testing functions defined in sort_joint_arrays.c.  
 * Returns EXIT_SUCCESS if all tests pass, else prints an error message 
 * and returns EXIT_FAILURE.
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
#include "test_includes.h"
#include <bebop/util/sort_joint_arrays.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * Assumes that `array' is a void*, an array of elements, each of which has size
 * `size'.  `idx' is the index into the array.  This macro is needed because
 * `array' is a void pointer, so array[idx] is an invalid expression (since the
 * compiler doesn't know how big the elements of `array' are).  So we do the 
 * pointer math ourselves.  
 *
 * @note The usual `AREF' macro dereferences the pointer; this macro does NOT.
 *
 * @warn This won't work if sizeof(char) != 1.  The reason we cast to char* is 
 * to get a T* such that sizeof(T) == 1.
 */
#define VOIDAREF( array, idx )  ((void*) ((char*) (array) + size*(idx)))


static void
test_sort_joint_arrays (void* A1, 
			void* A2, 
			const void* const sorted_A1, 
			const void* const sorted_A2, 
			size_t nmemb, 
			size_t size,
			int (*compar) (const void*, const void*, const void*, const void*))
{
  int i;

  sort_joint_arrays (A1, A2, nmemb, size, compar);

  for (i = 0; i < nmemb - 1; i++)
    {
      if (compar (VOIDAREF(A1,i),
		  VOIDAREF(A1,i+1),
		  VOIDAREF(A2,i),
		  VOIDAREF(A2,i+1)) > 0)
	{
	  bebop_log (0, "test_sort_joint_arrays: test FAILED because "
		     "element %d (out of %d elements) is out of order!",
		     i, nmemb);
	  bebop_exit (EXIT_FAILURE);
	}
    }

  for (i = 0; i < nmemb; i++)
    {
      if (compar (VOIDAREF(A1,i),
		  VOIDAREF(sorted_A1,i),
		  VOIDAREF(A2,i),
		  VOIDAREF(sorted_A2,i)) > 0)
	{
	  bebop_log (0, "test_sort_joint_arrays: test FAILED because "
		     "element %d (out of %d elements) doesn\'t match "
		     "the expected sorted output!", i, nmemb);
	  bebop_exit (EXIT_FAILURE);
	}
    }
}


int
main (int argc, char** argv)
{
  int A1[5] = {5, 3, 4, 1, 2};
  int sorted_A1[5] = {1, 2, 3, 4, 5};
  int A2[5] = {50, 40, 30, 10, 20};
  int sorted_A2[5] = {10, 20, 40, 30, 50};

  int A3[5] = {5, 3, 3, 1, 1};
  int sorted_A3[5] = {1, 1, 3, 3, 5};
  int A4[5] = {50, 30, 40, 20, 10};
  int sorted_A4[5] = {10, 20, 30, 40, 50};
  int info = 0;

  bebop_default_initialize (argc, argv, &info);
  assert (info == 0);
  bebop_log (1, "test_sort_joint_arrays: starting tests\n");

  test_sort_joint_arrays (A1, A2, sorted_A1, sorted_A2, (size_t) 5, sizeof(int), compare_by_first_array_increasing);
  test_sort_joint_arrays (A3, A4, sorted_A3, sorted_A4, (size_t) 5, sizeof(int), compare_lexicographically);

  bebop_log (1, "test_sort_joint_arrays: passed all tests\n");
  bebop_exit (EXIT_SUCCESS);
  return EXIT_SUCCESS;
}
