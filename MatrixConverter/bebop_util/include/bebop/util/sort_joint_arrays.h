#ifndef _sort_joint_array_h
#define _sort_joint_array_h
/**
 * @file sort_joint_arrays.h
 * @author Mark Hoemmen
 * @since 22 Feb 2005
 * @date Time-stamp: <2008-07-16 10:19:36 mhoemmen>
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
#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h> /* size_t */
#else
   typedef unsigned int size_t;
#endif


/**
 * Sorts the two arrays a1 and a2 jointly, using the comparison function 
 * compar, which takes two elements of the first array and two elements of 
 * the second array and does the desired comparison on them (e.g. sort in 
 * increasing order by the first array, or sort in increasing order by the 
 * second array).  Does not necessarily use a stable sort algorithm, but 
 * does sort in place using only 3*size + O(1) extra storage.
 *
 * @param a1 [IN/OUT]   First array to sort.  Has "nmemb" members, each of 
 *                      which has size "size".
 * @param a2 [IN/OUT]   Second array to sort.  Has "nmemb" members, each of 
 *                      which has size "size".
 * @param nmemb [IN]    Number of elements in a1.  Also number of elements 
 *                      in a2.
 * @param size [IN]     Size in bytes of each element in the arrays.
 * @param compar [IN]   Comparison function:  compar(&a1[i], &a1[j],
 *                                                   &a2[i], &a2[j])
 */
void
sort_joint_arrays (void* a1, void* a2, size_t nmemb, size_t size, 
		   int (*compar) (const void*, const void*, const void*, const void*));

/**
 * Comparison function for sorting two integer arrays jointly, in increasing 
 * order by the elements in the first array.
 */
int
compare_by_first_array_increasing (const void* a1x, const void* a1y, 
				   const void* a2x, const void* a2y);

/**
 * Comparison function for sorting two integer arrays jointly, in increasing 
 * order by the elements in the second array.
 */
int
compare_by_second_array_increasing (const void* a1x, const void* a1y, 
				    const void* a2x, const void* a2y);

/* Comparison function for sorting two integer arrays jointly, 
 * in lexicographical order */
int 
compare_lexicographically (const void* a1x, const void* a1y, 
			   const void* a2x, const void* a2y);

/* Comparison function for sorting two integer arrays jointly, 
 * in reverse lexicographical order (i.e. reversing the pairs) */
int 
compare_reverse_lexicographically (const void* a1x, const void* a1y, 
				   const void* a2x, const void* a2y);

#endif /* NOT _sort_joint_array_h */
