#ifndef _merge_sort_h
#define _merge_sort_h
/**
 * @file merge_sort.h
 * @author Mark Hoemmen
 * @since 06/10/04 10:26:02 PDT
 * @date Time-stamp: <2008-07-16 10:19:08 mhoemmen>
 * 
 * A merge sort function with the same interface as `qsort()'.
 *
 * @warn NOT THREAD-SAFE!!!
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
 * A merge sort function with the same interface as `qsort()'.  
 * Merge sort is a stable sort with guaranteed O(n lg n) runtime, 
 * unlike Quicksort.
 *
 * @param v [IN/OUT]   Array of elements to sort
 * @param nmemb [IN]   Number of elements in the array
 * @param size  [IN]   Size of each element in the array
 * @param compar [IN]  Sorting function (works like strcmp)
 */
void 
merge_sort (void* v, 
	    size_t nmemb,
	    size_t size,
	    int (*compar) (const void*, const void*));

#endif /* _merge_sort_h */
