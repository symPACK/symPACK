#ifndef _util_h
#define _util_h
/****************************************************************************
 * @file util.h
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-16 10:20:28 mhoemmen>
 *
 * Declaration of some utility functions.  (Used to be called smvm_util.h.)
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
 ****************************************************************************/
#include <bebop/util/config.h>
#include <bebop/util/init.h>

#include <stdlib.h>  /* exit */
#include <stdio.h>   /* fprintf */

#define WITH_DEBUG( OP )  do { if (bebop_debug_level() > 0) { OP; } } while(0)

#define WITH_DEBUG2( OP )  do { if (bebop_debug_level() > 1) { OP; } } while(0)

#define WITH_DEBUG3( OP )  do { if (bebop_debug_level() > 2) { OP; } } while(0)

#define WITH_DEBUG4( OP )  do { if (bebop_debug_level() > 3) { OP; } } while(0)

#ifdef WARN
#  define WITH_WARN( OP )   do { OP; } while(0)
#else
#  define WITH_WARN( OP )  
#endif /* WARN */

#ifndef MIN
/**
 * The standard macro for finding the minimum of two things. 
 *
 * \bug GCC supports the "typeof" operator that lets us define
 * MIN in a typesafe way and only evaluate the inputs once.  Figure
 * out how to play with the #defines so that "typeof" is used when
 * we are using GCC.
 */
#  define MIN(a,b)  (((a)<(b))? (a) : (b))
#endif /* NOT MIN */

#ifndef MAX
/**
 * The standard macro for finding the maximum of two things. 
 *
 * \bug GCC supports the "typeof" operator that lets us define
 * MIN in a typesafe way and only evaluate the inputs once.  Figure
 * out how to play with the #defines so that "typeof" is used when
 * we are using GCC.
 */
#  define MAX(a,b)  (((a)>(b))? (a) : (b))
#endif /* NOT MAX */


/**
 * A comparison function to pass into qsort; compares two doubles.
 *
 * @bug May not handle the NaN case; I need to check that (mfh).
 *
 * @return -1 if *d1 < *d2, 1 if >, 0 if ==.
 */
int 
bebop_compare_double (const void *d1, const void *d2);


/**
 * Finds the median, min and max of the given array.  
 * Sorts the array as a side effect.  Results undefined if n < 1.
 *
 * @param a [IN,OUT]     Array of doubles with n elements.
 * @param n [IN]         Number of elements in a[].
 * @param p_median [OUT] Median of the given array is written here.
 * @param p_min [OUT]    Minimum of the given array is written here.
 * @param p_max [OUT]    Maximum of the given array is written here.
 */
void 
bebop_get_median_min_max (double a[], int n, double *p_median, double *p_min, double *p_max);

/**
 * Returns the median of a given array of doubles.  If the 
 * array is empty, returns 0.0 just to do something.
 *
 * @param values     array of doubles
 * @param n          length of values[] array
 */
double 
bebop_median (double* values, int n);

/*
 * Initialize x[] (of length n) with random doubles in [-1,1].
 * No memory allocation here!
 */
void 
bebop_init_vector_rand (int n, double* x);

/*
 * Initialize x[] (of length n), setting each value to val.
 * No memory allocation here!
 */
void 
bebop_init_vector_val (int n, double* x, double val);


/*
 * Prints msg to stderr (appending an endline), along with the current file name
 * and line number, and calls `exit' with the given exitval.  If msg is NULL,
 * prints nothing.  When calling this macro, append it with a semicolon.  I know
 * macros are icky, but the only way to get the file name and line number in it
 * is to use a macro.
 */
#define die_with_message( msg, exitval )   do { \
  fprintf (stderr, "%s at %s line %d\n", msg, __FILE__, __LINE__);  \
  bebop_exit (exitval); } while(0)


/**
 * Prints the current file name and line number to stderr, and calls `exit' 
 * with the given exitval.  When calling this macro, append it with a 
 * semicolon.
 */
#define die( exitval )   do { \
  fprintf (stderr, "Aborting in %s at line %d\n", __FILE__, __LINE__); \
  bebop_exit (exitval); } while(0)


/**
 * Returns the number of bits in the binary representation of x.
 */
int
bebop_num_bits (const int x);


/**
 * @brief Positive integer product with overflow check.
 *
 * Positive integer product with overflow check.  If x*y overflows, returns 0,
 * without changing *result.  Otherwise returns nonzero, and sets *result = x*y.  
 *
 * @param x        [IN] positive integer
 * @param y        [IN] positive integer
 * @param result   [OUT] the product of x and y, if no overflow.
 *
 * @return         Zero if overflow, otherwise nonzero.
 */
int 
bebop_checked_positive_product (int* result, const int x, const int y);

/**
 * Returns NAN.  Hopefully this doesn't signal, if running under a 
 * Lisp's FFI (computing 0.0 / 0.0 may signal if the C function is
 * called by a Lisp's FFI, as the Lisp may change the floating-point
 * signaling properties).  We make this a function so we can call 
 * it under Lisp, Python, etc.
 */
double
make_nan ();

/**
 * Swaps x and y, each having sizeof == size, without using a temporary,
 * by using the XOR trick.  Swap granularity is 1 byte; perhaps some 
 * clever software pipelining / loop unrolling and use of vectorized
 * integer instructions could speed this up further, though if size is 
 * not very large, the required IF tests would slow things down.
 */
void
bebop_swap (void* x, void* y, size_t size);


#endif /* #ifndef _util_h */
