/**
 * @file util.c
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-20 19:43:32 mhoemmen>
 *
 * Implementation of functions declared in util.h.
 * 
 * @note Originally called "smvm_util.c".  Renamed on 11 Nov 2007 
 * in the namespace renaming that was started 09 Nov 2007.
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
 *************************************************************************/
#include <bebop/util/config.h>

#include <bebop/util/util.h>
#include <bebop/util/random_number.h>

#include <assert.h>
#include <limits.h>    /* INT_MAX */
#include <stdlib.h>    /* qsort, ... */


int 
bebop_compare_double (const void *d1, const void *d2)
{
  double v_d1 = *((double *)d1);
  double v_d2 = *((double *)d2);

  if (v_d1 < v_d2) return -1;
  if (v_d1 > v_d2) return 1;
  return 0;
}

void 
bebop_get_median_min_max (double a[], int n, double *p_median, double *p_min, double *p_max)
{
  if (n < 1) return;
    
  qsort(a, n, sizeof(double), bebop_compare_double);
  *p_min = a[0];
  *p_max = a[n-1];
  *p_median = (n%2)? a[n/2] : (a[n/2-1] + a[n/2])/2;
}


void
bebop_init_vector_rand (int n, double* x)
{
  int i;
  for (i = 0; i < n; i++)
    x[i] = bebop_random_double (-1.0, 1.0);
}


void
bebop_init_vector_val (int n, double* x, double val)
{
  int i;
  for (i = 0; i < n; i++)
    x[i] = val;
}


double 
bebop_median (double* values, int n)
{
  if (n < 1)
    return 0.0;
  if (n < 2)
    return values[0];

  qsort (values, n, sizeof(double), bebop_compare_double);
  return (n%2) ? values[n/2] : (values[n/2-1] + values[n/2])/2;
}


int 
bebop_num_bits (const int x)
{
  int t = x;
  int L;
  for (L = 0; t; L++, t >>= 1)
    ; /* do nothing */

  return L;
}


int 
bebop_checked_positive_product (int* result, const int x, const int y)
{
  const int lintmax = bebop_num_bits (INT_MAX);
  const int lx = bebop_num_bits (x);
  const int ly = bebop_num_bits (y);
  int p;

  /* x*y requires (lx + ly + 1) bits.  Since INT_MAX is $2^{lintmax + 1} - 1$,
     any positive integer requiring the same number of bits as INT_MAX does not
     overflow. */

  if (lx + ly + 1 > lintmax)
    return 0; /* overflowed */

  /* else */
  p = x * y;

  /* Sanity check to make sure that there was no overflow. */
  if (p < 0)
    return 0; /* overflowed */

  /* else */
  *result = p;
  return 1;
}

double
make_nan ()
{
  /* Note that this works even without C99 features */
  union {
    double d;
    unsigned int u[2];
  } x;
  x.u[0] = 0x7ff80000;
  x.u[1] = 0;
  return x.d; 
}

void
bebop_swap (void* x, void* y, size_t size)
{
  char* __x = (char*) x;
  char* __y = (char*) y;
  int i;

  for (i = 0; i < size; i++)
    {
      char xx = __x[i];
      char yy = __y[i];

      xx ^= yy;
      yy ^= xx;
      xx ^= yy;

      __x[i] = xx;
      __y[i] = yy;
    }
}    


