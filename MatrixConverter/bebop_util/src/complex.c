/**
 * @file complex.c
 * @author Mark Hoemmen
 * @since 02 June 2005
 * @date Time-stamp: <2008-07-20 20:29:27 mhoemmen>
 *
 * Replaces "double _Complex" for compilers not supporting the C99 standard.
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
#include <bebop/util/complex.h>
#include <bebop/util/util.h>
#include <math.h>

#ifndef HAVE_COMPLEX_H
const double_Complex double_Complex_I = {0.0, 1.0};
const double_Complex double_Complex_ZERO = {0.0, 0.0};

double_Complex 
double_Complex_add (const double_Complex a, const double_Complex b)
{
  double_Complex c;
  c.real = a.real + b.real;
  c.imag = a.imag + b.imag;
  return c;
}

double_Complex 
double_Complex_negate (const double_Complex a)
{
  double_Complex c;
  c.real = -(a.real);
  c.imag = -(a.imag);
  return c;
}

double_Complex 
double_Complex_subtract (const double_Complex a, const double_Complex b)
{
  return double_Complex_add (a, double_Complex_negate (b));
}

double_Complex 
double_Complex_multiply (const double_Complex a, const double_Complex b)
{
  double_Complex c;
  c.real = a.real * b.real - a.imag * b.imag;
  c.imag = a.real * b.imag + a.imag * b.real;
  return c;
}

double_Complex
double_Complex_divide (const double_Complex a, const double_Complex b)
{
  /*
   * We do some scaling to avoid overflow.  However, there still 
   * could be problems when scaling a.  I need to do a full analysis.
   *
   * real(a/b) = (a.real/M)*(b.real/M) + (a.imag/M)*(b.imag/M)
   *             ---------------------------------------------
   *             (b.real/M)*(b.real/M) + (b.imag/M)*(b.imag/M)
   */
  double_Complex c;
  double_Complex acopy = new_double_Complex (a.real, a.imag);
  double_Complex bcopy = new_double_Complex (b.real, b.imag);
  /* scaling factor */

  const double M = MAX (fabs (b.real), fabs (b.imag));
  double d;

  bcopy.real = bcopy.real / M;
  bcopy.imag = bcopy.imag / M;
  acopy.real = acopy.real / M;
  acopy.imag = acopy.imag / M;
  d = bcopy.real * bcopy.real + bcopy.imag * bcopy.imag;

  c.real = (acopy.real * bcopy.real + acopy.imag * bcopy.imag) / d;
  /* 
   * FIXME: this might produce NaN instead of Inf if a and b are 
   * finite (and not NaN) and a/b overflows.
   */
  c.imag = (acopy.imag * bcopy.real - bcopy.imag * acopy.real) / d; 

  /*
  c.real = (a.real * b.real + a.imag * b.imag) / d;
  c.imag = (a.imag * b.real - b.imag * a.real) / d; 
  */
  return c; 
}

int 
double_Complex_not_equal (const double_Complex a, const double_Complex b)
{
  return !(a.real == b.real && a.imag == b.imag);  
}

int 
double_Complex_equal (const double_Complex a, const double_Complex b)
{
  return (a.real == b.real && a.imag == b.imag);  
}

double_Complex
new_double_Complex (const double a, const double b)
{
  double_Complex z;
  z.real = a;
  z.imag = b;
  return z;
}

double 
double_Complex_real_part (const double_Complex a)
{
  return a.real;
}

double 
double_Complex_imag_part (const double_Complex a)
{
  return a.imag;
}

double_Complex 
double_Complex_conj (const double_Complex a)
{
  double_Complex c;
  c.real = a.real;
  c.imag = -(a.imag);
  return c;
}

double
double_Complex_cabs (const double_Complex a)
{
  /* We scale the computation to avoid overflow. */
  double re = double_Complex_real_part (a);
  double im = double_Complex_imag_part (a);
  double d;

  /* We don't use the MAX macro because that evaluates fabs() 
   * once unnecessarily */
  {
    const double t1 = fabs (re);
    const double t2 = fabs (im);
    d = (t1 > t2) ? t1 : t2;
  }

  if (d == 0.0)
    return d;
  else
    {
      re = re / d;
      im = im / d;
      return d * sqrt (re*re + im*im);
    }
}



#endif /* HAVE_COMPLEX_H */


/*
 * mfh 20 Jul 2008: we always define this as a function, even if we
 * have the C99 headers.
 */
int
double_Complex_isnan (const double_Complex a)
{
#ifdef HAVE_ISNAN
  const double x = creal(a);
  const double y = cimag(a);

  return isnan (x) || isnan (y);
#else
  const double x = a.real;
  const double y = a.imag;

  /*
   * FIXME (mfh 20 Jul 2008): some (bad) compilers may optimize away "x
   * != x" to false, even though it's a correct way to test if x is NaN
   * (NaN is not equal to any number, including itself).  Also, the test
   * could be slower than bit comparisons, if the hardware invokes a
   * software interrupt to handle the NaN comparison, rather than
   * treating the NaN as a raw bit string.  However, C99 was supposed to
   * fix these issues with its isnan() predicate, and I don't feel like
   * redoing all the work that the C99 implementers did -- you should go
   * get a C99-compliant system if you want that!
   */
  return (x != x) || (y != y);
#endif /* HAVE_ISNAN */
}

