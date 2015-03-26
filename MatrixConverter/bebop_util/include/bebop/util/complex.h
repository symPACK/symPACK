#ifndef __complex_h
#define __complex_h
/**
 * @file __complex.h
 * @author Mark Hoemmen
 * @since 2 June 2005
 * @date Time-stamp: <2008-07-20 20:30:30 mhoemmen>
 *
 * Replaces "double _Complex" for compilers not supporting the C99 
 * standard.  Defines an interface for simple complex arithmetic 
 * that lets us use the same code for C99-compliant compilers and 
 * non-compliant compiler.
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


#ifndef HAVE_COMPLEX_H
/* #  warning "You seem to lack the C99 header file complex.h.  I will provide a stub complex datatype in its place." */

typedef struct 
{
  double real;
  double imag;
} double_Complex;

const double_Complex double_Complex_I;
const double_Complex double_Complex_ZERO;

double_Complex 
double_Complex_add (const double_Complex a, const double_Complex b);

double_Complex 
double_Complex_subtract (const double_Complex a, const double_Complex b);

double_Complex 
double_Complex_negate (const double_Complex a);

double_Complex 
double_Complex_multiply (const double_Complex a, const double_Complex b);

int 
double_Complex_not_equal (const double_Complex a, const double_Complex b);

int 
double_Complex_equal (const double_Complex a, const double_Complex b);

double_Complex
new_double_Complex (const double a, const double b);

double 
double_Complex_real_part (const double_Complex a);

double 
double_Complex_imag_part (const double_Complex a);

double_Complex 
double_Complex_conj (const double_Complex a);

double
double_Complex_cabs (const double_Complex a);

double_Complex
double_Complex_divide (const double_Complex a, const double_Complex b);

#else /* we do have C99's complex.h */
#  include <complex.h>
#  include <math.h> /* needed for isnan() */

/*
 * The following typedef and defines allow us to produce source files
 * that are compatible with both the C99 - standard implementation of
 * complex arithmetic, and our own implementation, at no performance
 * cost in the C99 case.
 */
typedef double _Complex  double_Complex;

#  define double_Complex_add(a,b)        ((a) + (b))
#  define double_Complex_subtract(a,b)   ((a) - (b))
#  define double_Complex_negate(a)       (-(a))
#  define double_Complex_multiply(a,b)   ((a) * (b))
#  define double_Complex_not_equal(a,b)  ((a) != (b))
#  define double_Complex_equal(a,b)      ((a) == (b))
#  define new_double_Complex(a,b)        ((a) + I*(b))
#  define double_Complex_real_part(a)    (creal(a))
#  define double_Complex_imag_part(a)    (cimag(a))
#  define double_Complex_conj(a)         (conj(a))
#  define double_Complex_I               (I)
#  define double_Complex_ZERO            (0.0 + 0.0*I)
#  define double_Complex_cabs            (cabs(a))
#  define double_Complex_divide(a,b)     ((a) / (b))

int
double_Complex_isnan (const double_Complex a);


#endif /* HAVE_COMPLEX_H */

#endif /* __complex_h */
