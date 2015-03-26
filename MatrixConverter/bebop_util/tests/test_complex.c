/**
 * @file test_complex.c
 * @author Mark Hoemmen
 * @since 11 July 2006
 * @date Time-stamp: <2008-07-16 10:21:24 mhoemmen>
 * 
 * Correctness tests for complex numbers.  Executable returns 0 if all
 * the tests pass, otherwise it breaks with an assertion.
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
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "test_includes.h"
#include <bebop/util/complex.h>


#define ASSERT_EQUAL( __z, __a, __b )  do {     \
  seqnum++;                                     \
  if ((double_Complex_real_part(__z) != __a) || \
      (double_Complex_imag_part(__z) != __b))   \
    {                                           \
      bebop_error ("test:complex", "failed test number %d, an equality test", seqnum); \
      bebop_exit (EXIT_FAILURE); \
    } \
  } while(0)

#define ASSERT_NOT_EQUAL( __w, __z )  do {      \
  seqnum++;                                     \
  if (! double_Complex_not_equal ((__w), (__z))) \
    {                                            \
      bebop_error ("test:complex", "failed test number %d, an inequality test", seqnum); \
      bebop_exit (EXIT_FAILURE); \
    } \
  } while(0)

int
main (int argc, char** argv)
{
  double_Complex a, b, c;
  int seqnum = 0; /* current test number */
  int info = 0;

  bebop_default_initialize (argc, argv, &info);
  assert (info == 0);
  bebop_log (1, "test_complex: starting tests\n");

  a = new_double_Complex (1.0, 2.0);
  ASSERT_EQUAL( a, 1.0, 2.0 ); 

  b = new_double_Complex (3.0, -3.0);
  ASSERT_EQUAL( b, 3.0, -3.0 );

  c = double_Complex_add (a, b);
  ASSERT_EQUAL( c, 4.0, -1.0 );

  c = double_Complex_divide (c, new_double_Complex(1.0, 0.0));
  ASSERT_EQUAL( c, 4.0, -1.0 );

  c = double_Complex_divide (c, new_double_Complex(4.0, 0.0));
  ASSERT_EQUAL( c, 1.0, -0.25 );

  c = double_Complex_conj (c);
  ASSERT_EQUAL( c, 1.0, +0.25 );

  c = double_Complex_negate (c);
  ASSERT_EQUAL( c, -1.0, -0.25 );

  ASSERT_NOT_EQUAL( c, a );
  ASSERT_NOT_EQUAL( c, b );
  ASSERT_NOT_EQUAL( a, b );
  ASSERT_EQUAL( double_Complex_ZERO, 0.0, 0.0 );

  bebop_log (1, "test_complex: passed all tests\n");
  bebop_exit (EXIT_SUCCESS);
  return EXIT_SUCCESS;
}

