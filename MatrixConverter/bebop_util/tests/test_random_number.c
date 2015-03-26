/**
 * @file test_random_number.c
 * @author Mark Hoemmen
 * @since 7 Jul 2005
 * @date Time-stamp: <2008-07-16 10:21:58 mhoemmen>
 *
 * Driver program for testing functions defined in random_number.c.  Returns 
 * EXIT_SUCCESS if all tests pass, else prints an error message and returns 
 * EXIT_FAILURE.
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
#include <bebop/util/random_number.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void
test_random_number (int low, int high)
{
  const int total = high - low + 1;
  struct random_integer_from_range_without_replacement_generator_t* gen = NULL;
  int* returned = NULL;
  int errcode = 0;
  int p = 0;
  int i;

  if (low > high)
    {
      gen = create_random_integer_from_range_without_replacement_generator (low, high);
      if (gen->num_remaining != 0)
	{
	  bebop_fatal_error ("test:PRNG", "create_random_integer_from_range_"
			     "without_replacement_generator (%d, %d) should "
			     "return a generator with num_remaining == 0, bu"
			     "t num_remaining = %d instead", 
			     low, high, gen->num_remaining);
	}

      errcode = return_random_integer_from_range_without_replacement (&p, gen);
      if (errcode == 0)
	{
	  bebop_fatal_error ("test:PRNG", "return_random_integer_from_range_"
			     "without_replacement should always fail when "
			     "high < low, but did not fail");
	}

      destroy_random_integer_from_range_without_replacement_generator (gen);
      return;
    }

  returned = bebop_calloc (sizeof (int), total);
  gen = create_random_integer_from_range_without_replacement_generator (low, high);
  assert (gen != NULL);
  assert (gen->low == low);
  assert (gen->high == high);

  /* Try to get high-low+1 distinct integers within the range [low,high]. */

  for (i = 0; i < total; i++)
    {
      int j;

      errcode = return_random_integer_from_range_without_replacement (&p, gen);
      if (errcode != 0)
	{
	  bebop_fatal_error ("test:PRNG", "At iteration %d of %d, re"
		   "turn_random_integer_from_range_without_replacement should"
		   " have returned zero for success, but returned nonzero ins"
		   "tead", i+1, total);
	}
      if (p < low || p > high)
	{
	  bebop_fatal_error ("test:PRNG", "At iteration %d of %d, re"
		   "turn_random_integer_from_range_without_replacement return"
		   "ed a number %d outside the valid range [%d,%d]", 
		   i+1, total, p, gen->low, gen->high);
	}
      returned[i] = p;

      /* Make sure that p is unique (has not been chosen already). */
      for (j = 0; i < i-1; i++)
	{
	  if (p == returned[j])
	    {
	      bebop_fatal_error ("test:PRNG", "At iteration %d of %d"
		       ", return_random_integer_from_range_without_replacemen"
		       "t return ed a number %d that was returned previously "
		       "at iteration %d ***\n", i+1, total, p, j+1);
	    }
	}

      /* Make sure that all the remaining numbers are unique. */
      for (j = 1; j < gen->num_remaining; j++)
	{
	  int k;
	  for (k = 0; k < j; k++)
	    {
	      if (gen->remaining[j] == gen->remaining[k])
		{
		  bebop_fatal_error ("test:PRNG", "At iteration %d of "
			   "%d, gen->remaining has a duplicate entry!",
			   i+1, total);
		}
	    }
	}
    }

  if (gen->num_remaining != 0)
    {
      bebop_fatal_error ("test:PRNG", "After removing all %d numbers, "
			 "%d are still remaining, though 0 (zero) shou"
			 "ld remain!", total, gen->num_remaining);
    }

  /* Make sure that returning a random integer fails.  Try this twice. */

  errcode = return_random_integer_from_range_without_replacement (&p, gen);
  if (errcode == 0)
    {
      bebop_fatal_error ("test:PRNG", "return_random_integer_from_range_"
			 "without_replacement should fail after returning"
			 " all %d numbers, but did not fail!\n", 6);
    }

  errcode = return_random_integer_from_range_without_replacement (&p, gen);
  if (errcode == 0)
    {
      bebop_fatal_error ("test:PRNG", "return_random_integer_from_range_"
			 "without_replacement should fail after returning"
			 " all %d numbers, but did not fail!\n", 6);
    }

  destroy_random_integer_from_range_without_replacement_generator (gen);
  bebop_free (returned);
}


void
test_random_double (double low, double high)
{
  double d = bebop_random_double (low, high);

  bebop_log (2, "Random double in [%g,%g): %g\n", low, high, d);
  assert (d >= low);
  assert (d < high);
}


int 
main (int argc, char** argv)
{
  int low, high;
  int info = 0;

  bebop_default_initialize (argc, argv, &info); /* starts PRNG by default */
  assert (info == 0);
  bebop_log (1, "test_random_number: starting tests\n");

  if (argc > 1)
    {
      if (0 == strcmp (argv[1], "without-replacement"))
	{
	  low = 1;
	  high = 6;
	  test_random_number (low, high);

	  /* Now try a larger interval. */

	  low = -50;
	  high = 200;
	  test_random_number (low, high);

	  /* Now make a generator where high < low -- this one should
	     never return any values. */
	  low = 0;
	  high = -1;
	  test_random_number (low, high);
	}
      else if (0 == strcmp (argv[1], "random_double"))
	{
	  int i;

	  for (i = 0; i < 20; i++)
	    {
	      test_random_double (0.0, 1.0);
	      test_random_double (0.0, 1.0);
	      test_random_double (0.0, 1.0);
	    }

	  test_random_double (-10.0, 10.0);
	  test_random_double (3.14159, 7.654321);
	}
      else if (0 == strcmp (argv[1], "random_integer"))
	{
	  int i, die;
	  for (i = 0; i < 100; i++)
	    {
	      die = bebop_random_integer (1, 6);
	      assert (die >= 1);
	      assert (die <= 6);
	      bebop_log (2, "%d\n", die);
	    }
	}
    }
  else
    {
      low = 1;
      high = 6;
      test_random_number (low, high);

      /* Now try a larger interval. */

      low = -50;
      high = 200;
      test_random_number (low, high);

      /* Now make a generator where high < low -- this one should
	 never return any values. */
      low = 0;
      high = -1;
      test_random_number (low, high);

      test_random_double (0.0, 1.0);
      test_random_double (0.0, 1.0);
      test_random_double (0.0, 1.0);

      test_random_double (-10.0, 10.0);
      test_random_double (3.14159, 7.654321);
    }

  bebop_log (1, "test_random_number: passed all tests\n");
  bebop_exit (0);
  return 0; /* to pacify the compiler */
}

