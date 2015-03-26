/**
 * @file test_extstring.c
 * @author Mark Hoemmen
 * @since 05 May 2007
 * @date Time-stamp: <2008-07-16 10:21:32 mhoemmen>
 * 
 * Correctness tests for extstring_t.
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
#include <string.h>

#include "test_includes.h"
#include <bebop/util/extstring.h>
#include <bebop/util/random_number.h>


static inline char
random_char ()
{
  return 'a' + (char) bebop_random_integer (0, 25);
}


static void
test_pushchars (extstring_t* str, int numchars)
{
  char *s1, *s2;
  int i;

  bebop_log (1, "test_pushchars: starting tests\n");
  assert (numchars >= 0);
  s2 = bebop_malloc ((numchars + 1) * sizeof (char));

  extstring_clear (str);
  for (i = 0; i < numchars; i++)
    {
      char c = random_char ();
      extstring_pushchar (str, c);
      s2[i] = c;
    }
  s2[numchars] = '\0';

  /* Now compare the extstring with s2 */
  s1 = extstring_string (str);
  if (s1 == NULL)
    {
      bebop_error ("test:extstring:pushchars", "extstring_string() returned NULL unexpectedly");
      bebop_exit (EXIT_FAILURE);
    }
  bebop_log (2, "String 1: %s\n", s1);
  bebop_log (2, "String 2: %s\n", s2);
  if (0 != strcmp (s1, s2))
    {
      bebop_error ("test:extstring:pushchars", "s1 and s2 are not equal");
      bebop_exit (EXIT_FAILURE);
    }
  bebop_free (s2);
  bebop_log (1, "test_pushchars: passed all tests\n");
}

static void
iter_test_pushchars (int numiters)
{
  extstring_t s;
  int i;

  bebop_log (1, "iter_test_pushchars: starting tests\n");
  bebop_log (2, "numiters = %d\n", numiters);
  assert (numiters >= 0);
  extstring_init (&s, bebop_random_integer (0, 50));

  for (i = 0; i < numiters; i++)
    {
      int num_to_push = bebop_random_integer (0, 300);
      bebop_log (2, "\tPushing %d characters\n");
      test_pushchars (&s, num_to_push);
    }

  extstring_deinit (&s);
  bebop_log (1, "iter_test_pushchars: passed all tests\n");
}


int
main (int argc, char** argv)
{
  int info = 0;

  bebop_default_initialize (argc, argv, &info);
  assert (info == 0);

  bebop_log (1, "extstring_t: starting tests\n");
  iter_test_pushchars (100);
  bebop_log (1, "extstring_t: passed all tests\n");

  bebop_exit (EXIT_SUCCESS);
  return EXIT_SUCCESS;
}

