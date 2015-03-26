/**
 * @file test_split.c
 * @author Mark Hoemmen
 * @since 03 May 2007
 * @date Time-stamp: <2008-07-16 10:22:27 mhoemmen>
 * 
 * Correctness tests for split_t and manysplit_t.  This test requires
 * visual inspection of the results.
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
#include <bebop/util/extstring.h>
#include <bebop/util/split.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

void
test_manysplit (FILE* infile)
{
  const char delimiters[] = " \t";
  const int ndelimiters = 2;
  manysplit_t spl;
  extstring_t line;

  manysplit_init (&spl, 0);
  extstring_init (&line, 40);

  while (! feof (infile))
    {
      list_t list;
      int count = 0;
      char* theline = NULL;
      int linelen = 0;

      linelen = bebop_safe_getline (&line, infile);
      assert (linelen >= 0);
      theline = extstring_string (&line);
      if (linelen == 0)
	bebop_log (3, "test_manysplit: got empty line\n");
      else
	bebop_log (3, "test_manysplit: got length %d line: %s\n", 
		  linelen, theline);
      manysplit (&spl, theline, delimiters, ndelimiters);
      for (list = manysplit_tokens (&spl); 
	   ! list_empty_p (list); 
	   list = list_cdr (list), count++, printf("\n"))
	{
	  const char* str = (const char*) list_car (list);
	  assert (str != NULL);
	  printf ("Token %d: %s ", count, str);
	}
    }
  extstring_deinit (&line);
  manysplit_deinit (&spl);
}


int
main (int argc, char** argv)
{
  int info = 0;

  bebop_default_initialize (argc, argv, &info);
  assert (info == 0);

  if (argc < 2)
    return 0; /* trivial test */
  else
    {
      FILE* f = NULL;
      int i;

      /* Assume that the command-line arguments are input files.  Test
	 each one in turn by running manysplit on its lines. */
      for (i = 1; i < argc; i++)
	{
	  f = fopen (argv[i], "r");
	  if (f == NULL)
	    printf ("\nFile %s cannot be opened\n", argv[i]);
	  else
	    {
	      printf ("\nTesting manysplit() on file %s:\n", argv[i]);
	      test_manysplit (f);
	      (void) fclose (f);
	    }
	}
    }

  bebop_exit (EXIT_SUCCESS);
  return EXIT_SUCCESS;
}

