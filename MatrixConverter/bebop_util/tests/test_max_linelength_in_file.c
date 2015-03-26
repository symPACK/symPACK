/**
 * @file test_max_linelength_in_file.c
 * @author Mark Hoemmen
 * @since 28 Jul 2005
 * @date Time-stamp: <2008-07-16 10:21:52 mhoemmen>
 *
 * Driver program for testing max_linelength_in_file, in file.h.  Returns 
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
#include <bebop/util/file.h>
#include <bebop/util/random_number.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h> /* unlink */

static char digits[] = "0123456789";

static void
make_test_file_helper (FILE* f, const int maxlinelen, const int numlines)
{
  int i = 0;
  int special_line = -1;
  char* curline = bebop_malloc ((maxlinelen+1) * sizeof (char));
  assert (f != NULL);

  /* 
   * Make numlines lines of length a random number between 0 and maxlinelen. 
   * Fill them with random content that makes it easy to tell how long they
   * are by inspecting the file visually.  We also fill exactly one of the
   * lines with exactly maxlinelen characters so that we're guaranteed to 
   * have at least one line with that many characters.
   */
  special_line = bebop_random_integer (0, numlines - 1);
  assert (special_line >= 0 && special_line < numlines);
  for (i = 0; i < numlines; i++)
    {
      int curlen, j;

      if (i == special_line)
	curlen = maxlinelen;
      else
	curlen = bebop_random_integer (0, maxlinelen);

      assert (curlen >= 0 && curlen <= maxlinelen);
      
      for (j = 0; j < curlen; j++)
	curline[j] = digits[j % 10];

      curline[curlen] = '\0';
      fprintf (f, "%s\n", curline);
    }
  bebop_free (curline);
}

/*
 * Returns the name of the file
 */
static char*
make_test_file (const int maxlinelen, const int numlines)
{
  char filename[100];
  FILE* f = NULL;
  int file_desc = -1;

  /* tmpnam() opens a big security hole, so don't use it */
  strcpy (filename, "smc.tmp.XXXXXX");
  file_desc = mkstemp (filename);
  assert (file_desc != -1);
	
  bebop_log (2, "tempname: %s\n", filename);
  f = fdopen (file_desc, "w");
  assert (f != NULL);
  make_test_file_helper (f, maxlinelen, numlines);
  assert (0 == fclose (f));
  return strdup (filename);
}

static void
make_and_test_file (const int maxlinelen, const int numlines)
{
  char* filename = make_test_file (maxlinelen, numlines);
  FILE* f = NULL;
  
  bebop_log (2, "Got filename from make_test_file: %s\n", filename);

  f = fopen (filename, "r");
  assert (f != NULL);
  assert (max_linelength_in_file (f) == maxlinelen);
  assert (0 == fclose (f));
  /* delete the file */
  assert (0 == unlink (filename));

  if (filename != NULL) free (filename);
}



int 
main (int argc, char** argv)
{
  int info = 0;

  bebop_default_initialize (argc, argv, &info);
  assert (info == 0);
  bebop_log (1, "test_max_linelength_in_file: starting tests\n");

  make_and_test_file (10, 20);
  make_and_test_file (5, 50);
  make_and_test_file (15, 30);
  make_and_test_file (30, 15);
  make_and_test_file (40, 100);

  bebop_log (1, "test_max_linelength_in_file: passed all tests\n");
  bebop_exit (EXIT_SUCCESS);
  return EXIT_SUCCESS; /* to pacify the compiler */
}

