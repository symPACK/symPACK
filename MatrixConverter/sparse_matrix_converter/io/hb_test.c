/**
 * @file hb_test.c
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-16 10:59:00 mhoemmen>
 *
 * Preliminary, uncompleted new Harwell-Boeing file parser in C.
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


#include "hb_tok.h"
#include "hb_parser.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


static void
print_token_stream (const char s[], const size_t lenmax)
{
  hb_token_t tok;
  int error = 0;
  char* buf = NULL;
  size_t buflen = 0;

  for (tok = hb_token_start (s, lenmax, &error);
       ! hb_token_endp (tok);
       tok = next_hb_token (tok, &error))
    {
      if (hb_token_error (error))
	{
	  fprintf (stderr, "Got error %d after next_hb_token()\n", error);
	  goto print_token_stream_cleanup;
	}
      read_hb_token_with_resize (&buf, &buflen, tok, &error);
      if (hb_token_error (error))
	{
	  fprintf (stderr, "Got error %d after read_hb_token_with_"
		   "resize()\n", error);
	  goto print_token_stream_cleanup;
	}
      else
	printf ("Got token %s of type %s\n", buf, hb_token_type_string(tok));
    }

 print_token_stream_cleanup:
  if (buf != NULL)
    free (buf);
}


int
main (int argc, char** argv)
{
  hb_format_t fmt;
  int error = 0;
  char buf[100];
  size_t buflen = 99;

  /*  if (argc > 1)
    print_token_stream (argv[1], strlen(argv[1]));
  */

  parse_hb_format (&fmt, argv[1], strlen(argv[1]), buf, &buflen, &error);
  print_hb_format (stderr, fmt);
  

  return 0;
}
