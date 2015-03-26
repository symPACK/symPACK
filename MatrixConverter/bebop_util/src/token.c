/**
 * @file token.c
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-16 10:12:40 mhoemmen>
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
#include <bebop/util/string.h>

/*
 * Setup for using getline().  Must do this before including <stdio.h>.
 */
#ifdef HAVE_GETLINE
#  define _GNU_SOURCE
#endif /* HAVE_GETLINE */

#include <stdio.h>
#include <bebop/util/extstring.h>

static int
bebop_char_member (const char c, const char delims[], const int ndelims)
{
  int i;
  for (i = 0; i < ndelims; i++)
    if (delims[i] == c)
      return 1;
  return 0;
}

static void
bebop_ungetc (FILE* f, char c)
{
  if (EOF == ungetc (c, f))
    {
      bebop_error ("io:read", "bebop_ungetc: ungetc() failed");
      bebop_exit (EXIT_FAILURE);
    }
}

static int
bebop_pass_delimiters (FILE* f, const char delims[], const int ndelims)
{
  if (feof (f))
    return -1;

  for (c = getc (f);
       ! feof (f) && char_member (c, delims, ndelims);
       c = getc (f))
    ;
  if (feof (f))
    return -1;
  bebop_ungetc (f, c);
  return 0;
}

/**
 * Get the next token from the file, using the delimiters to separate
 * tokens.  Don't cross newline boundaries.  
 *
 * @return 0 if successfully got a token, else nonzero.
 */
static int
next_token (extstring_t* token, 
	    FILE* f, 
	    const char delimiters[], 
	    const int ndelimiters)
{
  int retval = 0;
  char c;

  if (feof (f))
    return -1;
  retval = bebop_pass_delimiters (f, delimiters, ndelimiters);
  
  


  while (tosplit[i] != '\0')
    {
      int len = 0;
      int start = i;

      /* Get past the current range of delimiters */
      while (tosplit[i] != '\0' && 
	     char_member (tosplit[i], delimiters, ndelimiters))
	i++;
      if (tosplit[i] == '\0')
	return spl;
      start = i;
      /* Now tosplit[start] is a non-delimiter character */

      /* Get to the end of the delimiter characters */
      while (tosplit[i] != '\0' &&
	     ! char_member (tosplit[i], delimiters, ndelimiters))
	i++;
      len = start + i; /* doesn't matter if tosplit[i] is '\0' */
      spl.tokens = list_append_item (spl.tokens, 
				     bebop_strdup_range (tosplit, start, len));
      spl.ntokens++;
    }
  return spl;

  ...
}

/**
 * Try to get "desired_ntokens" number of tokens from the given
 * filestream f.  Store them (destructively) in the given list
 * "tokens" of dynamically created extstring_t objects.  Store in
 * "gotten_ntokens" the actual number of tokens that we got. 
 */
int
bebop_get_tokens (list_t* tokens, 
		 int* gotten_ntokens, 
		 FILE* f,
		 const int desired_ntokens,
		 const char delimiters[],
		 const int ndelimiters)
{
  list_t current = *tokens;
  list_t outlist = *tokens;

  if (desired_ntokens < 0)
    {
      bebop_error ("bebop_get_tokens", 
		  "desired_ntokens = %d < 0", 
		  desired_ntokens);
      return -1;
    }
  for (*gotten_ntokens = 0; 
       *gotten_ntokens < desired_ntokens;
       (*gotten_ntokens)++)
    {
      extstring_t* token = NULL;
      int result = 0;

      if (list_empty_p (current))
	  /* We're out of nodes to reuse, so we have to create a new
	     extstring_t */
	token = extstring_create (0);
      else
	token = list_car (current);

      /* Try to get the next token */
      result = bebop_next_token (token, f, delimiters, ndelimiters);
      if (result != 0)
	{
	  extstring_destroy (token);
	  break;
	}
      else
	{
	  if (list_empty_p (current))
	    outlist = list_append_node (outlist, token);
	  else
	    {
	      list_set_car (current, token);
	      current = list_cdr (current);
	    }
	}
    }
  return 0;
}


static int
bebop_detect_newline (FILE* f)
{
  char c;

  if (feof (f))
    return -1; /* end-of-file */

  c = getc (f);
  if (c == EOF)
    return -1; /* end-of-file */
  else if (c == '\n')
    return 1; /* *nix - style endline */
  else if (c == '\r')
    {
      char next;
      
      /* Macs have "\r" as the endline character, whereas DOS (and
	 Windows) uses "\r\n".  We have to check the next character to
	 see if it's an "\n".  Luckily we're allowed at least one
	 pushback onto the stream, so we can put the character back if
	 it's not an '\n'.  We should first, of course, check for
	 end-of-file. */
      
      if (feof (f))
	return -1; /* end-of-file */
	  
      next = getc (infile);
      if (next == '\n')
	return 2; /* DOS/Windows endline */
      else
	{
	  /* It's a Mac endline; put back the character.  A single
	     putback isn't supposed to fail (according to the POSIX
	     standard), but we check anyway and report a fatal error
	     if it does fail. */
	  if (EOF == ungetc (next, f))
	    {
	      bebop_error ("io:read", "bebop_detect_newline:"
			  "ungetc() failed on only one call; not "
			  "supposed to happen!");
	      bebop_exit (EXIT_FAILURE);
	    }
	  return 3; /* Mac endline */
	}
    }
  return 0; /* no endline */
}		 

void
bebop_advance_past_newline (FILE* f)
{
  while (0 != bebop_detect_newline (f))
    ;
}
      







/*

ssize_t 
getline(char **lineptr, size_t *n, FILE *stream);

Format for ssize_t is %zu -- z refers specifically to a size_t or ssize_t argumet.  Obviously this is a GNU extension (???).
*/
