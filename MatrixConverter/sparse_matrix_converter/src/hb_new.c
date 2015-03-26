/**
 * @file hb_new.c
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-16 11:23:07 mhoemmen>
 *
 * A new and unfinished Harwell-Boeing file parser.
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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>





typedef
struct
{
  regex_t ifmt_regex;
  regex_t dfmt_regex;
} parser_data_t;

static parser_data_t parser_data;

static int parser_data_initialized_p = 0;

void
initialize_parser_data ()
{
  int errcode = 0;

  /* 
   * Matches the pattern "\((\d*)I(\d+)\)", in which \1 is the number
   * of integers per line (defaults to 1) and \2 is the (max) width of
   * the integer format.
   */
  errcode = regcomp (&(parser_data.ifmt_regex), 
		     "\([0-9]*[iI][0-9]+\)", 
		     REG_EXTENDED | REG_ICASE);
  assert (errcode == 0);
  errcode = regcomp (&(parser_data.dfmt_regex), 
		     "\([0-9]*[DEde][0-9]+\)", 
		     REG_EXTENDED | REG_ICASE);
  assert (errcode == 0);
  parser_data_initialized_p = 1;
}

void
free_parser_data ()
{
  if (parser_data_initialized_p)
    {
      regfree (&(parser_data.ifmt_regex));
      regfree (&(parser_data.dfmt_regex));
    }
  parser_data_initialized_p = 0;
}


int
parse_ifmt (int* num_elts_per_line, int* width, const char* fmt)
{
  int errcode = 0;
  regmatch_t pmatch[1];
  const size_t nmatch = (size_t) 1;
  char* ipos = NULL;
  char* ppos = NULL;
  size_t start, end;

  if (! parser_data_initialized_p)
    initialize_parser_data ();
  assert (parser_data_initialized_p);
  
  errcode = regexec (&(parser_data.ifmt_regex), fmt, nmatch, pmatch, 0);
  if (errcode == REG_NOMATCH)
    return -1;

  ipos = (char*) strchr ('I');
  if (ipos == NULL)
    {
      ipos == strchr ('i');
      if (ipos == NULL)
	{
	  fprintf (stderr, "*** parse_ifmt: Serious failure, failed to find "
		   "I or i after regex search said it was there! ***\n");
	  assert (0);
	}
    }
  assert ((unsigned int) ipos >= (unsigned int) fmt);
  end = (size_t) (ipos - fmt);
  start = (size_t) 1;

  if (start == end)
    {
      /* Number of elements per line not specified, so default to 1. */
      *num_elts_per_line = 1;
    }
  else
    {
      errcode = get_int_from_substring (num_elts_per_line, fmt, start, end);
      assert (errcode == 0);
      if (*num_elts_per_line < 1)
	{
	  fprintf (stderr, "*** parse_ifmt: Invalid number of elements per "
		   "line %d ***\n", *num_elts_per_line);
	  return -2;
	}
    }

  ppos = (char*) strchr (')');
  assert ((unsigned int) ppos > (unsigned int) ipos + 1);
  start = (size_t) (1 + ((unsigned int) ipos - (unsigned int) fmt));
  end = (size_t) ((unsigned int) ppos - (unsigned int) fmt) - 1;
  assert (end > start);
 
  errcode = get_int_from_substring (width, fmt, start, end);
  assert (errcode == 0);
  if (*width < 1)
    {
      fprintf (stderr, "*** parse_ifmt: Invalid width %d ***\n", *width);
      return -3;
    }

  return 0;
}







struct hb_header_t
{
  /* A72,A8 */
  char* title;
  char* key;

  /* 5I14 */
  int totcrd; /* total number of lines excluding header */
  int ptrcrd; /* number of lines for pointers */
  int indcrd; /* number of lines for row (or variable) indices */
  int valcrd; /* number of lines for numerical values */
  /* includes starting guesses and solution vectors if present */
  int rhscrd; /* number of lines for right-hand sides, including starting guesses and solution vectors if present; zero indicates no right-hand side data is present */

  /* A3,11X,4I14 */
  char* mxtype; /* matrix type */
  int nrow;  /* number of rows (or variables) */
  int ncol; /* number of columns or number of elements */
  int nnzero; /* number of row indices or variable indices (equal to number of entries for assembled matrix */
  int neltvl; /* number of elemental matrix entries; zero if assembled matrix */

  /* 2A16,2A20 */
  char* ptrfmt; /* format for pointers */
  char* indfmt; /* format for row (or variable) indices */
  char* valfmt; /* format for numerical values of coefficient matrix */
  char* rhsfmt; /* format for numerical values of right-hand sides */
  
  /* A3,11X,2I14 */ /* only present if there are RHSs present */
  char* rhstyp; /* rhstyp[0]: F for full storage or M for same format as matrix
		   rhstyp[1]: G if a starting vector(s) is supplied
		   rhstyp[2]: X if an exact solution vector(s) is supplied */
  int nrhs;     /* number of right-hand sides */
  int nrhsix;   /* number of row indices */
};

void
hb_header_clear (struct hb_header_t* header)
{
  if (header == NULL)
    return;

#define FREE_IF_VALID( p )  do { if((p) != NULL) free(p); } while(0)

  FREE_IF_VALID (header->title);
  FREE_IF_VALID (header->key);
  FREE_IF_VALID (header->mxtype);
  FREE_IF_VALID (header->ptrfmt);
  FREE_IF_VALID (header->indfmt);
  FREE_IF_VALID (header->valfmt);
  FREE_IF_VALID (header->rhsfmt);
  FREE_IF_VALID (header->rhstyp);
}


static int
get_int_from_substring (int* p_int, const char* buf, const int start, 
			const int end)
{
  char* s = get_substring (buf, start, end); 
  char* endptr = NULL; 
  int result = 0; 

  assert (s != NULL); 
  *p_int = strtol (s, &endptr, 10); 
  if (endptr == s) 
    { 
      free (s); 
      return errcode_on_convert;
    }
  free (s);
  return 0;
}


int
hb_header_read (struct hb_header_t* header, FILE* file)
{
#define VERIFIED_SUBSTR_GET( s, buf, start, end, errcode ) do { \
  s = get_substring(buf,start,end); \
  if (s == NULL) return errcode; \
} while(0)

#define VERIFIED_INT_SUBSTR_GET( p_int, buf, start, end, errcode ) do { \
  int e = get_int_from_substring ((p_int), buf, start, end); \
  if (e != 0) return errcode; \
} while(0)

  
  char buffer[82];
  int errcode = 0; 

  if (file == NULL)
    return -1;
  if (feof (file))
    return -1;

  if (fgets (buffer, 81, file) == NULL)
    return 1;

  /* extract title and key */
  VERIFIED_SUBSTR_GET (header->title, buffer, 0, 71, 1);
  VERIFIED_SUBSTR_GET (header->key, buffer, 73, 80, 1);

  /* extract line 2 data */
  if (fgets (buffer, 71, file) == NULL)
    return 2;

  VERIFIED_INT_SUBSTR_GET (&(header->ptrcrd), buffer, 28, 41, 2);
  VERIFIED_INT_SUBSTR_GET (&(header->indcrd), buffer, 42, 55, 2);
  VERIFIED_INT_SUBSTR_GET (&(header->valcrd), buffer, 56, 69, 2);
  
  /* extract line 3 data */
  if (fgets (buffer, 71, file) == NULL)
    return 3;
  
  VERIFIED_SUBSTR_GET (&(header->mxtype), buffer, 0, 2, 3);
  VERIFIED_INT_SUBSTR_GET (&(header->nrow), buffer, 14, 27, 3);
  VERIFIED_INT_SUBSTR_GET (&(header->ncol), buffer, 28, 41, 3);
  VERIFIED_INT_SUBSTR_GET (&(header->nnzero), buffer, 42, 55, 3);
  VERIFIED_INT_SUBSTR_GET (&(header-> neltvl), buffer, 56, 69, 3);

  /* extract line 4 data */
  if (fgets (buffer, 73) == NULL)
    return 4;

  VERIFIED_SUBSTR_GET (&(header->ptrfmt), buffer, 0, 15, 4);
  VERIFIED_SUBSTR_GET (&(header->indfmt), buffer, 16, 31, 4);
  VERIFIED_SUBSTR_GET (&(header->valfmt), buffer, 32, 51, 4);
  VERIFIED_SUBSTR_GET (&(header->rhsfmt), buffer, 52, 71, 4);

  /* Decide whether there is a line 5 in the header */
  if (header->rhscrd == 0)
    return;

  /* extract line 5 data */
  if (fgets (buffer, 43) == NULL)
    return 5;

  VERIFIED_SUBSTR_GET (&(header->rhstyp), buffer, 0, 2, 5);
  VERIFIED_INT_SUBSTR_GET (&(header->nrhs), buffer, 14, 27, 5);
  /* nrhsix, the number of row indices, is ignored in case of unassembled matrices */
  if (header->mxtype[2] == 'E' || header->mxtype[2] == 'e')
    {
      VERIFIED_INT_SUBSTR_GET (&(header->nrhsix), buffer, 28, 41, 5);
    }
  else
    header->nrhsix = 0;

  return 0;
}
