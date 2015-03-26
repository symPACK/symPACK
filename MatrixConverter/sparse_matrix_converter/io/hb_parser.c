/**
 * @file hb_parser.c
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-16 10:58:30 mhoemmen>
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

#include "hb_parser.h"
#include "hb_tok.h"

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>


static const int dbg = 1;
size_t indent = 0;

static void
do_indent ()
{
  size_t i;
  for (i = 0; i < indent; i++)
    fprintf (stderr, " ");
}

hb_format_t
make_new_hb_format ()
{
  hb_format_t fmt;

  if (dbg)
    {
      do_indent ();
      fprintf (stderr, "--- make_new_hb_format ---\n");
    }
  fmt.count = 0;
  fmt.type = 'A';
  fmt.params[0] = 0;
  fmt.params[1] = 0;
  return fmt;
}


static hb_token_t 
consume_token (hb_token_t token, hb_token_type_t type, int* error)
{
  if (dbg)
    {
      do_indent ();
      fprintf (stderr, "--- consume_token ---\n");
    }

  *error = 0;
  if (token.type != type)
    {
      *error = FMT_PARSE_ERROR;
      return token;
    }
  else
    {
      if (dbg)
	indent += 2;
      token = next_hb_token (token, error);
      if (dbg)
	indent -= 2;
      return token;
    }
}

static hb_token_t 
parse_uint (size_t* u, hb_token_t token, char* buf, size_t* buflen, int* error)
{
  size_t temp = -1; /* Helps with debugging */

  if (dbg)
    {
      do_indent ();
      fprintf (stderr, "--- parse_uint ---\n");
    }

  if (token.type != UINT)
    {
      *error = FMT_PARSE_ERROR;
      return token;
    }
  
  if (dbg)
    indent += 2;
  read_hb_token (buf, buflen, token, error);
  if (dbg)
    indent -= 2;
  if (*error != FMT_NO_ERROR)
    {
      if (dbg)
	{
	  do_indent ();
	  fprintf (stderr, "*** parse_uint: read_hb_token failed ***\n");
	}
      return token;
    }

  errno = 0;
  temp = strtoul (buf, NULL, 10);
  if (errno != 0)
    {
      errno = 0;
      *error = FMT_OTHER_ERROR;
      if (dbg)
	{
	  do_indent ();
	  fprintf (stderr, "*** parse_uint: strtoul() failed ***\n");
	}
      return token;
    }
  *u = temp;
  if (dbg)
    {
      do_indent ();
      fprintf (stderr, "parse_uint: got value %zu\n", *u);
    }
  *error = 0;
  if (dbg)
    indent += 2;
  token = next_hb_token (token, error);
  if (dbg)
    indent -= 2;
  return token;
}

static hb_token_t 
parse_letter (char* c, hb_token_t token, int* error)
{
  if (dbg)
    {
      do_indent ();
      fprintf (stderr, "--- parse_letter ---\n");
    }

  if (token.type != LETTER)
    {
      if (dbg)
	{
	  do_indent ();
	  fprintf (stderr, "*** parse_letter: token.type == %s instead of LETTER ***\n",
		   hb_token_type_string(token));
	}
      *error = FMT_PARSE_ERROR;
      return token;
    }

  *c = hb_token_first_char (token);
  *error = FMT_NO_ERROR;

  /* Canonicalize the letter by converting to uppercase */
  *c = toupper (*c);

  if (dbg)
    {
      do_indent ();
      fprintf (stderr, "Got letter %c\n", *c);
    }

  if (dbg)
    indent += 2;
  token = next_hb_token (token, error);
  if (dbg) 
    indent -= 2;
  return token;
}


static hb_token_t
parse_period_uint (size_t* u, hb_token_t token, 
		   char* buf, size_t* buflen, int* error)
{
  if (dbg)
    {
      do_indent ();
      fprintf (stderr, "--- parse_period_uint ---\n");
    }

  if (dbg)
    indent += 2;
  token = consume_token (token, PERIOD, error);
  if (dbg)
    indent -= 2;
  if (*error != FMT_NO_ERROR)
    return token;
  else
    {
      if (dbg)
	indent += 2;
      token = parse_uint (u, token, buf, buflen, error);
      if (dbg)
	indent -= 2;
      return token;
    }
  
  /*
  else if (token.type == UINT)
   return parse_uint (u, token, buf, buflen, error);
  else
    {
      *u = 0;
      return token;
    }
  */
}


static hb_token_t
parse_integer_format (hb_format_t* format, hb_token_t token, 
		      char* buf, size_t* buflen, int* error)
{
  if (dbg)
    {
      do_indent ();
      fprintf (stderr, "--- parse_integer_format ---\n");
    }

  /* Read field width specifier */
  if (dbg)
    indent += 2;
  token = parse_uint (&(format->field_width), token, buf, buflen, error);
  if (dbg)
    indent -= 2;
  if (*error != FMT_NO_ERROR)
    return token;

  /* If next token is a period, then get a minimum number of digits value */
  if (token.type == PERIOD)
    {
      if (dbg)
	indent += 2;
      token = parse_period_uint (&(format->params[0]), token,
				 buf, buflen, error);
      if (dbg)
	indent -= 2;
      return token;
    }
  else
    return token;
}

static hb_token_t
parse_fixed_format (hb_format_t* format, hb_token_t token, 
		    char* buf, size_t* buflen, int* error)
{
 /* Read field width specifier */
  token = parse_uint (&(format->field_width), token, buf, buflen, error);
  if (*error != FMT_NO_ERROR)
    return token;

  /* If next token is a period, then get a digits after decimal point value */
  if (token.type == PERIOD)
    return parse_period_uint (&(format->params[0]), token,
			      buf, buflen, error);
  else
    return token;
}


static hb_token_t
parse_float_format (hb_format_t* format, hb_token_t token, 
		    char* buf, size_t* buflen, int* error)
{
  if (dbg)
    {
      do_indent ();
      fprintf (stderr, "--- parse_integer_format ---\n");
    }

  /* Read field width specifier */
  if (dbg)
    indent += 2;
  token = parse_uint (&(format->field_width), token, buf, buflen, error);
  if (dbg)
    indent -= 2;
  if (*error != FMT_NO_ERROR)
    return token;

  /* You must have a period then a significand length value */
  if (dbg)
    indent += 2;
  token = parse_period_uint (&(format->params[0]), token,
			     buf, buflen, error);
  if (dbg)
    indent -= 2;
  if (*error != FMT_NO_ERROR)
    return token;

  /* You may optionally have an E and then an exponent length */
  if (token.type == LETTER && hb_token_first_char(token) == 'E')
    {
      if (dbg)
	indent += 2;
      token = consume_token (token, LETTER, error);
      if (dbg)
	indent -= 2;
      if (*error != FMT_NO_ERROR)
	return token;
      if (dbg)
	indent += 2;
      token = parse_uint (&(format->params[1]), token,
			  buf, buflen, error);
      if (dbg)
	indent -= 2;
    }

  return token;
}


static hb_token_t
parse_datatype_dispatch (hb_format_t* format, hb_token_t token, 
			 char* buf, size_t* buflen, int* error)
{
  char c;

  if (token.type != LETTER)
    {
      *error = FMT_PARSE_ERROR;
      return token;
    }
  c = hb_token_first_char (token);

  if (c == 'D' || c == 'E' || c == 'G')
    return parse_float_format (format, token, buf, buflen, error);
  else if (c == 'I')
    return parse_integer_format (format, token, buf, buflen, error);
  else if (c == 'F')
    return parse_fixed_format (format, token, buf, buflen, error);
  else
    {
      *error = FMT_PARSE_ERROR;
      return token;
    }
}


static hb_token_t
parse_format_string (hb_format_t* format, hb_token_t token, 
		     char* buf, size_t* buflen, int* error)
{
  if (dbg)
    {
      do_indent ();
      fprintf (stderr, "--- parse_format_string ---\n");
    }

  if (token.type == UINT)
    {
      if (dbg)
	indent += 2;
      token = parse_uint (&(format->count), token, buf, buflen, error);
      if (dbg)
	indent -= 2;
      if (*error != FMT_NO_ERROR)
	{
	  if (dbg)
	    {
	      do_indent ();
	      fprintf (stderr, "*** parse_format_string: failed after first parse_uint() ***\n");
	    }
	  return token;
	}
    }
  else
    {
      if (dbg)
	{
	  do_indent ();
	  fprintf (stderr, "Setting default count\n");
	}
      format->count = 1; /* Set default */
    }

  if (dbg)
    indent += 2;
  token = parse_letter (&(format->type), token, error);
  if (dbg)
    indent -= 2;
  if (*error != FMT_NO_ERROR)
    {
      if (dbg)
	{
	  do_indent ();
	  fprintf (stderr, "*** parse_format_string: failed after parse_letter() ***\n");
	}
      return token;
    }

  /* Dispatch to the right kind of format type */
  if (dbg)
    indent += 2;
  token = parse_datatype_dispatch (format, token, buf, buflen, error);
  if (dbg)
    indent -= 2;
  
  return token;
}


void
parse_hb_format (hb_format_t* format, 
		 const char s[],
		 const size_t lenmax,
		 char* buf, 
		 size_t* buflen, int* error)
{
  hb_token_t token;

  *format = make_new_hb_format ();

  if (dbg)
    indent += 2;
  token = hb_token_start (s, lenmax, error);
  if (dbg)
    indent -= 2;
  if (*error != FMT_NO_ERROR)
    return;

  if (token.type == LPAREN)
    {
      if (dbg)
	indent += 2;
      token = consume_token (token, LPAREN, error);
      if (dbg)
	indent -= 2;
      if (*error != FMT_NO_ERROR)
	{
	  if (dbg)
	    {
	      do_indent ();
	      fprintf (stderr, "*** parse_hb_format: failed after first consume_token() ***\n");
	    }
	  return;
	}
      if (dbg)
	indent += 2;
      token = parse_format_string (format, token, buf, buflen, error);
      if (dbg)
	indent -= 2;
      if (*error != FMT_NO_ERROR)
	{
	  if (dbg)
	    {
	      do_indent ();
	      fprintf (stderr, "*** parse_hb_format: failed after parse_format_string() ***\n");
	    }
	  return;
	}
      if (dbg)
	indent += 2;
      (void) consume_token (token, RPAREN, error);
      if (dbg) 
	indent -= 2;
    }
  else
    {
      if (dbg)
	indent += 2;
      (void) parse_format_string (format, token, buf, buflen, error);
      if (dbg)
	indent -= 2;
    }
}



static int
valid_format_letter_p (hb_format_t fmt)
{
  const char valid_letters[5] = "DEFGI";
  const size_t nvalid = 5;
  const size_t type = fmt.type;
  size_t i;

  for (i = 0; i < nvalid; i++)
    if (type == valid_letters[i])
      return 1;

  return 0;
}

void
print_hb_format (FILE* f, hb_format_t fmt)
{
  assert (valid_format_letter_p (fmt));

  fprintf (stderr, "(%zu%c", fmt.count, fmt.type);
  if (fmt.type == 'F' || fmt.type == 'I')
    fprintf (stderr, "%zu)", fmt.params[0]);
  else if (fmt.type == 'D' || fmt.type == 'E' || fmt.type == 'G')
    fprintf (stderr, "%zu.%zu)", fmt.params[0], fmt.params[1]);
  else
    {
      /* Shouldn't reach this point! */
      assert (0);
    }
}
