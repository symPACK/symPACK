/**
 * @file hb_tok.c
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-16 10:59:13 mhoemmen>
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

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h> 

static const int dbg = 1;
extern size_t indent; /* not thread-safe! */

static void 
do_indent ()
{
  size_t i;
  for (i = 0; i < indent; i++)
    fprintf (stderr, " ");
}


static inline int 
eos (const char c)
{
  return c == '\0';
}

static hb_token_t
skip_whitespace (hb_token_t tok);


/****************************************/

char
hb_token_first_char (hb_token_t token)
{
  return token.s[0]; /* NO SAFETY CHECKING! */
}

int
hb_token_endp (const hb_token_t token)
{
  if (dbg > 0)
    {
      do_indent();
      fprintf (stderr, "--- hb_token_endp ---\n");
    }

  if (token.type == NO_TOKEN)
    return 1;
  else if (token.type == END)
    return 1;
  else
    return 0;
}

static hb_token_t
uint_token (hb_token_t tok)
{
  const size_t lenmax = tok.lenmax;
  size_t end = -1;
  char c;
  
  if (tok.pos >= lenmax)
    {
      tok.type = END;
      return tok;
    }
  /* 
   * We now know that we can read the first character of tok.s.
   * 
   * We cut the parsing short once we reach lenmax.  Whatever is in s
   * after that isn't our business, so it's not an error to stop
   * there.
   */
  do {
    end++; /* will overflow to zero on first iteration */
    c = tok.s[end];
  } while (! eos(c) && isdigit(c) && end + tok.pos < lenmax);

  tok.end = end;
  if (end == 0)
    tok.type = NO_TOKEN;
  else
    tok.type = UINT;

  return tok;
}

hb_token_t 
next_hb_token (hb_token_t tok, int* error)
{
  char c = 'Z';

  if (dbg > 0)
    {
      do_indent ();
      fprintf (stderr, "--- next_hb_token ---\n");
    }

  if (tok.type == END)
    return tok;
  else if (tok.pos >= tok.lenmax)
    return tok;

  if (dbg > 0)
    {
      do_indent ();
      fprintf (stderr, "Past early exit\n");
      do_indent ();
      fprintf (stderr, "tok.pos = %zu, tok.end = %zu\n", tok.pos, tok.end);
    }
  
  /* 
   * Whether or not tok is in the START state, we need to advance s up
   * to the former location of tok.end.  Check carefully to make sure
   * that we haven't run over the max token stream length (lenmax).
   * 
   * Note that tok.pos < tok.lenmax (as we've checked above).
   */
  if (tok.pos + tok.end >= tok.lenmax)
    {
      if (dbg > 0)
	{
	  do_indent ();
	  fprintf (stderr, "next_hb_token: tok.pos %zu + tok.end %zu "
		   ">= tok.lenmax %zu, reached end of token stream\n", 
		   tok.pos, tok.end, tok.lenmax);
	}
      /* Restrict tok.end not to run over max token stream length */
      /* FIXME ??? */
      assert (0);

      tok.end = (tok.end > tok.lenmax) ? tok.lenmax : tok.end;
      /* Don't set tok to the END state yet; we'll check that after we
	 advance past whitespace */
    }
  else
    {
      if (dbg)
	{
	  do_indent ();
	  fprintf (stderr, "Advancing token stream\n");
	}
      tok.pos += tok.end;
      tok.s = &tok.s[tok.pos];
      tok.end = 0;
    }

  *error = FMT_NO_ERROR;
  if (dbg)
    indent += 2;
  tok = skip_whitespace (tok);
  if (dbg)
    indent -= 2;
  if (tok.type == NO_TOKEN)
    {
      if (dbg)
	{
	  do_indent();
	  fprintf (stderr, "skip_whitespace returned a NO_TOKEN\n");
	}
      *error = FMT_TOKEN_ERROR;
      return tok;
    }
  else if (tok.type == END)
    {
      if (dbg)
	{
	  do_indent();
	  fprintf (stderr, "skip_whitespace returned an END\n");
	}
      return tok;
    }
  /* We can have tok.end == 0 after skip_whitespace(), if there's no
     whitespace at the beginning of the stream */

  if (dbg > 0)
    {
      do_indent ();
      fprintf (stderr, "We can read first char\n");
    }

  c = hb_token_first_char (tok);
  if (dbg > 0)
    fprintf (stderr, "next_hb_token: got first char %c\n", c);

  if (c == '(')
    {
      tok.end = 1;
      tok.type = LPAREN;
    }
  else if (c == ')')
    {
      tok.end = 1;
      tok.type = RPAREN;
    }
  else if (c == '.')
    {
      tok.end = 1;
      tok.type = PERIOD;
    }
  else if (isalpha(c))
    {
      tok.end = 1;
      tok.type = LETTER;
    }
  else if (isdigit(c))
    {
      if (dbg > 0)
	indent += 2;
      tok = uint_token(tok);
      if (dbg > 0)
	indent -= 2;
    }
  else
    {
      if (dbg > 0)
	{
	  do_indent ();
	  fprintf (stderr, "bad first character \'%c\'\n", c);
	}
      tok.end = 0;
      tok.type = NO_TOKEN;
      *error = FMT_TOKEN_ERROR;
    }

  if (dbg > 0)
    {
      do_indent ();
      fprintf (stderr, "Got token of type %s\n", hb_token_type_string(tok));
    }
  return tok;
}

static hb_token_t
skip_whitespace (hb_token_t tok)
{
  if (dbg > 0)
    {
      do_indent ();
      fprintf (stderr, "--- skip_whitespace ---\n");
      do_indent ();
      fprintf (stderr, "tok.pos = %zu\n", tok.pos);
    }

  if (tok.type == NO_TOKEN)
    return tok;
  else if (tok.s == NULL)
    {
      tok.type = NO_TOKEN;
      return tok;
    }
  else if (tok.pos >= tok.lenmax)
    {
      tok.type = END;
      return tok;
    }
  else
    {
      const size_t lenmax = tok.lenmax;
      size_t start = (size_t) -1; /* will overflow on first iteration */
      char c;

      if (dbg)
	{
	  do_indent ();
	  fprintf (stderr, "Looking for end of whitespace\n");
	}

      /* We know that tok.pos < tok.lenmax, so we can read tok.s[0] */
      do {
	start++; /* will overflow to zero on first iteration */
	c = tok.s[start];
      } while (start < lenmax && ! eos(c) && isspace(c));

      if (start >= lenmax || eos(tok.s[start]))
	{
	  /* Whitespace isn't a token, so if we find end-of-string
	     after whitespace, we return "no token" */
	  if (dbg)
	    {
	      do_indent ();
	      fprintf (stderr, "There\'s nothing but whitespace: start == %zu\n", start);
	    }
	  tok.end = 0;
	  tok.type = NO_TOKEN;
	}
      else
	{
	  tok.end = start;
	  tok.pos += start;
	  tok.s = &tok.s[tok.pos];
	  if (dbg)
	    {
	      do_indent ();
	      fprintf (stderr, "Got to start of non-whitespace: first character is %c\n", c);
	    }
	}

      return tok;
    }
}

size_t
hb_token_length (hb_token_t tok)
{
  return tok.end;
}

void
read_hb_token_with_resize (char** buf, size_t* buflen, 
			   hb_token_t tok, int* error)
{
  *error = 0;
  read_hb_token (*buf, buflen, tok, error);
  if (*error == FMT_NOMEM)
    {
      /* We've already accounted for the '\0' in read_hb_token() */
      *buf = realloc (*buf, (*buflen) * sizeof(char));
      assert (*buf != NULL);
      *error = 0;
      read_hb_token (*buf, buflen, tok, error);
    }
}

void
read_hb_token (char buf[], size_t* buflen, hb_token_t current, int* error)
{
  size_t toklen;

  if (dbg > 0)
    {
      size_t i;
      for (i = 0; i < indent; i++)
	fprintf (stderr, " ");
      fprintf (stderr, "--- read_hb_token ---\n");
      indent += 2;
    }
  toklen = hb_token_length (current);
  if (dbg > 0)
    indent -= 2;
    
  if (*buflen < toklen + 1) /* Need room for '\0' */
    {
      *error = FMT_NOMEM;
      *buflen = toklen + 1;
    }
  else
    {
      size_t i;
      for (i = 0; i < toklen; i++)
	buf[i] = current.s[i];
      buf[toklen] = '\0';
      *error = FMT_NO_ERROR;
    }
}

hb_token_t
hb_token_start (const char* s, const size_t lenmax, int* error)
{
  hb_token_t token;

  if (dbg > 0)
    {
      do_indent();
      fprintf (stderr, "--- hb_token_start ---\n");
      do_indent();
      fprintf (stderr, "s == %s, lenmax = %zu\n", s, lenmax);
    }

  *error = 0;
  token.s = s;
  token.end = 0; /* to be changed later */
  token.pos = 0;
  token.lenmax = lenmax;
  if (lenmax == 0 || s == NULL || s[0] == '\0')
    {
      token.type = END;
      return token;
    }
  else
    {
      token.type = START;
      if (dbg)
	indent += 2;
      token = next_hb_token (token, error);
      if (dbg)
	indent -= 2;
      return token;
    }
}


const char*
hb_token_type_string (hb_token_t token)
{
  switch (token.type) 
    {
    case NO_TOKEN:
      return "NO_TOKEN";
    case START:
      return "START";
    case LPAREN:
      return "LPAREN";
    case UINT:
      return "UINT";
    case LETTER:
      return "LETTER";
    case PERIOD:
      return "PERIOD";
    case RPAREN:
      return "RPAREN";
    case END:
      return "END";
      /*
       * Should never reach this point, since there aren't any
       * other kinds of tokens.
       */
      assert (0); 
    }
  return ""; /* to pacify the compiler */
}

int
hb_token_error (int error)
{
  return error != FMT_NO_ERROR;
}

