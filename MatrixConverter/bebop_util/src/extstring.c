/**
 * @file extstring.c
 * @author Mark Hoemmen
 * @since 03 May 2007
 * @date Time-stamp: <2008-07-16 10:09:36 mhoemmen>
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
#include <string.h>

#include <bebop/util/extstring.h>
#include <bebop/util/init.h>
#include <bebop/util/log.h>
#include <bebop/util/malloc.h>
#include <bebop/util/util.h>

static inline int
extstring_valid_p (extstring_t* str)
{
  int validp = 1;

  if (str->s == NULL)
    {
      bebop_log (0, "*** extstring_valid_p: str->s == NULL ***\n");
      validp = 0;
    }
  else if (str->len < 0)
    {
      bebop_log (0, "*** extstring_valid_p: str->len = %d < 0 ***\n", 
		str->len);
      validp = 0;
    }
  else if (str->maxlen < str->len)
    {
      bebop_log (0, "*** extstring_valid_p: str->maxlen = %d "
		"< str->len = %d ***\n", str->maxlen, str->len);
      validp = 0;
    }
  return validp;
}

void
extstring_init (extstring_t* str, int maxlen)
{
  /* bebop_log (2, "=== extstring_init ===\n"); */
  if (maxlen < 0)
    {
      bebop_log (0, "*** extstring_init: maxlen = %d < 0; I silently "
		"corrected maxlen to zero ***\n", maxlen);
      maxlen = 0;
    }
  /* bebop_calloc() gives us automatic '\0'-termination */
  str->s = bebop_calloc ((maxlen + 1), sizeof (char));
  str->len = 0;
  str->maxlen = maxlen;

  assert (extstring_valid_p (str));
  /* bebop_log (2, "=== Done with extstring_init ===\n"); */
}

/**
 * Returns numeric true (nonzero) if the string needs resizing in
 * order to accomodate s[spot], otherwise returns zero.
 */
static inline int
extstring_need_resize_p (const extstring_t* str, const int spot)
{
  return (spot + 1 >= str->maxlen);
}

/**
 * Resize the string by doubling its max length.  
 */
static inline void
extstring_expand (extstring_t* str)
{
  /* bebop_log (2, "=== extstring_expand ===\n"); */
  assert (extstring_valid_p (str));

  if (str->maxlen == 0)
    str->maxlen = 4; /* better than 1 */
  else
    str->maxlen = 2 * str->maxlen;
  str->s = bebop_realloc (str->s, (str->maxlen+1) * sizeof (char));

  assert (extstring_valid_p (str));
  /* bebop_log (2, "=== Done with extstring_expand ===\n"); */
}

/**
 * Push c onto end of str, without len++.
 */
static void
extstring_pushchar_nolen (extstring_t* str, char c)
{
  int spot;

  assert (extstring_valid_p (str));
  spot = str->len;

  if (extstring_need_resize_p (str, spot))
    extstring_expand (str);

  str->s[spot] = c;
  assert (extstring_valid_p (str));
}

void
extstring_pushchar (extstring_t* str, char c)
{
  extstring_pushchar_nolen (str, c);
  str->len++;
  assert (extstring_valid_p (str));
}

void
extstring_terminate (extstring_t* str)
{
  /* bebop_log (2, "=== extstring_terminate ===\n"); */
  assert (extstring_valid_p (str));

  if (str->s[str->len] != '\0')
    /* Don't do len++, or the string's length will include the '\0' */
    extstring_pushchar_nolen (str, '\0');

  assert (extstring_valid_p (str));
  /* bebop_log (2, "=== Done with extstring_terminate ===\n"); */
}

char*
extstring_string (extstring_t* str)
{
  char* retval; 

  /* bebop_log (2, "=== extstring_string ===\n"); */
  assert (extstring_valid_p (str));
  extstring_terminate (str);
  assert (extstring_valid_p (str));
  retval = str->s;
  /* bebop_log (2, "=== Done with extstring_string ===\n"); */
  return retval;
}

void
extstring_deinit (extstring_t* str)
{
  /* bebop_log (2, "=== extstring_deinit ===\n"); */
  assert (extstring_valid_p (str));

  if (str->s != NULL)
    {
      bebop_free (str->s);
      str->s = NULL;
    }

  str->len = 0;
  str->maxlen = 0;

  /* str is no longer valid, since str->s is NULL */
  /* bebop_log (2, "=== Done with extstring_deinit ===\n"); */
}

void
extstring_clear (extstring_t* str)
{
  /* bebop_log (2, "=== extstring_clear ===\n"); */
  assert (extstring_valid_p (str));

  str->len = 0;
  str->s[0] = '\0';
  /* maxlen stays where it was, since we don't deallocate storage */

  assert (extstring_valid_p (str));
  /* bebop_log (2, "=== Done with extstring_clear ===\n"); */
}

/** 
 * A safe but potentially slow version of gets().  Reads characters
 * one by one, pushing them onto the output string stored in "line",
 * until it encounters EOF or an end-of-line character sequence (which
 * is OS-dependent).  Neither the EOF character or the end-of-line
 * sequence is stored in "line".
 *
 * @param line [OUT]  Where to write the line read from the file.
 *                    Assume initialized.
 * @param linelen [OUT]  Number of characters written to line, not
 *                       counting the '\0' (note that any '\n' are
 *                       not written)
 * @param linebuflen [IN/OUT]  Number of characters stored in line.
 *
 * @return Number of characters read, not counting EOF or any
 * end-of-line characters.
 */
int
bebop_safe_getline (extstring_t* line, FILE* infile)
{
  extstring_clear (line);

  if (feof (infile))
    return 0;
  else
    {
      while (! feof (infile))
	{
	  char c = fgetc (infile);
	  if (c == '\n')
	    break; /* Unix(tm) - style endline */
	  else if (c == EOF)
	    break; /* end of file, so stop reading */
	  else if (c == '\r')
	    {
	      char next;

	      /* Macs have "\r" as the endline character, whereas DOS
		 (and Windows) uses "\r\n".  We have to check the next
		 character to see if it's an "\n".  Luckily we're
		 allowed at least one pushback onto the stream, so we
		 can put the character back if it's not an '\n'.  We
		 should first, of course, check for end-of-file. */

	      if (feof (infile))
		break;

	      next = fgetc (infile);
	      if (next == '\n')
		break; /* It's a DOS/Windows endline */

	      /* It's a Mac endline; put back the character.  A single
		 putback isn't supposed to fail, but we should check
		 anyway. */
	      if (EOF == ungetc (next, infile))
	        bebop_fatal_error ("IO", "safe_getline: ungetc() failed (not supposed to happen)");
	    }
	  extstring_pushchar (line, c);
	}
      return line->len;
    }
}


