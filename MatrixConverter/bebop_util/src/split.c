/** 
 * @file split.c
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-16 10:12:12 mhoemmen>
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
#include <bebop/util/split.h>
#include <bebop/util/malloc.h>
#include <bebop/util/string.h>
#include <bebop/util/util.h>
#include <bebop/util/log.h>

#include <assert.h>
#include <string.h>

#if 0
static void
print_string_list (list_t L)
{
  list_t L2;
  for (L2 = L;
       ! list_empty_p (L2);
       L2 = list_cdr (L2))
    fprintf (stderr, "%s ", (char*) list_car (L2));
  fprintf (stderr, "\n");
}

static void
print_buf_contents (manysplit_t* ms)
{
  int i;
  for (i = 0; i < ms->buflen; i++)
    {
      char c = ms->buf[i];
      if (c != '\0')
	fprintf (stderr, "%c", ms->buf[i]);
      else
	fprintf (stderr, "-");
    }
  fprintf (stderr, "\n");
}
#endif /* 0 */


split_t 
new_split ()
{
  split_t s;
  s.tokens = list_create ();
  s.ntokens = 0;
  return s;
}

list_t
split_tokens (split_t s)
{
  return s.tokens;
}

int
split_ntokens (split_t s)
{
  return s.ntokens;
}

static int
char_member (const char c, const char elts[], const int nelts)
{
  int i;

  for (i = 0; i < nelts; i++)
    {
      if (c == elts[i])
	return 1;
    }
  return 0;
}

split_t
split (const char* tosplit, const char delimiters[], const int ndelimiters)
{
  split_t spl = new_split ();
  int i = 0;

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
}

void
split_destroy (split_t spl)
{
  /* FIXME: this only works because bebop_free() (which is a macro,
     otherwise we would use it here with list_destroy_custom()) calls
     free(). */
  spl.tokens = list_destroy (spl.tokens);
  spl.ntokens = 0;
}


#if 0
/**
 * A noop; destructor for list_t elements in the manysplit case.
 */
static void
nodel (void* v)
{
  return;
}
#endif /* 0 */

/**
 * Can alter this in order to tune the buffer recycling policy.
 */
static int
manysplit_recycle_buffer_p (const int linelen_estimate,
			    const int buflen)
{
  return linelen_estimate == buflen ||
    linelen_estimate == buflen + 1 ||
    linelen_estimate == buflen - 1;
}

void
manysplit_reset (manysplit_t* s, int linelen_estimate)
{
  bebop_log (2, "=== manysplit_reset ===\n");
  s->ntokens = 0;
  s->curpos = 0;

  if (linelen_estimate < 0)
    {
      bebop_log (0, "new_manysplit: linelen_estimate = %d "
		"which is < 0 !!!\n", linelen_estimate);
      /* Silently set it to zero */
      linelen_estimate = 0;
    }

  if (manysplit_recycle_buffer_p (linelen_estimate, s->buflen))
    {
      bebop_log (2, "recycling buffer\n");
    }
  else
    {
      bebop_log (2, "resizing buffer to %d\n", linelen_estimate + 1);
      s->buf = bebop_realloc (s->buf, (linelen_estimate + 1) * sizeof (char));
      s->buf[linelen_estimate] = '\0';
    }
  s->buflen = linelen_estimate + 1;
  bebop_log (2, "=== Done with manysplit_reset ===\n");
  return;
}

void
manysplit_init (manysplit_t* s, int linelen_estimate)
{
  bebop_log (2, "=== manysplit_init ===\n");
  s->ntokens = 0;
  s->curpos = 0;

  if (linelen_estimate < 0)
    {
      bebop_log (0, "new_manysplit: linelen_estimate = %d "
		"which is < 0 !!!\n", linelen_estimate);
      /* Silently set it to zero */
      linelen_estimate = 0;
    }
  
  /* fill with zeros (important!) */
  s->buf = bebop_calloc ((linelen_estimate + 1), sizeof (char));
  /* s->buf[linelen_estimate] = '\0'; */ /* <- don't need this */
  s->buflen = linelen_estimate + 1;

  bebop_log (2, "=== Done with manysplit_init ===\n");
  return;
}

list_t
manysplit_tokens (manysplit_t* s)
{ 
  /* not the buffer length, but the number of non-junk chars */
  const int nbufelts = s->curpos; 
  char* buf = s->buf;
  list_t tokens = list_create ();
  int token_count = 0;
  int start = 0;
  int i = 0;  

  while (i < nbufelts)
    {
      while (buf[i] == '\0' && start < nbufelts)
	i++;
      if (i >= nbufelts)
	break; /* ran out of non-junk data */
      start = i; /* Now buf[start] is a non-junk datum */

      /* Find the end of this non-junk chunk of data */
      while (buf[i] != '\0' && i < nbufelts)
	i++;
      /* the string MUST be null-terminated, so that we can append
	 just the pointer onto the list of tokens, rather than doing a
	 deep copy. */
      assert (buf[i] == '\0'); 
      tokens = list_append_item (tokens, &buf[start]);
      token_count++;
    }
  assert (token_count == s->ntokens);
  return tokens;
}

/**
 * Destructive version of manysplit_tokens(); conses only if
 * necessary.
 */
void
manysplit_dtokens (list_t* out, manysplit_t* s)
{
  list_t outlist = *out;
  list_t current = *out; /* current spot in outlist to overwrite */

  /* not the buffer length, but the number of non-junk chars */
  const int nbufelts = s->curpos; 
  char* buf = s->buf;
  int token_count = 0;
  int start = 0;
  int i = 0;  

  while (i < nbufelts)
    {
      while (buf[i] == '\0' && start < nbufelts)
	i++;
      if (i >= nbufelts)
	break; /* ran out of non-junk data */
      start = i; /* Now buf[start] is a non-junk datum */

      /* Find the end of this non-junk chunk of data */
      while (buf[i] != '\0' && i < nbufelts)
	i++;
      /* the string MUST be null-terminated, so that we can append
	 just the pointer onto the list of tokens, rather than doing a
	 deep copy. */
      assert (buf[i] == '\0'); 
      if (list_empty_p (current))
	{
	  /* We're out of nodes to overwrite, so we have to cons */
	  outlist = list_append_item (outlist, (void*) (&buf[start]));
	}
      else /* overwrite the current node, and move to the next one */
	{
	  list_set_car (current, &buf[start]);
	  current = list_cdr (current);
	}
      token_count++;
    }
  assert (token_count == s->ntokens);
  *out = outlist; /* necessary so that tail pointer is set correctly */
  return;
}


int
manysplit_ntokens (manysplit_t* s)
{
  return s->ntokens;
}

static char*
manysplit_new_tokenspot (manysplit_t* manysplit, const int toklen)
{
  char* buf = manysplit->buf;
  int buflen = manysplit->buflen;
  int curpos = manysplit->curpos;
  char* retval = NULL;

  bebop_log (2, "=== manysplit_new_tokenspot ===\n");
  bebop_log (3, "old curpos: %d\n", curpos);
  assert (curpos >= 0);
  assert (buflen >= 0);
  assert (buf != NULL);

  /* Expand buffer if necessary */
  if (curpos + toklen >= buflen)
    {
      buflen = MAX( 2 * buflen, curpos + toklen + 1 );
      buf = bebop_realloc (buf, buflen * sizeof (char));
    }

  retval = &(buf[curpos]);
  manysplit->curpos += toklen + 1; /* save a spot for the '\0' */
  manysplit->buf = buf;
  manysplit->buflen = buflen;

  bebop_log (3, "old curpos: %d\n", manysplit->curpos);
  bebop_log (2, "=== Done with manysplit_new_tokenspot ===\n");
  return retval;
}

static char*
manysplit_strdup (manysplit_t *ms, 
		  const char* token,
		  const int toklen)
{
  char* tokspot = NULL;
  bebop_log (2, "=== manysplit_strdup ===\n");

  if (token == NULL)
    {
      bebop_log (1, "manysplit_strdup: token is NULL\n");
    }
  else
    {
      int i;

      assert (toklen >= 0);
      tokspot = manysplit_new_tokenspot (ms, toklen);
      assert (tokspot != NULL);
      for (i = 0; i < toklen; i++)
	{
	  tokspot[i] = token[i];
	  fprintf (stderr, "i = %d\n", i);
	}
      /* memcpy (tokspot, token, toklen * sizeof (char)); */
      tokspot[toklen] = '\0';
    }

  bebop_log (2, "=== Done with manysplit_strdup\n");
  return tokspot;
}

void
manysplit (manysplit_t* ms,
	   const char* tosplit, 
	   const char delimiters[], 
	   const int ndelimiters)
{
  int i = 0, start = 0;

  bebop_log (2, "=== manysplit ===\n");

  /* Using buflen as the length estimate just recycles existing
     storage.  The buffer inside manysplit will be resized as 
     necessary. */
  manysplit_reset (ms, ms->buflen);

  if (tosplit == NULL)
    {
      bebop_log (1, "manysplit:  tosplit == NULL\n");
      bebop_log (2, "=== Done with manysplit ===\n");
      return;
    }

  while (tosplit[i] != '\0')
    {
      int toklen = 0; /* length of current token */
      char* newtoken = NULL;

      /* Get past the current range of delimiters */
      while (tosplit[i] != '\0' && 
	     tosplit[i] != EOF && 
	     char_member (tosplit[i], delimiters, ndelimiters))
	i++;
      if (tosplit[i] == '\0' || tosplit[i] == EOF)
	{
	  bebop_log (3, "manysplit: tosplit[%d] = null or EOF char\n", i);
	  bebop_log (2, "=== Done with manysplit ===\n");
	  return;
	}
      start = i;
      assert (! char_member (tosplit[start], delimiters, ndelimiters));
      /* Now tosplit[start] is not (a delimiter, EOF, or end-of-line) */

      /* Get to the end of the non-delimiter characters */
      while (tosplit[i] != '\0' 
	     && tosplit[i] != EOF &&
	     ! char_member (tosplit[i], delimiters, ndelimiters))
	i++;
      assert (i >= start);
      toklen = i - start;  /* don't include the last char */
#if 0
      if (bebop_debug_level() > 1)
	{
	  int j;
	  fprintf (stderr, "The token:\n");
	  for (j = 0; j < toklen; j++)
	    fprintf (stderr, "%c", tosplit[start + j]);
	  fprintf (stderr, "\n");
	}
#endif /* 0 */
      newtoken = manysplit_strdup (ms, &tosplit[start], toklen);
      
      if (bebop_debug_level() > 1)
	{
	  int success = 1;
	  int i;
	  for (i = 0; i < toklen; i++)
	    if (newtoken[i] == EOF || 
		newtoken[i] == '\0' || 
		char_member (newtoken[i], delimiters, ndelimiters))
	      {
		success = 0;
		bebop_log (1, "manysplit: newtoken[%d] is a delimiter"
			  " or other forbidden character\n", i);
	      }
	  assert (success);
	}

      bebop_log (2, "manysplit: new token: %s\n", newtoken);
#if 0
      if (bebop_debug_level() > 1)
	{
	  fprintf (stderr, "buffer contents: ");
	  print_buf_contents (ms);
	}
#endif /* 0 */
      ms->ntokens++;
    }
  bebop_log (2, "=== Done with manysplit ===\n");
  return;
}


void
manysplit_deinit (manysplit_t* spl)
{
  bebop_log (2, "=== manysplit_deinit ===\n");
  spl->ntokens = 0;
  bebop_free (spl->buf);
  spl->buf = NULL;
  spl->buflen = 0;
  spl->curpos = 0;
  bebop_log (2, "=== Done with manysplit_deinit ===\n");
}

