/**
 * @file string.c
 * @author Mark Hoemmen
 * @since 23 Nov 2007
 * @date Time-stamp: <2008-07-16 10:12:22 mhoemmen>
 *
 * Implementation of string utility functions declared in string.h.  
 * Extracted from util.c on 23 Nov 2007.
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
 *************************************************************************/
#include <bebop/util/config.h>
#include <bebop/util/string.h>
#include <bebop/util/malloc.h>

#include <assert.h>
#include <string.h>


char* 
get_substring (const char* s, int start, int end)
{
  int len, i;
  char* sub = NULL;
  int substring_length = 0;

  if (s == NULL)
    return NULL;

  len = strlen (s);
  if (s == 0)
    return NULL;
  if (end >= len || start < 0 || start > end)
    return NULL;

  substring_length = end - start + 1;
  assert (substring_length >= 0);
  sub = bebop_malloc ((substring_length + 1) * sizeof (char));

  for (i = 0; i < substring_length; i++)
    sub[i] = s[i + start];

  sub[substring_length] = (char) 0;
  return sub;
}

unsigned long
djb2_hash (unsigned char *str)
{
  unsigned long hash = 5381; /* magic seed number */
  int c;

  while ((c = *str++))
    hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

  return hash;
}

int
string_equal (char* s, char* t)
{
  return (0 == strcmp (s, t));
}

char*
bebop_strdup (const char* s)
{
  char* t;
  int len;
 
  if (s == NULL)
    return NULL;
  else 
    {
      len = strlen (s);
      t = bebop_malloc ((len+1) * sizeof (char));
      return strcpy (t, s);
    }
}

char*
bebop_strdup_range (const char* s, const int first, const int len)
{
  char* scopy; 

  if (s == NULL)
    return NULL;
   
  scopy = bebop_malloc ((len + 1) * sizeof (char));
  memcpy (scopy, &s[first], len * sizeof (char));
  scopy[first+len] = '\0';
  return scopy;
}


