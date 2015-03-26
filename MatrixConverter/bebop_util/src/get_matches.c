/** 
 * @file get_matches.c
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-16 10:09:52 mhoemmen>
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
#include <sys/types.h>
#include <regex.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bebop/util/list.h>


static list_node_t*
get_match (const regex_t* regex, const char* string, int eflags, int* new_offset)
{
  const size_t nmatch = 1;
  regmatch_t pmatch[1];
  int errcode = 0;
  int len = 0;
  char* match = NULL;
  int j;

  errcode = regexec (regex, string, nmatch, pmatch, eflags);
  if (errcode == REG_NOMATCH)
    {
      printf ("get_match: no match\n");
      return NULL;
    }

  assert (pmatch[0].rm_so != -1);
  len = pmatch[0].rm_eo - pmatch[0].rm_so;
  assert (len > 0);
  printf ("Match %d,%d, length:  %d\n", (int) (pmatch[0].rm_so), (int) (pmatch[0].rm_eo), len);
  match = malloc ((len+1) * sizeof (char));
  assert (match != NULL);

  for (j = 0; j < len; j++)
    match[j] = string[j + pmatch[0].rm_so];

  match[len] = (char) 0;
  printf ("Match: %s\n", match);

  *new_offset = pmatch[0].rm_eo;
  return list_node_create ((void*) match, NULL);
}


list_t
get_matches (const regex_t* regex, char* line, int eflags)
{
  list_t L = list_create ();
  list_node_t* node = NULL;
  int len = strlen (line);
  int new_offset = 0;
  int num_matches = 0;
  int max_num_matches = 10;
  char* string = line;

  printf ("Original string: %s\n", line);

  if (len == 0)
    return L;

  node = get_match (regex, line, eflags, &new_offset);
  while (node != NULL && num_matches < max_num_matches)
    {
      /* Got a match */
      num_matches++;
      L = list_append_node (L, node);
      if (new_offset == len)
	break;
      string = &string[new_offset];
      node = get_match (regex, string, eflags, &new_offset);
    }
  
  return L;
}


