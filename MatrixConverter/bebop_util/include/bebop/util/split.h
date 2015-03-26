#ifndef _split_h
#define _split_h
/**
 * @file split.h
 * @author mfh
 * @since 02 May 2007
 * @date Time-stamp: <2008-07-16 10:19:55 mhoemmen>
 * 
 * split_t and member functions: split a string into tokens based on
 * delimiters.  The tokens are deep copies of the relevant portions of
 * the original string (thus avoiding the thread-(un)safety issues of
 * strtok()).
 *
 * manysplit_t and member functions: optimization of split_t for
 * running splits on many strings in turn.  It keeps an internal
 * buffer for holding tokens, rather than allocating deep copies for
 * each tokenized string.  The tokens are still deep copies of the
 * original string, though.  The manysplit_t interface is different
 * because a manysplit_t object is conceptually "weightier" (and
 * therefore better handled as a pointer rather than a shallow-copied
 * struct).
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
#include <bebop/util/list.h>

/**
 * Struct representing the results of a split() operation on a string.
 */
typedef struct __split_t {
  list_t tokens;
  int ntokens;
} split_t;

split_t 
new_split ();

list_t
split_tokens (split_t s);

int
split_ntokens (split_t s);

split_t
split (const char* tosplit, 
       const char delimiters[], 
       const int ndelimiters);

void
split_destroy (split_t spl);


/** 
 * Object that stores the results of a manysplit operation.
 */
typedef struct __manysplit_t {
  list_t tokens;
  int ntokens;
  char* buf;
  int buflen;
  int curpos;
} manysplit_t;

/**
 * Initializes the given manysplit_t object, setting up the internal
 * string buffer to hold a line containing linelen_estimate
 * characters.  The string buffer will be expanded if necessary.
 */
void
manysplit_init (manysplit_t* s, int linelen_estimate);

/**
 * "Resets" a previously initialized manysplit_t object s -- releases
 * the stored tokens and reinitializes the internal string buffer to
 * hold linelen_estimate characters.
 */
void
manysplit_reset (manysplit_t* s, int linelen_estimate);

/**
 * Returns the current token list stored in s.  See the manysplit()
 * function for details.
 */
list_t
manysplit_tokens (manysplit_t* s);

/**
 * Destructive version of manysplit_tokens(); conses only if
 * necessary.
 */
void
manysplit_dtokens (list_t* out, manysplit_t* s);

/**
 * Returns the current number of tokens stored in the token list in s.
 */
int
manysplit_ntokens (manysplit_t* s);

/**
 * Given a pointer to an initialized manysplit_t buffer, splits the
 * string "tosplit" using the given delimiters, and stores the
 * resulting tokens in a list in manysplit.  The list can be extracted
 * by calling manysplit_tokens().
 */
void
manysplit (manysplit_t* manysplit,
	   const char* tosplit, 
	   const char delimiters[], 
	   const int ndelimiters);

/**
 * De-initializes the contents of spl -- frees the list of tokens and
 * the string buffer.  Does not free the spl pointer itself.
 */
void
manysplit_deinit (manysplit_t* spl);

#endif /* _split_h */
