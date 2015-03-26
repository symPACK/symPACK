#ifndef _hb_tok_h
#define _hb_tok_h
/**
 * @file hb_tok.h
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-16 10:59:31 mhoemmen>
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

#include <stdlib.h> /* size_t */


typedef enum __hb_format_error_t {
  FMT_NO_ERROR = 0,
  FMT_NOMEM,
  FMT_TOKEN_ERROR,
  FMT_PARSE_ERROR,
  FMT_OTHER_ERROR
} hb_format_error_t;

typedef enum __hb_token_type_t {
  START,
  LPAREN,
  RPAREN,
  UINT,
  LETTER,
  PERIOD,
  END,
  NO_TOKEN
} hb_token_type_t;

typedef struct __hb_token_t {
  const char* s;
  size_t end; /* end of current token */
  size_t lenmax; /* upper bound on length */
  size_t pos; /* index of s[0] in original string input */
  hb_token_type_t type;
} hb_token_t;

/**
 * Start tokenizing a given string: convert the string into a token
 * stream.
 * 
 * @param s [IN]: The string to turn into a token stream.  Modifying
 *                this string after calling hb_token_start() puts the
 *                token stream in an undefined state.
 * @param lenmax [IN]: Upper bound on length of s.  If (size_t) -1,
 *               then there is no upper bound.
 * @param error [OUT]: If no errors, set to FMT_NO_ERROR, else set to
 *              the corresponding error value.
 *
 * @return A token stream.
 */
hb_token_t
hb_token_start (const char* s, const size_t lenmax, int* error);

/**
 * Advance the input token stream TOK by one token.
 */
hb_token_t 
next_hb_token (hb_token_t tok, int* error);

void
read_hb_token (char buf[], size_t* buflen, hb_token_t current, int* error);

void
read_hb_token_with_resize (char** buf, size_t* buflen, 
			   hb_token_t tok, int* error);

int
hb_token_error (int error);

const char*
hb_token_type_string (hb_token_t token);

char
hb_token_first_char (hb_token_t token);

int
hb_token_endp (const hb_token_t token);

size_t
hb_token_length (hb_token_t tok);


#endif /* _hb_tok_h */
