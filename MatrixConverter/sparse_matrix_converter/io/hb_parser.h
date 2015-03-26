#ifndef _hb_parser_h
#define _hb_parser_h
/**
 * @file hb_parser.h
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-16 10:58:46 mhoemmen>
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
#include <stdio.h>
#include <stdlib.h>


/**
 * Every Harwell-Boeing format string has at least a count and a field width.
 * It may have at most two additional parameters after that:  
 *
 * For fixed formats, digits after decimal point;
 * For integer formats, minimum number of digits in the number 
 *   (defaults to zero, which means no particular number (???));
 * For float formats, significand width and exponent width.
 */
typedef struct __hb_format_t {
  size_t count;
  char type;
  size_t field_width;
  size_t params[2];
} hb_format_t;


/**
 * Parse the token stream given by TOKEN into a Fortran 77 format
 * specifier of the form acceptable for Harwell-Boeing sparse matrix
 * files.  Return the remaining token stream.
 *
 * @param format [OUT]:  The parsed format specifier.
 * @param s [IN]:        The string to parse
 * @param buf [OUT]:     String buffer for parsing.
 * @param buflen [IN/OUT]: Length of BUF array.  If not long enough,
 *                      then set to the required length (see below).
 * @param error [OUT]: If FMT_NOMEM, then BUF isn't long enough, and
 *              BUFLEN is set to the required value.  Otherwise, if
 *              not FMT_NO_ERROR, then something went wrong.
 */
void
parse_hb_format (hb_format_t* format, 
		 const char s[],
		 const size_t lenmax,
		 char* buf, 
		 size_t* buflen, 
		 int* error);

void
print_hb_format (FILE* f, hb_format_t fmt);

#endif /* _hb_parser_h */
