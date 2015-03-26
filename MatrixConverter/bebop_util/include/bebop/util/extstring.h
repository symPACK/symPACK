#ifndef _extstring_h
#define _extstring_h
/**
 * @file extstring.h
 * @author mfh
 * @since 03 May 2007
 * @date Time-stamp: <2008-07-16 10:16:43 mhoemmen>
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
#include <stdio.h> /* FILE* declaration */

/**
 * A string buffer holding len characters (not counting the '\n')
 * with space for maxlen characters (not counting the '\n').
 */
typedef struct __extstring_t {
  char* s;
  int len;
  int maxlen;
} extstring_t;

/**
 * Initialize str to hold maxlen max characters.
 */
void
extstring_init (extstring_t* str, int maxlen);

/**
 * De-initialize str (free the storage if necessary).
 */
void
extstring_deinit (extstring_t* str);

/**
 * Push a character onto the end of the string.  Doesn't guarantee
 * that the string is null-terminated; call extstring_terminate() for
 * that.
 */
void
extstring_pushchar (extstring_t* str, char c);

/**
 * "Complete" the string (if necessary) by ending it with a null
 * character.
 */
void
extstring_terminate (extstring_t* str);

/**
 * Get the string.  This automatically null-terminates it if it's not
 * already.
 */
char*
extstring_string (extstring_t* str);

/**
 * Clear current contents of the string and set the length to zero,
 * without freeing any extra buffer space that there may be.
 */
void
extstring_clear (extstring_t* str);

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
bebop_safe_getline (extstring_t* line, FILE* infile);

#endif /* _extstring_h */
