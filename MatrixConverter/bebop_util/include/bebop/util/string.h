#ifndef _string_h
#define _string_h
/****************************************************************************
 * @file string.h
 * @author Mark Hoemmen
 * @since 23 Nov 2007
 * @date Time-stamp: <2008-07-16 10:20:06 mhoemmen>
 *
 * Declaration of some string utility functions (extracted from util.h).  
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
 ****************************************************************************/

/**
 * Returns a deep copy (malloc'd) of s[start : end], if this is a valid
 * substring.  If not, returns NULL.
 */
char* 
get_substring (const char* s, int start, int end);

/**
 * A classic string hash function "djb2" reported by Dan Bernstein 
 * in comp.lang.c.  See e.g. http://www.cse.yorku.ca/~oz/hash.html.
 */
unsigned long
djb2_hash (unsigned char *str);

/**
 * Returns 1 if the two strings are equal, else returns 0.
 */
int
string_equal (char* s, char* t);

/**
 * Like strdup(), but calls bebop_malloc() internally instead of malloc().
 */
char*
bebop_strdup (const char* s);

/**
 * Like strdup(), but duplicates s[first:first+len-1] and appends a
 * '\0' to the duplicate.
 */
char*
bebop_strdup_range (const char* s, const int first, const int len);

#endif /* #ifndef _string_h */
