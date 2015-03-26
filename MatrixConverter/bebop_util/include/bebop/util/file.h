#ifndef _file_h
#define _file_h
/****************************************************************************
 * @file file.h
 * @author Mark Hoemmen
 * @since 23 Nov 2007
 * @date Time-stamp: <2008-07-16 10:16:54 mhoemmen>
 *
 * Declaration of filesystem utility functions.  
 *
 * @note Moved out of util.h and into this file on 23 Nov 2007.
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
#include <stdio.h>   /* FILE* declaration */

/**
 * Separates a complete pathname for a file into a parent directory, a 
 * namestem (with the .* extension removed) and an extension.  If there is no 
 * extension, *extn==NULL.  The returned char* are copies, allocated using 
 * bebop_malloc() (and not malloc(), so you have to use bebop_free() rather 
 * than free() to free them).
 *
 * @return Nonzero if error, else zero.
 */
void
split_pathname (char** parentdir, char** namestem, char** extn, const char* const path);

/**
 * Returns nonzero if path is a directory, else returns zero. 
 */
int
directory_p (const char* const path);

/**
 * Returns nonzero if path is a regular file, else returns zero.
 */
int
regular_file_p (const char* const path);

/**
 * Returns the number of characters in the longest-length line in the
 * given file (starting at the current file stream position).  Aborts
 * with a log message if that number of characters is larger than can
 * be represented by an unsigned int.  Resets the file stream position
 * to its original position (before calling this function).
 *
 * @param file [IN] Valid file stream
 *
 * @return Length of longest line in file
 */
unsigned int
max_linelength_in_file (FILE* file);


#endif /* #ifndef _file_h */
