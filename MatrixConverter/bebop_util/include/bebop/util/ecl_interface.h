#ifndef _ecl_interface_h
#define _ecl_interface_h
/**
 * @file ecl_interface.h
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-16 10:16:08 mhoemmen>
 * 
 * Interface to ECL (Embedded Common Lisp) functionality
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

#include <stdint.h>
#include <stdio.h>
#include <ecl/ecl.h>


/**
 * Return nonzero iff obj represents an error from the Lisp world.
 * Does not evaluate obj, just compares it to a known object
 * representing an error condition.  If errmsg != NULL, set it to
 * point to a string (which you shouldn't modify!) indicating the
 * error condition.
 */
int
lisp_error_p (cl_object obj, 
	      const char** errmsg);

/**
 * Start up the embedded Lisp system
 * 
 * @param argc [IN]  Same as in C's main() function
 * @param argv [IN]  Same as in C's main() function
 * @param info [OUT] Zero if successful, else nonzero (in which case, 
 *                   assume that startup failed and don't call any 
 *                   Lisp functions)
 */
void
lisp_boot (int argc, char** argv, int* const info);

/**
 * Return nonzero if x is non-NIL, else return zero.
 */
int
lisp_true (cl_object x);

/**
 * Return nonzero if x is NIL, else return zero.
 */
int
lisp_false (cl_object x);

/**
 * Convert a Lisp integer to a C int32_t, if possible.
 *
 * @param x [IN]     Lisp integer to convert to a C integer
 * @param info [OUT] Zero if no errors, positive on overflow (the Lisp
 *                   integer doesn't fit in an int32_t), negative on
 *                   conversion error
 * @return If no errors, a C int32_t which is equal to the given Lisp 
 *         integer x.  If *info != 0 on exit, then the return value is
 *         undefined.
 */
int32_t
lisp_object_to_int32 (cl_object x, int* const info);

/**
 * Return a C char* corresponding to the given Lisp string.  If the
 * conversion failed, set *info to nonzero, else set *info to zero.
 */
const char* const
lisp_string_to_C_string (cl_object obj, int* const info);

void
lisp_open_file (cl_object *stream,
		const char* const filename,
		int* const info);

void
lisp_close_file (cl_object stream,
		 int* const info);

/**
 * Wrap an open C FILE stream in a Lisp stream. 
 */
cl_object 
make_lisp_stream_from_C_FILE (FILE* f);


#endif /* _ecl_interface_h */
