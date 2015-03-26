/**
 * @file convert.c
 * @author Mark Hoemmen
 * @since 09 May 2007
 * @date Time-stamp: <2008-07-16 10:09:25 mhoemmen>
 * 
 * Functions for converting strings to various numeric types.
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
#include <bebop/util/config.h>
#include <bebop/util/convert.h>

#include <errno.h>
#include <stdlib.h>

int
string_to_double (double* d, const char* s)
{
  extern int errno;
  double __d = 0.0;

  errno = 0;
  __d = strtod (s, NULL);
  if (errno != 0)
    {
      int retval = errno;
      errno = 0;
      return retval;
    }
  else
    {
      *d = __d;
      return 0;
    }
}

int
string_to_int (int *i, const char* s)
{
  extern int errno;
  int __i = 0;

  errno = 0;
  __i = (int) strtol (s, NULL, 10);
  if (errno != 0)
    {
      int retval = errno;
      errno = 0;
      return retval;
    }
  else
    {
      *i = __i;
      return 0;
    }
}
