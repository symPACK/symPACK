/**
 * @file timer.c
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-16 10:12:30 mhoemmen>
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

#ifdef GETTIMEOFDAY_IN_SYS_TIME_H
#  include <sys/time.h>
#else
#  ifdef GETTIMEOFDAY_IN_TIME_H
#    include <time.h>
#  else
#    error "Which header file are you supposed to include for gettimeofday()?"
#  endif
#endif


void
deinit_timer()
{
  /* When using gettimeofday, there is no need to de-initialize the timer. */
}


double
get_seconds()
{
  struct timeval tp;
  int i;

  /* The POSIX standard (2003 edition) says that if the 2nd argument of
   * gettimeofday is not NULL, the behavior of the function is unspecified.
   */
  i = gettimeofday (&tp,0);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}


void
init_timer()
{
  int i;
  double s = 0.0;
  double t;
  double x = 1.0;

  /*
   * MFH 2004 Jan 24: Sometimes gettimeofday doesn't work right on the
   * first call.  So I decided to put a "warmup run" here.  We time
   * some convergent operation which the compiler won't optimize away.
   */
  t = get_seconds();
  for (i = 0; i < 100000; i++)
    {
      s = s + x / 2;
      x = x / 2.0;
    }
  t = get_seconds() - t;
}











