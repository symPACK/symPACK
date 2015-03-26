#ifndef _config_h
#define _config_h
/**
 * @file config.h
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-20 19:07:42 mhoemmen>
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

/* Define if you have the <sys/time.h> header file and gettimeofday() is in
   it. */
#define GETTIMEOFDAY_IN_SYS_TIME_H 1

/* Define if you have the <time.h> header file and gettimeofday() is in it. */
/* #undef GETTIMEOFDAY_IN_TIME_H */

/* 
 * mfh 20 Jul 2008: this automation process doesn't seem to be
 * working, so I got rid of it.  Now users have to set a variable in
 * bebop_make/options to 0 in order to turn off the requirement of C99
 * complex support.
 */
#if 0
/**
 * mfh 02 Aug 2006: automated the process of determining whether the 
 * complex data types are available, to make life easier for Cygwin users.
 */
#ifdef __STDC__
#  if (__STDC_VERSION__ >= 199901L)
#    define HAVE_C99 1
#  endif /* __STDC_VERSION__ >= 199901L */
#endif /* __STDC__ */

#ifdef HAVE_C99
/* Defined if you have the <complex.h> C99 standard header file */
#  ifdef HAVE_COMPLEX_H
   /* Defined if you have "float _Complex" and "double _Complex" data
      types (from C99) */
#    define HAVE__COMPLEX 1
#  endif
#endif /* HAVE_C99 */
#endif /* 0 */

/* Define to 1 if you have the <errno.h> header file. */
#define HAVE_ERRNO_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <libgen.h> standard POSIX header file */
#define HAVE_LIBGEN_H 1

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#define HAVE_MALLOC 1

/* Define to 1 if you have the <math.h> header file. */
#define HAVE_MATH_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if your system has a GNU libc compatible `realloc' function,
   and to 0 otherwise. */
#define HAVE_REALLOC 1

/* Define to 1 if you have the <stddef.h> header file. */
#define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if on AIX 3.
   System headers sometimes define this.
   We just want to avoid a redefinition error message.  */
#ifndef _ALL_SOURCE
/* # undef _ALL_SOURCE */
#endif

/* Define to 2 if the system does not provide POSIX.1 features except with
   this defined. */
/* #undef _POSIX_1_SOURCE */

/* Define to 1 if you need to in order for `stat' and other things to work. */
/* #undef _POSIX_SOURCE */


#endif /* _config_h */
