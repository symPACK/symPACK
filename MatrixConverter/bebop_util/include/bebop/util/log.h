#ifndef _log_h
#define _log_h
/**
 * @file log.h
 * @author Mark Hoemmen
 * @since 25 Jul 2006
 * @date Time-stamp: <2008-07-16 10:18:11 mhoemmen>
 *
 * Simple implementation of logging.  If USE_PTHREADS is defined, 
 * uses POSIX pthreads mutexes to protect accesses to global data.
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
#include <stdarg.h> /* lets us do printf-style variable # of args */

/**
 * Predicate for telling if logging is currently active.
 *
 * @warn If USE_PTHREADS is #defined and any of the mutex functions fail, 
 * calls bebop_exit().
 *
 * @returns Nonzero if logging is currently active, else returns zero.
 */
int
bebop_logging_p ();

/** 
 * Starts logging.
 *
 * If use_stream is nonzero, second argument is a FILE* to an open valid file
 * stream.  You are responsible for closing it if necessary after
 * bebop_stop_logging() is called.  If use_stream is zero, second argument
 * (char*) is a path to a log file; if third argument (int) is nonzero, we
 * append, otherwise we truncate.  We open and close the log file for you.
 *
 * @note If USE_PTHREADS is #defined, we use a "once only" technique to 
 * ensure that this function only takes effect once for all the threads,
 * so that we e.g. don't try to open the logfile more than once.
 *
 * @warn If USE_PTHREADS is #defined and any of the mutex or once only 
 * functions fail, calls bebop_exit().
 * 
 * @return zero if successful, nonzero if not.
 */
int 
bebop_start_logging (const int use_stream, ...);

/**
 * Stops logging.  Closes the logfile if necessary.
 *
 * @warn If USE_PTHREADS is #defined and any of the mutex functions fail, 
 * calls bebop_exit().
 */
void
bebop_stop_logging ();

/**
 * Stops logging, without calling bebop_exit() on failure of any of the mutex
 * functions.  This allows us to call bebop_emergency_stop_logging() in
 * bebop_exit() itself (otherwise there would be an infinite recursion).
 */
void
bebop_emergency_stop_logging ();

/**
 * Does the logging:  if the current debug level is >= min_debug_level, 
 * logs to the logfile (arguments including and after fmt work just like 
 * printf).
 *
 * @warn If USE_PTHREADS is #defined and any of the mutex functions fail, 
 * calls bebop_exit().
 */
void
bebop_log (int min_debug_level, const char* fmt, ...);

/**
 * Logs a warning of class "class".
 */
void
bebop_warn (const char* class, const char* fmt, ...);

/**
 * Logs an error of class "class".
 */
void
bebop_error (const char* class, const char* fmt, ...);

/**
 * Equivalent to "bebop_log (min_debug_level, "%s", s)".
 * Thus, very useful for calling from other programming languages, 
 * e.g. Common Lisp via an FFI, as other programming languages may
 * not support calling C functions with a variable number of args.
 */
void
bebop_log_string (int min_debug_level, const char* s);

/**
 * Macro for fatal errors.  It reports the error, as well as the 
 * function, line number, and filename in which the error occurred, 
 * and then invokes bebop_exit() to abort.
 * 
 * @param __errclass_string   The error class, as a string
 * @param __VA_ARGS__         format string for the error display,
 *                            and any further args (like printf)
 *
 * @note __VA_ARGS__ is part of the C99 standard.  Here, we avoid 
 * extensions particular to gcc (actually GNU CPP) that let you do
 * a very basic sort of pattern-matching (here, extracting the head
 * of the list) on the variable arguments.
 */
#define bebop_fatal_error( __errclass_string, ... )  do {         \
  extern void bebop_exit (const int errcode);                     \
  bebop_error ((__errclass_string), __VA_ARGS__);                   \
  bebop_log (0, "*** Fatal error in function %s at line %d of "   \
	     "file %s, aborting ***\n", __func__, __LINE__, __FILE__); \
  bebop_exit (EXIT_FAILURE);                                      \
} while(0)

#define bebop_log_fnstart( __level )   do {        \
  bebop_log ((__level), "=== %s ===\n", __func__); \
} while(0)
#define bebop_log_fnend( __level )   do {        \
  bebop_log ((__level), "=== Done with %s ===\n", __func__); \
} while(0)


#endif /* _log_h */
