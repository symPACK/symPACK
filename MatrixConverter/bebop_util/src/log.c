/**
 * @file log.c
 * @author Mark Hoemmen
 * @since 25 Jul 2006
 * @date Time-stamp: <2008-07-16 10:10:24 mhoemmen>
 *
 * Implementation of a thread-safe (serialized) logger.
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
#include <bebop/util/log.h>
#include <bebop/util/util.h>

#include <stdio.h>
#include <assert.h>
#include <time.h>

#ifdef USE_PTHREADS
#include <pthread.h>
#endif /* USE_PTHREADS */

/**
 * Pointer to the logfile.
 */
static FILE* bebop_logfile = NULL;

/**
 * If nonzero, then we are (the library is) responsible for closing the 
 * logfile when bebop_stop_logging() is called.  If zero, the logfile is 
 * not closed (this is useful e.g. if we want to log to stderr).
 */
static int bebop_must_close_logfile = 0;
/**
 * Nonzero if the logging initialization is given a stream (FILE*) to 
 * which to log, zero if it is given a log filename to which to log.
 */
static int bebop_use_stream = 0;
/**
 * If a log file is supplied (instead of a stream), this is the filename 
 * of the log file.
 */
static char* bebop_log_filename = NULL;
/** 
 * Nonzero if logging has started and not been stopped, else zero 
 */
static int bebop_logging = 0; 
/**
 * If using a logfile, nonzero if logging should append to any existing 
 * file and 0 if logging should truncate.
 */
static int bebop_log_append = 0;

#ifdef USE_PTHREADS 
/**
 * Ensures that the initialize() function is only called once.
 */
static pthread_once_t once_block = PTHREAD_ONCE_INIT;
/**
 * Lock to serialize logging.
 */
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
#endif /* USE_PTHREADS */

/**
 * Initialize logging.  It's a void->void function because of the 
 * restrictions of pthread_once_t.
 */
static void
initialize ()
{
  if (bebop_logging != 0)
    {
      /* This is just a warning message */
      fprintf (stderr, "*** log.c: initialize has been called twice "
	       "without stopping logging in between!  Don\'t do that!  "
	       "I\'ll ignore it this time but don\'t do it again! ***\n");
      return;
    }
  if (bebop_use_stream)
    {
      bebop_must_close_logfile = 0;
      /* assume that bebop_logfile is already set */
      if (ferror (bebop_logfile))
	{
	  fprintf (stderr, "*** log.c: Logfile stream is invalid!!! ***\n");
	  bebop_logging = 0;
	  return;
	}
    }
  else
    {
      bebop_must_close_logfile = 1;
      if (bebop_log_append)
	bebop_logfile = fopen (bebop_log_filename, "a");
      else 
	bebop_logfile = fopen (bebop_log_filename, "w");

      if (bebop_logfile == NULL || ferror (bebop_logfile))
	{
	  fprintf (stderr, "log.c: failed to open logfile %s ***\n", 
		   bebop_log_filename);
	  bebop_logging = 0;
	  return;
	}
    }
  /* Set to 1 to indicate that logging has started successfully. */
  bebop_logging = 1;
  /* 
   * Log an initial message to the beginning of the file, with the current time
   */
  {
    time_t tod; 
    if (-1 == time (&tod))
      {
        fprintf (stderr, "*** log.c initialize(): WARNING: failed to "
                "get current time via time() (in <time.h>) ***\n");
        bebop_log (1, "\n\n### BEGIN LOGGING (can\'t get current time)\n\n");
      }
    else
      {
        /*
         * It's OK to call ctime() (which is not thread-safe) because
         * this initialize() function is only called once by one thread.
         */
        bebop_log (1, "\n\n### BEGIN LOGGING %s\n\n", ctime (&tod));
      }
  } 
}

int
bebop_logging_p ()
{
  int retval = 0;

#ifdef USE_PTHREADS
  int error = pthread_mutex_lock (&lock);
  if (error != 0)
    {
      fprintf (stderr, "*** bebop_logging_p: failed to lock mutex! ***\n");
      bebop_exit (EXIT_FAILURE);
    }
#endif /* USE_PTHREADS */

  retval = bebop_logging;

#ifdef USE_PTHREADS
  error = pthread_mutex_unlock (&lock);
  if (error != 0)
    {
      fprintf (stderr, "*** bebop_logging_p: failed to unlock mutex! ***\n");
      bebop_exit (EXIT_FAILURE);
    }
#endif /* USE_PTHREADS */

  return retval;
}

int
bebop_start_logging (const int use_stream, ...)
{
  va_list ap;
#ifdef USE_PTHREADS
  int error = 0;
#endif /* USE_PTHREADS */

  /*
   * We use a lock to protect the "global" static variables that the 
   * initialize() function sees.  We have to protect them all as a 
   * group, not individually, so that one thread is responsible for
   * changing them all, and we won't get some values from some threads
   * and other values from other threads.  However, we don't care if 
   * one thread sets all those values and then some other thread calls
   * initialize(); the values should be the same for all threads.
   * (If they aren't, though, we're in trouble...)
   */

#ifdef USE_PTHREADS
  error = pthread_mutex_lock (&lock);
  if (error != 0)
    {
      fprintf (stderr, "*** bebop_start_logging: failed to lock mutex! ***\n");
      bebop_exit (EXIT_FAILURE);
    }
#endif /* USE_PTHREADS */

  va_start (ap, use_stream);

  bebop_use_stream = use_stream;
  if (use_stream)
    {
      bebop_logfile = va_arg (ap, FILE*);
      bebop_log_filename = NULL;
      bebop_log_append = 1; /* doesn't matter what this is */
    }
  else
    {
      bebop_logfile = NULL; /* will be set by initialize() */
      bebop_log_filename = va_arg (ap, char*);
      bebop_log_append = va_arg (ap, int);
    }
  va_end (ap);

#ifdef USE_PTHREADS
  error = pthread_mutex_unlock (&lock);
  if (error != 0)
    {
      fprintf (stderr, "*** bebop_start_logging: failed to unlock mutex! ***\n");
      bebop_exit (EXIT_FAILURE);
    }

  /* Call initialize() exactly once, no matter how many threads there are. */
  pthread_once (&once_block, initialize);
#else 
  /* Not using Pthreads, so just initialize the logger without thread safety */
  initialize ();
#endif /* USE_PTHREADS */

  if (bebop_logging)
    return 0;
  else
    return -1;
}

void
bebop_stop_logging ()
{
#ifdef USE_PTHREADS
  int error = 0;
#endif /* USE_PTHREADS */

#ifdef USE_PTHREADS
  error = pthread_mutex_lock (&lock);
  if (error != 0)
    {
      fprintf (stderr, "*** bebop_stop_logging: failed to lock mutex! ***\n");
      bebop_exit (EXIT_FAILURE);
    }
#endif /* USE_PTHREADS */

  if (bebop_logging != 0 && bebop_logfile != NULL && bebop_must_close_logfile)
    {
      fclose (bebop_logfile);
      bebop_logfile = NULL;
    }
  bebop_logging = 0;

#ifdef USE_PTHREADS
  error = pthread_mutex_unlock (&lock);
  if (error != 0)
    {
      fprintf (stderr, "*** bebop_stop_logging: failed to unlock mutex! ***\n");
      bebop_exit (EXIT_FAILURE);
    }
#endif /* USE_PTHREADS */
}

/**
 * Implementation of logging.  Not to be called by users.
 * Thread-safe.
 */
static void
bebop_log_helper (const char* fmt, va_list ap)
{
#ifdef USE_PTHREADS
  int error = pthread_mutex_lock (&lock);

  if (error != 0)
    {
      fprintf (stderr, "*** bebop_log_helper: failed to lock mutex! ***\n");
      bebop_exit (EXIT_FAILURE);
    }
#endif /* USE_PTHREADS */

  if (bebop_logging != 0)
    {
      assert (bebop_logfile != NULL && 0 == ferror (bebop_logfile));
      vfprintf (bebop_logfile, fmt, ap);
      /* Don't buffer output */
      fflush (bebop_logfile);
    }

#ifdef USE_PTHREADS
  error = pthread_mutex_unlock (&lock);
  if (error != 0)
    {
      fprintf (stderr, "*** bebop_log_helper: failed to unlock mutex! ***\n");
      bebop_exit (EXIT_FAILURE);
    }
#endif /* USE_PTHREADS */
}

static void
bebop_warn_helper (const char* class, const char* fmt, va_list ap)
{
#ifdef USE_PTHREADS
  int error = pthread_mutex_lock (&lock);

  if (error != 0)
    {
      fprintf (stderr, "*** bebop_warn_helper: failed to lock mutex! ***\n");
      bebop_exit (EXIT_FAILURE);
    }
#endif /* USE_PTHREADS */

  if (bebop_logging != 0)
    {
      assert (bebop_logfile != NULL && 0 == ferror (bebop_logfile));
      fprintf (bebop_logfile, "Warning %s: ", class);
      vfprintf (bebop_logfile, fmt, ap);
      /* Don't buffer output */
      fflush (bebop_logfile);
    }

#ifdef USE_PTHREADS
  error = pthread_mutex_unlock (&lock);
  if (error != 0)
    {
      fprintf (stderr, "*** bebop_warn_helper: failed to unlock mutex! ***\n");
      bebop_exit (EXIT_FAILURE);
    }
#endif /* USE_PTHREADS */
}

static void
bebop_error_helper (const char* class, const char* fmt, va_list ap)
{
#ifdef USE_PTHREADS
  int error = pthread_mutex_lock (&lock);

  if (error != 0)
    {
      fprintf (stderr, "*** bebop_error_helper: failed to lock mutex! ***\n");
      bebop_exit (EXIT_FAILURE);
    }
#endif /* USE_PTHREADS */

  if (bebop_logging != 0)
    {
      assert (bebop_logfile != NULL && 0 == ferror (bebop_logfile));
      fprintf (bebop_logfile, "*** Error %s: ", class);
      vfprintf (bebop_logfile, fmt, ap);
      fprintf (bebop_logfile, " ***\n");
      /* Don't buffer output */
      fflush (bebop_logfile);
    }

#ifdef USE_PTHREADS
  error = pthread_mutex_unlock (&lock);
  if (error != 0)
    {
      fprintf (stderr, "*** bebop_warn_helper: failed to unlock mutex! ***\n");
      bebop_exit (EXIT_FAILURE);
    }
#endif /* USE_PTHREADS */
}



void
bebop_warn (const char* class, const char* fmt, ...)
{
  va_list ap;

  va_start (ap, fmt);
  bebop_warn_helper (class, fmt, ap);
  va_end (ap);
}

void
bebop_error (const char* class, const char* fmt, ...)
{
  va_list ap;

  va_start (ap, fmt);
  bebop_error_helper (class, fmt, ap);
  va_end (ap);
}


void
bebop_log (int min_debug_level, const char* fmt, ...)
{
  va_list ap;

  va_start (ap, fmt);

  /**
   * FIXME (mfh 23 Nov 2007): This is not thread-safe, but it only 
   * reads the value and it only affects whether a log event is logged 
   * or not.  Thus, we consider it an acceptable tradeoff between
   * correctness and efficiency (since bebop_log() may be called 
   * frequently by different threads, whereas the debug level isn't 
   * usually set frequently at runtime).
   */
  if (bebop_debug_level () >= min_debug_level)
    bebop_log_helper (fmt, ap);

  va_end (ap);
}

void
bebop_emergency_stop_logging ()
{
  /*
   * Note that we don't use the once-only technique on this function,
   * but we do use a mutex to ensure that only one thread is in the 
   * body at a time, and we only close the logfile if it wasn't already
   * closed.  So effectively it works like a once-only function.
   */
#ifdef USE_PTHREADS
  int error = 0;
#endif /* USE_PTHREADS */

#ifdef USE_PTHREADS
  error = pthread_mutex_lock (&lock);
  if (error != 0)
    {
      fprintf (stderr, "*** bebop_emergency_stop_logging: failed to lock mutex! ***\n");
      /* Just ditch out; don't call bebop_exit(), else we'll get into an infinite loop,
       * because bebop_exit() calls bebop_emergency_stop_logging().  We'll just have to
       * deal with whatever errors result. */
      exit (EXIT_FAILURE);
    }
#endif /* USE_PTHREADS */

  if (bebop_logging != 0 && bebop_logfile != NULL && bebop_must_close_logfile)
    {
      fclose (bebop_logfile);
      bebop_logfile = NULL;
    }
  bebop_logging = 0;

#ifdef USE_PTHREADS
  error = pthread_mutex_unlock (&lock);
  if (error != 0)
    {
      fprintf (stderr, "*** bebop_emergency_stop_logging: failed to unlock mutex! ***\n");
      /* Just ditch out; don't call bebop_exit(), else we'll get into an infinite loop,
       * because bebop_exit() calls bebop_emergency_stop_logging(). */
      exit (EXIT_FAILURE);
    }
#endif /* USE_PTHREADS */
}


void
bebop_log_string (int min_debug_level, const char* s)
{
  bebop_log (min_debug_level, "%s", s);
}
