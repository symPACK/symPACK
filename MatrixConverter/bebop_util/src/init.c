/**
 * @file init.c
 * @author Mark Hoemmen
 * @since 23 Nov 2007
 * @date Time-stamp: <2008-08-19 09:34:14 mhoemmen>
 *
 * Implementation of BeBOP Utility Library initialization functions.
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
 *************************************************************************/
#include <bebop/util/config.h>

#include <bebop/util/init.h>
#include <bebop/util/log.h>
#include <bebop/util/timer.h>
#include <bebop/util/random_number.h>
#ifdef USE_ECL
#  include <bebop/util/ecl_interface.h>
#endif

#ifdef USE_PTHREADS
#  include <pthread.h>
#endif /* USE_PTHREADS */

#include <assert.h>
#include <errno.h>  
#include <limits.h>    /* INT_MAX */
#include <stdio.h>     /* the usual output stuff */
#include <stdlib.h>    /* strtol, getenv, ... */

/**
 * Current debug level.  Not to be read or modified directly by users.
 * The higher the number, the more events are printed.  Zero is the 
 * minimum.
 */
static int BEBOP_DEBUG_LEVEL = 0;

void
bebop_default_initialize (int argc, char** argv, int* const info)
{
  unsigned long seed; /* for the pseudorandom number generator */

  /* Read the BEBOP_DEBUG_LEVEL environment variable and use it to set
     the debug level */
  bebop_set_debug_level_from_environment ();

  /* Log to stderr; don't create a log file */
  bebop_start_logging (1, stderr);

  /* Start up the timer (low cost, medium usefulness) */
  init_timer ();

  /* Generate a seed for the custom PRNG */
  seed = bebop_random_seed ();

  /* FIXME: Start up the custom pseudorandom number generator */

#ifdef USE_ECL
  /* Start up ECL */
  lisp_boot (argc, argv, info);
#endif /* USE_ECL */
}

void
bebop_exit (const int errcode)
{
  /* Close the logfile if we are supposed to.  Use the "emergency stop" 
   * version so that we don't get into an infinite loop by calling 
   * bebop_exit() if the mutex in bebop_stop_logging() fails. */
  bebop_emergency_stop_logging ();

  /*
   * FIXME (mfh 23 Nov 2007): maybe it would be better to remove this 
   * call to exit(), and register bebop_exit() with atexit().  However, 
   * this registration isn't guaranteed to succeed, as atexit() is only
   * required to support up to 32 registered functions (on Darwin (MacOS 
   * 10.4.11) at least).  Then we would need to check whether atexit()
   * actually registered the function, and keep a flag around, etc...
   * That kind of defeats the purpose, no? ;-P
   */
  exit (errcode);
}

void
bebop_set_debug_level_from_environment ()
{
  char* env_var_name = NULL;
  char* valstring = NULL;

  env_var_name = "BEBOP_DEBUG_LEVEL";
  valstring = getenv (env_var_name);

  if (valstring == NULL)
    {
      env_var_name = "SPMV_DEBUG_LEVEL";
      valstring = getenv (env_var_name);
    }

  if (valstring == NULL)
    {
      env_var_name = "SMVM_DEBUG_LEVEL";
      valstring = getenv (env_var_name);
    }

  if (valstring == NULL) 
    {
      /* default */
      bebop_set_debug_level (0);
      return;
    }
  else
    {
      int value;
#ifdef HAVE_ERRNO_H
      errno = 0; /* Must always reset errno before checking it */

      /* Convert the "value" part of the "name=value" pair to an
       * int */

      value = strtol (valstring, NULL, 10);
      if (errno)
	{
	  fprintf (stderr, "*** bebop_set_debug_level_from_environment:"
		   " %s is set to %s, which is not an integer ***\n", 
		   env_var_name, valstring);
	  errno = 0; /* reset errno */
	  bebop_set_debug_level (0);
	  return;
	}
      else
	{
	  if (value < 0)
	    {
	      fprintf (stderr, "*** bebop_set_debug_level_from_environment: "
		       "%s is set to %d, which is negative and"
		       " therefore invalid ***\n", env_var_name, value);
	      bebop_set_debug_level (0);
	      return;
	    }
	  else
	    bebop_set_debug_level (value);
	}
#else /* No errno.h */
      char *endptr = NULL;

      value = strtol (valstring, &endptr, 10);
      /* Check out the GNU man page for strtol.   "If endptr is not NULL, 
       * strtol() stores the address of the first invalid character in 
       * *endptr.  If there were no digits at all, strtol() stores the 
       * original value of nptr in *endptr (and returns 0).  In particular, 
       * if *nptr is not ‘\0’ but **endptr is ‘\0’ on return, the entire 
       * string is valid."  */
      if (*valstring != '\0' && endptr != NULL && endptr == '\0')
	{
	  /* We got value successfully. */
	  if (value < 0)
	    {
	      fprintf (stderr, "*** bebop_set_debug_level_from_environment: "
		       "%s is set to %d, which is negative and"
		       " therefore invalid ***\n", env_var_name, value);
	      bebop_set_debug_level (0);
	      return;
	    }
	  else
	    bebop_set_debug_level (value);
	}
      else
	{
	  fprintf (stderr, "*** bebop_set_debug_level_from_environment:"
		   " %s is set to %s, which is not an integer ***\n", 
		   env_var_name, valstring);
	  bebop_set_debug_level (0);
	  return;
	}
#endif /* HAVE_ERRNO_H */
    }
}

int
bebop_debug_level ()
{
  return BEBOP_DEBUG_LEVEL;
}

#ifdef USE_PTHREADS
/**
 * Mutex for bebop_set_debug_level().
 */
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
#endif /* USE_PTHREADS */

void 
bebop_set_debug_level (const int level)
{
#ifdef USE_PTHREADS
  int error = pthread_mutex_lock (&lock);
  if (error != 0)
    {
      fprintf (stderr, "*** bebop_debug_level: failed to lock mutex! ***\n");
      bebop_exit (EXIT_FAILURE);
    }
#endif /* USE_PTHREADS */

  BEBOP_DEBUG_LEVEL = level;

#ifdef USE_PTHREADS
  error = pthread_mutex_unlock (&lock);
  if (error != 0)
    {
      fprintf (stderr, "*** bebop_debug_level: failed to unlock mutex! ***\n");
      bebop_exit (EXIT_FAILURE);
    }
#endif /* USE_PTHREADS */
}


