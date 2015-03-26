/**
 * @file random_number.c
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-16 10:11:03 mhoemmen>
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
 ******************************************************************/
#include <bebop/util/config.h>
#include <bebop/util/random_number.h>
#include <bebop/util/malloc.h>
#include <bebop/util/util.h>
#include <bebop/util/log.h>

#ifdef USE_PTHREADS
#  include <pthread.h>
#endif /* USE_PTHREADS */

#include <assert.h>
#include <limits.h>
#include <stdlib.h> 
#include <time.h>    /* time() */


/**
 * The random seed.
 */
static unsigned long the_seed = (unsigned long) 0;

#ifdef USE_PTHREADS
/**
 * Lock for protecting the random number generator's internal state.
 */
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
#endif /* USE_PTHREADS */


/* Initialize the Mersenne Twister PRNG. */
extern void
init_genrand (unsigned long s);

/* generates a random number on [0,0xffffffff]-interval */
extern unsigned long 
genrand_int32 (void);


/**
 * Mutex-protected genrand_int32
 */
static unsigned long
safe_genrand_int32 (void)
{
#ifdef USE_PTHREADS
  int error = 0;
  unsigned long retval = 0;

  error = pthread_mutex_lock (&lock);
  if (error != 0)
    bebop_fatal_error ("thread:mutex", "Failed to lock mutex!\n");
  retval = genrand_int32 ();  
  error = pthread_mutex_unlock (&lock);
  if (error != 0)
    bebop_fatal_error ("thread:mutex", "Failed to unlock mutex!\n");
  return retval;
#else
  return genrand_int32 ();
#endif /* USE_PTHREADS */
}

/* generates a random number on [0,0x7fffffff]-interval */
extern long 
genrand_int31 (void);

/* generates a random number on [0,1) with 53-bit resolution*/
extern double 
genrand_res53 (void);

/**
 * Mutex protected genrand_res53
 */
static double
safe_genrand_res53 (void)
{
#ifdef USE_PTHREADS
  int error = 0;
  double retval = 0.0;

  error = pthread_mutex_lock (&lock);
  if (error != 0)
    bebop_fatal_error ("thread:mutex", "Failed to lock mutex!\n");
  retval = genrand_res53 ();  
  error = pthread_mutex_unlock (&lock);
  if (error != 0)
    bebop_fatal_error ("thread:mutex", "Failed to unlock mutex!\n");
  return retval;
#else
  return genrand_res53 ();
#endif /* USE_PTHREADS */
}

#ifdef USE_PTHREADS
static pthread_once_t once_block = PTHREAD_ONCE_INIT;
/**
 * The POSIX standard provides for once-only initialization; the 
 * initialization routine must have the signature void (void).
 * Hence we make the seed (the argument to the MT initialization
 * routine init_genrand()) a static variable, "the_seed".
 */
static void
initialize ()
{
  init_genrand (the_seed);
}
#endif /* USE_PTHREADS */


void 
bebop_init_prng (unsigned long seed)
{
#ifdef USE_PTHREADS 
  int error = 0;
#endif /* USE_PTHREADS */

  bebop_log_fnstart (2);
  bebop_log (2, "seed = %x\n", seed);

#ifdef USE_PTHREADS
  error = pthread_mutex_lock (&lock);
  if (error != 0)
    bebop_fatal_error ("thread:mutex", "Failed to lock mutex!\n");
#endif /* USE_PTHREADS */

  the_seed = seed;

#ifdef USE_PTHREADS
  error = pthread_mutex_unlock (&lock);
  if (error != 0)
    bebop_fatal_error ("thread:mutex", "Failed to unlock mutex!\n");
#endif /* USE_PTHREADS */

  /*
   * If we use pthreads, we use the once-only initialization
   * technique so that the random number generator's state is
   * only initialized once.
   */
#ifdef USE_PTHREADS
  pthread_once (&once_block, initialize);
#else
  init_genrand (the_seed);
#endif /* USE_PTHREADS */

  bebop_log (2, "the_seed = %x\n", the_seed);
  bebop_log_fnend (2);
}

/**
 * Returns a random integer in the range 0 .. N-1, inclusive.
 */
static int
bebop_random_integer_in_range_zero_to_n_minus_one (int N)
{
  int num_trials = 0;
  long attempt = 0;
  int n;

  assert (N > 0);

  /*
   * See:  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/efaq.html
   *
   * 1. Find the min int n such that N <= 2^n.  INT_MAX is 2^31 - 1, so we 
   *    don't need any special cases in the algorithm below (these are signed
   *    ints).
   * 2. Take the most significant n bits of as many integer random numbers as
   *    we need to get n bits.  Since INT_MAX is 2^31 - 1, n will be at most
   *    32, so we only need 32 bits.
   * 3. If the random number is >= N, discard it and try again.  (Just in case
   *    the random number generator is broken, I have a loop that stops trying
   *    again if we don't have success after INT_MAX trials.)
   */

  /* Step 1:  Find n. */
  if (N >= INT_MAX)
    n = 31;
  else
    {
      register int D = 1;
      n = 0;
      while (N > D)
	{
	  D = D << 1;
	  n++;
	}
    }

  assert (n < 32);
  assert (n >= 0);

  while (num_trials < INT_MAX)
    {
      /* Step 2:  get 32 bits of randomness */
      unsigned long random_bits = safe_genrand_int32 ();

      /* Right-shift random_bits by (32 - n).  This fills in zeros for the 
       * sign bit, so we know that attempt >= 0. */
      if (n == 0)
	attempt = 0; /* the only possibility */
      else
        attempt = (long) (random_bits >> (32 - n));

      assert (attempt >= 0);

      /* Return attempt if it is in range. */
      if (attempt < N) 
	return attempt;

      num_trials++;
    }

  bebop_fatal_error ("PRNG", "The underlying pseudorandom number generator "
		     "is probably broken, since even after INT_MAX trials, "
		     "we were unable to get a random integer within the ran"
		     "ge [0,%d) ***\n", N);
  return 0; /* to pacify the compiler */
}


/***********************************************************************/
int 
bebop_random_integer (int low, int high)
{
  int k = 0;
  assert (high >= low);
  k = low + bebop_random_integer_in_range_zero_to_n_minus_one (high - low + 1);

#if 0
  assert (k >= low);
  assert (k <= high);
#endif /* 0 */
  return k;

#if 0
  /*
   * MFH 2004 Jan 15:  OK, this is too weird.  
   * double temp = 0.0; temp = (double) i / RAND_MAX;
   * causes temp = nan, but calculating (double) i / RAND_MAX again
   * after that makes everything work out.  You actually have to 
   * calculate (double) i / RAND_MAX _twice_  -- you can't just set
   * temp = 0.0.
   */
  /* rand() returns an int between 0 and RAND_MAX. */
  int i = rand ();
  /* double temp = (double) (i) / (double) (RAND_MAX); */
  /* double temp = (double) i / RAND_MAX; */
  double temp = 0.0;

  WITH_DEBUG4(fprintf(stderr, "### rand() = %d ###\n", i));

  temp = (double) i / RAND_MAX; 
  temp = (double) i / RAND_MAX; 

  /* Transform into interval [low,high], and round to nearest int.
   * This is essential -- if you just rely on the cast to int, the 
   * resulting truncation will throw off the probabilities.
   */
  return round_to_nearest_int ((high - low) * temp + low);
#endif
}


/***********************************************************************/
double
bebop_random_double (double low, double high)
{
  return low + (high - low) * safe_genrand_res53 ();
#if 0
  /* rand() returns an int between 0 and RAND_MAX. */
  double temp = (double) (rand ()) / (double) (RAND_MAX);

  /* Transform into interval [low, high]. */
  return (high - low) * temp + low;
#endif 
}


struct random_integer_from_range_without_replacement_generator_t*
create_random_integer_from_range_without_replacement_generator (const int low, const int high)
{
  struct random_integer_from_range_without_replacement_generator_t* gen = NULL;
  int i;

  bebop_log_fnstart (2);

  gen = bebop_calloc (sizeof (struct random_integer_from_range_without_replacement_generator_t), 1);
  gen->low = low;
  gen->high = high;

  /* It's legit for high < low -- that just means that there are no numbers 
   * left from which to pick. */
  if (high < low)
    gen->num_remaining = 0;
  else
    gen->num_remaining = high - low + 1;

  gen->remaining = bebop_malloc (gen->num_remaining * sizeof (int));
  for (i = 0; i < gen->num_remaining; i++)
    gen->remaining[i] = low + i;

  return gen;
}

void
destroy_random_integer_from_range_without_replacement_generator (struct random_integer_from_range_without_replacement_generator_t* gen)
{
  if (gen != NULL)
    {
      if (gen->remaining != NULL)
	bebop_free (gen->remaining);

      bebop_free (gen);
    }
}


int
return_random_integer_from_range_without_replacement (int* theint, struct random_integer_from_range_without_replacement_generator_t* gen)
{
  int p = 0;
  int i;

  if (gen->num_remaining < 1)
    return -1;

  if (bebop_debug_level() > 1)
    {
      int j;
      bebop_log (2, "### ");
      for (j = 0; j < gen->num_remaining; j++)
	bebop_log (2, "%d ", gen->remaining[j]);

      bebop_log (2, "\n");
    }

  /* Pick a random position in the array of remaining elements, and return 
   * the element there. */
  p = bebop_random_integer (0, gen->num_remaining - 1);
  *theint = gen->remaining[p];

  /* Remove the selected element from the array of remaining elements, by 
   * left-shifting over all the elements after it and decrementing 
   * num_remaining. */

  for (i = p; i < gen->num_remaining - 1; i++)
    gen->remaining[i] = gen->remaining[i+1];

  gen->num_remaining = gen->num_remaining - 1;

  if (bebop_debug_level() > 1)
    {
      int j;
      bebop_log (2, "### ");
      for (j = 0; j < gen->num_remaining; j++)
	bebop_log (2, "%d ", gen->remaining[j]);

      bebop_log (2, "\n");
    }

  return 0;
}


static inline int
__read_seed (unsigned long* seed, FILE* input)
{
  /* 
   * "char buffer[sizeof(unsigned long)]" compiles with "gcc -ansi
   * -pedantic" (C89 standard) as well as with "gcc -std=c99
   * -pedantic" (C99 standard).  I avoid assigning sizeof(unsigned
   * long) to a variable, as some compilers don't like stack array
   * declarations with a variable (even if the variable is constant).
   * "gcc -ansi -pedantic" rejects it; these options cause gcc to
   * reject any code that's not compliant with the ANSI C89 standard.
   */
  char buffer[sizeof(unsigned long)];
  const size_t nbytes_needed = sizeof (unsigned long);
  const size_t nbytes_read = fread (buffer, 
				    sizeof(char), 
				    nbytes_needed, 
				    input);
  if (nbytes_read != nbytes_needed)
    return -1;
  else
    {
      *seed = *((unsigned long*) buffer);
      return 0;
    }
}
static inline int
__read_seed_from_file (unsigned long* seed, const char filename[])
{
  FILE* f = fopen (filename, "r");
  if (f == NULL)
    return -1;
  else
    {
      const int errcode = __read_seed (seed, f);
      fclose (f);
      return errcode;
    }
}

#ifdef USE_PTHREADS
/**
 * Mutex for calling time(NULL) for computing the random seed.  We
 * only use time(NULL) as a fallback if other sources of entropy don't
 * exist, so the mutex may never be used.
 */
static pthread_mutex_t seed_compute_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif /* USE_PTHREADS */

static inline int
__read_seed_from_time (unsigned long* seed)
{
  time_t current_time;
#ifdef USE_PTHREADS
  int error;
#endif /* USE_PTHREADS */

#ifdef USE_PTHREADS
  error = pthread_mutex_lock (&seed_compute_mutex);
  if (error != 0)
    bebop_fatal_error ("thread:mutex", "Failed to lock mutex!\n");
#endif /* USE_PTHREADS */

  current_time = time (NULL);

#ifdef USE_PTHREADS
  error = pthread_mutex_unlock (&seed_compute_mutex);
  if (error != 0)
    bebop_fatal_error ("thread:mutex", "Failed to unlock mutex!\n");
#endif /* USE_PTHREADS */

  return (unsigned long) current_time;
}

unsigned long
bebop_random_seed ()
{
  /* The seed to return */
  unsigned long seed = (unsigned long) 0;
  /* Error code from various operations.  Zero is good. */
  int errcode = 0;

  /* 
   * Try to read some random bits from /dev/random.  The Darwin (MacOS
   * X) 10.4 documentation says that "[o]n Linux, /dev/urandom will
   * produce lower quality output if the entropy pool drains, while
   * /dev/random will prefer to block and wait for additional entropy
   * to be collected."  MacOS X allows (root) to write to /dev/random
   * in order to add what one thinks is entropy.
   */
  errcode = __read_seed_from_file (&seed, "/dev/random");
  if (errcode != 0)
    {
      /* Try to read some random bits from /dev/urandom */
      errcode = __read_seed_from_file (&seed, "/dev/random");
      if (errcode != 0)
        {
	  bebop_log (1, "Warning: Failed to open /dev/urandom or "
		     "/dev/urandom, so falling back on system clock for"
		     " random seed.  This is undesirable because the time"
		     " at which a program is run isn't usually random.");
	  errcode = __read_seed_from_time (&seed);
	  if (errcode != 0)
	    bebop_fatal_error ("PRNG", "Failed to read random seed "
			       "from time(NULL), which is the last "
			       "fallback source of entropy");
        }
    }
  return seed;
}

