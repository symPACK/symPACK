#ifndef _random_number_h
#define _random_number_h
/**
 * @file random_number.h
 * @author mfh
 * @date Time-stamp: <2008-07-16 10:19:21 mhoemmen>
 * 
 * Utility functions for creating random numbers.
 *
 * @note The functions here are (hopefully) thread-safe -- we wrap the
 * unlocked underlying random number generator library with pthreads
 * mutexes, and we wrap initialization with a once-only construct, so 
 * only the first thread to reach it gets to call it.
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
 ***************************************************************************/

/**
 * Generate a random seed to supply to bebop_init_prng(), in case you
 * don't want the same seed each time.  We get sizeof(unsigned long)
 * bytes from /dev/random if that (device) file exists, else we use the current
 * time as a fallback.
 *
 * @note This is NOT a cryptographically secure seed, but should be 
 * good enough for a Monte Carlo application.
 */
unsigned long
bebop_random_seed ();


/**
 * Call this exactly once in your program, before you generate any random 
 * numbers.  You have to specify the seed (as a 32-bit unsigned integer).  
 * Even if you want the same sequence of random numbers each time, you have 
 * to call this function.  Furthermore, this function may not be implemented
 * (in fact, is probably not implemented) using rand() and srand(), so you 
 * shouldn't count on srand() having been called.  But if you are, why are 
 * you using this library anyway?
 *
 * @param seed [IN]  Seed for the random number generator.
 */
void
bebop_init_prng (unsigned long seed);

/**
 * Returns a random integer in [low, high].  
 */
int  
bebop_random_integer (int low, int high);


/**
 * Returns a random double in the range [low, high).  This works by first 
 * generating a random double in the range [0,1), and then scaling it to fit.
 * The random number generator used ensures that the double in the range [0,1)
 * has 53 random bits, which is all that a double-precision number can provide.
 *
 * \warn NOTE THE CHANGE IN INTERFACE -- the range is no longer a closed set!!!
 *
 * @param [IN] low
 * @param [IN] high
 */
double
bebop_random_double (double low, double high);


/**
 * This data structure is used to generate random integers out of a range 
 * [low, high] without replacement.  You should consider it an opaque data 
 * structure, and not read from or write to its internals.
 */
struct 
random_integer_from_range_without_replacement_generator_t
{
  int low, high;  
  int* remaining;
  int num_remaining;
};




/**
 * If you want to choose random integers out of the range [low, high] without 
 * replacement, call this function and use the returned object in successive 
 * calls to extract a random integer out of the range without replacement.
 */
struct random_integer_from_range_without_replacement_generator_t*
create_random_integer_from_range_without_replacement_generator (const int low, const int high);

/**
 * Frees a random_integer_from_range_without_replacement_generator_t struct, 
 * once you are done with it.
 */
void
destroy_random_integer_from_range_without_replacement_generator (struct random_integer_from_range_without_replacement_generator_t* gen);

/**
 * Given an initialized
 * random_integer_from_range_without_replacement_generator_t object,
 * attempts to draw a random integer out of the integers that are
 * remaining.  If an integer remains to choose, it is chosen and
 * assigned to *theint, and zero is returned.  If no integers remain,
 * nothing is assigned to *theint, and nonzero is returned.
 */
int
return_random_integer_from_range_without_replacement (int* theint, struct random_integer_from_range_without_replacement_generator_t*);



#endif /* NOT _random_number_h */

