/**
 * @file ecl_interface.c
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-08-19 09:33:29 mhoemmen>
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
#include <bebop/util/config.h>
#include <bebop/util/ecl_interface.h>
#include <ecl/ecl.h>

#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h> /* atexit */



#ifdef USE_PTHREADS
#  include <pthread.h>
#endif /* USE_PTHREADS */


#ifdef USE_PTHREADS 
/**
 * Ensures that initialization is only called once.
 */
static pthread_once_t once_block = PTHREAD_ONCE_INIT;
/**
 * Protects the correct initialization flag (see below).
 */
static pthread_mutex_t init_flag_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif /* USE_PTHREADS */

/**
 * Set to 1 in the once_block's routine if initialization was done
 * correctly, else set to 0 there.  We protect this with a mutex (see
 * above) if USE_PTHREADS is defined, so that lisp_boot() is
 * reentrant.  The pthread_once() will still ensure that lisp_boot()
 * only initializes ECL once, even if lisp_boot() is called multiple
 * times.  That means initialized_correctly_flag will only be set
 * once.  However, two threads could call lisp_boot() and thus one
 * thread could try to read initialized_correctly_flag before the
 * other thread, still in the once block, gets to write to it.  Thus,
 * we protect it with a mutex.
 */
static int initialized_correctly_flag = 0;


/*
 * On 10 July 2008, we borrowed (slightly modified) a lot of code from
 * Dustin Long's "EclGui", a project that uses ECL as the back end of
 * a Windows Lisp REPL widget:
 * 
 * http://www.progmatism.com/hacks/eclgui/index.php
 * 
 * This code was released into the public domain without an explicit
 * license, but with the implicit intention to permit derivative
 * works.
 */

/*
 * Unique symbols to represent error conditions and other things
 */
static cl_object g_read_error_symbol;
static cl_object g_eval_error_symbol;
static cl_object g_read_from_string_symbol;
/**
 * A Boolean function taking two arguments (x,y).  x is a Lisp integer
 * which may or may not be a bignum, and y is a string representation
 * of an integer from the C world.  Returns true iff the two integers
 * are numerically equal (EQL).
 */
static cl_object g_compare_integer_and_string_fn;

/**
 * Construct a function call to perform the read operation safely
 */
static cl_object 
lisp_safe_read (const char* const string_cmd);

static int
lisp_string_p (cl_object obj);

/**
 * Attempt to guess the file mode (smm_input, smm_output, or smm_io)
 * of f, and return the result.  If successful, set *info to zero,
 * else set *info to nonzero.
 */
static enum ecl_smmode
ecl_file_mode (FILE* f, int* const info);



/* ****************************************************************** */
/* ****************************************************************** */
/* ****************************************************************** */



/**
 * Initialize unique symbols to represent error conditions, and also
 * perform Lisp-space computation that we only want to do once.
 * Called by lisp_boot().
 */
static void
lisp_init_globals ()
{
  g_read_error_symbol = cl_gensym (0);
  g_eval_error_symbol = cl_gensym (0);

  /* We only want to make these objects once, to avoid the overhead of
     calling into Lisp */
  g_read_from_string_symbol = c_string_to_object ("READ-FROM-STRING");
  g_compare_integer_and_string_fn = 
    lisp_safe_read("(lambda (possible-bignum integer-string) "
                     "(eql (let ((*read-eval* nil)) "
	                    "(read-from-string integer-string)) "
                          "possible-bignum))");

#ifdef USE_PTHREADS
  /* Protect access to initialized_correctly_flag.  If the lock
     doesn't succeed, we can't protect the flag in order to write to
     it so that we can notify users of the error, so just try to abort
     cleanly.  */
  if (pthread_mutex_lock (&init_flag_mutex) != 0)
    {
      /* We haven't registered this as an exit() handler yet (see
	 lisp_boot()), so call it now */
      /* fprintf (stderr, "Failed to lock in lisp_init_globals!\n"); */
      (void) cl_shutdown ();
      exit (EXIT_FAILURE);
    }
  else
    {
      /* Check whether the above two Lisp objects got created.
	 lisp_error_p() only uses g_*_error_symbol, so we can now call
	 it safely. */
      if (! lisp_error_p (g_read_from_string_symbol, NULL) &&
	  ! lisp_error_p (g_compare_integer_and_string_fn, NULL))
	initialized_correctly_flag = 1;
      else 
	initialized_correctly_flag = 0;

      if (pthread_mutex_unlock (&init_flag_mutex) != 0)
	{
	  /* We haven't registered this as an exit() handler yet (see
	     lisp_boot()), so call it now */
	  /* fprintf (stderr, "Failed to unlock in lisp_init_globals!\n"); */
	  (void) cl_shutdown ();
	  exit (EXIT_FAILURE);
	}
    }
#else /* NOT USE_PTHREADS */
  /* Check whether the above two Lisp objects got created.
     lisp_error_p() only uses g_*_error_symbol, so we can now call
     it safely. */
  if (! lisp_error_p (g_read_from_string_symbol, NULL) &&
      ! lisp_error_p (g_compare_integer_and_string_fn, NULL))
    initialized_correctly_flag = 1;
  else 
    initialized_correctly_flag = 0;
#endif /* USE_PTHREADS */
}


/**
 * Tell the Lisp system to shut down.
 */
static void
lisp_shutdown ()
{
  /* If you wanted to help the GC here, you could assign 0 to the g_*
     symbols, but generally we only call lisp_shutdown on full
     shutdown of the whole process, so a final GC wouldn't even be
     necessary. */

  /* cl_shutdown() returns int, but atexit() wants a void->void function */
  (void) cl_shutdown ();
}


void
lisp_boot (int argc, char** argv, int* const info)
{
  *info = 0;
  cl_boot (argc, argv);
  si_select_package (make_simple_base_string ("CL-USER"));

  /*
   * Make sure that the globals are initialized only once.
   */
#ifdef USE_PTHREADS
  pthread_once (&once_block, lisp_init_globals);

  /* Check whether the initialization in lisp_init_globals()
     succeeded.  Yay for C's short-circuiting sequential
     conditionals! */
  if (pthread_mutex_lock (&init_flag_mutex) != 0 ||
      ! initialized_correctly_flag ||
      pthread_mutex_unlock (&init_flag_mutex) != 0)
    *info = -1;
#else
  lisp_init_globals ();
  if (! initialized_correctly_flag)
    *info = -1;
#endif /* USE_PTHREADS */

  /* POSIX doesn't guarantee that atexit() will succeed, though it
     should always succeed on any reasonable implementation. */

  /* Register cl_shutdown() to call at exit.  This returns zero if
     successful, else nonzero. */
  if (atexit (&lisp_shutdown))
    *info = -2;
}



int
lisp_false (cl_object x)
{
  if (x == Cnil)
    return 1;
  else
    return 0;
}

int
lisp_true (cl_object x)
{
  if (x == Cnil)
    return 0;
  else
    return 1;
}


int32_t
lisp_object_to_int32 (cl_object x, int* const info)
{
  int32_t result = 0;
  char buf[22]; /* More than enough to hold an int32_t, including sign bit */
  cl_object str, samep;

  /* Straightforward conversion from Lisp integer to C int32_t */
  if (type_of (x) == t_fixnum)
    {
      result = (int32_t) fix(x);

      /* cl_fixnum is strictly contained in the C signed integer type
	 of the same number of bits, since fixnums use two of those
	 bits for type information.  So if x is a fixnum and it fits
	 in the return integer type (int32_t), then there was no
	 overflow and thus we can return without further checking. */
      if (sizeof(cl_fixnum) <= sizeof(int32_t))
	{
	  *info = 0;
	  return result;
	}
    }
  else if (type_of (x) == t_bignum)
    result = (int32_t) big_to_long (x);

  /* 
   * Above conversion may overflow without warning if x is a bignum,
   * and also may overflow if fixnums have a 64-bit representation (of
   * which ECL would use 62 bits for the integer itself).  Thus, it
   * doesn't suffice to check if x is a bignum.  We also can't pass
   * the result directly back into Lisp space as a fixnum, in case it
   * fits in an int32_t but not in a fixnum (say, if sizeof(cl_fixnum)
   * == sizeof(int32_t) and the result is big enough to overflow 30
   * bits).  Thus, we opt for an approach that is slow but safe and
   * relatively platform-independent: Print the int32_t result in C
   * space to a string, call Lisp's READ on it, and then compare the
   * result to x in Lisp space.  If the comparison fails, then the
   * conversion failed.
   */

  /* 
   * snprintf() returns >= its second argument when the output was
   * truncated, i.e., when the string was not long enough.  This
   * should never happen, but we check for it anyway. 
   */
  if (22 >= snprintf (buf, 22, "%d", result))
    {
      *info = -1;
      return result;
    }
  /* Convert to Lisp string */
  str = make_simple_base_string (buf);
  /* We use FUNCALL on an anonymous function so as not to pollute the
     Lisp namespace.  We've already constructed the function so it
     shouldn't be too expensive. */
  samep = cl_safe_eval (cl_list (4, 
				 c_string_to_object ("funcall"),
				 g_compare_integer_and_string_fn,
				 x,
				 str), 
			Cnil,
			g_eval_error_symbol);
  if (samep == g_eval_error_symbol)
    {
      *info = -2;
      return result;
    }
  else if (lisp_false (samep))
    {
      /* Overflow:  the C return value is not equal to the original Lisp integer */
      *info = +1;
      return result;
    }

  *info = 0;
  return result;
}


const char* const
lisp_string_to_C_string (cl_object obj, int* const info)
{
  /* "type_of(obj) == t_string" only works in a Unicode build.  */
  if (lisp_string_p (obj))
    {
      *info = 0;
      return ecl_base_string_pointer_safe(obj);
    }
  else
    {
      *info = -1;
      return NULL;
    }
}

/**
 * Broken because the C FILE* could have read ahead already and thus
 * consumed some data from the FILE*.
 */
#if 0
void
c_stream_to_lisp_stream (cl_object* lisp_stream,
			 FILE* c_stream,
			 int* const info)
{
  int fd = -1;

  /* Get integer file descriptor number for the given C stream.
     fileno must always return (under Linux at least), but may return
     -1 and set errno to EBADF if it can detect that the given stream
     (FILE*) is not valid. */
  errno = 0;
  fd = fileno (c_stream);
  if (errno == EBADF || fd == -1)
    {
      errno = 0;
      *info = -1;
      return;
    }   

  /* Make a Lisp stream from the given file descriptor.  In ECL, this
     calls fdopen() on the file descriptor; it's rather silly to go
     from a FILE* to an int and back to a FILE* again, but it gives
     the right level of abstraction from the details of how ECL works
     underneath.  Alternately, we could implement a function that
     creates an ECL stream directly from the C FILE* -- but that would
     make my code a derivative work of the ECL code, and thus force me
     to release the SMC under an LGPL license.  (I'd be happy to do
     so, except that it would displease some lawyers and manager
     types.)

     Note: it had better be an input stream of some kind, or you'll be
     sorry!!! */

  *lisp_stream = ecl_make_stream_from_fd (Cnil, fd, smm_input);
  *info = 0;
}
#endif /* 0 */	






static int
lisp_string_p (cl_object obj)
{
  cl_object form = cl_list (2, c_string_to_object("stringp"), obj);
  cl_object out = cl_safe_eval (form, Cnil, g_eval_error_symbol);
  if (out == Cnil || out == g_eval_error_symbol)
    return 0;
  else
    return 1;
}



void
lisp_open_file (cl_object *stream,
		const char* const filename,
		int* const info)
{
  cl_object form = cl_list (6, 
			    c_string_to_object("open"), 
			    make_base_string_copy((char*) filename), 
			    c_string_to_object(":direction"),
			    c_string_to_object(":input"),
			    c_string_to_object(":if-does-not-exist"),
			    c_string_to_object("nil"));
  assert (form != OBJNULL);

  /* FIXME: what if we don't have permission to read the file? could
     the funcall signal? */
  *stream = cl_safe_eval (form, Cnil, g_eval_error_symbol);
  if (*stream == g_eval_error_symbol)
    *info = -1;
  else
    *info = 0;
}



void
lisp_close_file (cl_object stream, int* const info)
{
  cl_object form = cl_list (2, c_string_to_object("close"), stream);
  cl_object output = cl_safe_eval (form, Cnil, g_eval_error_symbol);

  /* CLOSE returns T if successful, else implementation-dependent,
     according to the Common Lisp Hyperspec.  Thus we can't rely on
     the return value being not true on error. */

  if (output == g_eval_error_symbol)
    *info = -1;
  else
    *info = 0;
}





static cl_object 
lisp_safe_read (const char* const string_cmd)
{
  /* 
   * Convert string_cmd to a Lisp string (the cast is only to placate
   * ECL's API; ECL doesn't modify the string), create the read form,
   * and evaluate it, returning the global constant error symbol if
   * the read fails.
   */
  cl_object cmd = make_simple_base_string ((char*) string_cmd);
  return si_safe_eval (3, 
		       CONS (g_read_from_string_symbol, CONS (cmd, Cnil)),
		       Cnil, 
		       g_read_error_symbol);
}


/** 
 * Execute lisp code, with an option for safe execution 
 */
cl_object 
lisp_run (const char* const string_cmd, 
	  const int use_unsafe_reader)
{
  cl_object form;
  if (use_unsafe_reader)
    {
      /*
       * If we completely trust the lisp code, for example if we
       * generated it ourselves, it is okay to use the faster read
       * function. Bad input will lead to a crash.  The cast is only
       * to placate ECL's API; ECL doesn't modify the string.
       */
      form = c_string_to_object ((char*) string_cmd);
    }
  else
    {
      /*
       * If the lisp code has syntax errors (e.g., unbalanced
       * parentheses or misquoted strings), reading it normally will
       * crash ECL, so use this function instead.
       */
      form = lisp_safe_read (string_cmd);
      if (lisp_error_p (form, NULL))
	return form;
    }
  /* Evaluate the parsed lisp form and return result */
  return si_safe_eval (3, form, Cnil, g_eval_error_symbol);
}


static enum ecl_smmode
ecl_file_mode (FILE* f, int* const info)
{
  enum ecl_smmode smm;
  int fd = fileno (f);
  int flags = fcntl (fd, F_GETFL, 0);

  if (flags == -1)
    {
      *info = -1;
      return smm_input; /* have to return something */
    }

  *info = 0;
  if (O_ACCMODE && flags == O_RDONLY)
    return smm_input;
  else if (O_ACCMODE && flags == O_WRONLY)
    return smm_output;
  else if (O_ACCMODE && flags == O_RDWR)
    return smm_io;
  else
    {
      *info = -2;
      return smm_input; /* have to return something */
    }
}


cl_object 
make_lisp_stream_from_C_FILE (FILE* f)
{
  enum ecl_smmode smm;
  cl_object stream;
  int info = 0;

  smm = ecl_file_mode (f, &info);
  if (info != 0)
    return g_eval_error_symbol;
  else
    return ecl_make_stream_from_FILE (Cnil, f, smm);
}


int
lisp_error_p (cl_object obj, 
	      const char** errmsg)
{
  if (errmsg != NULL)
    {
      if (obj == g_read_error_symbol)
	*errmsg = "Error reading form";
      else if (obj == g_eval_error_symbol)
	*errmsg = "Error evaluating form";
    }
  return obj == g_read_error_symbol || obj == g_eval_error_symbol;
}


