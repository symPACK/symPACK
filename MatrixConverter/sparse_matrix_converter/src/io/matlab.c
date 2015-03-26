/**
 * @file matlab.c
 * @author Mark Hoemmen
 * @since 09 May 2007
 * @date Time-stamp: <2008-07-16 11:26:51 mhoemmen>
 *
 * A new and incomplete Matlab-format sparse matrix file parser
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
#include <errno.h>
#include <stdio.h>

#include <convert.h> /* for string conversions */
#include <list.h>
#include <get_matches.h> /* regex utility functions */

#ifdef read_assert 
#undef read_assert
#define read_assert( boolval )  do {		\
    if (! boolval) {							\
      smvm_error ("io:format", "Incorrect format at line %d", linenum); \
      return -1;							\
    } \ 
} while(0)
#endif /* read_assert */

#ifdef checked_fgets
#undef checked_fgets
#define checked_fgets( line, infile )  do {	\
    errno = 0; \
    fgets (line, 400, infile); \
    if (errno != 0) \
      { \
	smvm_error ("io:read", "From fgets(), errno = %d at " \
		    "line %d of file", errno, linenum); \
	errno = 0; \
	return NULL; \
      } \
    errno = 0; \ 
    linenum++;					\
} while(0)
#endif /* checked_fgets */

#ifdef store_pattern
#undef store_pattern
#define store_pattern( i, j ) do { \
      if (ind >= maxind) \
	{ \
	  maxind = MAX (4, 2*maxind); \
	  II = smvm_realloc (II, maxind * sizeof (int)); \
	  JJ = smvm_realloc (JJ, maxind * sizeof (int)); \
	} \
      II[ind] = i;		   \
      JJ[ind] = j;		   \
      ind++;                       \
} while(0)

#ifdef store_double
#undef store_double
#define store_double( i, j, d ) do {		\
      if (ind >= maxind) \
	{ \
	  maxind = MAX (4, 2*maxind); \
	  II = smvm_realloc (II, maxind * sizeof (int)); \
	  JJ = smvm_realloc (JJ, maxind * sizeof (int)); \
	  val = smvm_realloc (val, maxind * sizeof (double)); \
	} \
      II[ind] = i;		   \
      JJ[ind] = j;		   \
      ((double*) val)[ind] = d;    \
      ind++;                       \
} while(0)

#ifdef store_complex
#undef store_complex
#define store_complex( i, j, re, im ) do {		\
      if (ind >= maxind) \
	{ \
	  maxind = MAX (4, 2*maxind); \
	  II = smvm_realloc (II, maxind * sizeof (int)); \
	  JJ = smvm_realloc (JJ, maxind * sizeof (int)); \
	  val = smvm_realloc (val, maxind * sizeof (double_Complex)); \
	} \
      II[ind] = i;		   \
      JJ[ind] = j;		   \
      ((double_Complex*) (val))[ind] = new_double_Complex(re, im);	\
      ind++;								\
} while(0)


static int
read_pattern_tokens (int* i, int *j, list_t tokens, int ntokens)
{
  if (ntokens != 2)
    return -1;
  else
    {
      list_t head = tokens;
      if (0 != string_to_int (i, list_car (head)))
	return -2;
      head = list_cdr (head);
      if (0 != string_to_int (j, list_car (head)))
	return -3;
      return 0;
    }
}

static int
read_double_tokens (int* i, int *j, double* d, list_t tokens, int ntokens)
{
  if (ntokens != 3)
    return -1;
  else
    {
      list_t head = tokens;
      if (0 != string_to_int (i, list_car (head)))
	return -2;
      head = list_cdr (head);
      if (0 != string_to_int (j, list_car (head)))
	return -3;
      head = list_cdr (head);
      if (0 != string_to_double (d, list_car (head)))
	return -4;
      return 0;
    }
}

static int
read_complex_tokens (int* i, int *j, double* re, double* im,
		     list_t tokens, int ntokens)
{
  if (ntokens != 4)
    return -1;
  else
    {
      list_t head = tokens;
      if (0 != string_to_int (i, list_car (head)))
	return -2;
      head = list_cdr (head);
      if (0 != string_to_int (j, list_car (head)))
	return -3;
      head = list_cdr (head);
      if (0 != string_to_double (re, list_car (head)))
	return -4;
      if (0 != string_to_double (im, list_car (head)))
	return -5;
      return 0;
    }
}

/**
 * Determine number of (non-unique) entries in Matlab-format sparse
 * matrix file, as well as the value type.
 */
static int
__nelts_matlab_sparse_matrix (int* nelts, value_type_t* valtype, FILE* infile)
{
  manysplit_t ms;
  list_t toklist = list_create ();
  char line[400];
  int linenum = 0;
  int nfields = 0;
  int __nelts = 0;
  value_type_t value_type;

  manysplit_init (&ms, 0);

  /* Read the first line, and figure out how many fields of data there
     are.  This tells us the matrix type. */
  checked_fgets (line, 400, infile);
  manysplit (line, " \t", 2);
  manysplit_dtokens (&toklist, &ms);
  nfields = manysplit_ntokens (&ms);
  if (nfields == 2)
    {
      value_type = PATTERN;
      read_assert (0 == read_pattern_tokens (&i, &j, toklist, nfields));
    }
  else if (nfields == 3)
    {
      value_type = REAL;
      read_assert (0 == read_double_tokens (&i, &j, &d, toklist, nfields));
    }
  else if (nfields == 4)
    {
      value_type = COMPLEX;
      read_assert (0 == read_complex_tokens (&i, &j, &re, &im, 
					     toklist, nfields));
    }
  else
    {
      smvm_error ("io:format", "First line of Matlab-format sparse "
		  "matrix file has %d fields, but only 2 - 4 fields "
		  "(inclusive) are allowed", nfields);
      return -1;
    }
  __nelts++;

  if (value_type == PATTERN)
    {
      while (! feof (infile))
	{
	  checked_fgets (line, 400, infile);
	  manysplit (line, " \t", 2);
	  manysplit_dtokens (&toklist, &ms);
	  nfields = manysplit_ntokens (&ms);
	  read_assert (0 == read_pattern_tokens (&i, &j, toklist, nfields));
	  __nelts++;
	}
    }
  else if (value_type == REAL)
    {
      while (! feof (infile))
	{
	  checked_fgets (line, 400, infile);
	  manysplit (line, " \t", 2);
	  manysplit_dtokens (&toklist, &ms);
	  nfields = manysplit_ntokens (&ms);
	  read_assert (0 == read_double_tokens (&i, &j, &d, toklist, nfields));
	  __nelts++;
	}
    }
  else if (value_type == COMPLEX)
    {
      while (! feof (infile))
	{
	  checked_fgets (line, 400, infile);
	  manysplit (line, " \t", 2);
	  manysplit_dtokens (&toklist, &ms);
	  nfields = manysplit_ntokens (&ms);
	  read_assert (0 == read_complex_tokens (&i, &j, &re, &im, 
						 toklist, nfields));
	  __nelts++;
	}
    }
  
  list_destroy (&toklist);
  manysplit_deinit (&ms);

  *nelts = __nelts;
  *valtype = value_type;
  return 0;
}


/**
 * First call __nelts_matlab_sparse_matrix() to get the number of
 * elements in the matrix (nelts), and then call this function.
 */
static int
__read_matlab_sparse_matrix (int** __II, int** __JJ, void** __val, 
			     FILE* infile, const int nelts,
			     const enum value_type_t valtype)
{
  manysplit_t ms;
  list_t toklist = list_create ();
  char line[400];
  int linenum = 0;
  int nfields = 0;
  int* II = smvm_malloc (nelts * sizeof (int));
  int* JJ = smvm_malloc (nelts * sizeof (int));
  void* val = NULL;
  int ind = 0;

  if (valtype == PATTERN)
    val = NULL;
  else if (valtype == REAL)
    val = smvm_malloc (nelts * sizeof (double));
  else if (valtype == COMPLEX)
    val = smvm_malloc (nelts * sizeof (double_Complex));
  else
    {
      smvm_error ("enum:value", "In __read_matlab_sparse_matrix, valtype "
		  "(which is a value_type_t enum) has invalid value %d", 
		  (int) valtype);
      return -1;
    }
  manysplit_init (&ms, 0);

  /* We know the matrix data type, so just start reading! */

  if (value_type == PATTERN)
    {
      int i, j;
      while (! feof (infile))
	{
	  checked_fgets (line, 400, infile);
	  manysplit (line, " \t", 2);
	  manysplit_dtokens (&toklist, &ms);
	  nfields = manysplit_ntokens (&ms);
	  read_assert (0 == read_pattern_tokens (&i, &j, toklist, nfields));
	  store_pattern (i, j);
	}
    }
  else if (value_type == REAL)
    {
      int i, j;
      double d;
      while (! feof (infile))
	{
	  checked_fgets (line, 400, infile);
	  manysplit (line, " \t", 2);
	  manysplit_dtokens (&toklist, &ms);
	  nfields = manysplit_ntokens (&ms);
	  read_assert (0 == read_double_tokens (&i, &j, &d, toklist, nfields));
	  store_double (i, j, d);
	}
    }
  else if (value_type == COMPLEX)
    {
      int i, j;
      double re, im;
      while (! feof (infile))
	{
	  checked_fgets (line, 400, infile);
	  manysplit (line, " \t", 2);
	  manysplit_dtokens (&toklist, &ms);
	  nfields = manysplit_ntokens (&ms);
	  read_assert (0 == read_complex_tokens (&i, &j, &re, &im, 
						 toklist, nfields));
	  store_complex (i, j, re, im);
	}
    }
  
  list_destroy (&toklist);
  manysplit_deinit (&ms);

  *nelts = __nelts;
  return 0;
}


/**
 * Bare implementation of reading Matlab-format sparse matrices from
 * a file.  Allocates and expands matrix storage arrays as necessary.
 */
static int
__extread_matlab_sparse_matrix (int* nelts, value_type_t* valtype,
				int** __II, int** __JJ, 
			        void** __val, FILE* infile)
{
  manysplit_t ms;
  list_t toklist = list_create ();
  char line[400];
  int linenum = 0;
  int nfields = 0;
  int* II = NULL;
  int* JJ = NULL;
  void* val = NULL;
  int ind = 0;
  int i, j;

  manysplit_init (&ms, 0);

  /* Read the first line, and figure out how many fields of data there
     are.  This tells us the matrix type. */
  checked_fgets (line, 400, infile);
  manysplit (line, " \t", 2);
  manysplit_dtokens (&toklist, &ms);
  nfields = manysplit_ntokens (&ms);
  if (nfields == 2)
    {
      value_type = PATTERN;
      read_assert (0 == read_pattern_tokens (&i, &j, toklist, nfields));
      extstore_pattern (i, j);
    }
  else if (nfields == 3)
    {
      double d;
      value_type = REAL;
      read_assert (0 == read_double_tokens (&i, &j, &d, toklist, nfields));
      extstore_double (i, j, d);
    }
  else if (nfields == 4)
    {
      double re, im;
      value_type = COMPLEX;
      read_assert (0 == read_complex_tokens (&i, &j, &re, &im, 
					     toklist, nfields));
      extstore_complex (i, j, re, im);
    }
  else
    {
      smvm_error ("io:format", "First line of Matlab-format sparse "
		  "matrix file has %d fields, but only 2 - 4 fields "
		  "(inclusive) are allowed", nfields);
      return -1;
    }

  if (value_type == PATTERN)
    {
      while (! feof (infile))
	{
	  checked_fgets (line, 400, infile);
	  manysplit (line, " \t", 2);
	  manysplit_dtokens (&toklist, &ms);
	  nfields = manysplit_ntokens (&ms);
	  
	  read_assert (0 == read_pattern_tokens (&i, &j, toklist, nfields));
	  extstore_pattern (i, j);
	}
    }
  else if (value_type == REAL)
    {
      double d;
      while (! feof (infile))
	{
	  checked_fgets (line, 400, infile);
	  manysplit (line, " \t", 2);
	  manysplit_dtokens (&toklist, &ms);
	  nfields = manysplit_ntokens (&ms);
	  
	  read_assert (0 == read_double_tokens (&i, &j, &d, toklist, nfields));
	  extstore_double (i, j, d);
	}
    }
  else if (value_type == COMPLEX)
    {
      double re, im;
      while (! feof (infile))
	{
	  checked_fgets (line, 400, infile);
	  manysplit (line, " \t", 2);
	  manysplit_dtokens (&toklist, &ms);
	  nfields = manysplit_ntokens (&ms);
	  
	  read_assert (0 == read_complex_tokens (&i, &j, &re, &im, 
						 toklist, nfields));
	  extstore_complex (i, j, re, im);
	}
    }
  
  list_destroy (&toklist);
  manysplit_deinit (&ms);

  *__II = II;
  *__JJ = JJ;
  *__val = val;
  *nelts = ind;
  *valtype = value_type;
  return 0;
}


typedef struct __info_t {
  int nelts;
  value_type_t value_type;
  int save_memory_flag;
} info_t;

static int
set_info (info_t* info, const char* str)
{
  list_t matches;
  int result;

  info->nelts = -1;
  info->save_memory_flag = 0;

  result = regcomp (&match_nelts, 
		    "nelts[:space:]*=[:space:]*([:digit:]+)($|[^[:digit:]])", 
		    REG_EXTENDED | REG_ICASE);
  assert (result == 0);
  matches = get_matches (&match_nelts, str, 0); 
  if (! list_empty_p (matches))
    {
      result = string_to_int (&(info->nelts), 
			      (char*) (list_car (matches)));
      if (result != 0)
	{
	  smvm_warning ("io:format", "set_info: Unable to convert "
			"string %s to integer", (char*) (list_car (matches)));
	  return -1;
	}
    }
  list_destroy (matches);
  regfree (&match_nelts);

  result = regcomp (&match_nelts, 
		    "save_memory_flag[:space:]*=[:space:]*(T|True|Yes)",
		    REG_EXTENDED | REG_ICASE);
  assert (result == 0);
  matches = get_matches (&match_nelts, str, 0); 
  if (list_empty_p (matches))
    info->save_memory_flag = 0;
  else
    info->save_memory_flag = 1;

  list_destroy (matches);
  regfree (&match_nelts);
  return 0;
}

int
read_matlab_sparse_matrix (int* nelts, value_type_t* valtype,
			   int** __II, int** __JJ, void** __val, 
			   FILE* infile,
			   const char* info)
{
  regex_t match_nelts;
  info_t strinfo;
  int result = 0;

  result = set_info (&strinfo, info);
  if (result != 0)
    smvm_warning ("io:format", "read_matlab_sparse_matrix: set_info "
		  "failed with error code %d -- reverting to default"
		  " indications", result);

  if (strinfo.save_memory_flag)
    {
      if (strinfo.nelts == -1)
        {
          result = __nelts_matlab_sparse_matrix (nelts, valtype, infile);
	  if (result != 0)
	    {
	      smvm_error ("io:format", "Unable to determine number "
			  "of nonzeros in Matlab-format sparse matrix: "
			  "error code %d", result);
	      return result;
	    }
          if (*nelts != strinfo.nelts)
	    smvm_warning ("io:format", "read_matlab_sparse_matrix: "
			  "incorrect guess %d for nelts; correct ne"
			  "lts = %d\n", strinfo.nelts, nelts);
        }
      result = __read_matlab_sparse_matrix (__II, __JJ, __val, infile, 
				            *nelts, *valtype);
      if (result != 0)
	{
	  smvm_error ("io:format", "Error in reading Matlab-format "
		      "sparse matrix, though we were able to determine "
		      "number of nonzeros as %d.  Error code %d", 
		      *nelts, result);
	  return result;
	}
    }
  else
    {
      result = __extread_matlab_sparse_matrix (nelts, valtype, __II, __JJ,
					       __val, infile);
      if (result != 0)
	{
	  smvm_error ("io:format", "Error in reading Matlab-format "
		      "sparse matrix; error code %d", result);
	  return result;
	}
    }
  return 0;
}
