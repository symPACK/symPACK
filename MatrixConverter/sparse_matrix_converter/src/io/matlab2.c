/**
 * @file matlab2.c
 * @author Mark Hoemmen
 * @since 09 May 2007
 * @date Time-stamp: <2008-07-16 11:26:54 mhoemmen>
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

#include <stdlib.h>
#include <limits.h>
#include <errno.h>
#include <stdio.h>

typedef struct __matlab_spmat_t {
  void* values;
  void* col_indices;
  void* row_indices;
  void* nelts;      /* number of matrix elements */
  void* first_value;
  int   read_first_value_p;
  size_t index_size;
  size_t value_size;
  list_t token_list;
  value_type_t value_type;
  manysplit_t ms;
} matlab_spmat_t;

static void
set_defaults (matlab_spmat_t* M)
{
  M->index_size = sizeof (uint32_t);
  M->value_size = sizeof (double);
  M->value_type = DOUBLE_REAL;
}

void
init_matlab_spmat (matlab_spmat_t* M)
{
  M->values = NULL;
  M->col_indices = NULL;
  M->row_indices = NULL;
  M->nelts = NULL;

  M->first_value = NULL;
  M->read_first_value_p = 0;
  M->token_list = list_create ();
  manysplit_init (&ms, 4);

  set_defaults (M);
}

void
deinit_matlab_spmat (matlab_spmat_t* M)
{
#define free_if_not_null( ptr )  do {		\
    if ((ptr) != NULL) {			\
      smvm_free((ptr));				\
      (ptr) = NULL;				\
    }						\
  } while(0)

  free_if_not_null (M->values);
  free_if_not_null (M->col_indices);
  free_if_not_null (M->row_indices);
  free_if_not_null (M->nelts);
  free_if_not_null (M->first_value);
  M->read_first_value_p = 0;
}

void
matlab_spmat_append_double_real_nonzero (matlab_spmat_t* M, list_t tokens)
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

int
matlab_spmat_deduce_type (matlab_spmat_t* M, FILE* f)
{
  char line[400];
  int nfields;

  /* Read the first line, and count the number of fields of data.
     This tells us the matrix type.  Save the value and indices from
     the first line. */
  checked_fgets (line, 400, f);
  manysplit (line, " \t", 2);
  manysplit_dtokens (&(M->token_list), &(M->ms));
  nfields = manysplit_ntokens (&ms);

  if (nfields == 2)
    {
      value_type = PATTERN;
      read_assert (0 == read_pattern_tokens (&i, &j, toklist, nfields));
    }
  else if (nfields == 3)
    {
      value_type = DOUBLE_REAL;
      read_assert (0 == read_double_tokens (&i, &j, &d, toklist, nfields));
    }
  else if (nfields == 4)
    {
      value_type = DOUBLE_COMPLEX;
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

}
