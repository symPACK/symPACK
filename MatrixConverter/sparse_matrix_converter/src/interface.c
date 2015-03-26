/**
 * @file interface.c
 * @author Mark Hoemmen
 * @since 09 Jun 2006
 * @date Time-stamp: <2008-07-16 11:23:19 mhoemmen>
 *
 * Implementation of the interface defined in interface.h (which see).
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
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/sparse_matrix_ops.h>
#include <bebop/smc/csr_matrix.h>

#include <bebop/util/get_options.h>
#include <bebop/util/malloc.h>
#include <bebop/util/util.h>
#include <bebop/util/timer.h>

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>


struct sparse_matrix_t*
sp_load (const char* path, const char* fmt)
{
  struct sparse_matrix_t* A = NULL;
  double seconds = get_seconds ();
  A = load_sparse_matrix (sparse_matrix_file_format_string_to_enum (fmt), path);
  if (A == NULL)
    {
      fprintf (stderr, "*** sp_load: Failed to load sparse matrix! ***\n");
      return NULL;
    }
  seconds = get_seconds () - seconds;
  printf ("Loading the sparse matrix took %g seconds.\n", seconds); 
  return A;
}

int
sp_save (struct sparse_matrix_t* A, const char* path, 
      const char* fmt)
{
  double seconds = get_seconds ();
  save_sparse_matrix (path, A, sparse_matrix_file_format_string_to_enum (fmt));
  seconds = get_seconds () - seconds;
  printf ("Saving the sparse matrix took %g seconds.\n", seconds);
  return 0;
}

void
sp_format (struct sparse_matrix_t* A)
{
  if (A == NULL)
    printf ("Sparse matrix is NULL!\n");
  else
    printf ("%s\n", sparse_matrix_format_string (A));
}

int
sp_convert (struct sparse_matrix_t* A, const char* type)
{
  double seconds = get_seconds ();
  int errcode = sparse_matrix_convert (A, 
	sparse_matrix_storage_format_string_to_enum (type));
  seconds = get_seconds () - seconds;
  if (errcode != 0)
    {
      printf ("*** Failed to convert sparse matrix! ***\n");
      return errcode;
    }
  else
    {
      printf ("Converting the sparse matrix took %g seconds.\n", seconds);
      return 0;
    }
}


struct sparse_matrix_t* 
sp_mult (struct sparse_matrix_t* B, struct sparse_matrix_t* A)
{
  return sparse_matrix_matmatmult (B, A);
}


struct sparse_matrix_t* 
sp_triprod (struct sparse_matrix_t* RT, 
	    struct sparse_matrix_t* A, 
	    struct sparse_matrix_t* P)
{
  if (RT->format != CSR || A->format != CSR || P->format != CSR)
    return NULL;
 
  return create_sparse_matrix (CSR, (void*) csr_matrix_triple_product (RT->repr, A->repr, P->repr));
}

void
sp_destroy (struct sparse_matrix_t* A)
{
  destroy_sparse_matrix (A);
}

void
sp_print (struct sparse_matrix_t* A)
{
  print_sparse_matrix (stdout, A);
}
