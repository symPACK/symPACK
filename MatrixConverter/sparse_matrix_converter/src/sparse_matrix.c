/**
 * @file sparse_matrix.c
 * @author Mark Hoemmen
 * @since 26 May 2005
 * @date Time-stamp: <2008-07-16 11:24:50 mhoemmen>
 *
 * A wrapper representation of a general sparse matrix.
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
#include <bebop/smc/bcoo_matrix.h>
#include <bebop/smc/bcsr_matrix.h>
#include <bebop/smc/coo_matrix.h>
#include <bebop/smc/csc_matrix.h>
#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/jad_matrix.h>

#include <bebop/util/malloc.h>
#include <bebop/util/util.h>
#include <bebop/util/log.h>


/**
 * Returns nonzero if and only if "format" is one of the supported sparse matrix 
 * storage formats.
 */
static int
sparse_matrix_supported_storage_format_p (enum sparse_matrix_storage_format_t format)
{
  if (format == CSC || format == CSR || 
      format == COO || format == BCOO || 
      format == BCSR || format == JAD)
    return 1;
  else 
    return 0;
}


struct sparse_matrix_t*
create_sparse_matrix (enum sparse_matrix_storage_format_t format, 
		      void* repr)
{
  struct sparse_matrix_t* A = NULL;

  A = bebop_malloc (sizeof (struct sparse_matrix_t));
  A->format = format;
  A->repr = repr;
  return A;
}

int
destroy_sparse_matrix (struct sparse_matrix_t* A)
{
  if (A != NULL)
    {
      if (A->repr != NULL)
	{
	  if (A->format == CSC)
	    destroy_csc_matrix (A->repr);
	  else if (A->format == CSR)
	    destroy_csr_matrix (A->repr);
	  else if (A->format == COO)
	    destroy_coo_matrix (A->repr);
	  else if (A->format == BCOO)
	    destroy_bcoo_matrix (A->repr);
	  else if (A->format == BCSR)
	    destroy_bcsr_matrix (A->repr);
	  else if (A->format == JAD)
	    destroy_jad_matrix (A->repr);
	  else
	    {
	      bebop_log (0, "*** destroy_sparse_matrix: invalid matrix storage "
			"format %d ***\n", A->format);
	      return -1; 
	    }
	  A->repr = NULL;
	}
      bebop_free (A);
    }
  else
    {
      bebop_log (0, "*** destroy_sparse_matrix: input is NULL! ***\n");
      return -1;
    }
  return 0;
}


const char* 
sparse_matrix_format_string (struct sparse_matrix_t* A)
{
  enum sparse_matrix_storage_format_t format = A->format;

  if (format == CSC)
    return "CSC";
  else if (format == CSR)
    return "CSR";
  else if (format == COO)
    return "COO";
  else if (format == BCOO)
    return "BCOO";
  else if (format == BCSR)
    return "BCSR";
  else if (format == JAD)
    return "JAD";
  else 
    return "unknown";
}



struct sparse_matrix_t*
sparse_matrix_matmatmult (struct sparse_matrix_t* B, struct sparse_matrix_t* A)
{
  enum sparse_matrix_storage_format_t format = A->format;

  if (format != B->format)
    {
      bebop_log (0, "*** sparse_matrix_matmatmult: B and A do not have the "
		"same format! ***\n");
      return NULL;
    }
  else if (format != CSR)
    {
      bebop_log (0, "*** sparse_matrix_matmatmult: sparse matrix-matrix "
		"multiplication not supported for matrix format %s ***\n", 
		sparse_matrix_format_string (A));
      return NULL;
    }
  else 
    {
      struct csr_matrix_t* C = csr_matrix_matmatmult ((struct csr_matrix_t*) (B->repr), 
						      (struct csr_matrix_t*) (A->repr));
      if (C == NULL)
	return NULL;
      else 
	return create_sparse_matrix (CSR, C);
    }
}

