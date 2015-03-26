/**
 * @file sparse_matrix_ops.c
 * @author Mark Hoemmen
 * @since 26 May 2005
 * @date Time-stamp: <2008-07-16 11:25:08 mhoemmen>
 *
 * Load and save operations on general sparse matrices.
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
#include <bebop/smc/sparse_matrix_ops.h>
#include <bebop/smc/bcoo_matrix.h>
#include <bebop/smc/bcsr_matrix.h>
#include <bebop/smc/coo_matrix.h>
#include <bebop/smc/csc_matrix.h>
#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/jad_matrix.h>
#include <bebop/smc/read_mm.h>
#include <bebop/smc/sparse_matrix.h>

#include <bebop/util/malloc.h>
#include <bebop/util/util.h>
#include <bebop/util/log.h>

#include <stdio.h>
#include <string.h>

int
print_sparse_matrix (FILE* out, struct sparse_matrix_t* A)
{
  if (A->format == BCOO)
    return print_bcoo_matrix_in_matrix_market_format (out, A->repr);
  else if (A->format == BCSR)
    return print_bcsr_matrix_in_matrix_market_format (out, A->repr);
  else if (A->format == COO)
    return print_coo_matrix_in_matrix_market_format (out, A->repr);
  else if (A->format == CSC)
    return print_csc_matrix_in_matrix_market_format (out, A->repr);
  else if (A->format == CSR)
    return print_csr_matrix_in_matrix_market_format (out, A->repr);
  else if (A->format == JAD)
    {
      bebop_log (0, "*** save_sparse_matrix: I haven\'t yet implemented "
		"saving jagged diagonal format sparse matrices to "
		"MatrixMarket files!  Would you like to write the "
		"routine for me?  Submit an e-mail to mhoemmen at cs dot "
		"berkeley dot edu ***\n");
      return -1;
    }

  /* sentinel */
  bebop_log (0, "*** print_sparse_matrix: invalid format %d ***\n", A->format);
  return -1;
}


int
save_sparse_matrix (const char* const filename, 
		    struct sparse_matrix_t* A, 
		    enum sparse_matrix_file_format_t file_format)
{
  if (file_format == HARWELL_BOEING)
    {
      if (A->format == BCOO)
	{
	  bebop_log (0, "*** save_sparse_matrix: Harwell-Boeing output not"
		    " yet supported for BCOO format sparse matrices! ***\n");
	  return -1;
	}
      else if (A->format == COO)
	{
	  bebop_log (2, "Saving COO sparse matrix to file \"%s\""
		    " in Harwell-Boeing format\n", filename);
	  return save_coo_matrix_in_harwell_boeing_format (filename, A->repr);
	}
      else if (A->format == CSC)
	{
	  bebop_log (2, "Saving CSC sparse matrix to file \"%s\""
		    " in Harwell-Boeing format\n", filename);
	  return save_csc_matrix_in_harwell_boeing_format (filename, A->repr);
	}
      else if (A->format == CSR)
	{
	  bebop_log (2, "Saving CSR sparse matrix to file \"%s\""
		    " in Harwell-Boeing format\n", filename);
	  return save_csr_matrix_in_harwell_boeing_format (filename, A->repr);
	}
      else if (A->format == JAD)
	{
	  bebop_log (0, "*** save_sparse_matrix: I don\'t know how to save"
		    " jagged diagonal format sparse matrices to Harwell-Boeing"
		    " files! ***\n");
	  return -1;
	}
      else
	{
	  bebop_log (0, "*** save_sparse_matrix: Unsupported sparse matr"
		    "ix internal storage format %d ***\n", A->format);
	  return -1;
	}
    }
  else if (file_format == MATRIX_MARKET)
    {
      if (A->format == BCOO)
	{
	  bebop_log (2, "Saving BCOO sparse matrix to file \"%s\" in "
		    "MatrixMarket format\n", filename);
	  return save_bcoo_matrix_in_matrix_market_format (filename, A->repr);
	}
      else if (A->format == BCSR)
	{
	  bebop_log (2, "Saving BCSR sparse matrix to file \"%s\" in "
		    "MatrixMarket format\n", filename);
	  return save_bcsr_matrix_in_matrix_market_format (filename, A->repr);
	}
      else if (A->format == COO)
	{
	  bebop_log (2, "Saving COO sparse matrix to file \"%s\" in "
		    "MatrixMarket format\n", filename);
	  return save_coo_matrix_in_matrix_market_format (filename, A->repr);
	}
      else if (A->format == CSC)
	{
	  bebop_log (2, "Saving CSC sparse matrix to file \"%s\" in "
		    "MatrixMarket format\n", filename);
	  return save_csc_matrix_in_matrix_market_format (filename, A->repr);
	}
      else if (A->format == CSR)
	{
	  bebop_log (2, "Saving CSR sparse matrix to file \"%s\" in "
		    "MatrixMarket format\n", filename);
	  return save_csr_matrix_in_matrix_market_format (filename, A->repr);
	}
      else if (A->format == JAD)
	{
	  bebop_log (0, "*** save_sparse_matrix: I don\'t know how to save"
		    " jagged diagonal format sparse matrices to MatrixMarket"
		    " files! ***\n");
	  return -1;
	}
      else
	{
	  bebop_log (0, "*** save_sparse_matrix: Unsupported sparse matri"
		    "x internal storage format %d ***\n", A->format);
	  return -1;
	}
    }
  else if (file_format == MATLAB)
    {
      if (A->format == BCOO)
	{
	  bebop_log (0, "*** save_sparse_matrix: Matlab output not yet"
		    "supported for BCOO format sparse matrices! ***\n");
	  return -1;
	}
      else if (A->format == BCSR)
	{
	  bebop_log (0, "*** save_sparse_matrix: Matlab output not yet"
		   "supported for BCSR format sparse matrices! ***\n");
	  return -1;
	}
      else if (A->format == COO)
	{
	  bebop_log (2, "Saving COO sparse matrix to file \"%s\" in "
		    "Matlab format\n", filename);
	  return save_coo_matrix_in_matlab_format (filename, A->repr);
	}
      else if (A->format == CSC)
	{
	  bebop_log (2, "Saving CSC sparse matrix to file \"%s\" in "
		    "Matlab format\n", filename);
	  return save_csc_matrix_in_matlab_format (filename, A->repr);
	}
      else if (A->format == CSR)
	{
	  bebop_log (2, "Saving CSR sparse matrix to file \"%s\" in "
		    "Matlab format\n", filename);
	  return save_csr_matrix_in_matlab_format (filename, A->repr);
	}
      else if (A->format == JAD)
	{
	  bebop_log (2, "*** save_sparse_matrix: I don\'t know how to save"
		    " jagged diagonal format sparse matrices to Matlab files! "
		    "***\n");
	  return -1;
	}
      else
	{
	  bebop_log (0, "*** save_sparse_matrix: Unsupported sparse matri"
		    "x internal storage format %d ***\n", A->format);
	  return -1;
	}
    }
  else
    {
      bebop_log (0, "*** save_sparse_matrix: Invalid sparse matrix output"
		" file format %d ***\n", file_format);
      return -2;
    }

  bebop_log (0, "*** save_sparse_matrix: BUG: we should not have reached "
	    "this point! ***\n");
  return -3;
}


struct sparse_matrix_t*
load_sparse_matrix (enum sparse_matrix_file_format_t file_format, 
		    const char *const matrix_filename)
{
  struct sparse_matrix_t* A = NULL;
  int errcode = 0;

  if (file_format == HARWELL_BOEING)
    {
      struct csc_matrix_t* repr = bebop_malloc (sizeof (struct csc_matrix_t));

      bebop_log (2, "Loading sparse matrix from file \"%s\" "
		"in Harwell-Boeing format\n", matrix_filename);
      errcode = read_harwell_boeing_mat_double (matrix_filename, repr);

      if (errcode != 0)
	{
	  bebop_free (repr);
	  bebop_log (0, "*** load_sparse_matrix: Failed to read Harwell-Bo"
		    "eing format sparse matrix from file \"%s\", error code %d"
		    " ***\n", matrix_filename, errcode);
	  return NULL;
	}

      A = create_sparse_matrix (CSC, repr);
      if (A == NULL)
	{
	  bebop_log (0, "*** load_sparse_matrix: Unknown error initializin"
		    "g sparse matrix wrapper using Harwell-Boeing format spars"
		    "e matrix; probably a bug in create_sparse_matrix! ***\n");
	  destroy_csc_matrix (repr);
	  return NULL;
	}
    }
  else if (file_format == MATRIX_MARKET)
    {
      struct coo_matrix_t* repr = bebop_malloc (sizeof (struct coo_matrix_t));

      bebop_log (2, "Loading sparse matrix from file \"%s\" "
		"in MatrixMarket format\n", matrix_filename);

      /* Read COO-format matrix from file */
      errcode = read_matrix_market_sparse (matrix_filename, repr);
      if (errcode != 0)
	{
	  bebop_free (repr);
	  bebop_log (0, "*** load_sparse_matrix: Failed to read MatrixMark"
		    "et format sparse matrix from file \"%s\", error code %d *"
		    "**\n", matrix_filename, errcode);
	  return NULL;
	}

      A = create_sparse_matrix (COO, repr);
      if (A == NULL)
	{
	  bebop_log (0, "*** load_sparse_matrix: Unknown error initializing"
		    " sparse matrix wrapper using MatrixMarket format sparse ma"
		    "trix; probably a bug in create_sparse_matrix! ***\n");
	  destroy_coo_matrix (repr);
	  return NULL;
	}
    }
  else if (file_format == MATLAB)
    {
      struct coo_matrix_t* repr = NULL;

      bebop_log (2, "Loading sparse matrix from file \"%s\" "
		"in Matlab format\n", matrix_filename);

      /* Read COO-format matrix from file */
      repr = load_coo_matrix_in_matlab_format (matrix_filename);
      if (repr == NULL)
	{
	  bebop_log (0, "*** load_sparse_matrix: Failed to read Matlab"
		    "format sparse matrix from file \"%s\", error code %d *"
		    "**\n", matrix_filename, errcode);
	  return NULL;
	}

      A = create_sparse_matrix (COO, repr);
      if (A == NULL)
	{
	  bebop_log (0, "*** load_sparse_matrix: Unknown error initializing"
		   " sparse matrix wrapper using Matlab format sparse matrix; "
		   "probably a bug in create_sparse_matrix! ***\n");
	  destroy_coo_matrix (repr);
	  return NULL;
	}
    }
  else
    {
      bebop_log (0, "*** load_sparse_matrix: Unknown sparse matrix file for"
		"mat %d ***\n", file_format);
      return NULL;
    }

  return A;
}


int
sparse_matrix_expand_symmetric_storage (struct sparse_matrix_t* A)
{
  int errcode = 0;

  if (A->format == CSC)
    errcode = csc_matrix_expand_symmetric_storage ((struct csc_matrix_t*) A->repr);
  else if (A->format == COO)
    errcode = coo_matrix_expand_symmetric_storage ((struct coo_matrix_t*) A->repr);
  else if (A->format == CSR)
    errcode = csr_matrix_expand_symmetric_storage ((struct csr_matrix_t*) A->repr);
  else if (A->format == BCOO)
    {
      bebop_log (0, "*** sparse_matrix_expand_symmetric_storage: BCOO "
		"format not yet supported! ***\n");
      return -1;
    }
  else if (A->format == BCSR)
    {
      bebop_log (0, "*** sparse_matrix_expand_symmetric_storage: BCSR "
	       "format not yet supported! ***\n");
      return -1;
    }
  else if (A->format == JAD)
    {
      bebop_log (0, "*** sparse_matrix_expand_symmetric_storage: JAD "
	       "format not yet supported! ***\n");
      return -1;
    }
  else 
    {
      bebop_log (0, "*** sparse_matrix_expand_symmetric_storage: Invalid "
	       "matrix format %d ***\n", A->format);
      return -1;
    }

  return errcode;
}


int
sparse_matrix_convert (struct sparse_matrix_t* A, 
		       enum sparse_matrix_storage_format_t output_format)
{
  if (A->format == output_format)
    return 0; /* nothing to do */
  else if (A->format == BCOO)
    {
      if (output_format == BCSR)
	{
	  struct bcsr_matrix_t* B = bcoo_to_bcsr (A->repr);
	  if (B == NULL)
	    return -1;

	  destroy_bcoo_matrix (A->repr);
	  A->repr = B;
	  A->format = BCSR;
	  return 0;
	}
      else 
	{
	  bebop_log (0, "*** sparse_matrix_convert: unsupported conversion! ***\n");
	  return -1;
	}
    }
  else if (A->format == BCSR)
    {
      bebop_log (0, "*** sparse_matrix_convert: Conversions from BCSR to "
	       "other formats are not currently supported! ***\n");
      return -1;
    }
  else if (A->format == COO)
    {
      if (output_format == CSC)
	{
	  struct csc_matrix_t* B = coo_to_csc (A->repr);
	  if (B == NULL)
	    return -1;

	  destroy_coo_matrix (A->repr);
	  A->repr = B;
	  A->format = CSC;
	  return 0;
	}
      else if (output_format == CSR)
	{
	  struct csr_matrix_t* B = coo_to_csr (A->repr);
	  if (B == NULL)
	    return -1;

	  destroy_coo_matrix (A->repr);
	  A->repr = B;
	  A->format = CSR;
	  return 0;
	}
      else 
	{
	  bebop_log (0, "*** sparse_matrix_convert: unsupported conversion! ***\n");
	  return -1;
	}
    }
  else if (A->format == CSC)
    {
      if (output_format == CSR)
	{
	  struct csr_matrix_t* B = csc_to_csr (A->repr);
	  if (B == NULL)
	    return -1;

	  destroy_csc_matrix (A->repr);
	  A->repr = B;
	  A->format = CSR;
	  return 0;
	}
      else if (output_format == COO)
	{
	  struct coo_matrix_t* B = csc_to_coo (A->repr);
	  if (B == NULL)
	    return -1;

	  destroy_csc_matrix (A->repr);
	  A->repr = B;
	  A->format = COO;
	  return 0;
	}
      else 
	{
	  bebop_log (0, "*** sparse_matrix_convert: unsupported conversion! ***\n");
	  return -1;
	}
    }
  else if (A->format == CSR)
    {
      if (output_format == CSC)
	{
	  struct csc_matrix_t* B = csr_to_csc (A->repr);
	  if (B == NULL)
	    return -1;

	  destroy_csr_matrix (A->repr);
	  A->repr = B;
	  A->format = CSC;
	  return 0;
	}
      else if (output_format == COO)
	{
	  struct coo_matrix_t* B = csr_to_coo (A->repr);
	  if (B == NULL)
	    return -1;

	  destroy_csr_matrix (A->repr);
	  A->repr = B;
	  A->format = COO;
	  return 0;
	}
      else if (output_format == JAD)
    	{
	  struct jad_matrix_t* B = csr_to_jad (A->repr);
	  if (B == NULL)
	    return -1;

	  destroy_csr_matrix (A->repr);
	  A->repr = B;
	  A->format = JAD;
	  return 0;
	}
      else 
	{
	  bebop_log (0, "*** sparse_matrix_convert: unsupported conversion! ***\n");
	  return -1;
	}
    }

  return 0;
}

int
valid_sparse_matrix (struct sparse_matrix_t* A)
{
  if (A == NULL)
    {
      bebop_log (0, "*** valid_sparse_matrix: A is NULL ***\n");
      return 0;
    }
  else if (A->format == COO)
    return valid_coo_matrix_p (A->repr);
  else if (A->format == CSR)
    return valid_csr_matrix_p (A->repr);
  else if (A->format == CSC)
    return valid_csc_matrix_p (A->repr);
  else if (A->format == BCOO)
    return valid_bcoo_matrix_p (A->repr);
  else if (A->format == BCSR)
    return valid_bcsr_matrix_p (A->repr);
  else if (A->format == JAD)
    {
      bebop_log (0, "*** WARNING: valid_sparse_matrix: JAD format "
	       "not yet supported! ***\n");
      return 0;
    }
  else
    {
      bebop_log (0, "*** valid_sparse_matrix: invalid internal "
	       "storage format %d ***\n", A->format);
      return 0;
    }
}



enum sparse_matrix_storage_format_t 
sparse_matrix_storage_format_string_to_enum (const char* const s)
{
  if (strncmp (s, "CSC", 3) == 0)
    return CSC;
  else if (strncmp (s, "CSR", 3) == 0)
    return CSR;
  else if (strncmp (s, "COO", 3) == 0)
    return COO;
  else if (strncmp (s, "BCOO", 4) == 0)
    return BCOO;
  else if (strncmp (s, "BCSR", 4) == 0)
    return BCSR;
  else if (strncmp (s, "JAD", 3) == 0)
    return JAD;
  else 
    {
      bebop_log (0, "*** sparse_matrix_storage_format_string_to_enum: "
	       "invalid storage format %s ***\n", s);
      return -1;
    }
}

enum sparse_matrix_file_format_t
sparse_matrix_file_format_string_to_enum (const char* const s)
{
  if (strncmp (s, "HARWELL_BOEING", 14) == 0)
    return HARWELL_BOEING;
  else if (strncmp (s, "MATRIX_MARKET", 13) == 0)
    return MATRIX_MARKET;
  else if (strncmp (s, "MATLAB", 6) == 0)
    return MATLAB;
  else 
    {
      bebop_log (0, "*** sparse_matrix_file_format_string_to_enum: "
	       "invalid file format %s ***\n", s);
      return -1;
    }
}
