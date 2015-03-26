/**
 * @file bcoo_io.c
 * @author Mark Hoemmen
 * @since 07 Jul 2006
 * @date Time-stamp: <2008-07-16 11:13:56 mhoemmen>
 *
 * I/O routines taken from the original bcoo_matrix.c.
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

#include <bebop/util/complex.h>
#include <bebop/util/util.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 


int
bcoo_print_mm (FILE* out, const int nnzb, const int r, const int c, const int bm, const int bn, enum symmetry_type_t symmetry_type, enum value_type_t value_type, const int index_base, const int col_oriented_p, const int II[], const int JJ[], const void* __val)
{
  const int block_size = r * c;
  int i;
  char symmetry_type_label[20];
  char value_type_label[20];
  
  if (symmetry_type == UNSYMMETRIC)
    strncpy (symmetry_type_label, "general", 19);
  else if (symmetry_type == SYMMETRIC)
    strncpy (symmetry_type_label, "symmetric", 19);
  else if (symmetry_type == SKEW_SYMMETRIC)
    strncpy (symmetry_type_label, "skew-symmetric", 19);
  else if (symmetry_type == HERMITIAN)
    strncpy (symmetry_type_label, "hermitian", 19);
  else 
    {
      fprintf (stderr, "*** bcoo_print_mm: "
	       "Invalid symmetry type %d of the given BCOO format sparse "
	       "matrix! ***\n", A->symmetry_type);
      return -1;
    }

  if (value_type == REAL)
    strncpy (value_type_label, "real", 19);
  else if (value_type == COMPLEX)
    strncpy (value_type_label, "complex", 19);
  else if (value_type == PATTERN)
    strncpy (value_type_label, "pattern", 19);
  else
    {
      fprintf (stderr, "*** bcoo_print_mm: Unsupported value type %d ***\n", value_type);
      return -1;
    }

  fprintf (out, "%%%%MatrixMarket matrix coordinate %s %s\n", 
	   value_type_label, symmetry_type_label);
  fprintf (out, "%d %d %d\n", bm * r, bn * c, nnzb * block_size);

  if (value_type == REAL)
    {
      double* val = (double*) (__val);

      /* MatrixMarket format uses 1-based indices, so we have to convert
       * to 1-based when we print out the matrix.  */
      for (i = 0; i < nnzb; i++, val += block_size)
	{
	  int ii, jj;
	  const int rowindexbase = II[i] + (1 - index_base);
	  const int colindexbase = JJ[i] + (1 - index_base);

	  if (col_oriented_p)
	    {
	      for (ii = 0; ii < r; ii++)
		for (jj = 0; jj < c; jj++)
		  fprintf (out, "%d %d %.13e\n", rowindexbase + ii, 
			   colindexbase + jj, val[jj + ii*c]);
	    }
	  else 
	    {
	      for (ii = 0; ii < r; ii++)
		for (jj = 0; jj < c; jj++)
		  fprintf (out, "%d %d %.13e\n", rowindexbase + ii, 
			   colindexbase + jj, val[ii + jj*r]);
	    }
	}
    }
  else if (value_type == COMPLEX)
    {
      double_Complex* val = (double_Complex*) (__val);

      /* MatrixMarket format uses 1-based indices, so we have to convert
       * to 1-based when we print out the matrix.  */
      for (i = 0; i < nnzb; i++, val += block_size)
	{
	  int ii, jj;
	  const int rowindexbase = II[i] + (1 - index_base);
	  const int colindexbase = JJ[i] + (1 - index_base);

	  if (col_oriented_p)
	    {
	      for (ii = 0; ii < r; ii++)
		for (jj = 0; jj < c; jj++)
		  fprintf (out, "%d %d %.13e %.13e\n", rowindexbase + ii, 
			   colindexbase + jj, 
			   double_Complex_real_part (val[jj + ii*c]),
			   double_Complex_imag_part (val[jj + ii*c]));
	    }
	  else 
	    {
	      for (ii = 0; ii < r; ii++)
		for (jj = 0; jj < c; jj++)
		  fprintf (out, "%d %d %.13e %.13e\n", rowindexbase + ii, 
			   colindexbase + jj, 
			   double_Complex_real_part (val[block_size*i + ii + jj*r]),
			   double_Complex_imag_part (val[block_size*i + ii + jj*r]));
	    }
	}
    }
  else if (value_type == PATTERN)
    {
      /* MatrixMarket format uses 1-based indices, so we have to convert
       * to 1-based when we print out the matrix.  */
      for (i = 0; i < nnzb; i++)
	{
	  int ii, jj;
	  const int rowindexbase = II[i] + (1 - index_base);
	  const int colindexbase = JJ[i] + (1 - index_base);

	  /* It doesn't matter if the blocks are row-oriented or 
	   * column-oriented, since there are no values. */

	  for (ii = 0; ii < r; ii++)
	    for (jj = 0; jj < c; jj++)
	      fprintf (out, "%d %d\n", rowindexbase + ii, colindexbase + jj);
	}
    }

  return 0;
}



int
bcoo_save_mm (const char* const filename, const int nnzb, const int r, const int c, const int bm, const int bn, enum symmetry_type_t symmetry_type, enum value_type_t value_type, const int index_base, const int col_oriented_p, const int II[], const int JJ[], const void* __val)
{
  FILE* out = NULL;
  int errcode = 0;

  out = fopen (filename, "w");
  if (out == NULL)
    {
      fprintf (stderr, "*** bcoo_save_mm: failed to open output file \"%s\" ***\n", filename);
      return -1;
    }

  errcode = bcoo_print_mm (out, nnzb, r, c, bm, cn, symmetry_type, 
			   value_type, index_base, col_oriented_p, 
			   II, JJ, val);
  if (0 != fclose (out))
    {
      fprintf (stderr, "*** bcoo_save_mm: failed "
	       "to close output file \"%s\" ***\n", filename);
      return -1;
    }
  return errcode;
}





