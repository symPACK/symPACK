/**
 * @file mfh_matlib.c
 * @author Mark Hoemmen
 * @since 06/08/04 11:16:16 PDT
 * @date Time-stamp: <2008-07-16 11:24:11 mhoemmen>
 * @version 1.2
 *
 * Version 1.1 (mfh 1 July 2004): Moved csc_matrix_t struct definition and 
 * member functions into their own files.
 * Version 1.2 (mfh 26 May 2005): Moved coo_matrix_t and coord_elem_t into
 * their own files.
 * 
 * Adapted by mfh 25 Mar 2005.
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

#include <bebop/smc/mfh_matlib.h>
#include <bebop/smc/csc_matrix.h>

#include <bebop/util/avltree_intpair.h>
#include <bebop/util/malloc.h>
#include <bebop/util/util.h>
#include <bebop/util/log.h>

#include <string.h> /* strcmp */




/*======================================================================*/
void
ones (double* A, const int m, const int n)
{
  const int L = m*n;
  int i;
  for (i = 0; i < L; i++)
	A[i] = 1.0;
}


/*======================================================================*/
void
zeros (double* A, const int m, const int n)
{
  const int L = m*n;
  int i;
  for (i = 0; i < L; i++)
	A[i] = 0.0;
}


/*======================================================================*/
void
unit_vector (double* v, const int n, const int which)
{
  int i;
  if (which < 0 || which >= n)
	{
	  bebop_log (0, "*** unit_vector: which=%d out of valid range "
		    "[0,%d) ***\n", which, n);
	  bebop_exit (EXIT_FAILURE);
	}
  for (i = 0; i < n; i++) v[i] = 0.0;

  v[which] = 1.0;
}


/*======================================================================*/
void
print_dense_vector (FILE* out, const double* x, const int length)
{
  int i;
  for (i = 0; i < length; i++)
    fprintf (out, "%E\n", x[i]);
}


/*======================================================================*/
void
print_dense_matrix (FILE* out, const double* A, const int m, const int n, 
					const char* trans)
{
  int i, j;

  bebop_log (2, "=== print_dense_matrix ===");

  if (0 == strcmp (trans, "T"))
    {
      /* 
       * Lazy solution:  transpose the matrix and print the transposed 
       * version.  It's a human-readable-print routine, so the matrix 
       * probably isn't very big anyway.
       */
      double* B = transpose_dense_col_oriented_matrix (A, m, n);
      print_dense_matrix (out, B, n, m, "N");
      bebop_free (B);
    }
  else
    {
      for (i = 0; i < m; i++) /* for each row i */
	{
	  fprintf (out, "[");
	  if (n > 0)
	    {
	      for (j = 0; j < n - 1; j++) /* for each entry in the row */
		fprintf (out, "%E ", A[i + m*j]);

	      fprintf (out, "%E", A[i + m*j]);
	    }
	  fprintf (out, "]\n");
	}
    }

  bebop_log (2,  "=== Done with print_dense_matrix ===");
}


/*======================================================================*/
double* 
transpose_dense_col_oriented_matrix (const double* A, 
				     const int m, const int n)
{
  double* B = bebop_malloc (n * m * sizeof (double));
  int i, j;

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      B[j + n*i] = A[i + m*j];

  return B;
}


