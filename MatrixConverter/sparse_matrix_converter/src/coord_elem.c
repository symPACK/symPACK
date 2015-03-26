/**
 * @file coord_elem.c
 * @author Mark Hoemmen
 * @since 26 May 2005
 * @date Time-stamp: <2008-07-16 11:14:53 mhoemmen>
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
#include <bebop/smc/coord_elem.h>

#include <stdio.h>


void
print_coord_array_matrix (FILE* out, const struct coord_elem_t A[], const int nnz)
{
  int i;
  int max_col_idx = -1, max_row_idx = -1;

  /**
   * \bug (mfh 29 June 2004) Should really store m and n in the coordinate 
   * array data structure; otherwise it is necessary to have an explicit 
   * zero at (m-1,n-1) in order to get m and n.
   */
  /*
   * Find the number of rows and columns in the matrix, by finding the maximum
   * row and column indices.
   */
  if (nnz > 0)
    {
      for (i = 0; i < nnz; i++)
	{
	  const struct coord_elem_t curelem = A[i];

	  if (curelem.r > max_row_idx) max_row_idx = curelem.r;
	  if (curelem.r > max_col_idx) max_col_idx = curelem.c;
	}
    }

  fprintf (out, "%%%%MatrixMarket matrix coordinate real general\n");
  fprintf (out, "%d %d %d\n", max_row_idx + 1, max_col_idx + 1, nnz);

  for (i = 0; i < nnz; i++)
    fprintf (out, "%d %d %g\n", A[i].r, A[i].c, A[i].val);
}


