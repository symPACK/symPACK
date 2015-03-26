/**
 * @file csr_spmv.c
 * @author Mark Hoemmen
 * @since Jun 2006
 * @date Time-stamp: <2008-07-16 11:16:19 mhoemmen>
 *
 * Reference implementation of CSR SpMV.
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
#include <bebop/smc/csr_spmv.h>
#include <assert.h>


void
csr_spmv_transpose_double_real (double* y, const int m, const double beta, 
				const double alpha, const double* val, 
				const int* ind, const int* ptr, double* x, 
				const int n)
{
  int j, k;

  if (beta == 0.0)
    {
      for (j = 0; j < n; j++)
	y[j] = 0.0;
    }
  else 
    {
      /* Scale y by beta */
      for (j = 0; j < n; j++)
	y[j] = beta * y[j];
    }

  /* Test for easy exit */
  if (alpha == 0.0)
    return;

  /* for each column j of A^T, do a daxpy / scatter */
  for (j = 0; j < m; j++)
    {
      const double alpha_times_xj = alpha * x[j];
      const int start = ptr[j];
      const int end = ptr[j + 1];
      for (k = start; k < end; k++)
	{
	  const int row = ind[k];
	  y[row] = y[row] + val[k] * alpha_times_xj;
	}
    }
}

void
csr_spmv_transpose_double_complex (double_Complex* y, const int m, const double_Complex beta, 
				   const double_Complex alpha, const double_Complex* val, 
				   const int* ind, const int* ptr, double_Complex* x, 
				   const int n)
{
  int j, k;

  if (double_Complex_equal (beta, double_Complex_ZERO))
    {
      for (j = 0; j < n; j++)
	y[j] = double_Complex_ZERO;
    }
  else
    {
      /* Scale y by beta */
      for (j = 0; j < n; j++)
	y[j] = double_Complex_multiply (beta, y[j]);
    }

  /* Test for easy exit */
  if (double_Complex_equal (alpha, double_Complex_ZERO))
    return;

  /* for each column j of A^T, do a gather, daxpy and scatter */
  for (j = 0; j < m; j++)
    {
      const double_Complex alpha_times_xj = double_Complex_multiply (alpha, x[j]);
      const int start = ptr[j];
      const int end = ptr[j + 1];
      for (k = start; k < end; k++)
	{
	  const int row = ind[k];
	  y[row] = double_Complex_add (y[row], double_Complex_multiply (val[k], alpha_times_xj));
	}
    }
}


/**
 * y := alpha*A*x + beta*y
 */
void
csr_spmv_double_real (double* y, const int m, const double beta, 
		      const double alpha, const double* val, 
		      const int* ind, const int* ptr, double* x, 
		      const int n)
{
  int i, j, start, end;

  if (alpha == 0.0)
    {
      if (beta != 0.0)
	{
	  for (i = 0; i < m; i++)
	    y[i] = beta * y[i];
	}
      else 
	{
	  for (i = 0; i < m; i++)
	    y[i] = 0.0;
	}

      return;
    }

  if (beta != 0.0)
    {
      start = ptr[0];
      for (i = 0; i < m; i++)
	{
	  register double yi = beta * y[i];
	  end = ptr[i+1]; /* start is already set to ptr[i] from the previous iteration */
	  for (j = start; j < end; j++)
	    {
	      const int col = ind[j];
	      const double a = val[j];
	      yi += alpha * a * x[col];
	    }
	  y[i] = yi;
	  start = end;
	}
    }
  else
    {
      start = ptr[0];
      for (i = 0; i < m; i++)
	{
	  register double yi = 0.0;
	  end = ptr[i+1]; /* start is already set to ptr[i] from the previous iteration */
	  for (j = start; j < end; j++)
	    {
	      const int col = ind[j];
	      const double a = val[j];
	      yi += alpha * a * x[col];
	    }
	  y[i] = yi;
	  start = end;
	}
    }
}


void
csr_spmv_double_complex (double_Complex* y, const int m, const double_Complex beta, 
			 const double_Complex alpha, const double_Complex* val, 
			 const int* ind, const int* ptr, double_Complex* x, 
			 const int n)
{
  int i, j, start, end;

  if (double_Complex_equal (alpha, double_Complex_ZERO))
    {
      for (i = 0; i < m; i++)
	y[i] = double_Complex_multiply (beta, y[i]);

      return;
    }

  start = ptr[0];
  for (i = 0; i < m; i++)
    {
      double_Complex yi = double_Complex_multiply (beta, y[i]);
      end = ptr[i+1];
      for (j = start; j < end; j++)
	{
	  /* yi = yi + alpha * val[j] * x[ind[j]]; */
	  yi = double_Complex_add (yi, double_Complex_multiply (alpha, double_Complex_multiply (val[j], x[ind[j]])));
	}
      y[i] = yi;
      start = end;
    }
}
