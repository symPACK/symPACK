/**
 * @file dense_to_csr.c
 * @author Mark Hoemmen <mhoemmen@cs.berkeley.edu>
 * @since 04 Apr 2008
 * @date Time-stamp: <2008-07-16 11:22:35 mhoemmen>
 * 
 * Copyright (c) 2008, Regents of the University of California All
 * rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the
 * distribution.
 * 
 * * Neither the name of the University of California, Berkeley, nor
 * the names of its contributors may be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
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
#include <bebop/smc/dense_to_csr.h>
#include <bebop/util/log.h>
#include <stdio.h>

$IndType
dense_to_csr_unsym_count (const $ValType const A[],
			  const $IndType col_oriented, 
			  const $IndType m, 
			  const $IndType n, 
			  const $IndType pitch)

{
  $IndType count = 0;
  $IndType i, j;

  if (col_oriented)
    {
      for (j = 0; j < n; j++)
	for (i = 0; i < m; i++)
	  if (A[i + j * pitch] != ($ValType) 0.0)
	    count++;
    }
  else
    {
      for (i = 0; i < m; i++)
	for (j = 0; j < n; j++)
	  if (A[j + i * pitch] != ($ValType) 0.0)
	    count++;
    }
  return count;
}


int
dense_to_csr_unsym ($ValType A[], 
		    const $IndType col_oriented, 
		    const $IndType m, 
		    const $IndType m, 
		    const $IndType lda, 
		    const $ValType val[], 
		    const $IndType ind[], 
		    const $IndType ptr[],
		    const $IndType nnz) 
{
  $IndType cur = 0;
  $IndType i, j;

  if (col_oriented)
    {
      for (i = 0, cur = 0; i < m; i++)
	{
	  ptr[i] = cur;
	  for (j = 0; j < n; j++)
	    {
	      if (A[i + j * pitch] != ($ValType) (0))
		{
		  ind[cur] = j;
		  val[cur] = A[i + j * pitch];
		  cur++;
		}
	    }
	}
    }
  else
    {
      for (i = 0, cur = 0; i < m; i++)
	{
	  ptr[i] = cur;
	  for (j = 0; j < n; j++)
	    {
	      if (A[j + i * pitch] != ($ValType) (0))
		{
		  ind[cur] = j;
		  val[cur] = A[i + j * pitch];
		  cur++;
		}
	    }
	}
    }

  if (cur != nnz)
    {
      bebop_log (0, "There are %d nonzero entries in the input matrix,"
		 " but you said nnz = %d", (int) cur, (int) nnz);
      return -1;
    }
  ptr[m] = nnz;
  return 0;
}


