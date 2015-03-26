/**
 * @file csr_matmatmult.c
 * @author Mark Hoemmen <mhoemmen@cs.berkeley.edu>
 * @date 20 Jul 2008
 * 
 * Untuned implementation of sparse matrix-matrix multiplication.
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
#include <bebop/smc/csr_matmatmult.h>

#include <bebop/util/malloc.h>
#include <bebop/util/util.h>
#include <bebop/util/log.h>

#include <assert.h>
#include <stdio.h>
#include <math.h>



int
csr_matmatmult_double_real (int** pCptr, int** pCind, double** pCval, 
			    int* pCnnz, double alpha, int* Aptr, 
			    int* Aind, double* Aval, int* Bptr, 
			    int* Bind, double* Bval, 
			    const int m, const int p, const int n)
{
  /* Borrowed from the Hypre code hypre_CSRMatrixMultiply */
  int ia, ib, ic, ja, jb, num_nonzeros=0;
  int row_start, counter;
  double a_entry, b_entry;
  int* B_marker;
  int* Cptr;
  int* Cind;
  double* Cval;

  bebop_log (2, "=== csr_matmatmult_double_real ===\n");

  bebop_log (3, "\tAllocating B_marker and Cptr\n");
  B_marker = bebop_calloc (n, sizeof (int));
  Cptr = bebop_malloc ((m+1) * sizeof (int));
  for (ic = 1; ic < m; ic++)
    Cptr[ic] = -1; /* flag to detect errors */

  Cptr[0] = 0;
  
  for (ib = 0; ib < n; ib++)
    B_marker[ib] = -1;

  bebop_log (3, "\tCounting nonzeros in C and setting up Cptr\n");
  /* 
   * Count the number of nonzeros in each row of C, and use this
   * to set up the ptr array of C.
   */
  for (ic = 0; ic < m; ic++) /* for each row of C */
    {
      /*bebop_log (3, "\t\tRow %d of C\n", ic);*/
      /* For each row ic of A:
       *   For each entry (ic,ja) in row ic of A:
       *     Look in row ja of B (since entry (ic,ja) of A weights row ja of B):
       *     For each entry (ja,jb) in row ja of B:
       *       If we haven't seen an entry in column jb of B for the current row ic of A
       *         Mark it and increment nnz
       *       EndIf
       *     EndFor
       *   EndFor
       * EndFor
       */
      for (ia = Aptr[ic]; ia < Aptr[ic+1]; ia++) 
	{
	  ja = Aind[ia];
	  if (ja < 0 || ja >= p)
	    {
	      bebop_log (0, "*** csr_matmatmult_double_real: Element A(%d,%d) is "
			"out of range [0,%d] x [0,%d] ***\n", ic, ja, m-1, p-1);
	      bebop_free (Cptr);
	      bebop_free (B_marker);
	      return -1;
	    }
	  for (ib = Bptr[ja]; ib < Bptr[ja+1]; ib++)
	    {
	      jb = Bind[ib];
	      if (jb < 0 || jb >= n)
		{
		  bebop_log (0, "*** csr_matmatmult_double_real: Index jb = %d "
			    "is out of range [0,%d] ***\n", jb, n-1);
		  bebop_free (Cptr);
		  bebop_free (B_marker);
		  return -1;
		}
	      if (B_marker[jb] != ic)
		{
		  B_marker[jb] = ic;
		  num_nonzeros++;
		  /*bebop_log (4, "\t\t\tIdentified new nonzero (%d,%d)\n", ic, jb);*/
		}
	    }
	}
      Cptr[ic+1] = num_nonzeros;            
    }
  bebop_log (3, "\tThere are %d nonzeros in C.\n", num_nonzeros);
  
  bebop_log (3, "\tAllocating Cval and Cind\n");
  Cval = bebop_malloc (num_nonzeros * sizeof (double));
  Cind = bebop_malloc (num_nonzeros * sizeof (int));

  for (ic = 0; ic < num_nonzeros; ic++)
    {
      /* 
       * We set Cval[ic] and Cind[ic] to "invalid" values (NaN and -1, 
       * respectively) to make sure that we fill them all in.   Note 
       * that we can't use 0.0/0.0 to generate NaN because when running
       * this function via an FFI under some Lisps (SBCL in particular),
       * 0.0/0.0 raises a floating-point exception instead of being a 
       * silent NaN like we want it to be (mfh 28 July 2006).  That's 
       * why we need the C99 #defined value NAN and the fpclassify()
       * function in <math.h>.  We use the assert to make sure that the 
       * NAN really is a NAN.
       */
#ifdef USE_ISNAN
      Cval[ic] = NAN;
      assert (fpclassify (Cval[ic]) == FP_NAN);
#else
      /* 
       * We don't have C99 so we can't use NAN; initialize Cval[ic]
       * to some innocuous value instead.  Note that we can't just set
       * it to 0.0/0.0, as this may be evaluated; under plain C this 
       * typically doesn't signal, but under Lisp running the C function,
       * it may signal.
       */
      Cval[ic] = 0.0;
#endif
      Cind[ic] = -1;
    }

  for (ib = 0; ib < n; ib++)
    B_marker[ib] = -1;

  bebop_log (3, "\tComputing nonzero values in C\n");
  counter = 0;
  for (ic = 0; ic < m; ic++)
    {
      /*bebop_log (3, "\t\tRow %d of C\n", ic);*/
      row_start = Cptr[ic];
      /*bebop_log (3, "\t\trow_start = %d\n", row_start);*/
      for (ia = Aptr[ic]; ia < Aptr[ic+1]; ia++)
        {
	  ja = Aind[ia];
	  a_entry = Aval[ia];
	  /*bebop_log (4, "\t\t\tA(%d,%d) = %e\n", ic, ja, a_entry);*/
	  for (ib = Bptr[ja]; ib < Bptr[ja+1]; ib++)
	    {
	      jb = Bind[ib];
	      b_entry = Bval[ib];
	      /*bebop_log (4, "\t\t\t\tB(%d,%d) = %e\n", ja, jb, b_entry);*/
	      if (B_marker[jb] < row_start)
		{
		  B_marker[jb] = counter;
		  Cind[B_marker[jb]] = jb;
		  /*
		  bebop_log (4, "\t\t\t\tSetting Cval[%d] = %e\n", 
			    B_marker[jb], a_entry * b_entry);
		  */
		  Cval[B_marker[jb]] = a_entry * b_entry;
		  counter++;
		}
	      else
		{
		  Cval[B_marker[jb]] += a_entry * b_entry;
		  /*
		  bebop_log (4, "\t\t\t\tUpdated Cval[%d] to %e\n",
			    B_marker[jb], Cval[B_marker[jb]]);
	          */
		}
	    }
        }
    }
  bebop_free (B_marker);

  if (bebop_debug_level () > 0)
    {
      int failure = 0;
#ifdef HAVE_ISNAN
      int found_nan_in_A = 0;
      int found_nan_in_B = 0;
      int found_nan_in_C = 0;
#endif /* HAVE_ISNAN */
      int num_uninitialized_colindices = 0;

#ifdef HAVE_ISNAN
      bebop_log (3, "\tLooking for uninitialized values in Cval and Cind\n");
#else
      bebop_log (3, "\tLooking for uninitialized values in Cind\n");
#endif /* HAVE_ISNAN */

      for (ic = 0; ic < num_nonzeros; ic++)
	{
#ifdef HAVE_ISNAN
	  /* Note: Maybe Cval is supposed to contain a NaN, so don't report 
	   * failure if you find one.  Just log it. */
	  if (isnan (Cval[ic]))
	    {
	      found_nan_in_C = 1;
	    }
#endif /* HAVE_ISNAN */
	  if (-1 == Cind[ic])
	    {
	      failure = 1;
	      num_uninitialized_colindices++;
	    }
	}
      if (num_uninitialized_colindices > 0)
	bebop_log (0, "*** csr_matmatmult_double_real: C contains %d "
		  "uninitialized column indices, suggesting a bug in "
		  "CSR matrix-matrix multiplication! ***\n", 
		  num_uninitialized_colindices);
#ifdef HAVE_ISNAN
      if (found_nan_in_C)
	{
	  /* 
	   * Search the two input matrices for NaNs.  If they don't have any, 
	   * then there must be a bug, since C then shouldn't have any NaNs.
	   */
	  for (ic = 0; ic < Aptr[m]; ic++)
	    {
	      if (isnan (Aval[ic]))
		found_nan_in_A = 1;
	    }
	  for (ic = 0; ic < Bptr[p]; ic++)
	    {
	      if (isnan (Bval[ic]))
		found_nan_in_B = 1;
	    }
	  if (! found_nan_in_A && ! found_nan_in_B)
	    {
	      bebop_log (0, "*** csr_matmatmult_double_real: C contains NaNs "
			"but A and B (the two input matrices) do not, which "
			"suggests a bug in CSR matrix-matrix multiplication "
			"that left values of C uninitialized.  Before doing "
			"matrix-matrix multiplication, we fill Cval with "
			"NaNs to indicate uninitialized values, hence the "
			"detection of NaNs in Cval after matrix-matrix "
			"multiplication indicates that some of Cval wasn\'t "
			"filled in properly. ***\n");
	      failure = 1;
	    }
	  else
	    {
	      if (found_nan_in_A)
		{
		  if (found_nan_in_B)
		    bebop_log (0, "*** csr_matmatmult_double_real: we found "
			      "some NaNs in both the first and second input "
			      "matrices.  You may want to check your prior "
			      "computations. ***\n");
		  else
		    bebop_log (0, "*** csr_matmatmult_double_real: we found "
			      "some NaNs in the first input matrix but not "
			      "the second one.  You may want to check your "
			      "prior computations. ***\n");
		}
	      else if (found_nan_in_B)
		bebop_log (0, "*** csr_matmatmult_double_real: we found some "
			  "NaNs in the second input matrix but not the first "
			  "one.  You may want to check your prior computation"
			  "s. ***\n");
	    }
	}
#endif /* HAVE_ISNAN */
      if (failure)
	{
	  bebop_log (0, "*** csr_matmatmult_double_real: aborting due to fail"
		    "ures (which should be discussed in the error log above)"
		    " in matrix-matrix multiplication ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
    }

  if (bebop_debug_level () > 9)
    {
      bebop_log (10, "Cval:\n");
      for (ic = 0; ic < num_nonzeros; ic++)
        bebop_log (10, "%e\n", Cval[ic]);
    }

  /* Go back through and scale by alpha */
  bebop_log (3, "\tScaling by alpha\n");
  for (ic = 0; ic < m; ic++)
    {
      int k;
      for (k = Cptr[ic]; k < Cptr[ic+1]; k++)
	Cval[k] = alpha * Cval[k];
    }

  *pCnnz = num_nonzeros;
  *pCval = Cval;
  *pCptr = Cptr;
  *pCind = Cind;
  bebop_log (2, "=== Done with csr_matmatmult_double_real ===\n");
  return 0;
}  


int
csr_matmatmult_double_complex (int** pCptr, int** pCind, 
			       double_Complex** pCval, 
			       int* pCnnz, double_Complex alpha, int* Aptr, 
			       int* Aind, double_Complex* Aval, int* Bptr, 
			       int* Bind, double_Complex* Bval, 
			       const int m, const int p, const int n)
{
  /* Borrowed from the Hypre code hypre_CSRMatrixMultiply */
  int ia, ib, ic, ja, jb, num_nonzeros=0;
  int row_start, counter;
  double_Complex a_entry, b_entry;
  int* B_marker;
  int* Cptr;
  int* Cind;
  double_Complex* Cval;

  bebop_log (2, "=== csr_matmatmult_double_complex ===\n");

  bebop_log (3, "\tAllocating B_marker and Cptr\n");
  B_marker = bebop_calloc (n, sizeof (int));
  Cptr = bebop_malloc ((m+1) * sizeof (int));
  for (ic = 1; ic < m; ic++)
    Cptr[ic] = -1; /* flag to detect errors */

  Cptr[0] = 0;
  
  for (ib = 0; ib < n; ib++)
    B_marker[ib] = -1;

  bebop_log (3, "\tCounting nonzeros in C and setting up Cptr\n");
  /* 
   * Count the number of nonzeros in each row of C, and use this
   * to set up the ptr array of C.
   */
  for (ic = 0; ic < m; ic++) /* for each row of C */
    {
      /*bebop_log (3, "\t\tRow %d of C\n", ic);*/
      /* For each row ic of A:
       *   For each entry (ic,ja) in row ic of A:
       *     Look in row ja of B (since entry (ic,ja) of A weights row ja of B):
       *     For each entry (ja,jb) in row ja of B:
       *       If we haven't seen an entry in column jb of B for the current row ic of A
       *         Mark it and increment nnz
       *       EndIf
       *     EndFor
       *   EndFor
       * EndFor
       */
      for (ia = Aptr[ic]; ia < Aptr[ic+1]; ia++) 
	{
	  ja = Aind[ia];
	  if (ja < 0 || ja >= p)
	    {
	      bebop_log (0, "*** csr_matmatmult_double_complex: Element "
			"A(%d,%d) is out of range [0,%d] x [0,%d] ***\n", 
			ic, ja, m-1, p-1);
	      bebop_free (Cptr);
	      bebop_free (B_marker);
	      return -1;
	    }
	  for (ib = Bptr[ja]; ib < Bptr[ja+1]; ib++)
	    {
	      jb = Bind[ib];
	      if (jb < 0 || jb >= n)
		{
		  bebop_log (0, "*** csr_matmatmult_double_complex: Index "
			    "jb = %d is out of range [0,%d] ***\n", jb, n-1);
		  bebop_free (Cptr);
		  bebop_free (B_marker);
		  return -1;
		}
	      if (B_marker[jb] != ic)
		{
		  B_marker[jb] = ic;
		  num_nonzeros++;
		  /*bebop_log (4, "\t\t\tIdentified new nonzero (%d,%d)\n", ic, jb);*/
		}
	    }
	}
      Cptr[ic+1] = num_nonzeros;            
    }
  bebop_log (3, "\tThere are %d nonzeros in C.\n", num_nonzeros);
  
  bebop_log (3, "\tAllocating Cval and Cind\n");
  Cval = bebop_malloc (num_nonzeros * sizeof (double_Complex));
  Cind = bebop_malloc (num_nonzeros * sizeof (int));

  for (ic = 0; ic < num_nonzeros; ic++)
    {
      /* 
       * We set Cval[ic] and Cind[ic] to "invalid" values (NaN and -1, 
       * respectively) to make sure that we fill them all in.  NAN and
       * ilk only work when the compiler implements C99.
       */
#ifdef HAVE_ISNAN
      Cval[ic] = new_double_Complex (make_nan(), make_nan());
      assert (double_Complex_isnan (Cval[ic]));
#else
      Cval[ic] = new_double_Complex (0.0, 0.0);
#endif
      Cind[ic] = -1;
    }

  for (ib = 0; ib < n; ib++)
    B_marker[ib] = -1;

  bebop_log (3, "\tComputing nonzero values in C\n");
  counter = 0;
  for (ic = 0; ic < m; ic++)
    {
      /*bebop_log (3, "\t\tRow %d of C\n", ic);*/
      row_start = Cptr[ic];
      /*bebop_log (3, "\t\trow_start = %d\n", row_start);*/
      for (ia = Aptr[ic]; ia < Aptr[ic+1]; ia++)
        {
	  ja = Aind[ia];
	  a_entry = Aval[ia];
	  /*bebop_log (4, "\t\t\tA(%d,%d) = %e\n", ic, ja, a_entry);*/
	  for (ib = Bptr[ja]; ib < Bptr[ja+1]; ib++)
	    {
	      jb = Bind[ib];
	      b_entry = Bval[ib];
	      /*bebop_log (4, "\t\t\t\tB(%d,%d) = %e\n", ja, jb, b_entry);*/
	      if (B_marker[jb] < row_start)
		{
		  B_marker[jb] = counter;
		  Cind[B_marker[jb]] = jb;
		  /*
		  bebop_log (4, "\t\t\t\tSetting Cval[%d] = %e\n", 
			    B_marker[jb], a_entry * b_entry);
		  */
		  Cval[B_marker[jb]] = double_Complex_multiply (a_entry, b_entry);
		  counter++;
		}
	      else
		{
		  Cval[B_marker[jb]] = double_Complex_add (Cval[B_marker[jb]], double_Complex_multiply (a_entry, b_entry));
		  /*
		  bebop_log (4, "\t\t\t\tUpdated Cval[%d] to %e\n",
			    B_marker[jb], Cval[B_marker[jb]]);
	          */
		}
	    }
        }
    }
  bebop_free (B_marker);

  if (bebop_debug_level () > 0)
    {
      int failure = 0;
#ifdef HAVE_ISNAN
      int found_nan_in_A = 0;
      int found_nan_in_B = 0;
      int found_nan_in_C = 0;
#endif /* HAVE_ISNAN */
      int num_uninitialized_colindices = 0;

#ifdef HAVE_ISNAN
      bebop_log (3, "\tLooking for uninitialized values in Cval and Cind\n");
#else
      bebop_log (3, "\tLooking for uninitialized values in Cind\n");
#endif /* HAVE_ISNAN */

      for (ic = 0; ic < num_nonzeros; ic++)
	{
#ifdef HAVE_ISNAN
	  if (double_Complex_isnan (Cval[ic]))
	    {
	      found_nan_in_C = 1;
	    }
#endif /* HAVE_ISNAN */
	  if (-1 == Cind[ic])
	    {
	      failure = 1;
	      num_uninitialized_colindices++;
	    }
	}
      if (num_uninitialized_colindices > 0)
	bebop_log (0, "*** csr_matmatmult_double_complex: C contains %d "
		  "uninitialized column indices, suggesting a bug in "
		  "CSR matrix-matrix multiplication! ***\n", 
		  num_uninitialized_colindices);
#ifdef HAVE_ISNAN
      if (found_nan_in_C)
	{
	  /* 
	   * Search the two input matrices for NaNs.  If they don't have any, 
	   * then there must be a bug, since C then shouldn't have any NaNs.
	   */
	  for (ic = 0; ic < Aptr[m]; ic++)
	    {
	      if (double_Complex_isnan (Aval[ic]))
		found_nan_in_A = 1;
	    }
	  for (ic = 0; ic < Bptr[p]; ic++)
	    {
	      if (double_Complex_isnan (Bval[ic]))
		found_nan_in_B = 1;
	    }
	  if (! found_nan_in_A && ! found_nan_in_B)
	    {
	      bebop_log (0, "*** csr_matmatmult_double_complex: C contains NaNs "
			"but A and B (the two input matrices) do not, which "
			"suggests a bug in CSR matrix-matrix multiplication "
			"that left values of C uninitialized.  Before doing "
			"matrix-matrix multiplication, we fill Cval with "
			"NaNs to indicate uninitialized values, hence the "
			"detection of NaNs in Cval after matrix-matrix "
			"multiplication indicates that some of Cval wasn\'t "
			"filled in properly. ***\n");
	      failure = 1;
	    }
	  else
	    {
	      if (found_nan_in_A)
		{
		  if (found_nan_in_B)
		    bebop_log (0, "*** csr_matmatmult_double_complex: we found "
			      "some NaNs in both the first and second input "
			      "matrices.  You may want to check your prior "
			      "computations. ***\n");
		  else
		    bebop_log (0, "*** csr_matmatmult_double_complex: we found "
			      "some NaNs in the first input matrix but not "
			      "the second one.  You may want to check your "
			      "prior computations. ***\n");
		}
	      else if (found_nan_in_B)
		bebop_log (0, "*** csr_matmatmult_double_complex: we found some "
			  "NaNs in the second input matrix but not the first "
			  "one.  You may want to check your prior computation"
			  "s. ***\n");
	    }
	}
#endif /* HAVE_ISNAN */
      if (failure)
	{
	  bebop_log (0, "*** csr_matmatmult_double_complex: aborting due to fail"
		    "ures (which should be discussed in the error log above)"
		    " in matrix-matrix multiplication ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
    }

  if (bebop_debug_level () > 9)
    {
      bebop_log (10, "Cval:\n");
      for (ic = 0; ic < num_nonzeros; ic++)
        bebop_log (10, "%e+%eI\n", 
		  double_Complex_real_part (Cval[ic]), 
		  double_Complex_imag_part (Cval[ic]));
    }

  /* Go back through and scale by alpha */
  bebop_log (3, "\tScaling by alpha\n");
  for (ic = 0; ic < m; ic++)
    {
      int k;
      for (k = Cptr[ic]; k < Cptr[ic+1]; k++)
	Cval[k] = double_Complex_multiply (alpha, Cval[k]);
    }

  *pCnnz = num_nonzeros;
  *pCval = Cval;
  *pCptr = Cptr;
  *pCind = Cind;
  bebop_log (2, "=== Done with csr_matmatmult_double_complex ===\n");
  return 0;
}



int
csr_matmatmult_pattern (int** pCptr, int** pCind, int* pCnnz, int* Aptr, 
			int* Aind, int* Bptr, int* Bind, 
			const int m, const int p, const int n)
{
  /* Borrowed from the Hypre code hypre_CSRMatrixMultiply */
  int ia, ib, ic, ja, jb, num_nonzeros=0;
  int row_start, counter;
  int* B_marker;
  int* Cptr;
  int* Cind;

  bebop_log (2, "=== csr_matmatmult_pattern ===\n");

  bebop_log (3, "\tAllocating B_marker and Cptr\n");
  B_marker = bebop_calloc (n, sizeof (int));
  Cptr = bebop_malloc ((m+1) * sizeof (int));
  for (ic = 1; ic < m; ic++)
    Cptr[ic] = -1; /* flag to detect errors */

  Cptr[0] = 0;
  
  for (ib = 0; ib < n; ib++)
    B_marker[ib] = -1;

  bebop_log (3, "\tCounting nonzeros in C and setting up Cptr\n");
  /* 
   * Count the number of nonzeros in each row of C, and use this
   * to set up the ptr array of C.
   */
  for (ic = 0; ic < m; ic++) /* for each row of C */
    {
      /* For each row ic of A:
       *   For each entry (ic,ja) in row ic of A:
       *     Look in row ja of B (since entry (ic,ja) of A weights row ja of B):
       *     For each entry (ja,jb) in row ja of B:
       *       If we haven't seen an entry in column jb of B for the current row ic of A
       *         Mark it and increment nnz
       *       EndIf
       *     EndFor
       *   EndFor
       * EndFor
       */
      for (ia = Aptr[ic]; ia < Aptr[ic+1]; ia++) 
	{
	  ja = Aind[ia];
	  if (ja < 0 || ja >= p)
	    {
	      bebop_log (0, "*** csr_matmatmult_pattern: Element A(%d,%d) is "
			"out of range [0,%d] x [0,%d] ***\n", ic, ja, m-1, p-1);
	      bebop_free (Cptr);
	      bebop_free (B_marker);
	      return -1;
	    }
	  for (ib = Bptr[ja]; ib < Bptr[ja+1]; ib++)
	    {
	      jb = Bind[ib];
	      if (jb < 0 || jb >= n)
		{
		  bebop_log (0, "*** csr_matmatmult_pattern: Index jb = %d "
			    "is out of range [0,%d] ***\n", jb, n-1);
		  bebop_free (Cptr);
		  bebop_free (B_marker);
		  return -1;
		}
	      if (B_marker[jb] != ic)
		{
		  B_marker[jb] = ic;
		  num_nonzeros++;
		}
	    }
	}
      Cptr[ic+1] = num_nonzeros;            
    }
  bebop_log (3, "\tThere are %d nonzeros in C.\n", num_nonzeros);
  
  bebop_log (3, "\tAllocating Cind\n");
  Cind = bebop_malloc (num_nonzeros * sizeof (int));

  for (ic = 0; ic < num_nonzeros; ic++)
    {
      /* 
       * We set Cind[ic] to "invalid" values (-1) to make sure 
       * that we fill them all in.   
       */
      Cind[ic] = -1;
    }

  for (ib = 0; ib < n; ib++)
    B_marker[ib] = -1;

  bebop_log (3, "\tComputing nonzero values\' indices in C\n");
  counter = 0;
  for (ic = 0; ic < m; ic++)
    {
      row_start = Cptr[ic];
      for (ia = Aptr[ic]; ia < Aptr[ic+1]; ia++)
        {
	  ja = Aind[ia];
	  for (ib = Bptr[ja]; ib < Bptr[ja+1]; ib++)
	    {
	      jb = Bind[ib];
	      if (B_marker[jb] < row_start)
		{
		  B_marker[jb] = counter;
		  Cind[B_marker[jb]] = jb;
		  counter++;
		}
	    }
        }
    }
  bebop_free (B_marker);
  B_marker = NULL; /* "ground" the pointer */

  if (bebop_debug_level () > 0)
    {
      int failure = 0;
      int num_uninitialized_colindices = 0;

      bebop_log (3, "\tLooking for uninitialized values in Cind\n");
      for (ic = 0; ic < num_nonzeros; ic++)
	{
	  if (-1 == Cind[ic])
	    {
	      failure = 1;
	      num_uninitialized_colindices++;
	    }
	}
      if (num_uninitialized_colindices > 0)
	bebop_log (0, "*** csr_matmatmult_pattern: C contains %d "
		  "uninitialized column indices, suggesting a bug in "
		  "CSR matrix-matrix multiplication! ***\n", 
		  num_uninitialized_colindices);
      if (failure)
	{
	  bebop_log (0, "*** csr_matmatmult_pattern: aborting due to fail"
		    "ures (which should be discussed in the error log above)"
		    " in matrix-matrix multiplication ***\n");
	  bebop_exit (EXIT_FAILURE);
	}
    }

  *pCnnz = num_nonzeros;
  *pCptr = Cptr;
  *pCind = Cind;
  bebop_log (2, "=== Done with csr_matmatmult_pattern ===\n");
  return 0;
}  

      


#if 0
int
csr_matmatmult_double (int** pCptr, int** pCind, double** pCval, int* pCnnz,
		       double alpha,
		       int* Bptr, int* Bind, double* Bval,
		       int* Aptr, int* Aind, double* Aval,
		       const int m, const int p, const int n)
{
  int i, j, k, l, mm;
  int* Cptr, *Cind;
  double* Cval;

  int* marker;
  int* previous_index;
  int Cnnz = 0;

  return -1; /* BROKEN!!! */
  
  marker = (int*) bebop_malloc (n * sizeof (int));
  Cptr = (int*) bebop_malloc ((m+1) * sizeof (int));

  for (i = 0; i < n; i++)
    marker[i] = -1;
  
  Cptr[0] = Cnnz; /* = 0 */

  // for each row i of C
  for (i = 0; i < m; i++)
    {
      /* Go through the i-th row of B, using the index k.  This row
       * gives us the coefficients in the linear combination of the
       * rows of A.  Elements in the i-th row of B that are
       * (structurally) zero are skipped. */
      for (k = Bptr[i]; k < Bptr[i+1]; k++)
	{
	  /* Get the column index of the current coefficient
	   * (by which we will multiply the j-th row of A). */
	  j = Bind[k];
	  if (j < 0 || j >= p)
	    {
	      /* Column index out of range: clean up */
	      bebop_free (marker);
	      bebop_free (Cptr);
	      return -1; /* error flag */
	    }

	  /* Get the j-th row of A.  Any time we encounter a nonzero
	   * in this scan, there will be another nonzero in C. */
	  for (l = Aptr[j]; l < Aptr[j+1]; l++)
	    {
	      mm = Aind[l]; 
	      if (mm < 0 || mm >= n)
		{
		  /* Column index out of range: clean up */
		  bebop_free (marker);
		  bebop_free (Cptr);
		  return -2; /* error flag */
		}
	      /* Do this to avoid overcounting the number of nonzeros
	       * in C.  If we are currently computing row i of C, and
	       * we've already added C(i,mm) to C, then don't add it 
	       * again.  marker[mm] == i is true (on this i iteration)
	       * iff we have already added C(i,mm) to C. */
	      if (marker[mm] != i)
		{
		  marker[mm] = i; 
		  Cnnz++;
		}
	    }
	}
      Cptr[i+1] = Cnnz;
    }

  /* 
   * Now that we know how many nonzeros there will be in C, we can
   * allocate space for it. 
   */
  Cind = bebop_malloc (Cnnz * sizeof (int));
  Cval = bebop_calloc (Cnnz, sizeof (double));

  for (i = 0; i < n; i++)
    marker[i] = -1;

  previous_index = bebop_malloc (n * sizeof (int));
  for (i = 0; i < n; i++)
    previous_index[i] = -1;

  Cnnz = 0;
  for (i = 0; i < m; i++)
    {
      for (k = Bptr[i]; k < Bptr[i+1]; k++)
	{
	  const double b_coeff = Bval[k];
	  j = Bind[k];
	  for (l = Aptr[j]; l < Aptr[j+1]; l++)
	    {
	      mm = Aind[l];
	      if (marker[mm] != i)
		{
		  /* Remember that we've added a nonzero to C at (i,
		   * mm).  If while computing this i-th row of C, we
		   * encounter this column index mm again, then we
		   * don't need to increment Cnnz or store to Cind,
		   * because there is a structural nonzero already
		   * there.  (However, in that case we still need to
		   * increment the value of C(i,mm).  That's what the
		   * previous_index array is for.) */
		  marker[mm] = i;
		  /* Remember the location in the Cval array at which
		   * we stored this nonzero.  We'll need that location
		   * later, if while we are computing the i-th row of
		   * C, we encounter the same column index in a row of
		   * A. */
		  previous_index[mm] = Cnnz;
		  Cind[Cnnz] = mm;
		  Cval[Cnnz] = Cval[Cnnz] + alpha * b_coeff * Aval[l];
		  Cnnz++;
		}
	      else
		{
		  /* There's already a nonzero value in C at that
		   * location.  This means we need to add in the new
		   * nonzero value at that location, getting the
		   * location from the previous_index array. */
		  const int loc = previous_index[mm];
		  assert (loc != -1);
		  Cval[loc] = Cval[loc] + alpha * b_coeff * Aval[l];
		}
	    }
	}
    }

  /* Free up the temporary workspace */
  bebop_free (marker);
  bebop_free (previous_index);

  /* Return the matrix C through the output variables */
  *pCptr = Cptr;
  *pCind = Cind;
  *pCval = Cval;
  *pCnnz = Cnnz;

  /* No errors, so return zero */
  return 0;
}
#endif /* 0 */

#if 0
int
csr_matmatmult_complex (int** pCptr, int** pCind, 
			double_Complex** pCval, int* pCnnz,
			double_Complex alpha,
			int* Bptr, int* Bind, double_Complex* Bval,
			int* Aptr, int* Aind, double_Complex* Aval,
			const int m, const int p, const int n)
{
  int i, j, k, l, mm;
  int* Cptr, *Cind;
  double_Complex* Cval;

  int* marker;
  int* previous_index;
  int Cnnz = 0;

  return -1; /* BROKEN!!! */
  
  marker = (int*) bebop_malloc (n * sizeof (int));
  Cptr = (int*) bebop_malloc ((m+1) * sizeof (int));

  for (i = 0; i < n; i++)
    marker[i] = -1;
  
  Cptr[0] = Cnnz; /* = 0 */

  // for each row i of C
  for (i = 0; i < m; i++)
    {
      /* Go through the i-th row of B, using the index k.  This row
       * gives us the coefficients in the linear combination of the
       * rows of A.  Elements in the i-th row of B that are
       * (structurally) zero are skipped. */
      for (k = Bptr[i]; k < Bptr[i+1]; k++)
	{
	  /* Get the column index of the current coefficient
	   * (by which we will multiply the j-th row of A). */
	  j = Bind[k];
	  if (j < 0 || j >= p)
	    {
	      /* Column index out of range: clean up */
	      bebop_free (marker);
	      bebop_free (Cptr);
	      return -1; /* error flag */
	    }

	  /* Get the j-th row of A.  Any time we encounter a nonzero
	   * in this scan, there will be another nonzero in C. */
	  for (l = Aptr[j]; l < Aptr[j+1]; l++)
	    {
	      mm = Aind[l]; 
	      if (mm < 0 || mm >= n)
		{
		  /* Column index out of range: clean up */
		  bebop_free (marker);
		  bebop_free (Cptr);
		  return -2; /* error flag */
		}
	      /* Do this to avoid overcounting the number of nonzeros
	       * in C.  If we are currently computing row i of C, and
	       * we've already added C(i,mm) to C, then don't add it 
	       * again.  marker[mm] == i is true (on this i iteration)
	       * iff we have already added C(i,mm) to C. */
	      if (marker[mm] != i)
		{
		  marker[mm] = i; 
		  Cnnz++;
		}
	    }
	}
      Cptr[i+1] = Cnnz;
    }

  /* 
   * Now that we know how many nonzeros there will be in C, we can
   * allocate space for it. 
   */
  Cind = bebop_malloc (Cnnz * sizeof (int));
  Cval = bebop_calloc (Cnnz, sizeof (double_Complex));

  for (i = 0; i < n; i++)
    marker[i] = -1;

  previous_index = bebop_malloc (n * sizeof (int));
  for (i = 0; i < n; i++)
    previous_index[i] = -1;

  Cnnz = 0;
  for (i = 0; i < m; i++)
    {
      for (k = Bptr[i]; k < Bptr[i+1]; k++)
	{
	  const double_Complex b_coeff = Bval[k];
	  j = Bind[k];
	  for (l = Aptr[j]; l < Aptr[j+1]; l++)
	    {
	      mm = Aind[l];
	      if (marker[mm] != i)
		{
		  /* Remember that we've added a nonzero to C at (i,
		   * mm).  If while computing this i-th row of C, we
		   * encounter this column index mm again, then we
		   * don't need to increment Cnnz or store to Cind,
		   * because there is a structural nonzero already
		   * there.  (However, in that case we still need to
		   * increment the value of C(i,mm).  That's what the
		   * previous_index array is for.) */
		  marker[mm] = i;
		  /* Remember the location in the Cval array at which
		   * we stored this nonzero.  We'll need that location
		   * later, if while we are computing the i-th row of
		   * C, we encounter the same column index in a row of
		   * A. */
		  previous_index[mm] = Cnnz;
		  Cind[Cnnz] = mm;
		  Cval[Cnnz] = double_Complex_add (Cval[Cnnz], double_Complex_multiply (alpha, double_Complex_multiply (b_coeff, Aval[l])));
		  Cnnz++;
		}
	      else
		{
		  /* There's already a nonzero value in C at that
		   * location.  This means we need to add in the new
		   * nonzero value at that location, getting the
		   * location from the previous_index array. */
		  const int loc = previous_index[mm];
		  assert (loc != -1);
		  Cval[loc] = double_Complex_add (Cval[loc], double_Complex_multiply (alpha, double_Complex_multiply (b_coeff, Aval[l])));
		}
	    }
	}
    }

  /* Free up the temporary workspace */
  bebop_free (marker);
  bebop_free (previous_index);

  /* Return the matrix C through the output variables */
  *pCptr = Cptr;
  *pCind = Cind;
  *pCval = Cval;
  *pCnnz = Cnnz;

  /* No errors, so return zero */
  return 0;
}
#endif /* 0 */

#if 0
int
csr_matmatmult_pattern (int** pCptr, int** pCind, int* pCnnz,
			int* Bptr, int* Bind, 
			int* Aptr, int* Aind, 
			const int m, const int p, const int n)
{
  int i, j, k, l, mm;
  int* Cptr, *Cind;

  int* marker;
  int Cnnz = 0;
  
  return -1; /* BROKEN!!! */

  marker = (int*) bebop_malloc (n * sizeof (int));
  Cptr = (int*) bebop_malloc ((m+1) * sizeof (int));

  for (i = 0; i < n; i++)
    marker[i] = -1;
  
  Cptr[0] = Cnnz; /* = 0 */

  // for each row i of C
  for (i = 0; i < m; i++)
    {
      /* Go through the i-th row of B, using the index k.  This row
       * gives us the coefficients in the linear combination of the
       * rows of A.  Elements in the i-th row of B that are
       * (structurally) zero are skipped. */
      for (k = Bptr[i]; k < Bptr[i+1]; k++)
	{
	  /* Get the column index of the current coefficient
	   * (by which we will multiply the j-th row of A). */
	  j = Bind[k];
	  if (j < 0 || j >= p)
	    {
	      /* Column index out of range: clean up */
	      bebop_free (marker);
	      bebop_free (Cptr);
	      return -1; /* error flag */
	    }

	  /* Get the j-th row of A.  Any time we encounter a nonzero
	   * in this scan, there will be another nonzero in C. */
	  for (l = Aptr[j]; l < Aptr[j+1]; l++)
	    {
	      mm = Aind[l]; 
	      if (mm < 0 || mm >= n)
		{
		  /* Column index out of range: clean up */
		  bebop_free (marker);
		  bebop_free (Cptr);
		  return -2; /* error flag */
		}
	      /* Do this to avoid overcounting the number of nonzeros
	       * in C.  If we are currently computing row i of C, and
	       * we've already added C(i,mm) to C, then don't add it 
	       * again.  marker[mm] == i is true (on this i iteration)
	       * iff we have already added C(i,mm) to C. */
	      if (marker[mm] != i)
		{
		  marker[mm] = i; 
		  Cnnz++;
		}
	    }
	}
      Cptr[i+1] = Cnnz;
    }

  /* 
   * Now that we know how many nonzeros there will be in C, we can
   * allocate space for it. 
   */
  Cind = bebop_malloc (Cnnz * sizeof (int));

  for (i = 0; i < n; i++)
    marker[i] = -1;

  Cnnz = 0;
  for (i = 0; i < m; i++)
    {
      for (k = Bptr[i]; k < Bptr[i+1]; k++)
	{
	  j = Bind[k];
	  for (l = Aptr[j]; l < Aptr[j+1]; l++)
	    {
	      mm = Aind[l];
	      if (marker[mm] != i)
		{
		  /* Remember that we've added a nonzero to C at (i,
		   * mm).  If while computing this i-th row of C, we
		   * encounter this column index mm again, then we
		   * don't need to increment Cnnz or store to Cind,
		   * because there is a structural nonzero already
		   * there. */
		  marker[mm] = i;
		  Cind[Cnnz] = mm;
		  Cnnz++;
		}
	      /* No need for the else: this is a pattern matrix, so 
	       * we don't have to add to the nonzero value. */
	    }
	}
    }

  /* Free up the temporary workspace */
  bebop_free (marker);

  /* Return the matrix C through the output variables */
  *pCptr = Cptr;
  *pCind = Cind;
  *pCnnz = Cnnz;

  /* No errors, so return zero */
  return 0;
}
#endif /* 0 */
