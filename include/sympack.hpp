/*
	 Copyright (c) 2016 The Regents of the University of California,
	 through Lawrence Berkeley National Laboratory.  

   Author: Mathias Jacquelin
	 
   This file is part of symPACK. All rights reserved.

	 Redistribution and use in source and binary forms, with or without
	 modification, are permitted provided that the following conditions are met:

	 (1) Redistributions of source code must retain the above copyright notice, this
	 list of conditions and the following disclaimer.
	 (2) Redistributions in binary form must reproduce the above copyright notice,
	 this list of conditions and the following disclaimer in the documentation
	 and/or other materials provided with the distribution.
	 (3) Neither the name of the University of California, Lawrence Berkeley
	 National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
	 be used to endorse or promote products derived from this software without
	 specific prior written permission.

	 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
	 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
	 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
	 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
	 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

	 You are under no obligation whatsoever to provide any bug fixes, patches, or
	 upgrades to the features, functionality or performance of the source code
	 ("Enhancements") to anyone; however, if you choose to make your Enhancements
	 available either publicly, or directly to Lawrence Berkeley National
	 Laboratory, without imposing a separate written license agreement for such
	 Enhancements, then you hereby grant the following license: a non-exclusive,
	 royalty-free perpetual license to install, use, modify, prepare derivative
	 works, incorporate into other computer software, distribute, and sublicense
	 such enhancements or derivative works thereof, in binary and source code form.
*/
#ifndef _SYMPACK_HEADER_
#define _SYMPACK_HEADER_

#ifdef __cplusplus
#include  "sympack/Environment.hpp"
#include  "sympack/utility.hpp"
#include  "sympack/DistSparseMatrix.hpp"
#include  "sympack/symPACKMatrix.hpp"
#include  "sympack/Mapping.hpp"
#include  "sympack/CommTypes.hpp"
#include  "sympack/Ordering.hpp"
#include  "sympack/LogFile.hpp"
#endif

#include <mpi.h>


//namespace symPACK{
//
//extern template class symPACKMatrix<float>;
//extern template class symPACKMatrix<double>;
//extern template class symPACKMatrix<std::complex<float> >;
//extern template class symPACKMatrix<std::complex<double> >;
//
//extern template class supernodalTaskGraph<FBTask>;
//extern template class supernodalTaskGraph<CompTask>;
//}






#ifdef __cplusplus
extern "C" {
#endif


  int symPACK_Init(int *argc=NULL, char ***argv=NULL);
  int symPACK_Finalize();


  //returns an integer corresponding to a symPACK handle

  int symPACK_C_InitInstanceFloat(MPI_Comm ccomm);
  int symPACK_C_InitInstanceDouble(MPI_Comm ccomm);
  int symPACK_C_InitInstanceComplex(MPI_Comm ccomm);
  int symPACK_C_InitInstanceDoubleComplex(MPI_Comm ccomm);
  
  int symPACK_InitInstanceFloat(MPI_Fint * Fcomm);
  int symPACK_InitInstanceDouble(MPI_Fint * Fcomm);
  int symPACK_InitInstanceComplex(MPI_Fint * Fcomm);
  int symPACK_InitInstanceDoubleComplex(MPI_Fint * Fcomm);
  
  void symPACK_SymbolicFactorize(int * sp_handle, int * n, int * colptr , int * rowind);
  
  void symPACK_DistributeFloat(int * sp_handle, float * nzvals);
  void symPACK_DistributeDouble(int * sp_handle, double * nzvals);
  void symPACK_DistributeComplex(int * sp_handle, float * nzvals);
  void symPACK_DistributeDoubleComplex(int * sp_handle, double * nzvals);
  
  void symPACK_NumericalFactorize(int * sp_handle);
  
  void symPACK_NumericalSolveFloat(int * sp_handle, int * nrhs, float * rhs);
  void symPACK_NumericalSolveDouble(int * sp_handle, int * nrhs, double * rhs);
  void symPACK_NumericalSolveComplex(int * sp_handle, int * nrhs, float * rhs);
  void symPACK_NumericalSolveDoubleComplex(int * sp_handle, int * nrhs, double * rhs);
  
  void symPACK_FinalizeInstance(int * sp_handle);




#ifdef __cplusplus
}



#endif

#endif // _SYMPACK_HEADER_
