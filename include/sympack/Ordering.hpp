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
/******************************************************************************
    AMD, Copyright (c), 1996-2015, Timothy A. Davis,
    Patrick R. Amestoy, and Iain S. Duff.  All Rights Reserved.
    Used in symPACK under the BSD 3-clause license.
******************************************************************************/

#ifndef _ORDERING_HPP_ 
#define _ORDERING_HPP_

#include  <stdlib.h>

#include "sympack/Environment.hpp"
#include "sympack/DistSparseMatrixGraph.hpp"

namespace symPACK{
class SparseMatrixGraph;
class DistSparseMatrixGraph;
class Ordering;


#if RCMIDXSIZE==64
  typedef    int64_t   RCMInt;
#else
  typedef    int32_t   RCMInt;
#endif


#if MMDIDXSIZE==64
  typedef    int64_t   MMDInt;
#else
  typedef    int32_t   MMDInt;
#endif

#if AMDIDXSIZE==64
  typedef    int64_t   AMDInt;
#else
  typedef    int32_t   AMDInt;
#endif






class Ordering{
  public:
    std::vector<Int> perm;
    std::vector<Int> invp;
    int NpOrdering;

    Ordering():NpOrdering(0){};
    


    void RCM(const SparseMatrixGraph & g, MPI_Comm comm);
    void MMD(const SparseMatrixGraph & g, MPI_Comm comm);
    void AMD(const SparseMatrixGraph & g, MPI_Comm comm);
    void NDBOX(Int size, MPI_Comm comm);
    void NDGRID(Int size, MPI_Comm comm);

#ifdef USE_PARMETIS
  void PARMETIS(const DistSparseMatrixGraph & g);
#endif

#ifdef USE_PTSCOTCH
  void PTSCOTCH(const DistSparseMatrixGraph & g);
#endif

#ifdef USE_METIS
    void METIS(const SparseMatrixGraph & g, MPI_Comm comm);
#endif

#ifdef USE_SCOTCH
    void SCOTCH(const SparseMatrixGraph & g, MPI_Comm comm);
#endif
    void Compose(std::vector<Int> & invp2);
};

}
#endif // _ORDERING_HPP_
