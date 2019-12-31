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

    void GetRelativeInvp(const std::vector<Int> & frominvp, std::vector<Int> & relinvp);
    void Compose(std::vector<Int> & invp2);
};

}
#endif // _ORDERING_HPP_
