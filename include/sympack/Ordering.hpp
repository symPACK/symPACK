/******************************************************************************
    AMD, Copyright (c), 1996-2015, Timothy A. Davis,
    Patrick R. Amestoy, and Iain S. Duff.  All Rights Reserved.
    Used in symPACK under the BSD 3-clause license.
******************************************************************************/

#ifndef _ORDERING_HPP_ 
#define _ORDERING_HPP_

#include  <stdlib.h>

#include "sympack/Environment.hpp"
#include "sympack/SparseMatrixStructure.hpp"
#include "sympack/DistSparseMatrixGraph.hpp"

namespace SYMPACK{
class SparseMatrixGraph;
class DistSparseMatrixGraph;
class SparseMatrixStructure;
class Ordering;




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
    CommEnvironment * CommEnv_;
    SYMPACK::vector<Int> perm;
    SYMPACK::vector<Int> invp;

    Ordering():pStructure(NULL){};
    
    void SetCommEnvironment(CommEnvironment * CommEnv);

    SparseMatrixStructure * pStructure;
//    void SetStructure(SparseMatrixStructure & aGlobal);

    void MMD(const SparseMatrixGraph & g);
    void AMD(const SparseMatrixGraph & g);
    void NDBOX(Int size);
    void NDGRID(Int size);

#ifdef USE_PARMETIS
  void PARMETIS(const DistSparseMatrixGraph & g);
#endif

#ifdef USE_PTSCOTCH
  void PTSCOTCH(const DistSparseMatrixGraph & g);
#endif

#ifdef USE_METIS
    void METIS(const SparseMatrixGraph & g);
#endif

#ifdef USE_SCOTCH
    void SCOTCH(const SparseMatrixGraph & g);
#endif
    void Compose(SYMPACK::vector<Int> & invp2);
};

}
#endif // _ORDERING_HPP_
