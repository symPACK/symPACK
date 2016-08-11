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
class Ordering{
  public:
    CommEnvironment * CommEnv_;
    SYMPACK::vector<Int> perm;
    SYMPACK::vector<Int> invp;

    Ordering():pStructure(NULL){};
    
    void SetCommEnvironment(CommEnvironment * CommEnv);

    SparseMatrixStructure * pStructure;
    void SetStructure(SparseMatrixStructure & aGlobal);
    void AMD();
    void MMD();

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
