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
    SparseMatrixStructure * pStructure;

    Ordering():pStructure(NULL){};
    
    void SetCommEnvironment(CommEnvironment * CommEnv);
    void SetStructure(SparseMatrixStructure & aGlobal);
    void MMD();
    void AMD();
    void NDBOX();
    void NDGRID();

#ifdef USE_PARMETIS
  //FIXME: currently, colptr and rowind are not distributed
  void PARMETIS();
  void PARMETIS(const DistSparseMatrixGraph & g);
#endif

#ifdef USE_PTSCOTCH
  //FIXME: currently, colptr and rowind are not distributed
  void PTSCOTCH();
  void PTSCOTCH(const DistSparseMatrixGraph & g);
#endif

#ifdef USE_METIS
    void METIS();
    void METIS(const SparseMatrixGraph & g);
#endif

#ifdef USE_SCOTCH
    void SCOTCH();
    void SCOTCH(const SparseMatrixGraph & g);
#endif
    void Compose(SYMPACK::vector<Int> & invp2);
};

}
#endif // _ORDERING_HPP_
