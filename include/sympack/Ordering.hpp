#ifndef _ORDERING_HPP_ 
#define _ORDERING_HPP_

#include  <stdlib.h>

#include "sympack/Environment.hpp"
#include "sympack/SparseMatrixStructure.hpp"
#include "sympack/DistSparseMatrixGraph.hpp"

namespace SYMPACK{
class SparseMatrixStructure;
class Ordering{
  public:
    CommEnvironment * CommEnv_;
    vector<Int> perm;
    vector<Int> invp;
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
#endif

#ifdef USE_PTSCOTCH
  //FIXME: currently, colptr and rowind are not distributed
  void PTSCOTCH();
  void PTSCOTCH(const DistSparseMatrixGraph & g);
#endif

#ifdef USE_METIS
    void METIS();
#endif

#ifdef USE_SCOTCH
    void SCOTCH();
#endif
    void Compose(vector<Int> & invp2);
};

}
#endif // _ORDERING_HPP_
