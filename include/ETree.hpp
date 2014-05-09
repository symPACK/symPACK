/// @file ETree.hpp
/// @brief Various elimination tree related subroutines.
/// @author Mathias Jacquelin
/// @date 2013-08-31
#ifndef _ETREE_HPP_ 
#define _ETREE_HPP_

#include  <stdlib.h>

#include "Environment.hpp"
#include "NumVec.hpp"
#include "SparseMatrixStructure.hpp"

namespace LIBCHOLESKY{
  class DisjointSet{
    protected:
      NumVec<Int> pp_;
      NumVec<Int> root_;
    public:
      void Initialize(Int n);
      void Finalize();
      Int makeSet(Int i);
      Int link(Int s, Int t);
      void Union(Int s, Int t);
      Int find(Int i);
      inline Int & Root(Int i){return root_(i);};
  };

class ETree{

  friend std::ostream& operator<<( std::ostream& os, const ETree& tree);
  template <class F> friend class DistSparseMatrix;

protected:

  void BTreeToPO(IntNumVec & fson, IntNumVec & brother);

public:
  ETree();
  ETree(SparseMatrixStructure & aGlobal);

  void ConstructETree(SparseMatrixStructure & aGlobal);
  void ConstructETree2(SparseMatrixStructure & aGlobal);

  void PostOrderTree();

  ETree ToSupernodalETree(IntNumVec & aXsuper) const;

  inline bool IsPostOrdered() const { return bIsPostOrdered_;};
  inline Int n() const { return n_; };
  inline Int ToPostOrder(Int i) const { if(!bIsPostOrdered_){ throw std::logic_error("Tree must be postordered to use this function."); }  return postNumber_(i-1);};
  inline Int FromPostOrder(Int i) const  { if(!bIsPostOrdered_){ throw std::logic_error("Tree must be postordered to use this function."); }  return invPostNumber_(i-1);};
  inline IntNumVec ToPostOrder(IntNumVec & vec) const { if(!bIsPostOrdered_){ throw std::logic_error("Tree must be postordered to use this function."); } IntNumVec povec = vec; for(Int i=0;i<povec.m();i++){ if(vec[i]!=0){ povec[i]=postNumber_(vec[i]-1);} }   return povec;};
  inline Int PostParent(Int i) const { 
      if(!bIsPostOrdered_){ throw std::logic_error("Tree must be postordered to use this function."); } 
      Int node = invPostNumber_(i);
      Int parent = parent_(node-1);
      return parent==0?0:postNumber_(parent-1); 
  }
  inline Int Parent(Int i) const { return parent_(i); };
  inline Int Size() const { return parent_.m(); };


  IntNumVec SortChildren(IntNumVec & cc);


protected:
  Int n_;
  bool bIsPostOrdered_ = false;

  NumVec<Int> parent_;
  NumVec<Int> postNumber_;
  NumVec<Int> invPostNumber_;
//  NumVec<Int> postParent_;

};


}
#endif
