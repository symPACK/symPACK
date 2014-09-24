/// @file ETree.hpp
/// @brief Various elimination tree related subroutines.
/// @author Mathias Jacquelin
/// @date 2013-08-31
#ifndef _ETREE_HPP_ 
#define _ETREE_HPP_

#include  <stdlib.h>

#include "ngchol/Environment.hpp"
#include "ngchol/NumVec.hpp"
#include "ngchol/SparseMatrixStructure.hpp"
#include "ngchol/Ordering.hpp"

namespace LIBCHOLESKY{
  class DisjointSet{
    protected:
      NumVec<Int> pp_;
      NumVec<Int> root_;
    public:
      void Initialize(Int n);
      void Finalize();
      inline Int makeSet(Int i);
      inline Int link(Int s, Int t);
      inline void Union(Int s, Int t, Int root = -1);
      inline Int find(Int i);
      inline Int & Root(Int i){return root_(i);};
  };

class ETree{

  friend std::ostream& operator<<( std::ostream& os, const ETree& tree);
  template <class F> friend class DistSparseMatrix;

protected:

  void BTreeToPO(IntNumVec & fson, IntNumVec & brother, IntNumVec & perm);

public:
  ETree();
  ETree(SparseMatrixStructure & aGlobal, Ordering & aOrder);

  void ConstructETree(SparseMatrixStructure & aGlobal, Ordering & aOrder);
//  void ConstructETree(SparseMatrixStructure & aGlobal,IntNumVec & perm, IntNumVec & invp);
//  void ConstructETree2(SparseMatrixStructure & aGlobal);
//  void ConstructETree(int n, int * xadj, int * adj);

  void PostOrderTree(Ordering & aOrder);
  void DeepestFirst(Ordering & aOrder);

  ETree ToSupernodalETree(IntNumVec & aXsuper,IntNumVec & aSupMembership,Ordering & aOrder) const;

  inline bool IsPostOrdered() const { return bIsPostOrdered_;};
  inline Int n() const { return n_; };
//  inline Int ToPostOrder(Int i) const { if(!bIsPostOrdered_){ throw std::logic_error("Tree must be postordered to use this function."); }  return postNumber_(i-1);};
//  inline Int FromPostOrder(Int i) const  { if(!bIsPostOrdered_){ throw std::logic_error("Tree must be postordered to use this function."); }  return invPostNumber_(i-1);};
//  inline IntNumVec ToPostOrder(IntNumVec & vec) const { if(!bIsPostOrdered_){ throw std::logic_error("Tree must be postordered to use this function."); } IntNumVec povec = vec; for(Int i=0;i<povec.m();i++){ if(vec[i]!=0){ povec[i]=postNumber_(vec[i]-1);} }   return povec;};
  inline Int PostParent(Int i) const {
      return poparent_[i]; 
//      if(!bIsPostOrdered_){ throw std::logic_error("Tree must be postordered to use this function."); } 
//
//      Int node = invPostNumber_(i);
//      Int parent = parent_(node-1);
//      parent = parent==0?0:postNumber_(parent-1);
//      assert(parent == poparent_[i]);
//      return parent; 
  }
  inline Int Parent(Int i) const { return parent_(i); };
  inline Int Size() const { return parent_.m(); };

//  inline IntNumVec & Perm(){ return postNumber_; }
//  inline IntNumVec & Invp(){ return invPostNumber_; }

  void SortChildren(IntNumVec & cc, Ordering & aOrder);

//  void PermuteTree(IntNumVec & perm);

protected:
  Int n_;
  bool bIsPostOrdered_;

  NumVec<Int> parent_;
  NumVec<Int> poparent_;
//  NumVec<Int> postNumber_;
//  NumVec<Int> invPostNumber_;

};






  inline Int DisjointSet::makeSet(Int i){
    pp_(i-1)=i;
    return i;
  }

  inline Int DisjointSet::link(Int s, Int t){
    pp_(s-1)=t;

    return t;
  }

  inline Int DisjointSet::find(Int i){
    Int p, gp;

    p = pp_(i-1);
    gp=pp_(p-1);

    while(gp!=p){
      i = makeSet(gp);
      p = pp_(i-1);
      gp = pp_(p-1);
    }

    return p;
  }


  inline void DisjointSet::Union(Int s, Int t, Int root){
    Int tSet= find(t);
    Int sSet= find(s);
    sSet = link(sSet, tSet );
    root_(sSet-1) = root==-1?t:root;
  }














}
#endif
