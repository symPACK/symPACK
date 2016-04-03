/// @file ETree.hpp
/// @brief Various elimination tree related subroutines.
/// @author Mathias Jacquelin
/// @date 2013-08-31
#ifndef _ETREE_HPP_ 
#define _ETREE_HPP_

#include  <stdlib.h>

#include "sympack/Environment.hpp"
#include "sympack/DistSparseMatrixGraph.hpp"
#include "sympack/SparseMatrixStructure.hpp"
#include "sympack/Ordering.hpp"


namespace SYMPACK{
  class DisjointSet{
    protected:
      SYMPACK::vector<Int> pp_;
      SYMPACK::vector<Int> root_;
    public:
      void Initialize(Int n);
      void Finalize();
      inline Int makeSet(Int i);
      inline Int link(Int s, Int t);
      inline void Union(Int s, Int t, Int root = -1);
      inline Int find(Int i);
      inline Int & Root(Int i){return root_[i];};
  };

class ETree{

  friend std::ostream& operator<<( std::ostream& os, const ETree& tree);
  template <class F> friend class DistSparseMatrix;

protected:

  void BTreeToPO(SYMPACK::vector<Int> & fson, SYMPACK::vector<Int> & brother, SYMPACK::vector<Int> & perm);

public:
  ETree();
  ETree(SparseMatrixStructure & aGlobal, Ordering & aOrder);

  void ConstructETree(DistSparseMatrixGraph & aDistExp, Ordering & aOrder);
  void ConstructETree(SparseMatrixStructure & aGlobal, Ordering & aOrder);
  void PostOrderTree(Ordering & aOrder);
  void DeepestFirst(Ordering & aOrder);

  ETree ToSupernodalETree(SYMPACK::vector<Int> & aXsuper,SYMPACK::vector<Int> & aSupMembership,Ordering & aOrder) const;

  inline bool IsPostOrdered() const { return bIsPostOrdered_;};
  inline Int n() const { return n_; };
  inline Int PostParent(Int i) const { return poparent_[i]; }
  inline Int Parent(Int i) const { return parent_[i]; };
  inline Int Size() const { return parent_.size(); };
  void SortChildren(SYMPACK::vector<Int> & cc, Ordering & aOrder);

protected:
  Int n_;
  bool bIsPostOrdered_;

  SYMPACK::vector<Int> parent_;
  SYMPACK::vector<Int> poparent_;
};






  inline Int DisjointSet::makeSet(Int i){
    pp_[i-1]=i;
    return i;
  }

  inline Int DisjointSet::link(Int s, Int t){
    pp_[s-1]=t;

    return t;
  }

  inline Int DisjointSet::find(Int i){
    Int p, gp;

    p = pp_[i-1];
    gp=pp_[p-1];

    while(gp!=p){
      i = makeSet(gp);
      p = pp_[i-1];
      gp = pp_[p-1];
    }

    return p;
  }


  inline void DisjointSet::Union(Int s, Int t, Int root){
    Int tSet= find(t);
    Int sSet= find(s);
    sSet = link(sSet, tSet );
    root_[sSet-1] = root==-1?t:root;
  }














}
#endif
