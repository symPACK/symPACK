#ifndef _ETREE_HPP_ 
#define _ETREE_HPP_

#include  <stdlib.h>

#include "sympack/Environment.hpp"
#include "sympack/DistSparseMatrixGraph.hpp"
#include "sympack/Ordering.hpp"


namespace symPACK{
  class DisjointSet{
    protected:
      std::vector<Int> pp_;
      std::vector<Int> root_;
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

  void BTreeToPO(std::vector<Int> & fson, std::vector<Int> & brother, std::vector<Int> & perm);

public:
  ETree();

  //not working
  void ConstructETree(DistSparseMatrixGraph & aDistExp, Ordering & aOrder);
  
  void ConstructETree(SparseMatrixGraph & sgraph, Ordering & aOrder, MPI_Comm & aComm);
  void PostOrderTree(Ordering & aOrder,Int * relinvp = nullptr);
  void DeepestFirst(Ordering & aOrder);

  ETree ToSupernodalETree(std::vector<Int> & aXsuper,std::vector<Int> & aSupMembership,Ordering & aOrder) const;

  inline bool IsPostOrdered() const { return bIsPostOrdered_;};
  inline Int n() const { return n_; };
  inline Int PostParent(Int i) const { return poparent_[i]; }
  inline Int Parent(Int i) const { return parent_[i]; };
  inline Int Size() const { return parent_.size(); };
  void SortChildren(std::vector<Int> & cc, Ordering & aOrder);

protected:
  Int n_;
  bool bIsPostOrdered_;

  std::vector<Int> parent_;
  std::vector<Int> poparent_;
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
#endif //_ETREE_HPP_
