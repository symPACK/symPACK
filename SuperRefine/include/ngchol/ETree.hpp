/// @file ETree.hpp
/// @brief Various elimination tree related subroutines.
/// @author Mathias Jacquelin
/// @date 2013-08-31
#ifndef _ETREE_HPP_ 
#define _ETREE_HPP_

#include  <stdlib.h>

#include <vector>
//#include "ngchol/Environment.hpp"
//#include "ngchol/NumVec.hpp"
#include "ngchol/SparseMatrixStructure.hpp"
#include "ngchol/Ordering.hpp"


using namespace std;
namespace LIBCHOLESKY{
  class DisjointSet{
    protected:
      vector<int> pp_;
      vector<int> root_;
    public:
      void Initialize(int n);
      void Finalize();
      inline int makeSet(int i);
      inline int link(int s, int t);
      inline void Union(int s, int t, int root = -1);
      inline int find(int i);
      inline int & Root(int i){return root_[i];};
  };

class ETree{

  template <class F> friend class DistSparseMatrix;

protected:

  void BTreeToPO(vector<int> & fson, vector<int> & brother, vector<int> & perm);

public:
  ETree();
  ETree(SparseMatrixStructure & aGlobal, Ordering & aOrder);

  void ConstructETree(SparseMatrixStructure & aGlobal, Ordering & aOrder);

  void PostOrderTree(Ordering & aOrder);
  void DeepestFirst(Ordering & aOrder);

  ETree ToSupernodalETree(vector<int> & aXsuper,vector<int> & aSupMembership,Ordering & aOrder) const;

  inline bool IsPostOrdered() const { return bIsPostOrdered_;};
  inline int n() const { return n_; };
  inline int PostParent(int i) const {
      return poparent_[i]; 
  }
  inline int Parent(int i) const { return parent_[i]; };
  inline int Size() const { return parent_.size(); };

  void SortChildren(vector<int> & cc, Ordering & aOrder);


protected:
  int n_;
  bool bIsPostOrdered_;

  vector<int> parent_;
  vector<int> poparent_;

};






  inline int DisjointSet::makeSet(int i){
    pp_[i-1]=i;
    return i;
  }

  inline int DisjointSet::link(int s, int t){
    pp_[s-1]=t;

    return t;
  }

  inline int DisjointSet::find(int i){
    int p, gp;

    p = pp_[i-1];
    gp=pp_[p-1];

    while(gp!=p){
      i = makeSet(gp);
      p = pp_[i-1];
      gp = pp_[p-1];
    }

    return p;
  }


  inline void DisjointSet::Union(int s, int t, int root){
    int tSet= find(t);
    int sSet= find(s);
    sSet = link(sSet, tSet );
    root_[sSet-1] = root==-1?t:root;
  }














}
#endif
