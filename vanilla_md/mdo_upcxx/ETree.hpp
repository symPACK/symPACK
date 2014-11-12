/// @file ETree.hpp
/// @brief Various elimination tree related subroutines.
/// @author Mathias Jacquelin
/// @date 2013-08-31
#ifndef _ETREE_HPP_ 
#define _ETREE_HPP_

#include  <stdlib.h>
#include <vector>
#include <stdexcept>

#define I_ZERO 0

using namespace std;

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

protected:

  void BTreeToPO(vector<int> & fson, vector<int> & brother);

public:
  ETree();

  void ConstructETree(int n, int * xadj, int * adj);
  void PostOrderTree();
  void Dump();

  inline bool IsPostOrdered() const { return bIsPostOrdered_;};
  inline int n() const { return n_; };
  inline int ToPostOrder(int i) const { if(!bIsPostOrdered_){ throw std::logic_error("Tree must be postordered to use this function."); }  return postNumber_[i-1];};
  inline int FromPostOrder(int i) const  { if(!bIsPostOrdered_){ throw std::logic_error("Tree must be postordered to use this function."); }  return invPostNumber_[i-1];};
  inline vector<int> ToPostOrder(vector<int> & vec) const { if(!bIsPostOrdered_){ throw std::logic_error("Tree must be postordered to use this function."); } vector<int> povec = vec; for(int i=0;i<povec.size();i++){ if(vec[i]!=0){ povec[i]=postNumber_[vec[i]-1];} }   return povec;};
  inline int PostParent(int i) const { 
      if(!bIsPostOrdered_){ throw std::logic_error("Tree must be postordered to use this function."); } 
      int node = invPostNumber_[i];
      int parent = parent_[node-1];
      return parent==0?0:postNumber_[parent-1]; 
  }
  inline int Parent(int i) const { return parent_[i]; };
  inline int Size() const { return parent_.size(); };


  vector<int> SortChildren(vector<int> & cc);
  void PermuteTree(vector<int> & perm);

protected:
  int n_;
  bool bIsPostOrdered_;

  vector<int> parent_;
  vector<int> postNumber_;
  vector<int> invPostNumber_;
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














#endif
