/// @file ETree.hpp
/// @brief Various elimination tree related subroutines.
/// @author Mathias Jacquelin
/// @date 2013-08-31
#ifndef _ETREE_HPP_ 
#define _ETREE_HPP_

#include  <stdlib.h>

#include "Environment.hpp"
#include "NumVec.hpp"

namespace LIBCHOLESKY{

class ETree{

  friend std::ostream& operator<<( std::ostream& os, const ETree& tree);
  template <class F> friend class DistSparseMatrix;

protected:
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

public:
  ETree();
  void PostOrderTree();
  Int lca(Int u, Int v);

  inline Int n() const { return n_; };
  inline Int ToPostOrder(Int i) const { if(!isPostOrdered_){ throw std::logic_error("Tree must be postordered to use this function."); }  return postNumber_(i-1);};
  inline Int FromPostOrder(Int i) const  { if(!isPostOrdered_){ throw std::logic_error("Tree must be postordered to use this function."); }  return invPostNumber_(i-1);};
  inline Int PostParent(Int i) const { if(!isPostOrdered_){ throw std::logic_error("Tree must be postordered to use this function."); } return postNumber_(parent_(invPostNumber_(i-1)-1)-1); }
protected:
  Int n_;
  bool isPostOrdered_ = false;

  NumVec<Int> parent_;
  NumVec<Int> postNumber_;
  NumVec<Int> invPostNumber_;
//  NumVec<Int> postParent_;

};


}
#endif
