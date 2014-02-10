/// @file etree.cpp
/// @brief Implementation of the elimination-tree related algorithms.
/// @author Mathias Jacquelin
/// @date 2012-08-31
#include "ETree.hpp"
#include "utility.hpp"

namespace LIBCHOLESKY{

  void ETree::DisjointSet::Initialize(Int n){
    pp_.Resize(n);
    SetValue(pp_,I_ZERO);
    root_.Resize(n);
  }

  void ETree::DisjointSet::Finalize(){
    pp_.Resize(0);
  }

  Int ETree::DisjointSet::makeSet(Int i){
    pp_(i-1)=i;
    return i;
  }

  Int ETree::DisjointSet::link(Int s, Int t){
    pp_(s-1)=t;

    return t;
  }

  Int ETree::DisjointSet::find(Int i){
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


  void ETree::DisjointSet::Union(Int s, Int t){
      Int tSet= find(t);
      Int sSet= find(s);
      sSet = link(sSet, tSet );
      root_(sSet-1) = t;
  }



ETree::ETree(){

}


  void ETree::PostOrderTree(){
    if(n_>0 && !isPostOrdered_){


      IntNumVec fson(n_);
      SetValue(fson, I_ZERO);
      IntNumVec brother(n_);
      SetValue(brother, I_ZERO);


      //Get Binary tree representation
      for(Int vertex=n_-1; vertex>0; vertex--){
        Int curParent = parent_(vertex-1);
        brother(vertex-1) = fson(curParent-1);
        fson(curParent-1) = vertex;
      }

//      statusOFS<<"parent "<<parent_<<std::endl;
//      statusOFS<<"fson "<<fson<<std::endl;
//      statusOFS<<"brother "<<brother<<std::endl;

    
      //Do a depth first search to construct the postordered tree
      IntNumVec stack(n_);
      postNumber_.Resize(n_);
      invPostNumber_.Resize(n_);

      Int stacktop=0, vertex=n_,m=0;
      while( m<n_){
        while(vertex>0){
          stacktop++;
          stack(stacktop-1) = vertex;
          vertex = fson(vertex-1);
        }
       
        while(vertex==0 && stacktop>0){
          vertex = stack(stacktop-1);
          stacktop--;
          m++;

          postNumber_(vertex-1) = m;
          invPostNumber_(m-1) = vertex;

          vertex = brother(vertex-1);
        }

      }

//      postParent_.Resize(n_);
      //modify the parent list ?
      // node i is now node postNumber(i-1)
//      for(Int i=1; i<=n_;i++){
//        postParent_(postNumber_(i-1)-1)=postNumber_(parent_(i-1)-1);
//      }

      isPostOrdered_ = true;
    }

  }




}
