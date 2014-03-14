/// @file etree.cpp
/// @brief Implementation of the elimination-tree related algorithms.
/// @author Mathias Jacquelin
/// @date 2012-08-31
#include "ETree.hpp"
#include "utility.hpp"

namespace LIBCHOLESKY{

  void DisjointSet::Initialize(Int n){
    pp_.Resize(n);
    SetValue(pp_,I_ZERO);
    root_.Resize(n);
  }

  void DisjointSet::Finalize(){
    pp_.Resize(0);
  }

  Int DisjointSet::makeSet(Int i){
    pp_(i-1)=i;
    return i;
  }

  Int DisjointSet::link(Int s, Int t){
    pp_(s-1)=t;

    return t;
  }

  Int DisjointSet::find(Int i){
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


  void DisjointSet::Union(Int s, Int t){
    Int tSet= find(t);
    Int sSet= find(s);
    sSet = link(sSet, tSet );
    root_(sSet-1) = t;
  }



  ETree::ETree(){

  }

  ETree::ETree(SparseMatrixStructure & aGlobal){
    ConstructETree(aGlobal);
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

      //logfileptr->OFS()<<"postNumber: "<<postNumber_<<std::endl;
      //logfileptr->OFS()<<"invPostNumber: "<<invPostNumber_<<std::endl;

      isPostOrdered_ = true;
    }

  }




  void ETree::ConstructETree(SparseMatrixStructure & aGlobal){
    TIMER_START(ConstructETree);

logfileptr->OFS()<<"ALIVE"<<std::endl;
    //Expand to symmetric storage
    aGlobal.ExpandSymmetric();

logfileptr->OFS()<<"ALIVE"<<std::endl;
    n_ = aGlobal.size;

logfileptr->OFS()<<"ALIVE"<<std::endl;

    parent_.Resize(n_);
    SetValue(parent_,I_ZERO );

    DisjointSet sets;
    sets.Initialize(n_);

    Int cset,rset,rroot,row;
    for (Int col = 1; col <= n_; col++) {
      parent_(col-1)=col; //1 based indexes
      cset = sets.makeSet (col);
      sets.Root(cset-1) = col;
      parent_(col-1) = n_; 

#ifdef _DEBUG_
      logfileptr->OFS()<<"Examining col "<<col<<std::endl;
#endif
      for (Int p = aGlobal.expColptr(col-1); p < aGlobal.expColptr(col); p++) {
        row = aGlobal.expRowind(p-1);

#ifdef _DEBUG_
        logfileptr->OFS()<<"Row = "<<row<<" vs col = "<<col<<std::endl;
#endif


        if (row >= col) continue;

        rset = sets.find(row);
        rroot = sets.Root(rset-1);
#ifdef _DEBUG_
        logfileptr->OFS()<<"Row "<<row<<" is in set "<<rset<<" represented by "<<rroot<<std::endl;
#endif

        if (rroot != col) {
          parent_(rroot-1) = col;
          cset = sets.link(cset, rset);
          sets.Root(cset-1) = col;
#ifdef _DEBUG_
          logfileptr->OFS()<<"Parent of "<<rroot<<" is "<<col<<" which now represents set"<<cset<<std::endl;
#endif
        }
      }

    }


    parent_(n_-1) = 0;

    TIMER_STOP(ConstructETree);
  }


}
