/// @file etree.cpp
/// @brief Implementation of the elimination-tree related algorithms.
/// @author Mathias Jacquelin
/// @date 2012-08-31
#include "ngchol/ETree.hpp"
#include "ngchol/utility.hpp"

namespace LIBCHOLESKY{


  void DisjointSet::Initialize(Int n){
    pp_.Resize(n);
    SetValue(pp_,I_ZERO);
    root_.Resize(n);
  }

  void DisjointSet::Finalize(){
    pp_.Resize(0);
  }




  ETree::ETree(){

    bIsPostOrdered_=false;
  }

  ETree::ETree(SparseMatrixStructure & aGlobal, Ordering & aOrder){
    bIsPostOrdered_ = false;
    ConstructETree(aGlobal,aOrder);
  }



  void ETree::BTreeToPO(IntNumVec & fson, IntNumVec & brother, IntNumVec & invpos){
      //Do a depth first search to construct the postordered tree
      IntNumVec stack(n_);
      //postNumber_.Resize(n_);
      invpos.Resize(n_);

      Int stacktop=0, vertex=n_,m=0;
      bool exit = false;
      while( m<n_){
        do{
          stacktop++;
          stack(stacktop-1) = vertex;
          vertex = fson(vertex-1);
        }while(vertex>0);

        while(vertex==0){

          if(stacktop<=0){
            exit = true;
            break;
          }
          vertex = stack(stacktop-1);
          stacktop--;
          m++;

          invpos[vertex-1] = m;
//          perm[m-1] = vertex;

          vertex = brother(vertex-1);
        }

        if(exit){
          break;
        }
      }
}


  void ETree::PostOrderTree(Ordering & aOrder){
    if(n_>0 && !bIsPostOrdered_){

     TIMER_START(PostOrder);

      IntNumVec fson(n_);
      SetValue(fson, I_ZERO);
      IntNumVec & brother = poparent_;
      brother.Resize(n_);
      SetValue(brother, I_ZERO);

      Int lroot = n_;
      for(Int vertex=n_-1; vertex>0; vertex--){
        Int curParent = parent_(vertex-1);
        if(curParent==0 || curParent == vertex){
          brother(lroot-1) = vertex;
          lroot = vertex;
        }
        else{
          brother(vertex-1) = fson(curParent-1);
          fson(curParent-1) = vertex;
        }
      }


#ifdef _DEBUG_
      logfileptr->OFS()<<"parent "<<parent_<<std::endl;
      logfileptr->OFS()<<"fson "<<fson<<std::endl;
      logfileptr->OFS()<<"brother "<<brother<<std::endl;
#endif

      IntNumVec invpos;
      BTreeToPO(fson,brother,invpos);

      //modify the parent list ?
      // node i is now node invpos(i-1)
      for(Int i=1; i<=n_;i++){
        Int nunode = invpos[i-1];
        Int ndpar = parent_[i-1];
        if(ndpar>0){
          ndpar = invpos[ndpar-1];
        }
        poparent_[nunode-1] = ndpar;
      }

      //we need to compose aOrder.invp and invpos
      aOrder.Compose(invpos);
      


      bIsPostOrdered_ = true;


#ifdef _DEBUG_
      logfileptr->OFS()<<"new parent: "<<brother<<std::endl;
#endif

     TIMER_STOP(PostOrder);
    }

  }


  void ETree::SortChildren(IntNumVec & cc, Ordering & aOrder){
    if(!bIsPostOrdered_){
      this->PostOrderTree(aOrder);
    }

      IntNumVec fson(n_);
      SetValue(fson, I_ZERO);
      IntNumVec brother(n_);
      SetValue(brother, I_ZERO);
      IntNumVec lson(n_);
      SetValue(lson, I_ZERO);

      //Get Binary tree representation
      Int lroot = n_;
      for(Int vertex=n_-1; vertex>0; vertex--){
        Int curParent = PostParent(vertex-1);
        //Int curParent = parent_(vertex-1);
        if(curParent==0 || curParent == vertex){
          brother(lroot-1) = vertex;
          lroot = vertex;
        }
        else{
          Int ndlson = lson(curParent-1);
          if(ndlson > 0){
             if  ( cc(vertex-1) >= cc(ndlson-1) ) {
             //if  ( cc(ToPostOrder(vertex)-1) >= cc(ToPostOrder(ndlson)-1) ) {
                brother(vertex-1) = fson(curParent-1);
                fson(curParent-1) = vertex;
             }
             else{                                                                                                                                            
                brother(ndlson-1) = vertex;
                lson(curParent-1) = vertex;
             }                                                                                                                                                               
          }
          else{
             fson(curParent-1) = vertex;
             lson(curParent-1) = vertex;
          }
        }
      }
      brother(lroot-1)=0;


      IntNumVec invpos;
//      IntNumVec invperm;

      //Compute the parent permutation and update postNumber_
      //Do a depth first search to construct the postordered tree
      IntNumVec stack(n_);
      invpos.Resize(n_);
//      invperm.Resize(n_);

      Int stacktop=0, vertex=n_,m=0;
      bool exit = false;
      while( m<n_){
        do{
          stacktop++;
          stack(stacktop-1) = vertex;
          vertex = fson(vertex-1);
        }while(vertex>0);

        while(vertex==0){

          if(stacktop<=0){
            exit = true;
            break;
          }
          vertex = stack(stacktop-1);
          stacktop--;
          m++;

          invpos[vertex-1] = m;

          vertex = brother(vertex-1);
        }

        if(exit){
          break;
        }
      }


      //modify the parent list ?
      // node i is now node invpos(i-1)
      for(Int i=1; i<=n_;i++){
        Int nunode = invpos[i-1];
        Int ndpar = poparent_[i-1];
        if(ndpar>0){
          ndpar = invpos[ndpar-1];
        }
        brother[nunode-1] = ndpar;
      }
      poparent_ = brother;




      //Permute CC     
      for(Int node = 1; node <= n_; ++node){
        Int nunode = invpos[node-1];
        stack[nunode-1] = cc[node-1];
      }

      for(Int node = 1; node <= n_; ++node){
        cc[node-1] = stack[node-1];
      }

      //Compose the two permutations
      aOrder.Compose(invpos);

      //Sort the parent

////      IntNumVec poParent(n_+1);
////            for(Int i=1; i<=n_;i++){
////              Int nunode = perm_(i-1);
////              Int ndpar = parent_(i-1);
////              if(ndpar>0){
////                ndpar = postNumber_(ndpar-1);
////              }
////              poParent(nunode-1) = ndpar;
////            }
////#ifdef _DEBUG_
////      logfileptr->OFS()<<"ORDERED fson "<<fson<<std::endl;
////      logfileptr->OFS()<<"ORDERED brother "<<brother<<std::endl;
////
////      IntNumVec poParent(n_+1);
////            for(Int i=1; i<=n_;i++){
////              Int nunode = postNumber_(i-1);
////              Int ndpar = parent_(i-1);
////              if(ndpar>0){
////                ndpar = postNumber_(ndpar-1);
////              }
////              poParent(nunode-1) = ndpar;
////            }
////
////
////      logfileptr->OFS()<<"ORDERED new parent: "<<poParent<<std::endl;
////      logfileptr->OFS()<<"ORDERED postNumber: "<<postNumber_<<std::endl;
////      logfileptr->OFS()<<"ORDERED invPostNumber: "<<invPostNumber_<<std::endl;
////#endif
////
//////      logfileptr->OFS()<<"perm: "<<perm<<std::endl;
//////      logfileptr->OFS()<<"invperm: "<<invperm<<std::endl;



//    return perm;

  }


  void ETree::ConstructETree(SparseMatrixStructure & aGlobal, Ordering & aOrder){
    bIsPostOrdered_=false;
    n_ = aGlobal.size;

    //Expand to symmetric storage
    aGlobal.ExpandSymmetric();

TIMER_START(Construct_Etree_Classic);
    parent_.Resize(n_);
    SetValue(parent_,I_ZERO );


    IntNumVec ancstr(n_);



    for(Int i = 1; i<=n_; ++i){
            parent_(i-1) = 0;
            ancstr(i-1) = 0;
            Int node = aOrder.perm[i-1];

            Int jstrt = aGlobal.expColptr(node-1);
            Int jstop = aGlobal.expColptr(node) - 1;
            if  ( jstrt < jstop ){
              for(Int j = jstrt; j<=jstop; ++j){
                    Int nbr = aGlobal.expRowind(j-1);
                    nbr = aOrder.invp[nbr-1];
                    if  ( nbr < i ){
//                       -------------------------------------------
//                       for each nbr, find the root of its current
//                       elimination tree.  perform path compression
//                       as the subtree is traversed.
//                       -------------------------------------------
                      Int break_loop = 0;
                      if  ( ancstr(nbr-1) == i ){
                        break_loop = 1;
                      }
                      else{
                        while(ancstr(nbr-1) >0){
                          if  ( ancstr(nbr-1) == i ){
                            break_loop = 1;
                            break;
                          }
                          Int next = ancstr(nbr-1);
                          ancstr(nbr-1) = i;
                          nbr = next;
                        }
                        //                       --------------------------------------------
                        //                       now, nbr is the root of the subtree.  make i
                        //                       the parent node of this root.
                        //                       --------------------------------------------
                        if(!break_loop){
                          parent_(nbr-1) = i;
                          ancstr(nbr-1) = i;
                        }
                      }
                    }
              }
            }
  }

TIMER_STOP(Construct_Etree_Classic);

  }





  ETree ETree::ToSupernodalETree(IntNumVec & aXsuper,IntNumVec & aSupMembership,Ordering & aOrder) const{
    ETree newTree;
    newTree.n_ = aXsuper.m()-1;
    newTree.parent_.Resize(aXsuper.m()-1);
    

assert(bIsPostOrdered_);

    for(Int snode=1; snode<=newTree.n_; ++snode){
        Int fc = aXsuper[snode-1];
        Int lc = aXsuper[snode]-1;
        Int parent_col = this->PostParent(lc-1);
        Int parentSnode = ( parent_col == 0) ? 0:aSupMembership[parent_col-1];

          newTree.parent_(snode-1) = parentSnode;
#ifdef _DEBUG_
          logfileptr->OFS()<<"parent of curSnode "<<snode<<" is "<<parentSnode<<std::endl;
#endif
    } 

    newTree.poparent_ = newTree.parent_;
    newTree.bIsPostOrdered_ = true;


    return newTree;

  }





void ETree::DeepestFirst(Ordering & aOrder)
  {
  assert(bIsPostOrdered_);









    std::vector<Int> treesize(n_,0);
    std::vector<Int> depths(n_);
    //first, compute the depth of each node
    for(Int col=n_; col>=1; --col){
      Int parent = PostParent(col-1);
      if(parent==0){
        depths[col-1]=0;
      }
      else{
        treesize[parent-1]++;
        depths[col-1]=depths[parent-1]+1;
      }
    }


    for(Int col=n_; col>=1; --col){
      Int parent = PostParent(col-1);
      if(parent!=0){
        depths[parent-1]=max(depths[col-1],depths[parent-1]);
      }
    }




    //relabel the nodes within a subtree based on their depths

  }


}

