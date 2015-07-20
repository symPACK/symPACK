/// @file etree.cpp
/// @brief Implementation of the elimination-tree related algorithms.
/// @author Mathias Jacquelin
/// @date 2012-08-31
#include "ngchol/ETree.hpp"

namespace LIBCHOLESKY{


  void DisjointSet::Initialize(int n){
    pp_.resize(n,0);
    root_.resize(n);
  }

  void DisjointSet::Finalize(){
    pp_.resize(0);
  }




  ETree::ETree(){
    bIsPostOrdered_=false;
  }

  ETree::ETree(SparseMatrixStructure & aGlobal, Ordering & aOrder){
    bIsPostOrdered_ = false;
    ConstructETree(aGlobal,aOrder);
  }



  void ETree::BTreeToPO(vector<int> & fson, vector<int> & brother, vector<int> & invpos){
      //Do a depth first search to construct the postordered tree
      vector<int> stack(n_);
      invpos.resize(n_);

      int stacktop=0, vertex=n_,m=0;
      bool exit = false;
      while( m<n_){
        do{
          stacktop++;
          stack[stacktop-1] = vertex;
          vertex = fson[vertex-1];
        }while(vertex>0);

        while(vertex==0){

          if(stacktop<=0){
            exit = true;
            break;
          }
          vertex = stack[stacktop-1];
          stacktop--;
          m++;

          invpos[vertex-1] = m;

          vertex = brother[vertex-1];
        }

        if(exit){
          break;
        }
      }
}


  void ETree::PostOrderTree(Ordering & aOrder){
    if(n_>0 && !bIsPostOrdered_){


      vector<int> fson(n_,0);
      vector<int> & brother = poparent_;
      brother.assign(n_,0);

      int lroot = n_;
      for(int vertex=n_-1; vertex>0; vertex--){
        int curParent = parent_[vertex-1];
        if(curParent==0 || curParent == vertex){
          brother[lroot-1] = vertex;
          lroot = vertex;
        }
        else{
          brother[vertex-1] = fson[curParent-1];
          fson[curParent-1] = vertex;
        }
      }

      vector<int> invpos;
      BTreeToPO(fson,brother,invpos);

      //modify the parent list ?
      // node i is now node invpos(i-1)
      for(int i=1; i<=n_;i++){
        int nunode = invpos[i-1];
        int ndpar = parent_[i-1];
        if(ndpar>0){
          ndpar = invpos[ndpar-1];
        }
        poparent_[nunode-1] = ndpar;
      }

      //we need to compose aOrder.invp and invpos
      aOrder.Compose(invpos);
      


      bIsPostOrdered_ = true;



    }

  }


  void ETree::SortChildren(vector<int> & cc, Ordering & aOrder){
    if(!bIsPostOrdered_){
      this->PostOrderTree(aOrder);
    }

      vector<int> fson(n_,0);
      vector<int> brother(n_,0);
      vector<int> lson(n_,0);

      //Get Binary tree representation
      int lroot = n_;
      for(int vertex=n_-1; vertex>0; vertex--){
        int curParent = PostParent(vertex-1);
        //int curParent = parent_[vertex-1];
        if(curParent==0 || curParent == vertex){
          brother[lroot-1] = vertex;
          lroot = vertex;
        }
        else{
          int ndlson = lson[curParent-1];
          if(ndlson > 0){
             if  ( cc[vertex-1] >= cc[ndlson-1] ) {
             //if  ( cc(ToPostOrder(vertex)-1) >= cc(ToPostOrder(ndlson)-1) ) {
                brother[vertex-1] = fson[curParent-1];
                fson[curParent-1] = vertex;
             }
             else{                                                                                                                                            
                brother[ndlson-1] = vertex;
                lson[curParent-1] = vertex;
             }                                                                                                                                                               
          }
          else{
             fson[curParent-1] = vertex;
             lson[curParent-1] = vertex;
          }
        }
      }
      brother[lroot-1]=0;


      //Compute the parent permutation and update postNumber_
      //Do a depth first search to construct the postordered tree
      vector<int> invpos(n_);
      vector<int> stack(n_);

      int stacktop=0, vertex=n_,m=0;
      bool exit = false;
      while( m<n_){
        do{
          stacktop++;
          stack[stacktop-1] = vertex;
          vertex = fson[vertex-1];
        }while(vertex>0);

        while(vertex==0){

          if(stacktop<=0){
            exit = true;
            break;
          }
          vertex = stack[stacktop-1];
          stacktop--;
          m++;

          invpos[vertex-1] = m;

          vertex = brother[vertex-1];
        }

        if(exit){
          break;
        }
      }


      //modify the parent list ?
      // node i is now node invpos(i-1)
      for(int i=1; i<=n_;i++){
        int nunode = invpos[i-1];
        int ndpar = poparent_[i-1];
        if(ndpar>0){
          ndpar = invpos[ndpar-1];
        }
        brother[nunode-1] = ndpar;
      }
      poparent_ = brother;




      //Permute CC     
      for(int node = 1; node <= n_; ++node){
        int nunode = invpos[node-1];
        stack[nunode-1] = cc[node-1];
      }

      for(int node = 1; node <= n_; ++node){
        cc[node-1] = stack[node-1];
      }

      //Compose the two permutations
      aOrder.Compose(invpos);

  }


  void ETree::ConstructETree(SparseMatrixStructure & aGlobal, Ordering & aOrder){
    bIsPostOrdered_=false;
    n_ = aGlobal.size;

    //Expand to symmetric storage
    aGlobal.ExpandSymmetric();

    parent_.resize(n_,0);


    vector<int> ancstr(n_);



    for(int i = 1; i<=n_; ++i){
            parent_[i-1] = 0;
            ancstr[i-1] = 0;
            int node = aOrder.perm[i-1];

            int jstrt = aGlobal.expColptr[node-1];
            int jstop = aGlobal.expColptr[node] - 1;
            if  ( jstrt < jstop ){
              for(int j = jstrt; j<=jstop; ++j){
                    int nbr = aGlobal.expRowind[j-1];
                    nbr = aOrder.invp[nbr-1];
                    if  ( nbr < i ){
//                       -------------------------------------------
//                       for each nbr, find the root of its current
//                       elimination tree.  perform path compression
//                       as the subtree is traversed.
//                       -------------------------------------------
                      int break_loop = 0;
                      if  ( ancstr[nbr-1] == i ){
                        break_loop = 1;
                      }
                      else{
                        while(ancstr[nbr-1] >0){
                          if  ( ancstr[nbr-1] == i ){
                            break_loop = 1;
                            break;
                          }
                          int next = ancstr[nbr-1];
                          ancstr[nbr-1] = i;
                          nbr = next;
                        }
                        //                       --------------------------------------------
                        //                       now, nbr is the root of the subtree.  make i
                        //                       the parent node of this root.
                        //                       --------------------------------------------
                        if(!break_loop){
                          parent_[nbr-1] = i;
                          ancstr[nbr-1] = i;
                        }
                      }
                    }
              }
            }
  }

  }





  ETree ETree::ToSupernodalETree(vector<int> & aXsuper,vector<int> & aSupMembership,Ordering & aOrder) const{
    ETree newTree;
    newTree.n_ = aXsuper.size()-1;
    newTree.parent_.resize(aXsuper.size()-1);
    


    for(int snode=1; snode<=newTree.n_; ++snode){
        int lc = aXsuper[snode]-1;
        int parent_col = this->PostParent(lc-1);
        int parentSnode = ( parent_col == 0) ? 0:aSupMembership[parent_col-1];

          newTree.parent_[snode-1] = parentSnode;
    } 


    return newTree;

  }





void ETree::DeepestFirst(Ordering & aOrder)
  {
    std::vector<int> treesize(n_,0);
    std::vector<int> depths(n_);
    //first, compute the depth of each node
    for(int col=n_; col>=1; --col){
      int parent = PostParent(col-1);
      if(parent==0){
        depths[col-1]=0;
      }
      else{
        treesize[parent-1]++;
        depths[col-1]=depths[parent-1]+1;
      }
    }


    for(int col=n_; col>=1; --col){
      int parent = PostParent(col-1);
      if(parent!=0){
        depths[parent-1]=max(depths[col-1],depths[parent-1]);
      }
    }

    //relabel the nodes within a subtree based on their depths

  }
}

