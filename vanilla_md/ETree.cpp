/// @file etree.cpp
/// @brief Implementation of the elimination-tree related algorithms.
/// @author Mathias Jacquelin
/// @date 2012-08-31
#include "ETree.hpp"

#include <iostream>

using namespace std;

  void DisjointSet::Initialize(int n){
    pp_.resize(n,I_ZERO);
    root_.resize(n);
  }

  void DisjointSet::Finalize(){
    pp_.resize(0);
  }




  ETree::ETree(){
    bIsPostOrdered_=false;
  }

  void ETree::BTreeToPO(vector<int> & fson, vector<int> & brother){
      //Do a depth first search to construct the postordered tree
      vector<int> stack(n_);
      postNumber_.resize(n_);
      invPostNumber_.resize(n_);

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

          postNumber_[vertex-1] = m;
          invPostNumber_[m-1] = vertex;

          vertex = brother[vertex-1];
        }

        if(exit){
          break;
        }
      }
}


  void ETree::PostOrderTree(){
    if(n_>0 && !bIsPostOrdered_){

      vector<int> fson(n_,I_ZERO);
      vector<int> brother(n_,I_ZERO);

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


      BTreeToPO(fson,brother);

      //      postParent_.resize(n_);
      //modify the parent list ?
      // node i is now node postNumber(i-1)
            for(int i=1; i<=n_;i++){
              int nunode = postNumber_[i-1];
              int ndpar = parent_[i-1];
              if(ndpar>0){
                ndpar = postNumber_[ndpar-1];
              }
              brother[nunode-1] = ndpar;
            }


      bIsPostOrdered_ = true;
    }

  }


  vector<int> ETree::SortChildren(vector<int> & cc){
    if(!bIsPostOrdered_){
      this->PostOrderTree();
    }

      vector<int> fson(n_,I_ZERO);
      vector<int> brother(n_,I_ZERO);
      vector<int> lson(n_,I_ZERO);

      //Get Binary tree representation
      int lroot = n_;
      for(int vertex=n_-1; vertex>0; vertex--){
        int curParent = PostParent(vertex-1);
        //int curParent = parent_(vertex-1);
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


      vector<int> perm;
//      vector<int> invperm;

      //Compute the parent permutation and update postNumber_
      //Do a depth first search to construct the postordered tree
      vector<int> stack(n_);
      perm.resize(n_);
//      invperm.resize(n_);

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

          perm[vertex-1] = m;
//          invperm(m-1) = vertex;

          vertex = brother[vertex-1];
        }

        if(exit){
          break;
        }
      }

      //Permute CC     
      for(int node = 1; node <= n_; ++node){
        int nunode = perm[node-1];
        stack[nunode-1] = cc[node-1];
      }

      for(int node = 1; node <= n_; ++node){
        cc[node-1] = stack[node-1];
      }

    return perm;

  }

  void ETree::PermuteTree(vector<int> & perm){
      //Compose the two permutations
      for(int i = 1; i <= n_; ++i){
            int interm = postNumber_[i-1];
            postNumber_[i-1] = perm[interm-1];
      }
      for(int i = 1; i <= n_; ++i){
        int node = postNumber_[i-1];
        invPostNumber_[node-1] = i;
      } 
  }






  void ETree::ConstructETree(int n, int * xadj, int * adj){


//    n_ = n;
//
//
//    parent_.resize(n_,I_ZERO);
//
//    DisjointSet sets;
//    sets.Initialize(n_);
//
//    int cset,croot,rset,rroot,row;
//
//
//    for (int col = 1; col <= n_; col++) {
//      parent_[col-1]=col; //1 based indexes
//      cset = sets.makeSet (col);
//      sets.Root(cset-1) = col;
//      parent_[col-1] = 0; 
//    }
//
//
//
//
//
//    for (int col = 1; col <= n_; col++) {
//      parent_[col-1]=col; //1 based indexes
//      cset = sets.makeSet (col);
//      sets.Root(cset-1) = col;
//      parent_[col-1] = 0; 
//
//      for (int p = xadj[col-1]; p < xadj[col]; p++) {
//        row = adj[p-1];
//
//
//
//        if (row >= col) continue;
//
//        rset = sets.find(row);
//        rroot = sets.Root(rset-1);
//
//        if (rroot != col) {
//          parent_[rroot-1] = col;
//          cset = sets.link(cset, rset);
//          sets.Root(cset-1) = col;
//        }
//      }
//
//    }
//
//
//    parent_[n_-1] = 0;





    bIsPostOrdered_=false;
    n_ = n;


    parent_.assign(n_,I_ZERO);


    vector<int> ancstr(n_);



    for(int i = 1; i<=n_; ++i){
            parent_[i-1] = 0;
            ancstr[i-1] = 0;
            int node = i;

            int jstrt = xadj[node-1];
            int jstop = xadj[node] - 1;
            if  ( jstrt < jstop ){
              for(int j = jstrt; j<=jstop; ++j){
                    int nbr = adj[j-1];
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


  void ETree::Dump(){
    cout<<"Post ordered etree: ";
    for(int i = 1;i<=Size();++i){
      cout<<" "<<PostParent(i-1);
    }
    cout<<endl;
  }







