#include "cost.h"

#include <vector>
#include <iostream>
#include <algorithm>

#include "ETree.hpp"

using namespace std;

void GetLColRowCount(ETree & tree,const int * xadj, const int * adj, vector<int> & cc, vector<int> & rc);

  void GetLColRowCount(ETree & tree,const int * xadj, const int * adj, vector<int> & cc, vector<int> & rc){
     //The tree need to be postordered
    if(!tree.IsPostOrdered()){
      tree.PostOrderTree();
    }

    int size = tree.Size();
 
    cc.resize(size);
    rc.resize(size);

    vector<int> level(size+1);
    vector<int> weight(size+1);
    vector<int> fdesc(size+1);
    vector<int> nchild(size+1);
    vector<int> set(size);
    vector<int> prvlf(size);
    vector<int> prvnbr(size);


        int xsup = 1;
        level[0] = 0;
      for(int k = size; k>=1; --k){
            rc[k-1] = 1;
            cc[k-1] = 0;
            set[k-1] = k;
            prvlf[k-1] = 0;
            prvnbr[k-1] = 0;
            level[k] = level[tree.PostParent(k-1)] + 1;
            weight[k] = 1;
            fdesc[k] = k;
            nchild[k] = 0;
      }

      nchild[0] = 0;
      fdesc[0] = 0;
      for(int k =1; k<size; ++k){
            int parent = tree.PostParent(k-1);
            weight[parent] = 0;
            ++nchild[parent];
            int ifdesc = fdesc[k];
            if  ( ifdesc < fdesc[parent] ) {
                fdesc[parent] = ifdesc;
            }
      }







      for(int lownbr = 1; lownbr<=size; ++lownbr){
        int lflag = 0;
        int ifdesc = fdesc[lownbr];
        int oldnbr = tree.FromPostOrder(lownbr);
        int jstrt = xadj[oldnbr-1];
        int jstop = xadj[oldnbr] - 1;


        //           -----------------------------------------------
        //           for each ``high neighbor'', hinbr of lownbr ...
        //           -----------------------------------------------
        for(int j = jstrt; j<=jstop;++j){
          int hinbr = tree.ToPostOrder(adj[j-1]);
          if  ( hinbr > lownbr )  {
            if  ( ifdesc > prvnbr[hinbr-1] ) {
              //                       -------------------------
              //                       increment weight(lownbr).
              //                       -------------------------
              ++weight[lownbr];
              int pleaf = prvlf[hinbr-1];
              //                       -----------------------------------------
              //                       if hinbr has no previous ``low neighbor'' 
              //                       then ...
              //                       -----------------------------------------
              if  ( pleaf == 0 ) {
                //                           -----------------------------------------
                //                           ... accumulate lownbr-->hinbr path length 
                //                               in rowcnt(hinbr).
                //                           -----------------------------------------
                rc[hinbr-1] += level[lownbr] - level[hinbr];
              }
              else{
                //                           -----------------------------------------
                //                           ... otherwise, lca <-- find(pleaf), which 
                //                               is the least common ancestor of pleaf 
                //                               and lownbr.
                //                               (path halving.)
                //                           -----------------------------------------
                int last1 = pleaf;
                int last2 = set[last1-1];
                int lca = set[last2-1];
                while(lca != last2){
                  set[last1-1] = lca;
                  last1 = lca;
                  last2 = set[last1-1];
                  lca = set[last2-1];
                }
                //                           -------------------------------------
                //                           accumulate pleaf-->lca path length in 
                //                           rowcnt(hinbr).
                //                           decrement weight(lca).
                //                           -------------------------------------
                rc[hinbr-1] += level[lownbr] - level[lca];
                --weight[lca];
              }
              //                       ----------------------------------------------
              //                       lownbr now becomes ``previous leaf'' of hinbr.
              //                       ----------------------------------------------
              prvlf[hinbr-1] = lownbr;
              lflag = 1;
            }
            //                   --------------------------------------------------
            //                   lownbr now becomes ``previous neighbor'' of hinbr.
            //                   --------------------------------------------------
            prvnbr[hinbr-1] = lownbr;
          }
        }
        //           ----------------------------------------------------
        //           decrement weight ( parent(lownbr) ).
        //           set ( p(lownbr) ) <-- set ( p(lownbr) ) + set(xsup).
        //           ----------------------------------------------------
        int parent = tree.PostParent(lownbr-1);
        --weight[parent];


        //merge the sets
        if  ( lflag == 1  || nchild[lownbr] >= 2 ) {
          xsup = lownbr;
        }
        set[xsup-1] = parent;
      }

        for(int k = 1; k<=size; ++k){
            int temp = cc[k-1] + weight[k];
            cc[k-1] = temp;
            int parent = tree.PostParent(k-1);
            if  ( parent != 0 ) {
                cc[parent-1] += temp;
            }
        }

  }


//nnz must be adj.size()
double GetCost(int n, int nnz, int * xadj, int * adj,int * perm){
  int initEdgeCnt;

  vector<int> invperm(n);

  for(int step = 1; step<=n;++step){
    invperm[perm[step-1]-1] = step;
  }

//  cout<<"perm: ";
//  for(int step = 1; step<=n;++step){
//    cout<<" "<<perm[step-1];
//  }
//  cout<<endl;
//
//  cout<<"invperm: ";
//  for(int step = 1; step<=n;++step){
//    cout<<" "<<invperm[step-1];
//  }
//  cout<<endl;



  initEdgeCnt=0;
  //initialize nodes
  for(int i=0;i<n;++i){
    for(int idx = xadj[i]; idx <= xadj[i+1]-1;++idx){
      if(adj[idx-1]>i+1){
        initEdgeCnt++;
      }
    }
  }

//  cout<<"Initial edge count: "<<initEdgeCnt<<endl;

  //sort the adjacency structure according to the new labelling
  vector<int> newxadj(n+1);
  vector<int> newadj(nnz);

  newxadj[0] = 1;
  for(int step = 1; step<=n;++step){
    newxadj[step] = newxadj[step-1] + xadj[perm[step-1]]-xadj[perm[step-1]-1];
  }

//  cout<<"xadj: ";
//  for(int step = 1; step<=n+1;++step){
//    cout<<" "<<xadj[step-1];
//  }
//  cout<<endl;
//
//  cout<<"newxadj: ";
//  for(int step = 1; step<=n+1;++step){
//    cout<<" "<<newxadj[step-1];
//  }
//  cout<<endl;

  for(int step = 1; step<=n;++step){
    int fn = newxadj[step-1];
    int ln = newxadj[step]-1;

    int oldlabel = perm[step-1];
    int ofn = xadj[oldlabel-1];
    int oln = xadj[oldlabel]-1;
    for(int i =ofn;i<=oln;++i){
      newadj[fn+i-ofn-1] = invperm[adj[i-1]-1];
    }

    //now sort them
    sort(&newadj[fn-1],&newadj[ln]);

  }
//  newadj.back() = 0;
  
//  cout<<"adj: ";
//  for(int step = 1; step<=adj.size();++step){
//    cout<<" "<<adj[step-1];
//  }
//  cout<<endl;
//
//  cout<<"newadj: ";
//  for(int step = 1; step<=newadj.size();++step){
//    cout<<" "<<newadj[step-1];
//  }
//  cout<<endl;


  //Get the elimination tree
  ETree  tree;
  tree.ConstructETree(n,&newxadj[0],&newadj[0]);
//  tree.ConstructETree(n,&xadj[0],&adj[0]);

  vector<int> cc,rc;
  GetLColRowCount(tree,&newxadj[0],&newadj[0],cc,rc);
//  GetLColRowCount(tree,&xadj[0],&adj[0],cc,rc);

  //sum counts all the diagonal elements
  int sum = n;
  cout<<"Column count: ";
  for(int i =0; i<cc.size(); ++i){
    sum+= cc[i];
    cout<<" "<<cc[i];
  }
  cout<<endl;

//  cout<<"Sum is "<<sum<<endl;

  return (double)(sum - initEdgeCnt);

}
