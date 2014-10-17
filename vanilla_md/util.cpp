#include "util.h"


#include <iostream>
#include <algorithm>
#include <numeric>

#include <fstream>
#include <string>
#include <sstream>
#include <iterator>

#include <assert.h>

//#define _verbose_

void displayMatrix(vector<int> & xadj, vector<int> & adj){
  for(int i = 1; i< xadj.size();++i){
    int fi = xadj[i-1];
    int li = xadj[i]-1;

    int prevrow = 0;
    int diagpassed = 0;
    for(int rowi = fi;rowi<=li;++rowi){
      int row = adj[rowi-1];
      if(row == i)
        continue;

      if(row>i && !diagpassed){
        diagpassed=1;
        for(int prow = prevrow+1; prow<i; prow++){ 
          cout<<0<<" ";
        }
        cout<<1<<" ";
        prevrow = i;
      }

      for(int prow = prevrow+1; prow<row; prow++){ 
        cout<<0<<" ";
      }
      cout<<1<<" ";
      prevrow = row;
    }

    if(!diagpassed){
      for(int prow = prevrow+1; prow<i; prow++){ cout<<0<<" ";}
      cout<<1<<" ";
      for(int prow = i+1; prow<xadj.size(); prow++){ cout<<0<<" ";}
    }
    else{
      for(int prow = prevrow+1; prow<xadj.size(); prow++){ cout<<0<<" ";}
    }
    cout<<endl;
  }
}

void GetPermutedGraph(int n, int nnz, int * xadj, int * adj, int * perm, int * invp, int * newxadj, int * newadj){
  //sort the adjacency structure according to the new labelling

  vector<int> invperm;
  if(invp==NULL){
    invperm.resize(n);
    invp = &invperm[0];
  }

  for(int step = 1; step<=n;++step){
    invp[perm[step-1]-1] = step;
  }




  newxadj[0] = 1;
  for(int step = 1; step<=n;++step){
    newxadj[step] = newxadj[step-1] + xadj[perm[step-1]]-xadj[perm[step-1]-1];
  }

#ifdef _verbose_
  cout<<"xadj: ";
  for(int step = 1; step<=n+1;++step){
    cout<<" "<<xadj[step-1];
  }
  cout<<endl;

  cout<<"newxadj: ";
  for(int step = 1; step<=n+1;++step){
    cout<<" "<<newxadj[step-1];
  }
  cout<<endl;
#endif

  for(int step = 1; step<=n;++step){
    int fn = newxadj[step-1];
    int ln = newxadj[step]-1;

    int oldlabel = perm[step-1];
    int ofn = xadj[oldlabel-1];
    int oln = xadj[oldlabel]-1;
    for(int i =ofn;i<=oln;++i){
      newadj[fn+i-ofn-1] = invp[adj[i-1]-1];
    }

    //now sort them
    sort(&newadj[fn-1],&newadj[ln]);

  }
  //  newadj.back() = 0;




#ifdef _verbose_
  cout<<"adj: ";
  for(int step = 1; step<=nnz;++step){
    cout<<" "<<adj[step-1];
  }
  cout<<endl;

  cout<<"newadj: ";
  for(int step = 1; step<=nnz;++step){
    cout<<" "<<newadj[step-1];
  }
  cout<<endl;
#endif
}


void GetPermutedGraph(int n, int nnz, int * xadj, int * adj, int * perm, int * newxadj, int * newadj){
  GetPermutedGraph(n, nnz, xadj, adj, perm, NULL, newxadj, newadj);
}


void SymbolicFactorization(ETree& tree,const vector<int> & colptr,const vector<int> & rowind,const vector<int> & cc, vector<int> & xlindx, vector<int> & lindx){

  int size = tree.Size();
  int nsuper = tree.Size();

  int nzbeg = 0;
  //nzend points to the last used slot in lindx
  int nzend = 0;

  //tail is the end of list indicator (in rchlnk, not mrglnk)
  int tail = size +1;

  int head = 0;

  //Array of length nsuper containing the children of 
  //each supernode as a linked list
  vector<int> mrglnk(nsuper,0);

  //Array of length n+1 containing the current linked list 
  //of merged indices (the "reach" set)
  vector<int> rchlnk(size+1);

  //Array of length n used to mark indices as they are introduced
  // into each supernode's index set
  vector<int> marker(size,0);


  xlindx.resize(nsuper+1);

  //Compute the sum of the column count and resize lindx accordingly
  int nofsub = 1;
  for(int i =0; i<cc.size();++i){
    nofsub+=cc[i];
  }

  lindx.resize(nofsub);


  int point = 1;
  for(int ksup = 1; ksup<=nsuper; ++ksup){
    int fstcol = ksup;
    xlindx[ksup-1] = point;
    point += cc[fstcol-1]; 
  } 
  xlindx[nsuper] = point;

  for(int ksup = 1; ksup<=nsuper; ++ksup){
    int fstcol = ksup;
    int lstcol = ksup;
    int width = lstcol - fstcol +1;
    int length = cc[fstcol-1];
    int knz = 0;
    rchlnk[head] = tail;
    int jsup = mrglnk[ksup-1];

    //If ksup has children in the supernodal e-tree
    if(jsup>0){
      //copy the indices of the first child jsup into 
      //the linked list, and mark each with the value 
      //ksup.
      int jwidth = 1;
      int jnzbeg = xlindx[jsup-1] + jwidth;
      int jnzend = xlindx[jsup] -1;
      for(int jptr = jnzend; jptr>=jnzbeg; --jptr){
        int newi = lindx[jptr-1];
        ++knz;
        marker[newi-1] = ksup;
        rchlnk[newi] = rchlnk[head];
        rchlnk[head] = newi;
      }

      //for each subsequent child jsup of ksup ...
      jsup = mrglnk[jsup-1];
      while(jsup!=0 && knz < length){
        //merge the indices of jsup into the list,
        //and mark new indices with value ksup.

        jwidth = 1;
        jnzbeg = xlindx[jsup-1] + jwidth;
        jnzend = xlindx[jsup] -1;
        int nexti = head;
        for(int jptr = jnzbeg; jptr<=jnzend; ++jptr){
          int newi = lindx[jptr-1];
          int i;
          do{
            i = nexti;
            nexti = rchlnk[i];
          }while(newi > nexti);

          if(newi < nexti){
            ++knz;
            rchlnk[i] = newi;
            rchlnk[newi] = nexti;
            marker[newi-1] = ksup;
            nexti = newi;
          }
        }
        jsup = mrglnk[jsup-1];
      }
    }
    //structure of a(*,fstcol) has not been examined yet.  
    //"sort" its structure into the linked list,
    //inserting only those indices not already in the
    //list.
    if(knz < length){
      int node = tree.FromPostOrder(fstcol);
      int knzbeg = colptr[node-1];
      int knzend = colptr[node]-1;
      for(int kptr = knzbeg; kptr<=knzend;++kptr){
        int newi = rowind[kptr-1];
        newi = tree.ToPostOrder(newi);
        if(newi > fstcol && marker[newi-1] != ksup){
          //position and insert newi in list and
          // mark it with kcol
          int nexti = head;
          int i;
          do{
            i = nexti;
            nexti = rchlnk[i];
          }while(newi > nexti);
          ++knz;
          rchlnk[i] = newi;
          rchlnk[newi] = nexti;
          marker[newi-1] = ksup;
        }
      }
    }

    //if ksup has no children, insert fstcol into the linked list.
    if(rchlnk[head] != fstcol){
      rchlnk[fstcol] = rchlnk[head];
      rchlnk[head] = fstcol;
      ++knz;
    }

    //assert(knz == cc[fstcol-1]);


    //copy indices from linked list into lindx(*).
    nzbeg = nzend+1;
    nzend += knz;
    //      assert(nzend+1 == xlindx[ksup]);
    int i = head;
    for(int kptr = nzbeg; kptr<=nzend;++kptr){
      i = rchlnk[i];
      lindx[kptr-1] = i;
    } 

    //if ksup has a parent, insert ksup into its parent's 
    //"merge" list.
    if(length > width){
      int pcol = lindx[xlindx[ksup-1] + width -1];
      int psup = pcol;
      mrglnk[ksup-1] = mrglnk[psup-1];
      mrglnk[psup-1] = ksup;
    }
  }

  lindx.resize(nzend+1);

}


int ReadAdjacency(const char * pfilename, vector<int> & xadj, vector<int> & adj){
  string filename(pfilename);
  ifstream infile;
  infile.open(filename.c_str());

  adj.resize(0);
  xadj.resize(0);
  vector<int> ixadj;

  string line;
  //read xadj on the first line of the input file
  if(getline(infile, line))
  {
    istringstream iss(line);
    int i;
    while(iss>> i){
       ixadj.push_back(i);
    }
  }    
  else{
    return -2;
  }

  //read adj on the second line of the input file
  if(getline(infile, line))
  {
    istringstream iss(line);
    int i;
    int col = 1;
    int xpos = 0;
    int pos = 0;
    int offset = 0;
    int ifound =0;
    xadj.push_back(1);
    while(iss>> i){
       pos++;
       if(pos>=ixadj[col]){
        if(!ifound){
          adj.push_back(col);
        }
        col++;
        ifound=0; 
        xadj.push_back(adj.size()+1);
       }
       if(i==col){
        ifound=1;
       }
       adj.push_back(i);
    }
    if(!ifound){
      adj.push_back(col);
    }
    xadj.push_back(adj.size()+1);
  }    
  else{
    return -2;
  }


  infile.close();



#ifdef _verbose_
  cout<<"ixadj: ";
  for(int i = 0;i<ixadj.size();++i){
    cout<<" "<<ixadj[i];
  }
  cout<<endl;

  cout<<"xadj: ";
  for(int i = 0;i<xadj.size();++i){
    cout<<" "<<xadj[i];
  }
  cout<<endl;
  cout<<"adj: ";
  for(int i = 0;i<adj.size();++i){
    cout<<" "<<adj[i];
  }
  cout<<endl;
#endif

  return 0;

}



int ReadAdjacency(const char * pfilename, int ** pxadj, int ** padj, int * n , int * nnz){
  ifstream infile;
  infile.open(pfilename);

  


  vector<int> iixadj;
  vector<int> ixadj;
  vector<int> iadj;

  string line;
  //read xadj on the first line of the input file
  if(getline(infile, line))
  {
    istringstream iss(line);
    int pos = 0;
    int i;
    while(iss>> i){
       iixadj.push_back(i);
    }
  }    
  else{
    return -2;
  }

  //read adj on the second line of the input file
  if(getline(infile, line))
  {
    istringstream iss(line);
    int i;
    int col = 1;
    int xpos = 0;
    int pos = 0;
    int offset = 0;
    int ifound =0;
    ixadj.push_back(1);
    while(iss>> i){
       pos++;
       if(pos>=iixadj[col]){
        if(!ifound){
          iadj.push_back(col);
        }
        col++;
        ifound=0; 
        ixadj.push_back(iadj.size()+1);
       }
       if(i==col){
        ifound=1;
       }
       iadj.push_back(i);
    }
    if(!ifound){
      iadj.push_back(col);
    }
    ixadj.push_back(iadj.size()+1);
  }    
  else{
    return -2;
  }


  infile.close();

  *n = ixadj.size()-1;


  //expand to asymmetric storage
  vector<int> xadj;
  vector<int> adj;
  ExpandSymmetric(*n,&ixadj[0],&iadj[0], xadj, adj);


  *pxadj = (int*)malloc(xadj.size()*sizeof(int));
  *padj = (int*)malloc(adj.size()*sizeof(int));
  std::copy(xadj.begin(),xadj.end(),*pxadj);
  std::copy(adj.begin(),adj.end(),*padj);

  *nnz = adj.size();

#ifdef _verbose_
  cout<<"ixadj: ";
  for(int i = 0;i<ixadj.size();++i){
    cout<<" "<<ixadj[i];
  }
  cout<<endl;

  cout<<"xadj: ";
  for(int i = 0;i<xadj.size();++i){
    cout<<" "<<xadj[i];
  }
  cout<<endl;
  cout<<"adj: ";
  for(int i = 0;i<adj.size();++i){
    cout<<" "<<adj[i];
  }
  cout<<endl;
#endif

  return 0;


}


void ExpandSymmetric(int size,const int * colptr,const int * rowind, vector<int> & expColptr, vector<int> & expRowind){
  //code from sparsematrixconverter
  /* set-up */
  vector<int> cur_col_nnz(size);
  vector<int> new_col_nnz(size);

  /*
   * Scan A and count how many new non-zeros we'll need to create.
   *
   * Post:
   *   cur_col_nnz[i] == # of non-zeros in col i of the original symmetric 
   *                     matrix.
   *   new_col_nnz[i] == # of non-zeros to be stored in col i of the final 
   *                     expanded matrix.
   *   new_nnz == total # of non-zeros to be stored in the final expanded 
   *              matrix.
   */
  int new_nnz = 0; 
  int orig_nnz = 0; 
  for (int i = 0; i < size; i++) 
  { 
    //int * iptr = find(&rowind[colptr[i]-1],&rowind[colptr[i+1]-1],i+1);
    cur_col_nnz[i] = colptr[i+1] - colptr[i];
    new_col_nnz[i] = cur_col_nnz[i];
    new_nnz += new_col_nnz[i];
  }    

  orig_nnz = new_nnz;

  for (int i = 0; i < size; i++) 
  {    
    int k;
    for (k = colptr[i]; k < colptr[i+1]; k++) 
    {
      int j = rowind[k-1]-1;
      if(j<i){
        //the matrix is already in asymm form, copy and exit
        expColptr.resize(size+1);
        expRowind.resize(orig_nnz);
    
        for (int i = 0; i <= size; i++){
          expColptr[i] = colptr[i];
        }
        
        for (int i = 0; i < orig_nnz; i++){
          expRowind[i] = rowind[i];
        }

        return;
      }
      if (j > i)
      {
        new_col_nnz[j]++;
        new_nnz++;
      }
    }
  }

  /*
   *  Initialize row pointers in expanded matrix.
   *
   *  Post:
   *    new_colptr initialized to the correct, final values.
   *    new_col_nnz[i] reset to be equal to cur_col_nnz[i].
   */
  expColptr.resize(size+1);
  expColptr[0] = 1;
  for (int i = 1; i <= size; i++)
  {
    expColptr[i] = expColptr[i-1] + new_col_nnz[i-1];
    new_col_nnz[i-1] = cur_col_nnz[i-1];
  }
  expColptr[size] = new_nnz+1;

  expRowind.resize(new_nnz);

  /*
   *  Complete expansion of A to full storage.
   *
   *  Post:
   *    (new_colptr, new_rowind, new_values) is the full-storage equivalent of A.
   *    new_col_nnz[i] == # of non-zeros in col i of A.
   */




  for (int i = 0; i < size; i++)
  {
    int cur_nnz = cur_col_nnz[i];
    int k_cur   = colptr[i] -1;
    int k_new   = expColptr[i] -1;


    /* copy current non-zeros from old matrix to new matrix */
    std::copy(&rowind[k_cur], &rowind[k_cur + cur_nnz] , &expRowind[k_new]);

    /* fill in the symmetric "missing" values */
    while (k_cur < colptr[i+1]-1)
    {
      /* non-zero of original matrix */
      int j = rowind[k_cur]-1;

      if (j != i)  /* if not a non-diagonal element ... */
      {
        /* position of this transposed element in new matrix */
        k_new = expColptr[j]-1 + new_col_nnz[j];

        /* store */
        expRowind[k_new] = i+1;
        /*  update so next element stored at row j will appear
         *  at the right place.  */
        new_col_nnz[j]++;
      }

      k_cur++;
    }

    //sort
    std::sort(&expRowind[expColptr[i]-1],&expRowind[expColptr[i+1]-1]);
  }
}










void GetLColRowCount(ETree & tree,const int * pxadj, const int * padj, vector<int> & cc, vector<int> & rc);

void GetLColRowCount(ETree & tree,const int * pxadj, const int * padj, vector<int> & cc, vector<int> & rc){
  //The tree need to be postordered
  if(!tree.IsPostOrdered()){
    tree.PostOrderTree();
  }

  vector<int> xadj, adj;

  int size = tree.Size();
  ExpandSymmetric(size, pxadj,padj,xadj,adj);



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

void GetPrefixSum(int n, int * arr, int * arrout){
  partial_sum(arr,arr+n,arrout);
}

//nnz must be adj.size()
double GetCost(int n, int nnz, int * xadj, int * adj,int * perm){

  int initEdgeCnt;

  vector<int> newxadj(n+1);
  vector<int> newadj(nnz);
  GetPermutedGraph(n,nnz,xadj, adj, perm, NULL, &newxadj[0], &newadj[0]);


  initEdgeCnt=n;
  //initialize nodes
  for(int i=0;i<n;++i){
    for(int idx = newxadj[i]; idx <= newxadj[i+1]-1;++idx){
      if(newadj[idx-1]-1>i){
        initEdgeCnt++;
      }
    }
  }

#ifdef _verbose_
  cout<<"Initial edge count: "<<initEdgeCnt<<endl;
#endif


  //Get the elimination tree
  ETree  tree;
  tree.ConstructETree(n,&newxadj[0],&newadj[0]);
  tree.PostOrderTree();

  vector<int> cc,rc;
  GetLColRowCount(tree,&newxadj[0],&newadj[0],cc,rc);


  //sum counts all the diagonal elements
  int sum = 0;
    for(int i =0; i<cc.size(); ++i){
      sum+= cc[i];
    }

  #ifdef _verbose_
  cout<<"Column count: ";
  for(int i =0; i<cc.size(); ++i){
    cout<<" "<<cc[i];
  }
  cout<<endl;
  #endif

#ifdef _verbose_
  cout<<"Sum is "<<sum<<endl;
#endif

  return (double)(sum - initEdgeCnt);

}

#if 0
double GetCostPerCol(int n, int nnz, int * xadj, int * adj,int * perm, int * costc){
  int initEdgeCnt;



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



  vector<int> newxadj(n+1);
  vector<int> newadj(nnz);
  GetPermutedGraph(n,nnz,xadj, adj, perm, &newxadj[0], &newadj[0]);


  initEdgeCnt=n;
  //initialize nodes
  for(int i=0;i<n;++i){
    if(costc!=NULL){
      costc[i] = 1;
    }
    for(int idx = newxadj[i]; idx <= newxadj[i+1]-1;++idx){
      if(newadj[idx-1]-1>i){

        if(costc!=NULL){
          costc[i]++;
        }
        initEdgeCnt++;
      }
    }
  }

#ifdef _verbose_
  cout<<"Initial edge count: "<<initEdgeCnt<<endl;
#endif


  //Get the elimination tree
  ETree  tree;
  tree.ConstructETree(n,&newxadj[0],&newadj[0]);
  tree.PostOrderTree();
  vector<int> poperm(n);
  for(int i=1;i<=n;++i){poperm[i-1]=tree.ToPostOrder(i);}


  vector<int> cc2,rc;
  GetLColRowCount(tree,&newxadj[0],&newadj[0],cc2,rc);
#ifdef _verbose_
  cout<<"Column count (PO): ";
  for(int i =0; i<cc2.size(); ++i){
    cout<<" "<<cc2[i];
  }
  cout<<endl;
#endif

  vector<int> xlindx,lindx;
  SymbolicFactorization(tree,newxadj,newadj,cc2, xlindx, lindx);


  //Expand lindx and xlindx to assyn storage
  vector<int> expxlindx,explindx;
  ExpandSymmetric(xlindx.size()-1,&xlindx[0],&lindx[0], expxlindx, explindx);



#ifdef _verbose_
  for(int i = 1; i<=n; ++i){
    int fi = expxlindx[i-1];
    int li = expxlindx[i]-1;
    int check = 0;
    for(int rowi=fi;rowi<=li;++rowi){
      int row = explindx[rowi-1];
      check++;
    }
//    cerr<<check<<" vs "<<cc2[i-1]+rc[i-1]-1<<endl;
    assert(check==cc2[i-1]+rc[i-1]-1);
  }
#endif




#ifdef _verbose_
  cout<<endl;
  cout<<"Expanded L PO factor"<<endl;
  displayMatrix(expxlindx,explindx);
  cout<<endl;
  cout<<endl;

  cout<<"lindx: ";
  for(int i =0; i<lindx.size(); ++i){
    cout<<" "<<lindx[i];
  }
  cout<<endl;

  cout<<"explindx: ";
  for(int i =0; i<explindx.size(); ++i){
    cout<<" "<<explindx[i];
  }
  cout<<endl;
#endif


  //Permute expxlindx and explindx back to the original ordering

  vector<int> fxadj(expxlindx.size());
  vector<int> fadj(explindx.size());
  GetPermutedGraph(expxlindx.size()-1,explindx.size(),&expxlindx[0], &explindx[0], &poperm[0], &fxadj[0], &fadj[0]);

#ifdef _verbose_
  tree.Dump();
  cout<<"poperm: ";
  for(int i =0; i<poperm.size(); ++i){
    cout<<" "<<poperm[i];
  }
  cout<<endl;

  cout<<endl;
  cout<<endl;
  cout<<"Expanded reordered L factor"<<endl;
  displayMatrix(fxadj,fadj);
  cout<<endl;
  cout<<endl;

  cout<<endl;
  cout<<endl;
  cout<<"Expanded original matrix"<<endl;
  displayMatrix(newxadj,newadj);
  cout<<endl;
  cout<<endl;
#endif


  vector<int> cc(cc2.size());
  //count the sub diagonal non zeros for each column
  for(int i = 1; i<=n; ++i){
    int fi = fxadj[i-1];
    int li = fxadj[i]-1;
    cc[i-1]=0;
    int check = 0;
    for(int rowi=fi;rowi<=li;++rowi){
      int row = fadj[rowi-1];
      if(row>=i){
        ++cc[i-1];
      }
      check++;
    }
#ifdef _verbose_
    assert(check==cc2[tree.ToPostOrder(i)-1]+rc[tree.ToPostOrder(i)-1]-1);
#endif
  }


#ifdef _verbose_
  cout<<"Initial edge count: ";
  for(int i =0; i<cc.size(); ++i){
    cout<<" "<<costc[i];
  }
  cout<<endl;
#endif




  //sum counts all the diagonal elements
  int sum = 0;
  if(costc!=NULL){
    for(int i =0; i<cc.size(); ++i){
      costc[i] = cc[i] - costc[i];
      sum+= cc[i];
    }
  }
  else{
    for(int i =0; i<cc.size(); ++i){
      sum+= cc[i];
    }
  }

  #ifdef _verbose_
  cout<<"Column count: ";
  for(int i =0; i<cc.size(); ++i){
    cout<<" "<<cc[i];
  }
  cout<<endl;

  cout<<"Column costs: ";
  for(int i =0; i<cc.size(); ++i){
    cout<<" "<<costc[i];
  }
  cout<<endl;
  #endif

#ifdef _verbose_
  cout<<"Sum is "<<sum<<endl;
#endif

  return (double)(sum - initEdgeCnt);

}
#else
//In this version perm is replaced by its post ordered equivalent ordering
double GetCostPerCol(int n, int nnz, int * xadj, int * adj,int * perm, int * costc){
  int initEdgeCnt;

  vector<int> newxadj(n+1);
  vector<int> newadj(nnz);
  vector<int> invp(n);
  GetPermutedGraph(n,nnz,xadj, adj, perm,&invp[0], &newxadj[0], &newadj[0]);


  //Get the elimination tree
  ETree  tree;
  tree.ConstructETree(n,&newxadj[0],&newadj[0]);
  tree.PostOrderTree();


  //update perm
  vector<int> poinvp(n);
  for(int i=1;i<=n;++i){poinvp[i-1]=tree.ToPostOrder(i);}

#ifdef _verbose_
  tree.Dump();
  vector<int> poperm(n);
  for(int i=1;i<=n;++i){poperm[i-1]=tree.FromPostOrder(i);}
  cout<<"poperm: ";
  for(int i =0; i<n; ++i){
    cout<<" "<<poperm[i];
  }
  cout<<endl;

  cout<<"poinvp: ";
  for(int i =0; i<n; ++i){
    cout<<" "<<poinvp[i];
  }
  cout<<endl;
#endif


  //Compose the two permutations 
  for(int i = 1; i <= n; ++i){ 
    int interm = invp[i-1]; 
    invp[i-1] = poinvp[interm-1]; 
  } 
  for(int i = 1; i <= n; ++i){ 
    int node = invp[i-1]; 
    perm[node-1] = i; 
  } 



#ifdef _verbose_
  cout<<"perm: ";
  for(int i =0; i<n; ++i){
    cout<<" "<<perm[i];
  }
  cout<<endl;
#endif






  initEdgeCnt=n;
  //initialize nodes
  for(int i=0;i<n;++i){
    int col = tree.FromPostOrder(i+1);
    if(costc!=NULL){
      costc[i] = 1;
    }
    for(int idx = newxadj[col-1]; idx <= newxadj[col]-1;++idx){
      if(newadj[idx-1]>col){

        if(costc!=NULL){
          costc[i]++;
        }
        initEdgeCnt++;
      }
    }
  }

#ifdef _verbose_
  cout<<"Initial edge count: "<<initEdgeCnt<<endl;
#endif










  vector<int> cc,rc;
  GetLColRowCount(tree,&newxadj[0],&newadj[0],cc,rc);

  //sum counts all the diagonal elements
  int sum = 0;
  if(costc!=NULL){
    for(int i =0; i<cc.size(); ++i){
      costc[i] = cc[i] - costc[i];
      sum+= cc[i];
    }
  }
  else{
    for(int i =0; i<cc.size(); ++i){
      sum+= cc[i];
    }
  }

  #ifdef _verbose_
  cout<<"Column count: ";
  for(int i =0; i<cc.size(); ++i){
    cout<<" "<<cc[i];
  }
  cout<<endl;

  cout<<"Column costs: ";
  for(int i =0; i<cc.size(); ++i){
    cout<<" "<<costc[i];
  }
  cout<<endl;
  #endif

#ifdef _verbose_
  cout<<"Sum is "<<sum<<endl;
#endif

  return (double)(sum - initEdgeCnt);

}
#endif
