#include "ngchol/SparseMatrixStructure.hpp"
#include "ngchol/ETree.hpp"
#include "ngchol/Utility.hpp"

#include <limits>       // std::numeric_limits

#include <iostream>
#include <iterator>
#include <set>
#include <list>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cassert>

namespace LIBCHOLESKY{

  SparseMatrixStructure::SparseMatrixStructure(){
    bIsGlobal = false;
    bIsExpanded = false;
  }

  void SparseMatrixStructure::ClearExpandedSymmetric(){
    expColptr.resize(0);
    expRowind.resize(0);
    bIsExpanded=false;
  }

  void SparseMatrixStructure::ExpandSymmetric(){
    if(!bIsExpanded){
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
      for (int i = 0; i < size; i++) 
      {    
        cur_col_nnz[i] = colptr[i+1] - colptr[i];
        new_col_nnz[i] = cur_col_nnz[i];
        new_nnz += new_col_nnz[i];
      }    

      for (int i = 0; i < size; i++) 
      {    
        int k;
        for (k = colptr[i]; k < colptr[i+1]; k++) 
        {
          int j = rowind[k-1]-1;
          if (j != i)
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
      
      for (int i = 0; i < size; i++){
        new_col_nnz[i] = 0;
      }
      //expColptr[size] = new_nnz;

      expRowind.assign(new_nnz,-1);

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
        int k_new;



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


        /* position of this transposed element in new matrix */
        k_new = expColptr[i] -1 + new_col_nnz[i];

        k_cur   = colptr[i] -1;
        /* copy current non-zeros from old matrix to new matrix */
        std::copy(&rowind[0] + k_cur, &rowind[0] + k_cur + cur_nnz , &expRowind[0] + k_new);
      }

      bIsExpanded =true;
    }

  }


  void SparseMatrixStructure::ToGlobal(SparseMatrixStructure & pGlobal,MPI_Comm & comm){
    // test if structure hasn't been allocated yet
    if(bIsGlobal){
      pGlobal = *this;
    }
    else{


      int np;
      int iam;

      int ismpi=0;
      MPI_Initialized( &ismpi);

      //FIXME needs to be passed as an argument ?
      //MPI_Comm comm = MPI_COMM_WORLD;

      int isnull= (comm == MPI_COMM_NULL);
      //    logfileptr->OFS()<<ismpi<<std::endl;
      //    logfileptr->OFS()<<comm<<std::endl;
      //    logfileptr->OFS()<<MPI_COMM_NULL<<std::endl;
      //    logfileptr->OFS()<<isnull<<std::endl;

      if(ismpi && isnull==0){
        MPI_Comm_size(comm,&np);
        MPI_Comm_rank(comm, &iam);
      }
      else{
        //throw an exception
        throw std::logic_error("MPI needs to be available.");
      }


      pGlobal.bIsGlobal = true;
      pGlobal.bIsExpanded = false;
      pGlobal.size = size;
      pGlobal.colptr.resize(size+1);


      /* Allgatherv for row indices. */ 
      vector<int> prevnz(np);
      vector<int> rcounts(np);
      MPI_Allgather(&nnz, 1, MPI_INT, &rcounts[0], 1, MPI_INT, comm);

      prevnz[0] = 0;
      for (int i = 0; i < np-1; ++i) { prevnz[i+1] = prevnz[i] + rcounts[i]; } 

      pGlobal.nnz = 0;
      for (int i = 0; i < np; ++i) { pGlobal.nnz += rcounts[i]; } 
      pGlobal.rowind.resize(pGlobal.nnz);

      MPI_Allgatherv(&rowind[0], nnz, MPI_INT, &pGlobal.rowind[0],&rcounts[0], &prevnz[0], MPI_INT, comm); 


      /* Allgatherv for colptr */
      // Compute the number of columns on each processor
      int numColFirst = std::max(1,size / np);
      rcounts.assign(rcounts.size(),numColFirst);
      rcounts[np-1] = size - numColFirst * (np-1);  // Modify the last entry     


      vector<int> rdispls(np);
      rdispls[0] = 0;
      for (int i = 0; i < np-1; ++i) { rdispls[i+1] = rdispls[i] + rcounts[i]; } 


      MPI_Allgatherv(&colptr[0], colptr.size()-1, MPI_INT, &pGlobal.colptr[0],
          &rcounts[0], &rdispls[0], MPI_INT, comm);

      /* Recompute column pointers. */
      for (int p = 1; p < np; p++) {
        int idx = rdispls[p];
        for (int j = 0; j < rcounts[p]; ++j) pGlobal.colptr[idx++] += prevnz[p];
      }
      pGlobal.colptr[pGlobal.size]= pGlobal.nnz+1;



    }
  }

  void SparseMatrixStructure::GetLColRowCount(ETree & tree, Ordering & aOrder, vector<int> & cc, vector<int> & rc){
     //The tree need to be postordered
//    if(!tree.IsPostOrdered()){
//      tree.PostOrderTree(aOrder);
//    }

    ExpandSymmetric();
    
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
        int oldnbr = aOrder.perm[lownbr-1];
        int jstrt = expColptr[oldnbr-1];
        int jstop = expColptr[oldnbr] - 1;


        //           -----------------------------------------------
        //           for each ``high neighbor'', hinbr of lownbr ...
        //           -----------------------------------------------
        for(int j = jstrt; j<=jstop;++j){
          int hinbr = expRowind[j-1];
          hinbr = aOrder.invp[hinbr-1];
          if  ( hinbr > lownbr )  {
            if  ( ifdesc > prvnbr[hinbr-1] ) {
              //                       -------------------------
              //                       increment weight[lownbr].
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
                //                           decrement weight[lca].
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
        //           set ( p(lownbr) ) <-- set ( p(lownbr) ) + set[xsup].
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

  

  void SparseMatrixStructure::FindSupernodes(ETree& tree, Ordering & aOrder, vector<int> & cc,vector<int> & supMembership, vector<int> & xsuper, int maxSize ){

    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call FindSupernodes\n" );
    }

    supMembership.resize(size);

    int nsuper = 1;
    int supsize = 1;
    supMembership[0] = 1;

    for(int i =2; i<=size;i++){
      int prev_parent = tree.PostParent(i-2);
      if(prev_parent == i){
        if(cc[i-2] == cc[i-1]+1 ) {
          if(supsize<maxSize || maxSize==-1){
            ++supsize;
            supMembership[i-1] = nsuper;
            continue;
          }
        }
      }

        nsuper++;
      supsize = 1;
      supMembership[i-1] = nsuper;
    }

    xsuper.resize(nsuper+1);
    int lstsup = nsuper+1;
    for(int i = size; i>=1;--i){
      int ksup = supMembership[i-1];
      if(ksup!=lstsup){
       xsuper[lstsup-1] = i + 1; 
      }
      lstsup = ksup;
    }
    xsuper[0]=1;
  }



  void SparseMatrixStructure::FindFundamentalSupernodes(ETree& tree, Ordering & aOrder, vector<int> & cc,vector<int> & supMembership, vector<int> & xsuper, int maxSize ){

    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call FindSupernodes\n" );
    }

    supMembership.resize(size);

    int nsuper = 1;
    int supsize = 1;
    supMembership[0] = 1;




    vector<int> children(size,0);
    for(int I=1;I<=size;I++){
      int parent = tree.PostParent(I-1);
      if(parent!=0){
        ++children[parent-1];
      }
    }


    for(int i =2; i<=size;i++){
      int prev_parent = tree.PostParent(i-1-1);
      if(prev_parent == i){
        if(cc[i-1-1] == cc[i-1]+1 && children[i-1]==1) {
          if(supsize<maxSize || maxSize==-1){
            ++supsize;
            supMembership[i-1] = nsuper;
            continue;
          }
        }
      }

        nsuper++;
      supsize = 1;
      supMembership[i-1] = nsuper;
    }

    xsuper.resize(nsuper+1);
    int lstsup = nsuper+1;
    for(int i = size; i>=1;--i){
      int ksup = supMembership[i-1];
      if(ksup!=lstsup){
       xsuper[lstsup-1] = i + 1; 
      }
      lstsup = ksup;
    }
    xsuper[0]=1;
  }




//EXPERIMENTAL STUFF

  typedef std::set<int> nodeset;
  //typedef std::list<nodeset*> partitions;
  typedef std::vector<nodeset*> vecset;

  struct snode{
    int id;
    nodeset members;
  };

  typedef std::list<snode*> partitions;

void SparseMatrixStructure::RefineSupernodes(ETree& tree, Ordering & aOrder, vector<int> & supMembership, vector<int> & xsuper, vector<int64_t> & xlindx, vector<int32_t> & lindx, vector<int> & perm,vector<int> & origPerm){

#if 0
vector<int> newXsuper;
#endif

  perm.resize(size);
  origPerm.resize(size);

  int nsuper = xsuper.size()-1;
  //Using std datatypes first
  //adj[i-1] is the monotone adjacency set of supernode i 
  //(nodes numbered higher than the highest-numbered vertex in Supernode i)
  vecset adj(nsuper);

  //vecset snodes(nsuper);
  vector<snode*> snodes(nsuper);

  partitions L;
  //initialize L with curent supernodal partition
  int pos = 1;
  for(int i = 1; i<=nsuper; ++i){
    adj[i-1] = new nodeset();
    int64_t fi = xlindx[i-1];
    int64_t li = xlindx[i]-1;
    for(int64_t idx = fi; idx<=li;idx++){
       int32_t row = lindx[idx-1];
       if(row>xsuper[i]-1){
          adj[i-1]->insert(row);
       }
    }

    L.push_back(new snode());
    snode * curL = L.back();
    curL->id = i;
    snodes[i-1] = curL;

    int fc = xsuper[i-1];
    int lc = xsuper[i]-1;
    for(int node = fc; node<=lc;++node){
      curL->members.insert(node);
    }

//    L.push_back(new nodeset());
//    nodeset * curL = L.back();
//
//    snodes[i-1] = curL;
//
//    int fc = xsuper[i-1];
//    int lc = xsuper[i]-1;
//    for(int node = fc; node<=lc;++node){
//      curL->insert(node);
//    }

    for(int node = fc; node<=lc;++node){
      origPerm[node-1] = pos++;
    }
  }
 



    partitions::iterator it;

#if 1
    for(it = L.begin();it != L.end(); ++it){
      cout<<"[ ";
      for(nodeset::iterator nit = (*it)->members.begin();nit != (*it)->members.end(); ++nit){
        cout<<*nit<<" ";
      }
      cout<<"] ";
    }
    cout<<std::endl;
#endif







  partitions::reverse_iterator Kit;

  int prevK = nsuper+1;  
  for(Kit = L.rbegin();Kit != L.rend(); ++Kit){
    //assert( snodes[K-1] == *Kit);
    int K = (*Kit)->id;
    if(K<0 || K>=prevK){cout<<"K="<<K<<endl; continue;}

    cout<<"K_"<<K<<" is "<<(*Kit)->members<<std::endl;
    cout<<"Adj of K_"<<K<<" is "<<*adj[K-1]<<std::endl;

    partitions::reverse_iterator Jit;
//    partitions tmp;
    int newSets=0;
    auto KKit = Kit;
    for(Jit = L.rbegin(); Jit != KKit; ++Jit){
      cout<<"    L_"<<std::distance(L.begin(),Jit.base())<<" is "<<(*Jit)->members<<std::endl;
      //compute I = J inter M(Ki)
      nodeset * inter = new nodeset();
      std::set_intersection((*Jit)->members.begin(),(*Jit)->members.end(),
                              adj[K-1]->begin(),adj[K-1]->end(),
                                std::inserter(*inter, inter->begin()));

      if(inter->size()>0){
        //compute I' = J \ M(Ki)
        nodeset * diff = new nodeset();
        std::set_difference((*Jit)->members.begin(),(*Jit)->members.end(),
                              adj[K-1]->begin(),adj[K-1]->end(),
                                std::inserter(*diff, diff->begin()));

          //cout<<"        intersect is "<<*inter<<std::endl;
          //cout<<"        diff is "<<*diff<<std::endl;

        if(diff->size()>0){
          //increase Kit only if we are splitting an original set of L...
          int curJ =(*Jit)->id;
          int newJ = -abs(curJ);

          //replace J by I'
          (*Jit)->members.swap(*diff);
          (*Jit)->id=newJ;
          //Insert I before I'
          ++Jit; 
          L.insert(Jit.base(),new snode());
          (*Jit)->members.swap(*inter);
          (*Jit)->id=newJ;

          --Jit; 
          cout<<"B    L_"<<std::distance(L.begin(),Jit.base())<<" is now "<<(*Jit)->members<<std::endl;
          ++Jit; 
          cout<<"C    L_"<<std::distance(L.begin(),Jit.base())<<" is now "<<(*Jit)->members<<std::endl;

          KKit++;
        }
        delete diff;
        delete inter;
      }
      else{
        delete inter;
      }
    }
    prevK=K;
  }
/////
/////assert(K==0);
/////
///////    partitions::iterator it;
///////    for(it = L.begin();it != L.end(); ++it){
///////      //logfileptr->OFS()<<"[ ";
///////      for(nodeset::iterator nit = (*it)->begin();nit != (*it)->end(); ++nit){
///////      //  logfileptr->OFS()<<*nit<<" ";
///////      }
///////      //logfileptr->OFS()<<"] ";
///////    }
///////    //logfileptr->OFS()<<std::endl;
/////





#if 1
    for(it = L.begin();it != L.end(); ++it){
      cout<<"["<< (*it)->id <<": ";
      for(nodeset::iterator nit = (*it)->members.begin();nit != (*it)->members.end(); ++nit){
        cout<<*nit<<" ";
      }
      cout<<"] ";
    }
    cout<<std::endl;
#endif









  //construct perm
    pos = 1;
    for(it = L.begin();it != L.end(); ++it){
      for(nodeset::iterator nit = (*it)->members.begin();nit != (*it)->members.end(); ++nit){
        perm[*nit-1] = pos++;
      }
    }

    cout<<perm<<endl;


#if 0
    newXsuper.resize(L.size()+1);
    pos = 1;
    int I = 1;
    newXsuper[I-1] = 1;
    for(it = L.begin();it != L.end(); ++it){
      for(nodeset::iterator nit = (*it)->begin();nit != (*it)->end(); ++nit){
        pos++;
      }
      newXsuper[I]=pos;
      I++;
    }

    for(it = L.begin();it != L.end(); ++it){
      cerr<<"[ ";
      for(nodeset::iterator nit = (*it)->begin();nit != (*it)->end(); ++nit){
        cerr<<*nit<<" ";
      }
      cerr<<"] ";
    }
    cerr<<std::endl;



    for(int i =0;i<newXsuper.size();++i){
      cerr<<newXsuper[i]<<" ";
    }
    cerr<<endl;
#endif






  for(it = L.begin();it != L.end(); ++it){
    delete (*it);
  }

  for(int i = 1; i<=nsuper; ++i){
    delete adj[i-1];
  }

  aOrder.Compose(perm);

}

void SparseMatrixStructure::RelaxSupernodes(ETree& tree, vector<int> & cc,vector<int> & supMembership, vector<int> & xsuper, int maxSize ){

    int nsuper = xsuper.size()-1;

    DisjointSet sets;
    sets.Initialize(nsuper);
    vector<int> ncols(nsuper);
    vector<int> zeros(nsuper);
    vector<int> newCC(nsuper);
    for(int ksup=nsuper;ksup>=1;--ksup){
      int cset = sets.makeSet(ksup);
      sets.Root(cset-1)=ksup;
      
      int fstcol = xsuper[ksup-1];
      int lstcol = xsuper[ksup]-1;
      int width = lstcol - fstcol +1;
      int length = cc[fstcol-1];
      ncols[ksup-1] = width;
      zeros[ksup-1] = 0;
      newCC[ksup-1] = length;
    }


  //minsize
#if 0
  int nrelax0 = min(4,maxSize);
  int nrelax1 = min(16,maxSize);
  int nrelax2 = min(48,maxSize);
#else
  int nrelax0 = min(8,maxSize);
  int nrelax1 = min(32,maxSize);
  int nrelax2 = min(64,maxSize);




  double zrelax0 = 0.8;
  double zrelax1 = 0.1;
  double zrelax2 = 0.05;

  for(int ksup=nsuper;ksup>=1;--ksup){
      int fstcol = xsuper[ksup-1];
      int lstcol = xsuper[ksup]-1;
      int width = ncols[ksup-1];
      int length = cc[fstcol-1];

      int parent_fstcol = tree.PostParent(lstcol-1);
      if(parent_fstcol!=0){
        int parent_snode = supMembership[parent_fstcol-1];
        int pset = sets.find(parent_snode);
        parent_snode = sets.Root(pset-1);

        bool merge = (parent_snode == ksup+1);
 

        if(merge){
          int parent_width = ncols[parent_snode-1];

          int parent_fstcol = xsuper[parent_snode-1];
          int parent_lstcol = xsuper[parent_snode]-1;
          int totzeros = zeros[parent_snode-1];
          int fused_cols = width + parent_width;
          
          merge = false;
          if(fused_cols <= nrelax0){
            merge = true;
          }
          else if(fused_cols <=maxSize){
            double child_lnz = cc[fstcol-1];
            double parent_lnz = cc[parent_fstcol-1];
            double xnewzeros = width * (parent_lnz + width  - child_lnz);

            if(xnewzeros == 0){
              merge = true;
            }
            else{
              //all these values are the values corresponding to the merged snode
              double xtotzeros = (double)totzeros + xnewzeros;
              double xfused_cols = (double) fused_cols;
              //new number of nz
              double xtotsize = (xfused_cols * (xfused_cols+1)/2) + xfused_cols * (parent_lnz - parent_width);
              //percentage of explicit zeros
              double z = xtotzeros / xtotsize;

              int totsize = (fused_cols * (fused_cols+1)/2) + fused_cols * ((int)parent_lnz - parent_width);
              totzeros += (int)xnewzeros;

              merge = ((fused_cols <= nrelax1 && z < zrelax0) 
                          || (fused_cols <= nrelax2 && z < zrelax1)
                              || (z<zrelax2)) &&
                            (xtotsize < std::numeric_limits<int>::max() / sizeof(double));
            }

          }

          if(merge){
//              std::cout<<"merge "<<ksup<<" and "<<parent_snode<<std::endl;
            ncols[ksup-1] += ncols[parent_snode-1]; 
            zeros[ksup-1] = totzeros;
            newCC[ksup-1] = width + newCC[parent_snode-1];
            sets.Union(ksup,parent_snode,ksup);
          }
        } 

      }
  }

    vector<int> relXSuper(nsuper+1);
    int nrSuper = 0;
    for(int ksup=1;ksup<=nsuper;++ksup){
        int kset = sets.find(ksup);
        if(ksup == sets.Root(kset-1)){
          int fstcol = xsuper[ksup-1];
          relXSuper[nrSuper] = fstcol;
          newCC[nrSuper] = newCC[ksup-1];
          ++nrSuper;
        }
    }
    relXSuper[nrSuper] = xsuper[nsuper];
    relXSuper.resize(nrSuper+1);

    for(int ksup=1;ksup<=nrSuper;++ksup){
      int fstcol = relXSuper[ksup-1];
      int lstcol = relXSuper[ksup]-1;
      for(int col = fstcol; col<=lstcol;++col){
        supMembership[col-1] = ksup;
        cc[col-1] = newCC[ksup-1] + col-fstcol;
      }
    }
    
    xsuper = relXSuper;
///      //adjust the column counts
///      for(int col=i-2;col>=i-supsize;--col){
///        cc[col-1] = cc[col]+1;
///      }


}

  void SparseMatrixStructure::SymbolicFactorizationRelaxed(ETree& tree,Ordering & aOrder, const vector<int> & cc,const vector<int> & xsuper,const vector<int> & SupMembership, vector<int64_t> & xlindx, vector<int32_t> & lindx){


    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call SymbolicFactorization\n" );
    }

    int nsuper = xsuper.size()-1;




    int64_t nzbeg = 0;
    //nzend points to the last used slot in lindx
    int64_t nzend = 0;

    //tail is the end of list indicator (in rchlnk, not mrglnk)
    int32_t tail = size +1;

    int32_t head = 0;

    //Array of length nsuper containing the children of 
    //each supernode as a linked list
    vector<int32_t> mrglnk;
    mrglnk.assign(nsuper,0);

    //Array of length n+1 containing the current linked list 
    //of merged indices (the "reach" set)
    vector<int32_t> rchlnk(size+1);

    //Array of length n used to mark indices as they are introduced
    // into each supernode's index set
    vector<int> marker;
    marker.assign(size,0);


    xlindx.resize(nsuper+1);

    //Compute the sum of the column count and resize lindx accordingly
    int64_t nofsub = 1;
    for(int i =0; i<cc.size();++i){
      nofsub+=cc[i];
    }

    lindx.resize(nofsub);


    int64_t point = 1;
    for(int ksup = 1; ksup<=nsuper; ++ksup){
      int fstcol = xsuper[ksup-1];
      xlindx[ksup-1] = point;
      point += cc[fstcol-1]; 
    } 
    xlindx[nsuper] = point;



    for(int ksup = 1; ksup<=nsuper; ++ksup){
      int fstcol = xsuper[ksup-1];
      int lstcol = xsuper[ksup]-1;
      int width = lstcol - fstcol +1;
      int length = cc[fstcol-1];
      int64_t knz = 0;
      rchlnk[head] = tail;
      int jsup = mrglnk[ksup-1];

      //If ksup has children in the supernodal e-tree
      if(jsup>0){
        //copy the indices of the first child jsup into 
        //the linked list, and mark each with the value 
        //ksup.
        int jwidth = xsuper[jsup]-xsuper[jsup-1];
        int64_t jnzbeg = xlindx[jsup-1] + jwidth;
        int64_t jnzend = xlindx[jsup] -1;
        for(int64_t jptr = jnzend; jptr>=jnzbeg; --jptr){
          int32_t newi = lindx[jptr-1];
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

          jwidth = xsuper[jsup]-xsuper[jsup-1];
          jnzbeg = xlindx[jsup-1] + jwidth;
          jnzend = xlindx[jsup] -1;
          int nexti = head;
          for(int64_t jptr = jnzbeg; jptr<=jnzend; ++jptr){
            int32_t newi = lindx[jptr-1];
            int32_t i;
            do{
              i = nexti;
              nexti = rchlnk[i];
            }while(newi > nexti);

            if(newi < nexti){
#ifdef _DEBUG_
            logfileptr->OFS()<<jsup<<" is a child of "<<ksup<<" and "<<newi<<" is inserted in the structure of "<<ksup<<std::endl;
#endif
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
        for(int row = fstcol; row<=lstcol; ++row){
          int32_t newi = row;
          if(newi > fstcol && marker[newi-1] != ksup){
            //position and insert newi in list and
            // mark it with kcol
            int32_t nexti = head;
            int32_t i;
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


        for(int col = fstcol; col<=lstcol; ++col){
          int node = aOrder.perm[col-1];

          int64_t knzbeg = expColptr[node-1];
          int64_t knzend = expColptr[node]-1;
          for(int64_t kptr = knzbeg; kptr<=knzend;++kptr){
            int32_t newi = expRowind[kptr-1];
            newi = aOrder.invp[newi-1];
            
            if(newi > fstcol && marker[newi-1] != ksup){
              //position and insert newi in list and
              // mark it with kcol
              int32_t nexti = head;
              int32_t i;
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

      } 

      //if ksup has no children, insert fstcol into the linked list.
      if(rchlnk[head] != fstcol){
        rchlnk[fstcol] = rchlnk[head];
        rchlnk[head] = fstcol;
        ++knz;
      }

      assert(knz == cc[fstcol-1]);


      //copy indices from linked list into lindx(*).
      nzbeg = nzend+1;
      nzend += knz;
      assert(nzend+1 == xlindx[ksup]);
      int32_t i = head;
      for(int64_t kptr = nzbeg; kptr<=nzend;++kptr){
        i = rchlnk[i];
        lindx[kptr-1] = i;
      } 

      //if ksup has a parent, insert ksup into its parent's 
      //"merge" list.
      if(length > width){
        int32_t pcol = lindx[xlindx[ksup-1] + width -1];
        int psup = SupMembership[pcol-1];
        mrglnk[ksup-1] = mrglnk[psup-1];
        mrglnk[psup-1] = ksup;
      }
    }

    lindx.resize(nzend+1);

  }
#endif

  void SparseMatrixStructure::SymbolicFactorization(ETree& tree, Ordering & aOrder, const vector<int> & cc,const vector<int> & xsuper,const vector<int> & SupMembership, vector<int64_t> & xlindx, vector<int32_t> & lindx){


    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call SymbolicFactorization\n" );
    }

    int nsuper = xsuper.size()-1;

    int64_t nzbeg = 0;
    //nzend points to the last used slot in lindx
    int64_t nzend = 0;

    //tail is the end of list indicator (in rchlnk, not mrglnk)
    int32_t tail = size +1;

    int32_t head = 0;

    //Array of length nsuper containing the children of 
    //each supernode as a linked list
    vector<int32_t> mrglnk;
    mrglnk.assign(nsuper,0);

    //Array of length n+1 containing the current linked list 
    //of merged indices (the "reach" set)
    vector<int32_t> rchlnk;
    rchlnk.assign(size+1,0);

    //Array of length n used to mark indices as they are introduced
    // into each supernode's index set
    vector<int> marker(size,0);



    xlindx.resize(nsuper+1);

    //Compute the sum of the column count and resize lindx accordingly
    int64_t nofsub = 1;
    for(int i =0; i<cc.size();++i){
      nofsub+=cc[i];
    }

    lindx.resize(nofsub);


    int64_t point = 1;
    for(int ksup = 1; ksup<=nsuper; ++ksup){
      int fstcol = xsuper[ksup-1];
      xlindx[ksup-1] = point;
      point += cc[fstcol-1]; 
    } 
    xlindx[nsuper] = point;

    for(int ksup = 1; ksup<=nsuper; ++ksup){
      int fstcol = xsuper[ksup-1];
      int lstcol = xsuper[ksup]-1;
      int width = lstcol - fstcol +1;
      int length = cc[fstcol-1];
      int32_t knz = 0;
      rchlnk[head] = tail;
      int jsup = mrglnk[ksup-1];

      //If ksup has children in the supernodal e-tree
      if(jsup>0){
        //copy the indices of the first child jsup into 
        //the linked list, and mark each with the value 
        //ksup.
        int jwidth = xsuper[jsup]-xsuper[jsup-1];
        int64_t jnzbeg = xlindx[jsup-1] + jwidth;
        int64_t jnzend = xlindx[jsup] -1;
        for(int64_t jptr = jnzend; jptr>=jnzbeg; --jptr){
          int32_t newi = lindx[jptr-1];
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

          jwidth = xsuper[jsup]-xsuper[jsup-1];
          jnzbeg = xlindx[jsup-1] + jwidth;
          jnzend = xlindx[jsup] -1;
          int32_t nexti = head;
          for(int64_t jptr = jnzbeg; jptr<=jnzend; ++jptr){
            int32_t newi = lindx[jptr-1];
            int32_t i;
            do{
              i = nexti;
              nexti = rchlnk[i];
            }while(newi > nexti);

            if(newi < nexti){
#ifdef _DEBUG_
            logfileptr->OFS()<<jsup<<" is a child of "<<ksup<<" and "<<newi<<" is inserted in the structure of "<<ksup<<std::endl;
#endif
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
        int node = aOrder.perm[fstcol-1];
        int64_t knzbeg = expColptr[node-1];
        int64_t knzend = expColptr[node]-1;
        for(int64_t kptr = knzbeg; kptr<=knzend;++kptr){
          int32_t newi = expRowind[kptr-1];
          newi = aOrder.invp[newi-1];
          if(newi > fstcol && marker[newi-1] != ksup){
            //position and insert newi in list and
            // mark it with kcol
            int32_t nexti = head;
            int32_t i;
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

      assert(knz == cc[fstcol-1]);


      //copy indices from linked list into lindx(*).
      nzbeg = nzend+1;
      nzend += knz;
      assert(nzend+1 == xlindx[ksup]);
      int i = head;
      for(int kptr = nzbeg; kptr<=nzend;++kptr){
        i = rchlnk[i];
        lindx[kptr-1] = i;
      } 

      //if ksup has a parent, insert ksup into its parent's 
      //"merge" list.
      if(length > width){
        int pcol = lindx[xlindx[ksup-1] + width -1];
        int psup = SupMembership[pcol-1];
        mrglnk[ksup-1] = mrglnk[psup-1];
        mrglnk[psup-1] = ksup;
      }
    }

    lindx.resize(nzend+1);

  }



}
