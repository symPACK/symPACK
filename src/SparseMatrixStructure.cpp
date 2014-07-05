#include "SparseMatrixStructure.hpp"
#include "ETree.hpp"
#include "utility.hpp"
#include <limits>       // std::numeric_limits

#include <iterator>
#include <set>
#include <list>
#include <vector>
#include <algorithm>

namespace LIBCHOLESKY{
  void SparseMatrixStructure::ClearExpandedSymmetric(){
    expColptr.Resize(0);
    expRowind.Resize(0);
    bIsExpanded=false;
  }

  void SparseMatrixStructure::ExpandSymmetric(){
    if(!bIsExpanded){
      //code from sparsematrixconverter
      /* set-up */

      IntNumVec cur_col_nnz(size);
      IntNumVec new_col_nnz(size);

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
      Int new_nnz = 0; 
      for (Int i = 0; i < size; i++) 
      {    
        cur_col_nnz[i] = colptr[i+1] - colptr[i];
        new_col_nnz[i] = cur_col_nnz[i];
        new_nnz += new_col_nnz[i];
      }    

      for (Int i = 0; i < size; i++) 
      {    
        Int k;
        for (k = colptr[i]; k < colptr[i+1]; k++) 
        {
          Int j = rowind[k-1]-1;
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
      expColptr.Resize(size+1);
      expColptr[0] = 1;
      for (Int i = 1; i <= size; i++)
      {
        expColptr[i] = expColptr[i-1] + new_col_nnz[i-1];
        new_col_nnz[i-1] = cur_col_nnz[i-1];
      }
      expColptr[size] = new_nnz+1;

      expRowind.Resize(new_nnz);

      /*
       *  Complete expansion of A to full storage.
       *
       *  Post:
       *    (new_colptr, new_rowind, new_values) is the full-storage equivalent of A.
       *    new_col_nnz[i] == # of non-zeros in col i of A.
       */




      for (Int i = 0; i < size; i++)
      {
        Int cur_nnz = cur_col_nnz[i];
        Int k_cur   = colptr[i] -1;
        Int k_new   = expColptr[i] -1;


        /* copy current non-zeros from old matrix to new matrix */
        std::copy(rowind.Data() + k_cur, rowind.Data() + k_cur + cur_nnz , expRowind.Data() + k_new);

        /* fill in the symmetric "missing" values */
        while (k_cur < colptr[i+1]-1)
        {
          /* non-zero of original matrix */
          Int j = rowind[k_cur]-1;

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
      }
      bIsExpanded =true;
    }

  }


  void SparseMatrixStructure::ToGlobal(SparseMatrixStructure & pGlobal){

    TIMER_START(ToGlobalStructure);
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
      MPI_Comm comm = MPI_COMM_WORLD;

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


      pGlobal.bIsExpanded = false;
      pGlobal.size = size;
      pGlobal.colptr.Resize(size+1);


      /* Allgatherv for row indices. */ 
      IntNumVec prevnz(np);
      IntNumVec rcounts(np);
      MPI_Allgather(&nnz, 1, MPI_INT, rcounts.Data(), 1, MPI_INT, comm);

      prevnz[0] = 0;
      for (Int i = 0; i < np-1; ++i) { prevnz[i+1] = prevnz[i] + rcounts[i]; } 

      pGlobal.nnz = 0;
      for (Int i = 0; i < np; ++i) { pGlobal.nnz += rcounts[i]; } 
      pGlobal.rowind.Resize(pGlobal.nnz);


      //    logfileptr->OFS()<<"Global nnz is "<<pGlobal.nnz<<std::endl;

      MPI_Allgatherv(rowind.Data(), nnz, MPI_INT, pGlobal.rowind.Data(),rcounts.Data(), prevnz.Data(), MPI_INT, comm); 

      //    logfileptr->OFS()<<"Global rowind is "<<pGlobal.rowind<<std::endl;

      /* Allgatherv for colptr */
      // Compute the number of columns on each processor
      Int numColFirst = size / np;
      SetValue( rcounts, numColFirst );
      rcounts[np-1] = size - numColFirst * (np-1);  // Modify the last entry     


      IntNumVec rdispls(np);
      rdispls[0] = 0;
      for (Int i = 0; i < np-1; ++i) { rdispls[i+1] = rdispls[i] + rcounts[i]; } 


      MPI_Allgatherv(colptr.Data(), colptr.m()-1, MPI_INT, pGlobal.colptr.Data(),
          rcounts.Data(), rdispls.Data(), MPI_INT, comm);

      /* Recompute column pointers. */
      for (Int p = 1; p < np; p++) {
        Int idx = rdispls[p];
        for (Int j = 0; j < rcounts[p]; ++j) pGlobal.colptr[idx++] += prevnz[p];
      }
      pGlobal.colptr(pGlobal.size)= pGlobal.nnz+1;

      //    logfileptr->OFS()<<"Global colptr is "<<pGlobal.colptr<<std::endl;


      pGlobal.bIsGlobal = true;
    }

    TIMER_STOP(ToGlobalStructure);
  }

  void SparseMatrixStructure::GetLColRowCount2(ETree & tree, IntNumVec & cc, IntNumVec & rc){
     //The tree need to be postordered
    if(!tree.IsPostOrdered()){
      tree.PostOrderTree();
    }

    ExpandSymmetric();
    
    TIMER_START(GetColRowCount_Classic);
    cc.Resize(size);
    rc.Resize(size);

    IntNumVec level(size+1);
    IntNumVec weight(size+1);
    IntNumVec fdesc(size+1);
    IntNumVec nchild(size+1);
    IntNumVec set(size);
    IntNumVec prvlf(size);
    IntNumVec prvnbr(size);


        Int xsup = 1;
        level(0) = 0;
      for(Int k = size; k>=1; --k){
            rc(k-1) = 1;
            cc(k-1) = 0;
            set(k-1) = k;
            prvlf(k-1) = 0;
            prvnbr(k-1) = 0;
            level(k) = level(tree.PostParent(k-1)) + 1;
            weight(k) = 1;
            fdesc(k) = k;
            nchild(k) = 0;
      }

      nchild(0) = 0;
      fdesc(0) = 0;
      for(Int k =1; k<size; ++k){
            Int parent = tree.PostParent(k-1);
            weight(parent) = 0;
            ++nchild(parent);
            Int ifdesc = fdesc(k);
            if  ( ifdesc < fdesc(parent) ) {
                fdesc(parent) = ifdesc;
            }
      }







      for(Int lownbr = 1; lownbr<=size; ++lownbr){
        Int lflag = 0;
        Int ifdesc = fdesc(lownbr);
        Int oldnbr = tree.FromPostOrder(lownbr);
        Int jstrt = expColptr(oldnbr-1);
        Int jstop = expColptr(oldnbr) - 1;

//        if(nchild(lownbr)>=2){
//          weight(lownbr)--;
//        }

        //           -----------------------------------------------
        //           for each ``high neighbor'', hinbr of lownbr ...
        //           -----------------------------------------------
        for(Int j = jstrt; j<=jstop;++j){
          Int hinbr = tree.ToPostOrder(expRowind(j-1));
          if  ( hinbr > lownbr )  {
            if  ( ifdesc > prvnbr(hinbr-1) ) {
              //                       -------------------------
              //                       increment weight(lownbr).
              //                       -------------------------
              ++weight(lownbr);
              Int pleaf = prvlf(hinbr-1);
              //                       -----------------------------------------
              //                       if hinbr has no previous ``low neighbor'' 
              //                       then ...
              //                       -----------------------------------------
              if  ( pleaf == 0 ) {
                //                           -----------------------------------------
                //                           ... accumulate lownbr-->hinbr path length 
                //                               in rowcnt(hinbr).
                //                           -----------------------------------------
                rc(hinbr-1) += level(lownbr) - level(hinbr);
              }
              else{
                //                           -----------------------------------------
                //                           ... otherwise, lca <-- find(pleaf), which 
                //                               is the least common ancestor of pleaf 
                //                               and lownbr.
                //                               (path halving.)
                //                           -----------------------------------------
                Int last1 = pleaf;
                Int last2 = set(last1-1);
                Int lca = set(last2-1);
                while(lca != last2){
                  set(last1-1) = lca;
                  last1 = lca;
                  last2 = set(last1-1);
                  lca = set(last2-1);
                }
                //                           -------------------------------------
                //                           accumulate pleaf-->lca path length in 
                //                           rowcnt(hinbr).
                //                           decrement weight(lca).
                //                           -------------------------------------
                rc(hinbr-1) += level(lownbr) - level(lca);
                --weight(lca);
              }
              //                       ----------------------------------------------
              //                       lownbr now becomes ``previous leaf'' of hinbr.
              //                       ----------------------------------------------
              prvlf(hinbr-1) = lownbr;
              lflag = 1;
            }
            //                   --------------------------------------------------
            //                   lownbr now becomes ``previous neighbor'' of hinbr.
            //                   --------------------------------------------------
            prvnbr(hinbr-1) = lownbr;
          }
        }
        //           ----------------------------------------------------
        //           decrement weight ( parent(lownbr) ).
        //           set ( p(lownbr) ) <-- set ( p(lownbr) ) + set(xsup).
        //           ----------------------------------------------------
        Int parent = tree.PostParent(lownbr-1);
        --weight(parent);


        //merge the sets
        if  ( lflag == 1  || nchild(lownbr) >= 2 ) {
          xsup = lownbr;
        }
        set(xsup-1) = parent;
      }



#ifdef _DEBUG_
      logfileptr->OFS()<<"deltas "<<weight<<std::endl;
#endif

        for(Int k = 1; k<=size; ++k){
            Int temp = cc(k-1) + weight(k);
            cc(k-1) = temp;
            Int parent = tree.PostParent(k-1);
            if  ( parent != 0 ) {
                cc(parent-1) += temp;
            }
        }

//      logfileptr->OFS()<<"column counts "<<cc<<std::endl;

      TIMER_STOP(GetColRowCount_Classic);
  }

  void SparseMatrixStructure::GetLColRowCount(ETree & tree, IntNumVec & cc, IntNumVec & rc){


    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call GetLColRowCount\n" );
    }

    //The tree need to be postordered
    if(!tree.IsPostOrdered()){
      tree.PostOrderTree();
    }

    ExpandSymmetric();

    TIMER_START(GetColRowCount);
//    TIMER_START(Initialize_Data);
    //cc first contains the delta
    cc.Resize(size);
    //Compute size of subtrees
    IntNumVec treeSize(size);
    SetValue(treeSize,I_ONE);



    IntNumVec level(size);
    level(size-1)=1;
    for(Int vertex = 1; vertex<=size-1; vertex++){
      Int curParent = tree.PostParent(vertex-1);
      if(curParent!=0){
        treeSize(curParent-1)+=treeSize(vertex-1);
      }
      if(treeSize(vertex-1)==1){
        cc(vertex-1)=1;
      }
      else{
        cc(vertex-1)= 0;
      }
    }

    for(Int vertex = size-1; vertex>=1; --vertex){
      Int curParent = tree.PostParent(vertex-1);
      if(curParent!=0){
        level(vertex-1) = level(curParent-1)+1;
      }
      else{
        level(vertex-1) = 0;
      }
    }


    if(treeSize(size-1)==1){
      cc(size-1)=1;
    }
    else{
      cc(size-1)= 0 ;
    }



    IntNumVec prevLeaf(size);
    SetValue(prevLeaf,I_ZERO);
    IntNumVec prevNz(size);
    SetValue(prevNz,I_ZERO);

    rc.Resize(size);
    SetValue(rc,I_ONE);

    DisjointSet sets;
    sets.Initialize(size);
    for(Int vertex = 1; vertex<=size; vertex++){
      Int cset = sets.makeSet (vertex);
      sets.Root(cset-1)=vertex;
    }


//    TIMER_STOP(Initialize_Data);

//    TIMER_START(Compute_Col_Row_Count);
    for(Int col=1; col<size; col++){
      Int cset;

      Int colPar = tree.PostParent(col-1);
      if (col<size && colPar!=0){
        cc(colPar-1)--;
      }

      Int oCol = tree.FromPostOrder(col);
      for (Int i = expColptr(oCol-1); i < expColptr(oCol); i++) {
        Int row = tree.ToPostOrder(expRowind(i-1));
        if (row > col){
          Int k = prevNz(row-1);


#ifdef _DEBUG_
          logfileptr->OFS()<<"prevNz("<<row<<")="<<k<<" vs "<< col - treeSize(col-1) +1<<std::endl;
#endif
          if(k< col - treeSize(col-1) +1){
#ifdef _DEBUG_
            logfileptr->OFS()<<"Vertex "<<col<<" is a leaf of Tr["<<row<<"]"<<std::endl;
#endif
            cc(col-1)++;

            Int p = prevLeaf(row-1);
            if(p==0){
              rc(row-1)+=level(col-1)-level(row-1);
            }
            else {
//              TIMER_START(Get_LCA);
              Int pset = sets.find(p);
              Int q = sets.Root(pset-1);
//              TIMER_STOP(Get_LCA);

#ifdef _DEBUG_
              logfileptr->OFS()<<"Vertex "<<q<<" is the LCA of "<<p<<" and "<< col<<std::endl;
#endif
              rc(row-1)+= level(col-1) - level(q-1);
              cc(q-1)--;

            }
            prevLeaf(row-1)=col;
          }
#ifdef _DEBUG_
          else{
            logfileptr->OFS()<<"Vertex "<<col<<" is an internal vertex of Tr["<<row<<"]"<<std::endl;
          }
#endif
          prevNz(row-1)=col;
        }
      }

//      TIMER_START(Merge_LCA);
      //merge col and parent sets (for lca computation)
      if (colPar!=0){
        sets.Union(col,colPar);
      }
//      TIMER_STOP(Merge_LCA);

    }
//    TIMER_STOP(Compute_Col_Row_Count);


#ifdef _DEBUG_
        logfileptr->OFS()<<"Deltas "<<cc.m()<<std::endl;
        for(Int i = 0; i<cc.m();i++){
          logfileptr->OFS()<<cc(i)<<" ";
        }
        logfileptr->OFS()<<std::endl;
#endif





    //convert delta to col count
    for(Int col=1; col<size; col++){
      Int parent = tree.PostParent(col-1);
      if(parent!=0){
        cc(parent-1)+= cc(col-1);
      }
    }




    TIMER_STOP(GetColRowCount);
  }


  

  void SparseMatrixStructure::FindSupernodes(ETree& tree, IntNumVec & cc,IntNumVec & supMembership, IntNumVec & xsuper, Int maxSize ){
    TIMER_START(FindSupernodes);

    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call FindSupernodes\n" );
    }

    supMembership.Resize(size);

    Int nsuper = 1;
    Int supsize = 1;
    supMembership(0) = 1;

    for(Int i =2; i<=size;i++){
      Int prev_parent = tree.PostParent(i-2);
      if(prev_parent == i){
        if(cc(i-2) == cc(i-1)+1 ) {
          if(supsize<=maxSize || maxSize==-1){
            ++supsize;
            supMembership(i-1) = nsuper;
            continue;
          }
        }
      }

        nsuper++;
      supsize = 1;
      supMembership(i-1) = nsuper;
    }

    xsuper.Resize(nsuper+1);
    Int lstsup = nsuper+1;
    for(Int i = size; i>=1;--i){
      Int ksup = supMembership(i-1);
      if(ksup!=lstsup){
       xsuper(lstsup-1) = i + 1; 
      }
      lstsup = ksup;
    }
    xsuper(0)=1;
    TIMER_STOP(FindSupernodes);
  }



  typedef std::set<Int> nodeset;
  typedef std::list<nodeset*> partitions;
  typedef std::vector<nodeset*> vecset;
void SparseMatrixStructure::RefineSupernodes(ETree& tree, IntNumVec & supMembership, IntNumVec & xsuper, IntNumVec & xlindx, IntNumVec & lindx, IntNumVec & perm){

  perm.Resize(size);
//  for(Int i =1; i<=size; ++i){
//    perm[i-1]=i;
//  }

  Int nsuper = xsuper.m()-1;
  //Using std datatypes first

  vecset adj(nsuper);
  vecset snodes(nsuper);

  partitions L;

  IntNumVec origPerm(size);

  //init L with curent supernodal partition
  Int pos = 1;
  for(Int i = 1; i<=nsuper; ++i){
    adj[i-1] = new nodeset();
    Int fi = xlindx(i-1);
    Int li = xlindx(i)-1;
    for(Int idx = fi; idx<=li;idx++){
       adj[i-1]->insert(lindx(idx-1));
    }

    L.push_back(new nodeset());
    snodes[i-1] = L.back();
    nodeset * curL = L.back();
    Int fc = xsuper(i-1);
    Int lc = xsuper(i)-1;
    for(Int node = fc; node<=lc;++node){
      curL->insert(node);
      origPerm[node-1] = pos++;
    }
  }
 
  Int K = nsuper;
  partitions::reverse_iterator Kit;
  for(Kit = L.rbegin();Kit != L.rend(); ++Kit){
    logfileptr->OFS()<<"Looking at snode "<<K<<std::endl;

    assert( snodes[K-1] == *Kit);

//        logfileptr->OFS()<<"Adj is "<<*adj[K-1]<<std::endl;

    partitions::reverse_iterator Jit;
//    partitions tmp;
    Jit = L.rbegin();
    Int count = L.size() - K;
    while(count>0){
//        logfileptr->OFS()<<"L is "<<**Jit<<std::endl;
    //for(Jit = L.rbegin();Jit!=Kit;++Jit){
      nodeset * inter = new nodeset();

      std::set_intersection((*Jit)->begin(),(*Jit)->end(),
                              adj[K-1]->begin(),adj[K-1]->end(),
                                std::inserter(*inter, inter->begin()));
//        logfileptr->OFS()<<"Intersect is "<<*inter<<std::endl;

      if(inter->size()>0){
        nodeset * diff = new nodeset();
        std::set_difference((*Jit)->begin(),(*Jit)->end(),
                              adj[K-1]->begin(),adj[K-1]->end(),
                                std::inserter(*diff, diff->begin()));

//        logfileptr->OFS()<<"Diff is "<<*diff<<std::endl;

        if(diff->size()>0){
//          tmp.push_back(diff);
//          tmp.push_back(inter);
          //replace Jit by inter and diff
          (*Jit)->swap(*diff);
          L.insert(Jit.base(),inter);
          delete diff;
//          (*Jit)->swap(*inter);
//          L.insert(Jit.base(),diff);
//          delete inter;
          ++Jit; 
          ++Jit; 
        }
        else{
//          tmp.push_back(*Jit);
          delete diff;
        }
      }
      else{
//        tmp.push_back(*Jit);
        delete inter;
      }
    //}
      ++Jit;
      --count;
    }
    

    partitions::iterator it;
    for(it = L.begin();it != L.end(); ++it){
      logfileptr->OFS()<<"[ ";
      for(nodeset::iterator nit = (*it)->begin();nit != (*it)->end(); ++nit){
        logfileptr->OFS()<<*nit<<" ";
      }
      logfileptr->OFS()<<"] ";
    }
    logfileptr->OFS()<<std::endl;
    --K;
  }

    partitions::iterator it;
    for(it = L.begin();it != L.end(); ++it){
      logfileptr->OFS()<<"[ ";
      for(nodeset::iterator nit = (*it)->begin();nit != (*it)->end(); ++nit){
        logfileptr->OFS()<<*nit<<" ";
      }
      logfileptr->OFS()<<"] ";
    }
    logfileptr->OFS()<<std::endl;


  //construct perm
    pos = 1;
    for(it = L.begin();it != L.end(); ++it){
      for(nodeset::iterator nit = (*it)->begin();nit != (*it)->end(); ++nit){
        perm[*nit-1] = pos++;
      }
    }

    logfileptr->OFS()<<"Orig col order "<<origPerm<<std::endl;
    logfileptr->OFS()<<"Refined col order "<<perm<<std::endl;

//  IntNumVec lindxTemp = lindx;
//  //change lindx to reflect the new ordering
//  for(Int i = 1; i<=nsuper; ++i){
//    Int fi = xlindx(i-1);
//    Int li = xlindx(i)-1;
//    for(Int idx = fi; idx<=li;idx++){
//       lindx[idx-1] = perm[lindxTemp[idx-1]-1];
//    }
//  }
//
//    logfileptr->OFS()<<"Previous lindx "<<lindxTemp<<std::endl;
//    logfileptr->OFS()<<"Refined lindx "<<lindx<<std::endl;


  for(it = L.begin();it != L.end(); ++it){
    delete (*it);
  }

  for(Int i = 1; i<=nsuper; ++i){
    delete adj[i-1];
  }
}


















#ifdef RELAXED_SNODE
void SparseMatrixStructure::RelaxSupernodes(ETree& tree, IntNumVec & cc,IntNumVec & supMembership, IntNumVec & xsuper, Int maxSize ){

    Int nsuper = xsuper.m()-1;

    DisjointSet sets;
    sets.Initialize(nsuper);
    IntNumVec ncols(nsuper);
    IntNumVec zeros(nsuper);
    IntNumVec newCC(nsuper);
    for(Int ksup=nsuper;ksup>=1;--ksup){
      Int cset = sets.makeSet(ksup);
      sets.Root(cset-1)=ksup;
      
      Int fstcol = xsuper(ksup-1);
      Int lstcol = xsuper(ksup)-1;
      Int width = lstcol - fstcol +1;
      Int length = cc(fstcol-1);
      ncols[ksup-1] = width;
      zeros[ksup-1] = 0;
      newCC[ksup-1] = length;
    }


  //minsize
  Int nrelax0 = 4;
  Int nrelax1 = 16;
  Int nrelax2 = 48;

  double zrelax0 = 0.8;
  double zrelax1 = 0.1;
  double zrelax2 = 0.05;

  for(Int ksup=nsuper;ksup>=1;--ksup){
      Int fstcol = xsuper(ksup-1);
      Int lstcol = xsuper(ksup)-1;
      Int width = ncols[ksup-1];
      Int length = cc(fstcol-1);

      Int parent_fstcol = tree.PostParent(lstcol-1);
      if(parent_fstcol!=0){
        Int parent_snode = supMembership[parent_fstcol-1];
        Int pset = sets.find(parent_snode);
        parent_snode = sets.Root(pset-1);

        bool merge = (parent_snode == ksup+1);
 

        if(merge){
          Int parent_width = ncols[parent_snode-1];

          Int parent_fstcol = xsuper(parent_snode-1);
          Int parent_lstcol = xsuper(parent_snode)-1;
          Int totzeros = zeros[parent_snode-1];
          Int fused_cols = width + parent_width;
          
          merge = false;
          if(fused_cols <= nrelax0){
            merge = true;
          }
          else if(fused_cols <=maxSize){
            double child_lnz = cc(fstcol-1);
            double parent_lnz = cc(parent_fstcol-1);
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

              Int totsize = (fused_cols * (fused_cols+1)/2) + fused_cols * ((Int)parent_lnz - parent_width);
              totzeros += (Int)xnewzeros;

              merge = ((fused_cols <= nrelax1 && z < zrelax0) 
                          || (fused_cols <= nrelax2 && z < zrelax1)
                              || (z<zrelax2)) &&
                            (xtotsize < std::numeric_limits<Int>::max() / sizeof(double));
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

    IntNumVec relXSuper(nsuper+1);
    Int nrSuper = 0;
    for(Int ksup=1;ksup<=nsuper;++ksup){
        Int kset = sets.find(ksup);
        if(ksup == sets.Root(kset-1)){
          Int fstcol = xsuper[ksup-1];
          relXSuper[nrSuper] = fstcol;
          newCC[nrSuper] = newCC[ksup-1];
          ++nrSuper;
        }
    }
    relXSuper[nrSuper] = xsuper[nsuper];
    relXSuper.Resize(nrSuper+1);

    for(Int ksup=1;ksup<=nrSuper;++ksup){
      Int fstcol = relXSuper[ksup-1];
      Int lstcol = relXSuper[ksup]-1;
      for(Int col = fstcol; col<=lstcol;++col){
        supMembership[col-1] = ksup;
        cc[col-1] = newCC[ksup-1] + col-fstcol;
      }
    }
    
    xsuper = relXSuper;
///      //adjust the column counts
///      for(Int col=i-2;col>=i-supsize;--col){
///        cc[col-1] = cc[col]+1;
///      }


}

#endif



#ifdef RELAXED_SNODE
  void SparseMatrixStructure::SymbolicFactorizationRelaxed(ETree& tree,const IntNumVec & cc,const IntNumVec & xsuper,const IntNumVec & SupMembership, IntNumVec & xlindx, IntNumVec & lindx){
    TIMER_START(SymbolicFactorization);


    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call SymbolicFactorization\n" );
    }

    Int nsuper = xsuper.m()-1;




    Int nzbeg = 0;
    //nzend points to the last used slot in lindx
    Int nzend = 0;

    //tail is the end of list indicator (in rchlnk, not mrglnk)
    Int tail = size +1;

    Int head = 0;

    //Array of length nsuper containing the children of 
    //each supernode as a linked list
    IntNumVec mrglnk(nsuper);
    SetValue(mrglnk,0);

    //Array of length n+1 containing the current linked list 
    //of merged indices (the "reach" set)
    IntNumVec rchlnk(size+1);

    //Array of length n used to mark indices as they are introduced
    // into each supernode's index set
    IntNumVec marker(size);
    SetValue(marker,0);



    xlindx.Resize(nsuper+1);

    //Compute the sum of the column count and resize lindx accordingly
    Int nofsub = 1;
    for(Int i =0; i<cc.m();++i){
      nofsub+=cc(i);
    }

    lindx.Resize(nofsub);


    Int point = 1;
    for(Int ksup = 1; ksup<=nsuper; ++ksup){
      Int fstcol = xsuper(ksup-1);
      xlindx(ksup-1) = point;
      point += cc(fstcol-1); 
    } 
    xlindx(nsuper) = point;


//    Int point = 1;
//    for(Int ksup = 1; ksup<=nsuper; ++ksup){
//      Int fstcol = xsuper(ksup-1);
//      Int width = xsuper(ksup) - xsuper(ksup-1);
//      lindx(ksup-1) = point;
//      point += width*cc(fstcol-1); 
//    } 
//    lindx(nsuper) = point;
//
//    IntNumVec tmpLindx = lindx;
//    for(Int ksup = 1; ksup<=nsuper; ++ksup){
//
//    ETree superTree = tree.ToSupernodalETree(xsuper);
//      Int fstcol = xsuper(ksup-1);
//      Int lstcol = xsuper(ksup)-1;
//
//      for(Int row = fstcol; row <=lstcol; ++row){
////        lindx[tmpLindx[ksup-1]++-1]=row;
//      }
//
//
//      for(Int row = fstcol; row <=lstcol; ++row){
//        /* traverse the row subtree for each nonzero in A or AA' */
//          subtree (k, k, Ap, Ai, Anz, SuperMap, Sparent, mark,
//              Flag, Ls, Lpi2) ;
//
//
//
//
//
//
//    Int p, pend, i, si ;
//    p = Ap [j] ;
//    pend = (Anz == NULL) ? (Ap [j+1]) : (p + Anz [j]) ;
//    for ( ; p < pend ; p++)
//    {
//  Int row = Ai [p] ;
//  if (row < k)
//  {
//      /* (i,k) is an entry in the upper triangular part of A or A*F'.
//       * symmetric case:   A(i,k) is nonzero (j=k).
//       * unsymmetric case: A(i,j) and F(j,k) are both nonzero.
//       *
//       * Column i is in supernode si = SuperMap [i].  Follow path from si
//       * to root of supernodal etree, stopping at the first flagged
//       * supernode.  The root of the row subtree is supernode SuperMap[k],
//       * which is flagged already. This traversal will stop there, or it
//       * might stop earlier if supernodes have been flagged by previous
//       * calls to this routine for the same k. */
//      for (si = SuperMap [row-1] ; Flag [si] < mark ; si = Sparent [si])
//      {
//    ASSERT (si <= SuperMap [k]) ;
//    Ls [tmpLindx [si]++] = k ;
//    Flag [si] = mark ;
//      }
//  }
//    }
//      }
//    }



    for(Int ksup = 1; ksup<=nsuper; ++ksup){
      Int fstcol = xsuper(ksup-1);
      Int lstcol = xsuper(ksup)-1;
      Int width = lstcol - fstcol +1;
      Int length = cc(fstcol-1);
      Int knz = 0;
      rchlnk(head) = tail;
      Int jsup = mrglnk(ksup-1);

      //If ksup has children in the supernodal e-tree
      if(jsup>0){
        //copy the indices of the first child jsup into 
        //the linked list, and mark each with the value 
        //ksup.
        Int jwidth = xsuper(jsup)-xsuper(jsup-1);
        Int jnzbeg = xlindx(jsup-1) + jwidth;
        Int jnzend = xlindx(jsup) -1;
        for(Int jptr = jnzend; jptr>=jnzbeg; --jptr){
          Int newi = lindx(jptr-1);
          ++knz;
          marker(newi-1) = ksup;
          rchlnk(newi) = rchlnk(head);
          rchlnk(head) = newi;
        }

        //for each subsequent child jsup of ksup ...
        jsup = mrglnk(jsup-1);
        while(jsup!=0 && knz < length){
          //merge the indices of jsup into the list,
          //and mark new indices with value ksup.

          jwidth = xsuper(jsup)-xsuper(jsup-1);
          jnzbeg = xlindx(jsup-1) + jwidth;
          jnzend = xlindx(jsup) -1;
          Int nexti = head;
          for(Int jptr = jnzbeg; jptr<=jnzend; ++jptr){
            Int newi = lindx(jptr-1);
            Int i;
            do{
              i = nexti;
              nexti = rchlnk(i);
            }while(newi > nexti);

            if(newi < nexti){
#ifdef _DEBUG_
            logfileptr->OFS()<<jsup<<" is a child of "<<ksup<<" and "<<newi<<" is inserted in the structure of "<<ksup<<std::endl;
#endif
              ++knz;
              rchlnk(i) = newi;
              rchlnk(newi) = nexti;
              marker(newi-1) = ksup;
              nexti = newi;
            }
          }
          jsup = mrglnk(jsup-1);
        }
      }

      //structure of a(*,fstcol) has not been examined yet.  
      //"sort" its structure into the linked list,
      //inserting only those indices not already in the
      //list.
      if(knz < length){
        for(Int row = fstcol; row<=lstcol; ++row){
          Int newi = row;
          if(newi > fstcol && marker(newi-1) != ksup){
            //position and insert newi in list and
            // mark it with kcol
            Int nexti = head;
            Int i;
            do{
              i = nexti;
              nexti = rchlnk(i);
            }while(newi > nexti);
            ++knz;
            rchlnk(i) = newi;
            rchlnk(newi) = nexti;
            marker(newi-1) = ksup;
          }
        }


        for(Int col = fstcol; col<=lstcol; ++col){
          Int node = tree.FromPostOrder(col);
          Int knzbeg = colptr(node-1);
          Int knzend = colptr(node)-1;
          for(Int kptr = knzbeg; kptr<=knzend;++kptr){
            Int newi = rowind(kptr-1);
            newi = tree.ToPostOrder(newi);
            if(newi > fstcol && marker(newi-1) != ksup){
              //position and insert newi in list and
              // mark it with kcol
              Int nexti = head;
              Int i;
              do{
                i = nexti;
                nexti = rchlnk(i);
              }while(newi > nexti);
              ++knz;
              rchlnk(i) = newi;
              rchlnk(newi) = nexti;
              marker(newi-1) = ksup;
            }
          }
        }

      } 

      //if ksup has no children, insert fstcol into the linked list.
      if(rchlnk(head) != fstcol){
        rchlnk(fstcol) = rchlnk(head);
        rchlnk(head) = fstcol;
        ++knz;
      }

      assert(knz == cc(fstcol-1));


      //copy indices from linked list into lindx(*).
      nzbeg = nzend+1;
      nzend += knz;
      assert(nzend+1 == xlindx(ksup));
      Int i = head;
      for(Int kptr = nzbeg; kptr<=nzend;++kptr){
        i = rchlnk(i);
        lindx(kptr-1) = i;
      } 

      //if ksup has a parent, insert ksup into its parent's 
      //"merge" list.
      if(length > width){
        Int pcol = lindx(xlindx(ksup-1) + width -1);
        Int psup = SupMembership(pcol-1);
        mrglnk(ksup-1) = mrglnk(psup-1);
        mrglnk(psup-1) = ksup;
      }
    }

    lindx.Resize(nzend+1);

    TIMER_STOP(SymbolicFactorization);
  }
#endif


















  void SparseMatrixStructure::SymbolicFactorization2(ETree& tree,const IntNumVec & cc,const IntNumVec & xsuper,const IntNumVec & SupMembership, IntNumVec & xlindx, IntNumVec & lindx){
    TIMER_START(SymbolicFactorization);


    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call SymbolicFactorization\n" );
    }

    Int nsuper = xsuper.m()-1;

    Int nzbeg = 0;
    //nzend points to the last used slot in lindx
    Int nzend = 0;

    //tail is the end of list indicator (in rchlnk, not mrglnk)
    Int tail = size +1;

    Int head = 0;

    //Array of length nsuper containing the children of 
    //each supernode as a linked list
    IntNumVec mrglnk(nsuper);
    SetValue(mrglnk,0);

    //Array of length n+1 containing the current linked list 
    //of merged indices (the "reach" set)
    IntNumVec rchlnk(size+1);

    //Array of length n used to mark indices as they are introduced
    // into each supernode's index set
    IntNumVec marker(size);
    SetValue(marker,0);



    xlindx.Resize(nsuper+1);

    //Compute the sum of the column count and resize lindx accordingly
    Int nofsub = 1;
    for(Int i =0; i<cc.m();++i){
      nofsub+=cc(i);
    }

    lindx.Resize(nofsub);


    Int point = 1;
    for(Int ksup = 1; ksup<=nsuper; ++ksup){
      Int fstcol = xsuper(ksup-1);
      xlindx(ksup-1) = point;
      point += cc(fstcol-1); 
    } 
    xlindx(nsuper) = point;

    for(Int ksup = 1; ksup<=nsuper; ++ksup){
      Int fstcol = xsuper(ksup-1);
      Int lstcol = xsuper(ksup)-1;
      Int width = lstcol - fstcol +1;
      Int length = cc(fstcol-1);
      Int knz = 0;
      rchlnk(head) = tail;
      Int jsup = mrglnk(ksup-1);

      //If ksup has children in the supernodal e-tree
      if(jsup>0){
        //copy the indices of the first child jsup into 
        //the linked list, and mark each with the value 
        //ksup.
        Int jwidth = xsuper(jsup)-xsuper(jsup-1);
        Int jnzbeg = xlindx(jsup-1) + jwidth;
        Int jnzend = xlindx(jsup) -1;
        for(Int jptr = jnzend; jptr>=jnzbeg; --jptr){
          Int newi = lindx(jptr-1);
          ++knz;
          marker(newi-1) = ksup;
          rchlnk(newi) = rchlnk(head);
          rchlnk(head) = newi;
        }

        //for each subsequent child jsup of ksup ...
        jsup = mrglnk(jsup-1);
        while(jsup!=0 && knz < length){
          //merge the indices of jsup into the list,
          //and mark new indices with value ksup.

          jwidth = xsuper(jsup)-xsuper(jsup-1);
          jnzbeg = xlindx(jsup-1) + jwidth;
          jnzend = xlindx(jsup) -1;
          Int nexti = head;
          for(Int jptr = jnzbeg; jptr<=jnzend; ++jptr){
            Int newi = lindx(jptr-1);
            Int i;
            do{
              i = nexti;
              nexti = rchlnk(i);
            }while(newi > nexti);

            if(newi < nexti){
#ifdef _DEBUG_
            logfileptr->OFS()<<jsup<<" is a child of "<<ksup<<" and "<<newi<<" is inserted in the structure of "<<ksup<<std::endl;
#endif
              ++knz;
              rchlnk(i) = newi;
              rchlnk(newi) = nexti;
              marker(newi-1) = ksup;
              nexti = newi;
            }
          }
          jsup = mrglnk(jsup-1);
        }
      }

#ifdef RELAXED_SNODE
      //structure of a(*,fstcol) has not been examined yet.  
      //"sort" its structure into the linked list,
      //inserting only those indices not already in the
      //list.
      if(knz < length){
        for(Int row = fstcol; row<=lstcol; ++row){
          Int newi = row;
          if(newi > fstcol && marker(newi-1) != ksup){
            //position and insert newi in list and
            // mark it with kcol
            Int nexti = head;
            Int i;
            do{
              i = nexti;
              nexti = rchlnk(i);
            }while(newi > nexti);
            ++knz;
            rchlnk(i) = newi;
            rchlnk(newi) = nexti;
            marker(newi-1) = ksup;
          }
        }


        for(Int col = fstcol; col<=lstcol; ++col){
          Int node = tree.FromPostOrder(col);
          Int knzbeg = colptr(node-1);
          Int knzend = colptr(node)-1;
          for(Int kptr = knzbeg; kptr<=knzend;++kptr){
            Int newi = rowind(kptr-1);
            newi = tree.ToPostOrder(newi);
            if(newi > fstcol && marker(newi-1) != ksup){
              //position and insert newi in list and
              // mark it with kcol
              Int nexti = head;
              Int i;
              do{
                i = nexti;
                nexti = rchlnk(i);
              }while(newi > nexti);
              ++knz;
              rchlnk(i) = newi;
              rchlnk(newi) = nexti;
              marker(newi-1) = ksup;
            }
          }
        }

      } 
#else
      //structure of a(*,fstcol) has not been examined yet.  
      //"sort" its structure into the linked list,
      //inserting only those indices not already in the
      //list.
      if(knz < length){
        Int node = tree.FromPostOrder(fstcol);
        Int knzbeg = colptr(node-1);
        Int knzend = colptr(node)-1;
        for(Int kptr = knzbeg; kptr<=knzend;++kptr){
          Int newi = rowind(kptr-1);
          newi = tree.ToPostOrder(newi);
          if(newi > fstcol && marker(newi-1) != ksup){
            //position and insert newi in list and
            // mark it with kcol
            Int nexti = head;
            Int i;
            do{
              i = nexti;
              nexti = rchlnk(i);
            }while(newi > nexti);
            ++knz;
            rchlnk(i) = newi;
            rchlnk(newi) = nexti;
            marker(newi-1) = ksup;
          }
        }
      }
#endif

      //if ksup has no children, insert fstcol into the linked list.
      if(rchlnk(head) != fstcol){
        rchlnk(fstcol) = rchlnk(head);
        rchlnk(head) = fstcol;
        ++knz;
      }

      assert(knz == cc(fstcol-1));


      //copy indices from linked list into lindx(*).
      nzbeg = nzend+1;
      nzend += knz;
      assert(nzend+1 == xlindx(ksup));
      Int i = head;
      for(Int kptr = nzbeg; kptr<=nzend;++kptr){
        i = rchlnk(i);
        lindx(kptr-1) = i;
      } 

      //if ksup has a parent, insert ksup into its parent's 
      //"merge" list.
      if(length > width){
        Int pcol = lindx(xlindx(ksup-1) + width -1);
        Int psup = SupMembership(pcol-1);
        mrglnk(ksup-1) = mrglnk(psup-1);
        mrglnk(psup-1) = ksup;
      }
    }

    lindx.Resize(nzend+1);

    TIMER_STOP(SymbolicFactorization);
  }







//NOT WORKING
  void SparseMatrixStructure::SymbolicFactorization(ETree& tree,const IntNumVec & cc,const IntNumVec & xsuper, IntNumVec & xlindx, IntNumVec & lindx){
    TIMER_START(SymbolicFactorization);


    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call SymbolicFactorization\n" );
    }

    std::vector<std::set<Int> > sets;
    sets.resize(xsuper.m(),std::set<Int>());

    std::vector<IntNumVec > LIs;
    LIs.resize(xsuper.m());

    Int lindxCnt = 0;
    for(Int I=1;I<xsuper.m();I++){
      Int fi = tree.FromPostOrder(xsuper(I-1));
      Int width = xsuper(I)-xsuper(I-1);
      Int length = cc(fi-1);

      IntNumVec & LI = LIs[I-1];


      //Initialize LI with nnz struct of A_*fi
      Int begin = colptr(fi-1);
      Int end = colptr(fi);

      Int * start = &rowind(begin-1); 
      Int * stop = (end>=rowind.m())?&rowind(rowind.m()-1)+1:&rowind(end-1); 
      //find the diagonal block
      start=std::find(start,stop,fi);


      LI.Resize(stop-start);

      std::copy(start,stop,&LI(0));

      //     logfileptr->OFS()<<"L"<<I<<"<- A_*,fi: ";
      //     for(int i=0;i<LI.m();i++){logfileptr->OFS()<<LI(i)<< " ";}
      //     logfileptr->OFS()<<std::endl;

      LI = tree.ToPostOrder(LI);


      //     logfileptr->OFS()<<"PO L"<<I<<"<- A_*,fi: ";
      //     for(int i=0;i<LI.m();i++){logfileptr->OFS()<<LI(i)<< " ";}
      //     logfileptr->OFS()<<std::endl;




      std::set<Int> & SI = sets[I-1];
      for(std::set<Int>::iterator it = SI.begin(); it!=SI.end(); it++){
        Int K = *it;
        IntNumVec & LK = LIs[K-1];

        //        logfileptr->OFS()<<"merging "<<I<<" with "<<K<<std::endl;
        //LI = LI U LK \ K
        IntNumVec Ltmp(LI.m()+LK.m()-1);
        std::copy(&LI(0),&LI(LI.m()-1)+1,&Ltmp(0));


        if(LK.m()>1){

          //Be careful to not insert duplicates !
          Int firstidx =1;
          for(Int i =1;i<LK.m();i++){
            if(LK[i]>fi ){
              firstidx = i+1;
              break;
            }
          }
          Int * end = std::set_union(&LI(0),&LI(LI.m()-1)+1,&LK(firstidx-1),&LK(LK.m()-1)+1,&Ltmp(0));
          Ltmp.Resize(end - &Ltmp(0));
          LI = Ltmp;
        }


      }


      lindxCnt += LI.m();
#ifdef _DEBUG_
        logfileptr->OFS()<<"I:"<<I<<std::endl<<"LI:"<<LI<<std::endl;
#endif

      if(length>width  ){
#ifdef _DEBUG_
        logfileptr->OFS()<<"I:"<<I<<std::endl<<"LI:"<<LI<<std::endl;
#endif
        Int i = LI(width);
        //        logfileptr->OFS()<<I<<" : col "<<i<<" is the next to be examined. width="<<width<<" length="<<length<<std::endl;

        Int J = I+1;
        for(J = I+1;J<=xsuper.m()-1;J++){
          Int fc = xsuper(J-1);
          Int lc = xsuper(J)-1;
          //          logfileptr->OFS()<<"FC = "<<fc<<" vs "<<i<<std::endl;
          //          logfileptr->OFS()<<"LC = "<<lc<<" vs "<<i<<std::endl;
          if(fc <=i && lc >= i){
            //            logfileptr->OFS()<<I<<" : col "<<i<<" found in snode "<<J<<std::endl;
            break;
          }
        } 

        //        logfileptr->OFS()<<I<<" : col "<<i<<" is in snode "<<J<<std::endl;
        std::set<Int> & SJ = sets[J-1];
        SJ.insert(I);
        //        logfileptr->OFS()<<"S"<<J<<" U {"<<I<<"}"<<std::endl; 

      }

    }  

    Int nsuper = xsuper.m()-1;
    Int totNnz = 1;
    for(Int I=1;I<xsuper.m();I++){
      Int fc = xsuper(I-1);
      Int lc = xsuper(I)-1;
      for(Int i=fc;i<=lc;i++){
        totNnz+=cc(i-1);
      }

    }

    lindx.Resize(lindxCnt);
    xlindx.Resize(nsuper+1);
    Int head = 1;

    //    logfileptr->OFS()<<"There are "<<lindxCnt<<" slots in lindx"<<std::endl;
    for(Int I=1;I<=nsuper;I++){
      Int fi = tree.FromPostOrder(xsuper(I-1));
      IntNumVec & LI = LIs[I-1];
      xlindx(I-1)=head;


      //        logfileptr->OFS()<<"PO L"<<I<<":";
      //        for(int i=0;i<LI.m();i++){logfileptr->OFS()<<LI(i)<< " ";}
      //        logfileptr->OFS()<<std::endl;

      //      logfileptr->OFS()<<"Copying "<<LI.m()<<" elem into lindx("<<head-1<<")"<<std::endl;
      for(Int i=0;i<LI.m();++i){
        if(LI(i)!=0){
          lindx(head-1) = LI(i); 
          head++;
        }
      }
      //std::copy(&LI(0),LI.Data()+LI.m(),&(lindx(head-1)));
      //head+=LI.m();//cc(fi-1);
    }
    //    lindx.Resize(head-1);
    xlindx(nsuper) = head;

    TIMER_STOP(SymbolicFactorization);
  }

//FIXME correct these methods
//Return the row structure in the permuted matrix
void SparseMatrixStructure::GetARowStruct(const ETree & etree, const Int iPORow, std::vector<Int> & rowStruct){
  TIMER_START(SparseMatrixStructure::GetARowStruct);
    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call GetARowStruct\n" );
    }

  for(Int iPOCurCol = 1; iPOCurCol<iPORow;++iPOCurCol){
    Int iCurCol = etree.FromPostOrder(iPOCurCol);

    Int iFirstRowPtr = colptr(iCurCol-1);
    Int iLastRowPtr = colptr(iCurCol)-1;

//    logfileptr->OFS()<<iFirstRowPtr<<" to "<<iLastRowPtr<<std::endl;
    for(Int iCurRowPtr=iFirstRowPtr ;iCurRowPtr<=iLastRowPtr;++iCurRowPtr){
      Int iPOCurRow = etree.ToPostOrder(rowind(iCurRowPtr-1));
      if(iPOCurRow == iPORow){
        rowStruct.push_back(iPOCurCol);
//        logfileptr->OFS()<<"A("<<iPORow<<","<<iPOCurCol<<") is nz "<<std::endl;
      }

      if(iPOCurRow >= iPORow){
        break;
      }
    }
  }
  TIMER_STOP(SparseMatrixStructure::GetARowStruct);
}

void SparseMatrixStructure::GetLRowStruct(const ETree & etree, const Int iPORow, const std::vector<Int> & ARowStruct, std::set<Int> & LRowStruct){

  TIMER_START(SparseMatrixStructure::GetLRowStruct);
    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call GetLRowStruct\n" );
    }
//  LRowStruct.clear();
  for(Int i = 0; i<ARowStruct.size();++i){
    Int iCurNode = ARowStruct[i];
    //tracing from iCurRow to iRow;
//    logfileptr->OFS()<<"Starting from node "<<iCurNode<<std::endl;
    if(iCurNode==iPORow){
//      logfileptr->OFS()<<"Adding node "<<iCurNode<<std::endl;
      LRowStruct.insert(iCurNode);
    }
    else{
      while(iCurNode != iPORow && etree.PostParent(iCurNode-1) != 0){
//        logfileptr->OFS()<<"Adding node "<<iCurNode<<std::endl;
        LRowStruct.insert(iCurNode);
        iCurNode = etree.PostParent(iCurNode-1);
//        logfileptr->OFS()<<"Parent is "<<iCurNode<<std::endl;
      }
    }
  } 
  TIMER_STOP(SparseMatrixStructure::GetLRowStruct);
}



void SparseMatrixStructure::GetSuperARowStruct(const ETree & etree, const IntNumVec & Xsuper, const IntNumVec & SupMembership, const Int iSupNo, std::vector<Int> & SuperRowStruct){
  TIMER_START(SpStruct_GetSuperARowStruct);
//  TIMER_START(SparseMatrixStructure::GetSuperARowStruct);
    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call GetARowStruct\n" );
    }

  Int first_col = Xsuper[iSupNo-1];
  Int last_col = Xsuper[iSupNo]-1;

  for(Int iPOCurCol = 1; iPOCurCol<first_col;++iPOCurCol){
    Int iCurCol = etree.FromPostOrder(iPOCurCol);
    Int iCurSupno = SupMembership(iPOCurCol-1);
    

    Int iFirstRowPtr = colptr(iCurCol-1);
    Int iLastRowPtr = colptr(iCurCol)-1;

//    logfileptr->OFS()<<iFirstRowPtr<<" to "<<iLastRowPtr<<std::endl;
    for(Int iCurRowPtr=iFirstRowPtr ;iCurRowPtr<=iLastRowPtr;++iCurRowPtr){
      Int iPOCurRow = etree.ToPostOrder(rowind(iCurRowPtr-1));
      if(iPOCurRow >= first_col && iPOCurRow <= last_col){
//        SuperRowStruct.push_back(iPOCurCol);
        SuperRowStruct.push_back(iCurSupno);
//        logfileptr->OFS()<<"A("<<iPORow<<","<<iPOCurCol<<") is nz "<<std::endl;
      }

      if(iPOCurRow > last_col){
        break;
      }
    }
  }
//  TIMER_STOP(SparseMatrixStructure::GetSuperARowStruct);
  TIMER_STOP(SpStruct_GetSuperARowStruct);
}


void SparseMatrixStructure::GetSuperLRowStruct(const ETree & etree, const IntNumVec & Xsuper, const IntNumVec & SupMembership, const Int iSupNo, std::set<Int> & SuperLRowStruct){

  TIMER_START(SpStruct_GetSuperLRowStruct);
    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call GetSuperLRowStruct\n" );
    }

//  SuperLRowStruct.clear();

    //Get A row struct
    std::vector<Int> SuperARowStruct;
    GetSuperARowStruct(etree, Xsuper,SupMembership, iSupNo, SuperARowStruct);



#ifdef _DEBUG_
      logfileptr->OFS()<<"Row structure of A of Supernode "<<iSupNo<<" is ";
      for(std::vector<Int>::iterator it = SuperARowStruct.begin(); it != SuperARowStruct.end(); ++it){
        logfileptr->OFS()<<*it<<" ";
      }
      logfileptr->OFS()<<std::endl;
#endif

  Int first_col = Xsuper[iSupNo-1];
  Int last_col = Xsuper[iSupNo]-1;

  for(Int i = 0; i<SuperARowStruct.size();++i){
    Int iCurNode = SuperARowStruct[i];
    //tracing from iCurRow to iRow;
#ifdef _DEBUG_
    logfileptr->OFS()<<"Starting from node "<<iCurNode<<std::endl;
#endif
    //if(iCurNode==iPORow){
    if(iCurNode >= first_col && iCurNode <= last_col){
#ifdef _DEBUG_
      logfileptr->OFS()<<"Adding node "<<iCurNode<<std::endl;
#endif
      SuperLRowStruct.insert(iCurNode);
    }
    else{
      while( iCurNode != first_col && etree.PostParent(iCurNode-1) != 0){
//      while( !(iCurNode >= first_col && iCurNode <= last_col) && etree.PostParent(iCurNode-1) != 0){
#ifdef _DEBUG_
        logfileptr->OFS()<<"Adding node "<<iCurNode<<std::endl;
#endif
        SuperLRowStruct.insert(iCurNode);
        iCurNode = etree.PostParent(iCurNode-1);
#ifdef _DEBUG_
        logfileptr->OFS()<<"Parent is "<<iCurNode<<std::endl;
#endif
      }
    }
  } 
  TIMER_STOP(SpStruct_GetSuperLRowStruct);
}



}

