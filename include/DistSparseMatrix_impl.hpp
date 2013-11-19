/// @file sparse_matrix_impl.hpp
/// @brief Sparse matrix and Distributed sparse matrix in compressed
/// column format.
/// @author Lin Lin
/// @date 2012-11-10
#ifndef _DIST_SPARSE_MATRIX_IMPL_HPP_
#define _DIST_SPARSE_MATRIX_IMPL_HPP_

#include "Environment.hpp"
#include "DistSparseMatrix.hpp"
#include "utility.hpp"
#include "SuperNode.hpp"

#include <vector>

namespace LIBCHOLESKY{

  template <class F> void DistSparseMatrix<F>::ToGlobalStruct(){
    // test if structure hasn't been allocated yet
    if(globalAllocated){
      return;
    }

    Int numProcs; MPI_Comm_size(comm,&numProcs);
    Int mpirank;  MPI_Comm_rank(comm, &mpirank);


    Global_.rowind.Resize(nnz);
    Global_.colptr.Resize(size+1);


    /* Allgatherv for row indices. */ 
    IntNumVec prevnz(numProcs);
    IntNumVec rcounts(numProcs);
    MPI_Allgather(&Local_.nnz, 1, MPI_INT, rcounts.Data(), 1, MPI_INT, comm);

    prevnz[0] = 0;
    for (Int i = 0; i < numProcs-1; ++i) { prevnz[i+1] = prevnz[i] + rcounts[i]; } 

    MPI_Allgatherv(Local_.rowind.Data(), Local_.nnz, MPI_INT, Global_.rowind.Data(),rcounts.Data(), prevnz.Data(), MPI_INT, comm); 

    /* Allgatherv for colptr */
    // Compute the number of columns on each processor
    Int numColFirst = size / numProcs;
    SetValue( rcounts, numColFirst );
    rcounts[numProcs-1] = size - numColFirst * (numProcs-1);  // Modify the last entry	


    IntNumVec rdispls(numProcs);
    rdispls[0] = 0;
    for (Int i = 0; i < numProcs-1; ++i) { rdispls[i+1] = rdispls[i] + rcounts[i]; } 


    MPI_Allgatherv(Local_.colptr.Data(), Local_.colptr.m()-1, MPI_INT, Global_.colptr.Data(),
        rcounts.Data(), rdispls.Data(), MPI_INT, comm);


    /* Recompute column pointers. */
    for (Int p = 1; p < numProcs; p++) {
      Int idx = rdispls[p];
      for (Int j = 0; j < rcounts[p]; ++j) Global_.colptr[idx++] += prevnz[p];
    }

    Global_.colptr(size)=nnz+1;
    Global_.nnz = nnz;


    globalAllocated = true;



  }







  template <class F> void DistSparseMatrix<F>::ConstructETree(ETree & tree){
    TIMER_START(ConstructETree);

    if(!globalAllocated){
      ToGlobalStruct();
    }


    tree.n_ = size;
    tree.parent_.Resize(size);
    SetValue(tree.parent_,I_ZERO );

    ETree::DisjointSet sets;
    sets.Initialize(size);

    Int cset,rset,rroot,row;
    for (Int col = 1; col <= size; col++) {
      tree.parent_(col-1)=col; //1 based indexes
      cset = sets.makeSet (col);
      sets.Root(cset-1) = col;
      tree.parent_(col-1) = size; 

      for (Int p = Global_.colptr(col-1); p < Global_.colptr(col); p++) {
        row = Global_.rowind(p-1);
        if (row >= col) continue;

        rset = sets.find(row);
        rroot = sets.Root(rset-1);

        if (rroot != col) {
          tree.parent_(rroot-1) = col;
          cset = sets.link(cset, rset);
          sets.Root(cset-1) = col;
        }
      }

    }



    TIMER_STOP(ConstructETree);
  }




  template <class F> void DistSparseMatrix<F>::GetLColRowCount(ETree & tree, IntNumVec & cc, IntNumVec & rc){
    TIMER_START(GetColRowCount);
    //The tree need to be postordered
    if(!tree.isPostOrdered_){
      TIMER_START(PostOrder);
      tree.PostOrderTree();
      TIMER_START(PostOrder);
    }


    TIMER_START(Initialize_Data);
    //cc first contains the delta
    cc.Resize(size);
    //Compute size of subtrees
    IntNumVec treeSize(size);
    SetValue(treeSize,I_ONE);
    IntNumVec level(size);
    level(size-1)=1;
    for(Int vertex = 1; vertex<=size-1; vertex++){
//      Int curParent = tree.postParent_(vertex-1);
      Int curParent = tree.PostParent(vertex);
      treeSize(curParent-1)+=treeSize(vertex-1);
//      level(size-vertex-1) = level(tree.postParent_(size-vertex-1)-1)+1;
      level(size-vertex-1) = level(tree.PostParent(size-vertex)-1)+1;
      if(treeSize(vertex-1)==1){
        cc(vertex-1)=1;
      }
      else{
        cc(vertex-1)=0;
      }
    }

      if(treeSize(size-1)==1){
        cc(size-1)=1;
      }
      else{
        cc(size-1)=0;
      }



    IntNumVec prevLeaf(size);
    SetValue(prevLeaf,I_ZERO);
    IntNumVec prevNz(size);
    SetValue(prevNz,I_ZERO);

    rc.Resize(size);
    SetValue(rc,I_ONE);

    ETree::DisjointSet sets;
    sets.Initialize(size);
    for(Int vertex = 1; vertex<=size; vertex++){
      Int cset = sets.makeSet (vertex);
      sets.Root(cset-1)=vertex;
    }


    TIMER_STOP(Initialize_Data);

    TIMER_START(Compute_Col_Row_Count);
    for(Int col=1; col<=size; col++){
      Int cset;

//      Int colPar = tree.postParent_(col-1);
      Int colPar = tree.PostParent(col);
      if (col<size){
        cc(colPar-1)--;
      }

      Int oCol = tree.FromPostOrder(col);
      for (Int i = Global_.colptr(oCol-1); i < Global_.colptr(oCol); i++) {
        Int row = tree.ToPostOrder(Global_.rowind(i-1));
        if (row > col){
          Int k = prevNz(row-1);
          if(k< col - treeSize(col-1) +1){

            Int p = prevLeaf(row-1);
            cc(col-1)++;
            if(p==0){
              rc(row-1)+=level(col-1)-level(row-1);
            }
            else{
              TIMER_START(Get_LCA);
              Int pset = sets.find(p);
              Int q = sets.Root(pset-1);
              TIMER_STOP(Get_LCA);
              rc(row-1)+= level(col-1) - level(q-1);
              cc(q-1)--;
            }
            prevLeaf(row-1)=col;
          }
          prevNz(row-1)=col;
        }
      }

      TIMER_START(Merge_LCA);
      //merge col and parent sets (for lca computation)
      sets.Union(col,colPar);
      TIMER_STOP(Merge_LCA);

    }
    TIMER_STOP(Compute_Col_Row_Count);

    //convert delta to col count
    for(Int col=1; col<=size-1; col++){
//      Int parent = tree.postParent_(col-1);
      Int parent = tree.PostParent(col);
      cc(parent-1)+= cc(col-1);
    }

    TIMER_STOP(GetColRowCount);
  }



  template <class F> void DistSparseMatrix<F>::FindSupernodes(ETree& tree, IntNumVec & cc, IntNumVec & xsuper){
    TIMER_START(FindSupernodes);


    IntNumVec children(size);
    SetValue(children,I_ZERO);
    for(Int col=1; col<=size-1; col++){
//      children(tree.postParent_(col-1)-1)++;
      children(tree.PostParent(col)-1)++;
    }


    statusOFS<<"children "<<children.m()<<std::endl;
    for(Int i = 0; i<children.m();i++){
      statusOFS<<children(i)<<" ";
    }
    statusOFS<<std::endl;

    Int nsuper = 1;
    xsuper.Resize(2*size);
    SetValue(xsuper,I_ONE);
    for(Int i =2; i<=size;i++){
//      statusOFS<<"Column "<<i<<" has "<<children(i-1)<<" children, "<<cc(i-1)<<" nnz vs "<<cc(i-2)<<" in the prev col"<<std::endl;

      if(children(i-1)!=1 || cc(i-1) != (cc(i-2)-1)){
//        statusOFS<<"Col "<<i<<" and "<<i-1<<" are not in the same snode"<<std::endl; 
        nsuper++;
        xsuper(nsuper-1) = i;
      }
    }
      nsuper++;
      xsuper(nsuper-1) = size+1;

    xsuper.Resize(nsuper);


    

    TIMER_STOP(FindSupernodes);
  }






}


#endif // _DIST_SPARSE_MATRIX_IMPL_HPP_
