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


    logfileptr->OFS()<<"children "<<children.m()<<std::endl;
    for(Int i = 0; i<children.m();i++){
      logfileptr->OFS()<<children(i)<<" ";
    }
    logfileptr->OFS()<<std::endl;

    Int nsuper = 1;
    xsuper.Resize(2*size);
    xsuper(nsuper-1) = 1;
    //SetValue(xsuper,I_ONE);
    for(Int i =2; i<=size;i++){
      //logfileptr->OFS()<<"Column "<<i<<" has "<<children(i-1)<<" children, "<<cc(i-1)<<" nnz vs "<<cc(i-2)<<" in the prev col"<<std::endl;

      if(children(i-1)!=1 || cc(i-1) != (cc(i-2)-1)){
//        logfileptr->OFS()<<"Col "<<i<<" and "<<i-1<<" are not in the same snode"<<std::endl; 
        nsuper++;
        xsuper(nsuper-1) = i;
      }
//      else{
//        logfileptr->OFS()<<"Col "<<i<<" and "<<i-1<<" are in the same snode"<<std::endl; 
//      }
    }
      nsuper++;
      xsuper(nsuper-1) = size+1;

    xsuper.Resize(nsuper);


    

    TIMER_STOP(FindSupernodes);
  }


  template <class F> void DistSparseMatrix<F>::SymbolicFactorization(ETree& tree,const IntNumVec & cc,const IntNumVec & xsuper){
    TIMER_START(SymbolicFactorization);

logfileptr->OFS()<<"xsuper:"<<xsuper<<std::endl;


//    Int nsuper = xsuper.m();
//    IntNumVec snode(size+1);
//    IntNumVec marker(size+1);
//    IntNumVec mrglnk(nsuper+1);
//    IntNumVec rchlnk(size+1);
//
//    IntNumVec xlnz(size+1);
//    IntNumVec lindx(size*nsuper+1);
//    IntNumVec xlindx(nsuper+1);
//
//    Int nzbeg = 0;
//    Int nzend = 0;
//    Int head = 0;
//    Int tail = size + 1;
//    SetValue(marker,I_ZERO);
//    Int pointl = 1;
//    for(Int ksup=1;ksup<nsuper;ksup++){
//            Int fstcol = xsuper(ksup-1);
//            Int lstcol = xsuper(ksup) - 1;
//            for(Int jcol=fstcol;jcol<=lstcol;jcol++){
//                xlnz(jcol-1) = pointl;
//                pointl = pointl + cc(fstcol-1);
//            }
//    }
//    xlnz(size) = pointl;
//    Int point = 1;
//    for(Int ksup=1;ksup<nsuper;ksup++){
//            mrglnk(ksup-1) = 0;
//            Int fstcol = xsuper(ksup-1);
//            xlindx(ksup-1) = point;
//            point = point + cc(fstcol-1);
//    }
//    xlindx(nsuper) = point;
//
//
//    logfileptr->OFS()<<"xlnz"<<":";
//    for(int i=0;i<xlnz.m();i++){logfileptr->OFS()<<xlnz(i)<< " ";}
//    logfileptr->OFS()<<std::endl;
//
//        logfileptr->OFS()<<"xlindx"<<":";
//        for(int i=0;i<xlindx.m();i++){logfileptr->OFS()<<xlindx(i)<< " ";}
//        logfileptr->OFS()<<std::endl;
//
//
//    IntNumVec & xadj = this->Global_.colptr;
//    IntNumVec & adjncy = this->Global_.rowind;
//
//
//
//
//
//
//    for(Int ksup=1;ksup<nsuper;ksup++){
//            //initializations ...
//            //    fstcol : first column of supernode ksup.
//            //    lstcol : last column of supernode ksup.
//            //    knz    : will count the nonzeros of l in column kcol.
//            //    rchlnk : initialize empty index list for kcol.
//            Int fstcol = xsuper(ksup-1);
//            Int lstcol = xsuper(ksup) - 1;
//            Int width  = lstcol - fstcol + 1;
//            Int length = cc(fstcol-1);
//            Int knz = 0;
//            rchlnk(head) = tail;
//            Int jsup = mrglnk(ksup-1);
//
//            Int newi; Int nexti; Int i;
//            //if ksup has children in the supernodal e-tree ...
//            if  ( jsup > 0 ) {
//                //copy the indices of the first child jsup into 
//                //the linked list, and mark each with the value 
//                //ksup.
//                Int jwidth = xsuper(jsup) - xsuper(jsup-1);
//                Int jnzbeg = xlindx(jsup-1) + jwidth;
//                Int jnzend = xlindx(jsup) - 1;
//                for(Int jptr = jnzend;jptr>=jnzbeg;jptr--){
//                    newi = lindx(jptr-1);
//                    knz = knz+1;
//                    marker(newi-1) = ksup;
//                    rchlnk(newi) = rchlnk(head);
//                    rchlnk(head) = newi;
//                }
////  200           continue
//                //for each subsequent child jsup of ksup ...
//                jsup = mrglnk(jsup-1);
////  300           continue
//                if  ( jsup != 0  &&  knz < length ) {
//                    //merge the indices of jsup into the list,
//                    //and mark new indices with value ksup.
//                    Int jwidth = xsuper(jsup) - xsuper(jsup-1);
//                    Int jnzbeg = xlindx(jsup-1) + jwidth;
//                    Int jnzend = xlindx(jsup) - 1;
//                    nexti = head;
//                    for(Int jptr = jnzbeg;jptr<=jnzend;jptr++){
//                        newi = lindx(jptr-1);
//                            i = nexti;
//                            nexti = rchlnk(i);
//                            if  ( newi > nexti )  {continue;}
//                        if  ( newi < nexti ) {
//                            knz = knz+1;
//                            rchlnk(i) = newi;
//                            rchlnk(newi) = nexti;
//                            marker(newi-1) = ksup;
//                            nexti = newi;
//                        }
//                    }
//                    jsup = mrglnk(jsup-1);
//                    continue;//go to 300
//                }
//            }
//            //structure of a(*,fstcol) has not been examined yet.  
//            //"sort" its structure into the linked list,
//            //inserting only those indices not already in the
//            //list.
//            if  ( knz < length )  {
//                Int node = tree.ToPostOrder(fstcol-1);
//                Int knzbeg = xadj(node-1);
//                Int knzend = xadj(node) - 1;
//                for(Int kptr= knzbeg; kptr<=knzend;kptr++){
//                    Int newi = adjncy(kptr-1);
//                    newi = tree.FromPostOrder(newi-1);
//                    if  ( newi > fstcol  && marker(newi-1) != ksup ) {
//                        //position and insert newi in list
//                        //and mark it with kcol.
//                        Int nexti = head;
//                        Int i = nexti;
//                        nexti = rchlnk(i);
//                        if  ( newi > nexti ) {continue;}
//                        knz = knz + 1;
//                        rchlnk(i) = newi;
//                        rchlnk(newi) = nexti;
//                        marker(newi-1) = ksup;
//                    }
//              }
//            }
//            //if ksup has no children, insert fstcol into the linked list.
//            if  ( rchlnk(head) != fstcol ){
//                rchlnk(fstcol) = rchlnk(head);
//                rchlnk(head) = fstcol;
//                knz = knz + 1;
//            }
//
//            //copy indices from linked list into lindx(*).
//            nzbeg = nzend + 1;
//            nzend = nzend + knz;
//            if  ( nzend+1 != xlindx(ksup) ) {/*problem*/}
//            Int li = head;
//            for(Int kptr= nzbeg; kptr<=nzend;kptr++){
//                li = rchlnk(li);
//                lindx(kptr-1) = li;
//            }
//            //if ksup has a parent, insert ksup into its parent's 
//            //"merge" list.
//            if  ( length > width ) {
//                Int pcol = lindx ( xlindx(ksup-1) + width -1 );
//                Int psup = snode(pcol-1);
//                mrglnk(ksup-1) = mrglnk(psup-1);
//                mrglnk(psup-1) = ksup;
//            }
//      }
//
//
//        logfileptr->OFS()<<"xlindx"<<":";
//        for(int i=0;i<xlindx.m();i++){logfileptr->OFS()<<xlindx(i)<< " ";}
//        logfileptr->OFS()<<std::endl;
//
//    return;












    IntNumVec & rowind = this->Global_.rowind;
    IntNumVec & colptr = this->Global_.colptr;

//    IntNumVec marker(xsuper.m());
//    for(Int I=1;I<xsuper.m();I++){
//      marker(I-1)=I;
//    }
  


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
//      logfileptr->OFS()<<"fi="<<fi<<" begin="<<begin<<" end="<<end<<std::endl;
      


      Int * start = &rowind(begin-1); 
      Int * stop = (end-1<rowind.m())?&rowind(end-1):&rowind(rowind.m()-1)+1; 
      //find the diagonal block
      start=std::find(start,stop,fi);


      LI.Resize(stop-start);
      
      std::copy(start,stop,&LI(0));
      LI = tree.ToPostOrder(LI);

      std::set<Int> & SI = sets[I-1];
      for(std::set<Int>::iterator it = SI.begin(); it!=SI.end(); it++){
        Int K = *it;
        IntNumVec & LK = LIs[K-1];

        logfileptr->OFS()<<"merging "<<I<<" with "<<K<<std::endl;
        //LI = LI U LK \ K
        IntNumVec Ltmp(LI.m()+LK.m()-1);
        std::copy(&LI(0),&LI(LI.m()-1)+1,&Ltmp(0));

        
        if(LK.m()>1){
          Int head = LI.m();
          for(Int i =1;i<LK.m();i++){
            //insert element from LK
            if(LK[i]>I){
              Ltmp[head]=LK[i];
              head++;
            }
          }
          Ltmp.Resize(head);
          std::sort(&Ltmp(0),&Ltmp(Ltmp.m()-1)+1);
          
          //Int * end = std::set_union(&LI(0),&LI(LI.m()-1)+1,&LK(1),&LK(LK.m()-1)+1,&Ltmp(0));
          //Ltmp.Resize(end - &Ltmp(0));
          LI = Ltmp;
        }


      }
     
      
      lindxCnt += LI.m();

      if(length>width){
        logfileptr->OFS()<<"I:"<<I<<std::endl<<"LI:"<<LI<<std::endl;
        Int i = LI(width);
        logfileptr->OFS()<<I<<" : col "<<i<<" is the next to be examined"<<std::endl;

        Int J = I+1;
        for(J = I+1;J<=xsuper.m()-1;J++){
          Int fc = xsuper(J-1);
          Int lc = xsuper(J)-1;
          logfileptr->OFS()<<"FC = "<<fc<<" vs "<<i<<std::endl;
          logfileptr->OFS()<<"LC = "<<lc<<" vs "<<i<<std::endl;
          if(fc <=i && lc >= i){
            logfileptr->OFS()<<I<<" : col "<<i<<" found in snode "<<J<<std::endl;
            break;
          }
        } 

        logfileptr->OFS()<<I<<" : col "<<i<<" is in snode "<<J<<std::endl;
        std::set<Int> & SJ = sets[J-1];
        SJ.insert(I);

      }

    }  

    Int nsuper = xsuper.m()-1;
    IntNumVec xlnz(size+1);
    Int totNnz = 0;
    for(Int i=1;i<=cc.m();i++){xlnz(i-1)=totNnz+1; totNnz+=cc(i-1);}
    xlnz(size)=totNnz+1;

    DblNumVec lnz(totNnz+1);

    IntNumVec lindx(lindxCnt+1);
    IntNumVec xlindx(nsuper+1);
    Int head = 1;
    for(Int I=1;I<=nsuper;I++){
      Int fi = tree.FromPostOrder(xsuper(I-1));
      IntNumVec & LI = LIs[I-1];
      xlindx(I-1)=head;
      std::copy(&LI(0),&LI(LI.m()-1)+1,&(lindx(head-1)));
      head+=cc(fi-1);
    }
    xlindx(nsuper) = head;

 
        logfileptr->OFS()<<"xlindx"<<":";
        for(int i=0;i<xlindx.m();i++){logfileptr->OFS()<<xlindx(i)<< " ";}
        logfileptr->OFS()<<std::endl;

        logfileptr->OFS()<<"lindx"<<lindxCnt<<" :";
        for(int i=0;i<lindx.m();i++){logfileptr->OFS()<<lindx(i)<< " ";}
        logfileptr->OFS()<<std::endl;

        logfileptr->OFS()<<"xlnz"<<":";
        for(int i=0;i<xlnz.m();i++){logfileptr->OFS()<<xlnz(i)<< " ";}
        logfileptr->OFS()<<std::endl;


    //parsing the data structure
      
    for(Int I=1;I<xsuper.m();I++){
          Int fc = xsuper(I-1);
          Int lc = xsuper(I)-1;
          Int fi = xlindx(I-1);

          logfileptr->OFS()<<"FC = "<<fc<<std::endl;
          logfileptr->OFS()<<"LC = "<<lc<<std::endl;

          for(Int i = fc;i<=lc;i++){
            Int fnz = xlnz(i-1);
            Int lnz = xlnz(i)-1;

          logfileptr->OFS()<<"i = "<<i<<" FNZ = "<<fnz<<std::endl;
          logfileptr->OFS()<<"i = "<<i<<" LNZ = "<<lnz<<std::endl;
            //diag is lnz(fnz)
            logfileptr->OFS()<<"Diag element of "<<i<<" is "<<"lnz("<<fnz-1<<")"<<std::endl;
            Int idx = fi;
            for(Int s = fnz+1;s<=lnz;s++){
              idx++;

              //lnz(s) contains an off-diagonal nonzero entry in row lindx(idx)
              logfileptr->OFS()<<"L("<<lindx(idx-1)<<","<<i<<") = "<<"lnz("<<s-1<<")"<<std::endl;
            }
            fi++;

          }
    }




    TIMER_STOP(SymbolicFactorization);
  }



}


#endif // _DIST_SPARSE_MATRIX_IMPL_HPP_
