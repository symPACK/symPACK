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
#include "SparseMatrixStructure.hpp"

#include <vector>

#ifdef UPCXX
#include <upcxx.h>
#endif
namespace LIBCHOLESKY{






  template <class F> SparseMatrixStructure  DistSparseMatrix<F>::GetLocalStructure() const {
    return Local_;
  }

  template <class F> SparseMatrixStructure DistSparseMatrix<F>::GetGlobalStructure(){
    if(!globalAllocated){
      Local_.ToGlobal(Global_);
      globalAllocated = true;
    }
    return Global_;
  }




  template <class F> void DistSparseMatrix<F>::GetLColRowCount(ETree & tree, IntNumVec & cc, IntNumVec & rc){
    TIMER_START(GetColRowCount);
    //The tree need to be postordered
    if(!tree.IsPostOrdered()){
      TIMER_START(PostOrder);
      tree.PostOrderTree();
      TIMER_START(PostOrder);
    }



    IntNumVec & rowind = this->Global_.rowind;
    IntNumVec & colptr = this->Global_.colptr;

    TIMER_START(Initialize_Data);
    //cc first contains the delta
    cc.Resize(size);
    //Compute size of subtrees
    IntNumVec treeSize(size);
    SetValue(treeSize,I_ONE);



    IntNumVec level(size);
    level(size-1)=1;
    for(Int vertex = 1; vertex<=size-1; vertex++){
      Int curParent = tree.PostParent(vertex-1);
      treeSize(curParent-1)+=treeSize(vertex-1);
      level(size-vertex-1) = level(tree.PostParent(size-vertex-1)-1)+1;
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

    DisjointSet sets;
    sets.Initialize(size);
    for(Int vertex = 1; vertex<=size; vertex++){
      Int cset = sets.makeSet (vertex);
      sets.Root(cset-1)=vertex;
    }


    TIMER_STOP(Initialize_Data);

    TIMER_START(Compute_Col_Row_Count);
    for(Int col=1; col<size; col++){
      Int cset;

      Int colPar = tree.PostParent(col-1);
      if (col<size){
        cc(colPar-1)--;
      }

      Int oCol = tree.FromPostOrder(col);
      for (Int i = colptr(oCol-1); i < colptr(oCol); i++) {
        Int row = tree.ToPostOrder(rowind(i-1));
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



    //    logfileptr->OFS()<<"Deltas "<<cc.m()<<std::endl;
    //    for(Int i = 0; i<cc.m();i++){
    //      logfileptr->OFS()<<cc(i)<<" ";
    //    }
    //    logfileptr->OFS()<<std::endl;






    //convert delta to col count
    for(Int col=1; col<size; col++){
      Int parent = tree.PostParent(col-1);
      cc(parent-1)+= cc(col-1);
    }



    //    logfileptr->OFS()<<"colcnt "<<cc.m()<<std::endl;
    //    for(Int i = 0; i<cc.m();i++){
    //      logfileptr->OFS()<<cc(i)<<" ";
    //    }
    //    logfileptr->OFS()<<std::endl;



    TIMER_STOP(GetColRowCount);
  }



  template <class F> void DistSparseMatrix<F>::FindSupernodes(ETree& tree, IntNumVec & cc, IntNumVec & xsuper){
    TIMER_START(FindSupernodes);


    IntNumVec children(size);
    SetValue(children,I_ZERO);
    for(Int col=1; col<=size-1; col++){
      children(tree.PostParent(col-1)-1)++;
      //      children(tree.PostParent(col)-1)++;
    }


    //    logfileptr->OFS()<<"children "<<children.m()<<std::endl;
    //    for(Int i = 0; i<children.m();i++){
    //      logfileptr->OFS()<<children(i)<<" ";
    //    }
    //    logfileptr->OFS()<<std::endl;

    Int nsuper = 1;
    xsuper.Resize(2*size);
    xsuper(nsuper-1) = 1;
    //SetValue(xsuper,I_ONE);
    for(Int i =2; i<=size;i++){
      //      logfileptr->OFS()<<"Column "<<i<<" has "<<children(i-1)<<" children, "<<cc(i-1)<<" nnz vs "<<cc(i-2)<<" in the prev col"<<std::endl;

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


  template <class F> void DistSparseMatrix<F>::SymbolicFactorization(ETree& tree,const IntNumVec & cc,const IntNumVec & xsuper, IntNumVec & xlindx, IntNumVec & lindx){
    TIMER_START(SymbolicFactorization);


    IntNumVec & rowind = this->Global_.rowind;
    IntNumVec & colptr = this->Global_.colptr;

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
      Int * stop = (end-1<rowind.m())?&rowind(end-1):&rowind(rowind.m()-1)+1; 
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

      //     logfileptr->OFS()<<"Final L"<<I<<": ";
      //     for(int i=0;i<LI.m();i++){logfileptr->OFS()<<LI(i)<< " ";}
      //     logfileptr->OFS()<<std::endl;

      lindxCnt += LI.m();

      if(length>width){
        //logfileptr->OFS()<<"I:"<<I<<std::endl<<"LI:"<<LI<<std::endl;
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
      std::copy(&LI(0),LI.Data()+LI.m(),&(lindx(head-1)));
      head+=LI.m();//cc(fi-1);
    }
    //    lindx.Resize(head-1);
    xlindx(nsuper) = head;

    TIMER_STOP(SymbolicFactorization);
  }



  template <class F> void DistSparseMatrix<F>::CopyData(const csc_matrix_t * cscptr){
    int np;
    int iam;

      MPI_Comm_size(comm,&np);
      MPI_Comm_rank(comm, &iam);

    //fill global structure info as we have it directly
    this->size = cscptr->n; 
    this->nnz = cscptr->nnz; 
    this->Global_.size = this->size;
    this->Global_.nnz = this->nnz;
    this->Global_.colptr.Resize(cscptr->n+1);
    std::copy(cscptr->colptr,cscptr->colptr+cscptr->n+1,this->Global_.colptr.Data());
    this->Global_.rowind.Resize(cscptr->nnz+1);
    std::copy(cscptr->rowidx,cscptr->rowidx+cscptr->nnz+1,this->Global_.rowind.Data());

    //move to 1 based indices
    for(int i=0;i<this->Global_.colptr.m();i++){ ++this->Global_.colptr[i]; }
    for(int i=0;i<this->Global_.rowind.m();i++){ ++this->Global_.rowind[i]; }

    this->globalAllocated = true;

    //Compute local structure info
    // Compute the number of columns on each processor
    IntNumVec numColLocalVec(np);
    Int numColLocal, numColFirst;
    numColFirst = this->size / np;
    SetValue( numColLocalVec, numColFirst );
    numColLocalVec[np-1] = this->size - numColFirst * (np-1) ;  // Modify the last entry	
    numColLocal = numColLocalVec[iam];

    this->Local_.colptr.Resize( numColLocal + 1 );

    for( Int i = 0; i < numColLocal+1; i++ ){
      this->Local_.colptr[i] = this->Global_.colptr[iam * numColFirst+i] - this->Global_.colptr[iam * numColFirst] + 1;
    }

    this->Local_.size = this->size;
    // Calculate nnz_loc on each processor
    this->Local_.nnz = this->Local_.colptr[numColLocal] - this->Local_.colptr[0];

    // Resize rowind and nzval appropriately 
    this->Local_.rowind.Resize( this->Local_.nnz );
    this->nzvalLocal.Resize ( this->Local_.nnz );

    //Read my row indices
    Int prevRead = 0;
    Int numRead = 0;
    for( Int ip = 0; ip <iam; ip++ ){	
      prevRead += this->Global_.colptr[ip*numColFirst + numColLocalVec[ip]]
        - this->Global_.colptr[ip*numColFirst];
    }

    numRead = this->Global_.colptr[iam*numColFirst + numColLocalVec[iam]] - this->Global_.colptr[iam*numColFirst];
    std::copy(&this->Global_.rowind[prevRead],&this->Global_.rowind[prevRead+numRead],this->Local_.rowind.Data());

    //copy appropriate nnz values
    if(cscptr->value_type == REAL){
      std::copy(&((const double*)cscptr->values)[prevRead],&((const double*)cscptr->values)[prevRead+numRead],this->nzvalLocal.Data());
    }
    else if(cscptr->value_type == COMPLEX){
      std::copy(&((const double*)cscptr->values)[prevRead],&((const double*)cscptr->values)[prevRead+numRead],this->nzvalLocal.Data());
    }
  }


  template <class F> DistSparseMatrix<F>::DistSparseMatrix(const csc_matrix_t * cscptr,MPI_Comm oComm ):comm(oComm){
    this->CopyData(cscptr);
  }



}


#endif // _DIST_SPARSE_MATRIX_IMPL_HPP_
