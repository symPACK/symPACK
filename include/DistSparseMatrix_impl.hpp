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

#ifdef UPCXX
#include <upcxx.h>
#endif
namespace LIBCHOLESKY{













  template <class F> void DistSparseMatrix<F>::ToGlobalStruct(){
    // test if structure hasn't been allocated yet
    if(globalAllocated){
      return;
    }


    int np;
    int iam;

    int ismpi=0;
    MPI_Initialized( &ismpi);
    int isnull= (comm == MPI_COMM_NULL);
    logfileptr->OFS()<<ismpi<<std::endl;

    logfileptr->OFS()<<comm<<std::endl;
    logfileptr->OFS()<<MPI_COMM_NULL<<std::endl;
    logfileptr->OFS()<<isnull<<std::endl;

    if(ismpi && isnull==0){
      MPI_Comm_size(comm,&np);
      MPI_Comm_rank(comm, &iam);
    }
    else{
      //throw an exception
			throw std::logic_error("MPI needs to be available.");
    }


    Global_.rowind.Resize(nnz);
    Global_.colptr.Resize(size+1);


    /* Allgatherv for row indices. */ 
    IntNumVec prevnz(np);
    IntNumVec rcounts(np);
    MPI_Allgather(&Local_.nnz, 1, MPI_INT, rcounts.Data(), 1, MPI_INT, comm);

    prevnz[0] = 0;
    for (Int i = 0; i < np-1; ++i) { prevnz[i+1] = prevnz[i] + rcounts[i]; } 

    MPI_Allgatherv(Local_.rowind.Data(), Local_.nnz, MPI_INT, Global_.rowind.Data(),rcounts.Data(), prevnz.Data(), MPI_INT, comm); 

    /* Allgatherv for colptr */
    // Compute the number of columns on each processor
    Int numColFirst = size / np;
    SetValue( rcounts, numColFirst );
    rcounts[np-1] = size - numColFirst * (np-1);  // Modify the last entry	


    IntNumVec rdispls(np);
    rdispls[0] = 0;
    for (Int i = 0; i < np-1; ++i) { rdispls[i+1] = rdispls[i] + rcounts[i]; } 


    MPI_Allgatherv(Local_.colptr.Data(), Local_.colptr.m()-1, MPI_INT, Global_.colptr.Data(),
        rcounts.Data(), rdispls.Data(), MPI_INT, comm);


    /* Recompute column pointers. */
    for (Int p = 1; p < np; p++) {
      Int idx = rdispls[p];
      for (Int j = 0; j < rcounts[p]; ++j) Global_.colptr[idx++] += prevnz[p];
    }

    Global_.colptr(size)=nnz+1;
    Global_.nnz = nnz;


    globalAllocated = true;



  }


  template <class F> void DistSparseMatrix<F>::ConstructETreeBis(ETree & tree){
    TIMER_START(ConstructETree);

    if(!globalAllocated){
      ToGlobalStruct();
    }

    Global_.ExpandSymmetric();


    tree.n_ = size;
    tree.parent_.Resize(size);
    SetValue(tree.parent_,I_ZERO );

    IntNumVec ancestors(size);
    SetValue(ancestors,I_ZERO );

    Int k;
    //for each row
    for (Int i = 1; i <= size; i++) {

#ifdef _DEBUG_
      logfileptr->OFS()<<"Examining col "<<i<<std::endl;
#endif
      for (Int p = Global_.expColptr(i-1); p < Global_.expColptr(i); p++) {
        k = Global_.expRowind(p-1);


        if (k >= i) continue;

        Int r = k;
        Int t;
        while(ancestors(r-1)!=0 && ancestors(r-1)!=i){
          t = ancestors(r-1);
          ancestors(r-1)=i;
          r=t;
        }

        if(ancestors(r-1)==0){
          ancestors(r-1)=i;
          tree.parent_(r-1)=i;

#ifdef _DEBUG_
          logfileptr->OFS()<<"Parent of "<<r<<" is "<<i<<std::endl;
#endif
        }
      }

    }


    tree.parent_(size-1) = 0;

    TIMER_STOP(ConstructETree);
  }


  template <class F> void DistSparseMatrix<F>::ConstructETree(ETree & tree){
    TIMER_START(ConstructETree);

    if(!globalAllocated){
      ToGlobalStruct();
    }

    Global_.ExpandSymmetric();

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

#ifdef _DEBUG_
      logfileptr->OFS()<<"Examining col "<<col<<std::endl;
#endif
      for (Int p = Global_.expColptr(col-1); p < Global_.expColptr(col); p++) {
        row = Global_.expRowind(p-1);

#ifdef _DEBUG_
        logfileptr->OFS()<<"Row = "<<row<<" vs col = "<<col<<std::endl;
#endif


        if (row >= col) continue;

        rset = sets.find(row);
        rroot = sets.Root(rset-1);
#ifdef _DEBUG_
        logfileptr->OFS()<<"Row "<<row<<" is in set "<<rset<<" represented by "<<rroot<<std::endl;
#endif

        if (rroot != col) {
          tree.parent_(rroot-1) = col;
          cset = sets.link(cset, rset);
          sets.Root(cset-1) = col;
#ifdef _DEBUG_
          logfileptr->OFS()<<"Parent of "<<rroot<<" is "<<col<<" which now represents set"<<cset<<std::endl;
#endif
        }
      }

    }


    tree.parent_(size-1) = 0;

    TIMER_STOP(ConstructETree);
  }






//  template <class F> void DistSparseMatrix<F>::ConstructETree(ETree & tree){
//    TIMER_START(ConstructETree);
//
//    if(!globalAllocated){
//      ToGlobalStruct();
//    }
//
//    IntNumVec rowind, colptr;
//    Global_.ExpandSymmetric(colptr,rowind);
//
//    tree.n_ = size;
//    tree.parent_.Resize(size);
//    SetValue(tree.parent_,I_ZERO );
//
//    ETree::DisjointSet sets;
//    sets.Initialize(size);
//
//    Int cset,rset,rroot,row;
//    for (Int col = 1; col <= size; col++) {
//      tree.parent_(col-1)=col; //1 based indexes
//      cset = sets.makeSet (col);
//      sets.Root(cset-1) = col;
//      tree.parent_(col-1) = size; 
//
//#ifdef _DEBUG_
//      logfileptr->OFS()<<"Examining col "<<col<<std::endl;
//#endif
//      for (Int p = Global_.colptr(col-1); p < Global_.colptr(col); p++) {
//        row = Global_.rowind(p-1);
//
//#ifdef _DEBUG_
//        logfileptr->OFS()<<"Row = "<<row<<" vs col = "<<col<<std::endl;
//#endif
//
//
//        if (row >= col) continue;
//
//        rset = sets.find(row);
//        rroot = sets.Root(rset-1);
//#ifdef _DEBUG_
//        logfileptr->OFS()<<"Row "<<row<<" is in set "<<rset<<" represented by "<<rroot<<std::endl;
//#endif
//
//        if (rroot != col) {
//          tree.parent_(rroot-1) = col;
//          cset = sets.link(cset, rset);
//          sets.Root(cset-1) = col;
//#ifdef _DEBUG_
//          logfileptr->OFS()<<"Parent of "<<rroot<<" is "<<col<<" which now represents set"<<cset<<std::endl;
//#endif
//        }
//      }
//
//    }
//
//
//    tree.parent_(size-1) = 0;
//
//    TIMER_STOP(ConstructETree);
//  }




  template <class F> void DistSparseMatrix<F>::GetLColRowCount(ETree & tree, IntNumVec & cc, IntNumVec & rc){
    TIMER_START(GetColRowCount);
    //The tree need to be postordered
    if(!tree.isPostOrdered_){
      TIMER_START(PostOrder);
      tree.PostOrderTree();
      TIMER_START(PostOrder);
    }



    IntNumVec & rowind = this->Global_.rowind;
    IntNumVec & colptr = this->Global_.colptr;

//    Global_.ExpandSymmetric();
//    IntNumVec & rowind = this->Global_.expRowind;
//    IntNumVec & colptr = this->Global_.expColptr;

    TIMER_START(Initialize_Data);
    //cc first contains the delta
    cc.Resize(size);
    //Compute size of subtrees
    IntNumVec treeSize(size);
    SetValue(treeSize,I_ONE);

logfileptr->OFS()<<"computing levels"<<std::endl;
   

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

logfileptr->OFS()<<"levels computed"<<std::endl;

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

logfileptr->OFS()<<"data initialized"<<std::endl;

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



    logfileptr->OFS()<<"Deltas "<<cc.m()<<std::endl;
    for(Int i = 0; i<cc.m();i++){
      logfileptr->OFS()<<cc(i)<<" ";
    }
    logfileptr->OFS()<<std::endl;






    //convert delta to col count
    for(Int col=1; col<size; col++){
      Int parent = tree.PostParent(col-1);
      cc(parent-1)+= cc(col-1);
    }



    logfileptr->OFS()<<"colcnt "<<cc.m()<<std::endl;
    for(Int i = 0; i<cc.m();i++){
      logfileptr->OFS()<<cc(i)<<" ";
    }
    logfileptr->OFS()<<std::endl;



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


  template <class F> void DistSparseMatrix<F>::SymbolicFactorization(ETree& tree,const IntNumVec & cc,const IntNumVec & xsuper){
    TIMER_START(SymbolicFactorization);





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

//     logfileptr->OFS()<<"L"<<I<<"<- A_*,fi: ";
     for(int i=0;i<LI.m();i++){logfileptr->OFS()<<LI(i)<< " ";}
     logfileptr->OFS()<<std::endl;

      LI = tree.ToPostOrder(LI);


//     logfileptr->OFS()<<"PO L"<<I<<"<- A_*,fi: ";
     for(int i=0;i<LI.m();i++){logfileptr->OFS()<<LI(i)<< " ";}
     logfileptr->OFS()<<std::endl;




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
//          Int head = LI.m();
//          for(Int i =1;i<LK.m();i++){
//            //insert element from LK
//            if(LK[i]>I ){
//              Ltmp[head]=LK[i];
//              head++;
//            }
//          }
//          Ltmp.Resize(head);
//          std::sort(&Ltmp(0),&Ltmp(Ltmp.m()-1)+1);
          
          //Int * end = std::set_union(&LI(0),&LI(LI.m()-1)+1,&LK(1),&LK(LK.m()-1)+1,&Ltmp(0));
          //Ltmp.Resize(end - &Ltmp(0));
          LI = Ltmp;
        }


      }
     
//     logfileptr->OFS()<<"Final L"<<I<<": ";
     for(int i=0;i<LI.m();i++){logfileptr->OFS()<<LI(i)<< " ";}
     logfileptr->OFS()<<std::endl;
      
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
    IntNumVec xlnz(size+1);
    Int totNnz = 1;
    for(Int I=1;I<xsuper.m();I++){
          Int fc = xsuper(I-1);
          Int lc = xsuper(I)-1;
        for(Int i=fc;i<=lc;i++){
          xlnz(i-1)=totNnz;
          totNnz+=cc(fc-1);
        }

    }
    xlnz(size)=totNnz;

    DblNumVec lnz(totNnz+1);

//    IntNumVec lindx(size*nsuper+1);
    IntNumVec lindx(lindxCnt);
    IntNumVec xlindx(nsuper+1);
    Int head = 1;
    for(Int I=1;I<=nsuper;I++){
      Int fi = tree.FromPostOrder(xsuper(I-1));
      IntNumVec & LI = LIs[I-1];
      xlindx(I-1)=head;


        logfileptr->OFS()<<"PO L"<<I<<":";
        for(int i=0;i<LI.m();i++){logfileptr->OFS()<<LI(i)<< " ";}
        logfileptr->OFS()<<std::endl;

      std::copy(&LI(0),&LI(LI.m()-1)+1,&(lindx(head-1)));
      head+=cc(fi-1);
    }
//    lindx.Resize(head-1);
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
      
//    for(Int I=1;I<xsuper.m();I++){
//          Int fc = xsuper(I-1);
//          Int lc = xsuper(I)-1;
//          Int fi = xlindx(I-1);
//
//          logfileptr->OFS()<<"FC = "<<fc<<std::endl;
//          logfileptr->OFS()<<"LC = "<<lc<<std::endl;
//
//          for(Int i = fc;i<=lc;i++){
//            Int fnz = xlnz(i-1);
//            Int lnz = xlnz(i)-1;
//
//          logfileptr->OFS()<<"i = "<<i<<" FNZ = "<<fnz<<std::endl;
//          logfileptr->OFS()<<"i = "<<i<<" LNZ = "<<lnz<<std::endl;
//            //diag is lnz(fnz)
//            logfileptr->OFS()<<"Diag element of "<<i<<" is "<<"lnz("<<fnz-1<<")"<<std::endl;
//            Int idx = fi;
//            for(Int s = fnz+1;s<=lnz;s++){
//              idx++;
//
//              //lnz(s) contains an off-diagonal nonzero entry in row lindx(idx)
//              logfileptr->OFS()<<"L("<<lindx(idx-1)<<","<<i<<") = "<<"lnz("<<s-1<<")"<<std::endl;
//            }
//            fi++;
//
//          }
//    }




    TIMER_STOP(SymbolicFactorization);
  }



  template <class F> void DistSparseMatrix<F>::CopyData(const csc_matrix_t * cscptr){
    int np;
    int iam;

    int ismpi=0;
    MPI_Initialized( &ismpi);
    int isnull= (comm == MPI_COMM_NULL);
    logfileptr->OFS()<<ismpi<<std::endl;

    logfileptr->OFS()<<comm<<std::endl;
    logfileptr->OFS()<<MPI_COMM_NULL<<std::endl;
    logfileptr->OFS()<<isnull<<std::endl;

    if(ismpi && isnull==0){
      MPI_Comm_size(comm,&np);
      MPI_Comm_rank(comm, &iam);
    }
    else{
#ifdef UPCXX
      np = THREADS;
      iam = MYTHREAD;
#else
      //throw an exception
			throw std::logic_error("Either MPI OR UPCXX need to be available.");
#endif
    }

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

//logfileptr->OFS()<<this->Global_.rowind<<endl;

	  // Compute the number of columns on each processor
	  IntNumVec numColLocalVec(np);
	  Int numColLocal, numColFirst;
	  numColFirst = this->size / np;
    SetValue( numColLocalVec, numColFirst );
    numColLocalVec[np-1] = this->size - numColFirst * (np-1) ;  // Modify the last entry	
//logfileptr->OFS()<<numColLocalVec<<endl;
  	numColLocal = numColLocalVec[iam];

//logfileptr->OFS()<<"NumColFirst = "<<iam*numColFirst<<endl;
//logfileptr->OFS()<<"NumColLocal = "<<numColLocal<<endl;

	  this->Local_.colptr.Resize( numColLocal + 1 );



  	for( Int i = 0; i < numColLocal+1; i++ ){
      
//      logfileptr->OFS()<<"i = "<<i<<endl;
//      logfileptr->OFS()<<"this->Global_.colptr["<<iam * numColFirst+i<<"] = "<<this->Global_.colptr[iam * numColFirst+i]<<endl;
//      logfileptr->OFS()<<"this->Global_.colptr["<<iam * numColFirst<<"] = "<<this->Global_.colptr[iam * numColFirst]<<endl;
	  	this->Local_.colptr[i] = this->Global_.colptr[iam * numColFirst+i] - this->Global_.colptr[iam * numColFirst] + 1;


	  }

    //logfileptr->OFS()<<"Local colptr = "<<this->Local_.colptr<<endl;
    //logfileptr->OFS()<<"Global colptr = "<<this->Global_.colptr<<endl;


	  // Calculate nnz_loc on each processor
	  this->Local_.nnz = this->Local_.colptr[numColLocal] - this->Local_.colptr[0];

    //logfileptr->OFS()<<"Local NNZ = "<<this->Local_.nnz<<endl;

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
    //logfileptr->OFS()<<"PrevRead = "<<prevRead<<endl;

		numRead = this->Global_.colptr[iam*numColFirst + numColLocalVec[iam]] - this->Global_.colptr[iam*numColFirst];
    //logfileptr->OFS()<<"NumRead = "<<numRead<<endl;
    std::copy(&this->Global_.rowind[prevRead],&this->Global_.rowind[prevRead+numRead],this->Local_.rowind.Data());

    //logfileptr->OFS()<<"RowindLocal"<<this->Local_.rowind<<endl;


    //copy appropriate nnz values
    if(cscptr->value_type == REAL){
      std::copy(&((const double*)cscptr->values)[prevRead],&((const double*)cscptr->values)[prevRead+numRead],this->nzvalLocal.Data());
    }
    else if(cscptr->value_type == COMPLEX){
      std::copy(&((const double*)cscptr->values)[prevRead],&((const double*)cscptr->values)[prevRead+numRead],this->nzvalLocal.Data());
    }
  

    //logfileptr->OFS()<<"nzvalLocal"<<this->nzvalLocal<<endl;


  }


  template <class F> DistSparseMatrix<F>::DistSparseMatrix(const csc_matrix_t * cscptr){
    this->CopyData(cscptr);
  }
  template <class F> DistSparseMatrix<F>::DistSparseMatrix(MPI_Comm oComm , const csc_matrix_t * cscptr):comm(oComm){
    this->CopyData(cscptr);
  }



}


#endif // _DIST_SPARSE_MATRIX_IMPL_HPP_
