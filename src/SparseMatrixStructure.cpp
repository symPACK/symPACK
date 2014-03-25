#include "SparseMatrixStructure.hpp"
#include "ETree.hpp"
#include "utility.hpp"

namespace LIBCHOLESKY{
//  void SparseMatrixStructure::ClearExpandedSymmetric(){
//    expColptr.Resize(0);
//    expRowind.Resize(0);
//    expanded=false;
//  }
//
//  void SparseMatrixStructure::ExpandSymmetric(){
//    if(!expanded){
//      //code from sparsematrixconverter
//      /* set-up */
//
//      IntNumVec cur_col_nnz(size);
//      IntNumVec new_col_nnz(size);
//
//      /*
//       * Scan A and count how many new non-zeros we'll need to create.
//       *
//       * Post:
//       *   cur_col_nnz[i] == # of non-zeros in col i of the original symmetric 
//       *                     matrix.
//       *   new_col_nnz[i] == # of non-zeros to be stored in col i of the final 
//       *                     expanded matrix.
//       *   new_nnz == total # of non-zeros to be stored in the final expanded 
//       *              matrix.
//       */
//      Int new_nnz = 0; 
//      for (Int i = 0; i < size; i++) 
//      {    
//        cur_col_nnz[i] = colptr[i+1] - colptr[i];
//        new_col_nnz[i] = cur_col_nnz[i];
//        new_nnz += new_col_nnz[i];
//
//      }    
//
//
//
//
//
//      for (Int i = 0; i < size; i++) 
//      {    
//        Int k;
//        for (k = colptr[i]; k < colptr[i+1]; k++) 
//        {
//          Int j = rowind[k-1]-1;
//          if (j != i)
//          {
//            new_col_nnz[j]++;
//            new_nnz++;
//          }
//        }
//      }
//
//      /*
//       *  Initialize row pointers in expanded matrix.
//       *
//       *  Post:
//       *    new_colptr initialized to the correct, final values.
//       *    new_col_nnz[i] reset to be equal to cur_col_nnz[i].
//       */
//      expColptr.Resize(size+1);
//      expColptr[0] = 1;
//      for (Int i = 1; i <= size; i++)
//      {
//        expColptr[i] = expColptr[i-1] + new_col_nnz[i-1];
//        new_col_nnz[i-1] = cur_col_nnz[i-1];
//      }
//      expColptr[size] = new_nnz+1;
//
//      expRowind.Resize(new_nnz);
//
//      /*
//       *  Complete expansion of A to full storage.
//       *
//       *  Post:
//       *    (new_colptr, new_rowind, new_values) is the full-storage equivalent of A.
//       *    new_col_nnz[i] == # of non-zeros in col i of A.
//       */
//
//
//
//      for (Int i = 0; i < size; i++)
//      {
//        Int cur_nnz = cur_col_nnz[i];
//        Int k_cur   = colptr[i] -1;
//        Int k_new   = expColptr[i] -1;
//
//        /* copy current non-zeros from old matrix to new matrix */
//        std::copy(rowind.Data() + k_cur, rowind.Data() + k_cur + cur_nnz , expRowind.Data() + k_new);
//
//        /* fill in the symmetric "missing" values */
//        while (k_cur < colptr[i+1]-1)
//        {
//          /* non-zero of original matrix */
//          Int j = rowind[k_cur]-1;
//
//          if (j != i)  /* if not a non-diagonal element ... */
//          {
//            /* position of this transposed element in new matrix */
//            k_new = expColptr[j]-1 + new_col_nnz[j];
//
//            /* store */
//            expRowind[k_new] = i+1;
//            /*  update so next element stored at row j will appear
//             *  at the right place.  */
//            new_col_nnz[j]++;
//          }
//
//          k_cur++;
//        }
//      }
//      expanded =true;
//    }
//
//  }

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


  void SparseMatrixStructure::GetLColRowCount(ETree & tree, IntNumVec & cc, IntNumVec & rc){


    TIMER_START(GetColRowCount);
    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call GetLColRowCount\n" );
    }

    //The tree need to be postordered
    if(!tree.IsPostOrdered()){
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



  void SparseMatrixStructure::FindSupernodes(ETree& tree, IntNumVec & cc, IntNumVec & xsuper){
    TIMER_START(FindSupernodes);

    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call FindSupernodes\n" );
    }

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



void SparseMatrixStructure::GetSuperARowStruct(const ETree & etree, const IntNumVec & Xsuper, const Int iSupNo, std::vector<Int> & SuperRowStruct){
  TIMER_START(SpStruct_GetSuperARowStruct);
//  TIMER_START(SparseMatrixStructure::GetSuperARowStruct);
    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call GetARowStruct\n" );
    }

  Int first_col = Xsuper[iSupNo-1];
  Int last_col = Xsuper[iSupNo]-1;

  for(Int iPOCurCol = 1; iPOCurCol<first_col;++iPOCurCol){
    Int iCurCol = etree.FromPostOrder(iPOCurCol);

    Int iFirstRowPtr = colptr(iCurCol-1);
    Int iLastRowPtr = colptr(iCurCol)-1;

//    logfileptr->OFS()<<iFirstRowPtr<<" to "<<iLastRowPtr<<std::endl;
    for(Int iCurRowPtr=iFirstRowPtr ;iCurRowPtr<=iLastRowPtr;++iCurRowPtr){
      Int iPOCurRow = etree.ToPostOrder(rowind(iCurRowPtr-1));
      if(iPOCurRow >= first_col && iPOCurRow <= last_col){
        SuperRowStruct.push_back(iPOCurCol);
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


void SparseMatrixStructure::GetSuperLRowStruct(const ETree & etree, const IntNumVec & Xsuper, const Int iSupNo, std::set<Int> & SuperLRowStruct){

  TIMER_START(SpStruct_GetSuperLRowStruct);
    if(!bIsGlobal){
			throw std::logic_error( "SparseMatrixStructure must be global in order to call GetSuperLRowStruct\n" );
    }

//  SuperLRowStruct.clear();

    //Get A row struct
    std::vector<Int> SuperARowStruct;
    GetSuperARowStruct(etree, Xsuper, iSupNo, SuperARowStruct);

  Int first_col = Xsuper[iSupNo-1];
  Int last_col = Xsuper[iSupNo]-1;

  for(Int i = 0; i<SuperARowStruct.size();++i){
    Int iCurNode = SuperARowStruct[i];
    //tracing from iCurRow to iRow;
//    logfileptr->OFS()<<"Starting from node "<<iCurNode<<std::endl;
    //if(iCurNode==iPORow){
    if(iCurNode >= first_col && iCurNode <= last_col){
//      logfileptr->OFS()<<"Adding node "<<iCurNode<<std::endl;
      SuperLRowStruct.insert(iCurNode);
    }
    else{
      while( iCurNode != first_col && etree.PostParent(iCurNode-1) != 0){
//      while( !(iCurNode >= first_col && iCurNode <= last_col) && etree.PostParent(iCurNode-1) != 0){
//        logfileptr->OFS()<<"Adding node "<<iCurNode<<std::endl;
        SuperLRowStruct.insert(iCurNode);
        iCurNode = etree.PostParent(iCurNode-1);
//        logfileptr->OFS()<<"Parent is "<<iCurNode<<std::endl;
      }
    }
  } 
  TIMER_STOP(SpStruct_GetSuperLRowStruct);
}




}

