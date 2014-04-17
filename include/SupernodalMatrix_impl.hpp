#ifndef _SUPERNODAL_MATRIX_IMPL_HPP_
#define _SUPERNODAL_MATRIX_IMPL_HPP_

#include "SupernodalMatrix.hpp"

#define TAG_FACTOR 0

namespace LIBCHOLESKY{

  inline void gdb_lock(){
    int lock = 1;
    while (lock == 1){
      lock =1;
    }
  }


  template <typename T> SupernodalMatrix<T>::SupernodalMatrix(){
  }

  template <typename T> SupernodalMatrix<T>::SupernodalMatrix(const DistSparseMatrix<T> & pMat,MAPCLASS & pMapping, MPI_Comm & pComm ){
    Int iam;
    MPI_Comm_rank(pComm, &iam);
    Int np;
    MPI_Comm_size(pComm, &np);

    iSize_ = pMat.size;
    Local_ = pMat.GetLocalStructure();
    Local_.ToGlobal(Global_);

    ETree_.ConstructETree(Global_);


    IntNumVec cc,rc;
    Global_.GetLColRowCount(ETree_,cc,rc);
    logfileptr->OFS()<<"colcnt "<<cc<<std::endl;
    logfileptr->OFS()<<"rowcnt "<<rc<<std::endl;

    ETree_.SortChildren(cc);

    logfileptr->OFS()<<"colcnt "<<cc<<std::endl;

    Global_.FindSupernodes(ETree_,cc,SupMembership_,Xsuper_);

    logfileptr->OFS()<<"Membership list is "<<SupMembership_<<std::endl;
    logfileptr->OFS()<<"xsuper "<<Xsuper_<<std::endl;

    //  UpdatesCount.Resize(Xsuper_.Size());
    //  for(Int I = 1; I<Xsuper_.Size();++I){
    //    Int first_col = Xsuper_[I-1];
    //    Int last_col = Xsuper_[I]-1;
    //    for(Int 
    //  }


    Global_.SymbolicFactorization2(ETree_,cc,Xsuper_,SupMembership_,xlindx_,lindx_);
    //  Global_.SymbolicFactorization(ETree_,cc,Xsuper_,xlindx_,lindx_);
    logfileptr->OFS()<<"xlindx "<<xlindx_<<std::endl;
    logfileptr->OFS()<<"lindx "<<lindx_<<std::endl;


    GetUpdatingSupernodeCount(UpdateCount_,UpdateWidth_);
    logfileptr->OFS()<<"Supernodal counts are "<<UpdateCount_<<std::endl;
    logfileptr->OFS()<<"Supernodal update width are "<<UpdateWidth_<<std::endl;


    SupETree_ = ETree_.ToSupernodalETree(Xsuper_);
    logfileptr->OFS()<<"Supernodal Etree is "<<SupETree_<<std::endl;

    Mapping_ = pMapping;

    //copy
    for(Int I=1;I<Xsuper_.m();I++){
      Int fc = Xsuper_(I-1);
      Int lc = Xsuper_(I)-1;
      Int iWidth = lc-fc+1;
      Int fi = xlindx_(I-1);
      Int li = xlindx_(I)-1;
      Int iHeight = li-fi+1;

      Int iDest = pMapping.Map(I-1,I-1);

      //parse the first column to create the NZBlock
      if(iam==iDest){
        LocalSupernodes_.push_back( new SuperNode<T>(I,fc,lc,iHeight));
        SuperNode<T> & snode = *LocalSupernodes_.back();


        for(Int idx = fi; idx<=li;idx++){
          Int iStartRow = lindx_(idx-1);
          Int iPrevRow = iStartRow;
          Int iContiguousRows = 1;
          for(Int idx2 = idx+1; idx2<=li;idx2++){
            Int iCurRow = lindx_(idx2-1);
            if(iStartRow == lindx_(fi-1)){
              if(iCurRow>iStartRow+iWidth-1){
                //enforce the first block to be a square diagonal block
                break;
              }
            }

            if(iCurRow==iPrevRow+1){
              idx++;
              ++iContiguousRows;
              iPrevRow=iCurRow;
            }
            else{
              break;
            }
          }

          Int iCurNZcnt = iContiguousRows * iWidth;
#ifdef _DEBUG_
          logfileptr->OFS()<<"Creating a new NZBlock for rows "<<iStartRow<<" to "<<iStartRow + iContiguousRows-1<<" with "<<iCurNZcnt<<" nz."<<std::endl;
#endif
          snode.AddNZBlock(iContiguousRows , iWidth,iStartRow);
          if(I==Xsuper_.m()-1){
            logfileptr->OFS()<<"Creating a new NZBlock for rows "<<iStartRow<<" to "<<iStartRow + iContiguousRows-1<<" with "<<iCurNZcnt<<" nz."<<std::endl;
            logfileptr->OFS()<<snode.GetNZBlock(snode.NZBlockCnt()-1)<<endl;
          }

        }

        snode.Shrink();
      }

      //Distribute the data

      //look at the owner of the first column
      Int numColFirst = iSize_ / np;
      Int prevOwner = -1;
      DblNumVec aRemoteCol;
      T * pdNzVal = NULL;
      Int iStartIdxCopy = 0;
      //copy the data from A into this Block structure
      for(Int i = fc;i<=lc;i++){

        Int iOwner = std::min((i-1)/numColFirst,np-1);

        if(iOwner != prevOwner){

          prevOwner = iOwner;
          //we need to transfer the bulk of columns
          //first compute the number of cols we need to transfer
          Int iNzTransfered = 0;
          Int iLCTransfered = lc;

          for(Int j =i;j<=lc;++j){
            iOwner = std::min((j-1)/numColFirst,np-1);
            if(iOwner == prevOwner){
              Int nrows = Global_.colptr(j) - Global_.colptr(j-1);
              iNzTransfered+=nrows;
              //              logfileptr->OFS()<<"Looking at col "<<j<<" which is on P"<<iOwner<<std::endl;
            }
            else{
              iLCTransfered = j-1;
              break;
            }
          } 


          //#ifdef _DEBUG_
          //logfileptr->OFS()<<"Column "<<i<<" to "<<iLCTransfered<<" are owned by P"<<prevOwner<<" and should go on P"<<iDest<<std::endl;
          //logfileptr->OFS()<<"They contain "<<iNzTransfered<<" nz"<<std::endl;
          //#endif

          //if data needs to be transfered
          if(iDest!=prevOwner){
            if(iam == iDest){
              aRemoteCol.Resize(iNzTransfered);

              //MPI_Recv
              MPI_Recv(aRemoteCol.Data(),iNzTransfered*sizeof(T),MPI_BYTE,prevOwner,0,pComm,MPI_STATUS_IGNORE);
              iStartIdxCopy = 0;
            } 
            else if (iam == prevOwner){
              Int local_i = (i-(numColFirst)*iam);
              Int iColptrLoc = Local_.colptr(local_i-1);
              Int iRowIndLoc = Local_.rowind(iColptrLoc-1);
              Int iLastRowIndLoc = Local_.rowind(Local_.colptr(local_i)-1-1);
              const T * pdData = &pMat.nzvalLocal(iColptrLoc-1);

              //MPI_send
              MPI_Send(pdData,iNzTransfered*sizeof(T),MPI_BYTE,iDest,0,pComm);
            }
          }
        }

        //copy the data if I own it
        if(iam == iDest){
          SuperNode<T> & snode = *LocalSupernodes_.back();

          //isData transfered or local
          if(iam!=prevOwner){
            pdNzVal = aRemoteCol.Data();
            //logfileptr->OFS()<<"pdNzVal is the remote buffer"<<std::endl;
            //logfileptr->OFS()<<aRemoteCol<<std::endl;
          }
          else{
            Int local_i = (i-(numColFirst)*iam);
            Int iColptrLoc = Local_.colptr(local_i-1);
            pdNzVal = &pMat.nzvalLocal(iColptrLoc-1);
            //logfileptr->OFS()<<"pdNzVal is the local pMat"<<std::endl;
            iStartIdxCopy = 0;
          }

          //Copy the data from pdNzVal in the appropriate NZBlock

          Int iGcolptr = Global_.colptr(i-1);
          Int iRowind = Global_.rowind(iGcolptr-1);
          Int iNrows = Global_.colptr(i) - Global_.colptr(i-1);
          Int idxA = 0;
          Int iLRow = 0;
          Int firstrow = fi + i-fc;
          for(Int idx = firstrow; idx<=li;idx++){
            iLRow = lindx_(idx-1);
            //logfileptr->OFS()<<"Looking at L("<<iLRow<<","<<i<<")"<<std::endl;
            if( iLRow == iRowind){
              Int iNZBlockIdx = snode.FindBlockIdx(iLRow);

              T elem = pdNzVal[iStartIdxCopy + idxA];

              NZBlock<T> & pDestBlock = snode.GetNZBlock(iNZBlockIdx);
              //logfileptr->OFS()<<*pDestBlock<<std::endl;

              //find where we should put it
              Int localCol = i - fc;
              Int localRow = iLRow - pDestBlock.GIndex();
#ifdef _DEBUG_
              logfileptr->OFS()<<"Elem is A("<<iRowind<<","<<i<<") = "<<elem<<" at ("<<localRow<<","<<localCol<<") in "<<iNZBlockIdx<<"th NZBlock of L"<< std::endl;
#endif
              pDestBlock.Nzval(localRow,localCol) = elem;
              if(idxA<iNrows){
                idxA++;
                iRowind = Global_.rowind(iGcolptr+idxA-1-1);
              }
            }
          }

        }

      }
#ifdef _DEBUG_
      logfileptr->OFS()<<"--------------------------------------------------"<<std::endl;
#endif
    }

  }


  template <typename T> void SupernodalMatrix<T>::GetUpdatingSupernodeCount(IntNumVec & sc,IntNumVec & mw){
    sc.Resize(Xsuper_.m());
    SetValue(sc,I_ZERO);
    IntNumVec marker(Xsuper_.m());
    SetValue(marker,I_ZERO);
    mw.Resize(Xsuper_.m());
    SetValue(mw,I_ZERO);

    for(Int s = 1; s<Xsuper_.m(); ++s){
      Int first_col = Xsuper_(s-1);
      Int last_col = Xsuper_(s)-1;

      Int fi = xlindx_(s-1);
      Int li = xlindx_(s)-1;

      logfileptr->OFS()<<"Supernode "<<s<<" updates: ";

      for(Int row_idx = fi; row_idx<=li;++row_idx){
        Int row = lindx_(row_idx-1);
        Int supno = SupMembership_(row-1);

        if(marker(supno-1)!=s && supno!=s){

          logfileptr->OFS()<<supno<<" ";
          ++sc(supno-1);
          marker(supno-1) = s;

          mw(supno-1) = max(mw(supno-1),last_col - first_col+1);

        }
      }
      logfileptr->OFS()<<std::endl;
    }
  }


  template <typename T> SparseMatrixStructure SupernodalMatrix<T>::GetLocalStructure() const {
    return Local_;
  }

  template <typename T> SparseMatrixStructure SupernodalMatrix<T>::GetGlobalStructure(){
    if(!globalAllocated){
      Local_.ToGlobal(Global_);
      globalAllocated= true;
    }
    return Global_;
  }

  //FIXME write this function also in terms of SparseMatrixStructure and supno rather than Supernode object.
  template <typename T> inline bool SupernodalMatrix<T>::FindNextUpdate(Int src_snode_id, Int & src_first_row, Int & src_last_row, Int & tgt_snode_id){
    Int src_fc = Xsuper_(src_snode_id-1);
    Int src_lc = Xsuper_(src_snode_id)-1;

    //look in the global structure for the next nnz below row src_lc  
    Int first_row_ptr = xlindx_(src_snode_id-1);
    Int last_row_ptr = xlindx_(src_snode_id)-1;
    Int src_last_row_ptr = 0;
    Int src_first_row_ptr = 0;

    Int subdiag_row_cnt = last_row_ptr - first_row_ptr;

    Int first_row = lindx_(first_row_ptr-1);
    Int last_row = lindx_(last_row_ptr-1);

    //    logfileptr->OFS()<<"prev src_first_row = "<<src_first_row<<endl;
    //    logfileptr->OFS()<<"prev src_last_row = "<<src_last_row<<endl;
    //    logfileptr->OFS()<<"first_row = "<<first_row<<endl;
    //    logfileptr->OFS()<<"last_row = "<<last_row<<endl;
    //    logfileptr->OFS()<<"first_row_ptr = "<<first_row_ptr<<endl;
    //    logfileptr->OFS()<<"last_row_ptr = "<<last_row_ptr<<endl;
    //    logfileptr->OFS()<<"src_last_row_ptr = "<<src_last_row_ptr<<endl;
    //    logfileptr->OFS()<<"subdiag_row_cnt = "<<subdiag_row_cnt<<endl;


    //if tgt_snode_id == 0 , this is the first call to the function
    if(tgt_snode_id == 0){

      if(subdiag_row_cnt == 0 ){
        return false;
      }

      Int first_row = lindx_(first_row_ptr);
      tgt_snode_id = SupMembership_(first_row-1);
      src_first_row = first_row;
      src_last_row_ptr = first_row_ptr;
    }
    else{

      //find the block corresponding to src_last_row
      //src_nzblk_idx = src_snode.FindBlockIdx(src_last_row);

      //src_last_row_ptr = src_last_row - lindx_(first_row_ptr-1) + first_row_ptr;
      src_last_row_ptr = std::find(&lindx_(first_row_ptr-1),&lindx_(last_row_ptr), src_last_row) - &lindx_(first_row_ptr-1) + first_row_ptr;
      if(src_last_row_ptr == last_row_ptr){
        return false;
      }
      else{
        src_first_row_ptr = src_last_row_ptr+1;
        if(src_first_row_ptr > last_row_ptr){
          return false;
        }
        else{
          Int first_row = lindx_(src_first_row_ptr-1);
          tgt_snode_id = SupMembership_(first_row-1);
          src_first_row = first_row;
        }
      }
    }

    //Now we try to find src_last_row
    src_first_row_ptr = src_last_row_ptr+1;
    Int src_fr = src_first_row;

    //Find the last contiguous row
    Int src_lr = src_first_row;
    for(Int i = src_first_row_ptr+1; i<=last_row_ptr; ++i){
      if(src_lr+1 == lindx_(i-1)){ ++src_lr; }
      else{ break; }
    }

    Int tgt_snode_id_first = SupMembership_(src_fr-1);
    Int tgt_snode_id_last = SupMembership_(src_lr-1);
    if(tgt_snode_id_first == tgt_snode_id_last){
      //this can be a zero row in the src_snode
      tgt_snode_id = tgt_snode_id_first;

      //    src_last_row = min(src_lr,Xsuper_(tgt_snode_id_first)-1);

      //Find the last row in src_snode updating tgt_snode_id
      for(Int i = src_first_row_ptr +1; i<=last_row_ptr; ++i){
        Int row = lindx_(i-1);
        if(SupMembership_(row-1) != tgt_snode_id){ break; }
        else{ src_last_row = row; }
      }

    }
    else{
      src_last_row = Xsuper_(tgt_snode_id_first)-1;
      tgt_snode_id = tgt_snode_id_first;
    }

    return true;

  }


  template <typename T> inline bool SupernodalMatrix<T>::FindNextUpdate(SuperNode<T> & src_snode, Int & src_nzblk_idx, Int & src_first_row, Int & src_last_row, Int & tgt_snode_id){
    //src_nzblk_idx is the last nzblock index examined
    //src_first_row is the first row updating the supernode examined
    //src_last_row is the last row updating the supernode examined
    //if(src_snode.Id() == 25){gdb_lock();}

    //if tgt_snode_id == 0 , this is the first call to the function
    if(tgt_snode_id == 0){

      //find the first sub diagonal block
      src_nzblk_idx = -1;
      for(Int blkidx = 0; blkidx < src_snode.NZBlockCnt(); ++blkidx){
        if(src_snode.GetNZBlock(blkidx).GIndex() > src_snode.LastCol()){
          src_nzblk_idx = blkidx;
          break;
        }
      }

      if(src_nzblk_idx == -1 ){
        return false;
        //logfileptr->OFS()<<std::endl;
        //logfileptr->OFS()<<std::endl;
      }
      else{
        src_first_row = src_snode.GetNZBlock(src_nzblk_idx).GIndex();
      }
    }
    else{
      //find the block corresponding to src_last_row
      src_nzblk_idx = src_snode.FindBlockIdx(src_last_row);
      if(src_nzblk_idx==src_snode.NZBlockCnt()){
        return false;
        //logfileptr->OFS()<<std::endl;
        //logfileptr->OFS()<<std::endl;
      }
      else{
        NZBlock<T> & nzblk = src_snode.GetNZBlock(src_nzblk_idx);
        Int src_lr = nzblk.GIndex()+nzblk.NRows()-1;
        if(src_last_row == src_lr){
          src_nzblk_idx++;
          if(src_nzblk_idx==src_snode.NZBlockCnt()){
            return false;


            //logfileptr->OFS()<<std::endl;
            //logfileptr->OFS()<<std::endl;

          }
          else{
            src_first_row = src_snode.GetNZBlock(src_nzblk_idx).GIndex();
          }
        }
        else{
          src_first_row = src_last_row+1;
        }
      }
    }

    //Now we try to find src_last_row
    NZBlock<T> & nzblk = src_snode.GetNZBlock(src_nzblk_idx);

    assert(src_first_row >= nzblk.GIndex());

    Int src_fr = max(src_first_row,nzblk.GIndex());
    Int src_lr = nzblk.GIndex()+nzblk.NRows()-1;
    Int tgt_snode_id_first = SupMembership_(src_fr-1);
    Int tgt_snode_id_last = SupMembership_(src_lr-1);

    //logfileptr->OFS()<<"First row of Supernode "<< src_snode.Id() <<" updates Supernode "<<tgt_snode_id_first<<std::endl;
    //logfileptr->OFS()<<"NZ block is "<<nzblk<<std::endl;

    if(tgt_snode_id_first == tgt_snode_id_last){
      //this can be a zero row in the src_snode
      tgt_snode_id = tgt_snode_id_first;

      src_last_row = Xsuper_(tgt_snode_id_first)-1;

      //look into other nzblk

      //Find the last row in src_snode updating tgt_snode_id
      Int new_blk_idx = src_snode.FindBlockIdx(src_last_row);
      if(new_blk_idx>=src_snode.NZBlockCnt()){
        //src_last_row is the last row of the last nzblock
        new_blk_idx = src_snode.NZBlockCnt()-1; 
      }

      NZBlock<T> & last_block = src_snode.GetNZBlock(new_blk_idx); 
      src_last_row = min(src_last_row,last_block.GIndex() + last_block.NRows() - 1);
      assert(src_last_row<= Xsuper_(tgt_snode_id_first)-1);
    }
    else{
      src_last_row = Xsuper_(tgt_snode_id_first)-1;
      tgt_snode_id = tgt_snode_id_first;
    }


    assert(src_last_row>=src_first_row);
    assert(src_first_row>= Xsuper_(tgt_snode_id-1));
    assert(src_last_row<= Xsuper_(tgt_snode_id)-1);
    assert(src_first_row >= src_fr);

    return true;

  }





  template <typename T> inline bool SupernodalMatrix<T>::FindPivot(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode,Int & pivot_idx, Int & pivot_fr, Int & pivot_lr){
    //Determine which columns will be updated
    pivot_idx = 0;
    pivot_fr = 0;
    pivot_lr = 0;
    for(int blkidx=0;blkidx<src_snode.NZBlockCnt();++blkidx){
      NZBlock<T> & nzblk = src_snode.GetNZBlock(blkidx);
      Int src_fr = nzblk.GIndex();
      Int src_lr = nzblk.GIndex()+nzblk.NRows()-1;

      if(src_lr>=tgt_snode.FirstCol() && src_fr<=tgt_snode.FirstCol() ){
        //found the "pivot" nzblk

        pivot_fr = max(tgt_snode.FirstCol(),src_fr);
        pivot_lr = min(tgt_snode.LastCol(),src_lr);

        pivot_idx = blkidx;
        return true;
      }
    }
  }

  template <typename T> void SupernodalMatrix<T>::UpdateSuperNode(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode, Int &pivot_idx, Int  pivot_fr = I_ZERO){
    NZBlock<T> & pivot_nzblk = src_snode.GetNZBlock(pivot_idx);
    if(pivot_fr ==I_ZERO){
      pivot_fr = pivot_nzblk.GIndex();
    }
    assert(pivot_fr >= pivot_nzblk.GIndex());
    Int pivot_lr = min(pivot_nzblk.GIndex() +pivot_nzblk.NRows() -1, pivot_fr + tgt_snode.Size() -1);
    T * pivot = & pivot_nzblk.Nzval(pivot_fr-pivot_nzblk.GIndex(),0);

    Int tgt_updated_fc =  pivot_fr - tgt_snode.FirstCol();

    for(int src_idx=pivot_idx;src_idx<src_snode.NZBlockCnt();++src_idx){

      NZBlock<T> & src_nzblk = src_snode.GetNZBlock(src_idx);
      Int src_fr = max(pivot_fr, src_nzblk.GIndex()) ;
      Int src_lr = src_nzblk.GIndex()+src_nzblk.NRows()-1;

      do{
        Int tgt_idx = tgt_snode.FindBlockIdx(src_fr);
        NZBlock<T> & tgt_nzblk = tgt_snode.GetNZBlock(tgt_idx);
        Int tgt_fr = tgt_nzblk.GIndex();
        Int tgt_lr = tgt_nzblk.GIndex()+tgt_nzblk.NRows()-1;

        Int update_fr = max(tgt_fr, src_fr);
        Int update_lr = min(tgt_lr, src_lr);

        assert(update_fr >= tgt_fr); 
        assert(update_lr >= update_fr); 
        assert(update_lr <= tgt_lr); 

        //Update tgt_nzblk with src_nzblk
        //        logfileptr->OFS()<<"Updating SNODE "<<tgt_snode.Id()<<" Block "<<tgt_idx<<" ["<<update_fr<<".."<<update_lr<<"] with SNODE "<<src_snode.Id()<<" Block "<<src_idx<<" ["<<update_fr<<".."<<update_lr<<"]"<<endl;

        T * src = &src_nzblk.Nzval(update_fr - src_fr ,0);
        T * tgt = &tgt_nzblk.Nzval(update_fr - tgt_fr ,tgt_updated_fc);

#ifdef ROW_MAJOR
        blas::Gemm('T','N',pivot_lr-pivot_fr+1, update_lr - update_fr + 1,src_snode.Size(),1.0,pivot,pivot_nzblk.LDA(),src,src_nzblk.LDA(),1.0,tgt,tgt_nzblk.LDA());
#else
        blas::Gemm('N','T', update_lr - update_fr + 1,pivot_lr-pivot_fr+1,  src_snode.Size(),1.0,src,src_nzblk.LDA(),pivot,pivot_nzblk.LDA(),1.0,tgt,tgt_nzblk.LDA());
#endif

        if(tgt_idx+1<tgt_snode.NZBlockCnt()){
          NZBlock<T> & next_tgt_nzblk = tgt_snode.GetNZBlock(tgt_idx+1);
          Int next_tgt_fr = next_tgt_nzblk.GIndex();
          src_fr = next_tgt_fr;
        }
        else{
          break;
        }
      }while(src_fr<src_lr);
    }

  }




  template <typename T> void SupernodalMatrix<T>::Factorize( MPI_Comm & pComm ){
    Int iam,np;
    MPI_Comm_rank(pComm, &iam);
    MPI_Comm_size(pComm, &np);
    IntNumVec UpdatesToDo = UpdateCount_;
    std::vector<std::stack<Int> > LocalUpdates(LocalSupernodes_.size());

    //dummy right looking cholesky factorization
    for(Int I=1;I<Xsuper_.m();I++){
      Int src_first_col = Xsuper_(I-1);
      Int src_last_col = Xsuper_(I)-1;
      Int iOwner = Mapping_.Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam ){



        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> & src_snode = *LocalSupernodes_[iLocalI -1];
        logfileptr->OFS()<<"Supernode "<<I<<"("<<src_snode.Id()<<") is on P"<<iOwner<<" local index is "<<iLocalI<<std::endl; 

        assert(src_snode.Id() == I); 

        //      std::set<Int> SuperLRowStruct;
        //      Global_.GetSuperLRowStruct(ETree_, Xsuper_, I, SuperLRowStruct);
        //        
        //      logfileptr->OFS()<<"Row structure of Supernode "<<I<<" is ";
        //      for(std::set<Int>::iterator it = SuperLRowStruct.begin(); it != SuperLRowStruct.end(); ++it){
        //        logfileptr->OFS()<<*it<<" ";
        //      }
        //      logfileptr->OFS()<<endl;
        //      logfileptr->OFS()<<"Updates count for Supernode "<<I<<" is "<<UpdateCount_(I-1)<<endl;
        //
        //      assert(SuperLRowStruct.size()==UpdateCount_(I-1));


        //Build a "dense" representation of the factor
        NumMat<T> dense_snode = src_snode;
        Int nrhs = 5;
        NumMat<T> RHS(dense_snode.n(),nrhs);
        NumMat<T> XTrue(dense_snode.n(),nrhs);
        UniformRandom(XTrue);
        blas::Gemm('N','N',dense_snode.n(),nrhs,dense_snode.n(),1.0,&dense_snode(0,0),dense_snode.m(),&XTrue(0,0),XTrue.m(),0.0,&RHS(0,0),RHS.m());


        //Do all my updates (Local and remote)
        //Local updates
        while(!LocalUpdates[iLocalI-1].empty()){
          Int src_snode_id = LocalUpdates[iLocalI-1].top();
          LocalUpdates[iLocalI-1].pop();
          Int src_nzblk_idx = LocalUpdates[iLocalI-1].top();
          LocalUpdates[iLocalI-1].pop();
          Int src_first_row = LocalUpdates[iLocalI-1].top();
          LocalUpdates[iLocalI-1].pop();

          SuperNode<T> & local_src_snode = *LocalSupernodes_[(src_snode_id-1) / np];

          logfileptr->OFS()<<"LOCAL Supernode "<<I<<" is updated by Supernode "<<src_snode_id<<" from row "<<src_first_row<<" "<<src_nzblk_idx<<std::endl;

          UpdateSuperNode(local_src_snode,src_snode,src_nzblk_idx, src_first_row);

          --UpdatesToDo(I-1);
          logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left"<<endl;
        }


        //Remote updates
        std::vector<char> src_nzblocks;
        size_t num_bytes;
        while(UpdatesToDo(I-1)>0){
          logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left"<<endl;

          if(src_nzblocks.size()==0){
            //The upper bound must be of the width of the "largest" child
            logfileptr->OFS()<<"Maximum width is "<<UpdateWidth_(I-1)<<std::endl;

            num_bytes = sizeof(Int)+sizeof(NZBlock<T>);
            for(Int blkidx=0;blkidx<src_snode.NZBlockCnt();++blkidx){
              num_bytes += src_snode.GetNZBlock(blkidx).NRows()*NZBLOCK_ROW_SIZE<T>(UpdateWidth_(I-1));
            }
            logfileptr->OFS()<<"We allocate a buffer of size "<<num_bytes<<std::endl;
            src_nzblocks.resize(num_bytes);
          }


          //MPI_Recv
          MPI_Status recv_status;

          MPI_Probe(MPI_ANY_SOURCE,I,pComm,&recv_status);
          int bytes_to_receive = 0;
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_to_receive);
          //          if(src_nzblocks.size()< bytes_to_receive + sizeof(Int)+sizeof(NZBlock<T>)){
          //            num_bytes = bytes_to_receive+sizeof(Int)+sizeof(NZBlock<T>);
          //            src_nzblocks.resize(num_bytes);
          //          }

          assert(src_nzblocks.size()>=bytes_to_receive);




          char * recv_buf = &src_nzblocks[0]+ max(sizeof(Int),NZBLOCK_OBJ_SIZE<T>())-min(sizeof(Int),NZBLOCK_OBJ_SIZE<T>());
          MPI_Recv(recv_buf,&*src_nzblocks.end()-recv_buf,MPI_BYTE,MPI_ANY_SOURCE,I,pComm,&recv_status);
          //logfileptr->OFS()<<"Received something"<<endl;
          Int src_snode_id = *(Int*)recv_buf;


          //Resize the buffer to the actual number of bytes received
          int bytes_received = 0;
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
          src_nzblocks.resize(recv_buf+bytes_received - &src_nzblocks[0]);


          //Create the dummy supernode for that data
          SuperNode<T> dist_src_snode(src_snode_id,Xsuper_[src_snode_id-1],Xsuper_[src_snode_id]-1,&src_nzblocks);

          logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode.Id()<<std::endl;
          //      logfileptr->OFS()<<dist_src_snode<<endl;


          //Update everything I own with that factor
          //Update the ancestors
          Int tgt_snode_id = 0;
          Int src_first_row = 0;
          Int src_last_row = 0;
          Int src_nzblk_idx = 0;
          while(FindNextUpdate(dist_src_snode, src_nzblk_idx, src_first_row, src_last_row, tgt_snode_id)){ 
            Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
            if(iTarget == iam){
              logfileptr->OFS()<<"RECV Supernode "<<tgt_snode_id<<" is updated by Supernode "<<dist_src_snode.Id()<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;

              Int iLocalJ = (tgt_snode_id-1) / np +1 ;
              SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];

              UpdateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, src_first_row);

              --UpdatesToDo(tgt_snode_id-1);
              logfileptr->OFS()<<UpdatesToDo(tgt_snode_id-1)<<" updates left for Supernode "<<tgt_snode_id<<endl;
            }
          }



          //restore to its capacity
          src_nzblocks.resize(num_bytes);

        }
        //clear the buffer
        { vector<char>().swap(src_nzblocks);  }

        logfileptr->OFS()<<"  Factoring Supernode "<<I<<std::endl;

        //Factorize Diagonal block
        NZBlock<T> & diagonalBlock = src_snode.GetNZBlock(0);

#ifdef ROW_MAJOR
        lapack::Potrf( 'U', src_snode.Size(), diagonalBlock.Nzval(), diagonalBlock.LDA());
#else
        lapack::Potrf( 'L', src_snode.Size(), diagonalBlock.Nzval(), diagonalBlock.LDA());
#endif

        //      logfileptr->OFS()<<"    Diagonal block factored node "<<I<<std::endl;

        for(int blkidx=1;blkidx<src_snode.NZBlockCnt();++blkidx){
          NZBlock<T> & nzblk = src_snode.GetNZBlock(blkidx);

          //Update lower triangular blocks
#ifdef ROW_MAJOR
          blas::Trsm('R','U','T','N',nzblk.NCols(),src_snode.Size(), 1.0,  diagonalBlock.Nzval(), diagonalBlock.LDA(), nzblk.Nzval(), nzblk.LDA());
#else
          blas::Trsm('R','L','T','N',nzblk.NRows(),src_snode.Size(), 1.0,  diagonalBlock.Nzval(), diagonalBlock.LDA(), nzblk.Nzval(), nzblk.LDA());
#endif
          //          logfileptr->OFS()<<diagonalBlock<<std::endl;
          //          logfileptr->OFS()<<nzblk<<std::endl;
          //          logfileptr->OFS()<<"    "<<blkidx<<"th subdiagonal block updated node "<<I<<std::endl;
        }


        //Get a dense representation of the factor
        NumMat<T> dense_factor = src_snode;

        logfileptr->OFS()<<"Testing Supernode "<<I<<endl;
        //check the result
        double norm = 0;
        //do a solve
        NumMat<T> X = RHS;
        Int n = dense_factor.n();
        lapack::Potrs('L',n,nrhs,&dense_factor(0,0),dense_factor.m(),&X(0,0),X.m());
        blas::Axpy(n*nrhs,-1.0,&XTrue(0,0),1,&X(0,0),1);
        norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
        logfileptr->OFS()<<"Norm of residual for Supernode "<<I<<" is "<<norm<<std::endl;

        //      NumMat<T> B(Size(),nrhs);
        //      Trsm( ONE<T>(), src_snode, B, pComm );



        //Send my factor to my ancestors. 
        BolNumVec is_factor_sent(np);
        SetValue(is_factor_sent,false);

        Int tgt_snode_id = 0;
        Int src_nzblk_idx = 0;
        Int src_first_row = 0;
        Int src_last_row = 0;
        while(FindNextUpdate(src_snode, src_nzblk_idx, src_first_row, src_last_row, tgt_snode_id)){ 

          Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
          if(iTarget != iam && !is_factor_sent[iTarget]){
            is_factor_sent[iTarget] = true;

            logfileptr->OFS()<<"Supernode "<<tgt_snode_id<<" is updated by Supernode "<<I<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;

            Int J = tgt_snode_id;

            Int tgt_first_col = Xsuper_(J-1);
            Int tgt_last_col = Xsuper_(J)-1;

            //Send
            NZBlock<T> * pivot_nzblk = &src_snode.GetNZBlock(src_nzblk_idx);

#ifdef ROW_MAJOR
            int local_first_row = src_first_row-pivot_nzblk->GIndex();

            char * start_buf = reinterpret_cast<char*>(&pivot_nzblk->Nzval(local_first_row,0));
            size_t num_bytes = src_snode.End() - start_buf;

            //Create a header
            NZBlockHeader<T> * header = new NZBlockHeader<T>(pivot_nzblk->NRows()-local_first_row,pivot_nzblk->NCols(),src_first_row);

            MPI_Datatype type;
            int lens[3];
            MPI_Aint disps[3];
            MPI_Datatype types[3];
            Int err = 0;


            /* define a struct that describes all our data */
            lens[0] = sizeof(Int);
            lens[1] = sizeof(NZBlockHeader<T>);
            lens[2] = num_bytes;
            MPI_Address(&I, &disps[0]);
            MPI_Address(header, &disps[1]);
            MPI_Address(start_buf, &disps[2]);
            types[0] = MPI_BYTE;
            types[1] = MPI_BYTE;
            types[2] = MPI_BYTE;
            MPI_Type_struct(3, lens, disps, types, &type);
            MPI_Type_commit(&type);

            //logfileptr->OFS()<<src_snode<<endl;

            /* send to target */
            //FIXME : maybe this can be replaced by a scatterv ?
            MPI_Send(MPI_BOTTOM,1,type,iTarget,tgt_snode_id,pComm);

            MPI_Type_free(&type);

            delete header;
            logfileptr->OFS()<<"     Send factor "<<I<<" to node"<<J<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<" ("<<num_bytes+sizeof(*header)+sizeof(Int)<<" bytes)"<<std::endl;
#else
#ifdef SEND_WHOLE_BLOCK
#else
            int local_first_row = src_first_row-pivot_nzblk->GIndex();

            char * start_buf = reinterpret_cast<char*>(src_snode.NZBlockCnt()>src_nzblk_idx+1?&src_snode.GetNZBlock(src_nzblk_idx+1):src_snode.End());
            size_t num_bytes = src_snode.End() - start_buf;

            //Create a NZBlock
            Int nrows = pivot_nzblk->NRows()-local_first_row;
            Int ncols = pivot_nzblk->NCols();
            std::vector<char> tmpStorage(NZBLOCK_HEADER_SIZE<T>() + nrows*ncols*sizeof(T));
            NZBlock<T> * newBlock = new NZBlock<T>(nrows,ncols,src_first_row, &tmpStorage[0]);
            //Fill the block
            for(Int col = 0; col< pivot_nzblk->NCols(); ++col){
              for(Int row = local_first_row; row< pivot_nzblk->NRows(); ++row){
                newBlock->Nzval(row-local_first_row,col) = pivot_nzblk->Nzval(local_first_row,col);
              }
            }

            //          logfileptr->OFS()<<"New block is "<<newBlock<<endl;


            MPI_Datatype type;
            int lens[3];
            MPI_Aint disps[3];
            MPI_Datatype types[3];
            Int err = 0;


            /* define a struct that describes all our data */
            lens[0] = sizeof(Int);
            //don't send the object itself because it will be created again on the receiver side
            lens[1] = NZBLOCK_HEADER_SIZE<T>()+newBlock->ByteSize(); 
            lens[2] = num_bytes;
            MPI_Address(&I, &disps[0]);
            MPI_Address(&tmpStorage[0], &disps[1]);
            MPI_Address(start_buf, &disps[2]);
            types[0] = MPI_BYTE;
            types[1] = MPI_BYTE;
            types[2] = MPI_BYTE;
            MPI_Type_struct(3, lens, disps, types, &type);
            MPI_Type_commit(&type);

            MPI_Send(MPI_BOTTOM,1,type,iTarget,tgt_snode_id,pComm);

            delete newBlock;
            logfileptr->OFS()<<"     SEND Supernode "<<I<<" to Supernode "<<J<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<" ("<<lens[0] + lens[1] + lens[2]<<" bytes)"<<std::endl;
#endif // end ifdef SEND_WHOLE_BLOCK

#endif // end ifdef ROW_MAJOR


          }
        }

        //Update my ancestors right away. 
        tgt_snode_id = 0;
        while(FindNextUpdate(src_snode, src_nzblk_idx, src_first_row, src_last_row, tgt_snode_id)){ 

          Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
          if(iTarget == iam){
            logfileptr->OFS()<<"LOCAL Supernode "<<tgt_snode_id<<" is updated by Supernode "<<I<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;

            Int iLocalJ = (tgt_snode_id-1) / np +1 ;
            LocalUpdates[iLocalJ-1].push(src_first_row);
            LocalUpdates[iLocalJ-1].push(src_nzblk_idx);
            LocalUpdates[iLocalJ-1].push(I);

            //          SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];
            //
            //          UpdateSuperNode(src_snode,tgt_snode,src_nzblk_idx, src_first_row);
            //
            //          --UpdatesToDo(tgt_snode_id-1);
            //          logfileptr->OFS()<<UpdatesToDo(tgt_snode_id-1)<<" updates left"<<endl;
          }
        }
      }
      MPI_Barrier(pComm);
    }


    Int nrhs = 5;
    NumMat<T> B(Size(),nrhs);
    Solve(B,pComm);

  }






  template <typename T> void SupernodalMatrix<T>::Solve(NumMat<T> & B, MPI_Comm & pComm){
    Int iam,np;
    MPI_Comm_rank(pComm, &iam);
    MPI_Comm_size(pComm, &np);

    //forward-substitution phase
    //Sending contrib up the tree
    IntNumVec children(Xsuper_.m());
    SetValue(children,0);

    for(Int I=1;I<Xsuper_.m()-1;I++){
      Int fc = Xsuper_(I-1);
      Int lc = Xsuper_(I)-1;

      Int parent = ETree_.PostParent(lc-1);
      if(parent!=0){
        ++children(SupMembership_(parent-1)-1);
      }
    }

    logfileptr->OFS()<<"Children vector is"<<children<<std::endl;


    IntNumVec UpdatesToDo = children;
    std::vector<SuperNode<T> *> contributions(LocalSupernodes_.size());
    std::vector<std::stack<Int> > LocalUpdates(LocalSupernodes_.size());

    //Start from the leaves of the tree
    for(Int I=1;I<Xsuper_.m()-1;I++){
      Int iOwner = Mapping_.Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam ){


        //Create the blocks of my contrib with the same nz structure as L
        //and the same width as the final solution
        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];
        SuperNode<T> * contrib = new SuperNode<T>(I,1,B.n(),cur_snode->MaxHeight());
        contributions[iLocalI-1] = contrib;

        std::vector<Int> GlobToLocIndx(Size(),-1);


        for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){
          NZBlock<T> & cur_nzblk = cur_snode->GetNZBlock(blkidx);
          contrib->AddNZBlock(cur_nzblk.NRows(),B.n(),cur_nzblk.GIndex());
          //fill the indirect index array
          for(Int gi = cur_nzblk.GIndex(); gi<cur_nzblk.GIndex()+cur_nzblk.NRows();++gi){
            GlobToLocIndx[gi-1] = contrib->NZBlockCnt()-1;
          }
        }

        std::vector<Int> isBlockUpdated(contrib->NZBlockCnt(),0);




        //Do all my updates (Local and remote)
        //Local updates
        while(!LocalUpdates[iLocalI-1].empty()){
          Int contrib_snode_id = LocalUpdates[iLocalI-1].top();
          LocalUpdates[iLocalI-1].pop();

          SuperNode<T> & dist_contrib = *contributions[(contrib_snode_id-1) / np];

          logfileptr->OFS()<<"LOCAL Supernode "<<I<<" is updated by contrib of Supernode "<<contrib_snode_id<<std::endl;
          //skip the "diagonal" block
          for(Int blkidx = 1; blkidx<dist_contrib.NZBlockCnt();++blkidx){
            NZBlock<T> & dist_nzblk = dist_contrib.GetNZBlock(blkidx);
            Int local_blkidx = GlobToLocIndx[dist_nzblk.GIndex()-1];
            NZBlock<T> & local_nzblk = contrib->GetNZBlock(local_blkidx);


            Int src_local_fr = max(local_nzblk.GIndex() - dist_nzblk.GIndex(),0);
            Int src_lr = dist_nzblk.GIndex()+dist_nzblk.NRows()-1;

            Int tgt_local_fr = max(dist_nzblk.GIndex() - local_nzblk.GIndex(),0);
            Int tgt_lr = local_nzblk.GIndex()+local_nzblk.NRows()-1;
            Int tgt_local_lr = min(src_lr,tgt_lr) - local_nzblk.GIndex();


            if(tgt_local_fr>0){
              //copy B in the space preceding the updated block
              lapack::Lacpy('N',tgt_local_fr,local_nzblk.NCols(),
                  &B(local_nzblk.GIndex(),0),B.m(),
                  &local_nzblk.Nzval(0,0),local_nzblk.LDA());
            }

            //copy the contrib
            lapack::Lacpy('N',tgt_local_lr - tgt_local_fr +1,dist_nzblk.NCols(),
                &dist_nzblk.Nzval(src_local_fr,0),dist_nzblk.LDA(),
                &local_nzblk.Nzval(tgt_local_fr,0),local_nzblk.LDA());


            if(tgt_local_lr+1<local_nzblk.NRows()){
              //copy B in the space following the updated block
              lapack::Lacpy('N',local_nzblk.NRows()-(tgt_local_lr+1),local_nzblk.NCols(),
                  &B(local_nzblk.GIndex()+(tgt_local_lr+1),0),B.m(),
                  &local_nzblk.Nzval(tgt_local_lr+1,0),local_nzblk.LDA());
            }

            isBlockUpdated[local_blkidx] = 1;

            //this case happens only locally for the diagonal block
            if(src_lr>tgt_lr){
              assert(local_blkidx==0);
              ++local_blkidx;
              NZBlock<T> & local_nzblk = contrib->GetNZBlock(local_blkidx);


              Int src_local_fr = local_nzblk.GIndex() - dist_nzblk.GIndex();
              Int src_lr = dist_nzblk.GIndex()+dist_nzblk.NRows()-1;

              Int tgt_local_fr = max(dist_nzblk.GIndex() - local_nzblk.GIndex(),0);
              Int tgt_lr = local_nzblk.GIndex()+local_nzblk.NRows()-1;
              Int tgt_local_lr = min(src_lr,tgt_lr) - local_nzblk.GIndex();


              if(tgt_local_fr>0){
                //copy B in the space preceding the updated block
                lapack::Lacpy('N',tgt_local_fr,local_nzblk.NCols(),
                    &B(local_nzblk.GIndex(),0),B.m(),
                    &local_nzblk.Nzval(0,0),local_nzblk.LDA());
              }

              //copy the contrib
              lapack::Lacpy('N',tgt_local_lr - tgt_local_fr +1,dist_nzblk.NCols(),
                  &dist_nzblk.Nzval(src_local_fr,0),dist_nzblk.LDA(),
                  &local_nzblk.Nzval(tgt_local_fr,0),local_nzblk.LDA());


              if(tgt_local_lr+1<local_nzblk.NRows()){
                //copy B in the space following the updated block
                lapack::Lacpy('N',local_nzblk.NRows()-(tgt_local_lr+1),local_nzblk.NCols(),
                    &B(local_nzblk.GIndex()+(tgt_local_lr+1),0),B.m(),
                    &local_nzblk.Nzval(tgt_local_lr+1,0),local_nzblk.LDA());
              }

              //there won't be any other update to that block
              isBlockUpdated[local_blkidx] = 1;

            }



          }
            delete contributions[(contrib_snode_id-1) / np];
          --UpdatesToDo(I-1);
        }










        //do remote updates
        std::vector<char> src_nzblocks;
        size_t num_bytes;
        while(UpdatesToDo(I-1)>0){
          //receive children contrib
          logfileptr->OFS()<<UpdatesToDo(I-1)<<" contribs left"<<endl;

          if(src_nzblocks.size()==0){
            num_bytes = sizeof(Int)+sizeof(NZBlock<T>);
            for(Int blkidx=0;blkidx<cur_snode->NZBlockCnt();++blkidx){
              num_bytes += cur_snode->GetNZBlock(blkidx).NRows()*NZBLOCK_ROW_SIZE<T>(B.n());
            }
            logfileptr->OFS()<<"We allocate a buffer of size "<<num_bytes<<std::endl;
            src_nzblocks.resize(num_bytes);
          }


          //MPI_Recv
          MPI_Status recv_status;

          MPI_Probe(MPI_ANY_SOURCE,I,pComm,&recv_status);
          int bytes_to_receive = 0;
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_to_receive);

          assert(src_nzblocks.size()>=bytes_to_receive);


          char * recv_buf = &src_nzblocks[0]+ max(sizeof(Int),NZBLOCK_OBJ_SIZE<T>())-min(sizeof(Int),NZBLOCK_OBJ_SIZE<T>());
          MPI_Recv(recv_buf,&*src_nzblocks.end()-recv_buf,MPI_BYTE,MPI_ANY_SOURCE,I,pComm,&recv_status);
          Int src_snode_id = *(Int*)recv_buf;

          //Resize the buffer to the actual number of bytes received
          int bytes_received = 0;
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
          src_nzblocks.resize(recv_buf+bytes_received - &src_nzblocks[0]);

          //Create the dummy supernode for that data
          SuperNode<T> dist_contrib(src_snode_id,1,B.n(),&src_nzblocks);

          logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib.Id()<<std::endl;

          for(Int blkidx = 0; blkidx<dist_contrib.NZBlockCnt();++blkidx){
            NZBlock<T> & dist_nzblk = dist_contrib.GetNZBlock(blkidx);

            Int local_blkidx = GlobToLocIndx[dist_nzblk.GIndex()-1];
            NZBlock<T> & local_nzblk = contrib->GetNZBlock(local_blkidx);


            Int src_local_fr = max(local_nzblk.GIndex() - dist_nzblk.GIndex(),0);
            Int src_lr = dist_nzblk.GIndex()+dist_nzblk.NRows()-1;

            Int tgt_local_fr = max(dist_nzblk.GIndex() - local_nzblk.GIndex(),0);
            Int tgt_lr = local_nzblk.GIndex()+local_nzblk.NRows()-1;
            Int tgt_local_lr = min(src_lr,tgt_lr) - local_nzblk.GIndex();


            if(tgt_local_fr>0){
              //copy B in the space preceding the updated block
              lapack::Lacpy('N',tgt_local_fr,local_nzblk.NCols(),
                  &B(local_nzblk.GIndex(),0),B.m(),
                  &local_nzblk.Nzval(0,0),local_nzblk.LDA());
            }

            //copy the contrib
            lapack::Lacpy('N',tgt_local_lr - tgt_local_fr +1,dist_nzblk.NCols(),
                &dist_nzblk.Nzval(src_local_fr,0),dist_nzblk.LDA(),
                &local_nzblk.Nzval(tgt_local_fr,0),local_nzblk.LDA());


            if(tgt_local_lr+1<local_nzblk.NRows()){
              //copy B in the space following the updated block
              lapack::Lacpy('N',local_nzblk.NRows()-(tgt_local_lr+1),local_nzblk.NCols(),
                  &B(local_nzblk.GIndex()+(tgt_local_lr+1),0),B.m(),
                  &local_nzblk.Nzval(tgt_local_lr+1,0),local_nzblk.LDA());
            }

            isBlockUpdated[local_blkidx] = 1;
          }
          --UpdatesToDo(I-1);

          src_nzblocks.resize(num_bytes);
        }


        if(UpdatesToDo(I-1)==0){
          for(Int blkidx = 0; blkidx<contrib->NZBlockCnt();++blkidx){
            NZBlock<T> & cur_nzblk = contrib->GetNZBlock(blkidx);
            NZBlock<T> & chol_nzblk = cur_snode->GetNZBlock(blkidx);
            if(!isBlockUpdated[blkidx]){
              //copy B in it
              lapack::Lacpy('N',cur_nzblk.NRows(),cur_nzblk.NCols(),
                  &B(cur_nzblk.GIndex(),0),B.m(),
                  cur_nzblk.Nzval(),cur_nzblk.LDA());
            }

            //compute my contribution
            if(blkidx==0){
              //if we are processing the "pivot" block
              for(Int k = 0; k<cur_nzblk.NRows();++k){
                for(Int j = 0; j<cur_nzblk.NCols();++j){
                  cur_nzblk.Nzval(k,j) /= chol_nzblk.Nzval(k,k);
                }
              }

//              blas::Trsm('R','L','N','N',cur_nzblk.NRows(),chol_nzblk.NCols(),
//                  ONE<T>(),chol_nzblk.Nzval(),chol_nzblk.LDA(),
//                  cur_nzblk.Nzval(),cur_nzblk.LDA());
            }
            else{
              NZBlock<T> & pivot_nzblk = contrib->GetNZBlock(0);
              blas::Gemm('N','N',chol_nzblk.NRows(),pivot_nzblk.NCols(),
                  pivot_nzblk.NRows(),MINUS_ONE<T>(), 
                  chol_nzblk.Nzval(), chol_nzblk.LDA(),
                  pivot_nzblk.Nzval(), pivot_nzblk.LDA(),
                  ONE<T>(), cur_nzblk.Nzval(), cur_nzblk.LDA());
            } 
          }

          //send to my parent
          Int parent = ETree_.PostParent(cur_snode->LastCol()-1);
          Int parent_snode_id = SupMembership_[parent-1];

          Int iTarget = Mapping_.Map(parent_snode_id-1,parent_snode_id-1);

          if(iTarget!=iam){
            logfileptr->OFS()<<"Supernode "<<parent_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
            Int J = parent_snode_id;

            Int tgt_first_col = Xsuper_(J-1);
            Int tgt_last_col = Xsuper_(J)-1;


#ifdef ROW_MAJOR
            //          int local_first_row = src_first_row-pivot_nzblk->GIndex();
            //
            //          char * start_buf = reinterpret_cast<char*>(&pivot_nzblk->Nzval(local_first_row,0));
            //          size_t num_bytes = src_snode.End() - start_buf;
            //
            //          //Create a header
            //          NZBlockHeader<T> * header = new NZBlockHeader<T>(pivot_nzblk->NRows()-local_first_row,pivot_nzblk->NCols(),src_first_row);
            //
            //          MPI_Datatype type;
            //          int lens[3];
            //          MPI_Aint disps[3];
            //          MPI_Datatype types[3];
            //          Int err = 0;
            //
            //
            //          /* define a struct that describes all our data */
            //          lens[0] = sizeof(Int);
            //          lens[1] = sizeof(NZBlockHeader<T>);
            //          lens[2] = num_bytes;
            //          MPI_Address(&I, &disps[0]);
            //          MPI_Address(header, &disps[1]);
            //          MPI_Address(start_buf, &disps[2]);
            //          types[0] = MPI_BYTE;
            //          types[1] = MPI_BYTE;
            //          types[2] = MPI_BYTE;
            //          MPI_Type_struct(3, lens, disps, types, &type);
            //          MPI_Type_commit(&type);
            //
            //          //logfileptr->OFS()<<src_snode<<endl;
            //          
            //          /* send to target */
            //          //FIXME : maybe this can be replaced by a scatterv ?
            //          MPI_Send(MPI_BOTTOM,1,type,iTarget,tgt_snode_id,pComm);
            //
            //          MPI_Type_free(&type);
            //
            //          delete header;
            //          logfileptr->OFS()<<"     Send factor "<<I<<" to node"<<J<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<" ("<<num_bytes+sizeof(*header)+sizeof(Int)<<" bytes)"<<std::endl;
#else
#ifdef SEND_WHOLE_BLOCK
#else

            char * start_buf = reinterpret_cast<char*>(&contrib->GetNZBlock(1));
            size_t num_bytes = contrib->End() - start_buf;


            MPI_Datatype type;
            int lens[2];
            MPI_Aint disps[2];
            MPI_Datatype types[2];
            Int err = 0;


            /* define a struct that describes all our data */
            lens[0] = sizeof(Int);
            lens[1] = num_bytes;
            MPI_Address(&I, &disps[0]);
            MPI_Address(start_buf, &disps[1]);
            types[0] = MPI_BYTE;
            types[1] = MPI_BYTE;
            MPI_Type_struct(2, lens, disps, types, &type);
            MPI_Type_commit(&type);

            MPI_Send(MPI_BOTTOM,1,type,iTarget,parent_snode_id,pComm);

            logfileptr->OFS()<<"     SEND Supernode "<<I<<" to Supernode "<<J<<" on P"<<iTarget<<" ("<<lens[0] + lens[1] <<" bytes)"<<std::endl;
            MPI_Type_free(&type);
#endif // end ifdef SEND_WHOLE_BLOCK

#endif // end ifdef ROW_MAJOR

            //delete it
            delete contributions[iLocalI-1];

          }
          else{
            Int iLocalJ = (parent_snode_id-1) / np +1 ;
            LocalUpdates[iLocalJ-1].push(I);
          }

        }


      }
    }

    //Back-substitution phase

  }


} // namespace LIBCHOLESKY


#endif 
