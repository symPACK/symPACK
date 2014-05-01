#ifndef _SUPERNODAL_MATRIX_IMPL_HPP_
#define _SUPERNODAL_MATRIX_IMPL_HPP_

#include "SupernodalMatrix.hpp"

#define _DEBUG_

#define TAG_FACTOR 0

namespace LIBCHOLESKY{


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

    //ETree_.SortChildren(cc);

    logfileptr->OFS()<<"colcnt "<<cc<<std::endl;

    Global_.FindSupernodes(ETree_,cc,SupMembership_,Xsuper_,1);

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
        LocalSupernodes_.push_back( new SuperNode<T>(I,fc,lc,Size(),iHeight));
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
          Int iNextColptr = Global_.colptr(i);
          Int iRowind = Global_.rowind(iGcolptr-1);
          Int iNrows = Global_.colptr(i) - Global_.colptr(i-1);
          Int idxA = 0;
          Int iLRow = 0;
          Int firstrow = fi + i-fc;
          for(Int idx = firstrow; idx<=li;idx++){
            iLRow = lindx_(idx-1);
            //logfileptr->OFS()<<"Looking at L("<<iLRow<<","<<i<<")"<<std::endl;
            if( iLRow == iRowind){
              //TODO we can use the Globaltolocal index here also
              Int iNZBlockIdx = snode.FindBlockIdx(iLRow);

              T elem = pdNzVal[iStartIdxCopy + idxA];

                NZBlockDesc & desc = snode.GetNZBlockDesc(iNZBlockIdx);
                T * dest = snode.GetNZval(desc.Offset);

                //find where we should put it
                Int localCol = i - fc;
                Int localRow = iLRow - desc.GIndex;

                //row major storage
                dest[localRow*snode.Size()+localCol] = elem;
  
//              NZBlock<T> & pDestBlock = snode.GetNZBlock(iNZBlockIdx);
//              //logfileptr->OFS()<<*pDestBlock<<std::endl;
//
//              //find where we should put it
//              Int localCol = i - fc;
//              Int localRow = iLRow - pDestBlock.GIndex();
//#ifdef _DEBUG_
//              logfileptr->OFS()<<"Elem is A("<<iRowind<<","<<i<<") = "<<elem<<" at ("<<localRow<<","<<localCol<<") in "<<iNZBlockIdx<<"th NZBlock of L"<< std::endl;
//#endif
//              pDestBlock.Nzval(localRow,localCol) = elem;


              if(iGcolptr+idxA+1<iNextColptr){
                idxA++;
                iRowind = Global_.rowind(iGcolptr+idxA-1);
              }
              else{
                break;
              }
            }
          }

        }

      }
#ifdef _DEBUG_
      logfileptr->OFS()<<"--------------------------------------------------"<<std::endl;
#endif
    }

    for(Int I=1;I<Xsuper_.m();I++){
      Int src_first_col = Xsuper_(I-1);
      Int src_last_col = Xsuper_(I)-1;
      Int iOwner = Mapping_.Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam ){



        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> & src_snode = *LocalSupernodes_[iLocalI -1];


        logfileptr->OFS()<<"+++++++++++++"<<I<<"++++++++"<<std::endl;
        for(int blkidx=0;blkidx<src_snode.NZBlockCnt();++blkidx){

          NZBlockDesc & desc = src_snode.GetNZBlockDesc(blkidx);
          T * val = src_snode.GetNZval(desc.Offset);
          Int nRows = src_snode.NRows(blkidx);


          for(Int i = 0; i< nRows; ++i){
            for(Int j = 0; j< src_snode.Size(); ++j){
              logfileptr->OFS()<<val[i*src_snode.Size()+j]<<" ";
            }
            logfileptr->OFS()<<std::endl;
          }

//        logfileptr->OFS()<<"................................"<<std::endl;
//
//          NZBlock<T> & nzblk = src_snode.GetNZBlock(blkidx);
//          for(Int i = 0; i< nzblk.NRows(); ++i){
//            for(Int j = 0; j< nzblk.NCols(); ++j){
//              logfileptr->OFS()<<nzblk.Nzval(i,j)<<" ";
//            }
//            logfileptr->OFS()<<std::endl;
//          }
        logfileptr->OFS()<<"_______________________________"<<std::endl;
        }

//        Int nrhs = 5;
//        NumMat<T> RHS(dense_snode.n(),nrhs);
//        NumMat<T> XTrue(dense_snode.n(),nrhs);
//        UniformRandom(XTrue);
//        blas::Gemm('N','N',dense_snode.n(),nrhs,dense_snode.n(),1.0,&dense_snode(0,0),dense_snode.m(),&XTrue(0,0),XTrue.m(),0.0,&RHS(0,0),RHS.m());




      }
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
/*  template <typename T> inline bool SupernodalMatrix<T>::FindNextUpdate(Int src_snode_id, Int & src_first_row, Int & src_last_row, Int & tgt_snode_id){
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
*/

  template <typename T> inline bool SupernodalMatrix<T>::FindNextUpdate(SuperNode<T> & src_snode, Int & src_nzblk_idx, Int & src_first_row, Int & src_last_row, Int & tgt_snode_id){
    //src_nzblk_idx is the last nzblock index examined
    //src_first_row is the first row updating the supernode examined
    //src_last_row is the last row updating the supernode examined

    

    //if tgt_snode_id == 0 , this is the first call to the function
    if(tgt_snode_id == 0){

      //find the first sub diagonal block
      src_nzblk_idx = -1;
      for(Int blkidx = 0; blkidx < src_snode.NZBlockCnt(); ++blkidx){
        if(src_snode.GetNZBlockDesc(blkidx).GIndex > src_snode.LastCol()){
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
            assert(src_nzblk_idx<src_snode.NZBlockCnt());
        src_first_row = src_snode.GetNZBlockDesc(src_nzblk_idx).GIndex;
      }
    }
    else{
      //find the block corresponding to src_last_row
      src_nzblk_idx = src_snode.FindBlockIdx(src_last_row);
      if(src_nzblk_idx==-1/*src_snode.NZBlockCnt()*/){
        return false;
      }
      else{
            assert(src_nzblk_idx<src_snode.NZBlockCnt());
        NZBlockDesc & desc = src_snode.GetNZBlockDesc(src_nzblk_idx);
        Int src_lr = desc.GIndex + src_snode.NRows(src_nzblk_idx)-1;

//        NZBlock<T> & nzblk = src_snode.GetNZBlock(src_nzblk_idx);
//        Int src_lr = nzblk.GIndex()+nzblk.NRows()-1;
        if(src_last_row == src_lr){
          src_nzblk_idx++;
          if(src_nzblk_idx==src_snode.NZBlockCnt()){
            return false;
          }
          else{
            assert(src_nzblk_idx<src_snode.NZBlockCnt());
            src_first_row = src_snode.GetNZBlockDesc(src_nzblk_idx).GIndex;
          }
        }
        else{
          src_first_row = src_last_row+1;
        }
      }
    }

    //Now we try to find src_last_row
    NZBlockDesc & desc = src_snode.GetNZBlockDesc(src_nzblk_idx);
    assert(src_first_row >= desc.GIndex);
    Int src_fr = max(src_first_row,desc.GIndex);
    Int src_lr = desc.GIndex+src_snode.NRows(src_nzblk_idx)-1;

//    NZBlock<T> & nzblk = src_snode.GetNZBlock(src_nzblk_idx);
//    assert(src_first_row >= nzblk.GIndex());
//    Int src_fr = max(src_first_row,nzblk.GIndex());
//    Int src_lr = nzblk.GIndex()+nzblk.NRows()-1;


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
      if(new_blk_idx==-1/*>=src_snode.NZBlockCnt()*/){
        //src_last_row is the last row of the last nzblock
        new_blk_idx = src_snode.NZBlockCnt()-1; 
      }

      NZBlockDesc & last_desc = src_snode.GetNZBlockDesc(new_blk_idx); 
      src_last_row = min(src_last_row,last_desc.GIndex + src_snode.NRows(new_blk_idx) - 1);

//      NZBlock<T> & last_block = src_snode.GetNZBlock(new_blk_idx); 
//      src_last_row = min(src_last_row,last_block.GIndex() + last_block.NRows() - 1);
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



  template <typename T> void SupernodalMatrix<T>::UpdateSuperNode(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode, Int &pivot_idx, Int  pivot_fr = I_ZERO){
    NZBlockDesc & pivot_desc = src_snode.GetNZBlockDesc(pivot_idx);
    if(pivot_fr ==I_ZERO){
      pivot_fr = pivot_desc.GIndex;
    }
    assert(pivot_fr >= pivot_desc.GIndex);
    Int pivot_nrows = src_snode.NRows(pivot_idx);
    Int pivot_lr = min(pivot_desc.GIndex + pivot_nrows -1, pivot_fr + tgt_snode.Size() -1);
    T * pivot = &(src_snode.GetNZval(pivot_desc.Offset)[(pivot_fr-pivot_desc.GIndex)*src_snode.Size()]);

//    NZBlock<T> & pivot_nzblk = src_snode.GetNZBlock(pivot_idx);
//    if(pivot_fr ==I_ZERO){
//      pivot_fr = pivot_nzblk.GIndex();
//    }
//    assert(pivot_fr >= pivot_nzblk.GIndex());
//    Int pivot_lr = min(pivot_nzblk.GIndex() +pivot_nzblk.NRows() -1, pivot_fr + tgt_snode.Size() -1);
//    T * pivot = & pivot_nzblk.Nzval(pivot_fr-pivot_nzblk.GIndex(),0);

    Int tgt_updated_fc =  pivot_fr - tgt_snode.FirstCol();

    for(int src_idx=pivot_idx;src_idx<src_snode.NZBlockCnt();++src_idx){

      NZBlockDesc & src_desc = src_snode.GetNZBlockDesc(src_idx);
      Int src_fr = max(pivot_fr, src_desc.GIndex) ;
      Int src_nrows = src_snode.NRows(src_idx);
      Int src_lr = src_desc.GIndex+src_nrows-1;

//      NZBlock<T> & src_nzblk = src_snode.GetNZBlock(src_idx);
//      Int src_fr = max(pivot_fr, src_nzblk.GIndex()) ;
//      Int src_lr = src_nzblk.GIndex()+src_nzblk.NRows()-1;

      do{
        //TODO Need to be replaced by GlobToLoc index
        Int tgt_idx = tgt_snode.FindBlockIdx(src_fr);

        assert(tgt_idx!=-1 && tgt_idx<tgt_snode.NZBlockCnt());

        NZBlockDesc & tgt_desc = tgt_snode.GetNZBlockDesc(tgt_idx);
        Int tgt_nrows = tgt_snode.NRows(tgt_idx);
        Int tgt_fr = tgt_desc.GIndex;
        Int tgt_lr = tgt_desc.GIndex+tgt_nrows-1;

//        NZBlock<T> & tgt_nzblk = tgt_snode.GetNZBlock(tgt_idx);
//        Int tgt_fr = tgt_nzblk.GIndex();
//        Int tgt_lr = tgt_nzblk.GIndex()+tgt_nzblk.NRows()-1;

        Int update_fr = max(tgt_fr, src_fr);
        Int update_lr = min(tgt_lr, src_lr);

        assert(update_fr >= tgt_fr); 
        assert(update_lr >= update_fr); 
        assert(update_lr <= tgt_lr); 

        //Update tgt_nzblk with src_nzblk
        //        logfileptr->OFS()<<"Updating SNODE "<<tgt_snode.Id()<<" Block "<<tgt_idx<<" ["<<update_fr<<".."<<update_lr<<"] with SNODE "<<src_snode.Id()<<" Block "<<src_idx<<" ["<<update_fr<<".."<<update_lr<<"]"<<endl;

        T * src = &(src_snode.GetNZval(src_desc.Offset)[(update_fr - src_desc.GIndex)*src_snode.Size()]);
        T * tgt = &(tgt_snode.GetNZval(tgt_desc.Offset)[(update_fr - tgt_desc.GIndex)*tgt_snode.Size()+tgt_updated_fc]);
//        T * src = &src_nzblk.Nzval(update_fr - src_nzblk.GIndex() ,0);
//        T * tgt = &tgt_nzblk.Nzval(update_fr - tgt_nzblk.GIndex() ,tgt_updated_fc);

#ifdef _DEBUG_
          logfileptr->OFS()<<"| "<<tgt[0]<<" | = | "<<tgt[0]<<" | - |"<<src[0]<<"| * | ";
          for(Int i =0;i<pivot_lr-pivot_fr+1 ;++i){
            logfileptr->OFS()<<pivot[i]<<" ";
          }
          logfileptr->OFS()<<"|"<<std::endl;
          for(Int i =1;i<update_lr - update_fr + 1 ;++i){
            logfileptr->OFS()<<"| "<<tgt[i]<<" |   | "<<tgt[i]<<" | - |"<<src[i]<<"|"<<std::endl;
          }
#endif

        blas::Gemm('T','N',pivot_lr-pivot_fr+1, update_lr - update_fr + 1,src_snode.Size(),MINUS_ONE<T>(),pivot,src_snode.Size(),src,src_snode.Size(),ONE<T>(),tgt,tgt_snode.Size());

#ifdef ROW_MAJOR
//        blas::Gemm('T','N',pivot_lr-pivot_fr+1, update_lr - update_fr + 1,src_snode.Size(),1.0,pivot,pivot_nzblk.LDA(),src,src_nzblk.LDA(),1.0,tgt,tgt_nzblk.LDA());
#else
//        blas::Gemm('N','T', update_lr - update_fr + 1,pivot_lr-pivot_fr+1,  src_snode.Size(),MINUS_ONE<T>(),src,src_nzblk.LDA(),pivot,pivot_nzblk.LDA(),ONE<T>(),tgt,tgt_nzblk.LDA());
#endif

//        if(tgt_snode.Id()==46 && tgt_idx==1){
//          logfileptr->OFS()<<*tgt<<std::endl;
//        }
//        //special case for diagonal block
//        if(src_lr > tgt_lr){
//          assert(tgt_idx==0);
//          //keep the same source block, just change first line
//          src_fr = tgt_lr+1;
//        }
        
        if(tgt_idx+1<tgt_snode.NZBlockCnt()){
          NZBlockDesc & next_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_idx+1);
          src_fr = next_tgt_desc.GIndex;

//          NZBlock<T> & next_tgt_nzblk = tgt_snode.GetNZBlock(tgt_idx+1);
//          Int next_tgt_fr = next_tgt_nzblk.GIndex();
//          src_fr = next_tgt_fr;
        }
        else{
          break;
        }
      }while(src_fr<=src_lr);
    }

  }




  template <typename T> void SupernodalMatrix<T>::Factorize( MPI_Comm & pComm ){
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
        assert(src_snode.Id() == I); 


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
          
#ifdef _DEBUG_
          logfileptr->OFS()<<"LOCAL Supernode "<<I<<" is updated by Supernode "<<src_snode_id<<" from row "<<src_first_row<<" "<<src_nzblk_idx<<std::endl;
#endif

          UpdateSuperNode(local_src_snode,src_snode,src_nzblk_idx, src_first_row);

          --UpdatesToDo(I-1);
#ifdef _DEBUG_
          logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left"<<endl;
#endif
        }


//        //Remote updates
//        std::vector<char> src_nzblocks;
//
//        size_t num_bytes;
//        while(UpdatesToDo(I-1)>0){
//          logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left"<<endl;
//
//          if(src_nzblocks.size()==0){
//            //The upper bound must be of the width of the "largest" child
//            logfileptr->OFS()<<"Maximum width is "<<UpdateWidth_(I-1)<<std::endl;
//
//            num_bytes = sizeof(Int)+sizeof(NZBlock<T>);
//            for(Int blkidx=0;blkidx<src_snode.NZBlockCnt();++blkidx){
//              num_bytes += src_snode.GetNZBlock(blkidx).NRows()*NZBLOCK_ROW_SIZE<T>(UpdateWidth_(I-1));
//            }
//            logfileptr->OFS()<<"We allocate a buffer of size "<<num_bytes<<std::endl;
//            src_nzblocks.resize(num_bytes);
//          }
//
//
//          //receive the index array
//
//
//          //MPI_Recv
//          MPI_Status recv_status;
//
//          MPI_Probe(MPI_ANY_SOURCE,I,pComm,&recv_status);
//          int bytes_to_receive = 0;
//          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_to_receive);
//          assert(src_nzblocks.size()>=bytes_to_receive);
//
//
//
//
//          char * recv_buf = &src_nzblocks[0]+ max(sizeof(Int),NZBLOCK_OBJ_SIZE<T>())-min(sizeof(Int),NZBLOCK_OBJ_SIZE<T>());
//          MPI_Recv(recv_buf,&*src_nzblocks.end()-recv_buf,MPI_BYTE,MPI_ANY_SOURCE,I,pComm,&recv_status);
//          //logfileptr->OFS()<<"Received something"<<endl;
//          Int src_snode_id = *(Int*)recv_buf;
//
//
//          //Resize the buffer to the actual number of bytes received
//          int bytes_received = 0;
//          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
//          src_nzblocks.resize(recv_buf+bytes_received - &src_nzblocks[0]);
//
//
//          //Create the dummy supernode for that data
//          SuperNode<T> dist_src_snode(src_snode_id,Xsuper_[src_snode_id-1],Xsuper_[src_snode_id]-1,&src_nzblocks);
//
//          logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode.Id()<<std::endl;
//          //      logfileptr->OFS()<<dist_src_snode<<endl;
//
//
//          //Update everything I own with that factor
//          //Update the ancestors
//          Int tgt_snode_id = 0;
//          Int src_first_row = 0;
//          Int src_last_row = 0;
//          Int src_nzblk_idx = 0;
//          while(FindNextUpdate(dist_src_snode, src_nzblk_idx, src_first_row, src_last_row, tgt_snode_id)){ 
//            Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
//            if(iTarget == iam){
//              logfileptr->OFS()<<"RECV Supernode "<<tgt_snode_id<<" is updated by Supernode "<<dist_src_snode.Id()<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;
//
//              Int iLocalJ = (tgt_snode_id-1) / np +1 ;
//              SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];
//
//              UpdateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, src_first_row);
//
//              --UpdatesToDo(tgt_snode_id-1);
//              logfileptr->OFS()<<UpdatesToDo(tgt_snode_id-1)<<" updates left for Supernode "<<tgt_snode_id<<endl;
//            }
//          }
//
//
//
//          //restore to its capacity
//          src_nzblocks.resize(num_bytes);
//
//        }
//        //clear the buffer
//        { vector<char>().swap(src_nzblocks);  }
//



        //Remote updates
        std::vector<T> src_nzval;
        std::vector<char> src_blocks;

        Int nz_cnt;
        Int max_bytes;
        while(UpdatesToDo(I-1)>0){
#ifdef _DEBUG_
          logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left"<<endl;
#endif

          if(src_blocks.size()==0){
            //The upper bound must be of the width of the "largest" child
#ifdef _DEBUG_
            logfileptr->OFS()<<"Maximum width is "<<UpdateWidth_(I-1)<<std::endl;
#endif

            Int nrows = src_snode.NRowsBelowBlock(0);
            Int ncols = UpdateWidth_(I-1);
            nz_cnt = nrows * ncols;

            max_bytes = nrows*sizeof(NZBlockDesc) + sizeof(Int);
            src_blocks.resize(max_bytes);
            src_nzval.resize(nz_cnt);
          }

          

          //MPI_Recv
          MPI_Status recv_status;

          //receive the index array
          MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,MPI_ANY_SOURCE,I,pComm,&recv_status);
          int bytes_received = 0;
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
          //there an aditional integer storing the src id which needs to be removed
          Int src_nzblk_cnt = (bytes_received - sizeof(Int) ) / sizeof(NZBlockDesc);

          //receive the nzval array
          MPI_Recv(&src_nzval[0],nz_cnt*sizeof(T),MPI_BYTE,MPI_ANY_SOURCE,I,pComm,&recv_status);
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
          Int src_nzval_cnt = bytes_received / sizeof(T);


          Int src_snode_id = *(Int*)&src_blocks[0];
          NZBlockDesc * src_blocks_ptr = 
                    reinterpret_cast<NZBlockDesc*>(&src_blocks[sizeof(Int)]);
          T * src_nzval_ptr = &src_nzval[0];

          //Create the dummy supernode for that data
          SuperNode<T> dist_src_snode(src_snode_id,Xsuper_[src_snode_id-1],Xsuper_[src_snode_id]-1, Size(), src_blocks_ptr, src_nzblk_cnt, src_nzval_ptr, src_nzval_cnt);

#ifdef _DEBUG_
          logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode.Id()<<std::endl;
#endif
          //Update everything I own with that factor
          //Update the ancestors
          Int tgt_snode_id = 0;
          Int src_first_row = 0;
          Int src_last_row = 0;
          Int src_nzblk_idx = 0;


          while(FindNextUpdate(dist_src_snode, src_nzblk_idx, src_first_row, src_last_row, tgt_snode_id)){ 
            Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
            if(iTarget == iam){
#ifdef _DEBUG_
              logfileptr->OFS()<<"RECV Supernode "<<tgt_snode_id<<" is updated by Supernode "<<dist_src_snode.Id()<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;
#endif

              Int iLocalJ = (tgt_snode_id-1) / np +1 ;
              SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];

              UpdateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, src_first_row);

              --UpdatesToDo(tgt_snode_id-1);
#ifdef _DEBUG_
              logfileptr->OFS()<<UpdatesToDo(tgt_snode_id-1)<<" updates left for Supernode "<<tgt_snode_id<<endl;
#endif
            }
          }



          //restore to its capacity
//          src_nzblocks.resize(num_bytes);

        }
        //clear the buffer
//        { vector<char>().swap(src_nzblocks);  }








#ifdef _DEBUG_
        logfileptr->OFS()<<"  Factoring Supernode "<<I<<std::endl;
#endif

        //Factorize Diagonal block
        NZBlockDesc & diag_desc = src_snode.GetNZBlockDesc(0);
        T * diag_nzval = src_snode.GetNZval(diag_desc.Offset);
        lapack::Potrf( 'U', src_snode.Size(), diag_nzval, src_snode.Size());

        if(src_snode.NZBlockCnt()>1){
          NZBlockDesc & nzblk_desc = src_snode.GetNZBlockDesc(1);
          T * nzblk_nzval = src_snode.GetNZval(nzblk_desc.Offset);
          blas::Trsm('L','U','N','N',src_snode.Size(),  src_snode.NRowsBelowBlock(1), ONE<T>(),  diag_nzval, src_snode.Size(), nzblk_nzval, src_snode.Size());

//correct col-major
//          blas::Trsm('R','L','T','N',nzblk.NRows(),src_snode.Size(), ONE<T>(),  diagonalBlock.Nzval(), diagonalBlock.LDA(), nzblk.Nzval(), nzblk.LDA());
        }

//        NZBlock<T> & diagonalBlock = src_snode.GetNZBlock(0);
//
//#ifdef ROW_MAJOR
//        lapack::Potrf( 'U', src_snode.Size(), diagonalBlock.Nzval(), diagonalBlock.LDA());
//#else
//        lapack::Potrf( 'L', src_snode.Size(), diagonalBlock.Nzval(), diagonalBlock.LDA());
//#endif
//
//        for(int blkidx=1;blkidx<src_snode.NZBlockCnt();++blkidx){
//          NZBlock<T> & nzblk = src_snode.GetNZBlock(blkidx);
//
//          //Update lower triangular blocks
//#ifdef ROW_MAJOR
//          blas::Trsm('R','U','T','N',nzblk.NCols(),src_snode.Size(), ONE<T>(),  diagonalBlock.Nzval(), diagonalBlock.LDA(), nzblk.Nzval(), nzblk.LDA());
//#else
//          blas::Trsm('R','L','T','N',nzblk.NRows(),src_snode.Size(), ONE<T>(),  diagonalBlock.Nzval(), diagonalBlock.LDA(), nzblk.Nzval(), nzblk.LDA());
//#endif
//        }



        //Send my factor to my ancestors. 
        BolNumVec is_factor_sent(np);
        SetValue(is_factor_sent,false);

        std::vector<std::vector<char> *> blocks_sent;


        Int tgt_snode_id = 0;
        Int src_nzblk_idx = 0;
        Int src_first_row = 0;
        Int src_last_row = 0;
        while(FindNextUpdate(src_snode, src_nzblk_idx, src_first_row, src_last_row, tgt_snode_id)){ 

          Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
          if(iTarget != iam){
            if(!is_factor_sent[iTarget]){
              is_factor_sent[iTarget] = true;


#ifdef _DEBUG_
              logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated by Supernode "<<I<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;
#endif

              std::vector<char> * pNewDesc = new std::vector<char>();
              blocks_sent.push_back(pNewDesc);


              Int tgt_first_col = Xsuper_(tgt_snode_id-1);
              Int tgt_last_col = Xsuper_(tgt_snode_id)-1;

              //Send
              NZBlockDesc & pivot_desc = src_snode.GetNZBlockDesc(src_nzblk_idx);

              Int local_first_row = src_first_row-pivot_desc.GIndex;
              Int nzblk_cnt = src_snode.NZBlockCnt()-src_nzblk_idx;
              pNewDesc->resize(sizeof(Int) + nzblk_cnt*sizeof(NZBlockDesc));

              char * send_ptr = &(*pNewDesc)[0];

              Int * id_ptr = reinterpret_cast<Int *>(&(*pNewDesc)[0]);
              NZBlockDesc * block_desc_ptr = 
                reinterpret_cast<NZBlockDesc *>(&(*pNewDesc)[sizeof(Int)]);

              *id_ptr = src_snode.Id();
              //copy the block descriptors
              std::copy(&pivot_desc,&pivot_desc+nzblk_cnt, block_desc_ptr);
              //change the first one
              block_desc_ptr->Offset += (src_first_row - block_desc_ptr->GIndex)*src_snode.Size();
              block_desc_ptr->GIndex = src_first_row;

              //send the block descriptors
              MPI_Send(send_ptr,nzblk_cnt*sizeof(NZBlockDesc) + sizeof(Int),
                  MPI_BYTE,iTarget,tgt_snode_id,pComm);

              T * nzval_ptr = src_snode.GetNZval(pivot_desc.Offset
                  +local_first_row*src_snode.Size());

              Int nz_cnt = (src_snode.NRowsBelowBlock(src_nzblk_idx)
                  - local_first_row )*src_snode.Size();


              MPI_Send(nzval_ptr,nz_cnt*sizeof(T),
                  MPI_BYTE,iTarget,tgt_snode_id,pComm);

              delete blocks_sent.back();

#ifdef _DEBUG_            
              logfileptr->OFS()<<"     Send factor "<<I<<" to node"<<tgt_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
              logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
#endif


            }
          }
          else{
            Int iLocalJ = (tgt_snode_id-1) / np +1 ;

            assert(LocalSupernodes_[iLocalJ-1]->Id()==tgt_snode_id);

            LocalUpdates[iLocalJ-1].push(src_first_row);
            LocalUpdates[iLocalJ-1].push(src_nzblk_idx);
            LocalUpdates[iLocalJ-1].push(I);

#ifdef _DEBUG_
            logfileptr->OFS()<<"Local Supernode "<<tgt_snode_id<<" is updated by Supernode "<<I<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;
#endif

          }
        }
      }
      MPI_Barrier(pComm);
    }

    NumMat<T> fullMatrix;
    if(iam==0){
      fullMatrix.Resize(this->Size(),this->Size());
      SetValue(fullMatrix,ZERO<T>());
    }

    //output L
    for(Int I=1;I<Xsuper_.m();I++){
      Int src_first_col = Xsuper_(I-1);
      Int src_last_col = Xsuper_(I)-1;
      Int iOwner = Mapping_.Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam  && iam != 0){



        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> & src_snode = *LocalSupernodes_[iLocalI -1];

        logfileptr->OFS()<<"Supernode "<<I<<"("<<src_snode.Id()<<") is on P"<<iOwner<<" local index is "<<iLocalI<<std::endl; 
        logfileptr->OFS()<<src_snode<<std::endl;

        NZBlockDesc * nzblk_desc = &src_snode.GetNZBlockDesc(0);
        Int size_blocks = src_snode.NZBlockCnt();
        MPI_Send(&size_blocks,sizeof(Int),MPI_BYTE,0,I,pComm);
        MPI_Send(nzblk_desc,size_blocks*sizeof(NZBlockDesc),MPI_BYTE,0,I,pComm);

        T * nzblk_nzval = src_snode.GetNZval(0);
        Int size_nzval = src_snode.NRowsBelowBlock(0)*src_snode.Size();
        MPI_Send(&size_nzval,sizeof(Int),MPI_BYTE,0,I,pComm);
        MPI_Send(nzblk_nzval,size_nzval*sizeof(T),MPI_BYTE,0,I,pComm);



      }

      if(iam==0){
        if(iOwner != iam){ 
          Int snode_size = Xsuper_[I] - Xsuper_[I-1];

          Int size_blocks, size_nzval;
          std::vector<NZBlockDesc> blocks;
          std::vector<T> nzval;

          MPI_Recv(&size_blocks,sizeof(Int),MPI_BYTE,iOwner,I,pComm,MPI_STATUS_IGNORE);
          blocks.resize(size_blocks); 
          MPI_Recv(&blocks[0],size_blocks*sizeof(NZBlockDesc),MPI_BYTE,iOwner,I,pComm,MPI_STATUS_IGNORE);


          MPI_Recv(&size_nzval,sizeof(Int),MPI_BYTE,iOwner,I,pComm,MPI_STATUS_IGNORE);
          nzval.resize(size_nzval); 
          MPI_Recv(&nzval[0],size_nzval*sizeof(T),MPI_BYTE,iOwner,I,pComm,MPI_STATUS_IGNORE);

          SuperNode<T> src_snode(I,Xsuper_[I-1],Xsuper_[I]-1,Size(),&blocks[0],size_blocks,&nzval[0],size_nzval);
          logfileptr->OFS()<<"RECV "<<src_snode<<std::endl;

          for(int blkidx=0;blkidx<src_snode.NZBlockCnt();++blkidx){
            NZBlockDesc & nzblk_desc = src_snode.GetNZBlockDesc(blkidx);
            Int nRows = src_snode.NRows(blkidx);
            T * nzblk_nzval = src_snode.GetNZval(nzblk_desc.Offset);

            T * dest = fullMatrix.Data();
            for(Int i = 0; i<nRows;++i){
              for(Int j = 0; j<src_snode.Size();++j){
                dest[nzblk_desc.GIndex -1 + i + (src_snode.FirstCol()-1+j)*fullMatrix.m()] = nzblk_nzval[i * src_snode.Size() + j];
              }
            }
          } 
        }
        else{

          Int iLocalI = (I-1) / np +1 ;
          SuperNode<T> & src_snode = *LocalSupernodes_[iLocalI -1];
          for(int blkidx=0;blkidx<src_snode.NZBlockCnt();++blkidx){
            NZBlockDesc & nzblk_desc = src_snode.GetNZBlockDesc(blkidx);
            Int nRows = src_snode.NRows(blkidx);
            T * nzblk_nzval = src_snode.GetNZval(nzblk_desc.Offset);

            T * dest = fullMatrix.Data();
            for(Int i = 0; i<nRows;++i){
              for(Int j = 0; j<src_snode.Size();++j){
                dest[nzblk_desc.GIndex -1 + i + (src_snode.FirstCol()-1+j)*fullMatrix.m()] = nzblk_nzval[i * src_snode.Size() + j];
              }
            }
          }
        }
      }
    }

    
      if(iam==0){
    logfileptr->OFS()<<fullMatrix<<std::endl;
     } 
  }


  template <typename T> void SupernodalMatrix<T>::forward_update(SuperNode<T> * src_contrib,SuperNode<T> * tgt_contrib){


    Int iOwner = Mapping_.Map(src_contrib->Id()-1,src_contrib->Id()-1);
    Int src_ncols = src_contrib->Size();
    Int tgt_ncols = tgt_contrib->Size();

    Int startBlock = (iam==iOwner)?1:0;
    Int tgt_blkidx = -1;

    Int src_blkidx = startBlock;
    while(src_blkidx<src_contrib->NZBlockCnt()){
      NZBlockDesc & src_desc = src_contrib->GetNZBlockDesc(src_blkidx);
      Int src_nrows = src_contrib->NRows(src_blkidx);
  
      if(tgt_blkidx==-1){
          tgt_blkidx = tgt_contrib->FindBlockIdx(src_desc.GIndex);
      }

      NZBlockDesc & tgt_desc = tgt_contrib->GetNZBlockDesc(tgt_blkidx);
      Int tgt_nrows = tgt_contrib->NRows(tgt_blkidx);

      Int src_local_fr = max(tgt_desc.GIndex - src_desc.GIndex,0);
      Int src_lr = src_desc.GIndex+src_nrows-1;

      Int tgt_local_fr = max(src_desc.GIndex - tgt_desc.GIndex,0);
      Int tgt_lr = tgt_desc.GIndex+tgt_nrows-1;
      Int tgt_local_lr = min(src_lr,tgt_lr) - tgt_desc.GIndex;

      T * src = &src_contrib->GetNZval(src_desc.Offset)[src_local_fr*src_ncols];
      T * tgt = &tgt_contrib->GetNZval(tgt_desc.Offset)[tgt_local_fr*tgt_ncols];

      blas::Axpy((tgt_local_lr - tgt_local_fr +1)*src_ncols,
          ONE<T>(),src,1,tgt,1);

      if(src_lr>tgt_lr){
        //the src block hasn't been completely used and is
        // updating some lines in the nz block just below the diagonal block
        //this case happens only locally for the diagonal block
        assert(tgt_blkidx==0);
        //skip to the next tgt nz block
        ++tgt_blkidx;
      }
      else{
        //skip to the next src nz block
        ++src_blkidx;
        if(src_blkidx<src_contrib->NZBlockCnt()){
          NZBlockDesc & next_src_desc = src_contrib->GetNZBlockDesc(src_blkidx);
          tgt_blkidx = tgt_contrib->FindBlockIdx(next_src_desc.GIndex);
        }
      }
    }
  }


  template <typename T> void SupernodalMatrix<T>::back_update(SuperNode<T> * src_contrib,SuperNode<T> * tgt_contrib){
    Int nrhs = tgt_contrib->Size();
    for(Int blkidx = 1; blkidx<tgt_contrib->NZBlockCnt();++blkidx){
      NZBlockDesc & tgt_desc = tgt_contrib->GetNZBlockDesc(blkidx);
      Int tgt_nrows = tgt_contrib->NRows(blkidx);

      Int src_nzblk_idx = src_contrib->FindBlockIdx(tgt_desc.GIndex);
      Int tgt_lr = tgt_desc.GIndex+tgt_nrows-1;
      Int src_lr; 
      do
      {
        NZBlockDesc & src_desc = src_contrib->GetNZBlockDesc(src_nzblk_idx);
        Int src_nrows = src_contrib->NRows(src_nzblk_idx);
        src_lr = src_desc.GIndex+src_nrows-1;

        Int src_local_fr = max(tgt_desc.GIndex - src_desc.GIndex,0);

        Int tgt_local_fr = max(src_desc.GIndex - tgt_desc.GIndex,0);
        Int tgt_local_lr = min(src_lr,tgt_lr) - tgt_desc.GIndex;

        T * src = &src_contrib->GetNZval(src_desc.Offset)[src_local_fr*nrhs];
        T * tgt = &tgt_contrib->GetNZval(tgt_desc.Offset)[tgt_local_fr*nrhs];

        std::copy(src,src+(tgt_local_lr - tgt_local_fr +1)*nrhs,tgt);
        //              lapack::Lacpy('N',nrhs,(tgt_local_lr - tgt_local_fr +1),
        //                  src,nrhs,  
        //                  tgt,nrhs);

        //do other block
        if(tgt_lr>src_lr){
          assert(src_nzblk_idx==0);
          src_nzblk_idx++;
        }

      } while(tgt_lr>src_lr);
    }
  }

  template <typename T> Int SupernodalMatrix<T>::forward_update(NZBlock<T> & dist_nzblk, std::vector<Int> & GlobToLocIndx,SuperNode<T> * contrib,NumMat<T> & B, Int local_blkidx = -1){

////            if(local_blkidx==-1){
////              local_blkidx = GlobToLocIndx[dist_nzblk.GIndex()-1];
////            }
////
////            NZBlock<T> & local_nzblk = contrib->GetNZBlock(local_blkidx);
////            Int src_local_fr = max(local_nzblk.GIndex() - dist_nzblk.GIndex(),0);
////            Int src_lr = dist_nzblk.GIndex()+dist_nzblk.NRows()-1;
////
////            Int tgt_local_fr = max(dist_nzblk.GIndex() - local_nzblk.GIndex(),0);
////            Int tgt_lr = local_nzblk.GIndex()+local_nzblk.NRows()-1;
////            Int tgt_local_lr = min(src_lr,tgt_lr) - local_nzblk.GIndex();
////
//////            if(tgt_local_fr>0){
//////              //copy B in the space preceding the updated block
//////              lapack::Lacpy('N',tgt_local_fr,local_nzblk.NCols(),
//////                  &B(local_nzblk.GIndex()-1,0),B.m(),
//////                  &local_nzblk.Nzval(0,0),local_nzblk.LDA());
//////            }
//////
//////            //copy the contrib
//////            lapack::Lacpy('N',tgt_local_lr - tgt_local_fr +1,dist_nzblk.NCols(),
//////                &dist_nzblk.Nzval(src_local_fr,0),dist_nzblk.LDA(),
//////                &local_nzblk.Nzval(tgt_local_fr,0),local_nzblk.LDA());
//////            if(tgt_local_lr+1<local_nzblk.NRows()){
//////              //copy B in the space following the updated block
//////              lapack::Lacpy('N',local_nzblk.NRows()-(tgt_local_lr+1),local_nzblk.NCols(),
//////                  &B(local_nzblk.GIndex()-1+(tgt_local_lr+1),0),B.m(),
//////                  &local_nzblk.Nzval(tgt_local_lr+1,0),local_nzblk.LDA());
//////            }
////
////
////            //add the contrib
////#ifdef ROW_MAJOR
////            blas::Axpy((tgt_local_lr - tgt_local_fr +1)*dist_nzblk.NCols(),
////                          ONE<T>(),&dist_nzblk.Nzval(src_local_fr,0),1,
////                            &local_nzblk.Nzval(tgt_local_fr,0),1);
////#else
////            for(Int j = 0; j< dist_nzblk.NCols();++j){
////
////              blas::Axpy((tgt_local_lr - tgt_local_fr +1),
////                          ONE<T>(),&dist_nzblk.Nzval(src_local_fr,j),1,
////                             &local_nzblk.Nzval(tgt_local_fr,j),1);
////            }
////#endif
////
////      return local_blkidx;

}



  template <typename T> void SupernodalMatrix<T>::Solve(NumMat<T> * RHS, MPI_Comm & pComm, NumMat<T> * Xptr=NULL){


    NumMat<T> & B = *RHS;
    Int nrhs = B.n();

//    NumMat<T> TMP(B.m(),B.n());
//    NumMat<T> STEPS_SP(B.m(),B.m());
//    SetValue(STEPS_SP,ZERO<T>()); 
//    SetValue(TMP,ZERO<T>()); 

    Int iam,np;
    MPI_Comm_rank(pComm, &iam);
    MPI_Comm_size(pComm, &np);

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

#ifdef _DEBUG_
    logfileptr->OFS()<<"Children vector is"<<children<<std::endl;
#endif


    IntNumVec UpdatesToDo = children;
    std::vector<SuperNode<T> *> contributions(LocalSupernodes_.size());
    std::vector<std::stack<Int> > LocalUpdates(LocalSupernodes_.size());

    //forward-substitution phase
    //Sending contrib up the tree
    //Start from the leaves of the tree

    //This corresponds to the k loop in dtrsm
    for(Int I=1;I<Xsuper_.m();I++){
      Int iOwner = Mapping_.Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam ){
        //Create the blocks of my contrib with the same nz structure as L
        //and the same width as the final solution
        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];

        Int parent = ETree_.PostParent(cur_snode->LastCol()-1);

        SuperNode<T> * contrib = new SuperNode<T>(I,1,nrhs, Size(), cur_snode->NRowsBelowBlock(0) );
        contributions[iLocalI-1] = contrib;


        for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){
          NZBlockDesc & cur_desc = cur_snode->GetNZBlockDesc(blkidx);
          contrib->AddNZBlock(cur_snode->NRows(blkidx),nrhs,cur_desc.GIndex);
//          NZBlock<T> & contrib_nzblk = contrib->GetNZBlock(blkidx);
//          contrib_nzblk.Zero();
//          SetValue(contrib_nzblk.Nzval(),ZERO<T>());
          
        }

        Int oldsize =contrib->StorageSize();
        Int newsize =contrib->Shrink();
//        logfileptr->OFS()<<"new size of contrib is "<<newsize<<" vs "<<oldsize<<std::endl;


        //Do all my updates (Local and remote)
        //Local updates
        while(!LocalUpdates[iLocalI-1].empty()){
          Int contrib_snode_id = LocalUpdates[iLocalI-1].top();
          LocalUpdates[iLocalI-1].pop();

          SuperNode<T> * dist_contrib = contributions[(contrib_snode_id-1) / np];

          logfileptr->OFS()<<"LOCAL Supernode "<<I<<" is updated by contrib of Supernode "<<contrib_snode_id<<std::endl;
          forward_update(dist_contrib,contrib);
            //delete contributions[(contrib_snode_id-1) / np];
          --UpdatesToDo(I-1);
        }

        //do remote updates
        std::vector<char> src_blocks;
        std::vector<T> src_nzval;
        size_t max_bytes;
        Int nz_cnt;
        while(UpdatesToDo(I-1)>0){
          //receive children contrib
          logfileptr->OFS()<<UpdatesToDo(I-1)<<" contribs left"<<endl;




          if(src_blocks.size()==0){
            Int nrows = cur_snode->NRowsBelowBlock(0);
            Int ncols = nrhs;
            nz_cnt = nrows * ncols;

            max_bytes = nrows*sizeof(NZBlockDesc) + sizeof(Int);
            src_blocks.resize(max_bytes);
            src_nzval.resize(nz_cnt);
          }

          

          //MPI_Recv
          MPI_Status recv_status;

          //receive the index array
          MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,MPI_ANY_SOURCE,I,pComm,&recv_status);
          int bytes_received = 0;
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
          //there an aditional integer storing the src id which needs to be removed
          Int src_nzblk_cnt = (bytes_received - sizeof(Int) ) / sizeof(NZBlockDesc);

          //receive the nzval array
          MPI_Recv(&src_nzval[0],nz_cnt*sizeof(T),MPI_BYTE,MPI_ANY_SOURCE,I,pComm,&recv_status);
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
          Int src_nzval_cnt = bytes_received / sizeof(T);

          Int src_snode_id = *(Int*)&src_blocks[0];
          NZBlockDesc * src_blocks_ptr = 
                    reinterpret_cast<NZBlockDesc*>(&src_blocks[sizeof(Int)]);
          T * src_nzval_ptr = &src_nzval[0];

          //Create the dummy supernode for that data
          SuperNode<T> dist_contrib(src_snode_id,1,nrhs, Size(), src_blocks_ptr, src_nzblk_cnt, src_nzval_ptr, src_nzval_cnt);

#ifdef _DEBUG_
          logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib.Id()<<std::endl;
#endif

          forward_update(&dist_contrib,contrib);

          --UpdatesToDo(I-1);

        }


        if(UpdatesToDo(I-1)==0){
            //This corresponds to the i loop in dtrsm
            for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){

              NZBlockDesc & cur_desc = contrib->GetNZBlockDesc(blkidx);
              NZBlockDesc & chol_desc = cur_snode->GetNZBlockDesc(blkidx);
              NZBlockDesc & diag_desc = contrib->GetNZBlockDesc(0);
              
              Int cur_nrows = contrib->NRows(blkidx);
              Int chol_nrows = cur_snode->NRows(blkidx);
              Int diag_nrows = contrib->NRows(0);

              T * cur_nzval = contrib->GetNZval(cur_desc.Offset);
              T * chol_nzval = cur_snode->GetNZval(chol_desc.Offset);
              T * diag_nzval = contrib->GetNZval(diag_desc.Offset);

              //compute my contribution
              //Handle the diagonal block
              if(blkidx==0){
                //TODO That's where we can use the selective inversion
                //if we are processing the "pivot" block
                for(Int kk = 0; kk<cur_snode->Size(); ++kk){
                  for(Int j = 0; j<nrhs;++j){
                    diag_nzval[kk*nrhs+j] = (B(diag_desc.GIndex-1,j) + diag_nzval[kk*nrhs+j]) / chol_nzval[kk*cur_snode->Size()+kk];
                    for(Int i = kk+1; i<cur_nrows;++i){
                      diag_nzval[i*nrhs+j] += -diag_nzval[kk*nrhs+j]*chol_nzval[i*cur_snode->Size()+kk];
                    }
                  }
                }
              }
              else{
                for(Int kk = 0; kk<cur_snode->Size(); ++kk){
                  for(Int j = 0; j<nrhs;++j){
                    for(Int i = 0; i<cur_nrows;++i){
                      cur_nzval[i*nrhs+j] += -diag_nzval[kk*nrhs+j]*chol_nzval[i*cur_snode->Size()+kk];
                    }
                  }
                }
              }
            }
          

          //send to my parent
          if(parent!=0){
            Int parent_snode_id = SupMembership_[parent-1];

            Int iTarget = Mapping_.Map(parent_snode_id-1,parent_snode_id-1);

            if(iTarget!=iam){
#ifdef _DEBUG_
              logfileptr->OFS()<<"Remote Supernode "<<parent_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif



              std::vector<char> * pNewDesc = new std::vector<char>();

              Int tgt_first_col = Xsuper_(parent_snode_id-1);
              Int tgt_last_col = Xsuper_(parent_snode_id)-1;

              //Send
              Int src_nzblk_idx = 1;
              NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);

              Int src_first_row = pivot_desc.GIndex;
              Int local_first_row = 0;
              Int nzblk_cnt = contrib->NZBlockCnt()-src_nzblk_idx;
              pNewDesc->resize(sizeof(Int) + nzblk_cnt*sizeof(NZBlockDesc));

              char * send_ptr = &(*pNewDesc)[0];

              Int * id_ptr = reinterpret_cast<Int *>(&(*pNewDesc)[0]);
              NZBlockDesc * block_desc_ptr = 
                reinterpret_cast<NZBlockDesc *>(&(*pNewDesc)[sizeof(Int)]);

              *id_ptr = contrib->Id();
              //copy the block descriptors
              std::copy(&pivot_desc,&pivot_desc+nzblk_cnt, block_desc_ptr);
              //change the first one
              block_desc_ptr->Offset += (src_first_row - block_desc_ptr->GIndex)*nrhs;
              block_desc_ptr->GIndex = src_first_row;

              //send the block descriptors
              MPI_Send(send_ptr,nzblk_cnt*sizeof(NZBlockDesc) + sizeof(Int),
                  MPI_BYTE,iTarget,parent_snode_id,pComm);

              T * nzval_ptr = contrib->GetNZval(pivot_desc.Offset
                  +local_first_row*nrhs);

              Int nz_cnt = (contrib->NRowsBelowBlock(src_nzblk_idx)
                  - local_first_row )*nrhs;

              MPI_Send(nzval_ptr,nz_cnt*sizeof(T),
                  MPI_BYTE,iTarget,parent_snode_id,pComm);

              delete pNewDesc;

#ifdef _DEBUG_            
              logfileptr->OFS()<<"     Send contribution "<<I<<" to Supernode "<<parent_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
              logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
#endif
            }
            else{
#ifdef _DEBUG_
              logfileptr->OFS()<<"Local Supernode "<<parent_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif
              Int iLocalJ = (parent_snode_id-1) / np +1 ;
              LocalUpdates[iLocalJ-1].push(I);
            }
          }
        }
      }
    }

//////    logfileptr->OFS()<<"Solution after each step is "<<STEPS_SP<<std::endl;
//////    logfileptr->OFS()<<"Solution is after forward substitution is "<<TMP<<std::endl;
////
//////    SetValue(STEPS_SP,ZERO<T>());
    //Back-substitution phase

    //start from the root of the tree
    for(Int I=Xsuper_.m()-1;I>=1;--I){
      Int iOwner = Mapping_.Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam ){
        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];
        SuperNode<T> * contrib = contributions[iLocalI-1];

        Int parent = ETree_.PostParent(cur_snode->LastCol()-1);
        //Extend the contribution.

        std::vector<Int> isBlockUpdated(contrib->NZBlockCnt(),0);

        if(parent!=0){
          Int parent_snode_id = SupMembership_[parent-1];
          Int iTarget = Mapping_.Map(parent_snode_id-1,parent_snode_id-1);
          //Do all my updates (Local and remote)
          //Local updates
          SuperNode<T> * dist_contrib;
          if(!LocalUpdates[iLocalI-1].empty()){
            Int contrib_snode_id = LocalUpdates[iLocalI-1].top();
            LocalUpdates[iLocalI-1].pop();

            dist_contrib = contributions[(contrib_snode_id-1) / np];
            back_update(dist_contrib,contrib);
          }
          else{

            //Receive parent contrib
            std::vector<char> src_blocks;
            std::vector<T> src_nzval;


            //MPI_Recv
            MPI_Status recv_status;

            //Receive the size of the blocks array
            Int max_bytes = 0;
            MPI_Recv(&max_bytes,sizeof(Int),MPI_BYTE,MPI_ANY_SOURCE,I,pComm,&recv_status);
            src_blocks.resize(max_bytes);


            //receive the index array
            MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,MPI_ANY_SOURCE,I,pComm,&recv_status);
            //there an aditional integer storing the src id which needs to be removed
            Int src_nzblk_cnt = (max_bytes - sizeof(Int) ) / sizeof(NZBlockDesc);



            //Receive the size of the nzval array
            Int nz_cnt = 0;
            MPI_Recv(&nz_cnt,sizeof(Int),MPI_BYTE,MPI_ANY_SOURCE,I,pComm,&recv_status);
            src_nzval.resize(nz_cnt);

            //receive the nzval array
            MPI_Recv(&src_nzval[0],nz_cnt*sizeof(T),MPI_BYTE,MPI_ANY_SOURCE,I,pComm,&recv_status);

            Int src_snode_id = *(Int*)&src_blocks[0];
            NZBlockDesc * src_blocks_ptr = 
              reinterpret_cast<NZBlockDesc*>(&src_blocks[sizeof(Int)]);
            T * src_nzval_ptr = &src_nzval[0];

            //Create the dummy supernode for that data
            dist_contrib = new SuperNode<T>(src_snode_id,1,nrhs, Size(), src_blocks_ptr, src_nzblk_cnt, src_nzval_ptr, nz_cnt);

#ifdef _DEBUG_
            logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib->Id()<<std::endl;
#endif

            back_update(dist_contrib,contrib);
            delete dist_contrib;
          }
        }

        //now compute MY contribution

        NZBlockDesc & diag_desc = cur_snode->GetNZBlockDesc(0);
        NZBlockDesc & tgt_desc = contrib->GetNZBlockDesc(0);

        T* diag_nzval = cur_snode->GetNZval(diag_desc.Offset);
        T* tgt_nzval = contrib->GetNZval(tgt_desc.Offset);

        for(Int j = 0; j<nrhs;++j){
          for(Int ii = cur_snode->Size()-1; ii>=0; --ii){
            T temp = tgt_nzval[ii*nrhs+j];

            //This corresponds to the k loop in dtrsm
            for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){
              NZBlockDesc & chol_desc = cur_snode->GetNZBlockDesc(blkidx);
              Int chol_nrows = cur_snode->NRows(blkidx);

              Int src_blkidx = contrib->FindBlockIdx(chol_desc.GIndex);
              NZBlockDesc & cur_desc = contrib->GetNZBlockDesc(src_blkidx);
              Int cur_nrows = contrib->NRows(src_blkidx);

              T* chol_nzval = cur_snode->GetNZval(chol_desc.Offset);
              T* cur_nzval = contrib->GetNZval(cur_desc.Offset);

              for(Int kk = 0; kk< chol_nrows; ++kk){
                if(chol_desc.GIndex+kk>cur_snode->FirstCol()+ii){
                  Int src_row = chol_desc.GIndex - cur_desc.GIndex +kk;
                  if(src_row< cur_nrows){
                    temp += -chol_nzval[kk*cur_snode->Size()+ii]*cur_nzval[src_row*nrhs+j];

                    //                      if(j==0){
                    //                        STEPS_SP(ii,I-1)+= -chol_nzblk.Nzval(kk,ii)*cur_nzblk->Nzval(src_row,j);
                    //                      }
                  }
                  //                    else{
                  //                      src_blkidx++;
                  //                      cur_nzblk = &contrib->GetNZBlock(src_blkidx);
                  //                      src_row = chol_nzblk.GIndex() - cur_nzblk->GIndex() +kk;
                  //                      temp += -chol_nzblk.Nzval(kk,ii)*cur_nzblk->Nzval(src_row,j);
                  //
                  //                      if(j==0){
                  //                        STEPS_SP(ii,I-1+ii)+= -chol_nzblk.Nzval(kk,ii)*cur_nzblk->Nzval(src_row,j);
                  //                      }
                  //                    }
                }
              }
            }

            temp = temp / diag_nzval[ii*cur_snode->Size()+ii];
            tgt_nzval[ii*nrhs+j] = temp;
          }
        }



//        //TODO DO BETTER THAN THIS
//        //copy the updated block in B
//        lapack::Lacpy('N',tgt_nzblk.NRows(),tgt_nzblk.NCols(),
//            &tgt_nzblk.Nzval(0,0),tgt_nzblk.LDA(),
//            &B(tgt_nzblk.GIndex()-1,0),B.m());




        //send to my children
        Int colIdx = cur_snode->FirstCol()-1;
        if(colIdx>0){
          Int children_found = 0;
          while(children_found<children(I-1)){
            Int child_snode_id = SupMembership_[colIdx-1];
            if(ETree_.PostParent(child_snode_id-1)==cur_snode->FirstCol()){
              Int iTarget = Mapping_.Map(child_snode_id-1,child_snode_id-1);

              if(iTarget!=iam){



#ifdef _DEBUG_
                logfileptr->OFS()<<"Remote Supernode "<<child_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif



                std::vector<char> * pNewDesc = new std::vector<char>();

                Int tgt_first_col = Xsuper_(child_snode_id-1);
                Int tgt_last_col = Xsuper_(child_snode_id)-1;

                //Send
                Int src_nzblk_idx = 0;
                NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);

                Int src_first_row = pivot_desc.GIndex;
                Int local_first_row = 0;
                Int nzblk_cnt = contrib->NZBlockCnt()-src_nzblk_idx;
                pNewDesc->resize(sizeof(Int) + nzblk_cnt*sizeof(NZBlockDesc));

                char * send_ptr = &(*pNewDesc)[0];

                Int * id_ptr = reinterpret_cast<Int *>(&(*pNewDesc)[0]);
                NZBlockDesc * block_desc_ptr = 
                  reinterpret_cast<NZBlockDesc *>(&(*pNewDesc)[sizeof(Int)]);

                *id_ptr = contrib->Id();
                //copy the block descriptors
                std::copy(&pivot_desc,&pivot_desc+nzblk_cnt, block_desc_ptr);
                //change the first one
                block_desc_ptr->Offset += (src_first_row - block_desc_ptr->GIndex)*nrhs;
                block_desc_ptr->GIndex = src_first_row;

                Int bytes_size = nzblk_cnt*sizeof(NZBlockDesc)+sizeof(Int);
                MPI_Send(&bytes_size,sizeof(bytes_size),MPI_BYTE,iTarget,child_snode_id,pComm);

                //send the block descriptors
                MPI_Send(send_ptr,bytes_size, MPI_BYTE,iTarget,child_snode_id,pComm);

                T * nzval_ptr = contrib->GetNZval(pivot_desc.Offset
                    +local_first_row*nrhs);

                Int nz_cnt = (contrib->NRowsBelowBlock(src_nzblk_idx)
                    - local_first_row )*nrhs;

                MPI_Send(&nz_cnt,sizeof(nz_cnt),MPI_BYTE,iTarget,child_snode_id,pComm);

                MPI_Send(nzval_ptr,nz_cnt*sizeof(T),MPI_BYTE,iTarget,child_snode_id,pComm);

                delete pNewDesc;

#ifdef _DEBUG_            
                logfileptr->OFS()<<"     Send contribution "<<I<<" to Supernode "<<child_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
                logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
#endif
              }
              else{
                Int iLocalJ = (child_snode_id-1) / np +1 ;
                LocalUpdates[iLocalJ-1].push(I);
              }
              children_found++;
            }
            //last column of the prev supernode
            colIdx = Xsuper_[child_snode_id-1]-1;
            if(colIdx==0){
              break;
            }
          }
        }


      }
    }



    //Gather B from everybody
    std::vector<T> tmp_nzval;
    for(Int I=1; I<Xsuper_.m();++I){
      Int iOwner = Mapping_.Map(I-1,I-1);
      T * data;
      Int snode_size = Xsuper_[I] - Xsuper_[I-1];
      Int nzcnt = snode_size * nrhs;
      tmp_nzval.resize(nzcnt);

      if( iOwner == iam ){
        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> * contrib = contributions[iLocalI-1];
        data = contrib->GetNZval(0);
      }
      else{
        data = &tmp_nzval[0];
      }

      MPI_Bcast(data,nzcnt*sizeof(T),MPI_BYTE,iOwner,pComm);

      for(Int i = 0; i<snode_size; ++i){ 
        for(Int j = 0; j<nrhs; ++j){
          B(Xsuper_[I-1] -1 + i, j) = data[i*nrhs + j];
        }
      }
    }

////
////
//////    logfileptr->OFS()<<"Solution after each step is "<<STEPS_SP<<std::endl;
////    //Print B
    logfileptr->OFS()<<"Solution is "<<B<<std::endl;
////


















////    NumMat<T> & B = *RHS;
////
//////    NumMat<T> TMP(B.m(),B.n());
//////    NumMat<T> STEPS_SP(B.m(),B.m());
//////    SetValue(STEPS_SP,ZERO<T>()); 
//////    SetValue(TMP,ZERO<T>()); 
////
////    Int iam,np;
////    MPI_Comm_rank(pComm, &iam);
////    MPI_Comm_size(pComm, &np);
////
////    IntNumVec children(Xsuper_.m());
////    SetValue(children,0);
////
////    for(Int I=1;I<Xsuper_.m()-1;I++){
////      Int fc = Xsuper_(I-1);
////      Int lc = Xsuper_(I)-1;
////
////      Int parent = ETree_.PostParent(lc-1);
////      if(parent!=0){
////        ++children(SupMembership_(parent-1)-1);
////      }
////    }
////
////    logfileptr->OFS()<<"Children vector is"<<children<<std::endl;
////
////
////    IntNumVec UpdatesToDo = children;
////    std::vector<SuperNode<T> *> contributions(LocalSupernodes_.size());
////    std::vector<std::stack<Int> > LocalUpdates(LocalSupernodes_.size());
////
////    //forward-substitution phase
////    //Sending contrib up the tree
////    //Start from the leaves of the tree
////
////    //This corresponds to the k loop in dtrsm
////    for(Int I=1;I<Xsuper_.m();I++){
////      Int iOwner = Mapping_.Map(I-1,I-1);
////      //If I own the column, factor it
////      if( iOwner == iam ){
////        //Create the blocks of my contrib with the same nz structure as L
////        //and the same width as the final solution
////        Int iLocalI = (I-1) / np +1 ;
////        SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];
////
////        Int parent = ETree_.PostParent(cur_snode->LastCol()-1);
////
////        SuperNode<T> * contrib = new SuperNode<T>(I,1,B.n(),cur_snode->NRows());
////        contributions[iLocalI-1] = contrib;
////
////        std::vector<Int> GlobToLocIndx(Size(),-1);
////
////        for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){
////          NZBlock<T> & cur_nzblk = cur_snode->GetNZBlock(blkidx);
////          contrib->AddNZBlock(cur_nzblk.NRows(),B.n(),cur_nzblk.GIndex());
////          NZBlock<T> & contrib_nzblk = contrib->GetNZBlock(blkidx);
////          contrib_nzblk.Zero();
//////          SetValue(contrib_nzblk.Nzval(),ZERO<T>());
////          
////          //fill the indirect index array
////          for(Int gi = cur_nzblk.GIndex(); gi<cur_nzblk.GIndex()+cur_nzblk.NRows();++gi){
////            GlobToLocIndx[gi-1] = blkidx;
////          }
////        }
////
////        Int oldsize =contrib->StorageSize();
////        Int newsize =contrib->Shrink();
//////        logfileptr->OFS()<<"new size of contrib is "<<newsize<<" vs "<<oldsize<<std::endl;
////
////
////        //Do all my updates (Local and remote)
////        //Local updates
////        while(!LocalUpdates[iLocalI-1].empty()){
////          Int contrib_snode_id = LocalUpdates[iLocalI-1].top();
////          LocalUpdates[iLocalI-1].pop();
////
////          SuperNode<T> & dist_contrib = *contributions[(contrib_snode_id-1) / np];
////
////          logfileptr->OFS()<<"LOCAL Supernode "<<I<<" is updated by contrib of Supernode "<<contrib_snode_id<<std::endl;
////          //skip the "diagonal" block
////          for(Int blkidx = 1; blkidx<dist_contrib.NZBlockCnt();++blkidx){
////            NZBlock<T> & dist_nzblk = dist_contrib.GetNZBlock(blkidx);
////            Int local_blkidx = forward_update(dist_nzblk,GlobToLocIndx,contrib,B);
////
////            NZBlock<T> & local_nzblk = contrib->GetNZBlock(local_blkidx);
////            Int src_lr = dist_nzblk.GIndex()+dist_nzblk.NRows()-1;
////            Int tgt_lr = local_nzblk.GIndex()+local_nzblk.NRows()-1;
////            //this case happens only locally for the diagonal block
////            if(src_lr>tgt_lr){
////              assert(local_blkidx==0);
////              ++local_blkidx;
////
////              local_blkidx = forward_update(dist_nzblk,GlobToLocIndx,contrib,B,local_blkidx);
////              //there won't be any other update to that block
////            }
////          }
////            //delete contributions[(contrib_snode_id-1) / np];
////          --UpdatesToDo(I-1);
////        }
////
////        //do remote updates
////        std::vector<char> src_nzblocks;
////        size_t num_bytes;
////        while(UpdatesToDo(I-1)>0){
////          //receive children contrib
////          logfileptr->OFS()<<UpdatesToDo(I-1)<<" contribs left"<<endl;
////
////          if(src_nzblocks.size()==0){
////            num_bytes = sizeof(Int)+sizeof(NZBlock<T>);
////            for(Int blkidx=0;blkidx<cur_snode->NZBlockCnt();++blkidx){
////              num_bytes += cur_snode->GetNZBlock(blkidx).NRows()*NZBLOCK_ROW_SIZE<T>(B.n());
////            }
////            logfileptr->OFS()<<"We allocate a buffer of size "<<num_bytes<<std::endl;
////            src_nzblocks.resize(num_bytes);
////          }
////
////
////          //MPI_Recv
////          MPI_Status recv_status;
////
////          MPI_Probe(MPI_ANY_SOURCE,I,pComm,&recv_status);
////          int bytes_to_receive = 0;
////          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_to_receive);
////
////          assert(src_nzblocks.size()>=bytes_to_receive);
////
////
////          char * recv_buf = &src_nzblocks[0]+ max(sizeof(Int),NZBLOCK_OBJ_SIZE<T>())-min(sizeof(Int),NZBLOCK_OBJ_SIZE<T>());
////          MPI_Recv(recv_buf,&*src_nzblocks.end()-recv_buf,MPI_BYTE,MPI_ANY_SOURCE,I,pComm,&recv_status);
////          Int src_snode_id = *(Int*)recv_buf;
////
////          //Resize the buffer to the actual number of bytes received
////          int bytes_received = 0;
////          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
////          src_nzblocks.resize(recv_buf+bytes_received - &src_nzblocks[0]);
////
////
////          //Create the dummy supernode for that data
////          SuperNode<T> dist_contrib(src_snode_id,1,B.n(),&src_nzblocks);
////
////          logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib.Id()<<std::endl;
////
////          for(Int blkidx = 0; blkidx<dist_contrib.NZBlockCnt();++blkidx){
////            NZBlock<T> & dist_nzblk = dist_contrib.GetNZBlock(blkidx);
////            Int local_blkidx = forward_update(dist_nzblk,GlobToLocIndx,contrib,B);
////          }
////          --UpdatesToDo(I-1);
////
////          src_nzblocks.resize(num_bytes);
////        }
////
////
////        if(UpdatesToDo(I-1)==0){
////            //This corresponds to the i loop in dtrsm
////            for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){
////              NZBlock<T> & cur_nzblk = contrib->GetNZBlock(blkidx);
////              NZBlock<T> & chol_nzblk = cur_snode->GetNZBlock(blkidx);
////            NZBlock<T> & diag_nzblk = contrib->GetNZBlock(0);
////              //compute my contribution
////              //Handle the diagonal block
////              if(blkidx==0){
////                //TODO That's where we can use the selective inversion
////                //if we are processing the "pivot" block
////                for(Int kk = 0; kk<cur_snode->Size(); ++kk){
////                  for(Int j = 0; j<diag_nzblk.NCols();++j){
////                    diag_nzblk.Nzval(kk,j) = (B(diag_nzblk.GIndex()-1,j) + diag_nzblk.Nzval(kk,j)) / chol_nzblk.Nzval(kk,kk);
////                    for(Int i = kk+1; i<cur_nzblk.NRows();++i){
////                      diag_nzblk.Nzval(i,j) += -diag_nzblk.Nzval(kk,j)*chol_nzblk.Nzval(i,kk);
////                    }
////                  }
////                }
////              }
////              else{
////                for(Int kk = 0; kk<cur_snode->Size(); ++kk){
////                  for(Int j = 0; j<cur_nzblk.NCols();++j){
////                    for(Int i = 0; i<cur_nzblk.NRows();++i){
////                      cur_nzblk.Nzval(i,j) += -diag_nzblk.Nzval(kk,j)*chol_nzblk.Nzval(i,kk);
////                    }
////                  }
////                }
////              }
////            }
////          
////
////
//////          {
//////              NZBlock<T> & cur_nzblk = contrib->GetNZBlock(0);
//////                lapack::Lacpy('N',cur_nzblk.NRows(),cur_nzblk.NCols(),
//////                    cur_nzblk.Nzval(),cur_nzblk.LDA(),
//////                    &TMP(cur_nzblk.GIndex()-1,0),TMP.m()
//////                    );
//////          }
////
////
//////          for(Int blkidx=0; blkidx<contrib->NZBlockCnt();++blkidx)
//////          {
//////              NZBlock<T> & cur_nzblk = contrib->GetNZBlock(blkidx);
//////                lapack::Lacpy('N',cur_nzblk.NRows(),1,
//////                    &cur_nzblk.Nzval(0,0),cur_nzblk.LDA(),
//////                    &STEPS_SP(cur_nzblk.GIndex()-1,I-1),STEPS_SP.m()
//////                    );
//////
////// 
//////              Int start = blkidx==0?1:0;
//////              for(Int i =start; i<cur_nzblk.NRows();++i){
//////                //Add B to it
//////                STEPS_SP(cur_nzblk.GIndex()-1+i,I-1)*=-1;
//////              }
//////              
//////          }
////
////
////
////          //send to my parent
////          if(parent!=0){
////            Int parent_snode_id = SupMembership_[parent-1];
////
////            Int iTarget = Mapping_.Map(parent_snode_id-1,parent_snode_id-1);
////
////            if(iTarget!=iam){
////              logfileptr->OFS()<<"Supernode "<<parent_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
////              Int J = parent_snode_id;
////
////              Int tgt_first_col = Xsuper_(J-1);
////              Int tgt_last_col = Xsuper_(J)-1;
////
////
////#ifdef ROW_MAJOR
////#else
////#ifdef SEND_WHOLE_BLOCK
////#else
////              char * start_buf = reinterpret_cast<char*>(&contrib->GetNZBlock(1));
////              size_t num_bytes = contrib->End() - start_buf;
////
////              MPI_Datatype type;
////              int lens[2];
////              MPI_Aint disps[2];
////              MPI_Datatype types[2];
////              Int err = 0;
////
////
////              /* define a struct that describes all our data */
////              lens[0] = sizeof(Int);
////              lens[1] = num_bytes;
////              MPI_Address(&I, &disps[0]);
////              MPI_Address(start_buf, &disps[1]);
////              types[0] = MPI_BYTE;
////              types[1] = MPI_BYTE;
////              MPI_Type_struct(2, lens, disps, types, &type);
////              MPI_Type_commit(&type);
////
////              MPI_Send(MPI_BOTTOM,1,type,iTarget,parent_snode_id,pComm);
////
////              logfileptr->OFS()<<"     SEND Supernode "<<I<<" to Supernode "<<J<<" on P"<<iTarget<<" ("<<lens[0] + lens[1] <<" bytes)"<<std::endl;
////              MPI_Type_free(&type);
////#endif // end ifdef SEND_WHOLE_BLOCK
////
////#endif // end ifdef ROW_MAJOR
////
////              //delete it
////              //delete contributions[iLocalI-1];
////
////            }
////            else{
////              Int iLocalJ = (parent_snode_id-1) / np +1 ;
////              LocalUpdates[iLocalJ-1].push(I);
////            }
////          }
////        }
////      }
////    }
////
//////    logfileptr->OFS()<<"Solution after each step is "<<STEPS_SP<<std::endl;
//////    logfileptr->OFS()<<"Solution is after forward substitution is "<<TMP<<std::endl;
////
//////    SetValue(STEPS_SP,ZERO<T>());
////    //Back-substitution phase
////
////    //start from the root of the tree
////    for(Int I=Xsuper_.m()-1;I>=1;--I){
////      Int iOwner = Mapping_.Map(I-1,I-1);
////      //If I own the column, factor it
////      if( iOwner == iam ){
////        Int iLocalI = (I-1) / np +1 ;
////        SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];
////        SuperNode<T> * contrib = contributions[iLocalI-1];
////
////        Int parent = ETree_.PostParent(cur_snode->LastCol()-1);
////        //Extend the contribution.
////
////
////        std::vector<Int> GlobToLocIndx(Size(),-1);
////
////
////        for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){
////          NZBlock<T> & cur_nzblk = cur_snode->GetNZBlock(blkidx);
////          //fill the indirect index array
////          for(Int gi = cur_nzblk.GIndex(); gi<cur_nzblk.GIndex()+cur_nzblk.NRows();++gi){
////            GlobToLocIndx[gi-1] = blkidx;
////          }
////        }
////
////        std::vector<Int> isBlockUpdated(contrib->NZBlockCnt(),0);
////
////        if(parent!=0){
////          Int parent_snode_id = SupMembership_[parent-1];
////          Int iTarget = Mapping_.Map(parent_snode_id-1,parent_snode_id-1);
////          //Do all my updates (Local and remote)
////          //Local updates
////          SuperNode<T> * dist_contrib;
////          if(!LocalUpdates[iLocalI-1].empty()){
////            Int contrib_snode_id = LocalUpdates[iLocalI-1].top();
////            LocalUpdates[iLocalI-1].pop();
////
////            dist_contrib = contributions[(contrib_snode_id-1) / np];
////          }
////          else{
////
////            //do remote updates
////            std::vector<char> src_nzblocks;
////            size_t num_bytes;
////            //receive parent contrib
////
////            //MPI_Recv
////            MPI_Status recv_status;
////
////
////            MPI_Probe(iTarget,I,pComm,&recv_status);
////            int bytes_to_receive = 0;
////            MPI_Get_count(&recv_status, MPI_BYTE, &bytes_to_receive);
////            src_nzblocks.resize(bytes_to_receive+NZBLOCK_OBJ_SIZE<T>());
////
////
////            char * recv_buf = &src_nzblocks[0]+ max(sizeof(Int),NZBLOCK_OBJ_SIZE<T>())-min(sizeof(Int),NZBLOCK_OBJ_SIZE<T>());
////            MPI_Recv(recv_buf,&*src_nzblocks.end()-recv_buf,MPI_BYTE,iTarget,I,pComm,&recv_status);
////            Int src_snode_id = *(Int*)recv_buf;
////
////            //Resize the buffer to the actual number of bytes received
////            //            int bytes_received = 0;
////            //            MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
////            //            src_nzblocks.resize(recv_buf+bytes_received - &src_nzblocks[0]);
////
////            //Create the dummy supernode for that data
////            dist_contrib = new SuperNode<T>(src_snode_id,1,B.n(),&src_nzblocks);
////
////            logfileptr->OFS()<<"RECV contrib of Supernode "<<dist_contrib->Id()<<std::endl;
////          }
////
////
////          for(Int blkidx = 1; blkidx<contrib->NZBlockCnt();++blkidx){
////            NZBlock<T> & tgt_nzblk = contrib->GetNZBlock(blkidx);
////
////
////            //TODO this could be replaced if I had the GlobIndex array of the parent
////            Int src_nzblk_idx = dist_contrib->FindBlockIdx(tgt_nzblk.GIndex());
////            Int tgt_lr = tgt_nzblk.GIndex()+tgt_nzblk.NRows()-1;
////            Int src_lr; 
////            do
////            {
////              NZBlock<T> * src_nzblk = &dist_contrib->GetNZBlock(src_nzblk_idx);
////              src_lr = src_nzblk->GIndex()+src_nzblk->NRows()-1;
////
////
////              Int src_local_fr = max(tgt_nzblk.GIndex() - src_nzblk->GIndex(),0);
////
////              Int tgt_local_fr = max(src_nzblk->GIndex() - tgt_nzblk.GIndex(),0);
////              Int tgt_local_lr = min(src_lr,tgt_lr) - tgt_nzblk.GIndex();
////
////
////              lapack::Lacpy('N',(tgt_local_lr - tgt_local_fr +1),src_nzblk->NCols(),
////                  &src_nzblk->Nzval(src_local_fr,0),src_nzblk->LDA(),  
////                  &tgt_nzblk.Nzval(tgt_local_fr,0),tgt_nzblk.LDA());
////
////              //do other block
////              if(tgt_lr>src_lr){
////                assert(src_nzblk_idx==0);
////                src_nzblk_idx++;
////              }
////
////            } while(tgt_lr>src_lr);
////          }
////
////          if( iTarget != iam ){
////            delete dist_contrib;
////          }
////
////
////
////        }
////
////        //now compute MY contribution
////
////          NZBlock<T> & diag_nzblk = cur_snode->GetNZBlock(0);
////          NZBlock<T> & tgt_nzblk = contrib->GetNZBlock(0);
////          for(Int j = 0; j<tgt_nzblk.NCols();++j){
////            for(Int ii = cur_snode->Size()-1; ii>=0; --ii){
////              T temp = tgt_nzblk.Nzval(ii,j);
////
////              //This corresponds to the k loop in dtrsm
////              for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){
////                NZBlock<T> & chol_nzblk = cur_snode->GetNZBlock(blkidx);
////
////                Int src_blkidx = contrib->FindBlockIdx(chol_nzblk.GIndex());
////                NZBlock<T> * cur_nzblk = &contrib->GetNZBlock(src_blkidx);
////
////                for(Int kk = 0; kk< chol_nzblk.NRows(); ++kk){
////                  if(chol_nzblk.GIndex()+kk>cur_snode->FirstCol()+ii){
////                    Int src_row = chol_nzblk.GIndex() - cur_nzblk->GIndex() +kk;
////                    if(src_row< cur_nzblk->NRows()){
////                      temp += -chol_nzblk.Nzval(kk,ii)*cur_nzblk->Nzval(src_row,j);
////
//////                      if(j==0){
//////                        STEPS_SP(ii,I-1)+= -chol_nzblk.Nzval(kk,ii)*cur_nzblk->Nzval(src_row,j);
//////                      }
////                    }
//////                    else{
//////                      src_blkidx++;
//////                      cur_nzblk = &contrib->GetNZBlock(src_blkidx);
//////                      src_row = chol_nzblk.GIndex() - cur_nzblk->GIndex() +kk;
//////                      temp += -chol_nzblk.Nzval(kk,ii)*cur_nzblk->Nzval(src_row,j);
//////
//////                      if(j==0){
//////                        STEPS_SP(ii,I-1+ii)+= -chol_nzblk.Nzval(kk,ii)*cur_nzblk->Nzval(src_row,j);
//////                      }
//////                    }
////                  }
////                }
////              }
////
////              temp = temp / diag_nzblk.Nzval(ii,ii);    
////              tgt_nzblk.Nzval(ii,j) = temp;
////            }
////          }
////
////
////
////          //TODO DO BETTER THAN THIS
////          //copy the updated block in B
////              lapack::Lacpy('N',tgt_nzblk.NRows(),tgt_nzblk.NCols(),
////                  &tgt_nzblk.Nzval(0,0),tgt_nzblk.LDA(),
////                  &B(tgt_nzblk.GIndex()-1,0),B.m());
////
////
////
////
////          //send to my children
////          Int colIdx = cur_snode->FirstCol()-1;
////          if(colIdx>0){
////            Int children_found = 0;
////            while(children_found<children(I-1)){
////                Int nodeIdx = SupMembership_[colIdx-1];
////              if(ETree_.PostParent(colIdx-1)==cur_snode->FirstCol()){
////                Int iTarget = Mapping_.Map(nodeIdx-1,nodeIdx-1);
////
////                if(iTarget!=iam){
////                  logfileptr->OFS()<<"Supernode "<<nodeIdx<<" gets the contribution of Supernode "<<I<<std::endl;
////                  Int J = nodeIdx;
////
////                  Int tgt_first_col = Xsuper_(J-1);
////                  Int tgt_last_col = Xsuper_(J)-1;
////
////
////#ifdef ROW_MAJOR
////#else
////#ifdef SEND_WHOLE_BLOCK
////#else
////
////                  char * start_buf = reinterpret_cast<char*>(&contrib->GetNZBlock(0));
////                  size_t num_bytes = contrib->End() - start_buf;
////
////
////                  MPI_Datatype type;
////                  int lens[2];
////                  MPI_Aint disps[2];
////                  MPI_Datatype types[2];
////                  Int err = 0;
////
////
////                  lens[0] = sizeof(Int);
////                  lens[1] = num_bytes;
////                  MPI_Address(&I, &disps[0]);
////                  MPI_Address(start_buf, &disps[1]);
////                  types[0] = MPI_BYTE;
////                  types[1] = MPI_BYTE;
////                  MPI_Type_struct(2, lens, disps, types, &type);
////                  MPI_Type_commit(&type);
////
////                  MPI_Send(MPI_BOTTOM,1,type,iTarget,nodeIdx,pComm);
////
////                  logfileptr->OFS()<<"     SEND Supernode "<<I<<" to Supernode "<<J<<" on P"<<iTarget<<" ("<<lens[0] + lens[1] <<" bytes)"<<std::endl;
////                  MPI_Type_free(&type);
////#endif // end ifdef SEND_WHOLE_BLOCK
////#endif // end ifdef ROW_MAJOR
////
////
////                }
////                else{
////                  Int iLocalJ = (nodeIdx-1) / np +1 ;
////                  LocalUpdates[iLocalJ-1].push(I);
////                }
////                children_found++;
////              }
////              //last column of the prev supernode
////              colIdx = Xsuper_[nodeIdx-1]-1;
////              if(colIdx==0){
////                break;
////              }
////            }
////          }
////
////          
////      }
////    }
////
////
//////    //Gather B from everybody
//////    for(Int I=Xsuper_.m()-1;I>=1;--I){
//////      Int iOwner = Mapping_.Map(I-1,I-1);
//////      //If I own the column, factor it
//////      if( iOwner == iam ){
//////        Int iLocalI = (I-1) / np +1 ;
//////        SuperNode<T> * contrib = contributions[iLocalI-1];
//////
//////
//////
//////      }
//////    }
////
////
//////    logfileptr->OFS()<<"Solution after each step is "<<STEPS_SP<<std::endl;
////    //Print B
////    logfileptr->OFS()<<"Solution is "<<B<<std::endl;
////




  }


} // namespace LIBCHOLESKY


#endif 
