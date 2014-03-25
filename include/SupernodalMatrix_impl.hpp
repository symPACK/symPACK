#ifndef _SUPERNODAL_MATRIX_IMPL_HPP_
#define _SUPERNODAL_MATRIX_IMPL_HPP_

#include "SupernodalMatrix.hpp"

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


  Global_.FindSupernodes(ETree_,cc,Xsuper_);


//  UpdatesCount.Resize(Xsuper_.Size());
//  for(Int I = 1; I<Xsuper_.Size();++I){
//    Int first_col = Xsuper_[I-1];
//    Int last_col = Xsuper_[I]-1;
//    for(Int 
//  }


  IntNumVec xlindx,lindx;
  Global_.SymbolicFactorization(ETree_,cc,Xsuper_,xlindx,lindx);

  //build supernode membership
  SupMembership_.Resize(iSize_);
  Int cur_snode_idx = 1;
  for(Int i = 1; i<=iSize_;++i){
    SupMembership_(i-1) = cur_snode_idx;
    if(Xsuper_(cur_snode_idx) == i+1){
      cur_snode_idx++;
    } 
  }
  logfileptr->OFS()<<"Membership list is "<<SupMembership_<<std::endl;



  SupETree_ = ETree_.ToSupernodalETree(Xsuper_);
  logfileptr->OFS()<<"Supernodal Etree is "<<SupETree_<<std::endl;

  Mapping_ = pMapping;

  //copy
    for(Int I=1;I<Xsuper_.m();I++){
      Int fc = Xsuper_(I-1);
      Int lc = Xsuper_(I)-1;
      Int iWidth = lc-fc+1;
      Int fi = xlindx(I-1);
      Int li = xlindx(I)-1;
      Int iHeight = li-fi+1;

      Int iDest = pMapping.Map(I-1,I-1);

      //parse the first column to create the NZBlock
      if(iam==iDest){
        LocalSupernodes_.push_back( new SuperNode2<T>(I,fc,lc,iHeight));
        SuperNode2<T> & snode = *LocalSupernodes_.back();


        for(Int idx = fi; idx<=li;idx++){
          Int iStartRow = lindx(idx-1);
          Int iPrevRow = iStartRow;
          Int iContiguousRows = 1;
          for(Int idx2 = idx+1; idx2<=li;idx2++){
            Int iCurRow = lindx(idx2-1);
            if(iStartRow == lindx(fi-1)){
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

        //snode.Shrink();
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
          SuperNode2<T> & snode = *LocalSupernodes_.back();

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
            iLRow = lindx(idx-1);
            //logfileptr->OFS()<<"Looking at L("<<iLRow<<","<<i<<")"<<std::endl;
            if( iLRow == iRowind){
              Int iNZBlockIdx = snode.FindBlockIdx(iLRow);

              T elem = pdNzVal[iStartIdxCopy + idxA];

              NZBlock2<T> & pDestBlock = snode.GetNZBlock(iNZBlockIdx);
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


template <typename T> inline bool SupernodalMatrix<T>::FindNextUpdate(SuperNode2<T> & src_snode, Int & src_nzblk_idx, Int & src_first_row, Int & src_last_row, Int & tgt_snode_id){
  //src_nzblk_idx is the last nzblock index examined
  //src_first_row is the first row updating the supernode examined
  //src_last_row is the last row updating the supernode examined

  //if tgt_snode_id == 0 , this is the first call to the function
  if(tgt_snode_id == 0){

    if(src_snode.NZBlockCnt() == 1 ){
      return false;
    }
    src_nzblk_idx = 1;
    src_first_row = src_snode.GetNZBlock(src_nzblk_idx).GIndex();
  }
  else{
    //find the block corresponding to src_last_row
    src_nzblk_idx = src_snode.FindBlockIdx(src_last_row);
    if(src_nzblk_idx==src_snode.NZBlockCnt()){
      return false;
    }
    else{
      NZBlock2<T> & nzblk = src_snode.GetNZBlock(src_nzblk_idx);
      Int src_lr = nzblk.GIndex()+nzblk.NRows()-1;
      if(src_last_row == src_lr){
        src_nzblk_idx++;
        if(src_nzblk_idx==src_snode.NZBlockCnt()){
          return false;
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
  NZBlock2<T> & nzblk = src_snode.GetNZBlock(src_nzblk_idx);
  Int src_fr = max(src_first_row,nzblk.GIndex());
  Int src_lr = nzblk.GIndex()+nzblk.NRows()-1;

  Int tgt_snode_id_first = SupMembership_(src_fr-1);
  Int tgt_snode_id_last = SupMembership_(src_lr-1);
  if(tgt_snode_id_first == tgt_snode_id_last){
    //this can be a zero row in the src_snode
    src_last_row = Xsuper_(tgt_snode_id_first)-1;
    tgt_snode_id = tgt_snode_id_first;

    //look into other nzblk
    Int new_blk_idx = src_snode.FindBlockIdx(src_last_row);
    new_blk_idx = (new_blk_idx>src_snode.NZBlockCnt())?new_blk_idx-1:new_blk_idx+1;
    NZBlock2<T> & lb = src_snode.GetNZBlock(new_blk_idx-1); 
    src_last_row = lb.GIndex() + lb.NRows() - 1;
  }
  else{
    src_last_row = Xsuper_(tgt_snode_id_first)-1;
    tgt_snode_id = tgt_snode_id_first;
  }

  return true;

}





template <typename T> inline bool SupernodalMatrix<T>::FindPivot(SuperNode2<T> & src_snode, SuperNode2<T> & tgt_snode,Int & pivot_idx, Int & pivot_fr, Int & pivot_lr){
  //Determine which columns will be updated
  pivot_idx = 0;
  pivot_fr = 0;
  pivot_lr = 0;
  for(int blkidx=0;blkidx<src_snode.NZBlockCnt();++blkidx){
    NZBlock2<T> & nzblk = src_snode.GetNZBlock(blkidx);
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

template <typename T> void SupernodalMatrix<T>::UpdateSuperNode(SuperNode2<T> & src_snode, SuperNode2<T> & tgt_snode, Int &pivot_idx, Int  pivot_fr = I_ZERO){
    NZBlock2<T> & pivot_nzblk = src_snode.GetNZBlock(pivot_idx);
    if(pivot_fr ==I_ZERO){
      pivot_fr = pivot_nzblk.GIndex();
    }
    assert(pivot_fr >= pivot_nzblk.GIndex());
    Int pivot_lr = min(pivot_nzblk.GIndex() +pivot_nzblk.NRows() -1, pivot_fr + tgt_snode.Size() -1);
    T * pivot = & pivot_nzblk.Nzval(pivot_fr-pivot_nzblk.GIndex(),0);

    Int tgt_updated_fc =  pivot_fr - tgt_snode.FirstCol();
    Int tgt_updated_lc =  pivot_fr - tgt_snode.FirstCol();
    logfileptr->OFS()<<"pivotidx = "<<pivot_idx<<" Cols "<<tgt_snode.FirstCol()+tgt_updated_fc<<" to "<< tgt_snode.FirstCol()+tgt_updated_lc <<" will be updated"<<std::endl;

    for(int src_idx=pivot_idx;src_idx<src_snode.NZBlockCnt();++src_idx){

      NZBlock2<T> & src_nzblk = src_snode.GetNZBlock(src_idx);
      Int src_fr = max(pivot_fr, src_nzblk.GIndex()) ;
      Int src_lr = src_nzblk.GIndex()+src_nzblk.NRows()-1;
      
      do{
        Int tgt_idx = tgt_snode.FindBlockIdx(src_fr);
        NZBlock2<T> & tgt_nzblk = tgt_snode.GetNZBlock(tgt_idx);
        Int tgt_fr = tgt_nzblk.GIndex();
        Int tgt_lr = tgt_nzblk.GIndex()+tgt_nzblk.NRows()-1;

        Int update_fr = max(tgt_fr, src_fr);
        Int update_lr = min(tgt_lr, src_lr);
        

        //Update tgt_nzblk with src_nzblk
        logfileptr->OFS()<<"Updating SNODE "<<tgt_snode.Id()<<" Block "<<tgt_idx<<" ["<<update_fr<<".."<<update_lr<<"] with SNODE "<<src_snode.Id()<<" Block "<<src_idx<<" ["<<update_fr<<".."<<update_lr<<"]"<<endl;

        if(tgt_idx+1<tgt_snode.NZBlockCnt()){
          NZBlock2<T> & next_tgt_nzblk = tgt_snode.GetNZBlock(tgt_idx+1);
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
  Int iam;
  MPI_Comm_rank(pComm, &iam);
  Int np;
  MPI_Comm_size(pComm, &np);

 
  //dummy right looking cholesky factorization
  for(Int I=1;I<Xsuper_.m();I++){
    Int src_first_col = Xsuper_(I-1);
    Int src_last_col = Xsuper_(I)-1;
    Int iOwner = Mapping_.Map(I-1,I-1);
    //If I own the column, factor it
    if( iOwner == iam ){
      Int iLocalI = (I-1) / np +1 ;
      SuperNode2<T> & src_snode = *LocalSupernodes_[iLocalI -1];
      logfileptr->OFS()<<"Supernode "<<I<<"("<<src_snode.Id()<<") is on P"<<iOwner<<" local index is "<<iLocalI<<std::endl; 
     
      assert(src_snode.Id() == I); 


      //Do all my updates (Local and remote)
      //Get Row structure
      //FIXME

      std::set<Int> LRowStruct;
      std::vector<Int> supLStruct;

      LRowStruct.clear();
      Global_.GetSuperLRowStruct(ETree_, Xsuper_, I , LRowStruct);
      supLStruct.clear();
      for (std::set<Int>::iterator it=LRowStruct.begin(); it!=LRowStruct.end(); ++it){
        Int sup = SupMembership_(*it-1);
          if(*it == Xsuper_(sup-1)){
            supLStruct.push_back(sup);
          }
      }
 
      logfileptr->OFS()<<"  Supernodal Row structure is "<<supLStruct<<std::endl;


      std::vector<char> src_nzblocks;

      //Parse the row structure to look for remote factors to Recv
      for(int Jidx = 0; Jidx < supLStruct.size(); ++Jidx){
        Int src_snode_id = supLStruct[Jidx];
        Int iSource = Mapping_.Map(src_snode_id-1,src_snode_id-1);
//        if(iSource == iam){
//          logfileptr->OFS()<<"Supernode "<<src_snode_id<<" is updating Supernode "<<I<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;
//          Int iLocalJ = (tgt_snode_id-1) / np +1 ;
//          SuperNode2<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];
//          UpdateSuperNode(src_snode,tgt_snode,src_nzblk_idx, src_first_row);
//        }
//        else{
        if(iSource != iam){

          if(src_nzblocks.size()==0){
            //Create an upper bound buffer  (we have the space to store the first header)
            size_t num_bytes = 0;
            for(Int blkidx=0;blkidx<src_snode.NZBlockCnt();++blkidx){
              num_bytes += src_snode.GetNZBlock(blkidx).NRows()*(src_snode.Size()*sizeof(T) + sizeof(NZBlockHeader<T>)+ sizeof(NZBlock2<T>));
            }
            logfileptr->OFS()<<"We allocate a buffer of size "<<num_bytes<<std::endl;
            src_nzblocks.resize(num_bytes);
          }


          //MPI_Recv
          MPI_Status recv_status;
          char * recv_buf = &src_nzblocks[0]+sizeof(NZBlock2<T>);
          MPI_Recv(recv_buf,src_nzblocks.size(),MPI_BYTE,iSource,I,pComm,&recv_status);

          //Resize the buffer to the actual number of bytes received
          int bytes_received = 0;
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
          src_nzblocks.resize(bytes_received+sizeof(NZBlock2<T>));



          logfileptr->OFS()<<"Supernode "<<I<<" is updated by Supernode "<<src_snode_id<<", received "<<bytes_received<<" bytes"<<std::endl;
          //Create the dummy supernode for that data
          SuperNode2<T> dist_src_snode(src_snode_id,Xsuper_[src_snode_id-1],Xsuper_[src_snode_id]-1,&src_nzblocks);

          Int idx = 0;
          UpdateSuperNode(dist_src_snode,src_snode,idx);


          //restore to its capacity
          src_nzblocks.resize(src_nzblocks.capacity());


        }


      }
      //clear the buffer
      { vector<char>().swap(src_nzblocks);  }

      logfileptr->OFS()<<"  Factoring node "<<I<<std::endl;

      //Factorize Diagonal block
      NZBlock2<T> & diagonalBlock = src_snode.GetNZBlock(0);
      lapack::Potrf( 'L', src_snode.Size(), diagonalBlock.Nzval(), diagonalBlock.LDA());
      logfileptr->OFS()<<"    Diagonal block factored node "<<I<<std::endl;

      for(int blkidx=1;blkidx<src_snode.NZBlockCnt();++blkidx){
        NZBlock2<T> & nzblk = src_snode.GetNZBlock(blkidx);

          //Update lower triangular blocks
          blas::Trsm('R','L','T','N',nzblk.NRows(),src_snode.Size(), 1.0,  diagonalBlock.Nzval(), diagonalBlock.LDA(), nzblk.Nzval(), nzblk.LDA());
          logfileptr->OFS()<<"    "<<blkidx<<"th subdiagonal block updated node "<<I<<std::endl;
      }
      //Send my factor to my ancestors. 
      Int tgt_snode_id = 0;
      Int src_nzblk_idx = 0;
      Int src_first_row = 0;
      Int src_last_row = 0;
      while(FindNextUpdate(src_snode, src_nzblk_idx, src_first_row, src_last_row, tgt_snode_id)){ 

        Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
        if(iTarget != iam){
          logfileptr->OFS()<<"Supernode "<<tgt_snode_id<<" is updated by Supernode "<<I<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;

          Int J = tgt_snode_id;

          Int tgt_first_col = Xsuper_(J-1);
          Int tgt_last_col = Xsuper_(J)-1;

          //Send
          NZBlock2<T> * pivot_nzblk = &src_snode.GetNZBlock(src_nzblk_idx);

          int local_first_row = src_first_row-pivot_nzblk->GIndex();

          char * start_buf = reinterpret_cast<char*>(&pivot_nzblk->Nzval(local_first_row,0));
          size_t num_bytes = src_snode.End() - start_buf;

          //Create a header
          NZBlockHeader<T> * header = new NZBlockHeader<T>(pivot_nzblk->NRows()-local_first_row,pivot_nzblk->NCols(),src_first_row,1);


          MPI_Datatype type;
          int lens[2];
          MPI_Aint disps[2];
          MPI_Datatype types[2];
          Int err = 0;


          /* define a struct that describes all our data */
          lens[0] = sizeof(NZBlockHeader<T>);
          lens[1] = num_bytes;
          MPI_Address(header, &disps[0]);
          MPI_Address(start_buf, &disps[1]);
          types[0] = MPI_BYTE;
          types[1] = MPI_BYTE;
          MPI_Type_struct(2, lens, disps, types, &type);
          MPI_Type_commit(&type);

          /* send to target */
          //FIXME : maybe this can be replaced by a scatterv ?
          MPI_Send(MPI_BOTTOM,1,type,iTarget,tgt_snode_id,pComm);

          MPI_Type_free(&type);

          delete header;

          logfileptr->OFS()<<"     Send factor "<<I<<" to node"<<J<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<" ("<<num_bytes+sizeof(*header)<<" bytes)"<<std::endl;

        }
      }

      //Update my ancestors right away. 
      tgt_snode_id = 0;
      while(FindNextUpdate(src_snode, src_nzblk_idx, src_first_row, src_last_row, tgt_snode_id)){ 

        Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
        if(iTarget == iam){
          logfileptr->OFS()<<"Supernode "<<tgt_snode_id<<" is updated by Supernode "<<I<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;

          Int iLocalJ = (tgt_snode_id-1) / np +1 ;
          SuperNode2<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];

          UpdateSuperNode(src_snode,tgt_snode,src_nzblk_idx, src_first_row);
        }
      }
    }

///
///
///    else{
///      //If I have something updated by this column, update it
///      // i.e do I have ancestors of I ? 
///
///
///
///
///
/////Not possible
/////      //Update my ancestors. 
/////      tgt_snode_id = 0;
/////      while(FindNextUpdate(src_snode, src_nzblk_idx, src_first_row, src_last_row, tgt_snode_id)){ 
/////
/////        Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
/////        if(iTarget == iam){
/////          logfileptr->OFS()<<"Supernode "<<tgt_snode_id<<" is updated by Supernode "<<I<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;
/////
/////          Int iLocalJ = (tgt_snode_id-1) / np +1 ;
/////          SuperNode2<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];
/////
/////          UpdateSuperNode(src_snode,tgt_snode,src_nzblk_idx, src_first_row);
/////        }
/////      }
///
///
///
///
///
///
///  
///      Int J = SupETree_.Parent(I-1);
///      while(J != 0){
///
///        //Is it updated by I ? Is I nnz in LRowStruct ?
///
///        Int iAOwner = Mapping_.Map(J-1,J-1);
///        //If I own the ancestor supernode, update it
///        if( iAOwner == iam ){
///          Int iLocalJ = (J-1) / np +1 ;
///          SuperNode2<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];
///          //upper bound, can be replaced by something tighter
///          size_t num_bytes = tgt_snode.BytesToEnd(0);
///          std::vector<char> src_nzblocks(num_bytes);
///          //Recv
///          MPI_Status recv_status;
///          MPI_Recv(&src_nzblocks[0],num_bytes,MPI_BYTE,iOwner,TAG_FACTOR,pComm,&recv_status);
///          int bytes_received = 0;
///          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);
///          src_nzblocks.resize(bytes_received);
///  
///
///          SuperNode2<T> src_snode(I,src_first_col,src_last_col,&src_nzblocks);
///
///          logfileptr->OFS()<<"     Recv factor "<<I<<" from P"<<iOwner<<" to update snode"<<J<<" on P"<<iam<<std::endl;       
///          //Update
///
///          UpdateSuperNode(src_snode,tgt_snode,iLocalJ);
///        }
///        J = SupETree_.Parent(J-1);
///      }
///    }
///

    upcxx::barrier();
  }
}


} // namespace LIBCHOLESKY


#endif 
