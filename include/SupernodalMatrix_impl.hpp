#ifndef _SUPERNODAL_MATRIX_IMPL_HPP_
#define _SUPERNODAL_MATRIX_IMPL_HPP_

#include "SupernodalMatrix.hpp"

#include <queue>


#define BLOCKSIZE src_snode.Size()

#define TAG_INDEX 0
#define TAG_NZVAL 1
#define TAG_COUNT 2

#define PACKING


namespace LIBCHOLESKY{

  template<typename T> void SupernodalMatrix<T>::AddOutgoingComm(AsyncComms & outgoingSend, Int src_snode_id, Int src_snode_size, Int src_first_row, NZBlockDesc & pivot_desc, Int nzblk_cnt, T * nzval_ptr, Int nz_cnt){

    outgoingSend.push_back(new Icomm(3*sizeof(Int) + nzblk_cnt*sizeof(NZBlockDesc) + nz_cnt*sizeof(T),MPI_REQUEST_NULL));

    *outgoingSend.back()<<src_snode_id;
    *outgoingSend.back()<<nzblk_cnt;
    NZBlockDesc * dest_blocks_ptr = reinterpret_cast<NZBlockDesc *>(outgoingSend.back()->back());
    Serialize(*outgoingSend.back(),&pivot_desc,nzblk_cnt);
    dest_blocks_ptr->GIndex = src_first_row;
    dest_blocks_ptr->Offset = pivot_desc.Offset + (src_first_row - pivot_desc.GIndex)*src_snode_size;
    *outgoingSend.back()<<nz_cnt;
    Serialize(*outgoingSend.back(),nzval_ptr,nz_cnt);
  }



  template <typename T> SupernodalMatrix<T>::SupernodalMatrix(){
  }

  template <typename T> SupernodalMatrix<T>::SupernodalMatrix(const DistSparseMatrix<T> & pMat, Int maxSnode,MAPCLASS & pMapping, Int maxIsend, Int maxIrecv, MPI_Comm & pComm ){
    this->pComm = pComm;
    Int iam;
    MPI_Comm_rank(pComm, &iam);
    Int np;
    MPI_Comm_size(pComm, &np);

    maxIsend_ = maxIsend;
    maxIrecv_ = maxIrecv;

    iSize_ = pMat.size;
    Local_ = pMat.GetLocalStructure();
    Local_.ToGlobal(Global_);

//    ETree * tmp = new ETree();
//    tmp->ConstructETree2(Global_);
    
    ETree_.ConstructETree(Global_);

    IntNumVec cc,rc;
    Global_.GetLColRowCount2(ETree_,cc,rc);
    IntNumVec perm2 = ETree_.SortChildren(cc);

    IntNumVec poperm(Size());
      for(Int i =0; i<poperm.m();++i){
        poperm[i] = ETree_.FromPostOrder(i+1);
      }
  
    ETree_.PermuteTree(perm2);

    IntNumVec poperm2(Size());
      for(Int i =0; i<poperm2.m();++i){
        poperm2[i] = ETree_.FromPostOrder(i+1);
      }
#ifdef _DEBUG_
    logfileptr->OFS()<<"colcnt "<<cc<<std::endl;
    logfileptr->OFS()<<"rowcnt "<<rc<<std::endl;
#endif


#ifdef _DEBUG_
    logfileptr->OFS()<<"colcnt "<<cc<<std::endl;
#endif

    Global_.FindSupernodes(ETree_,cc,SupMembership_,Xsuper_,maxSnode);

#ifdef _DEBUG_
    logfileptr->OFS()<<"Membership list is "<<SupMembership_<<std::endl;
    logfileptr->OFS()<<"xsuper "<<Xsuper_<<std::endl;
#endif

#ifdef RELAXED_SNODE
    Global_.SymbolicFactorizationRelaxed(ETree_,cc,Xsuper_,SupMembership_,xlindx_,lindx_);


    logfileptr->OFS()<<"Membership list is "<<SupMembership_<<std::endl;
    logfileptr->OFS()<<"xsuper "<<Xsuper_<<std::endl;
    logfileptr->OFS()<<"colcnt "<<cc<<std::endl;
    logfileptr->OFS()<<"xlindx "<<xlindx_<<std::endl;
    logfileptr->OFS()<<"lindx "<<lindx_<<std::endl;

    ETree tmp2 = ETree_.ToSupernodalETree(Xsuper_);
    logfileptr->OFS()<<"Supernodal ETree "<<tmp2<<std::endl;


#else
    Global_.SymbolicFactorization2(ETree_,cc,Xsuper_,SupMembership_,xlindx_,lindx_);
    IntNumVec perm3;
    IntNumVec newPerm(Size());
if(0){
    Global_.RefineSupernodes(ETree_, SupMembership_, Xsuper_, xlindx_, lindx_, perm3);


      logfileptr->OFS()<<"po perm: "<<poperm<<std::endl;
      logfileptr->OFS()<<"po perm2: "<<poperm2<<std::endl;
      logfileptr->OFS()<<"perm2: "<<perm2<<std::endl;
      logfileptr->OFS()<<"perm3: "<<perm3<<std::endl;
newPerm = perm3;
//      for(Int i =0; i<perm3.m();++i){
//        newPerm[perm2[i]-1] = poperm[i];
//      }
//      logfileptr->OFS()<<"new perm == po perm 2 ?: "<<newPerm<<std::endl;
//      poperm2 = newPerm;
//
//      for(Int i =0; i<perm3.m();++i){
//        newPerm[perm3[i]-1] = perm2[i];
//      }
//      logfileptr->OFS()<<"perm2 and perm3: "<<newPerm<<std::endl;
//
//      perm3 = newPerm;
}

      for(Int i =0; i<newPerm.m();++i){
        newPerm[i] = ETree_.FromPostOrder(i+1);
      }
//      logfileptr->OFS()<<"new perm: "<<newPerm<<std::endl;
   perm_ = newPerm;

#endif

    GetUpdatingSupernodeCount(UpdateCount_,UpdateWidth_);

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

#ifdef _DEBUG_
      logfileptr->OFS()<<"Supernode "<<I<<" is owned by P"<<iDest<<std::endl;
#endif

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
        

        //snode.DumpITree();
      }

      //Distribute the data

      //look at the owner of the first column
      Int numColFirst = iSize_ / np;
      Int prevOwner = -1;
      DblNumVec aRemoteCol;
      T * pdNzVal = NULL;
      Int iStartIdxCopy = 0;

      //boolean to check if the prevOwner has sent all his columns or not
      bool allsend = false;
      //copy the data from A into this Block structure
      for(Int i = fc;i<=lc;i++){
        //corresponding column in the unsorted matrix A
//        Int orig_i = ETree_.FromPostOrder(perm(i-1));
        Int orig_i = perm_[i-1];
        Int iOwner = std::min((orig_i-1)/numColFirst,np-1);

        if(iOwner != prevOwner || !allsend){

          prevOwner = iOwner;
          //we need to transfer the bulk of columns
          //first compute the number of cols we need to transfer
          Int iNzTransfered = 0;
          Int iLCTransfered = lc;

          Int prevcol = orig_i-1;
          for(Int j =i;j<=lc;++j){
            //corresponding column in the unsorted matrix A
//            Int orig_j = ETree_.FromPostOrder(perm(j-1));
            Int orig_j = perm_[j-1];
       
            iOwner = std::min((orig_j-1)/numColFirst,np-1);
            //check if the column is owned by the same processor 
            //as the previous column and that they are contiguous
            //in the postordered matrix
            if(iOwner == prevOwner && prevcol+1==orig_j){
              Int nrows = Global_.colptr(orig_j) - Global_.colptr(orig_j-1);
              iNzTransfered+=nrows;

              iLCTransfered =orig_j;
              prevcol = orig_j;
            }
            else{
              //check if we have sent everything (po columns are contiguous) or not
              if(prevcol+1!=orig_j){
                allsend = false;;
              }
              else{
                allsend = true;
              }
              break;
            }
          } 

#ifdef _DEBUG_
logfileptr->OFS()<<fc<<" Col "<<orig_i<<" to "<<iLCTransfered<<" are owned by P"<<prevOwner<<std::endl;
#endif

          //if data needs to be transfered
          if(iDest!=prevOwner){
            if(iam == iDest){
              aRemoteCol.Resize(iNzTransfered);

              //MPI_Recv
              MPI_Recv(aRemoteCol.Data(),iNzTransfered*sizeof(T),MPI_BYTE,prevOwner,0,pComm,MPI_STATUS_IGNORE);
              iStartIdxCopy = 0;
            } 
            else if (iam == prevOwner){
              //USE THE PERM OBTAINED AFTER SORTING THE CHILDREN
              Int local_i = (orig_i-(numColFirst)*iam);
              Int iColptrLoc = Local_.colptr(local_i-1);
              Int iRowIndLoc = Local_.rowind(iColptrLoc-1);
              Int iLastRowIndLoc = Local_.rowind(Local_.colptr(local_i)-1-1);
              const T * pdData = &pMat.nzvalLocal(iColptrLoc-1);

#ifdef _DEBUG_
          logfileptr->OFS()<<"Sent elem are "<<std::endl;
          for(int i = 0; i<iNzTransfered;++i){
            logfileptr->OFS()<<pdData[i]<<std::endl;
          }
#endif
              //MPI_send        
#ifdef _DEBUG_
              assert(iDest<np);
#endif
              MPI_Send((void*)pdData,iNzTransfered*sizeof(T),MPI_BYTE,iDest,0,pComm);
            }
          }
        }

        //copy the data if I own it
        if(iam == iDest){
          SuperNode<T> & snode = *LocalSupernodes_.back();

          //isData transfered or local
          if(iam!=prevOwner){
            pdNzVal = aRemoteCol.Data();

#ifdef _DEBUG_
logfileptr->OFS()<<aRemoteCol<<std::endl;
#endif
          }
          else{
            //USE THE PERM OBTAINED AFTER SORTING THE CHILDREN
            Int local_i = (orig_i-(numColFirst)*iam);
            Int iColptrLoc = Local_.colptr(local_i-1);
            pdNzVal = (T*)(&pMat.nzvalLocal(iColptrLoc-1));
            //logfileptr->OFS()<<"pdNzVal is the local pMat"<<std::endl;
            iStartIdxCopy = 0;
          }

          //Copy the data from pdNzVal in the appropriate NZBlock
          Int iGcolptr = Global_.colptr(orig_i-1);
          Int iNextColptr = Global_.colptr(orig_i);
          Int iRowind = Global_.rowind(iGcolptr-1);
          Int iNrows = iNextColptr - iGcolptr;
          Int idxA = 0;
          Int iLRow = 0;
          Int firstrow = fi + i-fc;
          for(Int idx = firstrow; idx<=li;idx++){
            iLRow = lindx_(idx-1);
            // Original row index in the unsorted matrix A
//            Int orig_iLRow = ETree_.FromPostOrder(perm(lindx_(idx-1)-1));
            Int orig_iLRow = perm_[lindx_[idx-1]-1];
            if( orig_iLRow == iRowind){
              Int iNZBlockIdx = snode.FindBlockIdx(iLRow);
#ifdef _DEBUG_
              if(!(iNZBlockIdx>=0 && iNZBlockIdx<snode.NZBlockCnt())){
                snode.DumpITree();
                logfileptr->OFS()<<iNZBlockIdx<<endl;
              }
              assert(iNZBlockIdx>=0 && iNZBlockIdx<snode.NZBlockCnt());
#endif
              T elem = pdNzVal[iStartIdxCopy + idxA];

                NZBlockDesc & desc = snode.GetNZBlockDesc(iNZBlockIdx);
                T * dest = snode.GetNZval(desc.Offset);

                //find where we should put it
                Int localCol = i - fc;
                Int localRow = iLRow - desc.GIndex;

                dest[localRow*snode.Size()+localCol] = elem;

              if(iGcolptr+idxA+1<iNextColptr){
                idxA++;
                iRowind = Global_.rowind(iGcolptr+idxA-1);
              }
              else{
                idxA++;
                break;
              }
            }
          }
          //if(iam!=prevOwner){
            iStartIdxCopy+=iNrows;
          //}
        }

      }
#ifdef _DEBUG_
      logfileptr->OFS()<<"--------------------------------------------------"<<std::endl;
#endif
    }

#ifdef _DEBUG_
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

        logfileptr->OFS()<<"_______________________________"<<std::endl;
        }
      }
    }
#endif




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

#ifndef _DEBUG_
  #define nodebugtmp
  #define _DEBUG_
#endif



#ifdef nodebugtmp
  #undef _DEBUG_
#endif



#ifdef _DEBUG_
      logfileptr->OFS()<<"Supernode "<<s<<" updates: ";
#endif

      for(Int row_idx = fi; row_idx<=li;++row_idx){
        Int row = lindx_(row_idx-1);
        Int supno = SupMembership_(row-1);

        if(marker(supno-1)!=s && supno!=s){

#ifdef _DEBUG_
          logfileptr->OFS()<<supno<<" ";
#endif
          ++sc(supno-1);
          marker(supno-1) = s;

          mw(supno-1) = max(mw(supno-1),last_col - first_col+1);

        }
      }

#ifdef _DEBUG_
      logfileptr->OFS()<<std::endl;
#endif

#ifdef nodebugtmp
  #undef _DEBUG_
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

struct SnodeUpdate{
  Int tgt_snode_id;
  Int src_fr;
  SnodeUpdate(Int aSnodeId, Int aSrcFr):tgt_snode_id(aSnodeId),src_fr(aSrcFr){};
};

template<typename T> inline void SupernodalMatrix<T>::FindUpdates(SuperNode<T> & src_snode, std::list<SnodeUpdate> & updates  ){
  Int src_owner = Mapping_.Map(src_snode.Id()-1,src_snode.Id()-1);
  Int src_nzblk_idx = src_owner == iam ? 1 : 0;
  updates.clear();

  Int last_snode_id = -1;
  for(src_nzblk_idx; src_nzblk_idx< src_snode.NZBlockCnt(); ++src_nzblk_idx){
    NZBlockDesc & src_blk = src_snode.GetNZBlockDesc(src_nzblk_idx);
     Int src_fr = src_blk.GIndex;
     Int src_nrows = src_snode.NRows(src_nzblk_idx);
     Int src_lr = src_fr + src_nrows - 1; 
     for(Int row = src_fr; row<= src_lr; ++row){
      Int tgt_snode_id = SupMembership_[row-1];
      if(last_snode_id!= tgt_snode_id){
        updates.push_back(SnodeUpdate(tgt_snode_id,row));
        last_snode_id = tgt_snode_id;
      }
    }
  }
}






  template <typename T> inline bool SupernodalMatrix<T>::FindNextUpdate2(SuperNode<T> & src_snode, Int & tgt_snode_id, Int & f_ur, Int & f_ub, Int & n_ur, Int & n_ub){

//    TIMER_START(FIND_UPDATE);

    if(tgt_snode_id==0){   
      Int iOwner = Mapping_.Map(src_snode.Id()-1,src_snode.Id()-1);
      f_ub = iOwner==iam?1:0;
      n_ub = f_ub;
      n_ur = 0;
    }
    else{
      f_ub = n_ub;
    }

    if(src_snode.NZBlockCnt()>f_ub){
      NZBlockDesc * cur_desc = &src_snode.GetNZBlockDesc(f_ub); 
      f_ur = max(n_ur,cur_desc->GIndex); 
      //find the snode updated by that row
      tgt_snode_id = SupMembership_[f_ur-1];
      Int tgt_lc = Xsuper_[tgt_snode_id]-1;

      //or use FindBlockIdx
//      if(f_ub<src_snode.NZBlockCnt()-1){
        Int src_tgt_lb = src_snode.FindBlockIdx(tgt_lc);
        //if tgt_lc not found in the current column we need to resume at the next block 
        if(src_tgt_lb==-1){
          for(n_ub;n_ub<src_snode.NZBlockCnt();++n_ub){
            cur_desc = &src_snode.GetNZBlockDesc(n_ub);
            if(cur_desc->GIndex > tgt_lc){
              break;
            }
          }
          if(n_ub<src_snode.NZBlockCnt()){
            n_ur = cur_desc->GIndex;
          }
          else{
            n_ur = -1;
          }
        }
        else{
            n_ub = src_tgt_lb;
            cur_desc = &src_snode.GetNZBlockDesc(n_ub);
            if(cur_desc->GIndex + src_snode.NRows(n_ub)-1>tgt_lc){
              n_ur = tgt_lc+1;
            }
            else{
              ++n_ub;
              n_ur = (n_ub<src_snode.NZBlockCnt())?src_snode.GetNZBlockDesc(n_ub).GIndex:-1;
            }
        }


      

//      }
//      else{
//           
//          //tgt_lc is in current block or not
//          if(
//      }

      //src_snode updates tgt_snode_id. Then we need to look from row n_ur and block l_ub
      return true;
    }
    else{
      return false;
    }
  }
























  template <typename T> inline bool SupernodalMatrix<T>::FindNextUpdate(SuperNode<T> & src_snode, Int & src_nzblk_idx, Int & src_first_row, Int & src_last_row, Int & tgt_snode_id){
    //src_nzblk_idx is the last nzblock index examined
    //src_first_row is the first row updating the supernode examined
    //src_last_row is the last row updating the supernode examined

//    TIMER_START(FIND_UPDATE);
    
    bool returnval = true;

    //if tgt_snode_id == 0 , this is the first call to the function
    if(tgt_snode_id == 0){

      //find the first sub diagonal block
      Int iOwner = Mapping_.Map(src_snode.Id()-1,src_snode.Id()-1);
      src_nzblk_idx= iOwner==iam?1:0;

      if(src_nzblk_idx < 0  || src_nzblk_idx>=src_snode.NZBlockCnt()){
        returnval = false;
      }
      else{
#ifdef _DEBUG_
            assert(src_nzblk_idx<src_snode.NZBlockCnt());
#endif
        src_first_row = src_snode.GetNZBlockDesc(src_nzblk_idx).GIndex;
      }
    }
    else{
      //find the block corresponding to src_last_row
      src_nzblk_idx = src_snode.FindBlockIdx(src_last_row);
      
#ifdef _DEBUG_
      assert(src_nzblk_idx!=-1);
       assert(src_nzblk_idx<src_snode.NZBlockCnt());
#endif
        NZBlockDesc & desc = src_snode.GetNZBlockDesc(src_nzblk_idx);
        Int src_lr = desc.GIndex + src_snode.NRows(src_nzblk_idx)-1;

        if(src_last_row == src_lr){
          src_nzblk_idx++;
          if(src_nzblk_idx==src_snode.NZBlockCnt()){
            returnval = false;
          }
          else{
#ifdef _DEBUG_
            assert(src_nzblk_idx<src_snode.NZBlockCnt());
#endif
            src_first_row = src_snode.GetNZBlockDesc(src_nzblk_idx).GIndex;
          }
        }
        else{
          src_first_row = src_last_row+1;
        }
    }

    if(returnval){
    //Now we try to find src_last_row
    NZBlockDesc & desc = src_snode.GetNZBlockDesc(src_nzblk_idx);
#ifdef _DEBUG_
    assert(src_first_row >= desc.GIndex);
#endif
    Int src_fr = max(src_first_row,desc.GIndex);
    Int src_lr = desc.GIndex+src_snode.NRows(src_nzblk_idx)-1;

    Int tgt_snode_id_first = SupMembership_(src_fr-1);
    Int tgt_snode_id_last = SupMembership_(src_lr-1);

    if(tgt_snode_id_first == tgt_snode_id_last){
      //this can be a zero row in the src_snode
      tgt_snode_id = tgt_snode_id_first;

      Int last_tgt_col = Xsuper_(tgt_snode_id_first)-1;

      //We might have another block updating the same snode
      if(src_lr < last_tgt_col){

#ifdef FAST_INDEX_SEARCH
        src_last_row = src_lr;
        Int blkidx = src_snode.FindBlockIdx(last_tgt_col);
        if(blkidx<0){
          if(blkidx==-(iSize_+1)){
            blkidx = src_snode.NZBlockCnt()-1;
          }
          else{
            blkidx = src_snode.FindBlockIdx(-blkidx)-1;
          }
        }
        src_last_row = min(last_tgt_col,src_snode.GetNZBlockDesc(blkidx).GIndex + src_snode.NRows(blkidx) -1);
#else
        src_last_row = src_lr;
        for(Int blkidx = src_nzblk_idx; blkidx<src_snode.NZBlockCnt();++blkidx){
          if(src_snode.GetNZBlockDesc(blkidx).GIndex <= last_tgt_col){
            src_last_row = min(last_tgt_col,src_snode.GetNZBlockDesc(blkidx).GIndex + src_snode.NRows(blkidx) -1);
          }
          else{
            break;
          }
        }
#endif
      }
      else{
        src_last_row = last_tgt_col;
      }
#ifdef _DEBUG_
      assert(src_last_row<= Xsuper_(tgt_snode_id_first)-1);
#endif
    }
    else{
      src_last_row = Xsuper_(tgt_snode_id_first)-1;
      tgt_snode_id = tgt_snode_id_first;
    }


#ifdef _DEBUG_
    assert(src_last_row>=src_first_row);
    assert(src_first_row>= Xsuper_(tgt_snode_id-1));
    assert(src_last_row<= Xsuper_(tgt_snode_id)-1);
    assert(src_first_row >= src_fr);
#endif
    }

    return returnval;
  }


#ifdef SINGLE_BLAS
  template <typename T> inline void SupernodalMatrix<T>::UpdateSuperNode(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode, Int &pivot_idx, NumMat<T> & tmpBuf,IntNumVec & src_colindx, IntNumVec & src_rowindx, IntNumVec & src_to_tgt_offset
, Int  pivot_fr)
#else
  template <typename T> inline void SupernodalMatrix<T>::UpdateSuperNode(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode, Int &pivot_idx, Int  pivot_fr)
#endif
{

    TIMER_START(UPDATE_SNODE);

#ifdef SINGLE_BLAS

    TIMER_START(UPDATE_SNODE_FIND_INDEX);
    Int first_pivot_idx = -1;
    Int tgt_fc = pivot_fr;
    if(tgt_fc ==I_ZERO ){
#ifdef FAST_INDEX_SEARCH
      Int tgt_fc = tgt_snode.FirstCol();
      first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
      if(first_pivot_idx<0){
        tgt_fc = -first_pivot_idx;
      }
#else
      tgt_fc = tgt_snode.FirstCol();
      //find the pivot idx
      do {first_pivot_idx = src_snode.FindBlockIdx(tgt_fc); tgt_fc++;}
      while(first_pivot_idx<0 && tgt_fc<=tgt_snode.LastCol());
      tgt_fc--;
#endif
    }
    else{
      first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
    }
    NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);

#ifdef FAST_INDEX_SEARCH
    Int tgt_lc = tgt_snode.LastCol();
    Int last_pivot_idx = src_snode.FindBlockIdx(tgt_lc);
    if(last_pivot_idx<0){
      if(last_pivot_idx == -(iSize_+1)){
        last_pivot_idx = src_snode.NZBlockCnt()-1;
      }
      else{
        last_pivot_idx = src_snode.FindBlockIdx(-last_pivot_idx)-1;
      }
      assert(last_pivot_idx>=0);
    }
    NZBlockDesc & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);
    tgt_lc = min(tgt_lc,last_pivot_desc.GIndex + src_snode.NRows(last_pivot_idx)-1);
#else
    Int tgt_lc = tgt_snode.LastCol();
    Int last_pivot_idx = -1;
    //find the pivot idx
    do {last_pivot_idx = src_snode.FindBlockIdx(tgt_lc); tgt_lc--;}
    while(last_pivot_idx<0 && tgt_lc>=tgt_fc);
    tgt_lc++;
    NZBlockDesc & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);
#endif

    TIMER_STOP(UPDATE_SNODE_FIND_INDEX);

    //determine the first column that will be updated in the target supernode
    Int tgt_local_fc =  tgt_fc - tgt_snode.FirstCol();
    Int tgt_local_lc =  tgt_lc - tgt_snode.FirstCol();

    Int src_nrows = src_snode.NRowsBelowBlock(first_pivot_idx) - (tgt_fc - first_pivot_desc.GIndex);
    Int src_lr = tgt_fc+src_nrows-1;
    src_nrows = src_lr - tgt_fc + 1;
  

#ifdef ZEROCOLUMNS

    Int src_snode_size = src_snode.Size();
    Int tgt_snode_size = tgt_snode.Size();
    Int tgt_width = tgt_lc - tgt_fc+1;

    Int tgt_nrows = tgt_snode.NRowsBelowBlock(0);

    T * pivot = &(src_snode.GetNZval(first_pivot_desc.Offset)[(tgt_fc-first_pivot_desc.GIndex)*src_snode.Size()]);
    T * tgt = tgt_snode.GetNZval(0);
    if(src_nrows == tgt_nrows){
      //everything is in row-major
      TIMER_START(UPDATE_SNODE_GEMM);
      blas::Gemm('T','N',tgt_width, src_nrows,src_snode.Size(),MINUS_ONE<T>(),pivot,src_snode.Size(),pivot,src_snode.Size(),ONE<T>(),tgt,tgt_width);
      TIMER_STOP(UPDATE_SNODE_GEMM);
    }
    else{
      //Compute the update in a temporary buffer
      //we need to put the pivot in a "dense" buffer with zeros
      T * pivot_zero = tmpBuf.Data();
      Int pivot_nnz = src_snode_size*tgt_width;
      T * buf = tmpBuf.Data() + pivot_nnz;
      
      if(tgt_width>1){

//        if(first_pivot_idx != last_pivot_idx){gdb_lock();}

      std::fill(pivot_zero,pivot_zero+pivot_nnz,ZERO<T>());
//      if(first_pivot_idx != last_pivot_idx){
//        logfileptr->OFS()<<"INIT Pivot with zeros"<<endl;
//        for(Int row=0;row<tgt_width;++row){
//          for(Int col=0;col<src_snode_size;++col){
//            logfileptr->OFS()<<" "<<pivot_zero[row*src_snode_size+col];
//          }
//          logfileptr->OFS()<<endl;
//        }
//      }

        Int offset = 0;
        Int prevRow = max(tgt_fc, first_pivot_desc.GIndex);
        for(Int blkidx = first_pivot_idx; blkidx<= last_pivot_idx;++blkidx){
          NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
          Int cur_src_nrows = src_snode.NRows(blkidx);
          Int cur_src_lr = min(tgt_lc,cur_block_desc.GIndex + cur_src_nrows -1);
          Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
          cur_src_nrows = cur_src_lr - cur_src_fr +1;
          T * pivot_src = src_snode.GetNZval(cur_block_desc.Offset)+
            (cur_src_fr - cur_block_desc.GIndex)*src_snode.Size();

          offset += (cur_src_fr - prevRow)*src_snode_size ;
          std::copy(pivot_src,pivot_src+cur_src_nrows*src_snode.Size(),pivot_zero+offset);
          prevRow = cur_src_fr;
        }
      }
      else{
        pivot_zero = pivot;
      }

      //Do the Gemm
      TIMER_START(UPDATE_SNODE_GEMM);
      blas::Gemm('T','N',tgt_width, src_nrows,src_snode.Size(),MINUS_ONE<T>(),pivot_zero,src_snode.Size(),pivot,src_snode.Size(),ZERO<T>(),buf,tgt_width);
      TIMER_STOP(UPDATE_SNODE_GEMM);

      //Now do the assembly updating a CONTIGUOUS set of columns
      if(tgt_snode.Id()==14304 || tgt_snode.Id()==12872){
        logfileptr->OFS()<<"Pivot with zeros"<<endl;
        for(Int row=0;row<tgt_width;++row){
          for(Int col=0;col<src_snode_size;++col){
            logfileptr->OFS()<<" "<<pivot_zero[row*src_snode_size+col];
          }
          logfileptr->OFS()<<endl;
        }
        logfileptr->OFS()<<"Update with zeros"<<endl;
        for(Int row=0;row<src_nrows;++row){
          for(Int col=0;col<tgt_width;++col){
            logfileptr->OFS()<<" "<<buf[row*tgt_width+col];
          }
          logfileptr->OFS()<<endl;
        }
          logfileptr->OFS()<<endl;
      }

      //Index mapping
      src_to_tgt_offset.Resize(src_nrows);
      //SetValue(src_to_tgt_offset,-1);
      Int rowidx = 0;
      for(Int blkidx = first_pivot_idx; blkidx < src_snode.NZBlockCnt(); ++blkidx){
        NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
        Int cur_src_nrows = src_snode.NRows(blkidx);
        Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
        Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
        cur_src_nrows = cur_src_lr - cur_src_fr +1;

        //Except for the last pivot block which MIGHT be splitted onto multiple blocks in the target
        //The other one MUST reside into a single block in the target
        Int row = cur_src_fr;
        while(row<=cur_src_lr){
          Int tgt_blk_idx = tgt_snode.FindBlockIdx(row);
          NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
          Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + tgt_snode.NRows(tgt_blk_idx)-1);
          Int tgtOffset = cur_tgt_desc.Offset + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
          for(Int cr = row ;cr<=lr;++cr){
            //assert(tgtOffset + (cr - row)*tgt_snode_size <= tgt_nrows*tgt_snode_size);
            src_to_tgt_offset[rowidx] = tgtOffset + (cr - row)*tgt_snode_size;
            rowidx++;
          }
          row += (lr-row+1);
        }
      }

//      logfileptr->OFS()<<"rowidx "<<rowidx<<" vs "<<src_nrows<<std::endl;
//      assert(rowidx==src_nrows);

//      logfileptr->OFS()<<"Index map tgt :"<<src_to_tgt_offset<<std::endl;

      Int tgt_offset = (tgt_fc - tgt_snode.FirstCol());
      for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
        blas::Axpy(tgt_width,ONE<T>(),&buf[rowidx*tgt_width],1,&tgt[src_to_tgt_offset[rowidx] + tgt_offset],1);

//                      for(Int colidx = 0; colidx< tgt_width;++colidx){
//                        Int tgt_colidx = tgt_offset + colidx;
//                        tgt[src_to_tgt_offset[rowidx] + tgt_colidx] += buf[rowidx*tgt_width+colidx]; 
//                      }
      }


    }
#else
    Int tgt_width = src_snode.NRowsBelowBlock(first_pivot_idx) - (tgt_fc - first_pivot_desc.GIndex) - src_snode.NRowsBelowBlock(last_pivot_idx) + (tgt_lc - last_pivot_desc.GIndex)+1;


    Int tgt_nrows = tgt_snode.NRowsBelowBlock(0);

    T * pivot = &(src_snode.GetNZval(first_pivot_desc.Offset)[(tgt_fc-first_pivot_desc.GIndex)*src_snode.Size()]);
    T * tgt = tgt_snode.GetNZval(0);

    //If the target supernode has the same structure,
    //The GEMM is directly done in place
    if(src_nrows == tgt_nrows){
    //if(src_snode.NRowsBelowBlock(last_pivot_idx+1) == tgt_snode.NRowsBelowBlock(1) && first_pivot_idx == last_pivot_idx){

      Int tgt_offset = (tgt_fc - tgt_snode.FirstCol());
      tgt = &tgt[tgt_offset];

      //everything is in row-major
      TIMER_START(UPDATE_SNODE_GEMM);
      blas::Gemm('T','N',tgt_width, src_nrows,src_snode.Size(),MINUS_ONE<T>(),pivot,src_snode.Size(),pivot,src_snode.Size(),ONE<T>(),tgt,tgt_width);
      TIMER_STOP(UPDATE_SNODE_GEMM);


    }
    else{
      //Compute the update in a temporary buffer
#ifdef _DEBUG_
      tmpBuf.Resize(tgt_width,src_nrows);
#endif

      T * buf = tmpBuf.Data();

      //everything is in row-major
      TIMER_START(UPDATE_SNODE_GEMM);
      blas::Gemm('T','N',tgt_width, src_nrows,src_snode.Size(),MINUS_ONE<T>(),pivot,src_snode.Size(),pivot,src_snode.Size(),ZERO<T>(),buf,tgt_width);
      TIMER_STOP(UPDATE_SNODE_GEMM);


#ifdef _DEBUG_
      logfileptr->OFS()<<"tmpBuf is "<<tmpBuf<<std::endl;
#endif

      //now add the update to the target supernode


      if(1){

        TIMER_START(UPDATE_SNODE_INDEX_MAP);
        Int src_snode_size = src_snode.Size();
        Int tgt_snode_size = tgt_snode.Size();



        if(tgt_snode_size==1){

          Int rowidx = 0;
          for(Int blkidx = first_pivot_idx; blkidx < src_snode.NZBlockCnt(); ++blkidx){
            NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
            Int cur_src_nrows = src_snode.NRows(blkidx);
            Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
            Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);

            //NOTE: Except for the last pivot block which MIGHT 
            //      be splitted onto multiple blocks in the target
            //      The others MUST reside into single target block
            Int row = cur_src_fr;
            while(row<=cur_src_lr){
              Int tgt_blk_idx = tgt_snode.FindBlockIdx(row);
              NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
              Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + tgt_snode.NRows(tgt_blk_idx)-1);
              Int tgtOffset = cur_tgt_desc.Offset + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
              for(Int cr = row ;cr<=lr;++cr){
                tgt[tgtOffset + (cr - row)*tgt_snode_size] += buf[rowidx]; 
                rowidx++;
              }
              row += (lr-row+1);
            }
          }

        }
        else{

          //  IntNumVec src_colindx(tgt_width);
          //  IntNumVec src_rowindx(src_nrows);
          //  IntNumVec src_to_tgt_offset(src_nrows);

          src_colindx.Resize(tgt_width);
          src_to_tgt_offset.Resize(src_nrows);

          Int colidx = 0;
          Int rowidx = 0;
          Int offset = 0;


          for(Int blkidx = first_pivot_idx; blkidx < src_snode.NZBlockCnt(); ++blkidx){
            NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
            Int cur_src_nrows = src_snode.NRows(blkidx);
            Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
            Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
            cur_src_nrows = cur_src_lr - cur_src_fr +1;

            //Except for the last pivot block which MIGHT be splitted onto multiple blocks in the target
            //The other one MUST reside into a single block in the target
            //    if(blkidx==last_pivot_idx){
            Int row = cur_src_fr;
            while(row<=cur_src_lr){
              Int tgt_blk_idx = tgt_snode.FindBlockIdx(row);
              NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
              Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + tgt_snode.NRows(tgt_blk_idx)-1);
              Int tgtOffset = cur_tgt_desc.Offset + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
              for(Int cr = row ;cr<=lr;++cr){
                if(cr<=tgt_lc){
                  src_colindx[colidx++] = cr;
                }

                offset+=tgt_width;

                src_to_tgt_offset[rowidx] = tgtOffset + (cr - row)*tgt_snode_size;
                rowidx++;
              }
              row += (lr-row+1);
            }
          }


          ////  for(Int blkidx = first_pivot_idx; blkidx < src_snode.NZBlockCnt(); ++blkidx){
          ////    NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
          ////    Int cur_src_nrows = src_snode.NRows(blkidx);
          ////    Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
          ////    Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
          ////    cur_src_nrows = cur_src_lr - cur_src_fr +1;
          ////
          ////    for(Int row = cur_src_fr; row<= cur_src_lr;++row){
          ////      if(row<=tgt_lc){
          ////        src_colindx[colidx++] = row;
          ////      }
          ////      
          ////      src_rowindx[rowidx] = row;
          ////      offset+=tgt_width;
          ////
          ////      Int tgt_blk_idx = tgt_snode.FindBlockIdx(row);
          ////      NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
          ////      src_to_tgt_offset[rowidx] = cur_tgt_desc.Offset + (row - cur_tgt_desc.GIndex)*tgt_snode_size; 
          ////      rowidx++;
          ////    }
          ////  }




          //Multiple cases to consider
          // same structure between src and tgt : src_nrows == tgt_nrows
          // tgt has only one column 
          // single pivot block first_pivot idx == last_pivot_idx updating contiguous columns
          // full sparse case (done right now)







          TIMER_STOP(UPDATE_SNODE_INDEX_MAP);
          /////
          /////#ifdef _DEBUG_ 
          /////logfileptr->OFS()<<"src_rowindx :"<<src_rowindx<<std::endl;
          /////logfileptr->OFS()<<"src_colindx :"<<src_colindx<<std::endl;
          /////logfileptr->OFS()<<"Index map tgt :"<<src_to_tgt_offset<<std::endl;
          /////logfileptr->OFS()<<"Index map src :"<<src_offset<<std::endl;
          /////#endif
          /////
          /////    TIMER_START(UPDATE_SNODE_ASSEMBLY);

            if(first_pivot_idx==last_pivot_idx){
              Int tgt_offset = (tgt_fc - tgt_snode.FirstCol());
              for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
                //      for(Int colidx = 0; colidx< tgt_width;++colidx){
                //        Int tgt_colidx = tgt_offset + colidx;
                //        tgt[src_to_tgt_offset[rowidx] + tgt_colidx] += buf[rowidx*tgt_width+colidx]; 
                //      }

                blas::Axpy(tgt_width,ONE<T>(),&buf[rowidx*tgt_width],1,&tgt[src_to_tgt_offset[rowidx] + tgt_offset],1);
              }

            }
            else{
              for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
                for(Int colidx = 0; colidx< src_colindx.m();++colidx){
                  Int col = src_colindx[colidx];
                  Int tgt_colidx = col - tgt_snode.FirstCol();
                  tgt[src_to_tgt_offset[rowidx] + tgt_colidx] += buf[rowidx*tgt_width+colidx]; 
                }
              }
            }

          /////    TIMER_STOP(UPDATE_SNODE_ASSEMBLY);
          ///////logfileptr->OFS()<<"After "<<std::endl<<tgt_snode<<std::endl;
        }
      }
      else{






        //logfileptr->OFS()<<"TmpBuf is "<<tmpBuf<<std::endl; 

        //Now we can add the content into tgt_snode taking care of the indices
        Int cur_local_fc = 0;
        for(Int src_col_blk_idx = first_pivot_idx; src_col_blk_idx <= last_pivot_idx; ++src_col_blk_idx){
          Int Offset = 0;//cur_local_fc*tgt_width;

          NZBlockDesc & cur_col_desc = src_snode.GetNZBlockDesc(src_col_blk_idx);
          Int cur_updated_fc = max(tgt_fc,cur_col_desc.GIndex);
          Int cur_updated_lc = min(tgt_lc,cur_col_desc.GIndex+src_snode.NRows(src_col_blk_idx)-1);
          Int cur_updated_width = cur_updated_lc - cur_updated_fc +1;

          for(Int blkidx = first_pivot_idx; blkidx < src_snode.NZBlockCnt(); ++blkidx){
            NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
            Int cur_src_nrows = src_snode.NRows(blkidx);
            Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
            Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
            cur_src_nrows = cur_src_lr - cur_src_fr +1;

            Int tgt_blk_idx = tgt_snode.FindBlockIdx(cur_src_fr);
            Int last_tgt_blk_idx = tgt_snode.FindBlockIdx(cur_src_lr);

            //    assert(tgt_blk_idx != -1);
            //    assert(last_tgt_blk_idx != -1);

            for( ; tgt_blk_idx <= last_tgt_blk_idx; ++tgt_blk_idx){
              NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
              Int cur_tgt_nrows = tgt_snode.NRows(tgt_blk_idx);
              Int cur_tgt_fr = max(cur_src_fr, cur_tgt_desc.GIndex);
              Int cur_tgt_lr = min(cur_tgt_desc.GIndex + cur_tgt_nrows -1,cur_src_lr);
              cur_tgt_nrows = cur_tgt_lr - cur_tgt_fr +1;

              //        Int updated_nrows = min(tgt_nrows,cur_nrows); 


              TIMER_START(UPDATE_SNODE_AXPY);
              T * cur_src = &buf[ Offset + cur_local_fc];
              T * cur_tgt = &tgt_snode.GetNZval(cur_tgt_desc.Offset)[(cur_tgt_fr - cur_tgt_desc.GIndex )*tgt_snode.Size() + cur_updated_fc - tgt_snode.FirstCol() ];

#pragma loop unroll 
              for(Int row = 0; row< cur_tgt_nrows; ++row){
                //        blas::Axpy(cur_updated_width,ONE<T>(),&cur_src[row*tgt_width],1,&cur_tgt[row*tgt_snode.Size() ],1);
                for(Int col = 0; col< cur_updated_width; ++col){
                  cur_tgt[row*tgt_snode.Size() + col ] += cur_src[row*tgt_width + col];
                }
              }
              Offset += cur_tgt_nrows*tgt_width;
              TIMER_STOP(UPDATE_SNODE_AXPY);


              //      logfileptr->OFS()<<"After blkidx "<<blkidx<<std::endl<<tgt_snode<<std::endl;

            }
          }

          cur_local_fc +=cur_updated_width;
        }

      }

    }
#endif




#else
    NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(pivot_idx);
    Int first_pivot_fr = pivot_fr;
      if(first_pivot_fr ==I_ZERO ){
        first_pivot_fr = first_pivot_desc.GIndex;
      }

    //start with the first pivot
    for(int cur_piv_idx=pivot_idx;cur_piv_idx<src_snode.NZBlockCnt();++cur_piv_idx){


      NZBlockDesc & pivot_desc = src_snode.GetNZBlockDesc(cur_piv_idx);

      if(pivot_fr ==I_ZERO || cur_piv_idx != pivot_idx){
        pivot_fr = pivot_desc.GIndex;
      }

#ifdef _DEBUG_
      assert(pivot_fr >= pivot_desc.GIndex);
#endif

      if(pivot_fr>tgt_snode.LastCol()){
        break;
      }

      Int pivot_nrows = src_snode.NRows(cur_piv_idx);
      Int pivot_lr = min(pivot_desc.GIndex + pivot_nrows -1, tgt_snode.LastCol());
      T * pivot = &(src_snode.GetNZval(pivot_desc.Offset)[(pivot_fr-pivot_desc.GIndex)*src_snode.Size()]);

      //determine the first column that will be updated in the target supernode
      Int tgt_updated_fc =  pivot_fr - tgt_snode.FirstCol();
      Int tgt_updated_lc =  pivot_lr - tgt_snode.FirstCol();

      for(int src_idx=pivot_idx;src_idx<src_snode.NZBlockCnt();++src_idx){

        NZBlockDesc & src_desc = src_snode.GetNZBlockDesc(src_idx);

        Int src_fr = max(first_pivot_fr, src_desc.GIndex) ;
        //Int src_fr = src_desc.GIndex;
        Int src_nrows = src_snode.NRows(src_idx);
        Int src_lr = src_desc.GIndex+src_nrows-1;

        do{
          //TODO Need to be replaced by GlobToLoc index
          Int tgt_idx = tgt_snode.FindBlockIdx(src_fr);

          if(tgt_idx<0){
            break;
          }

#ifdef _DEBUG_
          assert(tgt_idx!=-1 && tgt_idx<tgt_snode.NZBlockCnt());
#endif

          NZBlockDesc & tgt_desc = tgt_snode.GetNZBlockDesc(tgt_idx);
          Int tgt_nrows = tgt_snode.NRows(tgt_idx);
          Int tgt_fr = tgt_desc.GIndex;
          Int tgt_lr = tgt_desc.GIndex+tgt_nrows-1;

          Int update_fr = max(tgt_fr, src_fr);
          Int update_lr = min(tgt_lr, src_lr);

#ifdef _DEBUG_
          assert(update_fr >= tgt_fr); 
          assert(update_lr >= update_fr); 
          assert(update_lr <= tgt_lr); 
#endif


#ifdef _DEBUG_
          logfileptr->OFS()<<"L("<<update_fr<<".."<<update_lr<<","<<tgt_snode.FirstCol()+tgt_updated_fc<<".."<< tgt_snode.FirstCol()+tgt_updated_lc <<") -= L("<<update_fr<<".."<<update_lr<<",:) * L("<<pivot_fr<<".."<<pivot_lr<<",:)'"<<endl;
#endif



          //Update tgt_nzblk with src_nzblk

          T * src = &(src_snode.GetNZval(src_desc.Offset)[(update_fr - src_desc.GIndex)*src_snode.Size()]);
          T * tgt = &(tgt_snode.GetNZval(tgt_desc.Offset)[(update_fr - tgt_desc.GIndex)*tgt_snode.Size()+tgt_updated_fc]);

          //everything is in row-major
    TIMER_START(UPDATE_SNODE_GEMM);
          blas::Gemm('T','N',pivot_lr-pivot_fr+1, update_lr - update_fr + 1,src_snode.Size(),MINUS_ONE<T>(),pivot,src_snode.Size(),src,src_snode.Size(),ONE<T>(),tgt,tgt_snode.Size());
    TIMER_STOP(UPDATE_SNODE_GEMM);

          if(tgt_idx+1<tgt_snode.NZBlockCnt()){
            NZBlockDesc & next_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_idx+1);
            src_fr = next_tgt_desc.GIndex;
          }
          else{
            break;
          }
        }while(src_fr<=src_lr);
      }
    } //end for pivots

#endif
    TIMER_STOP(UPDATE_SNODE);
  }


  template <typename T> inline void SupernodalMatrix<T>::AggregateSuperNode(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode, Int &pivot_idx, Int  pivot_fr)
{

    TIMER_START(AGGREGATE_SNODE);

#ifdef SINGLE_BLAS

    TIMER_START(AGGREGATE_SNODE_FIND_INDEX);
    Int first_pivot_idx = -1;
    Int tgt_fc = pivot_fr;
    if(tgt_fc ==I_ZERO ){
#ifdef FAST_INDEX_SEARCH
      Int tgt_fc = tgt_snode.FirstCol();
      first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
      if(first_pivot_idx<0){
        tgt_fc = -first_pivot_idx;
      }
#else
      tgt_fc = tgt_snode.FirstCol();
      //find the pivot idx
      do {first_pivot_idx = src_snode.FindBlockIdx(tgt_fc); tgt_fc++;}
      while(first_pivot_idx<0 && tgt_fc<=tgt_snode.LastCol());
      tgt_fc--;
#endif
    }
    else{
      first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
    }
    NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);

#ifdef FAST_INDEX_SEARCH
//    TIMER_START(UPDATE_SNODE_FIND_INDEX_LAST2);
    Int tgt_lc = tgt_snode.LastCol();
    Int last_pivot_idx = src_snode.FindBlockIdx(tgt_lc);
    if(last_pivot_idx<0){
      if(last_pivot_idx == -(iSize_+1)){
        last_pivot_idx = src_snode.NZBlockCnt()-1;
      }
      else{
        last_pivot_idx = src_snode.FindBlockIdx(-last_pivot_idx)-1;
      }
      assert(last_pivot_idx>=0);
    }
    NZBlockDesc & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);
    tgt_lc = min(tgt_lc,last_pivot_desc.GIndex + src_snode.NRows(last_pivot_idx)-1);
//    TIMER_STOP(UPDATE_SNODE_FIND_INDEX_LAST2);
#else
//    TIMER_START(UPDATE_SNODE_FIND_INDEX_LAST);
    Int tgt_lc = tgt_snode.LastCol();
    Int last_pivot_idx = -1;
    //find the pivot idx
    do {last_pivot_idx = src_snode.FindBlockIdx(tgt_lc); tgt_lc--;}
    while(last_pivot_idx<0 && tgt_lc>=tgt_fc);
    tgt_lc++;
    NZBlockDesc & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);
//    TIMER_STOP(UPDATE_SNODE_FIND_INDEX_LAST);
#endif

    TIMER_STOP(AGGREGATE_SNODE_FIND_INDEX);

    //determine the first column that will be updated in the target supernode
    Int tgt_local_fc =  tgt_fc - tgt_snode.FirstCol();
    Int tgt_local_lc =  tgt_lc - tgt_snode.FirstCol();

    Int src_nrows = src_snode.NRowsBelowBlock(first_pivot_idx) - (tgt_fc - first_pivot_desc.GIndex);
    Int src_lr = tgt_fc+src_nrows-1;
    src_nrows = src_lr - tgt_fc + 1;
   
    Int tgt_width = src_snode.NRowsBelowBlock(first_pivot_idx) - (tgt_fc - first_pivot_desc.GIndex) - src_snode.NRowsBelowBlock(last_pivot_idx) + (tgt_lc - last_pivot_desc.GIndex)+1;
    //tmpBuf.Resize(tgt_width,src_nrows);

    T * pivot = &(src_snode.GetNZval(first_pivot_desc.Offset)[(tgt_fc-first_pivot_desc.GIndex)*src_snode.Size()]);

if(1){

    TIMER_START(AGGREGATE_SNODE_INDEX_MAP);
  Int src_snode_size = src_snode.Size();
  Int tgt_snode_size = tgt_snode.Size();
  IntNumVec src_colindx(tgt_width);
  IntNumVec src_rowindx(src_nrows);
  IntNumVec src_to_tgt_offset(src_nrows);
  IntNumVec src_offset(src_nrows);
  SetValue(src_offset,-1);
  SetValue(src_to_tgt_offset,-1);

  Int colidx = 0;
  Int rowidx = 0;
  Int offset = 0;
  for(Int blkidx = first_pivot_idx; blkidx < src_snode.NZBlockCnt(); ++blkidx){
    NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
    Int cur_src_nrows = src_snode.NRows(blkidx);
    Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
    Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
    cur_src_nrows = cur_src_lr - cur_src_fr +1;

    for(Int row = cur_src_fr; row<= cur_src_lr;++row){
      if(row<=tgt_lc){
        src_colindx[colidx++] = row;
      }
      src_rowindx[rowidx] = row;
      src_offset[rowidx] = offset;
//cur_block_desc.Offset - (first_pivot_desc.Offset + (tgt_fc - first_pivot_desc.GIndex)*src_snode_size ) + (row - cur_block_desc.GIndex)*src_snode_size;
      offset+=src_snode_size;

      Int tgt_blk_idx = tgt_snode.FindBlockIdx(row);
      NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
      src_to_tgt_offset[rowidx] = cur_tgt_desc.Offset + (row - cur_tgt_desc.GIndex)*tgt_snode_size; 
      rowidx++;
    }
  }
    TIMER_STOP(AGGREGATE_SNODE_INDEX_MAP);

#ifdef _DEBUG_ 
logfileptr->OFS()<<"src_rowindx :"<<src_rowindx<<std::endl;
logfileptr->OFS()<<"src_colindx :"<<src_colindx<<std::endl;
logfileptr->OFS()<<"Index map tgt :"<<src_to_tgt_offset<<std::endl;
logfileptr->OFS()<<"Index map src :"<<src_offset<<std::endl;
#endif

    TIMER_START(AGGREGATE_SNODE_ASSEMBLY);
T* tgt = tgt_snode.GetNZval(0);
for(Int rowidx = 0; rowidx < src_rowindx.m(); ++rowidx){
  Int row = src_rowindx[rowidx];
  for(Int colidx = 0; colidx< src_colindx.m();++colidx){
    Int col = src_colindx[colidx];
    Int tgt_colidx = col - tgt_snode.FirstCol();
    Int src_colidx = col - src_snode.FirstCol();
      tgt[src_to_tgt_offset[rowidx] + tgt_colidx] += pivot[src_offset[rowidx]+src_colidx]; 
  }
}
    TIMER_STOP(AGGREGATE_SNODE_ASSEMBLY);
//logfileptr->OFS()<<"After "<<std::endl<<tgt_snode<<std::endl;



}








#else
    NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(pivot_idx);
    Int first_pivot_fr = pivot_fr;
      if(first_pivot_fr ==I_ZERO ){
        first_pivot_fr = first_pivot_desc.GIndex;
      }

    //start with the first pivot
    for(int cur_piv_idx=pivot_idx;cur_piv_idx<src_snode.NZBlockCnt();++cur_piv_idx){


      NZBlockDesc & pivot_desc = src_snode.GetNZBlockDesc(cur_piv_idx);

      if(pivot_fr ==I_ZERO || cur_piv_idx != pivot_idx){
        pivot_fr = pivot_desc.GIndex;
      }

#ifdef _DEBUG_
      assert(pivot_fr >= pivot_desc.GIndex);
#endif

      if(pivot_fr>tgt_snode.LastCol()){
        break;
      }

      Int pivot_nrows = src_snode.NRows(cur_piv_idx);
      Int pivot_lr = min(pivot_desc.GIndex + pivot_nrows -1, tgt_snode.LastCol());
      T * pivot = &(src_snode.GetNZval(pivot_desc.Offset)[(pivot_fr-pivot_desc.GIndex)*src_snode.Size()]);

      //determine the first column that will be updated in the target supernode
      Int tgt_updated_fc =  pivot_fr - tgt_snode.FirstCol();
      Int tgt_updated_lc =  pivot_lr - tgt_snode.FirstCol();

      for(int src_idx=pivot_idx;src_idx<src_snode.NZBlockCnt();++src_idx){

        NZBlockDesc & src_desc = src_snode.GetNZBlockDesc(src_idx);

        Int src_fr = max(first_pivot_fr, src_desc.GIndex) ;
        //Int src_fr = src_desc.GIndex;
        Int src_nrows = src_snode.NRows(src_idx);
        Int src_lr = src_desc.GIndex+src_nrows-1;

        do{
          //TODO Need to be replaced by GlobToLoc index
          Int tgt_idx = tgt_snode.FindBlockIdx(src_fr);

          if(tgt_idx<0){
            break;
          }

#ifdef _DEBUG_
          assert(tgt_idx!=-1 && tgt_idx<tgt_snode.NZBlockCnt());
#endif

          NZBlockDesc & tgt_desc = tgt_snode.GetNZBlockDesc(tgt_idx);
          Int tgt_nrows = tgt_snode.NRows(tgt_idx);
          Int tgt_fr = tgt_desc.GIndex;
          Int tgt_lr = tgt_desc.GIndex+tgt_nrows-1;

          Int update_fr = max(tgt_fr, src_fr);
          Int update_lr = min(tgt_lr, src_lr);

#ifdef _DEBUG_
          assert(update_fr >= tgt_fr); 
          assert(update_lr >= update_fr); 
          assert(update_lr <= tgt_lr); 
#endif


#ifdef _DEBUG_
          logfileptr->OFS()<<"L("<<update_fr<<".."<<update_lr<<","<<tgt_snode.FirstCol()+tgt_updated_fc<<".."<< tgt_snode.FirstCol()+tgt_updated_lc <<") -= L("<<update_fr<<".."<<update_lr<<",:) * L("<<pivot_fr<<".."<<pivot_lr<<",:)'"<<endl;
#endif



          //Update tgt_nzblk with src_nzblk

          T * src = &(src_snode.GetNZval(src_desc.Offset)[(update_fr - src_desc.GIndex)*src_snode.Size()]);
          T * tgt = &(tgt_snode.GetNZval(tgt_desc.Offset)[(update_fr - tgt_desc.GIndex)*tgt_snode.Size()+tgt_updated_fc]);

          //everything is in row-major
    TIMER_START(UPDATE_SNODE_GEMM);
          blas::Gemm('T','N',pivot_lr-pivot_fr+1, update_lr - update_fr + 1,src_snode.Size(),MINUS_ONE<T>(),pivot,src_snode.Size(),src,src_snode.Size(),ONE<T>(),tgt,tgt_snode.Size());
    TIMER_STOP(UPDATE_SNODE_GEMM);

          if(tgt_idx+1<tgt_snode.NZBlockCnt()){
            NZBlockDesc & next_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_idx+1);
            src_fr = next_tgt_desc.GIndex;
          }
          else{
            break;
          }
        }while(src_fr<=src_lr);
      }
    } //end for pivots

#endif
    TIMER_STOP(AGGREGATE_SNODE);
  }





    template <typename T> void SupernodalMatrix<T>::SendDelayedMessages(Int cur_snode_id, CommList & MsgToSend, AsyncComms & OutgoingSend){
    if(!MsgToSend.empty()){

     CommList::iterator it = MsgToSend.begin();
      while( it != MsgToSend.end()){
        Int src_snode_id = it->src_snode_id;
        Int tgt_id = it->tgt_snode_id;
        Int src_nzblk_idx = it->src_nzblk_idx;
        Int src_first_row = it->src_first_row;

        Int iLocalSrc = (src_snode_id-1) / np +1 ;
        SuperNode<T> & prev_src_snode = *LocalSupernodes_[iLocalSrc -1];

        assert(prev_src_snode.Id()==src_snode_id);

        Int last_local_id = LocalSupernodes_.back()->Id();
        Int iLocalI = (cur_snode_id-1) / np +1 ;
        Int local_snode_id = cur_snode_id>=last_local_id?-1:
                                LocalSupernodes_[iLocalI-1]->Id();

        bool is_sendable = ((cur_snode_id<last_local_id 
                                && tgt_id < local_snode_id)
                                  || (cur_snode_id>=last_local_id));

        if(is_sendable){
          //this can be sent now
          Int tgt_snode_id = tgt_id;

          Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
          if(iTarget != iam){
#ifdef _DEBUG_
            logfileptr->OFS()<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode.Id()<<" to Supernode "<<tgt_snode_id<<endl;
            cout<<"P"<<iam<<" has sent update from Supernode "<<prev_src_snode.Id()<<" to Supernode "<<tgt_snode_id<<endl;
#endif

#ifdef _DEBUG_
            logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated by Supernode "<<prev_src_snode.Id()<<" rows "<<src_first_row/*<<" to "<<src_last_row*/<<" "<<src_nzblk_idx<<std::endl;
#endif

            //Send
            Int tgt_first_col = Xsuper_(tgt_snode_id-1);
            Int tgt_last_col = Xsuper_(tgt_snode_id)-1;
            NZBlockDesc & pivot_desc = prev_src_snode.GetNZBlockDesc(src_nzblk_idx);
            Int local_first_row = src_first_row-pivot_desc.GIndex;
            Int nzblk_cnt = prev_src_snode.NZBlockCnt()-src_nzblk_idx;
            Int nz_cnt = (prev_src_snode.NRowsBelowBlock(src_nzblk_idx)
                - local_first_row )*prev_src_snode.Size();
            assert(nz_cnt>0);

            T * nzval_ptr = prev_src_snode.GetNZval(pivot_desc.Offset
                +local_first_row*prev_src_snode.Size());

            TIMER_START(SEND_MALLOC);
              AddOutgoingComm(OutgoingSend, prev_src_snode.Id(), prev_src_snode.Size(), src_first_row, pivot_desc, nzblk_cnt, nzval_ptr, nz_cnt);
            TIMER_STOP(SEND_MALLOC);

            if(OutgoingSend.size() > maxIsend_){
              TIMER_START(SEND_MPI);
              MPI_Send(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,pComm);
              TIMER_STOP(SEND_MPI);

#ifdef _DEBUG_            
            logfileptr->OFS()<<"DELAYED Sending "<<nz_cnt*sizeof(T)<<" bytes to P"<<iTarget<<std::endl;
            logfileptr->OFS()<<"DELAYED     Send factor "<<prev_src_snode.Id()<<" to node"<<tgt_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
            logfileptr->OFS()<<"DELAYED Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
            logfileptr->OFS()<<"DELAYED Sending "<<OutgoingSend.back()->size()<<" bytes to P"<<iTarget<<std::endl;
#endif
              OutgoingSend.pop_back();

            }
            else{
              TIMER_START(SEND_MPI);
              MPI_Isend(OutgoingSend.back()->front(),OutgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,pComm,&OutgoingSend.back()->Request);
              TIMER_STOP(SEND_MPI);

#ifdef _DEBUG_            
            logfileptr->OFS()<<"DELAYED Sending "<<nz_cnt*sizeof(T)<<" bytes to P"<<iTarget<<std::endl;
            logfileptr->OFS()<<"DELAYED     Send factor "<<prev_src_snode.Id()<<" to node"<<tgt_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
            logfileptr->OFS()<<"DELAYED Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
            logfileptr->OFS()<<"DELAYED Sending "<<OutgoingSend.back()->size()<<" bytes to P"<<iTarget<<std::endl;
#endif
            }






          }
        }

        if(is_sendable){
          //remove from the list
          it = MsgToSend.erase(it);
        }
        else{
          it++;
        }
      }


    }
    }



template <typename T> void SupernodalMatrix<T>::AdvanceOutgoing(AsyncComms & outgoingSend){
  //Check for completion of outgoing communication
  if(!outgoingSend.empty()){
    AsyncComms::iterator it = outgoingSend.begin();
    while(it != outgoingSend.end()){
      int flag = 0;
      MPI_Test(&(*it)->Request,&flag,MPI_STATUS_IGNORE);
      if(flag){
        it = outgoingSend.erase(it);
      }
      else{
        it++;
      }
    }
  }
}



template <typename T> inline AsyncComms::iterator SupernodalMatrix<T>::WaitIncomingFactors(AsyncComms & cur_incomingRecv, MPI_Status & recv_status, AsyncComms & outgoingSend) {
//        vector<MPI_Request> reqs;
//        reqs.reserve(cur_incomingRecv.size());
//      for(AsyncComms::iterator it = cur_incomingRecv.begin(); it!=cur_incomingRecv.end();++it){
//        reqs.push_back((*it)->Request);
//      }
//
//          TIMER_START(IRECV_MPI2);
//
//          AsyncComms::iterator it = cur_incomingRecv.begin();
//          if(cur_incomingRecv.size()==0){
//            it = cur_incomingRecv.end();
//          }
//          else{
//          int index=-1;
//          MPI_Waitany(reqs.size(),&reqs[0],&index,&recv_status);
//          advance(it,index);
//          }
//          TIMER_STOP(IRECV_MPI2);
//          return it;
//



          TIMER_START(IRECV_MPI);
        if(cur_incomingRecv.size()==0){
          TIMER_STOP(IRECV_MPI);
          return cur_incomingRecv.end();
        }
        else{
          Int done = 0;
          while(cur_incomingRecv.size()>0){
            for(AsyncComms::iterator it = cur_incomingRecv.begin(); it!=cur_incomingRecv.end();++it){
              Icomm * curComm = *it;
              MPI_Request req = (curComm->Request);
              //MPI_Test(&(curComm->Request),&done,&recv_status);
              MPI_Test(&req,&done,&recv_status);
              //Test if comm is done
              if(done==1){
          TIMER_STOP(IRECV_MPI);
                return it;
              }
            }

         //   AdvanceOutgoing(outgoingSend);
          }
        }
}





template<typename T> void SupernodalMatrix<T>::AsyncRecvFactors(Int iLocalI, std::vector<AsyncComms> & incomingRecvArr,IntNumVec & FactorsToRecv,IntNumVec & UpdatesToDo){

        for(Int nextLocalI = iLocalI;nextLocalI<=LocalSupernodes_.size();++nextLocalI){
          SuperNode<T> * next_src_snode = LocalSupernodes_[nextLocalI-1];
          AsyncComms & incomingRecv = incomingRecvArr[nextLocalI-1]; 
          
            assert(UpdateCount_[next_src_snode->Id()-1] >= FactorsToRecv[nextLocalI-1]);

          Int IrecvCnt = 0; 
          Int maxRecvCnt = min(UpdatesToDo[nextLocalI-1],FactorsToRecv[nextLocalI-1]);

          for(Int idx =0; idx<maxRecvCnt && incomingRecvCnt_ + IrecvCnt < maxIrecv_;
                 ++idx){
            Int max_bytes = 3*sizeof(Int); 
            //The upper bound must be of the width of the "largest" child
            Int nrows = next_src_snode->NRowsBelowBlock(0);
            Int ncols = UpdateWidth_(next_src_snode->Id()-1);
            Int nz_cnt = nrows * ncols;
            max_bytes += (std::max((Int)ceil(nrows/2)+1,next_src_snode->NZBlockCnt()))*sizeof(NZBlockDesc);
            max_bytes += nz_cnt*sizeof(T);


            assert(UpdateCount_[next_src_snode->Id()-1] >= incomingRecv.size()+1);

            incomingRecv.push_back(new Icomm(max_bytes,MPI_REQUEST_NULL));
            Icomm & Irecv = *incomingRecv.back();
            MPI_Irecv(Irecv.front(),Irecv.size(),MPI_BYTE,MPI_ANY_SOURCE,next_src_snode->Id(),pComm,&Irecv.Request);
            ++IrecvCnt;
          }

//  logfileptr->OFS()<<"Posting "<<IrecvCnt<<" IRECV for Supernode "<<next_src_snode->Id()<<std::endl;
          incomingRecvCnt_+=IrecvCnt;
          FactorsToRecv[nextLocalI-1] = maxRecvCnt - IrecvCnt;

          assert(FactorsToRecv[nextLocalI-1]>=0);

          if( incomingRecvCnt_ >= maxIrecv_){
            break;
          }
        }
#ifdef _DEBUG_
  logfileptr->OFS()<<"FactorsToRecv: "<<FactorsToRecv<<std::endl;
#endif

}












template <typename T> void SupernodalMatrix<T>::FanOut( MPI_Comm & pComm ){

  TIMER_START(FACTORIZATION_FO);

  Real timeSta, timeEnd;
  timeSta =  get_time( );
  

  MPI_Comm_rank(pComm, &iam);
  MPI_Comm_size(pComm, &np);
  IntNumVec UpdatesToDo = UpdateCount_;

  CommList FactorsToSend; 
  AsyncComms outgoingSend;

  std::vector<AsyncComms> incomingRecvArr(LocalSupernodes_.size());
  incomingRecvCnt_ = 0;
  IntNumVec FactorsToRecv(LocalSupernodes_.size());



  std::vector<std::queue<LocalUpdate> > LocalUpdates(LocalSupernodes_.size());

  for(Int i = LocalSupernodes_.size()-1; i>=0; --i){
    SuperNode<T> & cur_snode = *LocalSupernodes_[i];

    Int tgt_snode_id = 0;
    Int src_first_row = 0;
    Int src_last_row = 0;
    Int src_nzblk_idx = 0;
    Int src_next_nzblk_idx = 0;

    while(FindNextUpdate2(cur_snode, tgt_snode_id, src_first_row, src_nzblk_idx, src_last_row, src_next_nzblk_idx)){ 
    //while(FindNextUpdate(cur_snode, src_nzblk_idx, src_first_row, src_last_row, tgt_snode_id)){ 
      Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
      if(iTarget == iam){
        Int iLocalJ = (tgt_snode_id-1) / np +1 ;
        LocalUpdates[iLocalJ-1].push(LocalUpdate(cur_snode.Id(),src_nzblk_idx,src_first_row));

#ifdef _DEBUG_
          logfileptr->OFS()<<"FUTURE LOCAL Supernode "<<tgt_snode_id<<" is going to be updated by Supernode "<<cur_snode.Id()<<" from row "<<src_first_row<<" "<<src_nzblk_idx<<std::endl;
#endif
//        FactorsToRecv[iLocalJ-1]-=UpdatesToDo(cur_snode.Id()-1)+1;
      }
    }

    FactorsToRecv[i] = UpdatesToDo(cur_snode.Id()-1) % np;//UpdatesToDo(cur_snode.Id()-1) - LocalUpdates[i].size();
  }


#ifdef _DEBUG_
  logfileptr->OFS()<<"FactorsToRecv: "<<FactorsToRecv<<std::endl;
#endif



#ifdef UPDATE_LIST
  std::list<SnodeUpdate> updates;
#endif

  Int maxwidth = 0;
  for(Int i = 1; i<Xsuper_.m(); ++i){
    Int width =Xsuper_(i) - Xsuper_(i-1);
    if(width>=maxwidth){
      maxwidth = width;
    }
  }

  NumMat<T> tmpBuf(iSize_,maxwidth);
  IntNumVec src_colindx(maxwidth);
  IntNumVec src_rowindx(iSize_);
  IntNumVec src_to_tgt_offset(iSize_);


  //dummy right looking cholesky factorization
  Int I =1;
  while(I<Xsuper_.m() || !FactorsToSend.empty() || !outgoingSend.empty()){

    //Check for completion of outgoing communication
    AdvanceOutgoing(outgoingSend);

    //process some of the delayed send
    SendDelayedMessages(I,FactorsToSend,outgoingSend);


    if(I<Xsuper_.m()){
      Int src_first_col = Xsuper_(I-1);
      Int src_last_col = Xsuper_(I)-1;
      Int iOwner = Mapping_.Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam ){

#ifdef _DEBUG_
        logfileptr->OFS()<<"Processing Supernode "<<I<<std::endl;
#endif

        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> & src_snode = *LocalSupernodes_[iLocalI -1];


//if(I==1){gdb_lock();}


        //Launch Irecv for subsequent local supernodes if I can
        AsyncRecvFactors(iLocalI,incomingRecvArr,FactorsToRecv,UpdatesToDo);

        //Do all my updates (Local and remote)
        //Local updates
        while(!LocalUpdates[iLocalI-1].empty()){
          LocalUpdate & cur_update = LocalUpdates[iLocalI-1].front();
          Int src_snode_id  = cur_update.src_snode_id;
          Int src_nzblk_idx = cur_update.src_nzblk_idx;
          Int src_first_row = cur_update.src_first_row;
          LocalUpdates[iLocalI-1].pop();

          SuperNode<T> & local_src_snode = *LocalSupernodes_[(src_snode_id-1) / np];

#ifdef _DEBUG_
          logfileptr->OFS()<<"LOCAL Supernode "<<src_snode.Id()<<" is updated by Supernode "<<src_snode_id<<" from row "<<src_first_row<<" "<<src_nzblk_idx<<std::endl;
#endif

#ifdef SINGLE_BLAS
          UpdateSuperNode(local_src_snode,src_snode,src_nzblk_idx,tmpBuf,src_colindx,src_rowindx,src_to_tgt_offset, src_first_row);
#else
          UpdateSuperNode(local_src_snode,src_snode,src_nzblk_idx, src_first_row);
#endif
          //        logfileptr->OFS()<<"After "<<src_snode<<std::endl;

          --UpdatesToDo(I-1);
#ifdef _DEBUG_
          logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left"<<endl;
#endif
        }

        //Remote updates

        //first wait for the Irecv
        AsyncComms & cur_incomingRecv = incomingRecvArr[iLocalI-1];
        MPI_Status recv_status;
        AsyncComms::iterator it = WaitIncomingFactors(cur_incomingRecv, recv_status,outgoingSend);

        while( it != cur_incomingRecv.end() ){
          Icomm * curComm = *it;

          Int bytes_received = 0;
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);

          std::vector<char> & src_blocks = *curComm->pSrcBlocks;
          Int src_snode_id = *(Int*)&src_blocks[0];
          Int src_nzblk_cnt = *(((Int*)&src_blocks[0])+1);
          NZBlockDesc * src_blocks_ptr = 
            reinterpret_cast<NZBlockDesc*>(&src_blocks[2*sizeof(Int)]);
          Int src_nzval_cnt = *(Int*)(src_blocks_ptr + src_nzblk_cnt);
          T * src_nzval_ptr = (T*)((Int*)(src_blocks_ptr + src_nzblk_cnt)+1);

          //Create the dummy supernode for that data
          SuperNode<T> dist_src_snode(src_snode_id,Xsuper_[src_snode_id-1],Xsuper_[src_snode_id]-1, Size(), src_blocks_ptr, src_nzblk_cnt, src_nzval_ptr, src_nzval_cnt);


          //              logfileptr->OFS()<<"IRECV Supernode "<<dist_src_snode.Id()<<std::endl;
#ifdef _DEBUG_
          logfileptr->OFS()<<"IRECV Supernode "<<dist_src_snode.Id()<<std::endl;
#endif

#ifdef UPDATE_LIST
          TIMER_START(UPDATE_ANCESTORS);
          FindUpdates(dist_src_snode,updates);
          //now traverse the list
          for(std::list<SnodeUpdate>::iterator it = updates.begin(); it!=updates.end();it++){
            Int tgt_snode_id = it->tgt_snode_id;
            Int src_first_row = it->src_fr;
            Int src_nzblk_idx = dist_src_snode.FindBlockIdx(src_first_row);
            Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);

            if(iTarget == iam){
#ifdef _DEBUG_
              logfileptr->OFS()<<"IRECV Supernode "<<tgt_snode_id<<" is updated by Supernode "<<dist_src_snode.Id()<<" rows "<<src_first_row<<" to "<<src_last_row<<std::endl;
#endif

              Int iLocalJ = (tgt_snode_id-1) / np +1 ;
              SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];



#ifdef SINGLE_BLAS
              UpdateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, tmpBuf,src_colindx,src_rowindx,src_to_tgt_offset,src_first_row);
#else
              UpdateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, src_first_row);
#endif


              --UpdatesToDo(tgt_snode_id-1);

#ifdef _DEBUG_
              logfileptr->OFS()<<UpdatesToDo(tgt_snode_id-1)<<" updates left for Supernode "<<tgt_snode_id<<endl;
#endif
            }
          }
          TIMER_STOP(UPDATE_ANCESTORS);
#else
          //Update everything I own with that factor
          //Update the ancestors
          Int tgt_snode_id = 0;
          Int src_first_row = 0;
          Int src_last_row = 0;
          Int src_nzblk_idx = 0;

          TIMER_START(UPDATE_ANCESTORS);

          Int src_next_nzblk_idx = 0;
          while(FindNextUpdate2(dist_src_snode, tgt_snode_id, src_first_row, src_nzblk_idx, src_last_row, src_next_nzblk_idx)){ 
//          while(FindNextUpdate(dist_src_snode, src_nzblk_idx, src_first_row, src_last_row, tgt_snode_id)){ 
            Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
            if(iTarget == iam){
#ifdef _DEBUG_
              logfileptr->OFS()<<"IRECV Supernode "<<tgt_snode_id<<" is updated by Supernode "<<dist_src_snode.Id()<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;
#endif

              Int iLocalJ = (tgt_snode_id-1) / np +1 ;
              SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];



#ifdef SINGLE_BLAS
              UpdateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, tmpBuf,src_colindx,src_rowindx,src_to_tgt_offset,src_first_row);
#else
              UpdateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, src_first_row);
#endif

              --UpdatesToDo(tgt_snode_id-1);



#ifdef _DEBUG_
              logfileptr->OFS()<<UpdatesToDo(tgt_snode_id-1)<<" updates left for Supernode "<<tgt_snode_id<<endl;
#endif

            }
          }
          TIMER_STOP(UPDATE_ANCESTORS);
#endif





          //delete the request from the list
          cur_incomingRecv.erase(it);
          --incomingRecvCnt_;

          //              logfileptr->OFS()<<cur_incomingRecv.size()<<" async recv to do for Supernode "<<I<<endl;
          //              logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left for Supernode "<<I<<endl;

          if(UpdatesToDo(src_snode.Id()-1)==0){

            //cancel all requests
            Int flag;
            AsyncComms::iterator it2 = cur_incomingRecv.begin();
            while( it2!=cur_incomingRecv.end()){
              Icomm * curComm = *it2;
              MPI_Status recv_status;
              MPI_Test(&(curComm->Request),&flag,&recv_status);
              //Test if comm is done
              assert(flag==0);

              MPI_Cancel(&(curComm->Request));
              it2 = cur_incomingRecv.erase(it2);
              --incomingRecvCnt_; 


            }

            assert(cur_incomingRecv.size()==0);
          }


          it = WaitIncomingFactors(cur_incomingRecv,recv_status,outgoingSend);
        }




          //AdvanceOutgoing(outgoingSend);


        std::vector<T> src_nzval;
        std::vector<char> src_blocks;

        Int nz_cnt;
        Int max_bytes;
        while(UpdatesToDo(I-1)>0){
#ifdef _DEBUG_
          logfileptr->OFS()<<UpdatesToDo(I-1)<<" updates left"<<endl;
#endif


          TIMER_START(RECV_MALLOC);
          if(src_blocks.size()==0){
            max_bytes = 3*sizeof(Int); 
            //The upper bound must be of the width of the "largest" child
#ifdef _DEBUG_
            logfileptr->OFS()<<"Maximum width is "<<UpdateWidth_(I-1)<<std::endl;
#endif

            Int nrows = src_snode.NRowsBelowBlock(0);
            Int ncols = UpdateWidth_(I-1);
            nz_cnt = nrows * ncols;

            max_bytes += (std::max((Int)ceil(nrows/2)+1,src_snode.NZBlockCnt()))*sizeof(NZBlockDesc);
            max_bytes += nz_cnt*sizeof(T); 

            src_blocks.resize(max_bytes);
          }
          TIMER_STOP(RECV_MALLOC);

          TIMER_START(RECV_MPI);
          MPI_Status recv_status;
          int bytes_received = 0;

#ifdef PROBE_FIRST
          MPI_Probe(MPI_ANY_SOURCE,I,pComm,&recv_status);
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);

          bool doabort = false;
          int prev_size = 0;
          if(src_blocks.size()<bytes_received){
            prev_size = src_blocks.size();
            doabort = true;
            //receive anyway
            src_blocks.resize(bytes_received);
          }
#endif


#ifdef PROBE_FIRST
          MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,recv_status.MPI_SOURCE,I,pComm,&recv_status);
#else
          //receive the index array
          MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,MPI_ANY_SOURCE,I,pComm,&recv_status);
#endif
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);

          Int src_snode_id = *(Int*)&src_blocks[0];
          Int src_nzblk_cnt = *(((Int*)&src_blocks[0])+1);
          NZBlockDesc * src_blocks_ptr = 
            reinterpret_cast<NZBlockDesc*>(&src_blocks[2*sizeof(Int)]);
          Int src_nzval_cnt = *(Int*)(src_blocks_ptr + src_nzblk_cnt);
          T * src_nzval_ptr = (T*)((Int*)(src_blocks_ptr + src_nzblk_cnt)+1);
          TIMER_STOP(RECV_MPI);
          //Create the dummy supernode for that data
          SuperNode<T> dist_src_snode(src_snode_id,Xsuper_[src_snode_id-1],Xsuper_[src_snode_id]-1, Size(), src_blocks_ptr, src_nzblk_cnt, src_nzval_ptr, src_nzval_cnt);

#ifdef PROBE_FIRST
          if(doabort){
            cout<<"We have a problem !!!! on P"<<iam<<"\n";
            gdb_lock();

            abort();
          }
#endif
#ifdef _DEBUG_
          logfileptr->OFS()<<"RECV Supernode "<<dist_src_snode.Id()<<std::endl;
#endif

#ifdef UPDATE_LIST
          TIMER_START(UPDATE_ANCESTORS);
          FindUpdates(dist_src_snode,updates);
          //now traverse the list
          for(std::list<SnodeUpdate>::iterator it = updates.begin(); it!=updates.end();it++){
            Int tgt_snode_id = it->tgt_snode_id;
            Int src_first_row = it->src_fr;
            Int src_nzblk_idx = dist_src_snode.FindBlockIdx(src_first_row);
            Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);

            if(iTarget == iam){
#ifdef _DEBUG_
              logfileptr->OFS()<<"RECV Supernode "<<tgt_snode_id<<" is updated by Supernode "<<dist_src_snode.Id()<<" from row "<<src_first_row<<" "<<src_nzblk_idx<<std::endl;
#endif

              Int iLocalJ = (tgt_snode_id-1) / np +1 ;
              SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];



#ifdef SINGLE_BLAS
              UpdateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, tmpBuf,src_colindx,src_rowindx,src_to_tgt_offset,src_first_row);
#else
              UpdateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, src_first_row);
#endif


              --UpdatesToDo(tgt_snode_id-1);
#ifdef _DEBUG_
              logfileptr->OFS()<<UpdatesToDo(tgt_snode_id-1)<<" updates left for Supernode "<<tgt_snode_id<<endl;
#endif
            }
          }
          TIMER_STOP(UPDATE_ANCESTORS);
#else
          //Update everything I own with that factor
          //Update the ancestors
          Int tgt_snode_id = 0;
          Int src_first_row = 0;
          Int src_last_row = 0;
          Int src_nzblk_idx = 0;


          TIMER_START(UPDATE_ANCESTORS);
          Int src_next_nzblk_idx = 0;
          while(FindNextUpdate2(dist_src_snode, tgt_snode_id, src_first_row, src_nzblk_idx, src_last_row, src_next_nzblk_idx)){ 
//          while(FindNextUpdate(dist_src_snode, src_nzblk_idx, src_first_row, src_last_row, tgt_snode_id)){ 
            Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);
            if(iTarget == iam){
#ifdef _DEBUG_
              logfileptr->OFS()<<"RECV Supernode "<<tgt_snode_id<<" is updated by Supernode "<<dist_src_snode.Id()<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;
#endif

              Int iLocalJ = (tgt_snode_id-1) / np +1 ;
              SuperNode<T> & tgt_snode = *LocalSupernodes_[iLocalJ -1];



#ifdef SINGLE_BLAS
              UpdateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, tmpBuf,src_colindx,src_rowindx,src_to_tgt_offset,src_first_row);
#else
              UpdateSuperNode(dist_src_snode,tgt_snode,src_nzblk_idx, src_first_row);
#endif

              --UpdatesToDo(tgt_snode_id-1);
#ifdef _DEBUG_
              logfileptr->OFS()<<UpdatesToDo(tgt_snode_id-1)<<" updates left for Supernode "<<tgt_snode_id<<endl;
#endif

            

            }
          }
          TIMER_STOP(UPDATE_ANCESTORS);
#endif


        }
        //clear the buffer
        { vector<char>().swap(src_blocks);  }
        { vector<T>().swap(src_nzval);  }


        timeEnd =  get_time( );
#ifdef _DEBUG_
        if(UpdatesToDo(src_snode.Id()-1)!=0){gdb_lock();}
        assert(UpdatesToDo(src_snode.Id()-1)==0);
        logfileptr->OFS()<<"  Factoring Supernode "<<I<<" at "<<timeEnd-timeSta<<" s"<<std::endl;
#endif
//if((I*100/(Xsuper_.m()-1)) % 10 == 0){
//        logfileptr->OFS()<<"  Factoring Supernode "<<I<<" at "<<timeEnd-timeSta<<" s"<<std::endl;
//}

        TIMER_START(FACTOR_PANEL);
        //Factorize Diagonal block
        NZBlockDesc & diag_desc = src_snode.GetNZBlockDesc(0);
        for(Int col = 0; col<src_snode.Size();col+=BLOCKSIZE){
          Int bw = min(BLOCKSIZE,src_snode.Size()-col);
          T * diag_nzval = &src_snode.GetNZval(diag_desc.Offset)[col+col*src_snode.Size()];
          lapack::Potrf( 'U', bw, diag_nzval, src_snode.Size());
          T * nzblk_nzval = &src_snode.GetNZval(diag_desc.Offset)[col+(col+bw)*src_snode.Size()];
          blas::Trsm('L','U','T','N',bw, src_snode.NRowsBelowBlock(0)-(col+bw), ONE<T>(),  diag_nzval, src_snode.Size(), nzblk_nzval, src_snode.Size());

          //update the rest !!! (next blocks columns)
          T * tgt_nzval = &src_snode.GetNZval(diag_desc.Offset)[col+bw+(col+bw)*src_snode.Size()];
          blas::Gemm('T','N',src_snode.Size()-(col+bw), src_snode.NRowsBelowBlock(0)-(col+bw),bw,MINUS_ONE<T>(),nzblk_nzval,src_snode.Size(),nzblk_nzval,src_snode.Size(),ONE<T>(),tgt_nzval,src_snode.Size());
        }


        //          T * diag_nzval = src_snode.GetNZval(diag_desc.Offset);
        //          lapack::Potrf( 'U', src_snode.Size(), diag_nzval, src_snode.Size());
        //        if(src_snode.NZBlockCnt()>1){
        //          NZBlockDesc & nzblk_desc = src_snode.GetNZBlockDesc(1);
        //          T * nzblk_nzval = src_snode.GetNZval(nzblk_desc.Offset);
        //          blas::Trsm('L','U','T','N',src_snode.Size(), src_snode.NRowsBelowBlock(1), ONE<T>(),  diag_nzval, src_snode.Size(), nzblk_nzval, src_snode.Size());
        //        }

#ifdef _DEBUG_
        //        logfileptr->OFS()<<src_snode<<std::endl;
#endif

        TIMER_STOP(FACTOR_PANEL);

        //Send my factor to my ancestors. 
#ifdef _DEBUG_
        IntNumVec is_factor_sent(np);
        SetValue(is_factor_sent,0);
#else
        BolNumVec is_factor_sent(np);
        SetValue(is_factor_sent,false);
#endif

        BolNumVec is_skipped(np);
        SetValue(is_skipped,false);

        Int tgt_snode_id = 0;
        Int src_nzblk_idx = 0;
        Int src_first_row = 0;
        Int src_last_row = 0;


#ifdef UPDATE_LIST
        TIMER_START(FIND_UPDATED_ANCESTORS);
        FindUpdates(src_snode,updates);
        //now traverse the list
        for(std::list<SnodeUpdate>::iterator it = updates.begin(); it!=updates.end();it++){
          Int tgt_snode_id = it->tgt_snode_id;
          Int src_first_row = it->src_fr;
          Int src_nzblk_idx = src_snode.FindBlockIdx(src_first_row);
          Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);

          if(iTarget != iam){

#ifdef _DEBUG_
            logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated by Supernode "<<I<<" rows "<<src_first_row<<std::endl;

            if(is_factor_sent[iTarget]){
              logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated Supernode "<<I<<" in factor "<<is_factor_sent[iTarget]<<std::endl;
            }
#endif

            if(!is_factor_sent[iTarget] && !is_skipped[iTarget] ){

#ifdef DELAY_SNODES
              //need a std::unordered_set to check whether 
              if(iLocalI < LocalSupernodes_.size()){
                if(LocalSupernodes_[iLocalI]->Id()< tgt_snode_id){
                  FactorsToSend.push_back(DelayedComm(src_snode.Id(),tgt_snode_id,src_nzblk_idx,src_first_row));
#ifdef _DEBUG_
                  logfileptr->OFS()<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<tgt_snode_id<<" from row "<<src_first_row<<endl;
                  cout<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<tgt_snode_id<<" from row "<<src_first_row<<endl;
#endif
                  is_skipped[iTarget] = true;
                  continue;
                }
              }
#endif

#ifdef _DEBUG_
              is_factor_sent[iTarget] = I;
#else
              is_factor_sent[iTarget] = true;
#endif

              Int tgt_first_col = Xsuper_(tgt_snode_id-1);
              Int tgt_last_col = Xsuper_(tgt_snode_id)-1;

              //Send
              NZBlockDesc & pivot_desc = src_snode.GetNZBlockDesc(src_nzblk_idx);

              Int local_first_row = src_first_row-pivot_desc.GIndex;
              Int nzblk_cnt = src_snode.NZBlockCnt()-src_nzblk_idx;
              Int nz_cnt = (src_snode.NRowsBelowBlock(src_nzblk_idx) - local_first_row )*src_snode.Size();
              assert(nz_cnt>0);

              T * nzval_ptr = src_snode.GetNZval(pivot_desc.Offset
                  +local_first_row*src_snode.Size());
              TIMER_START(SEND_MALLOC);
              AddOutgoingComm(outgoingSend, src_snode.Id(), src_snode.Size(), src_first_row, pivot_desc, nzblk_cnt, nzval_ptr, nz_cnt);
              TIMER_STOP(SEND_MALLOC);



              if(outgoingSend.size() > maxIsend_){
                TIMER_START(SEND_MPI);
                MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,pComm);
                TIMER_STOP(SEND_MPI);

                outgoingSend.pop_back();
              }
              else{
                TIMER_START(SEND_MPI);
                MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,pComm,&outgoingSend.back()->Request);
                TIMER_STOP(SEND_MPI);
              }



#ifdef _DEBUG_            
              logfileptr->OFS()<<"Sending "<<nz_cnt*sizeof(T)<<" bytes to P"<<iTarget<<std::endl;
              logfileptr->OFS()<<"     Sent factor "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
              logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
              logfileptr->OFS()<<"Sending "<<src_blocks.size()<<" bytes to P"<<iTarget<<std::endl;
#endif

            }
          }
//          else{
//            Int iLocalJ = (tgt_snode_id-1) / np +1 ;
//
//#ifdef _DEBUG_
//            assert(LocalSupernodes_[iLocalJ-1]->Id()==tgt_snode_id);
//#endif
//            LocalUpdates[iLocalJ-1].push(LocalUpdate(src_snode.Id(),src_nzblk_idx,src_first_row));
//          }
        }
        TIMER_STOP(FIND_UPDATED_ANCESTORS);
#else
        TIMER_START(FIND_UPDATED_ANCESTORS);
          Int src_next_nzblk_idx = 0;
          while(FindNextUpdate2(src_snode, tgt_snode_id, src_first_row, src_nzblk_idx, src_last_row, src_next_nzblk_idx)){ 
//        while(FindNextUpdate(src_snode, src_nzblk_idx, src_first_row, src_last_row, tgt_snode_id)){ 
          Int iTarget = Mapping_.Map(tgt_snode_id-1,tgt_snode_id-1);

          if(iTarget != iam){


#ifdef _DEBUG_
            logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated by Supernode "<<I<<" rows "<<src_first_row<<" to "<<src_last_row<<" "<<src_nzblk_idx<<std::endl;

            if(is_factor_sent[iTarget]){
              logfileptr->OFS()<<"Remote Supernode "<<tgt_snode_id<<" is updated Supernode "<<I<<" in factor "<<is_factor_sent[iTarget]<<std::endl;
            }
#endif

            if(!is_factor_sent[iTarget] && !is_skipped[iTarget] ){

#ifdef DELAY_SNODES
              //need a std::unordered_set to check whether 
              if(iLocalI < LocalSupernodes_.size()){
                if(LocalSupernodes_[iLocalI]->Id()< tgt_snode_id){
                  //need to push the prev src_last_row
                  FactorsToSend.push_back(DelayedComm(src_snode.Id(),tgt_snode_id,src_nzblk_idx,src_first_row));
#ifdef _DEBUG_
                  logfileptr->OFS()<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<tgt_snode_id<<" from row "<<src_first_row<<endl;
                  cout<<"P"<<iam<<" has delayed update from Supernode "<<I<<" to "<<tgt_snode_id<<" from row "<<src_first_row<<endl;
#endif
                  is_skipped[iTarget] = true;
                  continue;
                }
              }
#endif



#ifdef _DEBUG_
              is_factor_sent[iTarget] = I;
#else
              is_factor_sent[iTarget] = true;
#endif

              Int tgt_first_col = Xsuper_(tgt_snode_id-1);
              Int tgt_last_col = Xsuper_(tgt_snode_id)-1;

              //Send
              NZBlockDesc & pivot_desc = src_snode.GetNZBlockDesc(src_nzblk_idx);

              Int local_first_row = src_first_row-pivot_desc.GIndex;



              //              logfileptr->OFS()<<src_first_row<<std::endl;
              //              logfileptr->OFS()<<local_first_row<<std::endl;
              //              logfileptr->OFS()<<pivot_desc.GIndex<<std::endl;
              //              assert(src_first_row < pivot_desc.GIndex + src_snode.NRows(src_nzblk_idx));

              Int nzblk_cnt = src_snode.NZBlockCnt()-src_nzblk_idx;

              Int nz_cnt = (src_snode.NRowsBelowBlock(src_nzblk_idx) - local_first_row )*src_snode.Size();
              assert(nz_cnt>0);

              T * nzval_ptr = src_snode.GetNZval(pivot_desc.Offset
                  +local_first_row*src_snode.Size());


              TIMER_START(SEND_MALLOC);
              AddOutgoingComm(outgoingSend, src_snode.Id(), src_snode.Size(), src_first_row, pivot_desc, nzblk_cnt, nzval_ptr, nz_cnt);
              TIMER_STOP(SEND_MALLOC);

              if(outgoingSend.size() > maxIsend_){
                TIMER_START(SEND_MPI);
                MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,pComm);
                TIMER_STOP(SEND_MPI);

#ifdef _DEBUG_            
              logfileptr->OFS()<<"Sending "<<nz_cnt*sizeof(T)<<" bytes to P"<<iTarget<<std::endl;
              logfileptr->OFS()<<"     Sent factor "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
              logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
              logfileptr->OFS()<<"Sending "<<outgoingSend.back()->size()<<" bytes to P"<<iTarget<<std::endl;
#endif
                outgoingSend.pop_back();
              }
              else{
                TIMER_START(SEND_MPI);
                MPI_Isend(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,tgt_snode_id,pComm,&outgoingSend.back()->Request);
                TIMER_STOP(SEND_MPI);

#ifdef _DEBUG_            
              logfileptr->OFS()<<"Sending "<<nz_cnt*sizeof(T)<<" bytes to P"<<iTarget<<std::endl;
              logfileptr->OFS()<<"     Sent factor "<<src_snode.Id()<<" to Supernode "<<tgt_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
              logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
              logfileptr->OFS()<<"Sending "<<outgoingSend.back()->size()<<" bytes to P"<<iTarget<<std::endl;
#endif
              }



            }
          }
//          else{
//            Int iLocalJ = (tgt_snode_id-1) / np +1 ;
//
//            assert(LocalSupernodes_[iLocalJ-1]->Id()==tgt_snode_id);
//
//            LocalUpdates[iLocalJ-1].push(LocalUpdate(src_snode.Id(),src_nzblk_idx,src_first_row));
//          }

        }
        TIMER_STOP(FIND_UPDATED_ANCESTORS);
#endif


      }
      //      MPI_Barrier(pComm);


      //      {
      //      NumMat<T> tmp;
      //      GetFullFactors(tmp,pComm);
      //      }


    }
    I++;
  }


  for(Int idx = 0; idx <incomingRecvArr.size();++idx){
    assert(incomingRecvArr[idx].size()==0);
  } 

  MPI_Barrier(pComm);

  TIMER_STOP(FACTORIZATION_FO);
}








#include "SupernodalMatrix_impl_FB.hpp"











template <typename T> void SupernodalMatrix<T>::Factorize( MPI_Comm & pComm ){
  TIMER_START(FACTORIZATION);
  FanOut(pComm);
  TIMER_STOP(FACTORIZATION);
}




template <typename T> void SupernodalMatrix<T>::GetFullFactors( NumMat<T> & fullMatrix, MPI_Comm &pComm){
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

#ifdef _DEBUG_
        logfileptr->OFS()<<"Supernode "<<I<<"("<<src_snode.Id()<<") is on P"<<iOwner<<" local index is "<<iLocalI<<std::endl; 
        logfileptr->OFS()<<src_snode<<std::endl;
#endif

        NZBlockDesc * nzblk_desc = &src_snode.GetNZBlockDesc(0);
        Int size_blocks = src_snode.NZBlockCnt();
        MPI_Send((void*)&size_blocks,sizeof(Int),MPI_BYTE,0,I,pComm);
        MPI_Send((void*)nzblk_desc,size_blocks*sizeof(NZBlockDesc),MPI_BYTE,0,I,pComm);

        T * nzblk_nzval = src_snode.GetNZval(0);
        Int size_nzval = src_snode.NRowsBelowBlock(0)*src_snode.Size();
        MPI_Send((void*)&size_nzval,sizeof(Int),MPI_BYTE,0,I,pComm);
        MPI_Send((void*)nzblk_nzval,size_nzval*sizeof(T),MPI_BYTE,0,I,pComm);



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

          for(int blkidx=0;blkidx<src_snode.NZBlockCnt();++blkidx){
            NZBlockDesc & nzblk_desc = src_snode.GetNZBlockDesc(blkidx);
            Int nRows = src_snode.NRows(blkidx);
            T * nzblk_nzval = src_snode.GetNZval(nzblk_desc.Offset);

            T * dest = fullMatrix.Data();
            for(Int i = 0; i<nRows;++i){
              for(Int j = 0; j<src_snode.Size();++j){
                if(nzblk_desc.GIndex -1 + i >= src_snode.FirstCol()-1+j){
                  dest[nzblk_desc.GIndex -1 + i + (src_snode.FirstCol()-1+j)*fullMatrix.m()] = nzblk_nzval[i * src_snode.Size() + j];
                }
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
                if(nzblk_desc.GIndex -1 + i >= src_snode.FirstCol()-1+j){
                  dest[nzblk_desc.GIndex -1 + i + (src_snode.FirstCol()-1+j)*fullMatrix.m()] = nzblk_nzval[i * src_snode.Size() + j];
                }
              }
            }
          }
        }
      }
      MPI_Barrier(pComm);
    }

    
      if(iam==0){
#ifdef _DEBUG_
    logfileptr->OFS()<<fullMatrix<<std::endl;
#endif
     } 

      MPI_Barrier(pComm);
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
  
      if(tgt_blkidx<1){
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

#ifdef _CHECK_RESULT_SEQ_
  template <typename T> void SupernodalMatrix<T>::Solve(NumMat<T> * RHS, MPI_Comm & pComm, NumMat<T> & forwardSol, NumMat<T> * Xptr)
#else
  template <typename T> void SupernodalMatrix<T>::Solve(NumMat<T> * RHS, MPI_Comm & pComm, /*NumMat<T> & forwardSol,*/ NumMat<T> * Xptr)
#endif
{
    TIMER_START(SPARSE_SOLVE);

    NumMat<T> & B = *RHS;
    Int nrhs = B.n();

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
    Contributions_.resize(LocalSupernodes_.size());
    std::vector<std::stack<Int> > LocalUpdates(LocalSupernodes_.size());

    AsyncComms outgoingSend;

    //This corresponds to the k loop in dtrsm
    for(Int I=1;I<Xsuper_.m();I++){
      Int iOwner = Mapping_.Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam ){
        //Create the blocks of my contrib with the same nz structure as L
        //and the same width as the final solution
        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];
        SuperNode<T> * contrib = new SuperNode<T>(I,1,nrhs, Size(), cur_snode->NRowsBelowBlock(0) );
        Contributions_[iLocalI-1] = contrib;


        for(Int blkidx = 0; blkidx<cur_snode->NZBlockCnt();++blkidx){
          NZBlockDesc & cur_desc = cur_snode->GetNZBlockDesc(blkidx);
          contrib->AddNZBlock(cur_snode->NRows(blkidx),nrhs,cur_desc.GIndex);
        }

        contrib->Shrink();
      }
    }


    //forward-substitution phase
    //Sending contrib up the tree
    //Start from the leaves of the tree
    TIMER_START(SPARSE_FWD_SUBST);
    for(Int I=1;I<Xsuper_.m();I++){
      Int iOwner = Mapping_.Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam ){
        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];
        Int parent = ETree_.PostParent(cur_snode->LastCol()-1);
        SuperNode<T> * contrib = Contributions_[iLocalI-1];


        //Do all my updates (Local and remote)
        //Local updates
        while(!LocalUpdates[iLocalI-1].empty()){
          Int contrib_snode_id = LocalUpdates[iLocalI-1].top();
          LocalUpdates[iLocalI-1].pop();

          SuperNode<T> * dist_contrib = Contributions_[(contrib_snode_id-1) / np];

#ifdef _DEBUG_
          logfileptr->OFS()<<"LOCAL Supernode "<<I<<" is updated by contrib of Supernode "<<contrib_snode_id<<std::endl;
#endif

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
#ifdef _DEBUG_
          logfileptr->OFS()<<UpdatesToDo(I-1)<<" contribs left"<<endl;
#endif


#ifndef PACKING
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
#else

//    blocklens[0] = sizeof(Int);
//    blocklens[1] = sizeof(Int);
//    blocklens[2] = sizeof(NZBlockDesc);
//    blocklens[3] = (nzblk_cnt_-1)*sizeof(NZBlockDesc);
//    blocklens[4] = sizeof(Int);
//    blocklens[5] = nz_cnt * sizeof(T);

    TIMER_START(RECV_MALLOC);
          if(src_blocks.size()==0){
            max_bytes = 3*sizeof(Int); 
            Int nrows = cur_snode->NRowsBelowBlock(0);
            Int ncols = nrhs;
            nz_cnt = nrows * ncols;

            max_bytes += (std::max((Int)ceil(nrows/2)+1,cur_snode->NZBlockCnt()))*sizeof(NZBlockDesc);
//            max_bytes += (ceil(nrows/2)+1)*sizeof(NZBlockDesc);
            max_bytes += nz_cnt*sizeof(T); 

            src_blocks.resize(max_bytes);
          }
    TIMER_STOP(RECV_MALLOC);

    TIMER_START(RECV_MPI);
          MPI_Status recv_status;
          int bytes_received = 0;

#ifdef PROBE_FIRST
          MPI_Probe(MPI_ANY_SOURCE,I,pComm,&recv_status);
          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);

logfileptr->OFS()<<"Preparing to receive "<<bytes_received<<" bytes"<<endl;

          bool doabort = false;
          int prev_size = 0;
          if(src_blocks.size()<bytes_received){

            cout<<"We have a problem !!!! on P"<<iam<<"\n";
            gdb_lock();
            prev_size = src_blocks.size();
            doabort = true;
            //receive anyway
            src_blocks.resize(bytes_received);
          }
#endif

#ifdef PROBE_FIRST
logfileptr->OFS()<<"Receiving from P"<<recv_status.MPI_SOURCE<<endl;
          MPI_Recv(&src_blocks[0],src_blocks.size(),MPI_BYTE,recv_status.MPI_SOURCE,I,pComm,&recv_status);
#else
          MPI_Recv(&src_blocks[0],src_blocks.size(),MPI_BYTE,MPI_ANY_SOURCE,I,pComm,&recv_status);
#endif


          MPI_Get_count(&recv_status, MPI_BYTE, &bytes_received);

          Int src_snode_id = *(Int*)&src_blocks[0];
          Int src_nzblk_cnt = *(((Int*)&src_blocks[0])+1);
          NZBlockDesc * src_blocks_ptr = 
                    reinterpret_cast<NZBlockDesc*>(&src_blocks[2*sizeof(Int)]);
          Int src_nzval_cnt = *(Int*)(src_blocks_ptr + src_nzblk_cnt);
          T * src_nzval_ptr = (T*)((Int*)(src_blocks_ptr + src_nzblk_cnt)+1);
    TIMER_STOP(RECV_MPI);
#endif

          //Create the dummy supernode for that data
          SuperNode<T> dist_contrib(src_snode_id,1,nrhs, Size(), src_blocks_ptr, src_nzblk_cnt, src_nzval_ptr, src_nzval_cnt);

#ifdef PROBE_FIRST
          if(doabort){

            abort();
          }
#endif

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
                    diag_nzval[kk*nrhs+j] = (B(diag_desc.GIndex-1+kk,j) + diag_nzval[kk*nrhs+j]) / chol_nzval[kk*cur_snode->Size()+kk];
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


#ifndef PACKING
//              std::vector<char> * pNewDesc = new std::vector<char>();
//
//              Int tgt_first_col = Xsuper_(parent_snode_id-1);
//              Int tgt_last_col = Xsuper_(parent_snode_id)-1;
//
//              //Send
//              Int src_nzblk_idx = 1;
//              NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);
//
//              Int src_first_row = pivot_desc.GIndex;
//              Int local_first_row = 0;
//              Int nzblk_cnt = contrib->NZBlockCnt()-src_nzblk_idx;
//              pNewDesc->resize(sizeof(Int) + nzblk_cnt*sizeof(NZBlockDesc));
//
//              char * send_ptr = &(*pNewDesc)[0];
//
//              Int * id_ptr = reinterpret_cast<Int *>(&(*pNewDesc)[0]);
//              NZBlockDesc * block_desc_ptr = 
//                reinterpret_cast<NZBlockDesc *>(&(*pNewDesc)[sizeof(Int)]);
//
//              *id_ptr = contrib->Id();
//              //copy the block descriptors
//              std::copy(&pivot_desc,&pivot_desc+nzblk_cnt, block_desc_ptr);
//              //change the first one
//              block_desc_ptr->Offset += (src_first_row - block_desc_ptr->GIndex)*nrhs;
//              block_desc_ptr->GIndex = src_first_row;
//
//              //send the block descriptors
//              MPI_Send(send_ptr,nzblk_cnt*sizeof(NZBlockDesc) + sizeof(Int),
//                  MPI_BYTE,iTarget,parent_snode_id,pComm);
//
//              T * nzval_ptr = contrib->GetNZval(pivot_desc.Offset
//                  +local_first_row*nrhs);
//
//              Int nz_cnt = (contrib->NRowsBelowBlock(src_nzblk_idx)
//                  - local_first_row )*nrhs;
//
//              MPI_Send(nzval_ptr,nz_cnt*sizeof(T),
//                  MPI_BYTE,iTarget,parent_snode_id,pComm);
//
//              delete pNewDesc;
#else



        Int tgt_first_col = Xsuper_(parent_snode_id-1);
              Int tgt_last_col = Xsuper_(parent_snode_id)-1;
              Int src_nzblk_idx = 1;
              NZBlockDesc & pivot_desc = contrib->GetNZBlockDesc(src_nzblk_idx);

              Int src_first_row = pivot_desc.GIndex;
              Int local_first_row = 0;
              Int nzblk_cnt = contrib->NZBlockCnt()-src_nzblk_idx;

              T * nzval_ptr = contrib->GetNZval(pivot_desc.Offset
                  +local_first_row*nrhs);

              Int nz_cnt = (contrib->NRowsBelowBlock(src_nzblk_idx)
                  - local_first_row )*nrhs;





    TIMER_START(SEND_MALLOC);
              AddOutgoingComm(outgoingSend, contrib->Id(), nrhs, src_first_row, pivot_desc, nzblk_cnt, nzval_ptr, nz_cnt);
    TIMER_STOP(SEND_MALLOC);


    TIMER_START(SEND_MPI);
                  assert(iTarget<np);
    MPI_Send(outgoingSend.back()->front(),outgoingSend.back()->size(), MPI_BYTE,iTarget,parent_snode_id,pComm);
    TIMER_STOP(SEND_MPI);
      outgoingSend.pop_back();
#endif


#ifdef _DEBUG_            
              logfileptr->OFS()<<"     Send contribution "<<I<<" to Supernode "<<parent_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
              logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
              logfileptr->OFS()<<"Sending "<<nzblk_cnt*sizeof(NZBlockDesc)+nz_cnt*sizeof(T) + 3*sizeof(Int)<<" bytes"<< std::endl;
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

#ifdef _CHECK_RESULT_SEQ_
      {
      MPI_Barrier(pComm);
      NumMat<T> tmp = B;
      GetSolution(tmp, pComm);
      Int nrows = 0;
      for(Int ii=1; ii<=I;++ii){ nrows+= Xsuper_[ii] - Xsuper_[ii-1];}

      NumMat<T> tmp2 = tmp;
      blas::Axpy(tmp.m()*tmp.n(),-1.0,&forwardSol(0,0),1,&tmp2(0,0),1);
      double norm = lapack::Lange('F',nrows,tmp.n(),&tmp2(0,0),tmp.m());
      logfileptr->OFS()<<"Norm after SuperNode "<<I<<" is "<<norm<<std::endl; 

        if(abs(norm)>=1e-1){
          for(Int i = 0;i<tmp.m();++i){
            logfileptr->OFS()<<forwardSol(i,0)<<"       "<<tmp(i,0)<<std::endl;
          }
        }
      }
#endif

      //MPI_Barrier(pComm);
    }
    TIMER_STOP(SPARSE_FWD_SUBST);




    //Back-substitution phase

    TIMER_START(SPARSE_BACK_SUBST);

    //start from the root of the tree
    for(Int I=Xsuper_.m()-1;I>=1;--I){
      Int iOwner = Mapping_.Map(I-1,I-1);
      //If I own the column, factor it
      if( iOwner == iam ){
        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> * cur_snode = LocalSupernodes_[iLocalI-1];
        SuperNode<T> * contrib = Contributions_[iLocalI-1];

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

            dist_contrib = Contributions_[(contrib_snode_id-1) / np];
            back_update(dist_contrib,contrib);
          }
          else{

            assert(iTarget!=iam);

            //Receive parent contrib
            std::vector<char> src_blocks;
            std::vector<T> src_nzval;


            //MPI_Recv
            MPI_Status recv_status;

            //Receive the size of the blocks array
            Int max_bytes = 0;
            MPI_Recv(&max_bytes,sizeof(Int),MPI_BYTE,iTarget,I,pComm,&recv_status);
            src_blocks.resize(max_bytes);


            //receive the index array
            MPI_Recv(&src_blocks[0],max_bytes,MPI_BYTE,iTarget,I,pComm,&recv_status);
            //there an aditional integer storing the src id which needs to be removed
            Int src_nzblk_cnt = (max_bytes - sizeof(Int) ) / sizeof(NZBlockDesc);



            //Receive the size of the nzval array
            Int nz_cnt = 0;
            MPI_Recv(&nz_cnt,sizeof(Int),MPI_BYTE,iTarget,I,pComm,&recv_status);
            src_nzval.resize(nz_cnt);

            //receive the nzval array
            MPI_Recv(&src_nzval[0],nz_cnt*sizeof(T),MPI_BYTE,iTarget,I,pComm,&recv_status);

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
                  }
                }
              }
            }

            temp = temp / diag_nzval[ii*cur_snode->Size()+ii];
            tgt_nzval[ii*nrhs+j] = temp;
          }
        }


        //send to my children
        Int colIdx = cur_snode->FirstCol()-1;
        if(colIdx>0){
          Int children_found = 0;
          while(children_found<children(I-1)){
            Int child_snode_id = SupMembership_[colIdx-1];
//            Int firstCol = Xsuper_[child_snode_id-1];
//            for(Int col = colIdx; col>=firstCol; --col){
//            }

            Int parent = ETree_.PostParent(colIdx-1);
            if(parent!=0){
            if(SupMembership_[parent-1]==cur_snode->Id()){
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


                T * nzval_ptr = contrib->GetNZval(pivot_desc.Offset
                    +local_first_row*nrhs);

                Int nz_cnt = (contrib->NRowsBelowBlock(src_nzblk_idx)
                    - local_first_row )*nrhs;

                //send the block descriptors
                  assert(iTarget<np);
                MPI_Send(&bytes_size,sizeof(bytes_size),MPI_BYTE,iTarget,child_snode_id,pComm);
                MPI_Send(send_ptr,bytes_size, MPI_BYTE,iTarget,child_snode_id,pComm);
                //send the nzvals
                MPI_Send(&nz_cnt,sizeof(nz_cnt),MPI_BYTE,iTarget,child_snode_id,pComm);
                MPI_Send(nzval_ptr,nz_cnt*sizeof(T),MPI_BYTE,iTarget,child_snode_id,pComm);

                delete pNewDesc;

#ifdef _DEBUG_            
                logfileptr->OFS()<<"     Send contribution "<<I<<" to Supernode "<<child_snode_id<<" on P"<<iTarget<<" from blk "<<src_nzblk_idx<<std::endl;
                logfileptr->OFS()<<"Sending "<<nzblk_cnt<<" blocks containing "<<nz_cnt<<" nz"<<std::endl;
                logfileptr->OFS()<<"Sending "<<nzblk_cnt*sizeof(NZBlockDesc)<<" and "<<nz_cnt*sizeof(T)<<" bytes during BS"<<std::endl;
#endif
              }
              else{

#ifdef _DEBUG_
                logfileptr->OFS()<<"Local Supernode "<<child_snode_id<<" gets the contribution of Supernode "<<I<<std::endl;
#endif
                Int iLocalJ = (child_snode_id-1) / np +1 ;
                LocalUpdates[iLocalJ-1].push(I);
              }
              children_found++;
            }
            }
            //last column of the prev supernode
            colIdx = Xsuper_[child_snode_id-1]-1;
            if(colIdx==0){
              break;
            }
          }
        }


      }

//      MPI_Barrier(pComm);

#ifdef _CHECK_RESULT_SEQ_
      {
//      MPI_Barrier(pComm);
//      NumMat<T> tmp = B;
//      GetSolution(tmp, pComm);
//            logfileptr->OFS()<<tmp<<std::endl;
//      Int nrows = 0;
//      for(Int ii=1; ii<=I;++ii){ nrows+= Xsuper_[ii] - Xsuper_[ii-1];}
//
//      NumMat<T> tmp2 = tmp;
//      blas::Axpy(tmp.m()*tmp.n(),-1.0,&forwardSol(0,0),1,&tmp2(0,0),1);
//      double norm = lapack::Lange('F',nrows,tmp.n(),&tmp2(0,0),tmp.m());
//      logfileptr->OFS()<<"Norm after SuperNode "<<I<<" is "<<norm<<std::endl; 
//
//        if(abs(norm)>=1e-1){
//          for(Int i = 0;i<tmp.m();++i){
//            logfileptr->OFS()<<forwardSol(i,0)<<"       "<<tmp(i,0)<<std::endl;
//          }
//        }
      }
#endif







    }
    TIMER_STOP(SPARSE_BACK_SUBST);

    MPI_Barrier(pComm);
    TIMER_STOP(SPARSE_SOLVE);

}



template<typename T> void SupernodalMatrix<T>::GetSolution(NumMat<T> & B, MPI_Comm pComm){
    Int nrhs = B.n();
    //Gather B from everybody and put it in the original matrix order
    std::vector<T> tmp_nzval;
    for(Int I=1; I<Xsuper_.m();++I){
      Int iOwner = Mapping_.Map(I-1,I-1);
      T * data;
      Int snode_size = Xsuper_[I] - Xsuper_[I-1];
      Int nzcnt = snode_size * nrhs;
      tmp_nzval.resize(nzcnt);

      if( iOwner == iam ){
        Int iLocalI = (I-1) / np +1 ;
        SuperNode<T> * contrib = Contributions_[iLocalI-1];
        data = contrib->GetNZval(0);
      }
      else{
        data = &tmp_nzval[0];
      }

      MPI_Bcast(data,nzcnt*sizeof(T),MPI_BYTE,iOwner,pComm);

      for(Int i = 0; i<snode_size; ++i){ 
        for(Int j = 0; j<nrhs; ++j){
//          Int destRow = 
          B(Xsuper_[I-1] -1 + i, j) = data[i*nrhs + j];
        }
      }
    }

////
////
////    //Print B
#ifdef _DEBUG_
//    logfileptr->OFS()<<"Solution is "<<B<<std::endl;
#endif
////

      MPI_Barrier(pComm);
}






} // namespace LIBCHOLESKY






#endif 
