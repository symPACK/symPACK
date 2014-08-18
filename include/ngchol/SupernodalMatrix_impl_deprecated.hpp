#ifndef _SUPERNODAL_MATRIX_IMPL_DEPRECATED_HPP_
#define _SUPERNODAL_MATRIX_IMPL_DEPRECATED_HPP_

#if 0
#ifdef UPDATE_LIST
template<typename T> inline void SupernodalMatrix<T>::FindUpdates(SuperNode<T> & src_snode, std::list<SnodeUpdateOld> & updates  ){
    Int iam = CommEnv_->MPI_Rank();
    Int np  = CommEnv_->MPI_Size();
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
        updates.push_back(SnodeUpdateOld(tgt_snode_id,row));
        last_snode_id = tgt_snode_id;
      }
    }
  }
}
#endif

  template <typename T> inline bool SupernodalMatrix<T>::FindNextUpdateOld(SuperNode<T> & src_snode, Int & src_nzblk_idx, Int & src_first_row, Int & src_last_row, Int & tgt_snode_id){
    //src_nzblk_idx is the last nzblock index examined
    //src_first_row is the first row updating the supernode examined
    //src_last_row is the last row updating the supernode examined

    Int iam = CommEnv_->MPI_Rank();
    Int np  = CommEnv_->MPI_Size();
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

////  template <typename T> inline void SupernodalMatrix<T>::UpdateSuperNode(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode, Int &pivot_idx, NumMat<T> & tmpBuf,IntNumVec & src_colindx, IntNumVec & src_rowindx, IntNumVec & src_to_tgt_offset
////, Int  pivot_fr)
////{
////
////    TIMER_START(UPDATE_SNODE);
////
////    TIMER_START(UPDATE_SNODE_FIND_INDEX);
////    Int first_pivot_idx = -1;
////    Int tgt_fc = pivot_fr;
////    if(tgt_fc ==I_ZERO ){
////#ifdef FAST_INDEX_SEARCH
////      Int tgt_fc = tgt_snode.FirstCol();
////      first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
////      if(first_pivot_idx<0){
////        tgt_fc = -first_pivot_idx;
////      }
////#else
////      tgt_fc = tgt_snode.FirstCol();
////      //find the pivot idx
////      do {first_pivot_idx = src_snode.FindBlockIdx(tgt_fc); tgt_fc++;}
////      while(first_pivot_idx<0 && tgt_fc<=tgt_snode.LastCol());
////      tgt_fc--;
////#endif
////    }
////    else{
////      first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
////    }
////    NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);
////
////#ifdef FAST_INDEX_SEARCH
////    Int tgt_lc = tgt_snode.LastCol();
////    Int last_pivot_idx = src_snode.FindBlockIdx(tgt_lc);
////    if(last_pivot_idx<0){
////      if(last_pivot_idx == -(iSize_+1)){
////        last_pivot_idx = src_snode.NZBlockCnt()-1;
////      }
////      else{
////        last_pivot_idx = src_snode.FindBlockIdx(-last_pivot_idx)-1;
////      }
////      assert(last_pivot_idx>=0);
////    }
////    NZBlockDesc & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);
////    tgt_lc = min(tgt_lc,last_pivot_desc.GIndex + src_snode.NRows(last_pivot_idx)-1);
////#else
////    Int tgt_lc = tgt_snode.LastCol();
////    Int last_pivot_idx = -1;
////    //find the pivot idx
////    do {last_pivot_idx = src_snode.FindBlockIdx(tgt_lc); tgt_lc--;}
////    while(last_pivot_idx<0 && tgt_lc>=tgt_fc);
////    tgt_lc++;
////    NZBlockDesc & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);
////#endif
////
////    TIMER_STOP(UPDATE_SNODE_FIND_INDEX);
////
////    //determine the first column that will be updated in the target supernode
////    Int tgt_local_fc =  tgt_fc - tgt_snode.FirstCol();
////    Int tgt_local_lc =  tgt_lc - tgt_snode.FirstCol();
////
////    Int src_nrows = src_snode.NRowsBelowBlock(first_pivot_idx) - (tgt_fc - first_pivot_desc.GIndex);
////    Int src_lr = tgt_fc+src_nrows-1;
////    src_nrows = src_lr - tgt_fc + 1;
////  
////
////#ifdef ZEROCOLUMNS
////{
////    Int src_snode_size = src_snode.Size();
////    Int tgt_snode_size = tgt_snode.Size();
////    Int tgt_width = tgt_lc - tgt_fc+1;
////
////    Int tgt_nrows = tgt_snode.NRowsBelowBlock(0);
////
////    T * pivot = &(src_snode.GetNZval(first_pivot_desc.Offset)[(tgt_fc-first_pivot_desc.GIndex)*src_snode.Size()]);
////    T * tgt = tgt_snode.GetNZval(0);
////    if(src_nrows == tgt_nrows){
////      //everything is in row-major
////      TIMER_START(UPDATE_SNODE_GEMM);
////      blas::Gemm('T','N',tgt_width, src_nrows,src_snode.Size(),MINUS_ONE<T>(),pivot,src_snode.Size(),pivot,src_snode.Size(),ONE<T>(),tgt,tgt_width);
////      TIMER_STOP(UPDATE_SNODE_GEMM);
////    }
////    else{
////      //Compute the update in a temporary buffer
////      //we need to put the pivot in a "dense" buffer with zeros
////      T * pivot_zero = tmpBuf.Data();
////      Int pivot_nnz = src_snode_size*tgt_width;
////      T * buf = tmpBuf.Data() + pivot_nnz;
////      
////      if(tgt_width>1){
////
//////        if(first_pivot_idx != last_pivot_idx){gdb_lock();}
////
////      std::fill(pivot_zero,pivot_zero+pivot_nnz,ZERO<T>());
//////      if(first_pivot_idx != last_pivot_idx){
//////        logfileptr->OFS()<<"INIT Pivot with zeros"<<endl;
//////        for(Int row=0;row<tgt_width;++row){
//////          for(Int col=0;col<src_snode_size;++col){
//////            logfileptr->OFS()<<" "<<pivot_zero[row*src_snode_size+col];
//////          }
//////          logfileptr->OFS()<<endl;
//////        }
//////      }
////
////        Int offset = 0;
////        Int prevRow = max(tgt_fc, first_pivot_desc.GIndex);
////        for(Int blkidx = first_pivot_idx; blkidx<= last_pivot_idx;++blkidx){
////          NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
////          Int cur_src_nrows = src_snode.NRows(blkidx);
////          Int cur_src_lr = min(tgt_lc,cur_block_desc.GIndex + cur_src_nrows -1);
////          Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
////          cur_src_nrows = cur_src_lr - cur_src_fr +1;
////          T * pivot_src = src_snode.GetNZval(cur_block_desc.Offset)+
////            (cur_src_fr - cur_block_desc.GIndex)*src_snode.Size();
////
////          offset += (cur_src_fr - prevRow)*src_snode_size ;
////          std::copy(pivot_src,pivot_src+cur_src_nrows*src_snode.Size(),pivot_zero+offset);
////          prevRow = cur_src_fr;
////        }
////      }
////      else{
////        pivot_zero = pivot;
////      }
////
////      //Do the Gemm
////      TIMER_START(UPDATE_SNODE_GEMM);
////      blas::Gemm('T','N',tgt_width, src_nrows,src_snode.Size(),MINUS_ONE<T>(),pivot_zero,src_snode.Size(),pivot,src_snode.Size(),ZERO<T>(),buf,tgt_width);
////      TIMER_STOP(UPDATE_SNODE_GEMM);
////
////      //Now do the assembly updating a CONTIGUOUS set of columns
////      if(tgt_snode.Id()==14304 || tgt_snode.Id()==12872){
////        logfileptr->OFS()<<"Pivot with zeros"<<endl;
////        for(Int row=0;row<tgt_width;++row){
////          for(Int col=0;col<src_snode_size;++col){
////            logfileptr->OFS()<<" "<<pivot_zero[row*src_snode_size+col];
////          }
////          logfileptr->OFS()<<endl;
////        }
////        logfileptr->OFS()<<"Update with zeros"<<endl;
////        for(Int row=0;row<src_nrows;++row){
////          for(Int col=0;col<tgt_width;++col){
////            logfileptr->OFS()<<" "<<buf[row*tgt_width+col];
////          }
////          logfileptr->OFS()<<endl;
////        }
////          logfileptr->OFS()<<endl;
////      }
////
////      //Index mapping
////      src_to_tgt_offset.Resize(src_nrows);
////      //SetValue(src_to_tgt_offset,-1);
////      Int rowidx = 0;
////      for(Int blkidx = first_pivot_idx; blkidx < src_snode.NZBlockCnt(); ++blkidx){
////        NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
////        Int cur_src_nrows = src_snode.NRows(blkidx);
////        Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
////        Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
////        cur_src_nrows = cur_src_lr - cur_src_fr +1;
////
////        //Except for the last pivot block which MIGHT be splitted onto multiple blocks in the target
////        //The other one MUST reside into a single block in the target
////        Int row = cur_src_fr;
////        while(row<=cur_src_lr){
////          Int tgt_blk_idx = tgt_snode.FindBlockIdx(row);
////          NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
////          Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + tgt_snode.NRows(tgt_blk_idx)-1);
////          Int tgtOffset = cur_tgt_desc.Offset + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
////          for(Int cr = row ;cr<=lr;++cr){
////            //assert(tgtOffset + (cr - row)*tgt_snode_size <= tgt_nrows*tgt_snode_size);
////            src_to_tgt_offset[rowidx] = tgtOffset + (cr - row)*tgt_snode_size;
////            rowidx++;
////          }
////          row += (lr-row+1);
////        }
////      }
////
//////      logfileptr->OFS()<<"rowidx "<<rowidx<<" vs "<<src_nrows<<std::endl;
//////      assert(rowidx==src_nrows);
////
//////      logfileptr->OFS()<<"Index map tgt :"<<src_to_tgt_offset<<std::endl;
////
////      Int tgt_offset = (tgt_fc - tgt_snode.FirstCol());
////      for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
////        blas::Axpy(tgt_width,ONE<T>(),&buf[rowidx*tgt_width],1,&tgt[src_to_tgt_offset[rowidx] + tgt_offset],1);
////
//////                      for(Int colidx = 0; colidx< tgt_width;++colidx){
//////                        Int tgt_colidx = tgt_offset + colidx;
//////                        tgt[src_to_tgt_offset[rowidx] + tgt_colidx] += buf[rowidx*tgt_width+colidx]; 
//////                      }
////      }
////
////
////    }
////}
////#else
////    Int tgt_width = src_snode.NRowsBelowBlock(first_pivot_idx) - (tgt_fc - first_pivot_desc.GIndex) - src_snode.NRowsBelowBlock(last_pivot_idx) + (tgt_lc - last_pivot_desc.GIndex)+1;
////
////
////    Int tgt_nrows = tgt_snode.NRowsBelowBlock(0);
////
////    T * pivot = &(src_snode.GetNZval(first_pivot_desc.Offset)[(tgt_fc-first_pivot_desc.GIndex)*src_snode.Size()]);
////    T * tgt = tgt_snode.GetNZval(0);
////
////    //If the target supernode has the same structure,
////    //The GEMM is directly done in place
////    if(src_nrows == tgt_nrows){
////    //if(src_snode.NRowsBelowBlock(last_pivot_idx+1) == tgt_snode.NRowsBelowBlock(1) && first_pivot_idx == last_pivot_idx){
////
////      Int tgt_offset = (tgt_fc - tgt_snode.FirstCol());
////      tgt = &tgt[tgt_offset];
////
////      //everything is in row-major
////      TIMER_START(UPDATE_SNODE_GEMM);
////      blas::Gemm('T','N',tgt_width, src_nrows,src_snode.Size(),MINUS_ONE<T>(),pivot,src_snode.Size(),pivot,src_snode.Size(),ONE<T>(),tgt,tgt_width);
////      TIMER_STOP(UPDATE_SNODE_GEMM);
////
////
////    }
////    else{
////      //Compute the update in a temporary buffer
////#ifdef _DEBUG_
////      tmpBuf.Resize(tgt_width,src_nrows);
////#endif
////
////      T * buf = tmpBuf.Data();
////
////      //everything is in row-major
////      TIMER_START(UPDATE_SNODE_GEMM);
////      blas::Gemm('T','N',tgt_width, src_nrows,src_snode.Size(),MINUS_ONE<T>(),pivot,src_snode.Size(),pivot,src_snode.Size(),ZERO<T>(),buf,tgt_width);
////      TIMER_STOP(UPDATE_SNODE_GEMM);
////
////
////#ifdef _DEBUG_
////      logfileptr->OFS()<<"tmpBuf is "<<tmpBuf<<std::endl;
////#endif
////
////      //now add the update to the target supernode
////
////
////      if(1){
////
////        TIMER_START(UPDATE_SNODE_INDEX_MAP);
////        Int src_snode_size = src_snode.Size();
////        Int tgt_snode_size = tgt_snode.Size();
////
////
////
////        if(tgt_snode_size==1){
////
////          Int rowidx = 0;
////          for(Int blkidx = first_pivot_idx; blkidx < src_snode.NZBlockCnt(); ++blkidx){
////            NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
////            Int cur_src_nrows = src_snode.NRows(blkidx);
////            Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
////            Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
////
////            //NOTE: Except for the last pivot block which MIGHT 
////            //      be splitted onto multiple blocks in the target
////            //      The others MUST reside into single target block
////            Int row = cur_src_fr;
////            while(row<=cur_src_lr){
////              Int tgt_blk_idx = tgt_snode.FindBlockIdx(row);
////              NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
////              Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + tgt_snode.NRows(tgt_blk_idx)-1);
////              Int tgtOffset = cur_tgt_desc.Offset + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
////              for(Int cr = row ;cr<=lr;++cr){
////                tgt[tgtOffset + (cr - row)*tgt_snode_size] += buf[rowidx]; 
////                rowidx++;
////              }
////              row += (lr-row+1);
////            }
////          }
////
////        }
////        else{
////
////          //  IntNumVec src_colindx(tgt_width);
////          //  IntNumVec src_rowindx(src_nrows);
////          //  IntNumVec src_to_tgt_offset(src_nrows);
////
////          src_colindx.Resize(tgt_width);
////          src_to_tgt_offset.Resize(src_nrows);
////
////          Int colidx = 0;
////          Int rowidx = 0;
////          Int offset = 0;
////
////
////          for(Int blkidx = first_pivot_idx; blkidx < src_snode.NZBlockCnt(); ++blkidx){
////            NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
////            Int cur_src_nrows = src_snode.NRows(blkidx);
////            Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
////            Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
////            cur_src_nrows = cur_src_lr - cur_src_fr +1;
////
////            //Except for the last pivot block which MIGHT be splitted onto multiple blocks in the target
////            //The other one MUST reside into a single block in the target
////            //    if(blkidx==last_pivot_idx){
////            Int row = cur_src_fr;
////            while(row<=cur_src_lr){
////              Int tgt_blk_idx = tgt_snode.FindBlockIdx(row);
////              NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
////              Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + tgt_snode.NRows(tgt_blk_idx)-1);
////              Int tgtOffset = cur_tgt_desc.Offset + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
////              for(Int cr = row ;cr<=lr;++cr){
////                if(cr<=tgt_lc){
////                  src_colindx[colidx++] = cr;
////                }
////
////                offset+=tgt_width;
////
////                src_to_tgt_offset[rowidx] = tgtOffset + (cr - row)*tgt_snode_size;
////                rowidx++;
////              }
////              row += (lr-row+1);
////            }
////          }
////
////
////          ////  for(Int blkidx = first_pivot_idx; blkidx < src_snode.NZBlockCnt(); ++blkidx){
////          ////    NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
////          ////    Int cur_src_nrows = src_snode.NRows(blkidx);
////          ////    Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
////          ////    Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
////          ////    cur_src_nrows = cur_src_lr - cur_src_fr +1;
////          ////
////          ////    for(Int row = cur_src_fr; row<= cur_src_lr;++row){
////          ////      if(row<=tgt_lc){
////          ////        src_colindx[colidx++] = row;
////          ////      }
////          ////      
////          ////      src_rowindx[rowidx] = row;
////          ////      offset+=tgt_width;
////          ////
////          ////      Int tgt_blk_idx = tgt_snode.FindBlockIdx(row);
////          ////      NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
////          ////      src_to_tgt_offset[rowidx] = cur_tgt_desc.Offset + (row - cur_tgt_desc.GIndex)*tgt_snode_size; 
////          ////      rowidx++;
////          ////    }
////          ////  }
////
////
////
////
////          //Multiple cases to consider
////          // same structure between src and tgt : src_nrows == tgt_nrows
////          // tgt has only one column 
////          // single pivot block first_pivot idx == last_pivot_idx updating contiguous columns
////          // full sparse case (done right now)
////
////
////
////
////
////
////
////          TIMER_STOP(UPDATE_SNODE_INDEX_MAP);
////          /////
////          /////#ifdef _DEBUG_ 
////          /////logfileptr->OFS()<<"src_rowindx :"<<src_rowindx<<std::endl;
////          /////logfileptr->OFS()<<"src_colindx :"<<src_colindx<<std::endl;
////          /////logfileptr->OFS()<<"Index map tgt :"<<src_to_tgt_offset<<std::endl;
////          /////logfileptr->OFS()<<"Index map src :"<<src_offset<<std::endl;
////          /////#endif
////          /////
////          /////    TIMER_START(UPDATE_SNODE_ASSEMBLY);
////
////            if(first_pivot_idx==last_pivot_idx){
////              Int tgt_offset = (tgt_fc - tgt_snode.FirstCol());
////              for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
////                //      for(Int colidx = 0; colidx< tgt_width;++colidx){
////                //        Int tgt_colidx = tgt_offset + colidx;
////                //        tgt[src_to_tgt_offset[rowidx] + tgt_colidx] += buf[rowidx*tgt_width+colidx]; 
////                //      }
////
////                blas::Axpy(tgt_width,ONE<T>(),&buf[rowidx*tgt_width],1,&tgt[src_to_tgt_offset[rowidx] + tgt_offset],1);
////              }
////
////            }
////            else{
////              for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
////                for(Int colidx = 0; colidx< src_colindx.m();++colidx){
////                  Int col = src_colindx[colidx];
////                  Int tgt_colidx = col - tgt_snode.FirstCol();
////                  tgt[src_to_tgt_offset[rowidx] + tgt_colidx] += buf[rowidx*tgt_width+colidx]; 
////                }
////              }
////            }
////
////          /////    TIMER_STOP(UPDATE_SNODE_ASSEMBLY);
////          ///////logfileptr->OFS()<<"After "<<std::endl<<tgt_snode<<std::endl;
////        }
////      }
////      else{
////
////
////
////
////
////
////        //logfileptr->OFS()<<"TmpBuf is "<<tmpBuf<<std::endl; 
////
////        //Now we can add the content into tgt_snode taking care of the indices
////        Int cur_local_fc = 0;
////        for(Int src_col_blk_idx = first_pivot_idx; src_col_blk_idx <= last_pivot_idx; ++src_col_blk_idx){
////          Int Offset = 0;//cur_local_fc*tgt_width;
////
////          NZBlockDesc & cur_col_desc = src_snode.GetNZBlockDesc(src_col_blk_idx);
////          Int cur_updated_fc = max(tgt_fc,cur_col_desc.GIndex);
////          Int cur_updated_lc = min(tgt_lc,cur_col_desc.GIndex+src_snode.NRows(src_col_blk_idx)-1);
////          Int cur_updated_width = cur_updated_lc - cur_updated_fc +1;
////
////          for(Int blkidx = first_pivot_idx; blkidx < src_snode.NZBlockCnt(); ++blkidx){
////            NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
////            Int cur_src_nrows = src_snode.NRows(blkidx);
////            Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
////            Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
////            cur_src_nrows = cur_src_lr - cur_src_fr +1;
////
////            Int tgt_blk_idx = tgt_snode.FindBlockIdx(cur_src_fr);
////            Int last_tgt_blk_idx = tgt_snode.FindBlockIdx(cur_src_lr);
////
////            //    assert(tgt_blk_idx != -1);
////            //    assert(last_tgt_blk_idx != -1);
////
////            for( ; tgt_blk_idx <= last_tgt_blk_idx; ++tgt_blk_idx){
////              NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
////              Int cur_tgt_nrows = tgt_snode.NRows(tgt_blk_idx);
////              Int cur_tgt_fr = max(cur_src_fr, cur_tgt_desc.GIndex);
////              Int cur_tgt_lr = min(cur_tgt_desc.GIndex + cur_tgt_nrows -1,cur_src_lr);
////              cur_tgt_nrows = cur_tgt_lr - cur_tgt_fr +1;
////
////              //        Int updated_nrows = min(tgt_nrows,cur_nrows); 
////
////
////              TIMER_START(UPDATE_SNODE_AXPY);
////              T * cur_src = &buf[ Offset + cur_local_fc];
////              T * cur_tgt = &tgt_snode.GetNZval(cur_tgt_desc.Offset)[(cur_tgt_fr - cur_tgt_desc.GIndex )*tgt_snode.Size() + cur_updated_fc - tgt_snode.FirstCol() ];
////
////#pragma loop unroll 
////              for(Int row = 0; row< cur_tgt_nrows; ++row){
////                //        blas::Axpy(cur_updated_width,ONE<T>(),&cur_src[row*tgt_width],1,&cur_tgt[row*tgt_snode.Size() ],1);
////                for(Int col = 0; col< cur_updated_width; ++col){
////                  cur_tgt[row*tgt_snode.Size() + col ] += cur_src[row*tgt_width + col];
////                }
////              }
////              Offset += cur_tgt_nrows*tgt_width;
////              TIMER_STOP(UPDATE_SNODE_AXPY);
////
////
////              //      logfileptr->OFS()<<"After blkidx "<<blkidx<<std::endl<<tgt_snode<<std::endl;
////
////            }
////          }
////
////          cur_local_fc +=cur_updated_width;
////        }
////
////      }
////
////    }
////#endif
////    TIMER_STOP(UPDATE_SNODE);
////  }
#endif

#endif //_SUPERNODAL_MATRIX_IMPL_DEPRECATED_HPP_
