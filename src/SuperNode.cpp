#include "sympack/SuperNode.hpp"

#define SCALAR double
//namespace SYMPACK{
//
//    Int trc_Merge(SuperNode<SCALAR> & src_snode, SnodeUpdate &update,SuperNode<SCALAR> & tgt_snode){
//      TIMER_START(MERGE_SNODE_OBJ);
//
//      assert(tgt_snode.b_own_storage_);
//
//      Int src_snode_size = src_snode.Size();
//      Int tgt_snode_size = tgt_snode.Size();
//
//      Int & pivot_idx = update.blkidx;
//      Int & pivot_fr = update.src_first_row;
//
//////      //find the first row updated by src_snode
//////      TIMER_START(MERGE_SNODE_FIND_INDEX);
//////      Int first_pivot_idx = -1;
//////      Int tgt_fc = pivot_fr;
//////      if(tgt_fc ==I_ZERO ){
//////        tgt_fc = tgt_snode.FirstCol();
//////        //find the pivot idx
//////        do {first_pivot_idx = src_snode.FindBlockIdx(tgt_fc); tgt_fc++;}
//////        while(first_pivot_idx<0 && tgt_fc<=tgt_snode.LastCol());
//////        tgt_fc--;
//////      }
//////      else{
//////        first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
//////      }
//////      assert(first_pivot_idx>=0);
//////
//////      TIMER_STOP(MERGE_SNODE_FIND_INDEX);
//
//      Int tgt_fc;
//      Int first_pivot_idx;
//      tgt_snode.FindUpdatedFirstCol(src_snode, 0, tgt_fc, first_pivot_idx);
//      NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);
//
//      //parse src_snode
//      ITree::Interval overlap;
//      ITree::Interval curInter;
//      ITree::Interval newInter;
//      std::queue< ITree::Interval > toInsert;
//      for(Int blkidx = first_pivot_idx; blkidx<src_snode.NZBlockCnt(); ++blkidx){
//        NZBlockDesc & blk_desc = src_snode.GetNZBlockDesc(blkidx);
//        Int fr = max(tgt_snode.FirstCol(),blk_desc.GIndex);
//        Int lr = blk_desc.GIndex + src_snode.NRows(blkidx) -1;
//      
//        curInter.low = fr;
//        curInter.high = lr;
//        toInsert.push(curInter);
//
//        while(!toInsert.empty()){
//          curInter = toInsert.front();
//          toInsert.pop();
//          if(tgt_snode.FindBlockIdx(curInter.low,curInter.high,overlap)==-1){
//            //Add the full block
//            tgt_snode.AddNZBlock( curInter.high - curInter.low + 1, tgt_snode_size, curInter.low);
//          }
//          else{
//
//            //check the overlap
//            //fr is curInter.low and lr is curInter.high
//            //                l-----overlap------h
//            //            l---------block-------------h
//            //        l--block----h
//            //                              l-----block----h
//            //we have two segments to look for : [overlap.high+1 - lr] and [fr - overlap.low -1]         
//            if(overlap.low>curInter.low){
////gdb_lock();
//              newInter.low = curInter.low;
//              newInter.high = overlap.low-1;
//              toInsert.push(newInter);
//            }
//
//            if(overlap.high < curInter.high){
////gdb_lock();
//              newInter.low = overlap.high+1;
//              newInter.high = curInter.high;
//              toInsert.push(newInter);
//            }
//
//          }
//
//
//        }
//
//
////        if(tgt_snode.FindBlockIdx(fr,lr,overlap)==-1){
////          //Add the full block
////          tgt_snode.AddNZBlock( lr - fr + 1, tgt_snode_size, fr);
////        }
////        else{
////          
////          //check the overlap
////          //                l-----overlap------h
////          //            l---------block-------------h
////          //        l--block----h
////          //                              l-----block----h
////          //we have two segments to look for : [overlap.high+1 - lr] and [fr - overlap.low -1]         
////         
////          if(overlap.high < lr){
////            //we need to add from high+1 to lr 
////            tgt_snode.AddNZBlock( lr - overlap.high, tgt_snode_size, overlap.high+1);
////          }
////          
////          if(overlap.low>fr){
////            //we need to add fr to low-1 
////            tgt_snode.AddNZBlock( overlap.low - fr, tgt_snode_size, fr);
////          }
////        }
//
//
////        Int nrows = src_snode.NRows(blkidx);
////        for(Int rowidx = 0; rowidx<nrows; ++rowidx){
////          Int row = blk_desc.GIndex + rowidx;
////          //if the row is updating the target
////          if(row>=tgt_fc){
////            //check if the row is not already in the structure
////            if(FindBlockIdx(row)==-1){
////              //add a nzblock with a single row in it
////              AddNZBlock( 1, tgt_snode_size, row);
////            }
////          }
////        }
//      }
//
//      TIMER_STOP(MERGE_SNODE_OBJ);
//      return 0;
//    }
//
//
//
//   Int trc_Aggregate(SuperNode<SCALAR> & src_snode,SuperNode<SCALAR> & tgt_snode){
//      TIMER_START(AGGREG_SNODE_OBJ);
//      Int  pivot_idx = 0;
//      Int  pivot_fr = 0;//FirstCol();
//
//      Int src_snode_size = src_snode.Size();
//      Int tgt_snode_size = tgt_snode.Size();
//
/////      TIMER_START(AGGREG_SNODE_FIND_INDEX);
/////      Int first_pivot_idx = -1;
/////      Int tgt_fc = pivot_fr;
/////      if(tgt_fc ==I_ZERO ){
/////        tgt_fc = tgt_snode.FirstCol();
/////        //find the pivot idx
/////        do {first_pivot_idx = src_snode.FindBlockIdx(tgt_fc); tgt_fc++;}
/////        while(first_pivot_idx<0 && tgt_fc<=tgt_snode.LastCol());
/////        tgt_fc--;
/////      }
/////      else{
/////        first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
/////      }
/////      assert(first_pivot_idx>=0);
/////      NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);
/////
/////      TIMER_STOP(AGGREG_SNODE_FIND_INDEX);
//
//
//      Int first_pivot_idx = 0 ;
//      Int tgt_fc = tgt_snode.FirstCol();
//
//      //parse src_snode and add everything
//
//
//      for(Int blkidx = first_pivot_idx; blkidx<src_snode.NZBlockCnt(); ++blkidx){
//        NZBlockDesc & blk_desc = src_snode.GetNZBlockDesc(blkidx);
//        Int nrows = src_snode.NRows(blkidx);
//        for(Int rowidx = 0; rowidx<nrows; ++rowidx){
//          Int row = blk_desc.GIndex + rowidx;
//
//          if(row>=tgt_fc){
//            Int src_offset = blk_desc.Offset + (row - blk_desc.GIndex)*src_snode_size;
//
//            Int tgt_blkidx = tgt_snode.FindBlockIdx(row);
//            assert(tgt_blkidx!=-1);
//            NZBlockDesc & tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blkidx);
//            Int tgt_offset = tgt_desc.Offset + (row - tgt_desc.GIndex)*tgt_snode_size;
//
//            SCALAR * src = src_snode.GetNZval(src_offset);
//            SCALAR * tgt = tgt_snode.GetNZval(tgt_offset);
//
//            //blas::Axpy(tgt_snode_size,ONE<T>(),src,1,tgt,1);
//            for(Int i = 0; i< tgt_snode_size;i+=1){
//              tgt[i] += src[i];
//            }
//
//          }
//        }
//      }
//
//
//      TIMER_STOP(AGGREG_SNODE_OBJ);
//      return 0;
//    }
//
//
//
//
//
//
//
//
//
//  Int trc_UpdateAggregate(SuperNode<SCALAR> & src_snode, SnodeUpdate &update, TempUpdateBuffers<SCALAR> & tmpBuffers, Int iTarget,SuperNode<SCALAR> & tgt_snode){
//
//      if(iTarget != iam){
//
//        tgt_snode.Merge(src_snode, update);
//
//        Int & pivot_idx = update.blkidx;
//        Int & pivot_fr = update.src_first_row;
//
//        TIMER_START(UPDATE_SNODE_OBJ);
//        Int src_snode_size = src_snode.Size();
//        Int tgt_snode_size = tgt_snode.Size();
//
//////        //find the first row updated by src_snode
//////        TIMER_START(UPDATE_SNODE_FIND_INDEX);
//////        Int first_pivot_idx = -1;
//////        Int tgt_fc = pivot_fr;
//////        if(tgt_fc ==I_ZERO ){
//////          tgt_fc = tgt_snode.FirstCol();
//////          //find the pivot idx
//////          do {first_pivot_idx = src_snode.FindBlockIdx(tgt_fc); tgt_fc++;}
//////          while(first_pivot_idx<0 && tgt_fc<=tgt_snode.LastCol());
//////          tgt_fc--;
//////        }
//////        else{
//////          first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
//////        }
//////        assert(first_pivot_idx>=0);
//////        NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);
//////
//////        //find the last row updated by src_snode
//////        Int tgt_lc = tgt_snode.LastCol();
//////        Int last_pivot_idx = -1;
//////        //find the pivot idx
//////        do {last_pivot_idx = src_snode.FindBlockIdx(tgt_lc); tgt_lc--;}
//////        while(last_pivot_idx<0 && tgt_lc>=tgt_fc);
//////        tgt_lc++;
//////        assert(last_pivot_idx>=0);
//////        NZBlockDesc & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);
//////        TIMER_STOP(UPDATE_SNODE_FIND_INDEX);
//
//      Int tgt_fc,tgt_lc;
//      Int first_pivot_idx,last_pivot_idx;
//      tgt_snode.FindUpdatedFirstCol(src_snode, update.src_first_row, tgt_fc, first_pivot_idx);
//      tgt_snode.FindUpdatedLastCol(src_snode, tgt_fc, first_pivot_idx, tgt_lc, last_pivot_idx);
//
//      NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);
//      NZBlockDesc & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);
//
//
//
//
//
//        //determine the first column that will be updated in the target supernode
//        Int tgt_local_fc =  tgt_fc - tgt_snode.FirstCol();
//        Int tgt_local_lc =  tgt_lc - tgt_snode.FirstCol();
//
//        Int tgt_nrows = tgt_snode.NRowsBelowBlock(0);
//        Int src_nrows = src_snode.NRowsBelowBlock(first_pivot_idx)
//          - (tgt_fc - first_pivot_desc.GIndex);
//        Int src_lr = tgt_fc+src_nrows-1;
//        src_nrows = src_lr - tgt_fc + 1;
//
//        Int tgt_width = src_nrows - src_snode.NRowsBelowBlock(last_pivot_idx)
//          + (tgt_lc - last_pivot_desc.GIndex)+1;
//
//        SCALAR * pivot = src_snode.GetNZval(first_pivot_desc.Offset)
//          + (tgt_fc-first_pivot_desc.GIndex)*src_snode_size;
//        SCALAR * tgt = tgt_snode.GetNZval(0);
//
//        //Pointer to the output buffer of the GEMM
//        SCALAR * buf = NULL;
//        SCALAR beta = ZERO<SCALAR>();
//        //Compute the update in a temporary buffer
//#ifdef _DEBUG_
//        tmpBuffers.tmpBuf.Resize(tgt_width,src_nrows);
//#endif
//
//        buf = tmpBuffers.tmpBuf.Data();
//
//        //everything is in row-major
//        TIMER_START(UPDATE_SNODE_GEMM);
//        blas::Gemm('T','N',tgt_width, src_nrows,src_snode_size,
//            MINUS_ONE<SCALAR>(),pivot,src_snode_size,
//            pivot,src_snode_size,beta,buf,tgt_width);
//        TIMER_STOP(UPDATE_SNODE_GEMM);
//
//        //If the GEMM wasn't done in place we need to aggregate the update
//        //This is the assembly phase
//        if(1){
//#ifdef _DEBUG_
//          logfileptr->OFS()<<"tmpBuf is "<<tmpBuffers.tmpBuf<<std::endl;
//#endif
//
//          //now add the update to the target supernode
//          TIMER_START(UPDATE_SNODE_INDEX_MAP);
//          if(tgt_snode_size==1){
//            Int rowidx = 0;
//            Int src_blkcnt = src_snode.NZBlockCnt();
//            for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
//              NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
//              Int cur_src_nrows = src_snode.NRows(blkidx);
//              Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
//              Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
//
//              Int row = cur_src_fr;
//              while(row<=cur_src_lr){
//                Int tgt_blk_idx = tgt_snode.FindBlockIdx(row);
//                assert(tgt_blk_idx>=0);
//                NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
//                Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + tgt_snode.NRows(tgt_blk_idx)-1);
//                Int tgtOffset = cur_tgt_desc.Offset 
//                  + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
//                for(Int cr = row ;cr<=lr;++cr){
//                  tgt[tgtOffset + (cr - row)*tgt_snode_size] += buf[rowidx]; 
//                  rowidx++;
//                }
//                row += (lr-row+1);
//              }
//            }
//          }
//          else{
//            tmpBuffers.src_colindx.Resize(tgt_width);
//            tmpBuffers.src_to_tgt_offset.Resize(src_nrows);
//            Int colidx = 0;
//            Int rowidx = 0;
//            Int offset = 0;
//
//            Int src_blkcnt = src_snode.NZBlockCnt();
//            for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
//              NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
//              Int cur_src_nrows = src_snode.NRows(blkidx);
//              Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
//              Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
//              cur_src_nrows = cur_src_lr - cur_src_fr +1;
//
//              //The other one MUST reside into a single block in the target
//              Int row = cur_src_fr;
//              while(row<=cur_src_lr){
//                Int tgt_blk_idx = tgt_snode.FindBlockIdx(row);
//                assert(tgt_blk_idx>=0);
//                NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
//                Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + tgt_snode.NRows(tgt_blk_idx)-1);
//                Int tgtOffset = cur_tgt_desc.Offset 
//                  + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
//                for(Int cr = row ;cr<=lr;++cr){
//                  if(cr<=tgt_lc){
//                    tmpBuffers.src_colindx[colidx++] = cr;
//                  }
//                  offset+=tgt_width;
//                  tmpBuffers.src_to_tgt_offset[rowidx] = tgtOffset + (cr - row)*tgt_snode_size;
//                  rowidx++;
//                }
//                row += (lr-row+1);
//              }
//            }
//
//
//            //Multiple cases to consider
//            TIMER_STOP(UPDATE_SNODE_INDEX_MAP);
//
//            if(first_pivot_idx==last_pivot_idx){
//              // Updating contiguous columns
//              Int tgt_offset = (tgt_fc - tgt_snode.FirstCol());
//              for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
//                blas::Axpy(tgt_width,ONE<SCALAR>(),&buf[rowidx*tgt_width],1,
//                    &tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_offset],1);
//              }
//            }
//            else{
//              // full sparse case (done right now)
//              for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
//                for(Int colidx = 0; colidx< tmpBuffers.src_colindx.m();++colidx){
//                  Int col = tmpBuffers.src_colindx[colidx];
//                  Int tgt_colidx = col - tgt_snode.FirstCol();
//                  tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_colidx] 
//                    += buf[rowidx*tgt_width+colidx]; 
//                }
//              }
//            }
//          }
//        }
//        TIMER_STOP(UPDATE_SNODE_OBJ);
//        return 0;
//
//      }
//      else{
//        trc_Update(src_snode, update, tmpBuffers,tgt_snode);
//      }
//
//    }
//
//
//    Int trc_Update(SuperNode<SCALAR> & src_snode, SnodeUpdate &update, TempUpdateBuffers<SCALAR> & tmpBuffers,SuperNode<SCALAR> & tgt_snode){
//      Int & pivot_idx = update.blkidx;
//      Int & pivot_fr = update.src_first_row;
//
//      TIMER_START(UPDATE_SNODE_OBJ);
//      Int src_snode_size = src_snode.Size();
//      Int tgt_snode_size = tgt_snode.Size();
//
//      //find the first row updated by src_snode
//      Int tgt_fc,tgt_lc;
//      Int first_pivot_idx,last_pivot_idx;
//      tgt_snode.FindUpdatedFirstCol(src_snode, update.src_first_row, tgt_fc, first_pivot_idx);
//      tgt_snode.FindUpdatedLastCol(src_snode, tgt_fc, first_pivot_idx, tgt_lc, last_pivot_idx);
//
//      NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);
//      NZBlockDesc & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);
//
//      //determine the first column that will be updated in the target supernode
//      Int tgt_local_fc =  tgt_fc - tgt_snode.FirstCol();
//      Int tgt_local_lc =  tgt_lc - tgt_snode.FirstCol();
//
//      Int tgt_nrows = tgt_snode.NRowsBelowBlock(0);
//      Int src_nrows = src_snode.NRowsBelowBlock(first_pivot_idx)
//        - (tgt_fc - first_pivot_desc.GIndex);
//      Int src_lr = tgt_fc+src_nrows-1;
//      src_nrows = src_lr - tgt_fc + 1;
//
//      Int tgt_width = src_nrows - src_snode.NRowsBelowBlock(last_pivot_idx)
//        + (tgt_lc - last_pivot_desc.GIndex)+1;
//
//      SCALAR * pivot = src_snode.GetNZval(first_pivot_desc.Offset)
//        + (tgt_fc-first_pivot_desc.GIndex)*src_snode_size;
//      SCALAR * tgt = tgt_snode.GetNZval(0);
//
//      //Pointer to the output buffer of the GEMM
//      SCALAR * buf = NULL;
//      SCALAR beta = ZERO<SCALAR>();
//      //If the target supernode has the same structure,
//      //The GEMM is directly done in place
//      if(src_nrows == tgt_nrows){
//        Int tgt_offset = (tgt_fc - tgt_snode.FirstCol());
//        buf = &tgt[tgt_offset];
//        beta = ONE<SCALAR>();
//      }
//      else{
//        //Compute the update in a temporary buffer
//#ifdef _DEBUG_
//        tmpBuffers.tmpBuf.Resize(tgt_width,src_nrows);
//#endif
//        buf = tmpBuffers.tmpBuf.Data();
//      }
//
//      //everything is in row-major
//      TIMER_START(UPDATE_SNODE_GEMM);
//      blas::Gemm('T','N',tgt_width, src_nrows,src_snode_size,
//          MINUS_ONE<SCALAR>(),pivot,src_snode_size,
//          pivot,src_snode_size,beta,buf,tgt_width);
//      TIMER_STOP(UPDATE_SNODE_GEMM);
//
//      //If the GEMM wasn't done in place we need to aggregate the update
//      //This is the assembly phase
//      if(src_nrows != tgt_nrows){
//#ifdef _DEBUG_
//        logfileptr->OFS()<<"tmpBuf is "<<tmpBuffers.tmpBuf<<std::endl;
//#endif
//
//        //now add the update to the target supernode
//        TIMER_START(UPDATE_SNODE_INDEX_MAP);
//        if(tgt_snode_size==1){
//          Int rowidx = 0;
//          Int src_blkcnt = src_snode.NZBlockCnt();
//          for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
//            NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
//            Int cur_src_nrows = src_snode.NRows(blkidx);
//            Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
//            Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
//
//            Int row = cur_src_fr;
//            while(row<=cur_src_lr){
//              Int tgt_blk_idx = tgt_snode.FindBlockIdx(row);
//              assert(tgt_blk_idx>=0);
//              NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
//              Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + tgt_snode.NRows(tgt_blk_idx)-1);
//              Int tgtOffset = cur_tgt_desc.Offset 
//                + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
//              for(Int cr = row ;cr<=lr;++cr){
//                tgt[tgtOffset + (cr - row)*tgt_snode_size] += buf[rowidx]; 
//                rowidx++;
//              }
//              row += (lr-row+1);
//            }
//          }
//        }
//        else{
//          tmpBuffers.src_colindx.Resize(tgt_width);
//          tmpBuffers.src_to_tgt_offset.Resize(src_nrows);
//          Int colidx = 0;
//          Int rowidx = 0;
//          Int offset = 0;
//
//          Int src_blkcnt = src_snode.NZBlockCnt();
//          for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
//            NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
//            Int cur_src_nrows = src_snode.NRows(blkidx);
//            Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
//            Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
//            cur_src_nrows = cur_src_lr - cur_src_fr +1;
//
//            //The other one MUST reside into a single block in the target
//            Int row = cur_src_fr;
//            while(row<=cur_src_lr){
//              Int tgt_blk_idx = tgt_snode.FindBlockIdx(row);
//              assert(tgt_blk_idx>=0);
//              NZBlockDesc & cur_tgt_desc = tgt_snode.GetNZBlockDesc(tgt_blk_idx);
//              Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + tgt_snode.NRows(tgt_blk_idx)-1);
//              Int tgtOffset = cur_tgt_desc.Offset 
//                + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
//              for(Int cr = row ;cr<=lr;++cr){
//                if(cr<=tgt_lc){
//                  tmpBuffers.src_colindx[colidx++] = cr;
//                }
//                offset+=tgt_width;
//                tmpBuffers.src_to_tgt_offset[rowidx] = tgtOffset + (cr - row)*tgt_snode_size;
//                rowidx++;
//              }
//              row += (lr-row+1);
//            }
//          }
//          //Multiple cases to consider
//          TIMER_STOP(UPDATE_SNODE_INDEX_MAP);
//
//          if(first_pivot_idx==last_pivot_idx){
//            // Updating contiguous columns
//            Int tgt_offset = (tgt_fc - tgt_snode.FirstCol());
//            for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
//              blas::Axpy(tgt_width,ONE<SCALAR>(),&buf[rowidx*tgt_width],1,
//                  &tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_offset],1);
//            }
//          }
//          else{
//            // full sparse case (done right now)
//            for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
//              for(Int colidx = 0; colidx< tmpBuffers.src_colindx.m();++colidx){
//                Int col = tmpBuffers.src_colindx[colidx];
//                Int tgt_colidx = col - tgt_snode.FirstCol();
//                tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_colidx] 
//                  += buf[rowidx*tgt_width+colidx]; 
//              }
//            }
//          }
//        }
//      }
//      TIMER_STOP(UPDATE_SNODE_OBJ);
//      return 0;
//    }
//
//    Int trc_Factorize(SuperNode<SCALAR> & tgt_snode){
//      Int tgt_snode_size = tgt_snode.Size();
//      Int BLOCKSIZE = tgt_snode_size;
//      NZBlockDesc & diag_desc = tgt_snode.GetNZBlockDesc(0);
//      for(Int col = 0; col<tgt_snode_size;col+=BLOCKSIZE){
//        Int bw = min(BLOCKSIZE,tgt_snode_size-col);
//        SCALAR * diag_nzval = &tgt_snode.GetNZval(diag_desc.Offset)[col+col*tgt_snode_size];
//        lapack::Potrf( 'U', bw, diag_nzval, tgt_snode_size);
//        SCALAR * nzblk_nzval = &tgt_snode.GetNZval(diag_desc.Offset)[col+(col+bw)*tgt_snode_size];
//        blas::Trsm('L','U','T','N',bw, tgt_snode.NRowsBelowBlock(0)-(col+bw), ONE<SCALAR>(),  diag_nzval, tgt_snode_size, nzblk_nzval, tgt_snode_size);
//
//        //update the rest !!! (next blocks columns)
//        SCALAR * tgt_nzval = &tgt_snode.GetNZval(diag_desc.Offset)[col+bw+(col+bw)*tgt_snode_size];
//        blas::Gemm('T','N',tgt_snode_size-(col+bw), tgt_snode.NRowsBelowBlock(0)-(col+bw),bw,MINUS_ONE<SCALAR>(),nzblk_nzval,tgt_snode_size,nzblk_nzval,tgt_snode_size,ONE<SCALAR>(),tgt_nzval,tgt_snode_size);
//      }
//      return 0;
//    }
//
//
//    void trc_SNodeInit(Int aiId, Int aiFc, Int aiLc, NZBlockDesc * a_block_desc, Int a_desc_cnt, SCALAR * a_nzval, Int a_nzval_cnt, SuperNode<SCALAR> & tgt_snode, Int aiN) {
//      tgt_snode.iId_ = aiId;
//      tgt_snode.iN_ = aiN;
//      tgt_snode.iFirstCol_=aiFc;
//      tgt_snode.iLastCol_ =aiLc;
//      //compute supernode size / width
//      tgt_snode.iSize_ = tgt_snode.iLastCol_ - tgt_snode.iFirstCol_+1;
//
//      tgt_snode.b_own_storage_ = false;
//      tgt_snode.nzval_= a_nzval;
//      tgt_snode.nzval_cnt_ = a_nzval_cnt;
//
//      assert(tgt_snode.nzval_cnt_ % tgt_snode.iSize_ == 0);
//
//
//      tgt_snode.blocks_ = a_block_desc;
//      tgt_snode.blocks_cnt_ = a_desc_cnt;
//
//      //restore 0-based offsets and compute global_to_local indices
//      for(Int blkidx=tgt_snode.blocks_cnt_-1; blkidx>=0;--blkidx){
//        tgt_snode.blocks_[blkidx].Offset -= tgt_snode.blocks_[0].Offset;
//      }
//
//
//
//#ifndef ITREE
//      //tgt_snode.iLastRow_ = max(tgt_snode.iLastCol_,tgt_snode.iFirstCol_ + tgt_snode.nzval_cnt_/tgt_snode.iSize_ - 1);
//      //tgt_snode.globalToLocal_ = new SYMPACK::vector<Int>(tgt_snode.iLastRow_ - tgt_snode.iFirstCol_ + 1 );
//      tgt_snode.globalToLocal_ = new SYMPACK::vector<Int>(tgt_snode.iN_+1,-1);
//#else
//      tgt_snode.idxToBlk_ = tgt_snode.CreateITree();
//#endif
//
//
//      //  tgt_snode.initIdxToBlk_(true);
//
//
//    }
//}
