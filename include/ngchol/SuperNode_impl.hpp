#ifndef _SUPERNODE_IMPL_HPP_
#define _SUPERNODE_IMPL_HPP_



//namespace LIBCHOLESKY{




//SuperNode implementation
  template<typename T>
    SuperNode<T>::SuperNode() :iId_(-1),iFirstCol_(-1),iLastCol_(-1) { }

  template<typename T>
    SuperNode<T>::SuperNode(Int aiId, Int aiFc, Int aiLc, Int ai_num_rows) :iId_(aiId),iFirstCol_(aiFc),iLastCol_(aiLc) {
      //this is an upper bound
      assert(ai_num_rows>=0);

      //compute supernode size / width
      iSize_ = iLastCol_ - iFirstCol_+1;

      nzval_cnt_ = 0;
      nzval_container_.reserve(iSize_*ai_num_rows);
      nzval_ = NULL;

      blocks_container_.reserve(ai_num_rows);
      blocks_ = NULL;
      blocks_cnt_ = 0;

      idxToBlk_ = CreateITree();

      b_own_storage_ = true;


    }; 

  template<typename T>
    SuperNode<T>::SuperNode(Int aiId, Int aiFc, Int aiLc) :iId_(aiId),iFirstCol_(aiFc),iLastCol_(aiLc) {

      //compute supernode size / width
      iSize_ = iLastCol_ - iFirstCol_+1;

      nzval_cnt_ = 0;
      nzval_container_.reserve(iSize_*iSize_);
      nzval_ = NULL;

      blocks_container_.reserve(1);
      blocks_ = NULL;
      blocks_cnt_ = 0;

      idxToBlk_ = CreateITree();

      b_own_storage_ = true;

      //add the diagonal block
      AddNZBlock(iSize_,0,iFirstCol_);

    }; 

  template<typename T>
    SuperNode<T>::SuperNode(Int aiId, Int aiFc, Int aiLc, NZBlockDesc * a_block_desc, Int a_desc_cnt,
        T * a_nzval, Int a_nzval_cnt) {
      Init(aiId, aiFc, aiLc, a_block_desc, a_desc_cnt, a_nzval, a_nzval_cnt);
    }

  template<typename T>
    SuperNode<T>::SuperNode(Int aiId, Int aiFc, Int aiLc, IntNumVec & xlindx, IntNumVec & lindx) :iId_(aiId),iFirstCol_(aiFc),iLastCol_(aiLc) {

      //compute supernode size / width
      iSize_ = iLastCol_ - iFirstCol_+1;
      b_own_storage_ = true;

      std::list<NZBlockDesc> tmpBlockIndex;

      Int fi = xlindx(iId_-1);
      Int li = xlindx(iId_)-1;



      Int prevRow = 0;
      nzval_cnt_ = 0;
      blocks_cnt_ = 0;


      for(Int idx= fi ; idx<=li; ++idx){
        Int curRow = lindx[idx-1];
        //create a new block
        if(curRow==iFirstCol_ || curRow != prevRow+1 || (curRow>iLastCol_ && blocks_cnt_==1) ){

          tmpBlockIndex.push_back(NZBlockDesc(curRow,nzval_cnt_));
          ++blocks_cnt_;
        }

        prevRow = curRow;

        nzval_cnt_+=iSize_;
      }

      blocks_container_.resize(tmpBlockIndex.size());
      std::copy(tmpBlockIndex.begin(),tmpBlockIndex.end(),blocks_container_.begin());





      nzval_container_.resize(nzval_cnt_,ZERO<T>());

      nzval_ = &nzval_container_.front();
      blocks_ = &blocks_container_.front();


      idxToBlk_ = CreateITree();

    }

  template<typename T>
    SuperNode<T>::~SuperNode(){
      delete idxToBlk_;
    }

  template<typename T>
    void SuperNode<T>::Init(Int aiId, Int aiFc, Int aiLc, NZBlockDesc * a_block_desc, Int a_desc_cnt,
        T * a_nzval, Int a_nzval_cnt) {
#ifndef _TAU_TRACE_
      iId_ = aiId;
      iFirstCol_=aiFc;
      iLastCol_ =aiLc;
      //compute supernode size / width
      iSize_ = iLastCol_ - iFirstCol_+1;

      b_own_storage_ = false;
      nzval_= a_nzval;
      nzval_cnt_ = a_nzval_cnt;

      assert(nzval_cnt_ % iSize_ == 0);


      blocks_ = a_block_desc;
      blocks_cnt_ = a_desc_cnt;
      //restore 0-based offsets and compute global_to_local indices
      for(Int blkidx=blocks_cnt_-1; blkidx>=0;--blkidx){
        blocks_[blkidx].Offset -= blocks_[0].Offset;
      }

      //initIdxToBlk_(true);
      idxToBlk_ = CreateITree();
#else
      trc_SNodeInit(aiId, aiFc, aiLc, a_block_desc, a_desc_cnt, a_nzval, a_nzval_cnt, *this);
#endif
    }

  template<typename T>
    inline void SuperNode<T>::AddNZBlock(Int aiNRows, Int aiNCols, Int aiGIndex){

      //Resize the container if I own the storage
      if(b_own_storage_){

        Int cur_fr = aiGIndex;
        Int cur_lr = cur_fr + aiNRows -1;

//#ifdef INTERVAL_TREE
        
        ITree::Interval cur_interv = { cur_fr, cur_lr, blocks_cnt_};
        idxToBlk_->Insert(cur_interv);
//#endif

        blocks_container_.push_back(NZBlockDesc(aiGIndex,nzval_cnt_));

        blocks_cnt_++;
        nzval_cnt_+=aiNRows*iSize_;
        nzval_container_.resize(nzval_cnt_,ZERO<T>());
        nzval_ = &nzval_container_.front();
        blocks_ = &blocks_container_.front();
      }
    }


  template<typename T>
    inline Int SuperNode<T>::FindBlockIdx(Int aiGIndex){
      TIMER_START(FindBlockIdx);
//      if(Id()==3 && aiGIndex==46){gdb_lock(0);}
      //      ITree::Interval it = {aiGIndex, aiGIndex,0};
      //      ITree::Interval * res = idxToBlk_->IntervalSearch(it);
//#ifdef INTERVAL_TREE
#ifdef _LAZY_INIT_
    if(!ITreeInitialized()){
      InitIdxToBlk();
    }
#endif


      ITree::Interval * res = idxToBlk_->IntervalSearch(aiGIndex,aiGIndex);
      if (res == NULL){
        TIMER_STOP(FindBlockIdx);
        return -1;
      }
      else{
        TIMER_STOP(FindBlockIdx);
        return res->block_idx;
      }
//#endif
    }



  template<typename T>
  inline Int SuperNode<T>::FindBlockIdx(Int fr, Int lr, ITree::Interval & overlap){
      TIMER_START(FindBlockIdx);
      //      ITree::Interval it = {aiGIndex, aiGIndex,0};
      //      ITree::Interval * res = idxToBlk_->IntervalSearch(it);
#ifdef _LAZY_INIT_
    if(!ITreeInitialized()){
      InitIdxToBlk();
    }
#endif


      ITree::Interval * res = idxToBlk_->IntervalSearch(fr,lr);
      if (res == NULL){
        TIMER_STOP(FindBlockIdx);
        return -1;
      }
      else{
        TIMER_STOP(FindBlockIdx);
        overlap = *res;
        return res->block_idx;
      }
    }





  template<typename T>
    inline void SuperNode<T>::DumpITree(){
      logfileptr->OFS()<<"Number of blocks: "<<blocks_cnt_<<endl;
      logfileptr->OFS()<<"log2(Number of blocks): "<<log2(blocks_cnt_)<<endl;
      idxToBlk_->Dump();
    }

  template<typename T>
    inline Int SuperNode<T>::Shrink(){
      if(b_own_storage_){
        blocks_ = &blocks_container_.front();

        nzval_container_.resize(nzval_cnt_);
        nzval_ = &nzval_container_.front();
      }

      return StorageSize();
    }







#ifdef COMPACT_AGGREGATES
  template<typename T>
    inline Int SuperNode<T>::Merge(SuperNode<T> & src_snode, SnodeUpdate &update){
#ifndef _TAU_TRACE_
      TIMER_START(MERGE_SNODE);

      assert(b_own_storage_);

      Int src_snode_size = src_snode.Size();
      Int tgt_snode_size = Size();

      Int & pivot_idx = update.blkidx;
      Int & pivot_fr = update.src_first_row;

////      //find the first row updated by src_snode
////      TIMER_START(MERGE_SNODE_FIND_INDEX);
////      Int first_pivot_idx = -1;
////      Int tgt_fc = pivot_fr;
////      if(tgt_fc ==I_ZERO ){
////        tgt_fc = FirstCol();
////        //find the pivot idx
////        do {first_pivot_idx = src_snode.FindBlockIdx(tgt_fc); tgt_fc++;}
////        while(first_pivot_idx<0 && tgt_fc<=LastCol());
////        tgt_fc--;
////      }
////      else{
////        first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
////      }
////      assert(first_pivot_idx>=0);
////
////      TIMER_STOP(MERGE_SNODE_FIND_INDEX);


      Int tgt_fc;
      Int first_pivot_idx;
      FindUpdatedFirstCol(src_snode, 0, tgt_fc, first_pivot_idx);
      NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);

      //parse src_snode
      ITree::Interval overlap;
      for(Int blkidx = first_pivot_idx; blkidx<src_snode.NZBlockCnt(); ++blkidx){
        NZBlockDesc & blk_desc = src_snode.GetNZBlockDesc(blkidx);
        Int fr = blk_desc.GIndex;
        Int lr = blk_desc.GIndex + src_snode.NRows(blkidx) -1;
        
        if(FindBlockIdx(fr,lr,overlap)==-1){
          //Add the full block
          AddNZBlock( lr - fr + 1, tgt_snode_size, fr);
        }
        else{
          
          //check the overlap
          //                l-----overlap------h
          //            l---------block-------------h
          //        l--block----h
          //                              l-----block----h
          
          if(overlap.high < lr){
            //we need to add from high+1 to lr 
            AddNZBlock( lr - overlap.high, tgt_snode_size, overlap.high+1);
          }
          
          if(overlap.low>fr){
            //we need to add fr to low-1 
            AddNZBlock( overlap.low - fr, tgt_snode_size, fr);
          }
        }


//        Int nrows = src_snode.NRows(blkidx);
//        for(Int rowidx = 0; rowidx<nrows; ++rowidx){
//          Int row = blk_desc.GIndex + rowidx;
//          //if the row is updating the target
//          if(row>=tgt_fc){
//            //check if the row is not already in the structure
//            if(FindBlockIdx(row)==-1){
//              //add a nzblock with a single row in it
//              AddNZBlock( 1, tgt_snode_size, row);
//            }
//          }
//        }
      }

      TIMER_STOP(MERGE_SNODE);
      return 0;
#else
    return trc_Merge(src_snode, update,*this);
#endif
    }
#endif

#ifndef _TAU_TRACE_
#ifdef COMPACT_AGGREGATES
  template<typename T>
    inline Int SuperNode<T>::Aggregate(SuperNode<T> & src_snode){
      TIMER_START(AGGREG_SNODE);
      Int  pivot_idx = 0;
      Int  pivot_fr = 0;//FirstCol();

      Int src_snode_size = src_snode.Size();
      Int tgt_snode_size = Size();

////      TIMER_START(AGGREG_SNODE_FIND_INDEX);
////      Int first_pivot_idx = -1;
////      Int tgt_fc = pivot_fr;
////      if(tgt_fc ==I_ZERO ){
////        tgt_fc = FirstCol();
////        //find the pivot idx
////        do {first_pivot_idx = src_snode.FindBlockIdx(tgt_fc); tgt_fc++;}
////        while(first_pivot_idx<0 && tgt_fc<=LastCol());
////        tgt_fc--;
////      }
////      else{
////        first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
////      }
////      assert(first_pivot_idx>=0);
////      NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);
////
////      TIMER_STOP(AGGREG_SNODE_FIND_INDEX);

//      Int tgt_fc;
//      Int first_pivot_idx;
//      FindUpdatedFirstCol(src_snode, 0, tgt_fc, first_pivot_idx);
//      NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);



      //parse src_snode and add everything

      Int first_pivot_idx = 0 ;
      Int tgt_fc = FirstCol();

      for(Int blkidx = first_pivot_idx; blkidx<src_snode.NZBlockCnt(); ++blkidx){
        NZBlockDesc & blk_desc = src_snode.GetNZBlockDesc(blkidx);
        Int nrows = src_snode.NRows(blkidx);
        for(Int rowidx = 0; rowidx<nrows; ++rowidx){
          Int row = blk_desc.GIndex + rowidx;

          if(row>=tgt_fc){
            Int src_offset = blk_desc.Offset + (row - blk_desc.GIndex)*src_snode_size;

            Int tgt_blkidx = FindBlockIdx(row);
            assert(tgt_blkidx!=-1);
            NZBlockDesc & tgt_desc = GetNZBlockDesc(tgt_blkidx);
            Int tgt_offset = tgt_desc.Offset + (row - tgt_desc.GIndex)*tgt_snode_size;

            T * src = src_snode.GetNZval(src_offset);
            T * tgt = GetNZval(tgt_offset);

            //blas::Axpy(tgt_snode_size,ONE<T>(),src,1,tgt,1);
            for(Int i = 0; i< tgt_snode_size;i+=1){
              tgt[i] += src[i];
            }

          }
        }
      }


      TIMER_STOP(AGGREG_SNODE);
      return 0;
    }
#else
  template<typename T>
    inline Int SuperNode<T>::Aggregate(SuperNode<T> & src_snode){

      TIMER_START(AGGREGATE_SNODE);

      Int src_snode_size = src_snode.Size();
      Int tgt_snode_size = Size();


      T * src = src_snode.GetNZval(0);
      T * tgt = GetNZval(0);

      //blas::Axpy(tgt_snode_size*NRowsBelowBlock(0),ONE<T>(),src,1,tgt,1);
      for(Int i = 0; i< tgt_snode_size*NRowsBelowBlock(0);i+=1){
        tgt[i] += src[i];
      }



      TIMER_STOP(AGGREGATE_SNODE);

      return 0;
    }
#endif
#else
  template<typename T>
    inline Int SuperNode<T>::Aggregate(SuperNode<T> & src_snode){
    return trc_Aggregate(src_snode,*this);
  }
#endif


#ifndef _TAU_TRACE_
#ifdef COMPACT_AGGREGATES
  template<typename T>
    inline Int SuperNode<T>::UpdateAggregate(SuperNode<T> & src_snode, SnodeUpdate &update, 
        TempUpdateBuffers<T> & tmpBuffers, Int iTarget){

      if(iTarget != iam){

        Merge(src_snode, update);

        Int & pivot_idx = update.blkidx;
        Int & pivot_fr = update.src_first_row;

        TIMER_START(UPDATE_SNODE);
        Int src_snode_size = src_snode.Size();
        Int tgt_snode_size = Size();

////        //find the first row updated by src_snode
////        TIMER_START(UPDATE_SNODE_FIND_INDEX);
////        Int first_pivot_idx = -1;
////        Int tgt_fc = pivot_fr;
////        if(tgt_fc ==I_ZERO ){
////          tgt_fc = FirstCol();
////          //find the pivot idx
////          do {first_pivot_idx = src_snode.FindBlockIdx(tgt_fc); tgt_fc++;}
////          while(first_pivot_idx<0 && tgt_fc<=LastCol());
////          tgt_fc--;
////        }
////        else{
////          first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
////        }
////        assert(first_pivot_idx>=0);
////        NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);
////
////        //find the last row updated by src_snode
////        Int tgt_lc = LastCol();
////        Int last_pivot_idx = -1;
////        //find the pivot idx
////        do {last_pivot_idx = src_snode.FindBlockIdx(tgt_lc); tgt_lc--;}
////        while(last_pivot_idx<0 && tgt_lc>=tgt_fc);
////        tgt_lc++;
////        assert(last_pivot_idx>=0);
////        NZBlockDesc & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);
////        TIMER_STOP(UPDATE_SNODE_FIND_INDEX);

      Int tgt_fc,tgt_lc;
      Int first_pivot_idx,last_pivot_idx;
      FindUpdatedFirstCol(src_snode, update.src_first_row, tgt_fc, first_pivot_idx);
      FindUpdatedLastCol(src_snode, tgt_fc, first_pivot_idx, tgt_lc, last_pivot_idx);

      NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);
      NZBlockDesc & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);





        //determine the first column that will be updated in the target supernode
        Int tgt_local_fc =  tgt_fc - FirstCol();
        Int tgt_local_lc =  tgt_lc - FirstCol();

        Int tgt_nrows = NRowsBelowBlock(0);
        Int src_nrows = src_snode.NRowsBelowBlock(first_pivot_idx)
          - (tgt_fc - first_pivot_desc.GIndex);
        Int src_lr = tgt_fc+src_nrows-1;
        src_nrows = src_lr - tgt_fc + 1;

        Int tgt_width = src_nrows - src_snode.NRowsBelowBlock(last_pivot_idx)
          + (tgt_lc - last_pivot_desc.GIndex)+1;

        T * pivot = src_snode.GetNZval(first_pivot_desc.Offset)
          + (tgt_fc-first_pivot_desc.GIndex)*src_snode_size;
        T * tgt = GetNZval(0);

        //Pointer to the output buffer of the GEMM
        T * buf = NULL;
        T beta = ZERO<T>();
        //If the target supernode has the same structure,
        //The GEMM is directly done in place
        //    if(src_nrows == tgt_nrows){
        //      Int tgt_offset = (tgt_fc - FirstCol());
        //      buf = &tgt[tgt_offset];
        //      beta = ONE<T>();
        //    }
        //    else{
        //Compute the update in a temporary buffer
#ifdef _DEBUG_
        tmpBuffers.tmpBuf.Resize(tgt_width,src_nrows);
#endif

        buf = tmpBuffers.tmpBuf.Data();
        //    }

        //everything is in row-major
        TIMER_START(UPDATE_SNODE_GEMM);
        blas::Gemm('T','N',tgt_width, src_nrows,src_snode_size,
            MINUS_ONE<T>(),pivot,src_snode_size,
            pivot,src_snode_size,beta,buf,tgt_width);
        TIMER_STOP(UPDATE_SNODE_GEMM);

        //If the GEMM wasn't done in place we need to aggregate the update
        //This is the assembly phase
        if(1){
#ifdef _DEBUG_
          logfileptr->OFS()<<"tmpBuf is "<<tmpBuffers.tmpBuf<<std::endl;
#endif

          //now add the update to the target supernode
          TIMER_START(UPDATE_SNODE_INDEX_MAP);
          if(tgt_snode_size==1){
            Int rowidx = 0;
            Int src_blkcnt = src_snode.NZBlockCnt();
            for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
              NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
              Int cur_src_nrows = src_snode.NRows(blkidx);
              Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
              Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);

              Int row = cur_src_fr;
              while(row<=cur_src_lr){
                Int tgt_blk_idx = FindBlockIdx(row);
                assert(tgt_blk_idx>=0);
                NZBlockDesc & cur_tgt_desc = GetNZBlockDesc(tgt_blk_idx);
                Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + NRows(tgt_blk_idx)-1);
                Int tgtOffset = cur_tgt_desc.Offset 
                  + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
                for(Int cr = row ;cr<=lr;++cr){
                  tgt[tgtOffset + (cr - row)*tgt_snode_size] += buf[rowidx]; 
                  rowidx++;
                }
                row += (lr-row+1);
              }
            }
          }
          else{
            tmpBuffers.src_colindx.Resize(tgt_width);
            tmpBuffers.src_to_tgt_offset.Resize(src_nrows);

            Int colidx = 0;
            Int rowidx = 0;
            Int offset = 0;


            Int src_blkcnt = src_snode.NZBlockCnt();
            for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
              NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
              Int cur_src_nrows = src_snode.NRows(blkidx);
              Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
              Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
              cur_src_nrows = cur_src_lr - cur_src_fr +1;

              //The other one MUST reside into a single block in the target
              //NOT TRUE
              Int row = cur_src_fr;
              while(row<=cur_src_lr){
                Int tgt_blk_idx = FindBlockIdx(row);
                assert(tgt_blk_idx>=0);
                NZBlockDesc & cur_tgt_desc = GetNZBlockDesc(tgt_blk_idx);
                Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + NRows(tgt_blk_idx)-1);
                Int tgtOffset = cur_tgt_desc.Offset 
                  + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
                for(Int cr = row ;cr<=lr;++cr){
                  if(cr<=tgt_lc){
                    tmpBuffers.src_colindx[colidx++] = cr;
                  }
                  offset+=tgt_width;
                  tmpBuffers.src_to_tgt_offset[rowidx] = tgtOffset + (cr - row)*tgt_snode_size;
                  rowidx++;
                }
                row += (lr-row+1);
              }
            }


            //Multiple cases to consider
            TIMER_STOP(UPDATE_SNODE_INDEX_MAP);

            if(first_pivot_idx==last_pivot_idx){
              // Updating contiguous columns
              Int tgt_offset = (tgt_fc - FirstCol());
              for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
                blas::Axpy(tgt_width,ONE<T>(),&buf[rowidx*tgt_width],1,
                    &tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_offset],1);
              }
            }
            else{
              // full sparse case (done right now)
              for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
                for(Int colidx = 0; colidx< tmpBuffers.src_colindx.m();++colidx){
                  Int col = tmpBuffers.src_colindx[colidx];
                  Int tgt_colidx = col - FirstCol();
                  tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_colidx] 
                    += buf[rowidx*tgt_width+colidx]; 
                }
              }
            }
          }
        }
        TIMER_STOP(UPDATE_SNODE);
        return 0;

      }
      else{
        Update(src_snode, update, tmpBuffers);
      }

    }
#else
  template<typename T>
    inline Int SuperNode<T>::UpdateAggregate(SuperNode<T> & src_snode, SnodeUpdate &update, 
        TempUpdateBuffers<T> & tmpBuffers, Int iTarget){
    Update(src_snode, update, tmpBuffers);
  }
#endif
#else
  template<typename T>
    inline Int SuperNode<T>::UpdateAggregate(SuperNode<T> & src_snode, SnodeUpdate &update, 
        TempUpdateBuffers<T> & tmpBuffers, Int iTarget){
    return trc_UpdateAggregate(src_snode, update, tmpBuffers, iTarget,*this);
  }
#endif

  template<typename T>
    inline void SuperNode<T>::FindUpdatedFirstCol(SuperNode<T> & src_snode, Int pivot_fr, Int & tgt_fc, Int & first_pivot_idx){
//gdb_lock();
      //find the first row updated by src_snode
      TIMER_START(UPDATE_SNODE_FIND_INDEX);


#ifdef _LINEAR_SEARCH_FCLC_
    if(src_snode.ITreeInitialized()){
#endif 
      tgt_fc = pivot_fr;
      first_pivot_idx = -1;
      if(tgt_fc == I_ZERO ){
        tgt_fc = FirstCol();
        //find the pivot idx
        do {first_pivot_idx = src_snode.FindBlockIdx(tgt_fc); tgt_fc++;}
        while(first_pivot_idx<0 && tgt_fc<=LastCol());
        tgt_fc--;
      }
      else{
        first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
      }
#ifdef _LINEAR_SEARCH_FCLC_
    }
    else{
      tgt_fc = pivot_fr;
      first_pivot_idx = -1;
      if(tgt_fc == I_ZERO ){
        tgt_fc = FirstCol();

        //this assumes that blocks are sorted...
        for(Int blkidx = 0; blkidx<src_snode.NZBlockCnt();++blkidx){
          NZBlockDesc & cur_block = src_snode.GetNZBlockDesc(blkidx);
          Int cur_fr = cur_block.GIndex;
          Int cur_lr = cur_block.GIndex + src_snode.NRows(blkidx)-1;

          if(tgt_fc < cur_fr){
            tgt_fc = cur_fr;
            first_pivot_idx = blkidx;
            break;
          }
          else if(tgt_fc<=cur_lr){
            first_pivot_idx = blkidx;
            break;
          }
        }

      }
      else{
        for(Int blkidx = 0; blkidx<src_snode.NZBlockCnt();++blkidx){
          NZBlockDesc & cur_block = src_snode.GetNZBlockDesc(blkidx);
          Int cur_fr = cur_block.GIndex;
          Int cur_lr = cur_block.GIndex + src_snode.NRows(blkidx)-1;
          if(cur_fr<= tgt_fc && tgt_fc<=cur_lr){
            first_pivot_idx = blkidx;
            break;
          }
        }
      }
    }
#endif
      if(first_pivot_idx<0){gdb_lock();}
      assert(first_pivot_idx>=0);
      TIMER_STOP(UPDATE_SNODE_FIND_INDEX);
  }

  template<typename T>
    inline void SuperNode<T>::FindUpdatedLastCol(SuperNode<T> & src_snode, Int tgt_fc, Int first_pivot_idx , Int & tgt_lc,  Int & last_pivot_idx){
//gdb_lock();
      TIMER_START(UPDATE_SNODE_FIND_INDEX);
#ifdef _LINEAR_SEARCH_FCLC_
    if(src_snode.ITreeInitialized()){
#endif 
      //find the last row updated by src_snode
      tgt_lc = LastCol();
      last_pivot_idx = -1;
      //find the pivot idx
      do {last_pivot_idx = src_snode.FindBlockIdx(tgt_lc); tgt_lc--;}
      while(last_pivot_idx<0 && tgt_lc>=tgt_fc);
      tgt_lc++;
#ifdef _LINEAR_SEARCH_FCLC_
    }
    else{
        tgt_lc = tgt_fc;
        for(Int blkidx = first_pivot_idx; blkidx<src_snode.NZBlockCnt();++blkidx){
          NZBlockDesc & cur_block = src_snode.GetNZBlockDesc(blkidx);
          Int cur_fr = cur_block.GIndex;
          Int cur_lr = cur_block.GIndex + src_snode.NRows(blkidx)-1;
          
          
          if(cur_fr <= LastCol()){
             if(LastCol()<cur_lr){
              last_pivot_idx = blkidx;
              tgt_lc = LastCol();
             }
             else if(LastCol() == cur_lr){
              last_pivot_idx = blkidx;
              tgt_lc = LastCol();
              break;
             }
             else{
              last_pivot_idx = blkidx;
              tgt_lc = cur_lr;
            }
          }
          else{
            break;
          }
        }
    }
#endif
      assert(last_pivot_idx>=0);
      TIMER_STOP(UPDATE_SNODE_FIND_INDEX);
  }

#ifndef _TAU_TRACE_
  template<typename T>
    inline Int SuperNode<T>::Update(SuperNode<T> & src_snode, SnodeUpdate &update, 
        TempUpdateBuffers<T> & tmpBuffers){
      Int & pivot_idx = update.blkidx;
      Int & pivot_fr = update.src_first_row;

      TIMER_START(UPDATE_SNODE);
      Int src_snode_size = src_snode.Size();
      Int tgt_snode_size = Size();

      //find the first row updated by src_snode
////      TIMER_START(UPDATE_SNODE_FIND_INDEX);
////      Int first_pivot_idx = -1;
////      Int tgt_fc = pivot_fr;
////      if(tgt_fc ==I_ZERO ){
////        tgt_fc = FirstCol();
////        //find the pivot idx
////        do {first_pivot_idx = src_snode.FindBlockIdx(tgt_fc); tgt_fc++;}
////        while(first_pivot_idx<0 && tgt_fc<=LastCol());
////        tgt_fc--;
////      }
////      else{
////        first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
////      }
////      assert(first_pivot_idx>=0);
////
////      //find the last row updated by src_snode
////      Int tgt_lc = LastCol();
////      Int last_pivot_idx = -1;
////      //find the pivot idx
////      do {last_pivot_idx = src_snode.FindBlockIdx(tgt_lc); tgt_lc--;}
////      while(last_pivot_idx<0 && tgt_lc>=tgt_fc);
////      tgt_lc++;
////      assert(last_pivot_idx>=0);
////      TIMER_STOP(UPDATE_SNODE_FIND_INDEX);

      Int tgt_fc,tgt_lc;
      Int first_pivot_idx,last_pivot_idx;
      FindUpdatedFirstCol(src_snode, update.src_first_row, tgt_fc, first_pivot_idx);
      FindUpdatedLastCol(src_snode, tgt_fc, first_pivot_idx, tgt_lc, last_pivot_idx);

      NZBlockDesc & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);
      NZBlockDesc & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);



      //determine the first column that will be updated in the target supernode
      Int tgt_local_fc =  tgt_fc - FirstCol();
      Int tgt_local_lc =  tgt_lc - FirstCol();

      Int tgt_nrows = NRowsBelowBlock(0);
      Int src_nrows = src_snode.NRowsBelowBlock(first_pivot_idx)
        - (tgt_fc - first_pivot_desc.GIndex);
      Int src_lr = tgt_fc+src_nrows-1;
      src_nrows = src_lr - tgt_fc + 1;

      Int tgt_width = src_nrows - src_snode.NRowsBelowBlock(last_pivot_idx)
        + (tgt_lc - last_pivot_desc.GIndex)+1;

      T * pivot = src_snode.GetNZval(first_pivot_desc.Offset)
        + (tgt_fc-first_pivot_desc.GIndex)*src_snode_size;
      T * tgt = GetNZval(0);

      //Pointer to the output buffer of the GEMM
      T * buf = NULL;
      T beta = ZERO<T>();
      //If the target supernode has the same structure,
      //The GEMM is directly done in place
      if(src_nrows == tgt_nrows){
        Int tgt_offset = (tgt_fc - FirstCol());
        buf = &tgt[tgt_offset];
        beta = ONE<T>();
      }
      else{
        //Compute the update in a temporary buffer
#ifdef _DEBUG_
        tmpBuffers.tmpBuf.Resize(tgt_width,src_nrows);
#endif

        buf = tmpBuffers.tmpBuf.Data();
      }

      //everything is in row-major
      TIMER_START(UPDATE_SNODE_GEMM);
      blas::Gemm('T','N',tgt_width, src_nrows,src_snode_size,
          MINUS_ONE<T>(),pivot,src_snode_size,
          pivot,src_snode_size,beta,buf,tgt_width);
      TIMER_STOP(UPDATE_SNODE_GEMM);

      //If the GEMM wasn't done in place we need to aggregate the update
      //This is the assembly phase
      if(src_nrows != tgt_nrows){
#ifdef _DEBUG_
        logfileptr->OFS()<<"tmpBuf is "<<tmpBuffers.tmpBuf<<std::endl;
#endif

        //now add the update to the target supernode
        TIMER_START(UPDATE_SNODE_INDEX_MAP);
        if(tgt_snode_size==1){
          Int rowidx = 0;
          Int src_blkcnt = src_snode.NZBlockCnt();
          for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
            NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
            Int cur_src_nrows = src_snode.NRows(blkidx);
            Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
            Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);

            Int row = cur_src_fr;
            while(row<=cur_src_lr){
              Int tgt_blk_idx = FindBlockIdx(row);
              assert(tgt_blk_idx>=0);
              NZBlockDesc & cur_tgt_desc = GetNZBlockDesc(tgt_blk_idx);
              Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + NRows(tgt_blk_idx)-1);
              Int tgtOffset = cur_tgt_desc.Offset 
                + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
              for(Int cr = row ;cr<=lr;++cr){
                tgt[tgtOffset + (cr - row)*tgt_snode_size] += buf[rowidx]; 
                rowidx++;
              }
              row += (lr-row+1);
            }
          }
        }
        else{
          tmpBuffers.src_colindx.Resize(tgt_width);
          tmpBuffers.src_to_tgt_offset.Resize(src_nrows);

          Int colidx = 0;
          Int rowidx = 0;
          Int offset = 0;


          Int src_blkcnt = src_snode.NZBlockCnt();
          for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
            NZBlockDesc & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
            Int cur_src_nrows = src_snode.NRows(blkidx);
            Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
            Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
            cur_src_nrows = cur_src_lr - cur_src_fr +1;

            //The other one MUST reside into a single block in the target
            Int row = cur_src_fr;
            while(row<=cur_src_lr){
              Int tgt_blk_idx = FindBlockIdx(row);
              assert(tgt_blk_idx>=0);
              NZBlockDesc & cur_tgt_desc = GetNZBlockDesc(tgt_blk_idx);
              Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + NRows(tgt_blk_idx)-1);
              Int tgtOffset = cur_tgt_desc.Offset 
                + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
              for(Int cr = row ;cr<=lr;++cr){
                if(cr<=tgt_lc){
                  tmpBuffers.src_colindx[colidx++] = cr;
                }
                offset+=tgt_width;
                tmpBuffers.src_to_tgt_offset[rowidx] = tgtOffset + (cr - row)*tgt_snode_size;
                rowidx++;
              }
              row += (lr-row+1);
            }
          }


          //Multiple cases to consider
          TIMER_STOP(UPDATE_SNODE_INDEX_MAP);

          if(first_pivot_idx==last_pivot_idx){
            // Updating contiguous columns
            Int tgt_offset = (tgt_fc - FirstCol());
            for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
              blas::Axpy(tgt_width,ONE<T>(),&buf[rowidx*tgt_width],1,
                  &tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_offset],1);
            }
          }
          else{
            // full sparse case (done right now)
            for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
              for(Int colidx = 0; colidx< tmpBuffers.src_colindx.m();++colidx){
                Int col = tmpBuffers.src_colindx[colidx];
                Int tgt_colidx = col - FirstCol();
                tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_colidx] 
                  += buf[rowidx*tgt_width+colidx]; 
              }
            }
          }
        }
      }
      TIMER_STOP(UPDATE_SNODE);
      return 0;
    }
#else
  template<typename T>
    inline Int SuperNode<T>::Update(SuperNode<T> & src_snode, SnodeUpdate &update, 
        TempUpdateBuffers<T> & tmpBuffers){
    return trc_Update(src_snode, update, tmpBuffers, *this);
  }
#endif


  template<typename T>
    inline Int SuperNode<T>::Factorize(){
#ifndef _TAU_TRACE_
      Int BLOCKSIZE = Size();
      NZBlockDesc & diag_desc = GetNZBlockDesc(0);
      for(Int col = 0; col<Size();col+=BLOCKSIZE){
        Int bw = min(BLOCKSIZE,Size()-col);
        T * diag_nzval = &GetNZval(diag_desc.Offset)[col+col*Size()];
        lapack::Potrf( 'U', bw, diag_nzval, Size());
        T * nzblk_nzval = &GetNZval(diag_desc.Offset)[col+(col+bw)*Size()];
        blas::Trsm('L','U','T','N',bw, NRowsBelowBlock(0)-(col+bw), ONE<T>(),  diag_nzval, Size(), nzblk_nzval, Size());

        //update the rest !!! (next blocks columns)
        T * tgt_nzval = &GetNZval(diag_desc.Offset)[col+bw+(col+bw)*Size()];
        blas::Gemm('T','N',Size()-(col+bw), NRowsBelowBlock(0)-(col+bw),bw,MINUS_ONE<T>(),nzblk_nzval,Size(),nzblk_nzval,Size(),ONE<T>(),tgt_nzval,Size());
      }
      return 0;

      //          T * diag_nzval = GetNZval(diag_desc.Offset);
      //          lapack::Potrf( 'U', Size(), diag_nzval, Size());
      //        if(NZBlockCnt()>1){
      //          NZBlockDesc & nzblk_desc = GetNZBlockDesc(1);
      //          T * nzblk_nzval = GetNZval(nzblk_desc.Offset);
      //          blas::Trsm('L','U','T','N',Size(), NRowsBelowBlock(1), ONE<T>(),  diag_nzval, Size(), nzblk_nzval, Size());
      //        }
#else
  return trc_Factorize(*this);
#endif


    }




  template<typename T>
    bool SuperNode<T>::FindNextUpdate(SnodeUpdate & nextUpdate, const IntNumVec & Xsuper,  const IntNumVec & SupMembership, bool isLocal){
      Int & tgt_snode_id = nextUpdate.tgt_snode_id;
      Int & f_ur = nextUpdate.src_first_row;
      Int & f_ub = nextUpdate.blkidx;
      Int & n_ur = nextUpdate.src_next_row;
      Int & n_ub = nextUpdate.next_blkidx;


      if(tgt_snode_id==0){   
        f_ub = isLocal?1:0;
        n_ub = f_ub;
        n_ur = 0;
      }
      else{
        f_ub = n_ub;
      }

      if(NZBlockCnt()>f_ub){
        NZBlockDesc * cur_desc = &GetNZBlockDesc(f_ub); 
        f_ur = max(n_ur,cur_desc->GIndex); 
        //find the snode updated by that row
        tgt_snode_id = SupMembership[f_ur-1];
        Int tgt_lc = Xsuper[tgt_snode_id]-1;

        
        Int src_tgt_lb = -1;
#ifdef _LINEAR_SEARCH_FCLC_
    if(ITreeInitialized()){
#endif 
        //or use FindBlockIdx
        src_tgt_lb = FindBlockIdx(tgt_lc);
#ifdef _LINEAR_SEARCH_FCLC_
    }
    else{
        for(Int blkidx = n_ub; blkidx<NZBlockCnt();++blkidx){
          NZBlockDesc & cur_block = GetNZBlockDesc(blkidx);
          Int cur_fr = cur_block.GIndex;
          Int cur_lr = cur_block.GIndex + NRows(blkidx)-1;

          if(cur_fr<= tgt_lc && tgt_lc<=cur_lr){
            src_tgt_lb = blkidx;
            break;
          }
        }
    }
#endif
        //if tgt_lc not found in the current column we need to resume at the next block 
        if(src_tgt_lb==-1){
          for(n_ub;n_ub<NZBlockCnt();++n_ub){
            cur_desc = &GetNZBlockDesc(n_ub);
            if(cur_desc->GIndex > tgt_lc){
              break;
            }
          }
          if(n_ub<NZBlockCnt()){
            n_ur = cur_desc->GIndex;
          }
          else{
            n_ur = -1;
          }
        }
        else{
          n_ub = src_tgt_lb;
          cur_desc = &GetNZBlockDesc(n_ub);
          if(cur_desc->GIndex + NRows(n_ub)-1>tgt_lc){
            n_ur = tgt_lc+1;
          }
          else{
            ++n_ub;
            n_ur = (n_ub<NZBlockCnt())?GetNZBlockDesc(n_ub).GIndex:-1;
          }
        }

        //src_snode updates tgt_snode_id. Then we need to look from row n_ur and block l_ub
        nextUpdate.src_snode_id = Id();
        return true;
      }
      else{
        nextUpdate.src_snode_id = -1;
        return false;
      }
    }


  template<typename T>
  inline ITree * SuperNode<T>::CreateITree(){
    #if defined(_AVL_ITREE_)
      return new AVLITree();
    #elif defined(_DSW_ITREE_)
      return new DSWITree();
    #else
      return new ITree();
    #endif
  } 

  template<typename T>
  inline void SuperNode<T>::InitIdxToBlk(){
//#ifdef INTERVAL_TREE
        for(Int blkidx=0; blkidx<blocks_cnt_;++blkidx){
          Int cur_fr = blocks_[blkidx].GIndex;
          Int cur_lr = cur_fr + NRows(blkidx) -1;

          ITree::Interval cur_interv = { cur_fr, cur_lr, blkidx};
          idxToBlk_->Insert(cur_interv);
        }
      idxToBlk_->Rebalance();
//#endif
  }










  template <typename T> inline std::ostream& operator<<( std::ostream& os,  SuperNode<T>& snode){
    os<<"ooooooooooo   Supernode "<<snode.Id()<<" oooooooooooo"<<std::endl;
    os<<"     size = "<<snode.Size()<<std::endl;
    os<<"     fc   = "<<snode.FirstCol()<<std::endl;
    os<<"     lc   = "<<snode.LastCol()<<std::endl;
    for(Int blkidx =0; blkidx<snode.NZBlockCnt();++blkidx){
      NZBlockDesc & nzblk_desc = snode.GetNZBlockDesc(blkidx);
      T * nzblk_nzval = snode.GetNZval(nzblk_desc.Offset);
      os<<"--- NZBlock "<<nzblk_desc.GIndex<<" ---"<<std::endl;
      for(Int i = 0; i<snode.NRows(blkidx);++i){
        for(Int j = 0; j<snode.Size();++j){
          os<<nzblk_nzval[i*snode.Size()+j]<<" ";
        }
        os<<std::endl;
      }
    }
    os<<"oooooooooooooooooooooooooooooooooooooooo"<<std::endl;
    return os;
  }


  template <typename T> inline void Serialize(Icomm & buffer,SuperNode<T> & snode){
    Int nzblk_cnt = snode.NZBlockCnt();
    Int nzval_cnt = snode.Size()*snode.NRowsBelowBlock(0);
    T* nzval_ptr = snode.GetNZval(0);
    NZBlockDesc* nzblk_ptr = &snode.GetNZBlockDesc(0);
    Int snode_id = snode.Id();
    Int snode_fc = snode.FirstCol();
    Int snode_lc = snode.LastCol();

    buffer.clear();
    buffer.resize(5*sizeof(Int) + nzblk_cnt*sizeof(NZBlockDesc) + nzval_cnt*sizeof(T));

    buffer<<snode_id;
    buffer<<snode_fc;
    buffer<<snode_lc;
    buffer<<nzblk_cnt;
    Serialize(buffer,nzblk_ptr,nzblk_cnt);
    buffer<<nzval_cnt;
    Serialize(buffer,nzval_ptr,nzval_cnt);
  }

  template <typename T> inline void Serialize(Icomm & buffer,SuperNode<T> & snode, Int first_blkidx, Int first_row){
    NZBlockDesc* nzblk_ptr = &snode.GetNZBlockDesc(first_blkidx);
    Int local_first_row = first_row - nzblk_ptr->GIndex;
    Int nzblk_cnt = snode.NZBlockCnt() - first_blkidx;
    Int nzval_cnt = snode.Size()*(snode.NRowsBelowBlock(first_blkidx)-local_first_row);
    T* nzval_ptr = snode.GetNZval(nzblk_ptr->Offset) + local_first_row*snode.Size();
    Int snode_id = snode.Id();
    Int snode_fc = snode.FirstCol();
    Int snode_lc = snode.LastCol();

    buffer.clear();
    buffer.resize(5*sizeof(Int) + nzblk_cnt*sizeof(NZBlockDesc) + nzval_cnt*sizeof(T));

    buffer<<snode_id;
    buffer<<snode_fc;
    buffer<<snode_lc;
    buffer<<nzblk_cnt;
    NZBlockDesc * new_nzblk_ptr = reinterpret_cast<NZBlockDesc *>(buffer.back());
    Serialize(buffer,nzblk_ptr,nzblk_cnt);
    //replace the GIndex of the serialized block descriptor
    // by the new first row, and update the offset appropriately
    new_nzblk_ptr->GIndex = first_row;
    new_nzblk_ptr->Offset = nzblk_ptr->Offset + local_first_row*snode.Size();

    buffer<<nzval_cnt;
    Serialize(buffer,nzval_ptr,nzval_cnt);
  }


  template <typename T> inline void Serialize(Icomm & buffer,SuperNode<T> & snode, Int first_blkidx, Int first_row, size_t extra_bytespace){
    NZBlockDesc* nzblk_ptr = &snode.GetNZBlockDesc(first_blkidx);
    Int local_first_row = first_row - nzblk_ptr->GIndex;
    Int nzblk_cnt = snode.NZBlockCnt() - first_blkidx;
    Int nzval_cnt = snode.Size()*(snode.NRowsBelowBlock(first_blkidx)-local_first_row);
    T* nzval_ptr = snode.GetNZval(nzblk_ptr->Offset) + local_first_row*snode.Size();
    Int snode_id = snode.Id();
    Int snode_fc = snode.FirstCol();
    Int snode_lc = snode.LastCol();

    buffer.clear();
    buffer.resize(5*sizeof(Int) + nzblk_cnt*sizeof(NZBlockDesc) + nzval_cnt*sizeof(T) + extra_bytespace);

    buffer<<snode_id;
    buffer<<snode_fc;
    buffer<<snode_lc;
    buffer<<nzblk_cnt;
    NZBlockDesc * new_nzblk_ptr = reinterpret_cast<NZBlockDesc *>(buffer.back());
    Serialize(buffer,nzblk_ptr,nzblk_cnt);
    //replace the GIndex of the serialized block descriptor
    // by the new first row, and update the offset appropriately
    new_nzblk_ptr->GIndex = first_row;
    new_nzblk_ptr->Offset = nzblk_ptr->Offset + local_first_row*snode.Size();

    buffer<<nzval_cnt;
    Serialize(buffer,nzval_ptr,nzval_cnt);
  }




  template <typename T> inline size_t Deserialize(char * buffer, SuperNode<T> & snode){
    Int snode_id = *(Int*)&buffer[0];
    Int snode_fc = *(((Int*)&buffer[0])+1);
    Int snode_lc = *(((Int*)&buffer[0])+2);
    Int nzblk_cnt = *(((Int*)&buffer[0])+3);
    NZBlockDesc * blocks_ptr = 
      reinterpret_cast<NZBlockDesc*>(&buffer[4*sizeof(Int)]);
    Int nzval_cnt = *(Int*)(blocks_ptr + nzblk_cnt);
    T * nzval_ptr = (T*)((Int*)(blocks_ptr + nzblk_cnt)+1);
    char * last_ptr = (char *)(nzval_ptr+nzval_cnt);

    //Create the dummy supernode for that data
    snode.Init(snode_id,snode_fc,snode_lc, blocks_ptr, nzblk_cnt, nzval_ptr, nzval_cnt);

    return (last_ptr - buffer);
  }


//}

#endif // _SUPERNODE_IMPL_HPP_
