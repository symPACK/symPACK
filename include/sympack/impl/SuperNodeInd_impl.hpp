#ifndef _SUPERNODE_IND_IMPL_HPP_
#define _SUPERNODE_IND_IMPL_HPP_

#ifdef _USE_COREDUMPER_
#include <google/coredumper.h>
#endif



/****************************************/
/*            _________________         */
/*           |     nzval       |        */
/*           |                 |        */
/*           |_________________|        */
/*           |                 |        */
/*           |    updrows      |        */
/*           |                 |        */
/*           |_________________|        */
/*           |      Diag       |        */
/*           |_________________|        */
/*           |      Meta       |        */
/*           |_________________|        */
/*           |                 |        */
/*           |   block desc    |        */
/*           |                 |        */
/*           |_________________|        */
/*                                      */
/****************************************/

namespace symPACK{

  template<typename T, class Allocator>
    SuperNodeInd<T,Allocator>::SuperNodeInd() :SuperNode<T,Allocator>(), diag_(nullptr){ }

  template<typename T, class Allocator>
    SuperNodeInd<T,Allocator>::SuperNodeInd(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int ai_num_rows, Int aiN, Int aiNZBlkCnt, Int panel):SuperNodeInd<T,Allocator>() {
      //compute supernode size / width
      Int size = aiLc - aiFc +1;
      //compute maximum number of blocks, number of off-diagonal rows + 1
      Int num_blocks = std::max((Int)1,ai_num_rows-size + 1);
      if(aiNZBlkCnt!=-1){
        num_blocks=aiNZBlkCnt;
      }


      if ( panel > 0 ) {
        Int numPanels = std::ceil((double)size/(double)panel);
        Int lastPanel = std::max(size - (numPanels-1)*panel,0);
        this->storage_size_ = sizeof(T)*((numPanels-1)*ai_num_rows*panel - 
            - (numPanels-1)*(numPanels-2)*panel/2 
            + (ai_num_rows-(numPanels-1)*panel)*lastPanel);
      }
      else {
        this->storage_size_  = sizeof(T)*size*ai_num_rows; //size of nzvals
      }




      this->storage_size_ += num_blocks*sizeof(NZBlockDesc); //size of block desc
      this->storage_size_ += sizeof(SuperNodeDesc); //size of the supernode metadata
      this->storage_size_ += sizeof(T)*size; //extra space for diagonal (indefinite matrices)

      Int num_updrows = 0;
      if (aiFr == aiFc){
        num_updrows = num_blocks;
      }


      try{
        this->loc_storage_container_ = Allocator::allocate(this->storage_size_);
      }
      catch(const MemoryAllocationException & e){
        this->loc_storage_container_=nullptr;
        this->storage_size_ = 0;
        throw;
      }
      bassert(this->loc_storage_container_!=nullptr);

      this->nzval_ = (T*)&this->loc_storage_container_[0];
      this->diag_ = (T*)(this->nzval_+size*ai_num_rows);
      this->meta_ = (SuperNodeDesc*)(this->diag_ + size);
      char * last = this->loc_storage_container_+this->storage_size_-1 - (sizeof(NZBlockDesc) -1);
      this->blocks_ = (NZBlockDesc*) last;

      this->meta_->panel_sz_ = panel;
      this->meta_->iId_ = aiId;
      this->meta_->iFirstRow_ = aiFr;
      this->meta_->iFirstCol_ = aiFc;
      this->meta_->iLastCol_ = aiLc;
      this->meta_->iN_=aiN;
      this->meta_->iSize_ = size;
      this->meta_->nzval_cnt_ = 0;
      this->meta_->blocks_cnt_ = 0;
      this->meta_->b_own_storage_ = true;

#ifndef ITREE
      globalToLocal_ = new std::vector<Int>(aiN+1,-1);
#else
      this->idxToBlk_ = this->CreateITree();
#endif
    } 

  template<typename T, class Allocator>
    SuperNodeInd<T,Allocator>::SuperNodeInd(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices, Int panel):SuperNodeInd<T,Allocator>() {
      Init(aiId,aiFr,aiFc,aiLc,aiN,rowIndices,panel);
    } 

  template<typename T, class Allocator>
    SuperNodeInd<T,Allocator>::SuperNodeInd(char * storage_ptr,size_t storage_size, Int GIndex ):SuperNodeInd<T,Allocator>() {
      this->Init(storage_ptr,storage_size,GIndex);
    }

  //CHECKED ON 11-18-2014
  template<typename T, class Allocator>
    void SuperNodeInd<T,Allocator>::Init(char * storage_ptr,size_t storage_size, Int GIndex ) {


      //loop through the block descriptors
      char * last = (char*)(storage_ptr+storage_size-1) - (sizeof(NZBlockDesc) -1);

      this->blocks_ = (NZBlockDesc*) last;

      Int blkCnt = 0;
      NZBlockDesc * curBlockPtr = nullptr;
      do{
        curBlockPtr = &this->GetNZBlockDesc(blkCnt);
        ++blkCnt;
      }
      while(!curBlockPtr->Last);


      this->meta_= (SuperNodeDesc*)(this->blocks_ - blkCnt +1) - 1;
      this->diag_= (T*)(this->meta_) - this->Size();
      this->nzval_ = (T*) storage_ptr;
      //we now need to update the meta data
      this->meta_->b_own_storage_ = false;
      this->meta_->blocks_cnt_ = blkCnt;
      this->meta_->nzval_cnt_ = (T*)this->diag_ - this->nzval_;

      if(GIndex==-1){
        GIndex = this->GetNZBlockDesc(0).GIndex;
      }

      Int loc_fr = GIndex - this->GetNZBlockDesc(0).GIndex;
      Int offset = this->GetNZBlockDesc(0).Offset + loc_fr*this->Size();
      //restore 0-based offsets and compute global_to_local indices
      this->GetNZBlockDesc(0).Offset = 0;
      for(Int blkidx=1; blkidx<this->NZBlockCnt();++blkidx){
        curBlockPtr = &this->GetNZBlockDesc(blkidx);
        curBlockPtr->Offset -= offset;
      }



      this->blocks_->GIndex = GIndex;


#ifndef ITREE
      globalToLocal_ = new std::vector<Int>(this->meta_->iN_+1,-1);
#else
      this->idxToBlk_ = this->CreateITree();
#endif

    }


  template<typename T, class Allocator>
    void SuperNodeInd<T,Allocator>::Init(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices, Int panel){
      //compute supernode size / width
      Int size = aiLc - aiFc +1;

      Int num_blocks = 0;
      if(rowIndices.size()>0){
        //go through the set and count the number of nz blocks
        Idx prevRow = *rowIndices.begin();
        Idx firstRow = *rowIndices.begin();
        for(auto it = rowIndices.begin();it!=rowIndices.end();it++){
          Idx row = *it;

          if(row>prevRow+1){
            num_blocks++;
            firstRow = row;
          }
          prevRow = row;
        }
        num_blocks++;
      }

      assert(num_blocks>0);

      Int numRows = rowIndices.size();
      if ( panel > 0 ) {
        Int numPanels = std::ceil((double)size/(double)panel);
        Int lastPanel = std::max(size - (numPanels-1)*panel,0);
        this->storage_size_ = sizeof(T)*((numPanels-1)*numRows*panel - 
            - (numPanels-1)*(numPanels-2)*panel/2 
            + (numRows-(numPanels-1)*panel)*lastPanel);
        this->storage_size_ += num_blocks*sizeof(NZBlockDesc) + sizeof(SuperNodeDesc);
      }
      else {
        this->storage_size_ = sizeof(T)*size*numRows + num_blocks*sizeof(NZBlockDesc) + sizeof(SuperNodeDesc);
      }

      this->storage_size_ += sizeof(T)*size; //extra space for diagonal (indefinite matrices)

      Int num_updrows = 0;
      if (aiFr == aiFc){
        num_updrows = num_blocks;
      }


      try{
        this->loc_storage_container_ = Allocator::allocate(this->storage_size_);
      }
      catch(const MemoryAllocationException & e){
        this->loc_storage_container_=nullptr;
        this->storage_size_ = 0;
        throw;
      }

      this->nzval_ = (T*)&this->loc_storage_container_[0];
      this->diag_ = (T*)(this->nzval_+size*numRows);
      this->meta_ = (SuperNodeDesc*)(this->diag_ + size);
      char * last = this->loc_storage_container_+this->storage_size_-1 - (sizeof(NZBlockDesc) -1);
      this->blocks_ = (NZBlockDesc*) last;

      this->meta_->panel_sz_ = panel;
      this->meta_->iId_ = aiId;
      this->meta_->iFirstRow_ = aiFr;
      this->meta_->iFirstCol_ = aiFc;
      this->meta_->iLastCol_ = aiLc;
      this->meta_->iN_=aiN;
      this->meta_->iSize_ = size;
      this->meta_->nzval_cnt_ = 0;
      this->meta_->blocks_cnt_ = 0;
      this->meta_->b_own_storage_ = true;


#ifndef ITREE
      this->globalToLocal_ = new std::vector<Int>(aiN+1,-1);
#else
      this->idxToBlk_ = this->CreateITree();
#endif

      //now add the blocks 
      if(rowIndices.size()>0){
        //go through the set and count the number of nz blocks
        Idx prevRow = *rowIndices.begin();
        Idx firstRow = *rowIndices.begin();
        for(auto it = rowIndices.begin();it!=rowIndices.end();it++){
          Idx row = *it;

          if(row>prevRow+1){
            this->AddNZBlock( prevRow - firstRow + 1, firstRow);
            firstRow = row;
          }
          prevRow = row;
        }
        this->AddNZBlock( prevRow - firstRow + 1, firstRow);
      }
    } 

  template<typename T, class Allocator>
    inline void SuperNodeInd<T,Allocator>::AddNZBlock(Int aiNRows, Int aiGIndex){

      //Resize the container if I own the storage
      if(this->meta_->b_own_storage_){
        logfileptr->OFS()<<"Adding\n";
        scope_timer(a,RESIZE_SUPERNODE);

        Int cur_fr = aiGIndex;
        Int cur_lr = cur_fr + aiNRows -1;
        Int cur_nzval_cnt = aiNRows*this->meta_->iSize_;
        if ( cur_fr == this->FirstCol() && this->meta_->panel_sz_ > 0 ) {
          Int panel = this->meta_->panel_sz_;
          Int numPanels = std::ceil((double)this->meta_->iSize_/(double)panel);
          Int lastPanel = std::max(this->meta_->iSize_ - (numPanels-1)*panel,0);
          cur_nzval_cnt = (numPanels-1)*aiNRows*panel - 
            - (numPanels-1)*(numPanels-2)*panel/2 
            + (aiNRows-(numPanels-1)*panel)*lastPanel;
        }



#ifndef ITREE
        std::fill(&(*globalToLocal_)[cur_fr],&(*globalToLocal_)[cur_lr]+1,this->meta_->blocks_cnt_); 
#else
        ITree<Int>::Interval<Int> cur_interv;
        cur_interv.low = cur_fr;
        cur_interv.high = cur_lr;
        cur_interv.data = this->meta_->blocks_cnt_;
        this->idxToBlk_->Insert(cur_interv);
#endif

        bool ownDiagonal = this->OwnDiagonal();
        //if there is no more room for either nzval or blocks, extend
        Int block_space = (Int)(this->blocks_+1 - (NZBlockDesc*)(this->meta_ +1)) - this->meta_->blocks_cnt_;
        size_t nzval_space = ((size_t)((char*)this->diag_ - (char*)this->nzval_) - this->meta_->nzval_cnt_*sizeof(T) );

        if(block_space==0 || nzval_space<cur_nzval_cnt*sizeof(T)){
          //need to resize storage space. this is expensive !
          Int size = this->Size();
          size_t extra_nzvals_bytes = std::max((size_t)0,cur_nzval_cnt*sizeof(T) - nzval_space);
          Int extra_blocks = std::max((Int)0,1 - block_space);
          size_t new_size = this->storage_size_ + extra_nzvals_bytes + extra_blocks*sizeof(NZBlockDesc);

          size_t offset_diag = (char*)this->diag_ - (char*)this->nzval_;
          size_t offset_meta = (char*)this->meta_ - (char*)this->nzval_;
          size_t offset_block = (char*)this->blocks_ - (char*)this->nzval_;

          char * locTmpPtr = Allocator::allocate(new_size);
          bassert(locTmpPtr!=nullptr);

          std::copy(this->loc_storage_container_,this->loc_storage_container_+this->storage_size_,locTmpPtr);


          Allocator::deallocate(this->loc_storage_container_);
          this->loc_storage_container_=locTmpPtr;
          this->storage_size_=new_size;

          this->nzval_=(T*)&this->loc_storage_container_[0];
          //move the block descriptors if required
          char * cur_blocks_ptr = (char*)&this->loc_storage_container_[0] + offset_block;
          //move the meta data if required
          char * cur_meta_ptr = (char*)&this->loc_storage_container_[0] + offset_meta;

          //move the diag entries if required
          char * cur_diag_ptr = (char*)&this->loc_storage_container_[0] + offset_diag;

          this->meta_ = (SuperNodeDesc*) cur_meta_ptr;
          //we need to move everything, starting from the blocks, then meta
          //blocks need to be moved by extra_nzvals_bytes + extra_blocks
          char * new_blocks_ptr = cur_blocks_ptr + extra_nzvals_bytes +extra_blocks*sizeof(NZBlockDesc);
          std::copy_backward(cur_blocks_ptr - (this->meta_->blocks_cnt_-1)*sizeof(NZBlockDesc),cur_blocks_ptr+sizeof(NZBlockDesc),new_blocks_ptr+sizeof(NZBlockDesc));

          //now move the meta data by extra_nzvals_bytes
          char * new_meta_ptr = cur_meta_ptr + extra_nzvals_bytes;
          std::copy(cur_meta_ptr,cur_meta_ptr + sizeof(SuperNodeDesc),new_meta_ptr);

          //now move the diag entries by extra_nzvals_bytes
          char * new_diag_ptr = cur_diag_ptr + extra_nzvals_bytes;
          std::copy(cur_diag_ptr,cur_diag_ptr + sizeof(T)*size,new_diag_ptr);

          //update pointers
          this->diag_ = (T*) new_diag_ptr;
          this->meta_ = (SuperNodeDesc*) new_meta_ptr;
          this->blocks_ = (NZBlockDesc*) new_blocks_ptr;
        }

        this->GetNZBlockDesc(this->meta_->blocks_cnt_).GIndex = aiGIndex;
        this->GetNZBlockDesc(this->meta_->blocks_cnt_).Offset = this->meta_->nzval_cnt_;
        this->GetNZBlockDesc(this->meta_->blocks_cnt_).Last = true;
        if(this->meta_->blocks_cnt_>0){
          this->GetNZBlockDesc(this->meta_->blocks_cnt_-1).Last = false;
        }

        this->meta_->blocks_cnt_++;

        //fill the new block with zeros
        std::fill(this->nzval_+this->meta_->nzval_cnt_,this->nzval_+this->meta_->nzval_cnt_+cur_nzval_cnt,ZERO<T>());

        //update nzval count
        this->meta_->nzval_cnt_+=cur_nzval_cnt;

      }
    }

  template<typename T, class Allocator>
    inline Int SuperNodeInd<T,Allocator>::Shrink(){
      if(this->meta_->b_own_storage_){
        logfileptr->OFS()<<"Shrinking\n";
        //TODO make sure that we do not have any extra space anywhere.

        bool ownDiagonal = this->OwnDiagonal();
        //if there is too much room for either nzval or blocks, contract
        Int block_space = (Int)(this->blocks_+1 - (NZBlockDesc*)(this->meta_ +1)) - this->meta_->blocks_cnt_;
        size_t nzval_space = ((size_t)((char*)this->diag_ - (char*)this->nzval_) - this->meta_->nzval_cnt_*sizeof(T) );

        if(block_space >0 || nzval_space >0){

          size_t new_size = this->storage_size_ - nzval_space - block_space*sizeof(NZBlockDesc);

          size_t offset_diag = this->meta_->nzval_cnt_*sizeof(T);
          size_t offset_meta = (this->meta_->nzval_cnt_+this->meta_->iSize_)*sizeof(T);
          size_t offset_block = (this->meta_->nzval_cnt_+this->meta_->iSize_)*sizeof(T)+sizeof(SuperNodeDesc);

          char * locTmpPtr = Allocator::allocate(new_size);
          bassert(locTmpPtr!=nullptr);

          //copy nzvals
          std::copy(this->nzval_,this->nzval_+this->meta_->nzval_cnt_,(T*)locTmpPtr);


          //copy diag
          std::copy(this->diag_,this->diag_+this->Size(),(T*)(locTmpPtr+offset_diag));

          //copy meta
          std::copy(this->meta_,this->meta_+1,(SuperNodeDesc*)(locTmpPtr+offset_meta));
          //copy blocks
          std::copy(this->blocks_+1-this->meta_->blocks_cnt_,this->blocks_+1,(NZBlockDesc*)(locTmpPtr+offset_block));

          Allocator::deallocate(this->loc_storage_container_);
          this->loc_storage_container_=locTmpPtr;
          this->storage_size_=new_size;

          this->nzval_=(T*)&this->loc_storage_container_[0];
          //move the block descriptors if required
          char * new_blocks_ptr = this->loc_storage_container_+this->storage_size_-1 - (sizeof(NZBlockDesc) -1);

          //move the meta data if required
          char * new_meta_ptr = this->loc_storage_container_ + offset_meta;

          //now move the diag entries by extra_nzvals
          char * new_diag_ptr = this->loc_storage_container_ + offset_diag;
          //update pointers
          this->diag_ = (T*) new_diag_ptr;
          this->meta_ = (SuperNodeDesc*) new_meta_ptr;
          this->blocks_ = (NZBlockDesc*) new_blocks_ptr;
        }
      }

      return this->StorageSize();
    }

  template<typename T, class Allocator>
    inline Int SuperNodeInd<T,Allocator>::Factorize( TempUpdateBuffers<T> & tmpBuffers){
#if defined(_NO_COMPUTATION_)
      return 0;
#endif

      Int snodeSize = this->Size();
      NZBlockDesc & diag_desc = this->GetNZBlockDesc(0);
      T * diag_nzval = this->GetNZval(diag_desc.Offset);
      Idx totalNRows = this->NRowsBelowBlock(0);
      Int INFO = 0;
      SYMPACK_TIMER_START(FACT_POTRF_LDL);
      lapack::Potrf_LDL( "U", snodeSize, diag_nzval, snodeSize, tmpBuffers, INFO);
      SYMPACK_TIMER_STOP(FACT_POTRF_LDL);
      //copy the diagonal entries into the diag portion of the supernode
      SYMPACK_TIMER_START(FACT_DIAG_COPY);
#pragma unroll
      for(Int col = 0; col< snodeSize;col++){ this->diag_[col] = diag_nzval[col+ (col)*snodeSize]; }
      SYMPACK_TIMER_STOP(FACT_DIAG_COPY);


      SYMPACK_TIMER_START(FACT_TRSM);
      T * nzblk_nzval = &diag_nzval[snodeSize*snodeSize];
      blas::Trsm('L','U','T','U',snodeSize, totalNRows-snodeSize, T(1),  diag_nzval, snodeSize, nzblk_nzval, snodeSize);
      SYMPACK_TIMER_STOP(FACT_TRSM);

      SYMPACK_TIMER_START(FACT_SCALE);
      //scale column I
      for ( Idx I = 1; I<=snodeSize;I++) {
        blas::Scal( totalNRows-snodeSize, T(1.0)/this->diag_[I-1], &nzblk_nzval[I-1], snodeSize );
      }
      SYMPACK_TIMER_STOP(FACT_SCALE);

      return 0;

    }


  template<typename T, class Allocator>
    inline Int SuperNodeInd<T,Allocator>::Factorize(SuperNode<T,Allocator> * diag_snode, TempUpdateBuffers<T> & tmpBuffers){
      abort();
#if defined(_NO_COMPUTATION_)
      return 0;
#endif
      Int snodeSize = this->Size();
      NZBlockDesc & diag_desc = this->GetNZBlockDesc(0);
      T * diag_nzval = this->GetNZval(diag_desc.Offset);
      Idx totalNRows = this->NRowsBelowBlock(0);
      Int INFO = 0;
      lapack::Potrf_LDL( "U", snodeSize, diag_nzval, snodeSize, tmpBuffers, INFO);
      //copy the diagonal entries into the diag portion of the supernode
#pragma unroll
      for(Int col = 0; col< snodeSize;col++){ this->diag_[col] = diag_nzval[col+ (col)*snodeSize]; }


      T * nzblk_nzval = &diag_nzval[snodeSize*snodeSize];
      blas::Trsm('L','U','T','U',snodeSize, totalNRows-snodeSize, T(1),  diag_nzval, snodeSize, nzblk_nzval, snodeSize);
      //scale column I
      for ( Idx I = 1; I<=snodeSize;I++) {
        blas::Scal( totalNRows-snodeSize, T(1.0)/this->diag_[I-1], &nzblk_nzval[I-1], snodeSize );
      }
      return 0;
    }

  template<typename T, class Allocator>
    inline Int SuperNodeInd<T,Allocator>::UpdateAggregate(SuperNode<T,Allocator> * src_snode, SnodeUpdate &update, TempUpdateBuffers<T> & tmpBuffers, Int iTarget, Int iam ){
      bassert(this->meta_!=nullptr);
      scope_timer(a,UPDATE_AGGREGATE_SNODE);

#if defined(_NO_COMPUTATION_)
      return 0;
#endif
      if(iTarget != iam){
        this->Merge(src_snode, update);

        Int & pivot_idx = update.blkidx;
        Int & pivot_fr = update.src_first_row;

        Int src_snode_size = src_snode->Size();
        Int tgt_snode_size = this->Size();

        Int tgt_fc,tgt_lc;
        Int first_pivot_idx,last_pivot_idx;
        this->FindUpdatedFirstCol(src_snode, update.src_first_row, tgt_fc, first_pivot_idx);
        this->FindUpdatedLastCol(src_snode, tgt_fc, first_pivot_idx, tgt_lc, last_pivot_idx);

        NZBlockDesc & first_pivot_desc = src_snode->GetNZBlockDesc(first_pivot_idx);
        NZBlockDesc & last_pivot_desc = src_snode->GetNZBlockDesc(last_pivot_idx);

        //determine the first column that will be updated in the target supernode
        Int tgt_local_fc =  tgt_fc - this->FirstCol();
        Int tgt_local_lc =  tgt_lc - this->FirstCol();

        Int tgt_nrows = this->NRowsBelowBlock(0);
        Int src_nrows = src_snode->NRowsBelowBlock(first_pivot_idx)
          - (tgt_fc - first_pivot_desc.GIndex);
        Int src_lr = tgt_fc+src_nrows-1;
        src_nrows = src_lr - tgt_fc + 1;

        Int tgt_width = src_nrows - src_snode->NRowsBelowBlock(last_pivot_idx)
          + (tgt_lc - last_pivot_desc.GIndex)+1;

        T * pivot = src_snode->GetNZval(first_pivot_desc.Offset)
          + (tgt_fc-first_pivot_desc.GIndex)*src_snode_size;
        T * tgt = this->GetNZval(0);


        //Pointer to the output buffer of the GEMM
        T * buf = nullptr;
        T beta = ZERO<T>();
        //If the target supernode has the same structure,
        //The GEMM is directly done in place
#ifdef SP_THREADS
        tmpBuffers.tmpBuf.resize(src_snode_size*tgt_width+tgt_width*src_nrows);
#endif
        T * bufLDL = &tmpBuffers.tmpBuf[0];
        buf = bufLDL + src_snode_size*tgt_width;

        SYMPACK_TIMER_START(UPDATE_SNODE_GEMM);

        //everything is in row-major
        //First do W = DLT 
        T * diag = static_cast<SuperNodeInd<T,Allocator> * >(src_snode)->GetDiag();
        for(Int row = 0; row<src_snode_size; row++){
          for(Int col = 0; col<tgt_width; col++){
            bufLDL[col+row*tgt_width] = diag[row]*(pivot[row+col*src_snode_size]);
          }
        }

        //Then do -L*W (gemm)
        blas::Gemm('N','N',tgt_width,src_nrows,src_snode_size,
            T(-1.0),bufLDL,tgt_width,pivot,src_snode_size,beta,buf,tgt_width);
        SYMPACK_TIMER_STOP(UPDATE_SNODE_GEMM);

        //If the GEMM wasn't done in place we need to aggregate the update
        //This is the assembly phase
#ifdef _DEBUG_
        logfileptr->OFS()<<"tmpBuf is "<<tmpBuffers.tmpBuf<<std::endl;
#endif

        //now add the update to the target supernode
        SYMPACK_TIMER_START(UPDATE_SNODE_INDEX_MAP);
        if(tgt_snode_size==1){
          Int rowidx = 0;
          Int src_blkcnt = src_snode->NZBlockCnt();
          for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
            NZBlockDesc & cur_block_desc = src_snode->GetNZBlockDesc(blkidx);
            Int cur_src_nrows = src_snode->NRows(blkidx);
            Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
            Int cur_src_fr = std::max(tgt_fc, cur_block_desc.GIndex);

            Int row = cur_src_fr;
            while(row<=cur_src_lr){
              Int tgt_blk_idx = this->FindBlockIdx(row);
              assert(tgt_blk_idx>=0);
              NZBlockDesc & cur_tgt_desc = this->GetNZBlockDesc(tgt_blk_idx);
              Int lr = std::min(cur_src_lr,cur_tgt_desc.GIndex + this->NRows(tgt_blk_idx)-1);
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
          tmpBuffers.src_colindx.resize(tgt_width);
          tmpBuffers.src_to_tgt_offset.resize(src_nrows);
          Int colidx = 0;
          Int rowidx = 0;
          Int offset = 0;

          Int src_blkcnt = src_snode->NZBlockCnt();
          for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
            NZBlockDesc & cur_block_desc = src_snode->GetNZBlockDesc(blkidx);
            Int cur_src_nrows = src_snode->NRows(blkidx);
            Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
            Int cur_src_fr = std::max(tgt_fc, cur_block_desc.GIndex);
            cur_src_nrows = cur_src_lr - cur_src_fr +1;

            Int row = cur_src_fr;
            while(row<=cur_src_lr){
              Int tgt_blk_idx = this->FindBlockIdx(row);
              assert(tgt_blk_idx>=0);
              NZBlockDesc & cur_tgt_desc = this->GetNZBlockDesc(tgt_blk_idx);
              Int lr = std::min(cur_src_lr,cur_tgt_desc.GIndex + this->NRows(tgt_blk_idx)-1);
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
          SYMPACK_TIMER_STOP(UPDATE_SNODE_INDEX_MAP);


          //Multiple cases to consider
          if(first_pivot_idx==last_pivot_idx){
            // Updating contiguous columns
            Int tgt_offset = (tgt_fc - this->FirstCol());
            for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
              T * A = &buf[rowidx*tgt_width];
              T * B = &tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_offset];
#pragma unroll
              for(Int i = 0; i < tgt_width; ++i){
                B[i] += A[i];
              }
            }
          }
          else{
            // full sparse case (done right now)
            for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
              for(Int colidx = 0; colidx< tmpBuffers.src_colindx.size();++colidx){
                Int col = tmpBuffers.src_colindx[colidx];
                Int tgt_colidx = col - this->FirstCol();
                tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_colidx] 
                  += buf[rowidx*tgt_width+colidx]; 
              }
            }
          }
        }
        return 0;

      }
      else{
        return this->Update(src_snode, update, tmpBuffers);
      }





    }


  template<typename T, class Allocator>
    inline Int SuperNodeInd<T,Allocator>::Update(SuperNode<T,Allocator> * src_snode, SnodeUpdate &update, TempUpdateBuffers<T> & tmpBuffers){
      scope_timer(a,UPDATE_SNODE);
#if defined(_NO_COMPUTATION_)
      return 0;
#endif
      Int & pivot_idx = update.blkidx;
      Int & pivot_fr = update.src_first_row;

      Int src_snode_size = src_snode->Size();
      Int tgt_snode_size = this->Size();

      //find the first row updated by src_snode
      Int tgt_fc,tgt_lc;
      Int first_pivot_idx,last_pivot_idx;
      this->FindUpdatedFirstCol(src_snode, update.src_first_row, tgt_fc, first_pivot_idx);
      this->FindUpdatedLastCol(src_snode, tgt_fc, first_pivot_idx, tgt_lc, last_pivot_idx);

      NZBlockDesc & first_pivot_desc = src_snode->GetNZBlockDesc(first_pivot_idx);
      NZBlockDesc & last_pivot_desc = src_snode->GetNZBlockDesc(last_pivot_idx);

      //determine the first column that will be updated in the target supernode
      Int tgt_local_fc =  tgt_fc - this->FirstCol();
      Int tgt_local_lc =  tgt_lc - this->FirstCol();

      Int tgt_nrows = this->NRowsBelowBlock(0);
      Int src_nrows = src_snode->NRowsBelowBlock(first_pivot_idx)
        - (tgt_fc - first_pivot_desc.GIndex);
      Int src_lr = tgt_fc+src_nrows-1;
      src_nrows = src_lr - tgt_fc + 1;

      Int src_belowLast = src_snode->NRowsBelowBlock(last_pivot_idx);
      Int tgt_width = src_nrows - src_belowLast
        + (tgt_lc - last_pivot_desc.GIndex)+1;

      T * pivot = src_snode->GetNZval(first_pivot_desc.Offset)
        + (tgt_fc-first_pivot_desc.GIndex)*src_snode_size;
      T * tgt = this->GetNZval(0);

      //Pointer to the output buffer of the GEMM
      T * buf = nullptr;
      T beta = ZERO<T>();
      //If the target supernode has the same structure,
      //The GEMM is directly done in place
#ifdef SP_THREADS
      tmpBuffers.tmpBuf.resize(tgt_width*src_nrows + src_snode_size*tgt_width);
#endif


      T * bufLDL = &tmpBuffers.tmpBuf[0];
      if(src_nrows == tgt_nrows){
        Int tgt_offset = (tgt_fc - this->FirstCol());
        buf = &tgt[tgt_offset];
        beta = ONE<T>();
      }
      else{
        //Compute the update in a temporary buffer
        buf = bufLDL + src_snode_size*tgt_width;
      }

      SYMPACK_TIMER_START(UPDATE_SNODE_GEMM);

      //everything is in row-major
      //First do W = DLT 
      T * diag = static_cast<SuperNodeInd<T,Allocator> * >(src_snode)->GetDiag();
      for(Int row = 0; row<src_snode_size; row++){
        for(Int col = 0; col<tgt_width; col++){
          bufLDL[col+row*tgt_width] = diag[row]*(pivot[row+col*src_snode_size]);
        }
      }

      //Then do -L*W (gemm)
      blas::Gemm('N','N',tgt_width,src_nrows,src_snode_size,
          T(-1.0),bufLDL,tgt_width,pivot,src_snode_size,beta,buf,tgt_width);
      SYMPACK_TIMER_STOP(UPDATE_SNODE_GEMM);

      //If the GEMM wasn't done in place we need to aggregate the update
      //This is the assembly phase
      if(src_nrows != tgt_nrows){
#ifdef _DEBUG_
        logfileptr->OFS()<<"tmpBuf is "<<tmpBuffers.tmpBuf<<std::endl;
#endif

        //now add the update to the target supernode
        SYMPACK_TIMER_START(UPDATE_SNODE_INDEX_MAP);
        if(tgt_snode_size==1){
          Int rowidx = 0;
          Int src_blkcnt = src_snode->NZBlockCnt();
          {
            for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
              NZBlockDesc & cur_block_desc = src_snode->GetNZBlockDesc(blkidx);
              Int cur_src_nrows = src_snode->NRows(blkidx);
              Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
              Int cur_src_fr = std::max(tgt_fc, cur_block_desc.GIndex);

              Int row = cur_src_fr;
              while(row<=cur_src_lr){
                Int tgt_blk_idx = this->FindBlockIdx(row);
                if(tgt_blk_idx<0){src_snode->DumpITree(); this->DumpITree(); 

                  logfileptr->OFS()<<"LOCK 3: tgt_blk_idx<0"<<std::endl;
                  gdb_lock();}
                assert(tgt_blk_idx>=0);
                NZBlockDesc & cur_tgt_desc = this->GetNZBlockDesc(tgt_blk_idx);
                Int lr = std::min(cur_src_lr,cur_tgt_desc.GIndex + this->NRows(tgt_blk_idx)-1);
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
        }
        else{
          tmpBuffers.src_colindx.resize(tgt_width);
          tmpBuffers.src_to_tgt_offset.resize(src_nrows);
          Int colidx = 0;
          Int rowidx = 0;
          Int offset = 0;

          Int src_blkcnt = src_snode->NZBlockCnt();

          {
            for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
              NZBlockDesc & cur_block_desc = src_snode->GetNZBlockDesc(blkidx);
              Int cur_src_nrows = src_snode->NRows(blkidx);
              Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
              Int cur_src_fr = std::max(tgt_fc, cur_block_desc.GIndex);
              cur_src_nrows = cur_src_lr - cur_src_fr +1;

              //The other one MUST reside into a single block in the target
              Int row = cur_src_fr;
              while(row<=cur_src_lr){
                Int tgt_blk_idx = this->FindBlockIdx(row);
                if(tgt_blk_idx<0){
                  logfileptr->OFS()<<"LOCK 4: tgt_blk_idx<0"<<std::endl;
                  gdb_lock();}
                assert(tgt_blk_idx>=0);
                NZBlockDesc & cur_tgt_desc = this->GetNZBlockDesc(tgt_blk_idx);
                Int lr = std::min(cur_src_lr,cur_tgt_desc.GIndex + this->NRows(tgt_blk_idx)-1);
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
          }
          //Multiple cases to consider
          SYMPACK_TIMER_STOP(UPDATE_SNODE_INDEX_MAP);

          if(first_pivot_idx==last_pivot_idx){
            // Updating contiguous columns
            Int tgt_offset = (tgt_fc - this->FirstCol());

            {
              for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
                T * A = &buf[rowidx*tgt_width];
                T * B = &tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_offset];
#pragma unroll
                for(Int i = 0; i < tgt_width; ++i){
                  B[i] += A[i];
                }
              }
            }
          }
          else{
            // full sparse case
            {
              for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
                for(Int colidx = 0; colidx< tmpBuffers.src_colindx.size();++colidx){
                  Int col = tmpBuffers.src_colindx[colidx];
                  Int tgt_colidx = col - this->FirstCol();
                  tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_colidx] 
                    += buf[rowidx*tgt_width+colidx]; 
                }
              }
            }
          }
        }
      }
      return 0;
    }





  template <typename T, class Allocator> 
    inline void SuperNodeInd<T,Allocator>::back_update_contrib(SuperNode<T> * cur_snode, Int nrhsOffset, Int pnrhs){
      SuperNodeInd<T,Allocator> * contrib = this;
      SuperNodeInd<T> * cur_snodeInd = static_cast<SuperNodeInd<T>*>(cur_snode);

      Int nrhs = pnrhs;
      if(pnrhs==-1){
        nrhs = contrib->Size();
        assert(nrhsOffset==0);
      }
      else{
        nrhs = std::min(pnrhs,contrib->Size() - nrhsOffset);
      }

      Int ldsol = contrib->Size();
      Int ldfact = cur_snode->Size();

      Int nzblk_cnt = cur_snode->NZBlockCnt();
      if(nzblk_cnt>=1){
        //for each row (of the diagonal block)
        T * updated_nzval = contrib->GetNZval(0);

        Int snodeSize = cur_snode->Size();
        for (Int kk = snodeSize-1; kk>=0; --kk){
          for (Int blkidx = 0; blkidx<nzblk_cnt;++blkidx){
            Int offset = blkidx==0?1:0;
            Int rowK = blkidx==0?kk+offset:0;


            NZBlockDesc & cur_desc = contrib->GetNZBlockDesc(blkidx);
            NZBlockDesc & fact_desc = cur_snodeInd->GetNZBlockDesc(blkidx);

            Int cur_nrows = contrib->NRows(blkidx) - rowK;
            T * cur_nzval = contrib->GetNZval(cur_desc.Offset);
            T * fact_nzval = cur_snodeInd->GetNZval(fact_desc.Offset);

            //TODO CHECK THIS
            blas::Gemv( 'N', nrhs, cur_nrows, T(-1.0), &cur_nzval[(rowK)*ldsol+nrhsOffset],ldsol,
                &fact_nzval[rowK*ldfact+kk],ldfact, T(1.0), &updated_nzval[(kk)*ldsol+nrhsOffset], 1);
          }
        }
      }

    }


  template <typename T, class Allocator> 
    inline void SuperNodeInd<T,Allocator>::forward_update_contrib( SuperNode<T> * cur_snode, Int nrhsOffset, Int pnrhs){
      SuperNodeInd<T> * cur_snodeInd = static_cast<SuperNodeInd<T>*>(cur_snode);
      SuperNodeInd<T,Allocator> * contrib = this;
      Int snodeSize = cur_snodeInd->Size();

      Int nrhs = pnrhs;
      if(pnrhs==-1){
        nrhs = contrib->Size();
        assert(nrhsOffset==0);
      }
      else{
        nrhs = std::min(pnrhs,contrib->Size() - nrhsOffset);
      }

      Int ldsol = contrib->Size();
      Int ldfact = cur_snodeInd->Size();




      if(cur_snodeInd->NZBlockCnt()>=1){
        NZBlockDesc & diag_desc = contrib->GetNZBlockDesc(0);
        Int diag_nrows = contrib->NRows(0);
        T * diag_nzval = contrib->GetNZval(diag_desc.Offset);
        T * chol_diag = cur_snodeInd->GetDiag();

        for(Int blkidx = 0; blkidx < cur_snodeInd->NZBlockCnt() ; ++blkidx){
          NZBlockDesc & cur_desc = contrib->GetNZBlockDesc(blkidx);
          NZBlockDesc & chol_desc = cur_snodeInd->GetNZBlockDesc(blkidx);

          Int cur_nrows = contrib->NRows(blkidx);
          Int chol_nrows = cur_snodeInd->NRows(blkidx);

          T * cur_nzval = contrib->GetNZval(cur_desc.Offset);
          T * chol_nzval = cur_snodeInd->GetNZval(chol_desc.Offset);

          //compute my contribution
          //Handle the diagonal block
          if(blkidx==0){
            //TODO That's where we can use the selective inversion
            //if we are processing the "pivot" block
            for(Int kk = 0; kk<snodeSize; ++kk){
              //then compute the rank one update
              blas::Geru(nrhs,cur_nrows-kk-1, T(-1.0),  &diag_nzval[kk*ldsol+nrhsOffset], 1,&chol_nzval[(kk + 1)*ldfact+kk], ldfact, &cur_nzval[(kk + 1)*ldsol+nrhsOffset], ldsol );
            }
          }
          else{
            for(Int kk = 0; kk<snodeSize; ++kk){
              //compute the rank one update
              blas::Geru(nrhs,cur_nrows, T(-1.0), &diag_nzval[kk*ldsol+nrhsOffset], 1, &chol_nzval[kk], ldfact, &cur_nzval[nrhsOffset], ldsol );
            }
          }
        }

        for(Int kk = 0; kk<snodeSize; ++kk){
          //scale the kk-th row of the solution
          lapack::Scal(nrhs, T(1.0)/chol_diag[kk], &diag_nzval[kk*ldsol+nrhsOffset], 1);
        }



      }
    }


  template <typename T, class Allocator> 
    inline void SuperNodeInd<T,Allocator>::forward_update_contrib( T * RHS, SuperNode<T> * cur_snode, std::vector<Int> & perm){
      SuperNodeInd<T> * cur_snodeInd = static_cast<SuperNodeInd<T>*>(cur_snode);
      SuperNodeInd<T,Allocator> * contrib = this;
      Int n = perm.size();
      Int nrhs = this->Size();
      Int snodeSize = cur_snodeInd->Size();

      if(cur_snodeInd->NZBlockCnt()>=1){
        NZBlockDesc & diag_desc = contrib->GetNZBlockDesc(0);
        Int diag_nrows = contrib->NRows(0);
        T * diag_nzval = contrib->GetNZval(diag_desc.Offset);
        T * chol_diag = cur_snodeInd->GetDiag();

        for(Int blkidx = 0; blkidx < cur_snodeInd->NZBlockCnt() ; ++blkidx){
          NZBlockDesc & cur_desc = contrib->GetNZBlockDesc(blkidx);
          NZBlockDesc & chol_desc = cur_snodeInd->GetNZBlockDesc(blkidx);

          Int cur_nrows = contrib->NRows(blkidx);
          Int chol_nrows = cur_snodeInd->NRows(blkidx);

          T * cur_nzval = contrib->GetNZval(cur_desc.Offset);
          T * chol_nzval = cur_snodeInd->GetNZval(chol_desc.Offset);

          //compute my contribution
          //Handle the diagonal block
          if(blkidx==0){
            //TODO That's where we can use the selective inversion
            //if we are processing the "pivot" block
            for(Int kk = 0; kk<cur_snodeInd->Size(); ++kk){
              //First, add the RHS into the contribution
              Int srcRow = perm[diag_desc.GIndex+kk-1];
              for(Int j = 0; j<nrhs;++j){
                diag_nzval[kk*nrhs+j] += RHS[srcRow-1 + j*n];
              }

              //then compute the rank one update
              blas::Geru(nrhs,cur_nrows-kk-1, T(-1.0),  &diag_nzval[kk*nrhs], 1,&chol_nzval[(kk + 1)*snodeSize+kk], snodeSize, &cur_nzval[(kk + 1)*nrhs], nrhs );
            }
          }
          else{
            for(Int kk = 0; kk<cur_snodeInd->Size(); ++kk){
              //compute the rank one update
              blas::Geru(nrhs,cur_nrows, T(-1.0), &diag_nzval[kk*nrhs], 1, &chol_nzval[kk], snodeSize, &cur_nzval[0], nrhs );
            }

          }
        }


        for(Int kk = 0; kk<cur_snodeInd->Size(); ++kk){
          //scale the kk-th row of the solution
          lapack::Scal(nrhs, T(1.0)/chol_diag[kk], &diag_nzval[kk*nrhs], 1);
        }


      }
    }









  template <typename T, class Allocator> inline std::ostream& operator<<( std::ostream& os,  SuperNodeInd<T, Allocator>& snode){
    os<<"ooooooooooo   Supernode "<<snode.Id()<<" oooooooooooo"<<std::endl;
    os<<"     size = "<<snode.Size()<<std::endl;
    os<<"     fc   = "<<snode.FirstCol()<<std::endl;
    os<<"     lc   = "<<snode.LastCol()<<std::endl;
    os<<"     n    = "<<snode.N()<<std::endl;
    os<<"--- Diagonal ---"<<std::endl;
    for(Int j = 0; j<snode.Size();++j){
      os<<snode.GetDiag()[j]<<" ";
    }
    os<<std::endl;
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




  template <typename T,class Allocator> inline size_t Deserialize(char * buffer, size_t size, SuperNodeInd<T,Allocator> & snode){
    snode.Init(&buffer[0],size);
    snode.InitIdxToBlk();
    return 0;//(last_ptr - buffer);
  }



} // namespace SYMPACK
#endif // _SUPERNODE_IND_IMPL_HPP_
