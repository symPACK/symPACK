#ifndef _SUPERNODE_IMPL_HPP_
#define _SUPERNODE_IMPL_HPP_

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
/*           |      Meta       |        */
/*           |_________________|        */
/*           |                 |        */
/*           |   block desc    |        */
/*           |                 |        */
/*           |_________________|        */
/*                                      */
/****************************************/

namespace symPACK{

  //SuperNode implementation
  template<typename T, class Allocator>
    SuperNode<T,Allocator>::SuperNode() : meta_(nullptr), blocks_(nullptr), nzval_(nullptr), storage_size_(0),
    loc_storage_container_(nullptr), storage_container_(nullptr){
#ifndef ITREE
      globalToLocal_=nullptr;
#else
      idxToBlk_=nullptr;
#endif
    }

  //this allocate the supernode with the number of rows and number of nzblocks known a priori
  template<typename T, class Allocator>
    SuperNode<T,Allocator>::SuperNode(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int ai_num_rows, Int aiN, Int aiNZBlkCnt, Int panel):SuperNode<T,Allocator>() {
      assert(ai_num_rows>=0);

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
        storage_size_ = sizeof(T)*((numPanels-1)*ai_num_rows*panel - 
            - (numPanels-1)*(numPanels-2)*panel/2 
            + (ai_num_rows-(numPanels-1)*panel)*lastPanel);
        storage_size_ += num_blocks*sizeof(NZBlockDesc) + sizeof(SuperNodeDesc);
      }
      else {
        storage_size_ = sizeof(T)*size*ai_num_rows 
          + num_blocks*sizeof(NZBlockDesc) 
          + sizeof(SuperNodeDesc);
      }

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

      nzval_ = (T*)&loc_storage_container_[0];
      meta_ = (SuperNodeDesc*)(nzval_+size*ai_num_rows);
      char * last = loc_storage_container_+storage_size_-1 - (sizeof(NZBlockDesc) -1);
      blocks_ = (NZBlockDesc*) last;

      meta_->panel_sz_ = panel;
      meta_->iId_ = aiId;
      meta_->iFirstRow_ = aiFr;
      meta_->iFirstCol_ = aiFc;
      meta_->iLastCol_ = aiLc;
      meta_->iN_=aiN;
      meta_->iSize_ = size;
      meta_->nzval_cnt_ = 0;
      meta_->blocks_cnt_ = 0;
      meta_->b_own_storage_ = true;

#ifndef ITREE
      globalToLocal_ = new std::vector<Int>(aiN+1,-1);
#else
      idxToBlk_ = CreateITree();
#endif
    } 

  template<typename T, class Allocator>
    SuperNode<T,Allocator>::SuperNode(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices, Int panel):SuperNode<T,Allocator>() {
      Init(aiId,aiFr,aiFc,aiLc,aiN,rowIndices,panel);
    } 

  template<typename T, class Allocator>
    SuperNode<T,Allocator>::SuperNode(char * storage_ptr,size_t storage_size, Int GIndex ):SuperNode<T,Allocator>() {
      Init(storage_ptr,storage_size,GIndex);
    }

  template<typename T, class Allocator>
    SuperNode<T,Allocator>::~SuperNode(){
      if(loc_storage_container_!=nullptr && meta_!=nullptr){
        if(meta_->b_own_storage_){
          Allocator::deallocate(loc_storage_container_);
        }
#ifndef ITREE
        if(globalToLocal_!=nullptr){
          delete globalToLocal_;
        }
#else
        if(idxToBlk_!=nullptr){
          delete idxToBlk_;
        }
#endif
      }
    }

  //CHECKED ON 11-18-2014
  //
  template<typename T, class Allocator>
    void SuperNode<T,Allocator>::Init(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices, Int panel) {
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
      //TODO revise this to handle panel size
      if ( panel > 0 ) {
        Int numPanels = std::ceil((double)size/(double)panel);
        Int lastPanel = std::max(size - (numPanels-1)*panel,0);
        storage_size_ = sizeof(T)*((numPanels-1)*numRows*panel - 
            - (numPanels-1)*(numPanels-2)*panel/2 
            + (numRows-(numPanels-1)*panel)*lastPanel);
        storage_size_ += num_blocks*sizeof(NZBlockDesc) + sizeof(SuperNodeDesc);
      }
      else {
        storage_size_ = sizeof(T)*size*numRows + num_blocks*sizeof(NZBlockDesc) + sizeof(SuperNodeDesc);
      }

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



      nzval_ = (T*)&loc_storage_container_[0];
      meta_ = (SuperNodeDesc*)(nzval_+size*numRows);
      char * last = loc_storage_container_+storage_size_-1 - (sizeof(NZBlockDesc) -1);
      blocks_ = (NZBlockDesc*) last;

      meta_->panel_sz_ = panel;
      meta_->iId_ = aiId;
      meta_->iFirstRow_ = aiFr;
      meta_->iFirstCol_ = aiFc;
      meta_->iLastCol_ = aiLc;
      meta_->iN_=aiN;
      meta_->iSize_ = size;
      meta_->nzval_cnt_ = 0;
      meta_->blocks_cnt_ = 0;
      meta_->b_own_storage_ = true;

#ifndef ITREE
      globalToLocal_ = new std::vector<Int>(aiN+1,-1);
#else
      idxToBlk_ = CreateITree();
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
    void SuperNode<T,Allocator>::Init(char * storage_ptr,size_t storage_size, Int GIndex ) {
      //loop through the block descriptors
      char * last = (char*)(storage_ptr+storage_size-1) - (sizeof(NZBlockDesc) -1);

      blocks_ = (NZBlockDesc*) last;

      Int blkCnt = 0;
      NZBlockDesc * curBlockPtr = nullptr;
      do{
        curBlockPtr = &GetNZBlockDesc(blkCnt);
        ++blkCnt;
      }
      while(!curBlockPtr->Last);


      meta_= (SuperNodeDesc*)(blocks_ - blkCnt +1) - 1;
      nzval_ = (T*) storage_ptr;

      //we now need to update the meta data
      meta_->b_own_storage_ = false;
      meta_->blocks_cnt_ = blkCnt;
      meta_->nzval_cnt_ = (T*)meta_ - nzval_;

      if(GIndex==-1){
        GIndex = GetNZBlockDesc(0).GIndex;
      }

      meta_->iFirstRow_ = GIndex;

      Int loc_fr = GIndex - GetNZBlockDesc(0).GIndex;
      Int offset = GetNZBlockDesc(0).Offset + loc_fr*Size();
      //restore 0-based offsets and compute global_to_local indices
      GetNZBlockDesc(0).Offset = 0;
      for(Int blkidx=1; blkidx<NZBlockCnt();++blkidx){
        curBlockPtr = &GetNZBlockDesc(blkidx);
        curBlockPtr->Offset -= offset;
      }



      blocks_->GIndex = GIndex;


#ifndef ITREE
      globalToLocal_ = new std::vector<Int>(meta_->iN_+1,-1);
#else
      idxToBlk_ = CreateITree();
#endif

    }

  template<typename T, class Allocator>
    inline void SuperNode<T,Allocator>::AddNZBlock(Int aiNRows, Int aiGIndex){

      //Resize the container if I own the storage
      if(meta_->b_own_storage_){
        logfileptr->OFS()<<"Adding\n";
        scope_timer(a,RESIZE_SUPERNODE);

        Int cur_fr = aiGIndex;
        Int cur_lr = cur_fr + aiNRows -1;
        Int cur_nzval_cnt = aiNRows*meta_->iSize_;
        if ( cur_fr == FirstCol() && meta_->panel_sz_ > 0 ) {
          Int panel = meta_->panel_sz_;
          Int numPanels = std::ceil((double)meta_->iSize_/(double)panel);
          Int lastPanel = std::max(meta_->iSize_ - (numPanels-1)*panel,0);
          cur_nzval_cnt = (numPanels-1)*aiNRows*panel - 
            - (numPanels-1)*(numPanels-2)*panel/2 
            + (aiNRows-(numPanels-1)*panel)*lastPanel;
        }

#ifndef ITREE
        std::fill(&(*globalToLocal_)[cur_fr],&(*globalToLocal_)[cur_lr]+1,meta_->blocks_cnt_); 
#else
        ITree<Int>::Interval<Int> cur_interv;
        cur_interv.low = cur_fr;
        cur_interv.high = cur_lr;
        cur_interv.data = meta_->blocks_cnt_;
        idxToBlk_->Insert(cur_interv);
#endif

        bool ownDiagonal = OwnDiagonal();
        //if there is no more room for either nzval or blocks, extend
        Int block_space = (Int)(blocks_+1 - (NZBlockDesc*)(meta_ +1)) - meta_->blocks_cnt_;
        size_t nzval_space = ((size_t)((char*)meta_ - (char*)nzval_) - meta_->nzval_cnt_*sizeof(T) );

        if(block_space==0 || nzval_space<cur_nzval_cnt*sizeof(T)){
          //need to resize storage space. this is expensive !
          Int size = Size();
          size_t extra_nzvals_bytes = std::max((size_t)0,(cur_nzval_cnt*sizeof(T) - nzval_space));
          Int extra_blocks = std::max((Int)0,(Int)1 - block_space);
          size_t new_size = storage_size_ + extra_nzvals_bytes + extra_blocks*sizeof(NZBlockDesc);

          size_t offset_meta = (char*)meta_ - (char*)nzval_;
          size_t offset_block = (char*)blocks_ - (char*)nzval_;

          char * locTmpPtr = Allocator::allocate(new_size);
          bassert(locTmpPtr!=nullptr);
          std::copy(loc_storage_container_,loc_storage_container_+storage_size_,locTmpPtr);

          Allocator::deallocate(loc_storage_container_);
          loc_storage_container_=locTmpPtr;
          storage_size_=new_size;

          nzval_=(T*)&loc_storage_container_[0];
          //move the block descriptors if required
          char * cur_blocks_ptr = (char*)&loc_storage_container_[0] + offset_block;
          //move the meta data if required
          char * cur_meta_ptr = (char*)&loc_storage_container_[0] + offset_meta;
          //move the updated rows if required

          meta_ = (SuperNodeDesc*) cur_meta_ptr;
          //we need to move everything, starting from the blocks, then meta
          //blocks need to be moved by extra_nzvals_bytes + extra_blocks*(sizeof(NZBlockDesc) + sizeof(Idx))
          char * new_blocks_ptr = cur_blocks_ptr + extra_nzvals_bytes +extra_blocks*sizeof(NZBlockDesc);
          if(ownDiagonal){
            new_blocks_ptr += extra_blocks*sizeof(Idx);
          }
          std::copy_backward(cur_blocks_ptr - (meta_->blocks_cnt_-1)*sizeof(NZBlockDesc),cur_blocks_ptr+sizeof(NZBlockDesc),new_blocks_ptr+sizeof(NZBlockDesc));

          //now move the meta data by extra_nzvals_bytes + extra_blocks*sizeof(Idx)
          char * new_meta_ptr = cur_meta_ptr + extra_nzvals_bytes;
          if(ownDiagonal){
            new_meta_ptr += extra_blocks*sizeof(Idx);
          }
          std::copy(cur_meta_ptr,cur_meta_ptr + sizeof(SuperNodeDesc),new_meta_ptr);

          //now move the updated rows by extra_nzvals_bytes

          //update pointers
          meta_ = (SuperNodeDesc*) new_meta_ptr;
          blocks_ = (NZBlockDesc*) new_blocks_ptr;
        }

        GetNZBlockDesc(meta_->blocks_cnt_).GIndex = aiGIndex;
        GetNZBlockDesc(meta_->blocks_cnt_).Offset = meta_->nzval_cnt_;
        GetNZBlockDesc(meta_->blocks_cnt_).Last = true;
        if(meta_->blocks_cnt_>0){
          GetNZBlockDesc(meta_->blocks_cnt_-1).Last = false;
        }

        meta_->blocks_cnt_++;

        //fill the new block with zeros
        std::fill(nzval_+meta_->nzval_cnt_,nzval_+meta_->nzval_cnt_+cur_nzval_cnt,ZERO<T>());

        //update nzval count
        meta_->nzval_cnt_+=cur_nzval_cnt;

        //Handle device buffer

      }
    }

  template<typename T, class Allocator>
    inline Int SuperNode<T,Allocator>::FindBlockIdx(Int aiGIndex){
      scope_timer(a,FindBlockIdx);


      Int rval = -1;

#ifndef ITREE
      rval = globalToLocal_->at(aiGIndex);
#else

#ifdef _LAZY_INIT_
      if(!ITreeInitialized()){
        InitIdxToBlk();
      }
#endif


      ITree<Int>::Interval<Int> * res = idxToBlk_->IntervalSearch(aiGIndex,aiGIndex);
      if (res != nullptr){
        rval = res->data;
      }
#endif
      return rval;
    }

  template<typename T, class Allocator>
    inline Int SuperNode<T,Allocator>::FindBlockIdx(Int aiGIndex,Int & closestR, Int & closestL){
      scope_timer(a,FindBlockIdxRL);
      Int rval = -1;
      ITree<Int>::Interval<Int> * L = nullptr;
      ITree<Int>::Interval<Int> * R = nullptr;
      ITree<Int>::Interval<Int> * res = idxToBlk_->IntervalSearch(aiGIndex,aiGIndex,R,L);
      if(R!=nullptr){
        closestR = R->low;
      }

      if(L!=nullptr){
        closestL = L->high;
      }

      if (res != nullptr){
        rval = res->data;
      }
      return rval;
    }

  template<typename T, class Allocator>
    inline Int SuperNode<T,Allocator>::FindBlockIdx(Int fr, Int lr, ITree<Int>::Interval<Int> & overlap){
      scope_timer(a,FindBlockIdxOverlap);

      Int rval = -1;

      ITree<Int>::Interval<Int> * res = idxToBlk_->IntervalSearch(fr,lr);
      if (res != nullptr){
        overlap = *res;
        rval = res->data;
      }
      return rval;
    }

  template<typename T, class Allocator>
    inline void SuperNode<T,Allocator>::DumpITree(){
#ifdef ITREE
      logfileptr->OFS()<<"Number of blocks: "<<meta_->blocks_cnt_<<std::endl;
      logfileptr->OFS()<<"log2(Number of blocks): "<<log2(meta_->blocks_cnt_)<<std::endl;
      idxToBlk_->Dump();
#endif
    }

  template<typename T, class Allocator>
    inline Int SuperNode<T,Allocator>::Shrink(){
      if(meta_->b_own_storage_){
        //TODO make sure that we do not have any extra space anywhere.
        logfileptr->OFS()<<"Shrinking\n";

        bool ownDiagonal = OwnDiagonal();
        //if there is too much room for either nzval or blocks, contract
        Int block_space = (Int)(blocks_+1 - (NZBlockDesc*)(meta_ +1)) - meta_->blocks_cnt_;
        size_t nzval_space = ((size_t)((char*)meta_ - (char*)nzval_) - meta_->nzval_cnt_*sizeof(T) );

        if(block_space >0 || nzval_space >0){

          size_t new_size = storage_size_ - nzval_space - block_space*sizeof(NZBlockDesc);

          size_t offset_meta = meta_->nzval_cnt_*sizeof(T);
          size_t offset_block = offset_meta +sizeof(SuperNodeDesc);

          char * locTmpPtr = Allocator::allocate(new_size);
          bassert(locTmpPtr!=nullptr);

          //copy nzvals
          std::copy(nzval_,nzval_+meta_->nzval_cnt_,(T*)locTmpPtr);

          //copy meta
          std::copy(meta_,meta_+1,(SuperNodeDesc*)(locTmpPtr+offset_meta));
          //copy blocks
          std::copy(blocks_+1-meta_->blocks_cnt_,blocks_+1,(NZBlockDesc*)(locTmpPtr+offset_block));

          Allocator::deallocate(loc_storage_container_);
          loc_storage_container_=locTmpPtr;
          storage_size_=new_size;

          nzval_=(T*)&loc_storage_container_[0];
          //move the block descriptors if required
          char * new_blocks_ptr = loc_storage_container_+storage_size_-1 - (sizeof(NZBlockDesc) -1);
          //move the meta data if required
          char * new_meta_ptr = loc_storage_container_ + offset_meta;

          //update pointers
          meta_ = (SuperNodeDesc*) new_meta_ptr;
          blocks_ = (NZBlockDesc*) new_blocks_ptr;

        }
      }

      return StorageSize();
    }

  template<typename T, class Allocator>
    inline void SuperNode<T,Allocator>::FindUpdatedFirstCol(SuperNode<T,Allocator> * src_snode, Int pivot_fr, Int & tgt_fc, Int & first_pivot_idx){
      //find the first row updated by src_snode
      scope_timer(a,UPDATE_SNODE_FIND_INDEX);


#ifdef _LINEAR_SEARCH_FCLC_
      if(src_snode->ITreeInitialized()){
#endif 
        tgt_fc = pivot_fr;
        first_pivot_idx = -1;
        if(tgt_fc == I_ZERO ){
          tgt_fc = FirstCol();
          //find the pivot idx

          Int tgt_lc = LastCol();
#ifndef _BINARY_BLOCK_SEARCH_
          do {first_pivot_idx = src_snode->FindBlockIdx(tgt_fc); tgt_fc++;}
          while(first_pivot_idx<0 && tgt_fc<=tgt_lc);
          tgt_fc--;
#else
          Int closestR = -1;
          Int closestL = -1;
          do{
            first_pivot_idx = src_snode->FindBlockIdx(tgt_fc,closestR,closestL);
            if(closestR!=-1){
              logfileptr->OFS()<<"ClosestR of "<<tgt_fc<<" is "<<closestR<<std::endl;
              tgt_fc = closestR;
            }
          }
          while(first_pivot_idx<0);
#endif

        }
        else{
          first_pivot_idx = src_snode->FindBlockIdx(tgt_fc);
        }
#ifdef _LINEAR_SEARCH_FCLC_
      }
      else{
        tgt_fc = pivot_fr;
        first_pivot_idx = -1;
        if(tgt_fc == I_ZERO ){
          tgt_fc = FirstCol();

          //this assumes that blocks are sorted...
          for(Int blkidx = 0; blkidx<src_snode->NZBlockCnt();++blkidx){
            NZBlockDesc & cur_block = src_snode->GetNZBlockDesc(blkidx);
            Int cur_fr = cur_block.GIndex;
            Int cur_lr = cur_block.GIndex + src_snode->NRows(blkidx)-1;

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
          for(Int blkidx = 0; blkidx<src_snode->NZBlockCnt();++blkidx){
            NZBlockDesc & cur_block = src_snode->GetNZBlockDesc(blkidx);
            Int cur_fr = cur_block.GIndex;
            Int cur_lr = cur_block.GIndex + src_snode->NRows(blkidx)-1;
            if(cur_fr<= tgt_fc && tgt_fc<=cur_lr){
              first_pivot_idx = blkidx;
              break;
            }
          }
        }
      }
#endif
      if(first_pivot_idx<0){
        logfileptr->OFS()<<"LOCK 1: first_pivot_idx<0"<<std::endl;
        gdb_lock();
      }
      assert(first_pivot_idx>=0);
    }

  template<typename T, class Allocator>
    inline void SuperNode<T,Allocator>::FindUpdatedLastCol(SuperNode<T,Allocator> * src_snode, Int tgt_fc, Int first_pivot_idx , Int & tgt_lc,  Int & last_pivot_idx){
      scope_timer(a,UPDATE_SNODE_FIND_INDEX);
#ifdef _LINEAR_SEARCH_FCLC_
      if(src_snode->ITreeInitialized()){
#endif 
        //find the last row updated by src_snode
        tgt_lc = LastCol();
        last_pivot_idx = -1;
        //find the pivot idx
#ifndef _BINARY_BLOCK_SEARCH_
        do {last_pivot_idx = src_snode->FindBlockIdx(tgt_lc); tgt_lc--;}
        while(last_pivot_idx<0 && tgt_lc>=tgt_fc);
        tgt_lc++;
#else
        for(tgt_lc;tgt_lc>=tgt_fc;tgt_lc--){
          last_pivot_idx = src_snode->FindBlockIdx(tgt_lc);
          if(last_pivot_idx>=0){
            break;
          }
        }
#endif

#ifdef _LINEAR_SEARCH_FCLC_
      }
      else{
        tgt_lc = tgt_fc;
        for(Int blkidx = first_pivot_idx; blkidx<src_snode->NZBlockCnt();++blkidx){
          NZBlockDesc & cur_block = src_snode->GetNZBlockDesc(blkidx);
          Int cur_fr = cur_block.GIndex;
          Int cur_lr = cur_block.GIndex + src_snode->NRows(blkidx)-1;


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
    }

  template<typename T, class Allocator>
    inline Int SuperNode<T,Allocator>::Merge(SuperNode<T,Allocator> * src_snode, SnodeUpdate &update){
      scope_timer(a,MERGE_SNODE);

      bassert(meta_->b_own_storage_);

      Int src_snode_size = src_snode->Size();
      Int tgt_snode_size = Size();

      Int & pivot_idx = update.blkidx;
      Int & pivot_fr = update.src_first_row;

      Int tgt_fc;
      Int first_pivot_idx;
      FindUpdatedFirstCol(src_snode, 0, tgt_fc, first_pivot_idx);
      NZBlockDesc & first_pivot_desc = src_snode->GetNZBlockDesc(first_pivot_idx);

      //parse src_snode
      ITree<Int>::Interval<Int> overlap;
      ITree<Int>::Interval<Int> curInter;
      ITree<Int>::Interval<Int> newInter;
      std::queue< ITree<Int>::Interval<Int> > toInsert;
      for(Int blkidx = first_pivot_idx; blkidx<src_snode->NZBlockCnt(); ++blkidx){
        NZBlockDesc & blk_desc = src_snode->GetNZBlockDesc(blkidx);
        Int fr = std::max(FirstCol(),blk_desc.GIndex);
        Int lr = blk_desc.GIndex + src_snode->NRows(blkidx) -1;
        if(lr<fr){gdb_lock();}
        assert(lr>=fr);
        curInter.low = fr;
        curInter.high = lr;
        toInsert.push(curInter);

        while(!toInsert.empty()){
          curInter = toInsert.front();
          toInsert.pop();
          if(FindBlockIdx(curInter.low,curInter.high,overlap)==-1){
            //Add the full block
            AddNZBlock( curInter.high - curInter.low + 1, curInter.low);
          }
          else{

            //check the overlap
            //fr is curInter.low and lr is curInter.high
            //                l-----overlap------h
            //            l---------block-------------h
            //        l--block----h
            //                              l-----block----h
            //we have two segments to look for : [overlap.high+1 - lr] and [fr - overlap.low -1]         
            if(overlap.low>curInter.low){
              newInter.low = curInter.low;
              newInter.high = overlap.low-1;
              toInsert.push(newInter);
            }

            if(overlap.high < curInter.high){
              newInter.low = overlap.high+1;
              newInter.high = curInter.high;
              toInsert.push(newInter);
            }

          }
        }
      }

      return 0;
    }

  template<typename T, class Allocator>
    inline Int SuperNode<T,Allocator>::Aggregate(SuperNode<T,Allocator> * src_snode){
      logfileptr->OFS()<<"Aggregation normal\n";
      scope_timer(a,AGGREGATE_SNODE);

#if defined(_NO_COMPUTATION_)
      return 0;
#endif

      Int  pivot_idx = 0;
      Int  pivot_fr = 0;

      Int src_snode_size = src_snode->Size();
      Int tgt_snode_size = Size();

      //parse src_snode and add everything

      Int first_pivot_idx = 0 ;
      Int tgt_fc = FirstCol();

      assert(src_snode->NRowsBelowBlock(0) == NRowsBelowBlock(0) && src_snode->Size() == Size());
      if(src_snode->NRowsBelowBlock(0) == NRowsBelowBlock(0) && src_snode->Size() == Size()){
        T * src = src_snode->nzval_;
        T * tgt = nzval_;
#pragma unroll
        for(Int i = 0; i< meta_->nzval_cnt_ ;i++){ tgt[i] += src[i]; }
      }
      else{
        for(Int blkidx = first_pivot_idx; blkidx<src_snode->NZBlockCnt(); ++blkidx){
          NZBlockDesc & blk_desc = src_snode->GetNZBlockDesc(blkidx);
          Int nrows = src_snode->NRows(blkidx);
          for(Int rowidx = 0; rowidx<nrows; ++rowidx){
            Int row = blk_desc.GIndex + rowidx;

            if(row>=tgt_fc){
              Int src_offset = blk_desc.Offset + (row - blk_desc.GIndex)*src_snode_size;

              Int tgt_blkidx = FindBlockIdx(row);
              NZBlockDesc & tgt_desc = GetNZBlockDesc(tgt_blkidx);
              Int tgt_offset = tgt_desc.Offset + (row - tgt_desc.GIndex)*tgt_snode_size;

              T * src = src_snode->GetNZval(src_offset);
              T * tgt = GetNZval(tgt_offset);

#pragma unroll
              for(Int i = 0; i< tgt_snode_size;i+=1){ tgt[i] += src[i]; }

            }
          }
        }
      }

      return 0;
    }


  //CHECKED ON 11-18-2014
  template<typename T, class Allocator>
    inline Int SuperNode<T,Allocator>::UpdateAggregate(SuperNode<T,Allocator> * src_snode, SnodeUpdate &update, TempUpdateBuffers<T> & tmpBuffers, Int iTarget, Int iam){

      scope_timer(a,UPDATE_AGGREGATE_SNODE);
#if defined(_NO_COMPUTATION_)
      return 0;
#endif


      if(iTarget != iam){
        Merge(src_snode, update);

        Int & pivot_idx = update.blkidx;
        Int & pivot_fr = update.src_first_row;

        Int src_snode_size = src_snode->Size();
        Int tgt_snode_size = Size();

        Int tgt_fc,tgt_lc;
        Int first_pivot_idx,last_pivot_idx;
        FindUpdatedFirstCol(src_snode, update.src_first_row, tgt_fc, first_pivot_idx);
        FindUpdatedLastCol(src_snode, tgt_fc, first_pivot_idx, tgt_lc, last_pivot_idx);

        NZBlockDesc & first_pivot_desc = src_snode->GetNZBlockDesc(first_pivot_idx);
        NZBlockDesc & last_pivot_desc = src_snode->GetNZBlockDesc(last_pivot_idx);

        //determine the first column that will be updated in the target supernode
        Int tgt_local_fc =  tgt_fc - FirstCol();
        Int tgt_local_lc =  tgt_lc - FirstCol();

        Int tgt_nrows = NRowsBelowBlock(0);
        Int src_nrows = src_snode->NRowsBelowBlock(first_pivot_idx)
          - (tgt_fc - first_pivot_desc.GIndex);
        Int src_lr = tgt_fc+src_nrows-1;
        src_nrows = src_lr - tgt_fc + 1;

        Int tgt_width = src_nrows - src_snode->NRowsBelowBlock(last_pivot_idx)
          + (tgt_lc - last_pivot_desc.GIndex)+1;

        T * pivot = src_snode->GetNZval(first_pivot_desc.Offset)
          + (tgt_fc-first_pivot_desc.GIndex)*src_snode_size;
        T * tgt = GetNZval(0);

        //Pointer to the output buffer of the GEMM
        T * buf = nullptr;
        T beta = ZERO<T>();

#ifdef SP_THREADS
        tmpBuffers.tmpBuf.resize(tgt_width*src_nrows + src_snode_size*tgt_width);
#endif
        buf = &tmpBuffers.tmpBuf[0];

        //everything is in row-major
        SYMPACK_TIMER_START(UPDATE_SNODE_GEMM);
        blas::Gemm('T','N',tgt_width, src_nrows,src_snode_size,
            T(-1.0),pivot,src_snode_size,
            pivot,src_snode_size,beta,buf,tgt_width);
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
              Int tgt_blk_idx = FindBlockIdx(row);
              assert(tgt_blk_idx>=0);
              NZBlockDesc & cur_tgt_desc = GetNZBlockDesc(tgt_blk_idx);
              Int lr = std::min(cur_src_lr,cur_tgt_desc.GIndex + NRows(tgt_blk_idx)-1);
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
              Int tgt_blk_idx = FindBlockIdx(row);
              assert(tgt_blk_idx>=0);
              NZBlockDesc & cur_tgt_desc = GetNZBlockDesc(tgt_blk_idx);
              Int lr = std::min(cur_src_lr,cur_tgt_desc.GIndex + NRows(tgt_blk_idx)-1);
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
            Int tgt_offset = (tgt_fc - FirstCol());
            for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
              T * A = &buf[rowidx*tgt_width];
              T * B = &tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_offset];
#pragma unroll
              for(Int i = 0; i < tgt_width; ++i){ B[i] += A[i]; }
            }
          }
          else{
            // full sparse case (done right now)
            for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
              for(Int colidx = 0; colidx< tmpBuffers.src_colindx.size();++colidx){
                Int col = tmpBuffers.src_colindx[colidx];
                Int tgt_colidx = col - FirstCol();
                tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_colidx] 
                  += buf[rowidx*tgt_width+colidx]; 
              }
            }
          }
        }
        return 0;

      }
      else{
        return Update(src_snode, update, tmpBuffers);
      }
    }


  //CHECKED ON 11-18-2014
  template<typename T, class Allocator>
    inline Int SuperNode<T,Allocator>::Update(SuperNode<T,Allocator> * src_snode, SnodeUpdate &update, TempUpdateBuffers<T> & tmpBuffers){

      scope_timer(a,UPDATE_SNODE);
#if defined(_NO_COMPUTATION_)
      return 0;
#endif


      Int & pivot_idx = update.blkidx;
      Int & pivot_fr = update.src_first_row;

      Int src_snode_size = src_snode->Size();
      Int tgt_snode_size = Size();

      //find the first row updated by src_snode
      Int tgt_fc,tgt_lc;
      Int first_pivot_idx,last_pivot_idx;
      FindUpdatedFirstCol(src_snode, update.src_first_row, tgt_fc, first_pivot_idx);
      FindUpdatedLastCol(src_snode, tgt_fc, first_pivot_idx, tgt_lc, last_pivot_idx);

      NZBlockDesc & first_pivot_desc = src_snode->GetNZBlockDesc(first_pivot_idx);
      NZBlockDesc & last_pivot_desc = src_snode->GetNZBlockDesc(last_pivot_idx);

      //determine the first column that will be updated in the target supernode
      Int tgt_local_fc =  tgt_fc - FirstCol();
      Int tgt_local_lc =  tgt_lc - FirstCol();

      Int tgt_nrows = NRowsBelowBlock(0);
      Int src_nrows = src_snode->NRowsBelowBlock(first_pivot_idx)
        - (tgt_fc - first_pivot_desc.GIndex);
      Int src_lr = tgt_fc+src_nrows-1;
      src_nrows = src_lr - tgt_fc + 1;

      Int src_belowLast = src_snode->NRowsBelowBlock(last_pivot_idx);
      Int tgt_width = src_nrows - src_belowLast
        + (tgt_lc - last_pivot_desc.GIndex)+1;

      T * pivot = src_snode->GetNZval(first_pivot_desc.Offset)
        + (tgt_fc-first_pivot_desc.GIndex)*src_snode_size;
      bassert(  src_snode->GetNZval(0) + src_snode->NNZ() -  pivot >= 0 );
      T * tgt = GetNZval(0);

      //Pointer to the output buffer of the GEMM
      T * buf = nullptr;
      T beta = ZERO<T>();
#ifdef SP_THREADS
      tmpBuffers.tmpBuf.resize(tgt_width*src_nrows + src_snode_size*tgt_width);
#endif
      //If the target supernode has the same structure,
      //The GEMM is directly done in place
      if(src_nrows == tgt_nrows){
        Int tgt_offset = (tgt_fc - FirstCol());
        buf = &tgt[tgt_offset];
        beta = ONE<T>();
      }
      else{
        //Compute the update in a temporary buffer
        buf = &tmpBuffers.tmpBuf[0];
      }

      //everything is in row-major
      SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_GEMM);
      blas::Gemm('T','N',tgt_width, src_nrows,src_snode_size,
          T(-1.0),pivot,src_snode_size,
          pivot,src_snode_size,beta,buf,tgt_width);
      SYMPACK_TIMER_SPECIAL_STOP(UPDATE_SNODE_GEMM);

      //If the GEMM wasn't done in place we need to aggregate the update
      //This is the assembly phase
      if(src_nrows != tgt_nrows){
#ifdef _DEBUG_
        logfileptr->OFS()<<"tmpBuf is "<<tmpBuffers.tmpBuf<<std::endl;
#endif

        //now add the update to the target supernode
        SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_INDEX_MAP);
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
              Int tgt_blk_idx = FindBlockIdx(row);
              assert(tgt_blk_idx>=0);
              NZBlockDesc & cur_tgt_desc = GetNZBlockDesc(tgt_blk_idx);
              Int lr = std::min(cur_src_lr,cur_tgt_desc.GIndex + NRows(tgt_blk_idx)-1);
              Int tgtOffset = cur_tgt_desc.Offset 
                + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
              SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_ADD_SINGLE);
              for(Int cr = row ;cr<=lr;++cr){
                tgt[tgtOffset + (cr - row)*tgt_snode_size] += buf[rowidx]; 
                rowidx++;
              }
              SYMPACK_TIMER_SPECIAL_STOP(UPDATE_SNODE_ADD_SINGLE);
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

            //The other one MUST reside into a single block in the target
            Int row = cur_src_fr;
            while(row<=cur_src_lr){
              Int tgt_blk_idx = FindBlockIdx(row);
              bassert(tgt_blk_idx>=0);
              NZBlockDesc & cur_tgt_desc = GetNZBlockDesc(tgt_blk_idx);
              Int lr = std::min(cur_src_lr,cur_tgt_desc.GIndex + NRows(tgt_blk_idx)-1);
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
          SYMPACK_TIMER_SPECIAL_STOP(UPDATE_SNODE_INDEX_MAP);

          SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_ADD);
          if(first_pivot_idx==last_pivot_idx){
            // Updating contiguous columns
            Int tgt_offset = (tgt_fc - FirstCol());
            for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
              T * A = &buf[rowidx*tgt_width];
              T * B = &tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_offset];
#pragma unroll
              for(Int i = 0; i < tgt_width; ++i){ B[i] += A[i]; }
            }
          }
          else{
            // full sparse case
            for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
              for(Int colidx = 0; colidx< tmpBuffers.src_colindx.size();++colidx){
                Int col = tmpBuffers.src_colindx[colidx];
                Int tgt_colidx = col - FirstCol();
                tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_colidx] 
                  += buf[rowidx*tgt_width+colidx]; 
              }
            }
          }
          SYMPACK_TIMER_SPECIAL_STOP(UPDATE_SNODE_ADD);
        }
      }
      return 0;
    }


  //CHECKED ON 11-18-2014
  template<typename T, class Allocator>
    inline Int SuperNode<T,Allocator>::Factorize( TempUpdateBuffers<T> & tmpBuffers){
#if defined(_NO_COMPUTATION_)
      return 0;
#endif

      Int snode_size = Size();
      NZBlockDesc & diag_desc = GetNZBlockDesc(0);
      if ( diag_desc.GIndex == this->FirstCol() && this->meta_->panel_sz_ > 0 ) {
        gdb_lock();
        Int panel = meta_->panel_sz_;
        Int numPanels = std::ceil((double)meta_->iSize_/(double)panel);
        Int lastPanel = std::max(meta_->iSize_ - (numPanels-1)*panel,0);
        T * diag_nzval = &GetNZval(diag_desc.Offset)[0];
        for ( Int pan = 0; pan<numPanels-1; pan++ ) {
          Int ld_diag = std::min(meta_->iSize_,(pan+1)*panel);
          Int pansize = std::min(panel,meta_->iSize_ - pan*panel);
          //factor current panel's top diagonal block
          lapack::Potrf( 'U', pansize, diag_nzval, pansize);
          T * nzblk_nzval = diag_nzval + pansize*pansize;
          //update current panel
          for ( Int lpan = pan+1; lpan<numPanels; lpan++ ) {
            Int ld_lpan = std::min(meta_->iSize_,(lpan+1)*panel);
            Int lpansize = std::min(panel,meta_->iSize_ - lpan*panel);
            blas::Trsm('L','U','T','N',pansize, lpansize, T(1.0),  diag_nzval, ld_diag, nzblk_nzval, ld_lpan );
          }
          //TODO update trailing matrix
          diag_nzval += (pan+1)*panel*panel;
        }
      }
      else {
        T * diag_nzval = &GetNZval(diag_desc.Offset)[0];
        lapack::Potrf( 'U', snode_size, diag_nzval, snode_size);
        T * nzblk_nzval = &GetNZval(diag_desc.Offset)[(snode_size)*snode_size];
        blas::Trsm('L','U','T','N',snode_size, NRowsBelowBlock(0)-snode_size, T(1.0),  diag_nzval, snode_size, nzblk_nzval, snode_size);
      }
      return 0;

    }

  template<typename T, class Allocator>
    inline Int SuperNode<T,Allocator>::Factorize(SuperNode<T,Allocator> * diag_snode, TempUpdateBuffers<T> & tmpBuffers){
      abort();
#if defined(_NO_COMPUTATION_)
      return 0;
#endif

      Int snode_size = Size();
      NZBlockDesc & diag_desc = diag_snode->GetNZBlockDesc(0);
      T * diag_nzval = &diag_snode->GetNZval(diag_desc.Offset)[0];
      T * nzblk_nzval = &GetNZval(diag_desc.Offset)[(snode_size)*snode_size];
      blas::Trsm('L','U','T','N',snode_size, NRowsBelowBlock(0)-snode_size, T(1.0),  diag_nzval, snode_size, nzblk_nzval, snode_size);
      return 0;

    }





  template<typename T, class Allocator>
    bool SuperNode<T,Allocator>::FindNextUpdate(SnodeUpdate & nextUpdate, const std::vector<Int> & Xsuper,  const std::vector<Int> & SupMembership, bool isLocal){
      scope_timer(a,FIND_NEXT_UPDATE);
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
        f_ur = std::max(n_ur,cur_desc->GIndex);
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


  template<typename T, class Allocator>
    inline ITree<Int> * SuperNode<T,Allocator>::CreateITree(){
#if defined(_AVL_ITREE_)
      return new AVLITree<Int>();
#elif defined(_DSW_ITREE_)
      return new DSWITree<Int>();
#else
      return new ITree<Int>();
#endif
    } 

  template<typename T, class Allocator>
    inline void SuperNode<T,Allocator>::InitIdxToBlk(){
#ifndef ITREE
      gdb_lock();
      SYMPACK_TIMER_START(ARRAY_INSERT);
      for(Int blkidx=0; blkidx<meta_->blocks_cnt_;++blkidx){
        Int cur_fr = GetNZBlockDesc(blkidx).GIndex;
        Int cur_lr = cur_fr + NRows(blkidx) -1;

        std::fill(&(*globalToLocal_)[cur_fr],&(*globalToLocal_)[cur_lr]+1,blkidx); 
      }
      SYMPACK_TIMER_STOP(ARRAY_INSERT);
#else
      SYMPACK_TIMER_START(BST_INSERT);
      for(Int blkidx=0; blkidx<meta_->blocks_cnt_;++blkidx){
        Int cur_fr = GetNZBlockDesc(blkidx).GIndex;
        Int cur_lr = cur_fr + NRows(blkidx) -1;

        ITree<Int>::Interval<Int> cur_interv;
        cur_interv.low = cur_fr;
        cur_interv.high = cur_lr;
        cur_interv.data = blkidx;

        idxToBlk_->Insert(cur_interv);
      }
      SYMPACK_TIMER_STOP(BST_INSERT);
      idxToBlk_->Rebalance();
#endif
    }


  template<typename T, class Allocator>
    inline void SuperNode<T,Allocator>::Reserve(size_t storage_size){
    }




  template <typename T, class Allocator> 
    inline void SuperNode<T,Allocator>::Serialize(Icomm & buffer, Int first_blkidx, Idx first_row){
      logfileptr->OFS()<<"Serialize called\n";
      NZBlockDesc* nzblk_ptr = &this->GetNZBlockDesc(first_blkidx);
      if(first_row==0){
        first_row = nzblk_ptr->GIndex;
      }
      Idx local_first_row = first_row - nzblk_ptr->GIndex;

      Int nzblk_cnt = this->NZBlockCnt() - first_blkidx;
      T* nzval_ptr = this->GetNZval(nzblk_ptr->Offset) + local_first_row*this->Size();

      size_t size = (char*)(nzblk_ptr+1)-(char*)nzval_ptr;
      buffer.clear();
      buffer.resize(size);
      //copy the whole thing in the buffer
      symPACK::Serialize(buffer,(char*)nzval_ptr,size);
      //now we need to modify the first block data
      char * tail = buffer.front()+size; 
      NZBlockDesc* new_blk_ptr= (NZBlockDesc*)(tail-sizeof(NZBlockDesc));

      Int offset = new_blk_ptr->Offset + local_first_row*this->Size();

      if(offset!=0){
        Int blkCnt = 1;
        if(blkCnt<nzblk_cnt){
          NZBlockDesc * curBlockPtr = nullptr;
          do{
            curBlockPtr = new_blk_ptr - blkCnt;
            curBlockPtr->Offset -= offset;
            ++blkCnt;
          }
          while(!curBlockPtr->Last);
        }
      }
      new_blk_ptr->Offset = 0;
      new_blk_ptr->GIndex=first_row;
    }





  template <typename T, class Allocator> 
    inline void SuperNode<T,Allocator>::forward_update(SuperNode<T,Allocator> * src_contrib, Int iOwner,Int iam, Int nrhsOffset, Int pnrhs){

      SuperNode<T,Allocator> * tgt_contrib = this;

      Int nrhs = pnrhs;
      if(pnrhs==-1){
        nrhs = src_contrib->Size();
        assert(nrhsOffset==0);
      }
      else{
        nrhs = std::min(pnrhs,src_contrib->Size() - nrhsOffset);
      }

      Int ldsol = src_contrib->Size();





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

        Int src_local_fr = std::max(tgt_desc.GIndex - src_desc.GIndex,0);
        Int src_lr = src_desc.GIndex+src_nrows-1;

        Int tgt_local_fr = std::max(src_desc.GIndex - tgt_desc.GIndex,0);
        Int tgt_lr = tgt_desc.GIndex+tgt_nrows-1;
        Int tgt_local_lr = std::min(src_lr,tgt_lr) - tgt_desc.GIndex;

        T * src = &src_contrib->GetNZval(src_desc.Offset)[src_local_fr*ldsol+nrhsOffset];
        T * tgt = &tgt_contrib->GetNZval(tgt_desc.Offset)[tgt_local_fr*ldsol+nrhsOffset];

        if(nrhs!=ldsol){
          Int nrows = tgt_local_lr - tgt_local_fr +1;
          for(Int row=0;row<nrows;row++){
            for(Int col=0; col<nrhs;col++){
              tgt[row*ldsol+col] += src[row*ldsol+col];
            }
          }
        }
        else{
          Int nrows = tgt_local_lr - tgt_local_fr +1;
          for(Int row=0;row<nrows;row++){
            for(Int col=0; col<nrhs;col++){
              tgt[row*ldsol+col] += src[row*ldsol+col];
            }
          }
        }

        if(src_lr>tgt_lr){
          //the src block hasn't been completely used and is
          // updating some lines in the nz block just below the diagonal block
          //this case happens only locally for the diagonal block
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




  template <typename T, class Allocator> 
    inline void SuperNode<T,Allocator>::back_update(SuperNode<T,Allocator> * src_contrib, Int nrhsOffset, Int pnrhs){
      SuperNode<T,Allocator> * tgt_contrib = this;
      Int nrhs = pnrhs;
      if(pnrhs==-1){
        nrhs = tgt_contrib->Size();
        assert(nrhsOffset==0);
      }
      else{
        nrhs = std::min(pnrhs,tgt_contrib->Size() - nrhsOffset);
      }
      Int ldsol = tgt_contrib->Size();



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

          Int src_local_fr = std::max(tgt_desc.GIndex - src_desc.GIndex,(Int)0);
          Int tgt_local_fr = std::max(src_desc.GIndex - tgt_desc.GIndex,(Int)0);
          Int tgt_local_lr = std::min(src_lr,tgt_lr) - tgt_desc.GIndex;

          T * src = &src_contrib->GetNZval(src_desc.Offset)[src_local_fr*ldsol+nrhsOffset];
          T * tgt = &tgt_contrib->GetNZval(tgt_desc.Offset)[tgt_local_fr*ldsol+nrhsOffset];

          if(nrhs!=ldsol){
            Int nrows = tgt_local_lr - tgt_local_fr +1;
            for(Int row=0;row<nrows;row++){
              std::copy(&src[row*ldsol],&src[row*ldsol]+nrhs,&tgt[row*ldsol]);
            }
          }
          else{
            std::copy(src,src+(tgt_local_lr - tgt_local_fr +1)*nrhs,tgt);
          }

          //do other block
          if(tgt_lr>src_lr){
            src_nzblk_idx++;
          }

        } while(tgt_lr>src_lr);
      }
    }


  template <typename T, class Allocator> 
    inline void SuperNode<T,Allocator>::back_update_contrib(SuperNode<T> * cur_snode, Int nrhsOffset, Int pnrhs){
      SuperNode<T,Allocator> * contrib = this;
      Int nrhs = pnrhs;
      nrhs = contrib->Size();

      Int ldsol = contrib->Size();
      Int ldfact = cur_snode->Size();

      NZBlockDesc & diag_desc = cur_snode->GetNZBlockDesc(0);
      NZBlockDesc & tgt_desc = contrib->GetNZBlockDesc(0);

      T* diag_nzval = cur_snode->GetNZval(diag_desc.Offset);
      T* tgt_nzval = contrib->GetNZval(tgt_desc.Offset);

      for(Int j = 0; j<nrhs;++j){
        for(Int ii = ldfact-1; ii>=0; --ii){
          T temp = tgt_nzval[ii*ldsol+j];

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
                  temp += -(chol_nzval[kk*ldfact+ii])*cur_nzval[src_row*ldsol+j];
                }
              }
            }
          }

          //if nounit
          temp = temp / (diag_nzval[ii*ldfact+ii]);
          tgt_nzval[ii*ldsol+j] = temp;
        }
      }
    }

  template <typename T, class Allocator> 
    inline void SuperNode<T,Allocator>::forward_update_contrib(SuperNode<T> * cur_snode, Int nrhsOffset, Int pnrhs){
      SuperNode<T,Allocator> * contrib = this;

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

      nrhs += nrhsOffset;

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
          for(Int kk = 0; kk< ldfact; ++kk){
            for(Int j = nrhsOffset; j<nrhs;++j){
              //if non unit
              diag_nzval[kk*ldsol+j] = ( diag_nzval[kk*ldsol+j]) / (chol_nzval[kk*ldfact+kk]);

              for(Int i = kk+1; i<cur_nrows;++i){
                diag_nzval[i*ldsol+j] += -diag_nzval[kk*ldsol+j]*(chol_nzval[i*ldfact+kk]);
              }
            }
          }
        }
        else{
          for(Int kk = 0; kk<ldfact; ++kk){
            for(Int j = nrhsOffset; j<nrhs;++j){
              for(Int i = 0; i<cur_nrows;++i){
                cur_nzval[i*ldsol+j] += -diag_nzval[kk*ldsol+j]*(chol_nzval[i*ldfact+kk]);
              }
            }
          }
        }
      }
    }



  template <typename T, class Allocator> 
    inline void SuperNode<T,Allocator>::forward_update_contrib( T * RHS, SuperNode<T> * cur_snode, std::vector<Int> & perm){
      SuperNode<T,Allocator> * contrib = this;
      Int n = perm.size();
      Int nrhs = this->Size();

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
              //NOTE: RHS is stored in COLUMN major format, and is not permuted
              Int srcRow = perm[diag_desc.GIndex+kk-1];

              //if non unit
              diag_nzval[kk*nrhs+j] = (RHS[srcRow-1 + j*n] + diag_nzval[kk*nrhs+j]) / (chol_nzval[kk*cur_snode->Size()+kk]);

              for(Int i = kk+1; i<cur_nrows;++i){
                diag_nzval[i*nrhs+j] += -diag_nzval[kk*nrhs+j]*(chol_nzval[i*cur_snode->Size()+kk]);
              }
            }
          }
        }
        else{
          for(Int kk = 0; kk<cur_snode->Size(); ++kk){
            for(Int j = 0; j<nrhs;++j){
              for(Int i = 0; i<cur_nrows;++i){
                cur_nzval[i*nrhs+j] += -diag_nzval[kk*nrhs+j]*(chol_nzval[i*cur_snode->Size()+kk]);
              }
            }
          }
        }
      }
    }




  template <typename T,class Allocator> 
    inline size_t SuperNode<T,Allocator>::Deserialize(char * buffer, size_t size){
      this->Init(&buffer[0],size);
      this->InitIdxToBlk();
      return 0;
    }









  template <typename T, class Allocator> 
    inline void Serialize(Icomm & buffer,SuperNode<T, Allocator> & snode, Int first_blkidx, Idx first_row){
      NZBlockDesc* nzblk_ptr = &snode.GetNZBlockDesc(first_blkidx);
      if(first_row==0){
        first_row = nzblk_ptr->GIndex;
      }
      Idx local_first_row = first_row - nzblk_ptr->GIndex;

      Int nzblk_cnt = snode.NZBlockCnt() - first_blkidx;
      T* nzval_ptr = snode.GetNZval(nzblk_ptr->Offset) + local_first_row*snode.Size();

      size_t size = (char*)(nzblk_ptr+1)-(char*)nzval_ptr;
      buffer.clear();
      buffer.resize(size);
      //copy the whole thing in the buffer
      Serialize(buffer,(char*)nzval_ptr,size);
      //now we need to modify the first block data
      char * tail = buffer.front()+size; 
      NZBlockDesc* new_blk_ptr= (NZBlockDesc*)(tail-sizeof(NZBlockDesc));

      Int offset = new_blk_ptr->Offset + local_first_row*snode.Size();

      if(offset!=0){
        Int blkCnt = 1;
        if(blkCnt<nzblk_cnt){
          NZBlockDesc * curBlockPtr = nullptr;
          do{
            curBlockPtr = new_blk_ptr - blkCnt;
            curBlockPtr->Offset -= offset;
            ++blkCnt;
          }
          while(!curBlockPtr->Last);
        }
      }
      new_blk_ptr->Offset = 0;
      new_blk_ptr->GIndex=first_row;
    }





  template <typename T, class Allocator> inline std::ostream& operator<<( std::ostream& os,  SuperNode<T, Allocator>& snode){
    os<<"ooooooooooo   Supernode "<<snode.Id()<<" oooooooooooo"<<std::endl;
    os<<"     size = "<<snode.Size()<<std::endl;
    os<<"     fc   = "<<snode.FirstCol()<<std::endl;
    os<<"     lc   = "<<snode.LastCol()<<std::endl;
    os<<"     n    = "<<snode.N()<<std::endl;
    for(Int blkidx =0; blkidx<snode.NZBlockCnt();++blkidx){
      NZBlockDesc & nzblk_desc = snode.GetNZBlockDesc(blkidx);
      T * nzblk_nzval = snode.GetNZval(nzblk_desc.Offset);
      os<<"--- NZBlock "<<nzblk_desc.GIndex<<" ---"<<std::endl;
      std::streamsize p = logfileptr->OFS().precision();
      logfileptr->OFS().precision(std::numeric_limits< T >::max_digits10);
      for(Int i = 0; i<snode.NRows(blkidx);++i){
        for(Int j = 0; j<snode.Size();++j){
          os<<std::scientific<<ToMatlabScalar(nzblk_nzval[i*snode.Size()+j])<<" ";
        }
        os<<std::endl;
      }
      logfileptr->OFS().precision(p);
    }
    os<<"oooooooooooooooooooooooooooooooooooooooo"<<std::endl;
    return os;
  }


  inline std::ostream& operator<<( std::ostream& os,  NZBlockDesc& block){
    os<<"GIndex = "<<block.GIndex<<std::endl;
    os<<"Offset = "<<block.Offset<<std::endl;
    os<<"Last   = "<<(block.Last?std::string("true"):std::string("false"))<<std::endl;
    return os;
  }


  inline std::ostream& operator<<( std::ostream& os,  SuperNodeDesc& desc){
    os<<"ooooooooooo   Supernode "<<desc.iId_<<" oooooooooooo"<<std::endl;
    os<<"     size = "<<desc.iSize_<<std::endl;
    os<<"     fc   = "<<desc.iFirstCol_<<std::endl;
    os<<"     lc   = "<<desc.iLastCol_<<std::endl;
    os<<"     n    = "<<desc.iN_<<std::endl;
    os<<"   blocks = "<<desc.blocks_cnt_<<std::endl;
    os<<"   nzvals = "<<desc.nzval_cnt_<<std::endl;
    return os;
  }





  template <typename T,class Allocator> inline size_t Deserialize(char * buffer, size_t size, SuperNode<T,Allocator> & snode){
    snode.Init(&buffer[0],size);
    snode.InitIdxToBlk();
    return 0;
  }



} // namespace SYMPACK





#endif // _SUPERNODE_IMPL_HPP_
