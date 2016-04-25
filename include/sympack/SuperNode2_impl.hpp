#ifndef _SUPERNODE2_IMPL_HPP_
#define _SUPERNODE2_IMPL_HPP_

#ifdef _USE_COREDUMPER_
#include <google/coredumper.h>
#endif

//SuperNode2 implementation
template<typename T, class Allocator>
SuperNode2<T,Allocator>::SuperNode2() : meta_(NULL), blocks_(NULL), nzval_(NULL) { }

template<typename T, class Allocator>
SuperNode2<T,Allocator>::SuperNode2(Int aiId, Int aiFc, Int aiLc, Int ai_num_rows, Int aiN, Int aiNZBlkCnt) {

  //this is an upper bound
  assert(ai_num_rows>=0);

  //compute supernode size / width
  Int size = aiLc - aiFc +1;
  //compute maximum number of blocks, number of off-diagonal rows + 1
  Int num_blocks = ai_num_rows-size + 1;
  if(aiNZBlkCnt!=-1){
    num_blocks=aiNZBlkCnt;
  }

  storage_size_ = sizeof(T)*size*ai_num_rows + num_blocks*sizeof(NZBlockDesc2) + sizeof(SuperNodeDesc);
  loc_storage_container_ = Allocator::allocate(storage_size_);
  //storage_container_ = upcxx::allocate<char>(iam,storage_size_); 
  //loc_storage_container_ = (char *)storage_container_;
#ifdef _USE_COREDUMPER_
  if(loc_storage_container_==NULL){
    std::stringstream corename;
    corename << "core.sympack." << iam;
    WriteCoreDump(corename.str().c_str());
  }
#endif
  assert(loc_storage_container_!=NULL);

  nzval_ = (T*)&loc_storage_container_[0];
  meta_ = (SuperNodeDesc*)(nzval_+size*ai_num_rows);
  char * last = loc_storage_container_+storage_size_-1 - (sizeof(NZBlockDesc2) -1);
  blocks_ = (NZBlockDesc2*) last;

  meta_->iId_ = aiId;
  meta_->iFirstCol_ = aiFc;
  meta_->iLastCol_ = aiLc;
  meta_->iN_=aiN;
  meta_->iSize_ = size;
  meta_->nzval_cnt_ = 0;
  meta_->blocks_cnt_ = 0;
  meta_->b_own_storage_ = true;

#ifndef ITREE
  globalToLocal_ = new SYMPACK::vector<Int>(aiN+1,-1);
#else
  idxToBlk_ = CreateITree();
#endif
}; 

///  template<typename T>
///    SuperNode2<T>::SuperNode2(Int aiId, Int aiFc, Int aiLc, Int aiN) :iId_(aiId),iFirstCol_(aiFc),iLastCol_(aiLc) {
///
///      iN_=aiN;
///      //compute supernode size / width
///      iSize_ = iLastCol_ - iFirstCol_+1;
///
///      nzval_cnt_ = 0;
///      nzval_container_.reserve(iSize_*iSize_);
///      nzval_ = NULL;
///
///      blocks_container_.reserve(1);
///      blocks_ = NULL;
///      blocks_cnt_ = 0;
///
///
///#ifndef ITREE
///      //iLastRow_ = max(iLastCol_,iFirstCol_ + nzval_cnt_/iSize_ - 1);
///      //globalToLocal_ = new SYMPACK::vector<Int>(iLastRow_ - iFirstCol_ + 1 );
///      globalToLocal_ = new SYMPACK::vector<Int>(iN_+1,-1);
///#else
///      idxToBlk_ = CreateITree();
///#endif
///
///      b_own_storage_ = true;
///
///      //add the diagonal block
///      AddNZBlock(iSize_,0,iFirstCol_);
///
///    }; 

template<typename T, class Allocator>
SuperNode2<T,Allocator>::SuperNode2(char * storage_ptr,size_t storage_size, Int GIndex ) {
  //Init(aiId, aiFc, aiLc, a_block_desc, a_desc_cnt, a_nzval, a_nzval_cnt,aiN);
  Init(storage_ptr,storage_size,GIndex);
}

template<typename T, class Allocator>
SuperNode2<T,Allocator>::~SuperNode2(){
  if(meta_->b_own_storage_){
//    upcxx::deallocate(storage_container_);
    Allocator::deallocate(loc_storage_container_);
  }
#ifndef ITREE
  delete globalToLocal_;
#else
  delete idxToBlk_;
#endif
}

//CHECKED ON 11-18-2014
template<typename T, class Allocator>
void SuperNode2<T,Allocator>::Init(char * storage_ptr,size_t storage_size, Int GIndex ) {


  //loop through the block descriptors
  char * last = (char*)(storage_ptr+storage_size-1) - (sizeof(NZBlockDesc2) -1);

  blocks_ = (NZBlockDesc2*) last;

  Int blkCnt = 0;
  NZBlockDesc2 * curBlockPtr = NULL;
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
  globalToLocal_ = new SYMPACK::vector<Int>(meta_->iN_+1,-1);
#else
  idxToBlk_ = CreateITree();
#endif

}

template<typename T, class Allocator>
inline void SuperNode2<T,Allocator>::AddNZBlock(Int aiNRows, Int aiNCols, Int aiGIndex){

  //Resize the container if I own the storage
  if(meta_->b_own_storage_){
    scope_timer(a,RESIZE_SUPERNODE);

    Int cur_fr = aiGIndex;
    Int cur_lr = cur_fr + aiNRows -1;
    Int cur_nzval_cnt = aiNRows*meta_->iSize_;


#ifndef ITREE
    std::fill(&(*globalToLocal_)[cur_fr],&(*globalToLocal_)[cur_lr]+1,meta_->blocks_cnt_); 
#else
    ITree::Interval cur_interv = { cur_fr, cur_lr, meta_->blocks_cnt_};
    idxToBlk_->Insert(cur_interv);
#endif

    //if there is no more room for either nzval or blocks, extend
    Int block_space = (Int)(blocks_+1 - (NZBlockDesc2*)(meta_ +1)) - meta_->blocks_cnt_;
    Int nzval_space = (Int)((T*)meta_ - nzval_) - meta_->nzval_cnt_;

    if(block_space==0 || nzval_space<cur_nzval_cnt){
      //need to resize storage space. this is expensive !
      Int extra_nzvals = max(0,cur_nzval_cnt - nzval_space);
      Int extra_blocks = max(0,1 - block_space);
      size_t new_size = storage_size_ + extra_nzvals*sizeof(T) + extra_blocks*sizeof(NZBlockDesc2);
      size_t offset_meta = (char*)meta_ - (char*)nzval_;
      size_t offset_block = (char*)blocks_ - (char*)nzval_;

      //upcxx::global_ptr<char> tmpPtr = upcxx::allocate<char>(iam,new_size);
      //char * locTmpPtr = (char*)tmpPtr;
      char * locTmpPtr = Allocator::allocate(new_size);

//#ifdef _NO_MEMORY_PROGRESS_
//      while(locTmpPtr==NULL){
//        logfileptr->OFS()<<"No more memory, calling advance"<<endl;
//        upcxx::advance();
//        tmpPtr = upcxx::allocate<char>(iam,new_size);
//        locTmpPtr = (char*)tmpPtr;
//        break;
//      }
//#endif

#ifdef _USE_COREDUMPER_
      if(locTmpPtr==NULL){
        std::stringstream corename;
        corename << "core.sympack." << iam;
        WriteCoreDump(corename.str().c_str());
      }
#endif

      assert(locTmpPtr!=NULL);

      std::copy(loc_storage_container_,loc_storage_container_+storage_size_,locTmpPtr);
      //upcxx::deallocate(storage_container_);
      //storage_container_=tmpPtr;
      Allocator::deallocate(loc_storage_container_);
      loc_storage_container_=locTmpPtr;
      storage_size_=new_size;
      //storage_container_.resize(new_size);

      nzval_=(T*)&loc_storage_container_[0];
      //move the block descriptors if required
      char * cur_blocks_ptr = (char*)&loc_storage_container_[0] + offset_block;

      //move the meta data if required
      char * cur_meta_ptr = (char*)&loc_storage_container_[0] + offset_meta;
      meta_ = (SuperNodeDesc*) cur_meta_ptr;
      //we need to move everything, starting from the blocks, then meta
      //blocks need to be moved by extra_nzvals + extra_blocks
      char * new_blocks_ptr = cur_blocks_ptr + extra_nzvals*sizeof(T) +extra_blocks*sizeof(NZBlockDesc2);
      std::copy_backward(cur_blocks_ptr - (meta_->blocks_cnt_-1)*sizeof(NZBlockDesc2),cur_blocks_ptr+sizeof(NZBlockDesc2),new_blocks_ptr+sizeof(NZBlockDesc2));

      //now move the meta data by extra_nzvals
      char * new_meta_ptr = cur_meta_ptr + extra_nzvals*sizeof(T);
      std::copy(cur_meta_ptr,cur_meta_ptr + sizeof(SuperNodeDesc),new_meta_ptr);

      //update pointers
      meta_ = (SuperNodeDesc*) new_meta_ptr;
      blocks_ = (NZBlockDesc2*) new_blocks_ptr;
    }

    GetNZBlockDesc(meta_->blocks_cnt_).GIndex = aiGIndex;
    GetNZBlockDesc(meta_->blocks_cnt_).Offset = meta_->nzval_cnt_;
    GetNZBlockDesc(meta_->blocks_cnt_).Last = true;
    if(meta_->blocks_cnt_>0){
      GetNZBlockDesc(meta_->blocks_cnt_-1).Last = false;
    }

    //blocks_container_.push_back(NZBlockDesc2(aiGIndex,nzval_cnt_));

    meta_->blocks_cnt_++;

    //fill the new block with zeros
    std::fill(nzval_+meta_->nzval_cnt_,nzval_+meta_->nzval_cnt_+cur_nzval_cnt,ZERO<T>());

    //update nzval count
    meta_->nzval_cnt_+=cur_nzval_cnt;

  }
}


template<typename T, class Allocator>
inline Int SuperNode2<T,Allocator>::FindBlockIdx(Int aiGIndex){
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


  ITree::Interval * res = idxToBlk_->IntervalSearch(aiGIndex,aiGIndex);
  if (res != NULL){
    rval = res->block_idx;
  }
#endif
  return rval;
}


template<typename T, class Allocator>
inline Int SuperNode2<T,Allocator>::FindBlockIdx(Int aiGIndex,Int & closestR, Int & closestL){
  scope_timer(a,FindBlockIdx);
  Int rval = -1;
  ITree::Interval * L = NULL;
  ITree::Interval * R = NULL;
  ITree::Interval * res = idxToBlk_->IntervalSearch(aiGIndex,aiGIndex,R,L);
  if(R!=NULL){
    closestR = R->low;
  }

  if(L!=NULL){
    closestL = L->high;
  }

  if (res != NULL){
    rval = res->block_idx;
  }
  return rval;
}




template<typename T, class Allocator>
inline Int SuperNode2<T,Allocator>::FindBlockIdx(Int fr, Int lr, ITree::Interval & overlap){
  scope_timer(a,FindBlockIdx);

  Int rval = -1;

  ITree::Interval * res = idxToBlk_->IntervalSearch(fr,lr);
  if (res != NULL){
    overlap = *res;
    rval = res->block_idx;
  }
  return rval;
}





template<typename T, class Allocator>
inline void SuperNode2<T,Allocator>::DumpITree(){
#ifdef ITREE
  logfileptr->OFS()<<"Number of blocks: "<<meta_->blocks_cnt_<<endl;
  logfileptr->OFS()<<"log2(Number of blocks): "<<log2(meta_->blocks_cnt_)<<endl;
  idxToBlk_->Dump();
#endif
}

template<typename T, class Allocator>
inline Int SuperNode2<T,Allocator>::Shrink(){
  if(meta_->b_own_storage_){
    //TODO make sure that we do not have any extra space anywhere.

    //if there is too much room for either nzval or blocks, contract
    Int block_space = (Int)(blocks_+1 - (NZBlockDesc2*)(meta_ +1)) - meta_->blocks_cnt_;
    Int nzval_space = (Int)((T*)meta_ - nzval_) - meta_->nzval_cnt_;

    if(block_space >0 || nzval_space >0){

      size_t new_size = storage_size_ - nzval_space*sizeof(T) - block_space*sizeof(NZBlockDesc2);
      size_t offset_meta = meta_->nzval_cnt_*sizeof(T);
      size_t offset_block = meta_->nzval_cnt_*sizeof(T)+sizeof(SuperNodeDesc);//+meta_->blocks_cnt_*sizeof(NZBlockDesc2);

      //upcxx::global_ptr<char> tmpPtr = upcxx::allocate<char>(iam,new_size);
      //char * locTmpPtr = (char*)tmpPtr;
      char * locTmpPtr = Allocator::allocate(new_size);


#ifdef _USE_COREDUMPER_
      if(locTmpPtr==NULL){
        std::stringstream corename;
        corename << "core.sympack." << iam;
        WriteCoreDump(corename.str().c_str());
      }
#endif


      assert(locTmpPtr!=NULL);

      //copy nzvals
      std::copy(loc_storage_container_,loc_storage_container_+meta_->nzval_cnt_*sizeof(T),locTmpPtr);
      //copy meta
      std::copy(meta_,meta_+1,(SuperNodeDesc*)(locTmpPtr+offset_meta));
      //copy blocks
      std::copy(blocks_-meta_->blocks_cnt_+1,blocks_+1,(NZBlockDesc2*)(locTmpPtr+offset_block));

      //upcxx::deallocate(storage_container_);
      //storage_container_=tmpPtr;
      Allocator::deallocate(loc_storage_container_);
      loc_storage_container_=locTmpPtr;
      storage_size_=new_size;

      nzval_=(T*)&loc_storage_container_[0];
      //move the block descriptors if required
      char * new_blocks_ptr = loc_storage_container_+storage_size_-1 - (sizeof(NZBlockDesc2) -1);

      //move the meta data if required
      char * new_meta_ptr = loc_storage_container_ + offset_meta;

      //update pointers
      meta_ = (SuperNodeDesc*) new_meta_ptr;
      blocks_ = (NZBlockDesc2*) new_blocks_ptr;
    }




    //        blocks_ = &blocks_container_.front();
    //        nzval_container_.resize(nzval_cnt_);
    //        nzval_ = &nzval_container_.front();
  }

  return StorageSize();
}





template<typename T, class Allocator>
inline void SuperNode2<T,Allocator>::FindUpdatedFirstCol(SuperNode2<T> & src_snode, Int pivot_fr, Int & tgt_fc, Int & first_pivot_idx){
  //find the first row updated by src_snode
  scope_timer(a,UPDATE_SNODE_FIND_INDEX);


#ifdef _LINEAR_SEARCH_FCLC_
  if(src_snode.ITreeInitialized()){
#endif 
    tgt_fc = pivot_fr;
    first_pivot_idx = -1;
    if(tgt_fc == I_ZERO ){
      tgt_fc = FirstCol();
      //find the pivot idx

      Int tgt_lc = LastCol();
#ifndef _BINARY_BLOCK_SEARCH_
//    for(tgt_fc;tgt_fc<tgt_lc;tgt_fc++){
//      first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
//      if(first_pivot_idx>=0){
//        break;
//      }
//    }
      do {first_pivot_idx = src_snode.FindBlockIdx(tgt_fc); tgt_fc++;}
      while(first_pivot_idx<0 && tgt_fc<=tgt_lc);
      tgt_fc--;
#else
      Int closestR = -1;
      Int closestL = -1;
      do{
        first_pivot_idx = src_snode.FindBlockIdx(tgt_fc,closestR,closestL);
        if(closestR!=-1){
          logfileptr->OFS()<<"ClosestR of "<<tgt_fc<<" is "<<closestR<<endl;
          tgt_fc = closestR;
        }
      }
      while(first_pivot_idx<0);

//    Int tgt_cur;
//    Int count = tgt_lc - tgt_fc+1;
//    Int step;
//    while (count > 0) {
//        step = count / 2;
//        tgt_cur = tgt_fc + step;
//        first_pivot_idx = src_snode.FindBlockIdx(tgt_cur); 
//        if (first_pivot_idx>=0) {
//            tgt_fc = ++tgt_cur;
//            count -= step + 1;
//        }
//        else{
//            count = step;
//        }
//    }
#endif

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
        NZBlockDesc2 & cur_block = src_snode.GetNZBlockDesc(blkidx);
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
        NZBlockDesc2 & cur_block = src_snode.GetNZBlockDesc(blkidx);
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
  if(first_pivot_idx<0){
    logfileptr->OFS()<<"LOCK 1: first_pivot_idx<0"<<endl;
    gdb_lock();
  }
  assert(first_pivot_idx>=0);
}

template<typename T, class Allocator>
inline void SuperNode2<T,Allocator>::FindUpdatedLastCol(SuperNode2<T> & src_snode, Int tgt_fc, Int first_pivot_idx , Int & tgt_lc,  Int & last_pivot_idx){
  //gdb_lock();
  scope_timer(a,UPDATE_SNODE_FIND_INDEX);
#ifdef _LINEAR_SEARCH_FCLC_
  if(src_snode.ITreeInitialized()){
#endif 
    //find the last row updated by src_snode
    tgt_lc = LastCol();
    last_pivot_idx = -1;
    //find the pivot idx
#ifndef _BINARY_BLOCK_SEARCH_
//    for(tgt_lc;tgt_lc>=tgt_fc;tgt_lc--){
//      last_pivot_idx = src_snode.FindBlockIdx(tgt_lc);
//      if(last_pivot_idx>=0){
//        break;
//      }
//    }
    do {last_pivot_idx = src_snode.FindBlockIdx(tgt_lc); tgt_lc--;}
    while(last_pivot_idx<0 && tgt_lc>=tgt_fc);
    tgt_lc++;
#else
    for(tgt_lc;tgt_lc>=tgt_fc;tgt_lc--){
      last_pivot_idx = src_snode.FindBlockIdx(tgt_lc);
      if(last_pivot_idx>=0){
        break;
      }
    }

//    do {last_pivot_idx = src_snode.FindBlockIdx(tgt_lc); tgt_lc--;  logfileptr->OFS()<<"lpi: "<<last_pivot_idx<<" | "<<tgt_lc+1<<endl; }
//    while(last_pivot_idx<0 && tgt_lc>=tgt_fc);
//    tgt_lc++;
//    Int bak_tgt_lc = tgt_lc;
//    tgt_lc = LastCol();
//    Int tgt_cur;
//    Int count = tgt_lc - tgt_fc+1;
//    Int step;
//    //last_pivot_idx = src_snode.FindBlockIdx(tgt_lc);
//    //if(last_pivot_idx<0)
//    {
//      while (count > 0) {
//        step = count / 2;
//        tgt_cur = tgt_lc - step;
//        last_pivot_idx = src_snode.FindBlockIdx(tgt_cur); 
//        logfileptr->OFS()<<"lpi2: "<<last_pivot_idx<<" | "<<tgt_cur<<endl;
//        if (last_pivot_idx<0) {
//          tgt_lc = --tgt_cur;
//          count -= step - 1;
//        }
//        else{
//          count = step;
//        }
//      }
//    }
//    logfileptr->OFS()<<bak_tgt_lc<<" vs "<<tgt_lc<<endl;
//    logfileptr->OFS()<<src_snode.FindBlockIdx(tgt_lc+1)<<endl;
//    logfileptr->OFS()<<src_snode.FindBlockIdx(tgt_lc)<<endl;
//    assert(src_snode.FindBlockIdx(tgt_lc+1)<0 && src_snode.FindBlockIdx(tgt_lc)>=0);
//    assert(bak_tgt_lc == tgt_lc);
#endif






#ifdef _LINEAR_SEARCH_FCLC_
  }
  else{
    tgt_lc = tgt_fc;
    for(Int blkidx = first_pivot_idx; blkidx<src_snode.NZBlockCnt();++blkidx){
      NZBlockDesc2 & cur_block = src_snode.GetNZBlockDesc(blkidx);
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
}







//#ifdef COMPACT_AGGREGATES
template<typename T, class Allocator>
inline Int SuperNode2<T,Allocator>::Merge(SuperNode2<T> & src_snode, SnodeUpdate &update){
  //#ifndef _TAU_TRACE_
  scope_timer(a,MERGE_SNODE);

  assert(meta_->b_own_storage_);

  Int src_snode_size = src_snode.Size();
  Int tgt_snode_size = Size();

  Int & pivot_idx = update.blkidx;
  Int & pivot_fr = update.src_first_row;

  Int tgt_fc;
  Int first_pivot_idx;
  FindUpdatedFirstCol(src_snode, 0, tgt_fc, first_pivot_idx);
  NZBlockDesc2 & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);

  //parse src_snode
  ITree::Interval overlap;
  ITree::Interval curInter;
  ITree::Interval newInter;
  std::queue< ITree::Interval > toInsert;
  for(Int blkidx = first_pivot_idx; blkidx<src_snode.NZBlockCnt(); ++blkidx){
    NZBlockDesc2 & blk_desc = src_snode.GetNZBlockDesc(blkidx);
    Int fr = max(FirstCol(),blk_desc.GIndex);
    Int lr = blk_desc.GIndex + src_snode.NRows(blkidx) -1;

    curInter.low = fr;
    curInter.high = lr;
    toInsert.push(curInter);

    while(!toInsert.empty()){
      curInter = toInsert.front();
      toInsert.pop();
      if(FindBlockIdx(curInter.low,curInter.high,overlap)==-1){
        //Add the full block
        AddNZBlock( curInter.high - curInter.low + 1, tgt_snode_size, curInter.low);
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
//#endif

//#ifndef _TAU_TRACE_
//#ifdef COMPACT_AGGREGATES
template<typename T, class Allocator>
inline Int SuperNode2<T,Allocator>::Aggregate(SuperNode2<T> & src_snode){

  scope_timer(a,AGGREGATE_SNODE);

#if defined(_NO_COMPUTATION_)
  return 0;
#endif

  Int  pivot_idx = 0;
  Int  pivot_fr = 0;

  Int src_snode_size = src_snode.Size();
  Int tgt_snode_size = Size();

  //parse src_snode and add everything

  Int first_pivot_idx = 0 ;
  Int tgt_fc = FirstCol();

  for(Int blkidx = first_pivot_idx; blkidx<src_snode.NZBlockCnt(); ++blkidx){
    NZBlockDesc2 & blk_desc = src_snode.GetNZBlockDesc(blkidx);
    Int nrows = src_snode.NRows(blkidx);
    for(Int rowidx = 0; rowidx<nrows; ++rowidx){
      Int row = blk_desc.GIndex + rowidx;

      if(row>=tgt_fc){
        Int src_offset = blk_desc.Offset + (row - blk_desc.GIndex)*src_snode_size;

        Int tgt_blkidx = FindBlockIdx(row);
        if(tgt_blkidx==-1){
          logfileptr->OFS()<<"LOCK 2: tgt_blkidx=-1"<<endl;
          gdb_lock();
        }
        NZBlockDesc2 & tgt_desc = GetNZBlockDesc(tgt_blkidx);
        Int tgt_offset = tgt_desc.Offset + (row - tgt_desc.GIndex)*tgt_snode_size;

        T * src = src_snode.GetNZval(src_offset);
        T * tgt = GetNZval(tgt_offset);

        //blas::Axpy(tgt_snode_size,ONE<T>(),src,1,tgt,1);
        #pragma omp simd
        for(Int i = 0; i< tgt_snode_size;i+=1){ tgt[i] += src[i]; }

      }
    }
  }


  return 0;
}


//#ifdef COMPACT_AGGREGATES
//CHECKED ON 11-18-2014
template<typename T, class Allocator>
inline Int SuperNode2<T,Allocator>::UpdateAggregate(SuperNode2<T> & src_snode, SnodeUpdate &update, 
    TempUpdateBuffers<T> & tmpBuffers, Int iTarget){

  scope_timer(a,UPDATE_AGGREGATE_SNODE);
#if defined(_NO_COMPUTATION_)
  return 0;
#endif

  if(iTarget != iam){
    Merge(src_snode, update);

    Int & pivot_idx = update.blkidx;
    Int & pivot_fr = update.src_first_row;

    Int src_snode_size = src_snode.Size();
    Int tgt_snode_size = Size();

    Int tgt_fc,tgt_lc;
    Int first_pivot_idx,last_pivot_idx;
    FindUpdatedFirstCol(src_snode, update.src_first_row, tgt_fc, first_pivot_idx);
    FindUpdatedLastCol(src_snode, tgt_fc, first_pivot_idx, tgt_lc, last_pivot_idx);

    NZBlockDesc2 & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);
    NZBlockDesc2 & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);

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
#ifdef _DEBUG_
    tmpBuffers.tmpBuf.Resize(tgt_width,src_nrows);
#endif

    buf = tmpBuffers.tmpBuf.Data();

    //everything is in row-major
    TIMER_START(UPDATE_SNODE_GEMM);
    blas::Gemm('T','N',tgt_width, src_nrows,src_snode_size,
        MINUS_ONE<T>(),pivot,src_snode_size,
        pivot,src_snode_size,beta,buf,tgt_width);
    TIMER_STOP(UPDATE_SNODE_GEMM);

    //If the GEMM wasn't done in place we need to aggregate the update
    //This is the assembly phase
#ifdef _DEBUG_
    logfileptr->OFS()<<"tmpBuf is "<<tmpBuffers.tmpBuf<<std::endl;
#endif

    //now add the update to the target supernode
    TIMER_START(UPDATE_SNODE_INDEX_MAP);
    if(tgt_snode_size==1){
      Int rowidx = 0;
      Int src_blkcnt = src_snode.NZBlockCnt();
      for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
        NZBlockDesc2 & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
        Int cur_src_nrows = src_snode.NRows(blkidx);
        Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
        Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);

        Int row = cur_src_fr;
        while(row<=cur_src_lr){
          Int tgt_blk_idx = FindBlockIdx(row);
          assert(tgt_blk_idx>=0);
          NZBlockDesc2 & cur_tgt_desc = GetNZBlockDesc(tgt_blk_idx);
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
      tmpBuffers.src_colindx.resize(tgt_width);
      tmpBuffers.src_to_tgt_offset.resize(src_nrows);
      Int colidx = 0;
      Int rowidx = 0;
      Int offset = 0;

      Int src_blkcnt = src_snode.NZBlockCnt();
      for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
        NZBlockDesc2 & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
        Int cur_src_nrows = src_snode.NRows(blkidx);
        Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
        Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
        cur_src_nrows = cur_src_lr - cur_src_fr +1;

        Int row = cur_src_fr;
        while(row<=cur_src_lr){
          Int tgt_blk_idx = FindBlockIdx(row);
          assert(tgt_blk_idx>=0);
          NZBlockDesc2 & cur_tgt_desc = GetNZBlockDesc(tgt_blk_idx);
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
      TIMER_STOP(UPDATE_SNODE_INDEX_MAP);


      //Multiple cases to consider
      if(first_pivot_idx==last_pivot_idx){
        // Updating contiguous columns
        Int tgt_offset = (tgt_fc - FirstCol());
        for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
          T * A = &buf[rowidx*tgt_width];
          T * B = &tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_offset];
          #pragma omp simd
          for(Int i = 0; i < tgt_width; ++i){
            B[i] += A[i];
          }
//          blas::Axpy(tgt_width,ONE<T>(),&buf[rowidx*tgt_width],1,
//              &tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_offset],1);
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
    Update(src_snode, update, tmpBuffers);
  }





}


//CHECKED ON 11-18-2014
template<typename T, class Allocator>
inline Int SuperNode2<T,Allocator>::Update(SuperNode2<T> & src_snode, SnodeUpdate &update, 
    TempUpdateBuffers<T> & tmpBuffers){

  scope_timer(a,UPDATE_SNODE);
#if defined(_NO_COMPUTATION_)
  return 0;
#endif

  Int & pivot_idx = update.blkidx;
  Int & pivot_fr = update.src_first_row;

  Int src_snode_size = src_snode.Size();
  Int tgt_snode_size = Size();

  //find the first row updated by src_snode
  Int tgt_fc,tgt_lc;
  Int first_pivot_idx,last_pivot_idx;
  FindUpdatedFirstCol(src_snode, update.src_first_row, tgt_fc, first_pivot_idx);
  FindUpdatedLastCol(src_snode, tgt_fc, first_pivot_idx, tgt_lc, last_pivot_idx);

  NZBlockDesc2 & first_pivot_desc = src_snode.GetNZBlockDesc(first_pivot_idx);
  NZBlockDesc2 & last_pivot_desc = src_snode.GetNZBlockDesc(last_pivot_idx);

  //determine the first column that will be updated in the target supernode
  Int tgt_local_fc =  tgt_fc - FirstCol();
  Int tgt_local_lc =  tgt_lc - FirstCol();

  Int tgt_nrows = NRowsBelowBlock(0);
  Int src_nrows = src_snode.NRowsBelowBlock(first_pivot_idx)
    - (tgt_fc - first_pivot_desc.GIndex);
  Int src_lr = tgt_fc+src_nrows-1;
  src_nrows = src_lr - tgt_fc + 1;

  Int src_belowLast = src_snode.NRowsBelowBlock(last_pivot_idx);
  Int tgt_width = src_nrows - src_belowLast
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
        NZBlockDesc2 & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
        Int cur_src_nrows = src_snode.NRows(blkidx);
        Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
        Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);

        Int row = cur_src_fr;
        while(row<=cur_src_lr){
          Int tgt_blk_idx = FindBlockIdx(row);
          if(tgt_blk_idx<0){src_snode.DumpITree(); DumpITree(); 

            logfileptr->OFS()<<"LOCK 3: tgt_blk_idx<0"<<endl;
            gdb_lock();}
            assert(tgt_blk_idx>=0);
            NZBlockDesc2 & cur_tgt_desc = GetNZBlockDesc(tgt_blk_idx);
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
      tmpBuffers.src_colindx.resize(tgt_width);
      tmpBuffers.src_to_tgt_offset.resize(src_nrows);
      Int colidx = 0;
      Int rowidx = 0;
      Int offset = 0;

      Int src_blkcnt = src_snode.NZBlockCnt();
      for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
        NZBlockDesc2 & cur_block_desc = src_snode.GetNZBlockDesc(blkidx);
        Int cur_src_nrows = src_snode.NRows(blkidx);
        Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
        Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
        cur_src_nrows = cur_src_lr - cur_src_fr +1;

        //The other one MUST reside into a single block in the target
        Int row = cur_src_fr;
        while(row<=cur_src_lr){
          Int tgt_blk_idx = FindBlockIdx(row);
          if(tgt_blk_idx<0){
            logfileptr->OFS()<<"LOCK 4: tgt_blk_idx<0"<<endl;
            gdb_lock();}
            assert(tgt_blk_idx>=0);
            NZBlockDesc2 & cur_tgt_desc = GetNZBlockDesc(tgt_blk_idx);
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
          T * A = &buf[rowidx*tgt_width];
          T * B = &tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_offset];
          #pragma omp simd
          for(Int i = 0; i < tgt_width; ++i){
            B[i] += A[i];
          }
//          blas::Axpy(tgt_width,ONE<T>(),&buf[rowidx*tgt_width],1,
//              &tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_offset],1);
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
    }
  }
  return 0;
}


//CHECKED ON 11-18-2014
template<typename T, class Allocator>
  inline Int SuperNode2<T,Allocator>::Factorize(){
#if defined(_NO_COMPUTATION_)
    return 0;
#endif


    Int BLOCKSIZE = Size();
    NZBlockDesc2 & diag_desc = GetNZBlockDesc(0);
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

  }




template<typename T, class Allocator>
bool SuperNode2<T,Allocator>::FindNextUpdate(SnodeUpdate & nextUpdate, const SYMPACK::vector<Int> & Xsuper,  const SYMPACK::vector<Int> & SupMembership, bool isLocal){
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
    NZBlockDesc2 * cur_desc = &GetNZBlockDesc(f_ub); 
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
        NZBlockDesc2 & cur_block = GetNZBlockDesc(blkidx);
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
  inline ITree * SuperNode2<T,Allocator>::CreateITree(){
#if defined(_AVL_ITREE_)
    return new AVLITree();
#elif defined(_DSW_ITREE_)
    return new DSWITree();
#else
    return new ITree();
#endif
  } 

template<typename T, class Allocator>
inline void SuperNode2<T,Allocator>::InitIdxToBlk(){
#ifndef ITREE
  TIMER_START(ARRAY_INSERT);
  for(Int blkidx=0; blkidx<meta_->blocks_cnt_;++blkidx){
    Int cur_fr = GetNZBlockDesc(blkidx).GIndex;
    Int cur_lr = cur_fr + NRows(blkidx) -1;

    std::fill(&(*globalToLocal_)[cur_fr],&(*globalToLocal_)[cur_lr]+1,blkidx); 
  }
  TIMER_STOP(ARRAY_INSERT);
#else
  TIMER_START(BST_INSERT);
  for(Int blkidx=0; blkidx<meta_->blocks_cnt_;++blkidx){
    Int cur_fr = GetNZBlockDesc(blkidx).GIndex;
    Int cur_lr = cur_fr + NRows(blkidx) -1;

    ITree::Interval cur_interv = { cur_fr, cur_lr, blkidx};
    idxToBlk_->Insert(cur_interv);
  }
  TIMER_STOP(BST_INSERT);
  idxToBlk_->Rebalance();
#endif
}


template<typename T, class Allocator>
inline void SuperNode2<T,Allocator>::Reserve(size_t storage_size){
}

//  template <typename T> inline size_t Deserialize(char * buffer, SuperNode2<T> & snode){
//    Int snode_id = *(Int*)&buffer[0];
//    Int snode_fc = *(((Int*)&buffer[0])+1);
//    Int snode_lc = *(((Int*)&buffer[0])+2);
//    Int n = *(((Int*)&buffer[0])+3);
//    Int nzblk_cnt = *(((Int*)&buffer[0])+4);
//    NZBlockDesc * blocks_ptr = 
//      reinterpret_cast<NZBlockDesc*>(&buffer[5*sizeof(Int)]);
//    Int nzval_cnt_ = *(Int*)(blocks_ptr + nzblk_cnt);
//    T * nzval_ptr = (T*)((Int*)(blocks_ptr + nzblk_cnt)+1);
//    char * last_ptr = (char *)(nzval_ptr+nzval_cnt_);
//
//    //Create the dummy supernode for that data
//    snode.Init(snode_id,snode_fc,snode_lc, blocks_ptr, nzblk_cnt, nzval_ptr, nzval_cnt_,n);
//
//    return (last_ptr - buffer);
//  }


template <typename T, class Allocator> inline void Serialize(Icomm & buffer,SuperNode2<T, Allocator> & snode, Int first_blkidx, Int first_row){
  NZBlockDesc2* nzblk_ptr = &snode.GetNZBlockDesc(first_blkidx);
  Int local_first_row = first_row - nzblk_ptr->GIndex;

  Int nzblk_cnt = snode.NZBlockCnt() - first_blkidx;
  T* nzval_ptr = snode.GetNZval(nzblk_ptr->Offset) + local_first_row*snode.Size();
  if(local_first_row!=0){

    logfileptr->OFS()<<"LOCK 5: serialize"<<endl;
    gdb_lock();}

    size_t size = (char*)(nzblk_ptr+1)-(char*)nzval_ptr;
    buffer.clear();
    buffer.resize(size);
    //copy the whole thing in the buffer
    Serialize(buffer,(char*)nzval_ptr,size);
    //now we need to modify the first block data
    char * tail = buffer.front()+size; 
    NZBlockDesc2* new_blk_ptr= (NZBlockDesc2*)(tail-sizeof(NZBlockDesc2));

    Int offset = new_blk_ptr->Offset + local_first_row*snode.Size();

    Int blkCnt = 1;
    if(blkCnt<nzblk_cnt){
      NZBlockDesc2 * curBlockPtr = NULL;
      do{
        curBlockPtr = new_blk_ptr - blkCnt;
        curBlockPtr->Offset -= offset;
        ++blkCnt;
      }
      while(!curBlockPtr->Last);
    }
    new_blk_ptr->Offset = 0;
}





template <typename T, class Allocator> inline std::ostream& operator<<( std::ostream& os,  SuperNode2<T, Allocator>& snode){
  os<<"ooooooooooo   Supernode "<<snode.Id()<<" oooooooooooo"<<std::endl;
  os<<"     size = "<<snode.Size()<<std::endl;
  os<<"     fc   = "<<snode.FirstCol()<<std::endl;
  os<<"     lc   = "<<snode.LastCol()<<std::endl;
  os<<"     n    = "<<snode.N()<<std::endl;
  for(Int blkidx =0; blkidx<snode.NZBlockCnt();++blkidx){
    NZBlockDesc2 & nzblk_desc = snode.GetNZBlockDesc(blkidx);
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









#endif // _SUPERNODE2_IMPL_HPP_
