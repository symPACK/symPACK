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
/*           |      Meta       |        */
/*           |_________________|        */
/*           |                 |        */
/*           |   block desc    |        */
/*           |                 |        */
/*           |_________________|        */
/*                                      */
/****************************************/

namespace SYMPACK{

//SuperNode implementation
template<typename T, class Allocator>
SuperNode<T,Allocator>::SuperNode() : meta_(NULL), blocks_(NULL), nzval_(NULL) { }

template<typename T, class Allocator>
SuperNode<T,Allocator>::SuperNode(Int aiId, Int aiFc, Int aiLc, Int ai_num_rows, Int aiN, Int aiNZBlkCnt) {

  //this is an upper bound
  assert(ai_num_rows>=0);

  //compute supernode size / width
  Int size = aiLc - aiFc +1;
  //compute maximum number of blocks, number of off-diagonal rows + 1
  Int num_blocks = max(1,ai_num_rows-size + 1);
  if(aiNZBlkCnt!=-1){
    num_blocks=aiNZBlkCnt;
  }

  storage_size_ = sizeof(T)*size*ai_num_rows + num_blocks*sizeof(NZBlockDesc) + sizeof(SuperNodeDesc);

  loc_storage_container_ = Allocator::allocate(storage_size_);
  //storage_container_ = upcxx::allocate<char>(iam,storage_size_); 
  //loc_storage_container_ = (char *)storage_container_;
#ifdef _USE_COREDUMPER_
  if(loc_storage_container_==NULL){
gdb_lock();
    std::stringstream corename;
    corename << "core.sympack." << iam;
    WriteCoreDump(corename.str().c_str());
  }
#endif
  assert(loc_storage_container_!=NULL);

  nzval_ = (T*)&loc_storage_container_[0];
  meta_ = (SuperNodeDesc*)(nzval_+size*ai_num_rows);
  char * last = loc_storage_container_+storage_size_-1 - (sizeof(NZBlockDesc) -1);
  blocks_ = (NZBlockDesc*) last;

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






template<typename T, class Allocator>
SuperNode<T,Allocator>::SuperNode(Int aiId, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices) {


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
  storage_size_ = sizeof(T)*size*numRows + num_blocks*sizeof(NZBlockDesc) + sizeof(SuperNodeDesc);

  loc_storage_container_ = Allocator::allocate(storage_size_);
  assert(loc_storage_container_!=NULL);

  nzval_ = (T*)&loc_storage_container_[0];
  meta_ = (SuperNodeDesc*)(nzval_+size*numRows);
  char * last = loc_storage_container_+storage_size_-1 - (sizeof(NZBlockDesc) -1);
  blocks_ = (NZBlockDesc*) last;

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



  //now add the blocks 
  if(rowIndices.size()>0){
    //go through the set and count the number of nz blocks
    Idx prevRow = *rowIndices.begin();
    Idx firstRow = *rowIndices.begin();
    for(auto it = rowIndices.begin();it!=rowIndices.end();it++){
      Idx row = *it;

      if(row>prevRow+1){
        this->AddNZBlock( prevRow - firstRow + 1, size, firstRow);
        firstRow = row;
      }
      prevRow = row;
    }
    this->AddNZBlock( prevRow - firstRow + 1, size, firstRow);
  }










}; 

















template<typename T, class Allocator>
SuperNode<T,Allocator>::SuperNode(char * storage_ptr,size_t storage_size, Int GIndex ) {
  //Init(aiId, aiFc, aiLc, a_block_desc, a_desc_cnt, a_nzval, a_nzval_cnt,aiN);
  Init(storage_ptr,storage_size,GIndex);
}

template<typename T, class Allocator>
SuperNode<T,Allocator>::~SuperNode(){
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
void SuperNode<T,Allocator>::Init(char * storage_ptr,size_t storage_size, Int GIndex ) {


  //loop through the block descriptors
  char * last = (char*)(storage_ptr+storage_size-1) - (sizeof(NZBlockDesc) -1);

  blocks_ = (NZBlockDesc*) last;

  Int blkCnt = 0;
  NZBlockDesc * curBlockPtr = NULL;
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
inline void SuperNode<T,Allocator>::AddNZBlock(Int aiNRows, Int aiNCols, Int aiGIndex){

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
    Int block_space = (Int)(blocks_+1 - (NZBlockDesc*)(meta_ +1)) - meta_->blocks_cnt_;
    Int nzval_space = (Int)((T*)meta_ - nzval_) - meta_->nzval_cnt_;


    if(block_space==0 || nzval_space<cur_nzval_cnt){
      //need to resize storage space. this is expensive !
      Int size = Size();
      Int extra_nzvals = max(0,cur_nzval_cnt - nzval_space);
      Int extra_blocks = max(0,1 - block_space);
      size_t new_size = storage_size_ + extra_nzvals*sizeof(T) + extra_blocks*sizeof(NZBlockDesc);
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

      //std::copy(loc_storage_container_,loc_storage_container_+storage_size_,locTmpPtr);
      for(size_t i = 0;i<storage_size_;i++){locTmpPtr[i] = loc_storage_container_[i];}


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
      char * new_blocks_ptr = cur_blocks_ptr + extra_nzvals*sizeof(T) +extra_blocks*sizeof(NZBlockDesc);
      std::copy_backward(cur_blocks_ptr - (meta_->blocks_cnt_-1)*sizeof(NZBlockDesc),cur_blocks_ptr+sizeof(NZBlockDesc),new_blocks_ptr+sizeof(NZBlockDesc));

      //now move the meta data by extra_nzvals
      char * new_meta_ptr = cur_meta_ptr + extra_nzvals*sizeof(T);
      std::copy(cur_meta_ptr,cur_meta_ptr + sizeof(SuperNodeDesc),new_meta_ptr);


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

    //blocks_container_.push_back(NZBlockDesc(aiGIndex,nzval_cnt_));

    meta_->blocks_cnt_++;

    //fill the new block with zeros
    std::fill(nzval_+meta_->nzval_cnt_,nzval_+meta_->nzval_cnt_+cur_nzval_cnt,ZERO<T>());

    //update nzval count
    meta_->nzval_cnt_+=cur_nzval_cnt;

  }
}

template<typename T, class Allocator>
inline Int SuperNode<T,Allocator>::FindBlockIdx(Int aiGIndex){
  scope_timer_special(a,FindBlockIdx);


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
inline Int SuperNode<T,Allocator>::FindBlockIdx(Int aiGIndex,Int & closestR, Int & closestL){
  scope_timer_special(a,FindBlockIdxRL);
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
inline Int SuperNode<T,Allocator>::FindBlockIdx(Int fr, Int lr, ITree::Interval & overlap){
  scope_timer_special(a,FindBlockIdxOverlap);

  Int rval = -1;

  ITree::Interval * res = idxToBlk_->IntervalSearch(fr,lr);
  if (res != NULL){
    overlap = *res;
    rval = res->block_idx;
  }
  return rval;
}

template<typename T, class Allocator>
inline void SuperNode<T,Allocator>::DumpITree(){
#ifdef ITREE
  logfileptr->OFS()<<"Number of blocks: "<<meta_->blocks_cnt_<<endl;
  logfileptr->OFS()<<"log2(Number of blocks): "<<log2(meta_->blocks_cnt_)<<endl;
  idxToBlk_->Dump();
#endif
}

template<typename T, class Allocator>
inline Int SuperNode<T,Allocator>::Shrink(){
  if(meta_->b_own_storage_){
    //TODO make sure that we do not have any extra space anywhere.

    //if there is too much room for either nzval or blocks, contract
    Int block_space = (Int)(blocks_+1 - (NZBlockDesc*)(meta_ +1)) - meta_->blocks_cnt_;
    Int nzval_space = (Int)((T*)meta_ - nzval_) - meta_->nzval_cnt_;

    if(block_space >0 || nzval_space >0){

      size_t new_size = storage_size_ - nzval_space*sizeof(T) - block_space*sizeof(NZBlockDesc);

#if 1
      size_t offset_meta = meta_->nzval_cnt_*sizeof(T);
      size_t offset_block = meta_->nzval_cnt_*sizeof(T)+sizeof(SuperNodeDesc);
#else
      size_t offset_meta = (char*)meta_ - (char*)nzval_;
      size_t offset_block = (char*)blocks_ - (char*)nzval_;
#endif


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
      std::copy(nzval_,nzval_+meta_->nzval_cnt_,(T*)locTmpPtr);

      //copy meta
      std::copy(meta_,meta_+1,(SuperNodeDesc*)(locTmpPtr+offset_meta));
      //copy blocks
      std::copy(blocks_+1-meta_->blocks_cnt_,blocks_+1,(NZBlockDesc*)(locTmpPtr+offset_block));

      //upcxx::deallocate(storage_container_);
      //storage_container_=tmpPtr;
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




    //        blocks_ = &blocks_container_.front();
    //        nzval_container_.resize(nzval_cnt_);
    //        nzval_ = &nzval_container_.front();
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
//    for(tgt_fc;tgt_fc<tgt_lc;tgt_fc++){
//      first_pivot_idx = src_snode.FindBlockIdx(tgt_fc);
//      if(first_pivot_idx>=0){
//        break;
//      }
//    }
      do {first_pivot_idx = src_snode->FindBlockIdx(tgt_fc); tgt_fc++;}
      while(first_pivot_idx<0 && tgt_fc<=tgt_lc);
      tgt_fc--;
#else
      Int closestR = -1;
      Int closestL = -1;
      do{
        first_pivot_idx = src_snode->FindBlockIdx(tgt_fc,closestR,closestL);
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
    logfileptr->OFS()<<"LOCK 1: first_pivot_idx<0"<<endl;
    gdb_lock();
  }
  assert(first_pivot_idx>=0);
}

template<typename T, class Allocator>
inline void SuperNode<T,Allocator>::FindUpdatedLastCol(SuperNode<T,Allocator> * src_snode, Int tgt_fc, Int first_pivot_idx , Int & tgt_lc,  Int & last_pivot_idx){
  //gdb_lock();
  scope_timer(a,UPDATE_SNODE_FIND_INDEX);
#ifdef _LINEAR_SEARCH_FCLC_
  if(src_snode->ITreeInitialized()){
#endif 
    //find the last row updated by src_snode
    tgt_lc = LastCol();
    last_pivot_idx = -1;
    //find the pivot idx
#ifndef _BINARY_BLOCK_SEARCH_
//    for(tgt_lc;tgt_lc>=tgt_fc;tgt_lc--){
//      last_pivot_idx = src_snode->FindBlockIdx(tgt_lc);
//      if(last_pivot_idx>=0){
//        break;
//      }
//    }
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

//    do {last_pivot_idx = src_snode->FindBlockIdx(tgt_lc); tgt_lc--;  logfileptr->OFS()<<"lpi: "<<last_pivot_idx<<" | "<<tgt_lc+1<<endl; }
//    while(last_pivot_idx<0 && tgt_lc>=tgt_fc);
//    tgt_lc++;
//    Int bak_tgt_lc = tgt_lc;
//    tgt_lc = LastCol();
//    Int tgt_cur;
//    Int count = tgt_lc - tgt_fc+1;
//    Int step;
//    //last_pivot_idx = src_snode->FindBlockIdx(tgt_lc);
//    //if(last_pivot_idx<0)
//    {
//      while (count > 0) {
//        step = count / 2;
//        tgt_cur = tgt_lc - step;
//        last_pivot_idx = src_snode->FindBlockIdx(tgt_cur); 
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
//    logfileptr->OFS()<<src_snode->FindBlockIdx(tgt_lc+1)<<endl;
//    logfileptr->OFS()<<src_snode->FindBlockIdx(tgt_lc)<<endl;
//    assert(src_snode->FindBlockIdx(tgt_lc+1)<0 && src_snode->FindBlockIdx(tgt_lc)>=0);
//    assert(bak_tgt_lc == tgt_lc);
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

//#ifdef COMPACT_AGGREGATES
template<typename T, class Allocator>
inline Int SuperNode<T,Allocator>::Merge(SuperNode<T,Allocator> * src_snode, SnodeUpdate &update){
  scope_timer_special(a,MERGE_SNODE);

  assert(meta_->b_own_storage_);

  Int src_snode_size = src_snode->Size();
  Int tgt_snode_size = Size();

  Int & pivot_idx = update.blkidx;
  Int & pivot_fr = update.src_first_row;

  Int tgt_fc;
  Int first_pivot_idx;
  FindUpdatedFirstCol(src_snode, 0, tgt_fc, first_pivot_idx);
  NZBlockDesc & first_pivot_desc = src_snode->GetNZBlockDesc(first_pivot_idx);

  //parse src_snode
  ITree::Interval overlap;
  ITree::Interval curInter;
  ITree::Interval newInter;
  std::queue< ITree::Interval > toInsert;
  for(Int blkidx = first_pivot_idx; blkidx<src_snode->NZBlockCnt(); ++blkidx){
    NZBlockDesc & blk_desc = src_snode->GetNZBlockDesc(blkidx);
    Int fr = max(FirstCol(),blk_desc.GIndex);
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

//#ifdef COMPACT_AGGREGATES
template<typename T, class Allocator>
inline Int SuperNode<T,Allocator>::Aggregate(SuperNode<T,Allocator> * src_snode){

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

#if 0
  for(Int blkidx = first_pivot_idx; blkidx<src_snode->NZBlockCnt(); ++blkidx){
    NZBlockDesc & blk_desc = src_snode->GetNZBlockDesc(blkidx);
    Int nrows = src_snode->NRows(blkidx);
    for(Int rowidx = 0; rowidx<nrows; ++rowidx){
      Int row = blk_desc.GIndex + rowidx;

      if(row>=tgt_fc){
        Int src_offset = blk_desc.Offset + (row - blk_desc.GIndex)*src_snode_size;

        Int tgt_blkidx = FindBlockIdx(row);
        if(tgt_blkidx==-1){
          logfileptr->OFS()<<"LOCK 2: tgt_blkidx=-1"<<endl;
          gdb_lock();
        }
        NZBlockDesc & tgt_desc = GetNZBlockDesc(tgt_blkidx);
        Int tgt_offset = tgt_desc.Offset + (row - tgt_desc.GIndex)*tgt_snode_size;

        T * src = src_snode->GetNZval(src_offset);
        T * tgt = GetNZval(tgt_offset);

        //blas::Axpy(tgt_snode_size,ONE<T>(),src,1,tgt,1);
        #pragma omp simd
        for(Int i = 0; i< tgt_snode_size;i+=1){ tgt[i] += src[i]; }

      }
    }
  }
#else
  Int tgt_blkcnt = NZBlockCnt();
  Int tgt_blkidx = 0;
  NZBlockDesc tgt_desc = GetNZBlockDesc(tgt_blkidx);
  Int tgt_nrows = NRows(tgt_blkidx);

  for(Int blkidx = first_pivot_idx; blkidx<src_snode->NZBlockCnt(); ++blkidx){
    NZBlockDesc & blk_desc = src_snode->GetNZBlockDesc(blkidx);
    Int nrows = src_snode->NRows(blkidx);
    for(Int rowidx = 0; rowidx<nrows; ++rowidx){
      Int row = blk_desc.GIndex + rowidx;

      if(row>=tgt_fc){
        Int src_offset = blk_desc.Offset + (row - blk_desc.GIndex)*src_snode_size;

        while(row>tgt_desc.GIndex + tgt_nrows-1){
          tgt_blkidx++;
          tgt_desc = GetNZBlockDesc(tgt_blkidx);
          tgt_nrows = NRows(tgt_blkidx);
          if(tgt_blkidx==tgt_blkcnt){
            tgt_blkidx=-1;
            break;
          }
        }

        if(tgt_blkidx==-1){
          logfileptr->OFS()<<"LOCK 2: tgt_blkidx=-1"<<endl;
          gdb_lock();
        }
        Int tgt_offset = tgt_desc.Offset + (row - tgt_desc.GIndex)*tgt_snode_size;

        T * src = src_snode->GetNZval(src_offset);
        T * tgt = GetNZval(tgt_offset);

        //blas::Axpy(tgt_snode_size,ONE<T>(),src,1,tgt,1);
        #pragma omp simd
        for(Int i = 0; i< tgt_snode_size;i+=1){ tgt[i] += src[i]; }

      }
    }
  }
#endif


  return 0;
}


//#ifdef COMPACT_AGGREGATES
//CHECKED ON 11-18-2014
template<typename T, class Allocator>
inline Int SuperNode<T,Allocator>::UpdateAggregate(SuperNode<T,Allocator> * src_snode, SnodeUpdate &update, 
    TempUpdateBuffers<T> & tmpBuffers, Int iTarget){

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
    T * buf = NULL;
    T beta = ZERO<T>();
#ifdef _DEBUG_
    tmpBuffers.tmpBuf.Resize(tgt_width,src_nrows);
#endif

    buf = &tmpBuffers.tmpBuf[0];

    //everything is in row-major
    SYMPACK_TIMER_START(UPDATE_SNODE_GEMM);
    blas::Gemm('T','N',tgt_width, src_nrows,src_snode_size,
        MINUS_ONE<T>(),pivot,src_snode_size,
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
      tmpBuffers.src_colindx.resize(tgt_width);
      tmpBuffers.src_to_tgt_offset.resize(src_nrows);
      Int colidx = 0;
      Int rowidx = 0;
      Int offset = 0;

      Int src_blkcnt = src_snode->NZBlockCnt();
#if 0
      for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
        NZBlockDesc & cur_block_desc = src_snode->GetNZBlockDesc(blkidx);
        Int cur_src_nrows = src_snode->NRows(blkidx);
        Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
        Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
        cur_src_nrows = cur_src_lr - cur_src_fr +1;

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
#else

  Int tgt_blkcnt = NZBlockCnt();
  Int tgt_blkidx = 0;
  NZBlockDesc tgt_desc = GetNZBlockDesc(tgt_blkidx);
  Int tgt_nrows = NRows(tgt_blkidx);

      for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
        NZBlockDesc & cur_block_desc = src_snode->GetNZBlockDesc(blkidx);
        Int cur_src_nrows = src_snode->NRows(blkidx);
        Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
        Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
        cur_src_nrows = cur_src_lr - cur_src_fr +1;

        Int row = cur_src_fr;
        while(row<=cur_src_lr){


        while(row>tgt_desc.GIndex + tgt_nrows-1){
          tgt_blkidx++;
          tgt_desc = GetNZBlockDesc(tgt_blkidx);
          tgt_nrows = NRows(tgt_blkidx);
          if(tgt_blkidx==tgt_blkcnt){
            tgt_blkidx=-1;
            break;
          }
        }

          assert(tgt_blkidx>=0);
          Int lr = min(cur_src_lr,tgt_desc.GIndex + tgt_nrows-1);
          Int tgtOffset = tgt_desc.Offset 
            + (row - tgt_desc.GIndex)*tgt_snode_size;
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

#endif
      SYMPACK_TIMER_STOP(UPDATE_SNODE_INDEX_MAP);


      //Multiple cases to consider
      if(first_pivot_idx==last_pivot_idx){
        // Updating contiguous columns
        Int tgt_offset = (tgt_fc - FirstCol());
        for(Int rowidx = 0; rowidx < src_nrows; ++rowidx){
          T * A = &buf[rowidx*tgt_width];
          T * B = &tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_offset];
          #pragma omp simd
          for(Int i = 0; i < tgt_width; ++i){ B[i] += A[i]; }
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
inline Int SuperNode<T,Allocator>::Update(SuperNode<T,Allocator> * src_snode, SnodeUpdate &update, 
    TempUpdateBuffers<T> & tmpBuffers){

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
    buf = &tmpBuffers.tmpBuf[0];
  }

  //everything is in row-major
  SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_GEMM);
logfileptr->OFS()<<"GEMM ("<<tgt_width<<"-by-"<<src_snode_size<<") x ("<<src_snode_size<<"-by-"<<src_nrows<<")"<<endl;
  blas::Gemm('T','N',tgt_width, src_nrows,src_snode_size,
      MINUS_ONE<T>(),pivot,src_snode_size,
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
        Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);

        Int row = cur_src_fr;
        while(row<=cur_src_lr){
          Int tgt_blk_idx = FindBlockIdx(row);
            assert(tgt_blk_idx>=0);
            NZBlockDesc & cur_tgt_desc = GetNZBlockDesc(tgt_blk_idx);
            Int lr = min(cur_src_lr,cur_tgt_desc.GIndex + NRows(tgt_blk_idx)-1);
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
#if 0
      for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
        NZBlockDesc & cur_block_desc = src_snode->GetNZBlockDesc(blkidx);
        Int cur_src_nrows = src_snode->NRows(blkidx);
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
#else

  Int tgt_blkcnt = NZBlockCnt();
  Int tgt_blkidx = 0;
  NZBlockDesc tgt_desc = GetNZBlockDesc(tgt_blkidx);
  Int tgt_nrows = NRows(tgt_blkidx);

      for(Int blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
        NZBlockDesc & cur_block_desc = src_snode->GetNZBlockDesc(blkidx);
        Int cur_src_nrows = src_snode->NRows(blkidx);
        Int cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
        Int cur_src_fr = max(tgt_fc, cur_block_desc.GIndex);
        cur_src_nrows = cur_src_lr - cur_src_fr +1;

        //The other one MUST reside into a single block in the target
        Int row = cur_src_fr;
        while(row<=cur_src_lr){

        while(row>tgt_desc.GIndex + tgt_nrows-1){
          tgt_blkidx++;
          tgt_desc = GetNZBlockDesc(tgt_blkidx);
          tgt_nrows = NRows(tgt_blkidx);
          if(tgt_blkidx==tgt_blkcnt){
            tgt_blkidx=-1;
            break;
          }
        }


            assert(tgt_blkidx>=0);
            Int lr = min(cur_src_lr,tgt_desc.GIndex + tgt_nrows-1);
            Int tgtOffset = tgt_desc.Offset 
              + (row - tgt_desc.GIndex)*tgt_snode_size;
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
#endif
      //Multiple cases to consider
      SYMPACK_TIMER_SPECIAL_STOP(UPDATE_SNODE_INDEX_MAP);

      SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_ADD);
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
      SYMPACK_TIMER_SPECIAL_STOP(UPDATE_SNODE_ADD);
    }
  }
  return 0;
}


//CHECKED ON 11-18-2014
template<typename T, class Allocator>
  inline Int SuperNode<T,Allocator>::Factorize(TempUpdateBuffers<T> & tmpBuffers){
#if defined(_NO_COMPUTATION_)
    return 0;
#endif


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

  }





template<typename T, class Allocator>
bool SuperNode<T,Allocator>::FindNextUpdate(SnodeUpdate & nextUpdate, const SYMPACK::vector<Int> & Xsuper,  const SYMPACK::vector<Int> & SupMembership, bool isLocal){
  scope_timer_special(a,FIND_NEXT_UPDATE);
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


template<typename T, class Allocator>
  inline ITree * SuperNode<T,Allocator>::CreateITree(){
#if defined(_AVL_ITREE_)
    return new AVLITree();
#elif defined(_DSW_ITREE_)
    return new DSWITree();
#else
    return new ITree();
#endif
  } 

template<typename T, class Allocator>
inline void SuperNode<T,Allocator>::InitIdxToBlk(){
#ifndef ITREE
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

    ITree::Interval cur_interv = { cur_fr, cur_lr, blkidx};
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
    SYMPACK::Serialize(buffer,(char*)nzval_ptr,size);
    //now we need to modify the first block data
    char * tail = buffer.front()+size; 
    NZBlockDesc* new_blk_ptr= (NZBlockDesc*)(tail-sizeof(NZBlockDesc));

    Int offset = new_blk_ptr->Offset + local_first_row*this->Size();

    if(offset!=0){
      Int blkCnt = 1;
      if(blkCnt<nzblk_cnt){
        NZBlockDesc * curBlockPtr = NULL;
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
    inline void SuperNode<T,Allocator>::forward_update(SuperNode<T,Allocator> * src_contrib,Int iOwner){
      SuperNode<T,Allocator> * tgt_contrib = this;

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
          //        assert(tgt_blkidx==0);
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
    inline void SuperNode<T,Allocator>::back_update(SuperNode<T,Allocator> * src_contrib){
      SuperNode<T,Allocator> * tgt_contrib = this;
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
            //          assert(src_nzblk_idx==0);
            src_nzblk_idx++;
          }

        } while(tgt_lr>src_lr);
      }
    }


template <typename T, class Allocator> 
    inline void SuperNode<T,Allocator>::back_update_contrib(SuperNode<T> * cur_snode){
        SuperNode<T,Allocator> * contrib = this;
        Int nrhs = this->Size();

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
}


template <typename T, class Allocator> 
inline void SuperNode<T,Allocator>::forward_update_contrib( T * RHS, SuperNode<T> * cur_snode, SYMPACK::vector<Int> & perm){
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
          diag_nzval[kk*nrhs+j] = (RHS[srcRow-1 + j*n] + diag_nzval[kk*nrhs+j]) / chol_nzval[kk*cur_snode->Size()+kk];
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
        NZBlockDesc * curBlockPtr = NULL;
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


inline std::ostream& operator<<( std::ostream& os,  NZBlockDesc& block){
  os<<"GIndex = "<<block.GIndex<<std::endl;
  os<<"Offset = "<<block.Offset<<std::endl;
  os<<"Last   = "<<(block.Last?string("true"):string("false"))<<std::endl;
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




//  template <typename T> inline void Serialize(Icomm & buffer,SuperNode<T,Allocator> & snode){
//    Int nzblk_cnt = snode.NZBlockCnt();
//    Int nzval_cnt_ = snode.Size()*snode.NRowsBelowBlock(0);
//    T* nzval_ptr = snode.GetNZval(0);
//    NZBlockDesc* nzblk_ptr = &snode.GetNZBlockDesc(0);
//    Int snode_id = snode.Id();
//    Int snode_fc = snode.FirstCol();
//    Int snode_lc = snode.LastCol();
//
//    buffer.clear();
//    buffer.resize(6*sizeof(Int) + nzblk_cnt*sizeof(NZBlockDesc) + nzval_cnt_*sizeof(T));
//
//    buffer<<snode_id;
//    buffer<<snode_fc;
//    buffer<<snode_lc;
//    buffer<<snode.N();
//    buffer<<nzblk_cnt;
//    Serialize(buffer,nzblk_ptr,nzblk_cnt);
//    buffer<<nzval_cnt_;
//    Serialize(buffer,nzval_ptr,nzval_cnt_);
//  }
//
//  template <typename T> inline void Serialize(Icomm & buffer,SuperNode<T> & snode, Int first_blkidx, Int first_row){
//    NZBlockDesc* nzblk_ptr = &snode.GetNZBlockDesc(first_blkidx);
//    Int local_first_row = first_row - nzblk_ptr->GIndex;
//    Int nzblk_cnt = snode.NZBlockCnt() - first_blkidx;
//    Int nzval_cnt_ = snode.Size()*(snode.NRowsBelowBlock(first_blkidx)-local_first_row);
//    T* nzval_ptr = snode.GetNZval(nzblk_ptr->Offset) + local_first_row*snode.Size();
//    Int snode_id = snode.Id();
//    Int snode_fc = snode.FirstCol();
//    Int snode_lc = snode.LastCol();
//
//    buffer.clear();
//    buffer.resize(6*sizeof(Int) + nzblk_cnt*sizeof(NZBlockDesc) + nzval_cnt_*sizeof(T));
//
//    buffer<<snode_id;
//    buffer<<snode_fc;
//    buffer<<snode_lc;
//    buffer<<snode.N();
//    buffer<<nzblk_cnt;
//    NZBlockDesc * new_nzblk_ptr = reinterpret_cast<NZBlockDesc *>(buffer.back());
//    Serialize(buffer,nzblk_ptr,nzblk_cnt);
//    //replace the GIndex of the serialized block descriptor
//    // by the new first row, and update the offset appropriately
//    new_nzblk_ptr->GIndex = first_row;
//    new_nzblk_ptr->Offset = nzblk_ptr->Offset + local_first_row*snode.Size();
//
//    buffer<<nzval_cnt_;
//    Serialize(buffer,nzval_ptr,nzval_cnt_);
//  }
//
//
//  template <typename T> inline void Serialize(Icomm & buffer,SuperNode<T> & snode, Int first_blkidx, Int first_row, size_t extra_bytespace){
//    NZBlockDesc* nzblk_ptr = &snode.GetNZBlockDesc(first_blkidx);
//    Int local_first_row = first_row - nzblk_ptr->GIndex;
//    Int nzblk_cnt = snode.NZBlockCnt() - first_blkidx;
//    Int nzval_cnt_ = snode.Size()*(snode.NRowsBelowBlock(first_blkidx)-local_first_row);
//    T* nzval_ptr = snode.GetNZval(nzblk_ptr->Offset) + local_first_row*snode.Size();
//    Int snode_id = snode.Id();
//    Int snode_fc = snode.FirstCol();
//    Int snode_lc = snode.LastCol();
//
//    buffer.clear();
//    buffer.resize(6*sizeof(Int) + nzblk_cnt*sizeof(NZBlockDesc) + nzval_cnt_*sizeof(T) + extra_bytespace);
//
//    buffer<<snode_id;
//    buffer<<snode_fc;
//    buffer<<snode_lc;
//    buffer<<snode.N();
//    buffer<<nzblk_cnt;
//    NZBlockDesc * new_nzblk_ptr = reinterpret_cast<NZBlockDesc *>(buffer.back());
//    Serialize(buffer,nzblk_ptr,nzblk_cnt);
//    //replace the GIndex of the serialized block descriptor
//    // by the new first row, and update the offset appropriately
//    new_nzblk_ptr->GIndex = first_row;
//    new_nzblk_ptr->Offset = nzblk_ptr->Offset + local_first_row*snode.Size();
//
//    buffer<<nzval_cnt_;
//    Serialize(buffer,nzval_ptr,nzval_cnt_);
//  }




  template <typename T,class Allocator> inline size_t Deserialize(char * buffer, size_t size, SuperNode<T,Allocator> & snode){
    snode.Init(&buffer[0],size);
    snode.InitIdxToBlk();
    return 0;
  }



} // namespace SYMPACK





#endif // _SUPERNODE_IMPL_HPP_
