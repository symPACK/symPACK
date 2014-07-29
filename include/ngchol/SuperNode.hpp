/// @file SuperNode.hpp
/// @brief TODO
/// @author Mathias Jacquelin
/// @date 2013-09-10
#ifndef _SUPERNODE_FACT_HPP_
#define _SUPERNODE_FACT_HPP_

#include "ngchol/Environment.hpp"
#include "ngchol/NumVec.hpp"
#include "ngchol/IntervalTree.hpp"
#include "ngchol/CommTypes.hpp"

#include <list>

namespace LIBCHOLESKY{

struct NZBlockDesc{
    Int GIndex;
    size_t Offset;
    NZBlockDesc():GIndex(-1),Offset(-1){};
    NZBlockDesc(Int aGIndex, size_t aOffset):GIndex(aGIndex),Offset(aOffset){};
};



template<typename T>
class SuperNode{

  protected:

  ITree * idxToBlk_;

  std::vector<T> nzval_container_;
  Int nzval_cnt_;
  T * nzval_;
  std::vector<NZBlockDesc> blocks_container_;
  NZBlockDesc * blocks_;
  Int blocks_cnt_;


  bool b_own_storage_;
  Int iId_;
  Int iSize_;
  Int iFirstCol_;
  Int iLastCol_;

  public:


  inline Int & Id(){ return iId_;}
  inline Int FirstCol(){ return iFirstCol_;}
  inline Int LastCol(){ return iLastCol_;}
  inline Int Size(){ return iSize_;}
  inline Int NZBlockCnt(){ return blocks_cnt_;}

  inline NZBlockDesc & GetNZBlockDesc(Int aiLocIndex){ return blocks_[aiLocIndex];}
  inline T* GetNZval(size_t offset){ return &nzval_[offset];}

  inline Int NRows(Int blkidx){
      NZBlockDesc & desc = GetNZBlockDesc(blkidx);
      size_t end = (blkidx<NZBlockCnt()-1)?GetNZBlockDesc(blkidx+1).Offset:nzval_cnt_;
      Int nRows = (end-desc.Offset)/Size();
      return nRows;
  }

  inline Int NRowsBelowBlock(Int blkidx){
      NZBlockDesc & desc = GetNZBlockDesc(blkidx);
      size_t end = nzval_cnt_;
      Int nRows = (end-desc.Offset)/Size();
      return nRows;
  }


  inline Int StorageSize(){ return nzval_cnt_*sizeof(T) + blocks_cnt_*sizeof(NZBlockDesc);}

  SuperNode() :iId_(-1),iFirstCol_(-1),iLastCol_(-1) {
  }

  SuperNode(Int aiId, Int aiFc, Int aiLc, Int ai_num_rows) :iId_(aiId),iFirstCol_(aiFc),iLastCol_(aiLc) {
     //this is an upper bound
     assert(ai_num_rows>0);

     //compute supernode size / width
     iSize_ = iLastCol_ - iFirstCol_+1;

     nzval_cnt_ = 0;
     nzval_container_.reserve(iSize_*ai_num_rows);
     nzval_ = NULL;

     blocks_container_.reserve(ai_num_rows);
     blocks_ = NULL;
     blocks_cnt_ = 0;

     idxToBlk_ = new ITree();

     b_own_storage_ = true;


  }; 



  SuperNode(Int aiId, Int aiFc, Int aiLc, IntNumVec & xlindx, IntNumVec & lindx) :iId_(aiId),iFirstCol_(aiFc),iLastCol_(aiLc) {

    //compute supernode size / width
    iSize_ = iLastCol_ - iFirstCol_+1;
     idxToBlk_ = new ITree();
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

  }





  void Init(Int aiId, Int aiFc, Int aiLc, NZBlockDesc * a_block_desc, Int a_desc_cnt,
              T * a_nzval, Int a_nzval_cnt) {
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
     idxToBlk_ = new ITree();

    //restore 0-based offsets and compute global_to_local indices
    for(Int blkidx=blocks_cnt_-1; blkidx>=0;--blkidx){
      blocks_[blkidx].Offset -= blocks_[0].Offset;
    }

    for(Int blkidx=0; blkidx<blocks_cnt_;++blkidx){
      Int cur_fr = blocks_[blkidx].GIndex;
      Int cur_lr = cur_fr + NRows(blkidx) -1;

      ITree::Interval cur_interv = { cur_fr, cur_lr, blkidx};
      idxToBlk_->Insert(cur_interv);
    }


    
  }


  SuperNode(Int aiId, Int aiFc, Int aiLc, NZBlockDesc * a_block_desc, Int a_desc_cnt,
              T * a_nzval, Int a_nzval_cnt) {
     Init(aiId, aiFc, aiLc, a_block_desc, a_desc_cnt, a_nzval, a_nzval_cnt);
  }



  ~SuperNode(){
  delete idxToBlk_;
 
  }
    

 

  inline void AddNZBlock(Int aiNRows, Int aiNCols, Int aiGIndex){

    //Resize the container if I own the storage
    if(b_own_storage_){
      Int cur_fr = aiGIndex;
      Int cur_lr = cur_fr + aiNRows -1;

      ITree::Interval cur_interv = { cur_fr, cur_lr, blocks_cnt_};
      idxToBlk_->Insert(cur_interv);

      blocks_container_.push_back(NZBlockDesc(aiGIndex,nzval_cnt_));

      blocks_cnt_++;
      nzval_cnt_+=aiNRows*iSize_;
      nzval_container_.resize(nzval_cnt_,ZERO<T>());
      nzval_ = &nzval_container_.front();
      blocks_ = &blocks_container_.front();
    }
  }
 

  Int FindBlockIdx(Int aiGIndex){
//      ITree::Interval it = {aiGIndex, aiGIndex,0};
//      ITree::Interval * res = idxToBlk_->IntervalSearch(it);
      ITree::Interval * res = idxToBlk_->IntervalSearch(aiGIndex,aiGIndex);
      if (res == NULL){
         return -1;
      }
      else{
         return res->block_idx;
      }
  }

  void DumpITree(){
    logfileptr->OFS()<<"Number of blocks: "<<blocks_cnt_<<endl;
    logfileptr->OFS()<<"log2(Number of blocks): "<<log2(blocks_cnt_)<<endl;
    idxToBlk_->Dump();
  }

  Int Shrink(){
    if(b_own_storage_){
      //blocks_container_.shrink_to_fit();
      blocks_ = &blocks_container_.front();

      nzval_container_.resize(nzval_cnt_);
      nzval_ = &nzval_container_.front();
    }

    return StorageSize();
  }



};




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




template <typename T> inline void Deserialize(char * buffer, SuperNode<T> & snode){
            Int snode_id = *(Int*)&buffer[0];
            Int snode_fc = *(((Int*)&buffer[0])+1);
            Int snode_lc = *(((Int*)&buffer[0])+2);
            Int nzblk_cnt = *(((Int*)&buffer[0])+3);
            NZBlockDesc * blocks_ptr = 
              reinterpret_cast<NZBlockDesc*>(&buffer[4*sizeof(Int)]);
            Int nzval_cnt = *(Int*)(blocks_ptr + nzblk_cnt);
            T * nzval_ptr = (T*)((Int*)(blocks_ptr + nzblk_cnt)+1);

            //Create the dummy supernode for that data
            snode.Init(snode_id,snode_fc,snode_lc, blocks_ptr, nzblk_cnt, nzval_ptr, nzval_cnt);
  }




} // namespace LIBCHOLESKY


#endif // _SUPERNODE_FACT_HPP_
