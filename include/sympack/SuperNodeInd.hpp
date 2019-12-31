#ifndef _SUPERNODE_IND_FACT_HPP_
#define _SUPERNODE_IND_FACT_HPP_

#include "sympack/Environment.hpp"
#include "sympack/blas.hpp"
#include "sympack/lapack.hpp"
#include "sympack/IntervalTree.hpp"
#include "sympack/CommTypes.hpp"
#include "sympack/SuperNode.hpp"

#ifdef NO_INTRA_PROFILE
#if defined (SPROFILE)
#define SYMPACK_TIMER_START(a) 
#define SYMPACK_TIMER_STOP(a) 
#define scope_timer(b,a)
#endif
#endif


namespace symPACK{


////////////////////////////////////////
/// Class representing a supernode.
/// Class representing a supernode stored as a collection of 
/// blocks of contiguous rows in row-major format.
/////////////////////////////////////////
template<typename T, class Allocator = UpcxxAllocator>
class SuperNodeInd: public SuperNode<T,Allocator>{
  public:

  protected:
  T * diag_;
  
  public:
  inline T* GetDiag(){ return diag_;}

  SuperNodeInd();
  SuperNodeInd(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices, Int panel=-1);
  SuperNodeInd(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int ai_num_rows, Int aiN, Int aiNZBlkCnt=-1, Int panel=-1);
  SuperNodeInd(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int aiN, Int panel=-1);
  SuperNodeInd(char * storage_ptr,size_t storage_size, Int GIndex = -1);

    
  virtual void Init(char * storage_ptr,size_t storage_size, Int GIndex = -1);
  virtual void Init(Int aiId, Int aiFr, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices, Int panel=-1);

  virtual inline void AddNZBlock(Int aiNRows, Int aiGIndex);

  virtual inline Int Shrink();
  
  //Factorize the supernode
  virtual inline Int Factorize(TempUpdateBuffers<T> & tmpBuffers);
  virtual inline Int Factorize(SuperNode<T,Allocator> * diag_snode, TempUpdateBuffers<T> & tmpBuffers);

  
  virtual inline Int UpdateAggregate(SuperNode<T,Allocator> * src_snode, SnodeUpdate &update, 
              TempUpdateBuffers<T> & tmpBuffers,Int iTarget,Int iam);
  virtual inline Int Update(SuperNode<T,Allocator> * src_snode, SnodeUpdate &update, 
              TempUpdateBuffers<T> & tmpBuffers);

  //forward and backward solve phases
  virtual inline void forward_update_contrib( T * RHS, SuperNode<T> * cur_snode, std::vector<Int> & perm);
  virtual inline void forward_update_contrib( SuperNode<T> * cur_snode, Int nrhsOffset = 0,Int pnrhs=-1);
  virtual inline void back_update_contrib(SuperNode<T> * cur_snode, Int nrhsOffset = 0,Int pnrhs=-1);

  //clear function to reset the supernode (when redistributing another matrix for instance
  virtual inline void clear(){ 
    if(this->NNZ()>0){
      std::fill(this->nzval_,this->nzval_+this->NNZ(),T(0));
    }

    if(this->Size()>0){
      std::fill(this->diag_,this->diag_+this->Size(),T(0));
    }
  }


};



//template <typename T, class Allocator> inline void Serialize(Icomm & buffer,SuperNode<T,Allocator> & snode, Int first_blkidx=0, Idx first_row=0);
template <typename T, class Allocator> inline std::ostream& operator<<( std::ostream& os,  SuperNodeInd<T,Allocator>& snode);
template <typename T, class Allocator> inline size_t Deserialize(char * buffer, size_t size, SuperNodeInd<T,Allocator> & snode);




} // namespace SYMPACK


#include "sympack/impl/SuperNodeInd_impl.hpp"


#ifdef NO_INTRA_PROFILE
#if defined (SPROFILE)
#define SYMPACK_TIMER_START(a) SYMPACK_FSTART(a);
#define SYMPACK_TIMER_STOP(a) SYMPACK_FSTOP(a);
#define scope_timer(b,a) symPACK_scope_timer b(#a)
#endif
#endif



#endif // _SUPERNODE_IND_FACT_HPP_
