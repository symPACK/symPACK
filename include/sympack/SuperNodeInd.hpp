/// @file SuperNode.hpp
/// @brief TODO
/// @author Mathias Jacquelin
/// @date 2013-09-10
#ifndef _SUPERNODE_IND_FACT_HPP_
#define _SUPERNODE_IND_FACT_HPP_

#include "sympack/Environment.hpp"
#include "sympack/blas.hpp"
#include "sympack/lapack.hpp"
#include "sympack/IntervalTree.hpp"
#include "sympack/CommTypes.hpp"
#include "sympack/SuperNode.hpp"

#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) 
#define TIMER_STOP(a) 
#define scope_timer(b,a)
#endif
#endif


namespace SYMPACK{


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

  SuperNodeInd(Int aiId, Int aiFc, Int aiLc, Int ai_num_rows, Int aiN, Int aiNZBlkCnt=-1);
  SuperNodeInd(Int aiId, Int aiFc, Int aiLc, Int aiN);
  SuperNodeInd(char * storage_ptr,size_t storage_size, Int GIndex = -1);

    
  virtual void Init(char * storage_ptr,size_t storage_size, Int GIndex = -1);

  virtual inline void AddNZBlock(Int aiNRows, Int aiNCols, Int aiGIndex);

  virtual inline Int Shrink();
  
  //Factorize the supernode
  virtual inline Int Factorize(TempUpdateBuffers<T> & tmpBuffers);
  virtual inline Int UpdateAggregate(SuperNode<T,Allocator> * src_snode, SnodeUpdate &update, 
              TempUpdateBuffers<T> & tmpBuffers,Int iTarget);
  virtual inline Int Update(SuperNode<T,Allocator> * src_snode, SnodeUpdate &update, 
              TempUpdateBuffers<T> & tmpBuffers);

  //forward and backward solve phases
  virtual inline void forward_update_contrib( T * RHS, SuperNode<T> * cur_snode, SYMPACK::vector<Int> & perm);
  virtual inline void back_update_contrib(SuperNode<T> * cur_snode);

};



//template <typename T, class Allocator> inline void Serialize(Icomm & buffer,SuperNode<T,Allocator> & snode, Int first_blkidx=0, Idx first_row=0);
template <typename T, class Allocator> inline std::ostream& operator<<( std::ostream& os,  SuperNodeInd<T,Allocator>& snode);
template <typename T, class Allocator> inline size_t Deserialize(char * buffer, size_t size, SuperNodeInd<T,Allocator> & snode);




} // namespace SYMPACK


#include "sympack/SuperNodeInd_impl.hpp"


#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) TAU_FSTART(a);
#define TIMER_STOP(a) TAU_FSTOP(a);
#define scope_timer(b,a) CTF_scope_timer b(#a)
#endif
#endif



#endif // _SUPERNODE_IND_FACT_HPP_
