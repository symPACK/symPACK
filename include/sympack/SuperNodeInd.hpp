/*
	 Copyright (c) 2016 The Regents of the University of California,
	 through Lawrence Berkeley National Laboratory.  

   Author: Mathias Jacquelin
	 
   This file is part of symPACK. All rights reserved.

	 Redistribution and use in source and binary forms, with or without
	 modification, are permitted provided that the following conditions are met:

	 (1) Redistributions of source code must retain the above copyright notice, this
	 list of conditions and the following disclaimer.
	 (2) Redistributions in binary form must reproduce the above copyright notice,
	 this list of conditions and the following disclaimer in the documentation
	 and/or other materials provided with the distribution.
	 (3) Neither the name of the University of California, Lawrence Berkeley
	 National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
	 be used to endorse or promote products derived from this software without
	 specific prior written permission.

	 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
	 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
	 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
	 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
	 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

	 You are under no obligation whatsoever to provide any bug fixes, patches, or
	 upgrades to the features, functionality or performance of the source code
	 ("Enhancements") to anyone; however, if you choose to make your Enhancements
	 available either publicly, or directly to Lawrence Berkeley National
	 Laboratory, without imposing a separate written license agreement for such
	 Enhancements, then you hereby grant the following license: a non-exclusive,
	 royalty-free perpetual license to install, use, modify, prepare derivative
	 works, incorporate into other computer software, distribute, and sublicense
	 such enhancements or derivative works thereof, in binary and source code form.
*/
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
  SuperNodeInd(Int aiId, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices);
  SuperNodeInd(Int aiId, Int aiFc, Int aiLc, Int ai_num_rows, Int aiN, Int aiNZBlkCnt=-1);
  SuperNodeInd(Int aiId, Int aiFc, Int aiLc, Int aiN);
  SuperNodeInd(char * storage_ptr,size_t storage_size, Int GIndex = -1);

    
  virtual void Init(char * storage_ptr,size_t storage_size, Int GIndex = -1);
  virtual void Init(Int aiId, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices);

  virtual inline void AddNZBlock(Int aiNRows, Int aiNCols, Int aiGIndex);

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
