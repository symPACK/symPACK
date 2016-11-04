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
#ifndef _TYPES_DECL_HPP_
#define _TYPES_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/timer.hpp"
#include "sympack/CommTypes.hpp"

#include <string>

namespace symPACK{

  typedef vector<Ptr> PtrVec;
  typedef vector<Idx> IdxVec;


  struct RelaxationParameters{
    Int nrelax0;
    Int nrelax1;
    Int nrelax2;
    Int maxSize;

    double zrelax0 = 0.8;
    double zrelax1 = 0.1;
    double zrelax2 = 0.05;

    RelaxationParameters(){
      nrelax0 = 8 ;
      nrelax1 = 32;
      nrelax2 = 64;
    }
    RelaxationParameters(Int pmaxSize):RelaxationParameters(){
      maxSize = pmaxSize;
      if(maxSize>0){
        nrelax0 = std::min(nrelax0,maxSize);
        nrelax1 = std::min(nrelax1,maxSize);
        nrelax2 = std::min(nrelax2,maxSize);
      }

    }

    RelaxationParameters(Int pnrelax0, Int pnrelax1, Int pnrelax2, Int pmaxSize){
      maxSize = pmaxSize;
      if(maxSize>0){
        nrelax0 = std::min(pnrelax0,maxSize);
        nrelax1 = std::min(pnrelax1,maxSize);
        nrelax2 = std::min(pnrelax2,maxSize);
      }
      else{
        nrelax0 = pnrelax0;
        nrelax1 = pnrelax1;
        nrelax2 = pnrelax2;
      }

    }

    void SetMaxSize(Int pmaxSize){
      maxSize = pmaxSize;
      if(maxSize>0){
      nrelax0 = std::min(nrelax0,pmaxSize);
      nrelax1 = std::min(nrelax1,pmaxSize);
      nrelax2 = std::min(nrelax2,pmaxSize);
      }
    }
    void SetNrelax0(Int pnrelax0){
      if(maxSize>0){
        nrelax0 = std::min(pnrelax0,maxSize);
      }
      else{
        nrelax0 = pnrelax0;
      }
    }
    void SetNrelax1(Int pnrelax1){
      if(maxSize>0){
      nrelax1 = std::min(pnrelax1,maxSize);
      }
      else{
        nrelax1 = pnrelax1;
      }
    }
    void SetNrelax2(Int pnrelax2){
      if(maxSize>0){
      nrelax2 = std::min(pnrelax2,maxSize);
      }
      else{
        nrelax2 = pnrelax2;
      }
    }

  };

  enum DecompositionType {LL,LDL};

  enum MappingType {ROW2D,COL2D,MODWRAP2D,MODWRAP2DNS,WRAP2D,WRAP2DFORCED};
  enum FactorizationType {FANOUT,FANBOTH,FANBOTH_STATIC};
  enum LoadBalanceType {NOLB,NNZ,NCOLS,WORK,SUBCUBE,SUBCUBE_NNZ};
  enum OrderingType {NATURAL,RCM,MMD,AMD,NDBOX,NDGRID,SCOTCH,PTSCOTCH,METIS,PARMETIS,USER};
  enum SchedulerType {DL,MCT,PR,FIFO};
  //enum OrderRefinementType {BarryDL,MCT,PR,FIFO};
  class symPACKOptions{
    public:
      int NpOrdering;
      DecompositionType decomposition;
      MappingType mappingType;
      std::string mappingTypeStr;
      FactorizationType factorization;
      LoadBalanceType load_balance;
      std::string load_balance_str;
      std::string order_refinement_str;
      OrderingType ordering;
      std::string orderingStr;
      SchedulerType scheduler;
      Int maxIsend;
      Int maxIrecv;
      CommEnvironment * commEnv;
      RelaxationParameters relax;
      MPI_Comm MPIcomm;
      int verbose;
      int dumpPerm;
      int * perm;
    protected:
      bool isSqrtP(){
        bool val = false;
        if(mappingTypeStr == "MODWRAP2D" || mappingTypeStr == "MODWRAP2DNS" || mappingTypeStr == "WRAP2D" || mappingTypeStr == "WRAP2DFORCED"){
          val = true;
        }

        //        switch(mappingType){
        //          case MODWRAP2D: case MODWRAP2DNS: case WRAP2D: case WRAP2DFORCED:
        //            val = true;
        //            break;
        //          default:
        //            val = false;
        //            break;
        //        }
        return val;
      }

    public:
      symPACKOptions(){
        verbose = 0;
        NpOrdering = 0;
        decomposition = LL; 
        factorization = FANBOTH;
        mappingTypeStr = "ROW2D";
        load_balance_str = "SUBCUBE-FO";
        orderingStr = "MMD";
        order_refinement_str = "NONE";
        scheduler = DL;
        maxIsend = 0;
        maxIrecv=0;
        relax = RelaxationParameters(0);
        relax.SetMaxSize(150);

        commEnv = NULL;
        MPIcomm = MPI_COMM_NULL;

        perm = NULL;
        dumpPerm = 0;
//        mappingType = ROW2D;
///        ordering = MMD;
//        load_balance = SUBCUBE;
      }

      Int used_procs(Int np){

        if(isSqrtP()){
          Int nr = (Int)sqrt((double)np);
          Int nc = np/nr;
          //Int nc = nr;
          np= nr*nc;
        }
        return np;
      }

  };







  struct LocalUpdate{
    Int src_snode_id;
    Int src_nzblk_idx;
    Int src_first_row;
    LocalUpdate(Int snode_id,Int nzblk_idx,Int first_row):src_snode_id(snode_id),src_nzblk_idx(nzblk_idx),src_first_row(first_row){};
  };

  struct SnodeUpdate{
    Int src_snode_id;
    Int tgt_snode_id;
    Int src_first_row;
    Int src_next_row;
    Int blkidx;
    Int next_blkidx;
    //std::vector<bool> is_factor_sent;

    SnodeUpdate(){
      src_snode_id = 0;
      tgt_snode_id = 0;
      src_first_row = 0;
      src_next_row = 0;
      blkidx = 0;
      next_blkidx = 0;
      //is_factor_sent.assign(np,false);
    }
  };

  template<typename T>
    class TempUpdateBuffers{
      public:
        std::vector<T> tmpBuf;
        std::vector<Int> src_colindx;
        std::vector<Int> src_to_tgt_offset;

        void Resize(Int size, Int mw){
          if(size*mw > tmpBuf.size()){
            tmpBuf.resize(size*mw);
          }
          if(mw > src_colindx.size()){
            src_colindx.resize(mw);
          }
          if(size > src_to_tgt_offset.size()){
            src_to_tgt_offset.resize(size);
          }
        }

        void Clear(){
          tmpBuf.clear();
          src_colindx.clear();
          src_to_tgt_offset.clear();
        }


        TempUpdateBuffers(){
        }
        TempUpdateBuffers(Int size, Int mw){
          Resize(size,mw);
        }
    };

  struct duet{
            Idx row;
            Idx col;
          };

struct sortDuet {
  bool operator() (const duet & a, const duet & b){
    bool retval = a.row<b.row;
    if(a.row==b.row){
      retval = a.col<b.col;
    }
    return retval;
  }
};

struct sortDuetInv {
  bool operator() (const duet & a, const duet & b){
    bool retval = a.row>b.row;
    if(a.row==b.row){
      retval = a.col>b.col;
    }
    return retval;
  }
};





  template<typename T>
  struct triplet{
            Idx row;
            Idx col;
            T val;
          };

template<typename T>
struct sortTriplet {
  bool operator() (const triplet<T> & a, const triplet<T> & b){
    bool retval = a.row<b.row;
    if(a.row==b.row){
      retval = a.col<b.col;
    }
    return retval;
  }
};

template<typename T>
struct sortTripletInv {
  bool operator() (const triplet<T> & a, const triplet<T> & b){
    bool retval = a.row>b.row;
    if(a.row==b.row){
      retval = a.col>b.col;
    }
    return retval;
  }
};




}

#endif //_TYPES_DECL_HPP_

