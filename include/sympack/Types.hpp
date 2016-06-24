#ifndef _TYPES_DECL_HPP_
#define _TYPES_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/timer.hpp"
//#include "sympack/NumMat.hpp"
#include "sympack/CommTypes.hpp"

#include <string>

namespace SYMPACK{

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
        nrelax0 = min(nrelax0,maxSize);
        nrelax1 = min(nrelax1,maxSize);
        nrelax2 = min(nrelax2,maxSize);
      }

    }

    RelaxationParameters(Int pnrelax0, Int pnrelax1, Int pnrelax2, Int pmaxSize){
      maxSize = pmaxSize;
      if(maxSize>0){
        nrelax0 = min(pnrelax0,maxSize);
        nrelax1 = min(pnrelax1,maxSize);
        nrelax2 = min(pnrelax2,maxSize);
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
      nrelax0 = min(nrelax0,pmaxSize);
      nrelax1 = min(nrelax1,pmaxSize);
      nrelax2 = min(nrelax2,pmaxSize);
      }
    }
    void SetNrelax0(Int pnrelax0){
      if(maxSize>0){
        nrelax0 = min(pnrelax0,maxSize);
      }
      else{
        nrelax0 = pnrelax0;
      }
    }
    void SetNrelax1(Int pnrelax1){
      if(maxSize>0){
      nrelax1 = min(pnrelax1,maxSize);
      }
      else{
        nrelax1 = pnrelax1;
      }
    }
    void SetNrelax2(Int pnrelax2){
      if(maxSize>0){
      nrelax2 = min(pnrelax2,maxSize);
      }
      else{
        nrelax2 = pnrelax2;
      }
    }

  };


  enum MappingType {ROW2D,COL2D,MODWRAP2D,MODWRAP2DNS,WRAP2D,WRAP2DFORCED};
  enum FactorizationType {FANOUT,FANBOTH,FANBOTH_STATIC};
  enum LoadBalanceType {NOLB,NNZ,NCOLS,FLOPS,SUBCUBE,SUBCUBE_NNZ};
  enum OrderingType {NATURAL,MMD,AMD,NDBOX,NDGRID,SCOTCH,PTSCOTCH,METIS,PARMETIS,USER};
  enum SchedulerType {DL,MCT,PR,FIFO};
  //enum OrderRefinementType {BarryDL,MCT,PR,FIFO};
  class NGCholOptions{
    public:
      MappingType mappingType;
      std::string mappingTypeStr;
      FactorizationType factorization;
      LoadBalanceType load_balance;
      std::string load_balance_str;
      std::string order_refinement_str;
      OrderingType ordering;
      SchedulerType scheduler;
      Int maxIsend;
      Int maxIrecv;
      Int maxSnode;
      CommEnvironment * commEnv;
      RelaxationParameters relax;
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
      NGCholOptions(){
        mappingType = MODWRAP2D;
        mappingTypeStr = "MODWRAP2D";
        factorization = FANBOTH;
        load_balance = NNZ;
        load_balance_str = "SUBCUBE";
        scheduler = DL;
        maxIsend = 0;
        maxIrecv=0;
        maxSnode=-1;
        relax = RelaxationParameters(maxSnode);
        commEnv = NULL;
        ordering = MMD;
        order_refinement_str = "NONE";
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
    //SYMPACK::vector<bool> is_factor_sent;

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
        SYMPACK::vector<T> tmpBuf;
        SYMPACK::vector<Int> src_colindx;
        SYMPACK::vector<Int> src_to_tgt_offset;

        void Resize(Int size, Int mw){
          tmpBuf.resize(size*mw);
          src_colindx.resize(mw);
          src_to_tgt_offset.resize(size);
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










}

#endif
