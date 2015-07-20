#ifndef _TYPES_DECL_HPP_
#define _TYPES_DECL_HPP_

#include "ngchol/Environment.hpp"
#include "ngchol/timer.hpp"
#include "ngchol/NumVec.hpp"
#include "ngchol/NumMat.hpp"
#include "ngchol/CommTypes.hpp"

namespace LIBCHOLESKY{

  enum MappingType {ROW2D,COL2D,MODWRAP2D,MODWRAP2DNS,WRAP2D,WRAP2DFORCED};
  enum FactorizationType {FANOUT,FANBOTH};
  enum LoadBalanceType {NOLB,NNZ,NCOLS,FLOPS,SUBCUBE,SUBCUBE_NNZ};
  enum OrderingType {MMD,AMD};
  enum SchedulerType {DL,MCT};
  class NGCholOptions{
    public:
      MappingType mappingType;
      FactorizationType factorization;
      LoadBalanceType load_balance;
      OrderingType ordering;
      SchedulerType scheduler;
      Int maxIsend;
      Int maxIrecv;
      Int maxSnode;
      CommEnvironment * commEnv;

    protected:
      bool isSqrtP(){
        bool val = false;
        switch(mappingType){
          case MODWRAP2D: case MODWRAP2DNS: case WRAP2D: case WRAP2DFORCED:
            val = true;
            break;
          default:
            val = false;
            break;
        }
        return val;
      }

    public:
      NGCholOptions(){
        mappingType = MODWRAP2D;
        factorization = FANBOTH;
        load_balance = NNZ;
        scheduler = DL;
        maxIsend = 0;
        maxIrecv=0;
        maxSnode=-1;
        commEnv = NULL;
        ordering = MMD;
      }

      Int used_procs(Int np){
        if(isSqrtP()){
          Int nr = (Int)sqrt((double)np);
          Int nc = np/nr;
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
    NumMat<T> tmpBuf;
    IntNumVec src_colindx;
    IntNumVec src_to_tgt_offset;

    void Resize(Int size, Int mw){
      tmpBuf.Resize(size,mw);
      src_colindx.Resize(mw);
      src_to_tgt_offset.Resize(size);
    }

    void Clear(){
      tmpBuf.Clear();
      src_colindx.Clear();
      src_to_tgt_offset.Clear();
     }


    TempUpdateBuffers(){
    }
    TempUpdateBuffers(Int size, Int mw){
      Resize(size,mw);
    }
  };










}

#endif
