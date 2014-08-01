#ifndef _COMM_TYPES_DECL_HPP_
#define _COMM_TYPES_DECL_HPP_


#include <list>
#include <deque>
#include <queue>
#include <vector>
#include <mpi.h>

#include "ngchol/Environment.hpp"
#include "ngchol/NumMat.hpp"
//#include "ngchol/debug.hpp"

#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) 
#define TIMER_STOP(a) 
#endif
#endif




namespace LIBCHOLESKY{

  template<typename T> class NumMat;


  struct SnodeUpdate;
  struct LocalUpdate;
  struct DelayedComm;
  struct Icomm;
  struct DelayedCommCompare;
  struct DelayedCommReverseCompare;
  class AsyncComms;
  typedef std::priority_queue<DelayedComm,deque<DelayedComm>,DelayedCommCompare> CommList;
  typedef std::priority_queue<DelayedComm,deque<DelayedComm>,DelayedCommReverseCompare> DownCommList;

  template <typename T> inline void Serialize( Icomm& os,  const T * val, const Int count);
  template <typename T> inline Icomm& operator<<( Icomm& os,  const T & val);


  class CommEnvironment;



  struct SnodeUpdateOld{
    Int tgt_snode_id;
    Int src_fr;
    SnodeUpdateOld(Int aSnodeId, Int aSrcFr):tgt_snode_id(aSnodeId),src_fr(aSrcFr){};
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
    //IntNumVec src_rowindx;
    IntNumVec src_to_tgt_offset;
    TempUpdateBuffers(Int size, Int mw){
      tmpBuf.Resize(size,mw);
      src_colindx.Resize(mw);
      //src_rowindx.Resize(size);
      src_to_tgt_offset.Resize(size);
    }
  };



  struct DelayedComm{
    Int tgt_snode_id;
    Int src_snode_id;

    Int src_nzblk_idx;
    Int src_first_row;
    Int src_last_row;

    DelayedComm(Int a_src_snode_id, Int a_tgt_snode_id, Int a_src_nzblk_idx, Int a_src_first_row):src_snode_id(a_src_snode_id),tgt_snode_id(a_tgt_snode_id),src_first_row(a_src_first_row),src_nzblk_idx(a_src_nzblk_idx){};
  };

  struct Icomm{
    std::vector<char> * pSrcBlocks;
    Int head;
    MPI_Request Request;
    Icomm(){
      Request = MPI_REQUEST_NULL;
      TIMER_START(ICOMM_MALLOC);
      pSrcBlocks = new std::vector<char>();
      TIMER_STOP(ICOMM_MALLOC);
      head = 0;
    };

    Icomm(Int aSize, MPI_Request aRequest):Request(aRequest){
      TIMER_START(ICOMM_MALLOC);
      pSrcBlocks = new std::vector<char>(aSize);
      TIMER_STOP(ICOMM_MALLOC);
      head = 0;
    };
    ~Icomm(){  
      delete pSrcBlocks; 
      if(Request !=MPI_REQUEST_NULL){
        MPI_Status recv_status;
        int flag = 0;

//        set_mpi_handler(MPI_COMM_WORLD);
        int error_code = MPI_Test(&Request,&flag,&recv_status);
//        check_mpi_error(error_code, MPI_COMM_WORLD, true);
        if(!flag){
          MPI_Cancel(&Request);
        }
        MPI_Request_free(&Request);
      }
    };
    inline char * back(){ return &pSrcBlocks->at(head);}
    inline char * front(){ return &pSrcBlocks->front();}
    inline Int capacity(){ return pSrcBlocks->size();}
    inline Int size(){ return head;}
    inline void resize(Int size){
      head = 0; 
      pSrcBlocks->resize(size);
    }
    inline void setHead(Int phead){ 
//if(this == 0x1cb5e00){gdb_lock();}
      assert(phead <= pSrcBlocks->size()); 
      head = phead;
    }
    inline void clear(){ 
      head = 0; 
      pSrcBlocks->resize(0); 
//      pSrcBlocks->clear();
 
      if(Request !=MPI_REQUEST_NULL){
        MPI_Status recv_status;
        int flag = 0;
//        set_mpi_handler(MPI_COMM_WORLD);
        int error_code = MPI_Test(&Request,&flag,&recv_status);
//        check_mpi_error(error_code, MPI_COMM_WORLD, true);
        assert(!flag);
        MPI_Request_free(&Request);
      }
      Request = MPI_REQUEST_NULL;     
    }

  };

  template <typename T> inline void Serialize( Icomm& os,  const T * val, const Int count){
    TIMER_START(ICOMM_SERIALIZE);
    T* dest = reinterpret_cast<T*>(&os.pSrcBlocks->at(os.head));
    std::copy(val,val+count,dest);
    os.head += count*sizeof(T);
    TIMER_STOP(ICOMM_SERIALIZE);
  }

  template <typename T> inline Icomm& operator<<( Icomm& os,  const T & val){
    TIMER_START(ICOMM_SERIALIZE);
    T* dest = reinterpret_cast<T*>(&os.pSrcBlocks->at(os.head));
    std::copy(&val,&val+1,dest);
    os.head += sizeof(T);
    TIMER_STOP(ICOMM_SERIALIZE);

    return os;
  }


  struct DelayedCommCompare{
    bool operator()(const DelayedComm & a,const DelayedComm & b) const
    {
//      bool return_value =  a.tgt_snode_id<b.tgt_snode_id;
      if(a.tgt_snode_id>b.tgt_snode_id){
        return true;
      }
      else if(a.tgt_snode_id==b.tgt_snode_id){
       return a.src_snode_id>b.src_snode_id;
      }
      else{
        return false;
      }
      //return a.tgt_snode_id>b.tgt_snode_id;
    }
  };


  struct DelayedCommReverseCompare{
    bool operator()(const DelayedComm & a,const DelayedComm & b) const
    {
//      bool return_value =  a.tgt_snode_id<b.tgt_snode_id;
      if(a.tgt_snode_id<b.tgt_snode_id){
        return true;
      }
      else if(a.tgt_snode_id==b.tgt_snode_id){
       return a.src_snode_id<b.src_snode_id;
      }
      else{
        return false;
      }
      //return a.tgt_snode_id>b.tgt_snode_id;
    }
  };


  class AsyncComms{
    protected:
      std::list<Icomm *> list_;


    public:
      typedef std::list<Icomm *>::iterator iterator;
      //push_back
      //pop_back
      //back

      void push_back(Icomm * elem){
        list_.push_back(elem);
      }
      void pop_back(){
        Icomm * elem = list_.back();
        delete elem;
        list_.pop_back();
      }

      Icomm * & back(){
        return list_.back();
      }

      int size(){
        return list_.size();
      }

      bool empty(){
        return list_.empty();
      }

      iterator begin(){
        return list_.begin();
      }

      iterator end(){
        return list_.end();
      }

      void clear(){
        for(iterator it = list_.begin();it!=list_.end();++it){
          delete *it;
        }
        list_.clear();
      }

      iterator erase(iterator position){
        static int count = 0;
        count++;

        if(position!=end()){
          Icomm * elem = *position;
          delete elem;
        }
        return list_.erase(position);
      }

      friend std::ostream& operator<<(std::ostream& out, AsyncComms& comm) // output 
      {
#ifdef _DEBUG_DELAY_
        out <<"(";
        for(AsyncComms::iterator it = comm.list_.begin();it!=comm.list_.end();++it){
          Icomm * elem = *it;
          Int src_snode_id = *(Int*)elem->front();
          out <<" " << src_snode_id;
        }
        out <<" )";
#endif
        return out; 
      } 



  };



  class CommEnvironment{
    protected:
      bool        isMpi_;
      /// @brief MPI communicator
      MPI_Comm    pComm_;
      Int         MPI_rank_;
      Int         MPI_size_;
    public:
      CommEnvironment(){
          isMpi_ = false;
          //pComm_ = MPI_COMM_NULL;
          MPI_size_ = -1;
          MPI_rank_ = -1;
      }

      CommEnvironment(MPI_Comm & aComm){
        if( aComm != MPI_COMM_NULL){
          isMpi_ = true;
          MPI_Comm_dup(aComm,&pComm_);
          MPI_Comm_size(pComm_,&MPI_size_);
          MPI_Comm_rank(pComm_,&MPI_rank_);
        }
        else{
          isMpi_ = false;
          pComm_ = MPI_COMM_NULL;
          MPI_size_ = -1;
          MPI_rank_ = -1;
        }
      }

      inline bool IsMPI(){return isMpi_;}
      inline MPI_Comm & MPI_GetComm() {return pComm_;}
      inline Int MPI_Size(){return MPI_size_;}
      inline Int MPI_Rank(){return MPI_rank_;}
  };







}

#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) TAU_FSTART(a);
#define TIMER_STOP(a) TAU_FSTOP(a);
#endif
#endif




#endif //_COMM_TYPES_DECL_HPP_
