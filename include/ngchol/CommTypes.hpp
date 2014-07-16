#ifndef _COMM_TYPES_DECL_HPP_
#define _COMM_TYPES_DECL_HPP_


#include <list>
#include <deque>
#include <queue>
#include <vector>
#include <mpi.h>

#include "ngchol/Environment.hpp"

namespace LIBCHOLESKY{

  struct SnodeUpdate;
  struct LocalUpdate;
  struct DelayedComm;
  struct Icomm;
  struct DelayedCommCompare;
  class AsyncComms;
  typedef std::priority_queue<DelayedComm,deque<DelayedComm>,DelayedCommCompare> CommList;
  template <typename T> inline void Serialize( Icomm& os,  const T * val, const Int count);
  template <typename T> inline Icomm& operator<<( Icomm& os,  const T & val);


  class CommEnvironment;



  struct SnodeUpdate{
    Int tgt_snode_id;
    Int src_fr;
    SnodeUpdate(Int aSnodeId, Int aSrcFr):tgt_snode_id(aSnodeId),src_fr(aSrcFr){};
  };


  struct LocalUpdate{
    Int src_snode_id;
    Int src_nzblk_idx;
    Int src_first_row;
    LocalUpdate(Int snode_id,Int nzblk_idx,Int first_row):src_snode_id(snode_id),src_nzblk_idx(nzblk_idx),src_first_row(first_row){};
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
        MPI_Test(&Request,&flag,&recv_status);
        if(!flag){
          //        MPI_Wait(&Request,&recv_status);
          //          MPI_Test_cancelled(&recv_status,&flag);
          //          if(!flag){
          //          }
          //            cout<<"CANCEL"<<endl;
          MPI_Cancel(&Request);
        }
        MPI_Request_free(&Request);
      }
    };
    inline char * back(){ return &pSrcBlocks->at(head);}
    inline char * front(){ return &pSrcBlocks->front();}
    //inline Int size(){ return pSrcBlocks->size();}
    inline Int size(){ return head;}
    inline void resize(Int size){
      head = 0; 
//      pSrcBlocks->resize(0); 
//      pSrcBlocks->reserve(size);
      pSrcBlocks->resize(size);
    }
    inline void setHead(Int phead){ 
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
        MPI_Test(&Request,&flag,&recv_status);
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
      return a.tgt_snode_id>b.tgt_snode_id;
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

#endif //_COMM_TYPES_DECL_HPP_
