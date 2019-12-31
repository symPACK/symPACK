#ifndef _COMM_TYPES_DECL_HPP_
#define _COMM_TYPES_DECL_HPP_

#include <list>
#include <deque>
#include <queue>
#include <vector>
#include <mpi.h>

#include "sympack/Environment.hpp"
#include "sympack/Task.hpp"

#ifdef NO_INTRA_PROFILE
#if defined (SPROFILE)
#define SYMPACK_TIMER_START(a) 
#define SYMPACK_TIMER_STOP(a) 
#endif
#endif

namespace symPACK{

  template<typename T, class Allocator> class SuperNode;


  class FBDelayedComm;
  class FBDelayedCommCompare;
#ifdef _DEADLOCK_
  typedef std::queue<FBDelayedComm,std::deque<FBDelayedComm > > FBCommList;
#else
  typedef std::priority_queue<FBDelayedComm,std::deque<FBDelayedComm >,FBDelayedCommCompare > FBCommList;
#endif

  class FBDelayedComm{
    public:
    TaskType type;
    void * src_data;
    Int tgt_snode_id;
    Int src_snode_id;

    Int src_nzblk_idx;
    Int src_first_row;
    Int src_last_row;

    FBDelayedComm(TaskType a_type,void * a_src_data,Int a_src_snode_id, Int a_tgt_snode_id, Int a_src_nzblk_idx,
                   Int a_src_first_row, Int a_target, Int a_tag, Int a_count=0){
          type = a_type;
          src_data = a_src_data;
          tgt_snode_id = a_tgt_snode_id;
          src_snode_id = a_src_snode_id;
          src_first_row = a_src_first_row;
          src_nzblk_idx = a_src_nzblk_idx;
   }


  };


 

  class FBDelayedCommCompare{
    public:
    inline bool base_compare(const Int & a_src_snode_id, const Int & a_tgt_snode_id,
        const Int & b_src_snode_id, const Int & b_tgt_snode_id) const{
        if(a_tgt_snode_id>b_tgt_snode_id){
          return true;
        }
        else if(a_tgt_snode_id==b_tgt_snode_id){
          return a_src_snode_id>b_src_snode_id;
        }
        else{
          return false;
        }
    }


    //Logic is: does a go after b ?
    inline bool compare(const Int & a_src_snode_id, const Int & a_tgt_snode_id, const TaskType & a_type,
        const Int & b_src_snode_id, const Int & b_tgt_snode_id, const TaskType & b_type) const{




      //If they are the same type, sort by tgt id then by src_id
      if(a_type == b_type){
        return base_compare(a_src_snode_id, a_tgt_snode_id, b_src_snode_id, b_tgt_snode_id);
      }
      //case Fx,y   vs   Ax,y    =>   F comes first
      else if(a_type == FACTOR &&  a_tgt_snode_id == b_tgt_snode_id){
        return false;
      }
      //case Ax,*   vs   Fx,*    =>   F comes first
      else if(b_type == FACTOR &&  a_tgt_snode_id == b_tgt_snode_id){
        return true;
      }
      //case Fy,*   vs   A*,y    =>   A comes first
      else if(a_type == FACTOR && a_src_snode_id == b_tgt_snode_id){
        return true;
      }
      //case A*,y   vs   Fy,*    =>   A comes first
      else if(b_type == FACTOR && b_src_snode_id == a_tgt_snode_id){
        return false;
      }
      else{
        return base_compare(a_src_snode_id, a_tgt_snode_id, b_src_snode_id, b_tgt_snode_id);
      }
    }




    // b is the task
    inline bool compare_task(const Int & a_src_snode_id, const Int & a_tgt_snode_id, const TaskType & a_type,
        const Int & b_src_snode_id, const Int & b_tgt_snode_id, const TaskType & b_type) const{


      //If they are the same type, sort by tgt id then by src_id
      if(a_type == b_type){
        return base_compare(a_src_snode_id, a_tgt_snode_id, b_src_snode_id, b_tgt_snode_id);
      }
      //case Fx,y   vs   Ax,y    =>   F comes first
      else if(a_type == FACTOR &&  a_tgt_snode_id == b_tgt_snode_id){
        return false;
      }
      //case Ax,*   vs   Fx,*    =>   F comes first
      else if(b_type == FACTOR &&  a_tgt_snode_id == b_tgt_snode_id){
        return true;
      }
      //case Fy,*   vs   A*,y    =>   A comes first
      else if(a_type == FACTOR && a_src_snode_id == b_tgt_snode_id){
        return true;
      }
      //case A*,y   vs   Fy,*    =>   A comes first
      else if(b_type == FACTOR && abs(b_src_snode_id) == a_tgt_snode_id){
        return false;
      }
      else{
        return base_compare(a_src_snode_id, a_tgt_snode_id, abs(b_src_snode_id), b_tgt_snode_id);
      }
    }













    bool operator()(const FBDelayedComm & a,const FBDelayedComm & b) const
    {
      return compare(a.src_snode_id,a.tgt_snode_id,a.type,b.src_snode_id,b.tgt_snode_id,b.type);
    }
  };















  struct SnodeUpdate;
  struct LocalUpdate;
  struct DelayedComm;
  struct Icomm;
  struct DelayedCommCompare;
  struct DelayedCommReverseCompare;
  struct SnodeUpdateFB;
  struct SnodeUpdateFBCompare;

  class AsyncComms;
#ifdef _DEADLOCK_
  typedef std::queue<DelayedComm,std::deque<DelayedComm> > CommList;
  typedef std::queue<DelayedComm,std::deque<DelayedComm> > DownCommList;
  typedef std::queue<SnodeUpdateFB,std::deque<SnodeUpdateFB> > FBTasks;
#else
  typedef std::priority_queue<DelayedComm,std::deque<DelayedComm>,DelayedCommCompare> CommList;
  typedef std::priority_queue<DelayedComm,std::deque<DelayedComm>,DelayedCommReverseCompare> DownCommList;
#ifdef _TASKLIST_
  typedef std::list<SnodeUpdateFB> FBTasks;
#else
  typedef std::priority_queue<SnodeUpdateFB,std::deque<SnodeUpdateFB>,SnodeUpdateFBCompare> FBTasks;
#endif

#endif

  template <typename T> inline void Serialize( Icomm& os,  const T * val, const Int count);
  template <typename T> inline Icomm& operator<<( Icomm& os,  const T & val);


  class CommEnvironment;







  struct DelayedComm{
    Int tgt_snode_id;
    Int src_snode_id;

    //for FanBoth
    Int target;
    Int tag;
    void * src_data;
    Int count;


    Int src_nzblk_idx;
    Int src_first_row;
    Int src_last_row;

    DelayedComm(Int a_src_snode_id, Int a_tgt_snode_id, Int a_src_nzblk_idx, Int a_src_first_row):src_snode_id(a_src_snode_id),tgt_snode_id(a_tgt_snode_id),src_first_row(a_src_first_row),src_nzblk_idx(a_src_nzblk_idx){};
    DelayedComm(void * a_src_data, Int a_tgt_snode_id, Int a_src_nzblk_idx, Int a_src_first_row, Int a_target, Int a_tag, Int a_count=0):src_data(a_src_data),tgt_snode_id(a_tgt_snode_id),src_first_row(a_src_first_row),src_nzblk_idx(a_src_nzblk_idx),target(a_target),tag(a_tag),count(a_count){};
  };

  class Icomm{
    public:

#ifdef PREALLOC_IRECV
    list<Icomm*>::iterator bufferPos_;
    void SetPos(list<Icomm*>::iterator pos){
      bufferPos_=pos;
    }

    list<Icomm*>::iterator GetPos(){
      return bufferPos_;
    }
#endif

    std::vector<char> * pSrcBlocks;
    size_t head;
    MPI_Request Request;
    Icomm(){
      Request = MPI_REQUEST_NULL;
      SYMPACK_TIMER_START(ICOMM_MALLOC);
      pSrcBlocks = new std::vector<char>();
      SYMPACK_TIMER_STOP(ICOMM_MALLOC);
      head = 0;
    };

    Icomm(size_t aSize, MPI_Request aRequest):Request(aRequest){
      SYMPACK_TIMER_START(ICOMM_MALLOC);
      pSrcBlocks = new std::vector<char>(aSize);
      SYMPACK_TIMER_STOP(ICOMM_MALLOC);
      head = 0;
    };

    Icomm(size_t aSize){
      Request = MPI_REQUEST_NULL;
      SYMPACK_TIMER_START(ICOMM_MALLOC);
      pSrcBlocks = new std::vector<char>(aSize);
      SYMPACK_TIMER_STOP(ICOMM_MALLOC);
      head = 0;
    };



    ~Icomm(){  
      delete pSrcBlocks; 
      if(Request !=MPI_REQUEST_NULL){
        MPI_Status recv_status;
        int flag = 0;
        int error_code = MPI_Test(&Request,&flag,&recv_status);
        if(!flag){
          MPI_Cancel(&Request);
        }
        MPI_Request_free(&Request);
      }
    };
    inline char * back(){ return &pSrcBlocks->at(head);}
    inline char * front(){ return &pSrcBlocks->front();}
    inline size_t capacity(){ return pSrcBlocks->size();}
    inline size_t size(){ return head;}

    inline void resize(size_t size){
      bassert(size<(size_t)20*1024*1024*1024);
      pSrcBlocks->resize(size);
    }

    inline void setHead(Int phead){ 
      bassert(phead <= pSrcBlocks->size()); 
      head = phead;
    }
    inline void clear(){ 
      head = 0; 
      pSrcBlocks->resize(0); 
 
      if(Request !=MPI_REQUEST_NULL){
        MPI_Status recv_status;
        int flag = 0;
        int error_code = MPI_Test(&Request,&flag,&recv_status);
        assert(!flag);
        MPI_Request_free(&Request);
      }
      Request = MPI_REQUEST_NULL;     
    }

    inline void reset(){ 
      head = 0; 
 
      if(Request !=MPI_REQUEST_NULL){
        MPI_Status recv_status;
        int flag = 0;
        int error_code = MPI_Test(&Request,&flag,&recv_status);
        assert(!flag);
        MPI_Request_free(&Request);
      }
      Request = MPI_REQUEST_NULL;     
    }

    char& operator[](std::size_t idx){ return (*pSrcBlocks)[idx]; }
    const char& operator[](std::size_t idx) const { return (*pSrcBlocks)[idx]; }

  };

  template <typename T> inline void Serialize( Icomm& os,  const T * val, const Int count){
    T* dest = reinterpret_cast<T*>(&os.pSrcBlocks->at(os.head));
    std::copy(val,val+count,dest);
    os.head += count*sizeof(T);
  }

  template <typename T> inline Icomm& operator<<( Icomm& os,  const T & val){
    T* dest = reinterpret_cast<T*>(&os.pSrcBlocks->at(os.head));
    std::copy(&val,&val+1,dest);
    os.head += sizeof(T);
    return os;
  }

  template <typename T> inline Icomm& operator>>( Icomm& os,  T & val){
    T* dest = reinterpret_cast<T*>(&os.pSrcBlocks->at(os.head));
    std::copy(dest,dest+1,&val);
    os.head += sizeof(T);
    return os;
  }



  struct DelayedCommCompare{
    bool operator()(const DelayedComm & a,const DelayedComm & b) const
    {
      if(a.tgt_snode_id>b.tgt_snode_id){
        return true;
      }
      else if(a.tgt_snode_id==b.tgt_snode_id){
       return a.src_snode_id>b.src_snode_id;
      }
      else{
        return false;
      }
    }
  };


  struct DelayedCommReverseCompare{
    bool operator()(const DelayedComm & a,const DelayedComm & b) const
    {
      if(a.tgt_snode_id<b.tgt_snode_id){
        return true;
      }
      else if(a.tgt_snode_id==b.tgt_snode_id){
       return a.src_snode_id<b.src_snode_id;
      }
      else{
        return false;
      }
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

      Icomm * & back(){
        return list_.back();
      }

      size_t size(){
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

      void push_back(Icomm * elem){
        list_.push_back(elem);
      }

#ifdef PREALLOC_IRECV
      list<int>::iterator bufferPos_;


      iterator erase(iterator position){
        return list_.erase(position);
      }

      Icomm * pop_back(){
        Icomm * elem = list_.back();
        list_.pop_back();
        return elem;
      }

      void clear(){
        list_.clear();
      }


#else
      iterator erase(iterator position){
        if(position!=end()){
          Icomm * elem = *position;
          delete elem;
        }
        return list_.erase(position);
      }

      void pop_back(){
        Icomm * elem = list_.back();
        delete elem;
        list_.pop_back();
      }

      void clear(){
        for(iterator it = list_.begin();it!=list_.end();++it){
          delete *it;
        }
        list_.clear();
      }


#endif



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
      int         MPI_rank_;
      int         MPI_size_;
    public:
      CommEnvironment(){
          isMpi_ = false;
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

      CommEnvironment(CommEnvironment & C){
        if(C.isMpi_){
          isMpi_ = true;
          MPI_Comm_dup(C.pComm_,&pComm_);
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


      ~CommEnvironment(){
        if(isMpi_){
          MPI_Comm_free(&pComm_);
        }
        else{

        }
      }

      inline bool IsMPI(){return isMpi_;}
      inline MPI_Comm & MPI_GetComm() {return pComm_;}
      inline int MPI_Size(){return MPI_size_;}
      inline int MPI_Rank(){return MPI_rank_;}
  };




  struct SnodeUpdateOld{
    Int tgt_snode_id;
    Int src_fr;
    SnodeUpdateOld(Int aSnodeId, Int aSrcFr):tgt_snode_id(aSnodeId),src_fr(aSrcFr){};
  };





}

#ifdef NO_INTRA_PROFILE
#if defined (SPROFILE)
#define SYMPACK_TIMER_START(a) SYMPACK_FSTART(a);
#define SYMPACK_TIMER_STOP(a) SYMPACK_FSTOP(a);
#endif
#endif




#endif //_COMM_TYPES_DECL_HPP_
