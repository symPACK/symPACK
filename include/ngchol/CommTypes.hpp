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
  template<typename T> class SuperNode;

  enum TaskType {FACTOR, AGGREGATE, UPDATE};

  class FBDelayedComm;
  class FBDelayedCommCompare;
#ifdef _DEADLOCK_
  typedef std::queue<FBDelayedComm,deque<FBDelayedComm > > FBCommList;
#else
  typedef std::priority_queue<FBDelayedComm,deque<FBDelayedComm >,FBDelayedCommCompare > FBCommList;
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

//   ~FBDelayedComm(){
//        if(type==AGGREGATE){
//          delete src_data;
//        }
//    }

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

//    inline bool base_compare2(const Int & a_src_snode_id, const Int & a_tgt_snode_id, const TaskType & a_type,
//        const Int & b_src_snode_id, const Int & b_tgt_snode_id, const TaskType & b_type) const{
//        if(a_tgt_snode_id>b_tgt_snode_id){
//          return true;
//        }
//        else if(a_tgt_snode_id==b_tgt_snode_id){
//          if(a_src_snode_id==b_src_snode_id){
//            //factor comes first
//            //case Fx,y   vs   Ax,y    =>   F comes first
//            return !(a_type == FACTOR);
//          }
//          else{
//             return a_src_snode_id>b_src_snode_id;
//          }
//        }
//        else{
//          return false;
//        }
//    }


    //Logic is: does a go after b ?
    inline bool compare(const Int & a_src_snode_id, const Int & a_tgt_snode_id, const TaskType & a_type,
        const Int & b_src_snode_id, const Int & b_tgt_snode_id, const TaskType & b_type) const{


//      return base_compare2(a_src_snode_id, a_tgt_snode_id, a_type, b_src_snode_id, b_tgt_snode_id, b_type);


      //If they are the same type, sort by tgt id then by src_id
      if(a_type == b_type){
        return base_compare(a_src_snode_id, a_tgt_snode_id, b_src_snode_id, b_tgt_snode_id);
      }
      //case Fx,y   vs   Ax,y    =>   F comes first
      else if(a_type == FACTOR && /*a_src_snode_id == b_src_snode_id &&*/ a_tgt_snode_id == b_tgt_snode_id){
        return false;
      }
      //case Ax,*   vs   Fx,*    =>   F comes first
      else if(b_type == FACTOR && /*a_src_snode_id == b_src_snode_id &&*/ a_tgt_snode_id == b_tgt_snode_id){
        return true;
      }
//#ifndef _SEPARATE_COMM_
      //case Fy,*   vs   A*,y    =>   A comes first
      else if(a_type == FACTOR && a_src_snode_id == b_tgt_snode_id){
        return true;
      }
      //case A*,y   vs   Fy,*    =>   A comes first
      else if(b_type == FACTOR && b_src_snode_id == a_tgt_snode_id){
        return false;
      }
//#endif
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
//#ifndef _SEPARATE_COMM_
      //case Fx,y   vs   Ax,y    =>   F comes first
      else if(a_type == FACTOR && /*a_src_snode_id == abs(b_src_snode_id) &&*/ a_tgt_snode_id == b_tgt_snode_id){
        return false;
      }
      //case Ax,*   vs   Fx,*    =>   F comes first
      else if(b_type == FACTOR && /*a_src_snode_id == abs(b_src_snode_id) &&*/ a_tgt_snode_id == b_tgt_snode_id){
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
//#endif
      else{
        return base_compare(a_src_snode_id, a_tgt_snode_id, abs(b_src_snode_id), b_tgt_snode_id);
      }
    }













    bool operator()(const FBDelayedComm & a,const FBDelayedComm & b) const
    {
//      return compare(a.src_data->Id(),a.tgt_snode_id,a.type,b.src_data->Id(),b.tgt_snode_id,b.type);
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
  typedef std::queue<DelayedComm,deque<DelayedComm> > CommList;
  typedef std::queue<DelayedComm,deque<DelayedComm> > DownCommList;
  typedef std::queue<SnodeUpdateFB,deque<SnodeUpdateFB> > FBTasks;
#else
  typedef std::priority_queue<DelayedComm,deque<DelayedComm>,DelayedCommCompare> CommList;
  typedef std::priority_queue<DelayedComm,deque<DelayedComm>,DelayedCommReverseCompare> DownCommList;
  typedef std::priority_queue<SnodeUpdateFB,deque<SnodeUpdateFB>,SnodeUpdateFBCompare> FBTasks;
#endif

  template <typename T> inline void Serialize( Icomm& os,  const T * val, const Int count);
  template <typename T> inline Icomm& operator<<( Icomm& os,  const T & val);


  class CommEnvironment;


  struct SnodeUpdateFB{
    TaskType type;
    Int src_snode_id;
    Int tgt_snode_id;
    //dependencies
    Int in_deps;
    //unused but preparing for task scheduling priorities
    Int rank;
    SnodeUpdateFB():rank(-1),in_deps(0){}
  };


  struct SnodeUpdateFBCompare{
    bool operator()(const SnodeUpdateFB & a,const SnodeUpdateFB & b) const
    {

      bool b_factor = b.tgt_snode_id == b.src_snode_id;


      //use the ranks first
      if(a.rank>=0 && b.rank>=0){
        return a.rank<b.rank;
      }
    


      //check whether it is an update or a factorization
      bool a_factor = a.tgt_snode_id == a.src_snode_id;

      //if same type apply the usual priorities
//      if(a_factor == b_factor){

      //use the classic priorities otherwise
      if(a.tgt_snode_id>b.tgt_snode_id){
        return true;
      }
      else if(a.tgt_snode_id==b.tgt_snode_id){
        return a.src_snode_id>b.src_snode_id;
      }
      else{
        return false;
      }

//      }
//      else if (a_factor){
//        if(a.tgt_snode_id
//      }
//      else if (b_factor){
//
//      }
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

  struct Icomm{


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
      TIMER_START(ICOMM_MALLOC);
      pSrcBlocks = new std::vector<char>();
      TIMER_STOP(ICOMM_MALLOC);
      head = 0;
    };

    Icomm(size_t aSize, MPI_Request aRequest):Request(aRequest){
      TIMER_START(ICOMM_MALLOC);
      pSrcBlocks = new std::vector<char>(aSize);
      TIMER_STOP(ICOMM_MALLOC);
      head = 0;
    };

    Icomm(size_t aSize){
      Request = MPI_REQUEST_NULL;
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
      pSrcBlocks->resize(size);
    }
    inline void setHead(Int phead){ 
      assert(phead <= pSrcBlocks->size()); 
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

  template <typename T> inline Icomm& operator>>( Icomm& os,  T & val){
    TIMER_START(ICOMM_DESERIALIZE);
    T* dest = reinterpret_cast<T*>(&os.pSrcBlocks->at(os.head));
    std::copy(dest,dest+1,&val);
    os.head += sizeof(T);
    TIMER_STOP(ICOMM_DESERIALIZE);

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
      inline Int MPI_Size(){return MPI_size_;}
      inline Int MPI_Rank(){return MPI_rank_;}
  };




  struct SnodeUpdateOld{
    Int tgt_snode_id;
    Int src_fr;
    SnodeUpdateOld(Int aSnodeId, Int aSrcFr):tgt_snode_id(aSnodeId),src_fr(aSrcFr){};
  };





}

#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) TAU_FSTART(a);
#define TIMER_STOP(a) TAU_FSTOP(a);
#endif
#endif




#endif //_COMM_TYPES_DECL_HPP_
