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

  enum FBTaskType {FACTOR, AGGREGATE, COMPFACT,COMPAGG};

#ifndef _USE_TAU_
  template<typename T> class FBDelayedComm;
  template<typename T> class FBDelayedCommCompare;
#ifdef _DEADLOCK_
  template<typename T> using FBCommList = std::queue<FBDelayedComm<T>,deque<FBDelayedComm<T> > >;
#else
  template<typename T> using FBCommList = std::priority_queue<FBDelayedComm<T>,deque<FBDelayedComm<T> >,FBDelayedCommCompare<T> >;
#endif

  template<typename T>
  class FBDelayedComm{
    public:
    FBTaskType type;
    SuperNode<T> * src_data;
    Int tgt_snode_id;

    Int src_nzblk_idx;
    Int src_first_row;
    Int src_last_row;

    FBDelayedComm(FBTaskType a_type,SuperNode<T> * a_src_data,Int a_src_snode_id, Int a_tgt_snode_id, Int a_src_nzblk_idx,
                   Int a_src_first_row, Int a_target, Int a_tag, Int a_count=0){
          type = a_type;
          src_data = a_src_data;
          tgt_snode_id = a_tgt_snode_id;
//          src_snode_id = a_src_snode_id;
          src_first_row = a_src_first_row;
          src_nzblk_idx = a_src_nzblk_idx;
//          target = a_target;
//          tag = a_tag;
//          count = a_count;
   }

//   ~FBDelayedComm(){
//        if(type==AGGREGATE){
//          delete src_data;
//        }
//    }

  };


 

  template<typename T>
  class FBDelayedCommCompare{
    public:
    //Logic is: does a go after b ?
    inline bool compare(const Int & a_src_snode_id, const Int & a_tgt_snode_id, const FBTaskType & a_type,
        const Int & b_src_snode_id, const Int & b_tgt_snode_id, const FBTaskType & b_type) const{
      //If they are the same type, sort by tgt id then by src_id
      if(a_type == b_type){
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
      //If they are not of the same type
      else{
        if(a_type == FACTOR){
          //case Fx,* vs Ax,*
          if(a_src_snode_id == b_src_snode_id){
            return true; //a is factor ?
          }
          else{
            //case Fx,* vs A*,x
            if(a_src_snode_id == b_tgt_snode_id){
              return true;
            }
            else{
              //              return a_src_snode_id<b_tgt_snode_id;
              //              return b_src_snode_id<a_src_snode_id;
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
          }
        }
        else{

          //case Fx,* vs Ax,*
          if(b_src_snode_id == a_src_snode_id){
            return false; //a is factor ?
          }
          else{
            //case Fx,* vs A*,x
            if(b_src_snode_id == a_tgt_snode_id){
              return false;
            }
            else{
              //              return b_src_snode_id>=a_tgt_snode_id;
              //              return b_src_snode_id<a_src_snode_id;
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
          }
        }
      }
    }

    bool unitTest(){
        Int a_src_snode_id;
        Int a_tgt_snode_id;
        FBTaskType a_type;
        Int b_src_snode_id;
        Int b_tgt_snode_id;
        FBTaskType b_type;
        
        bool expected_result, passed;

        passed = true;


        //First case : F1-1 vs F2-2 : F1-1 before F2-2
        a_src_snode_id = 1;
        a_tgt_snode_id = 1;
        a_type = FACTOR;
        b_src_snode_id = 2;
        b_tgt_snode_id = 2;
        b_type = FACTOR;
        expected_result= false;
        passed = passed && (compare(a_src_snode_id,a_tgt_snode_id,a_type,b_src_snode_id,b_tgt_snode_id,b_type) == expected_result);
        expected_result= true;
        passed = passed && (compare(b_src_snode_id,b_tgt_snode_id,b_type,a_src_snode_id,a_tgt_snode_id,a_type) == expected_result);

        //Second case : A1-1 vs A2-2 : A1-1 before A2-2
        a_src_snode_id = 1;
        a_tgt_snode_id = 1;
        a_type = AGGREGATE;
        b_src_snode_id = 2;
        b_tgt_snode_id = 2;
        b_type = AGGREGATE;
        expected_result= false;
        passed = passed && (compare(a_src_snode_id,a_tgt_snode_id,a_type,b_src_snode_id,b_tgt_snode_id,b_type) == expected_result);
        expected_result= true;
        passed = passed && (compare(b_src_snode_id,b_tgt_snode_id,b_type,a_src_snode_id,a_tgt_snode_id,a_type) == expected_result);


        //Third case : A1-4 vs F4-6 : A1-4 before F4-6
        a_src_snode_id = 1;
        a_tgt_snode_id = 4;
        a_type = AGGREGATE;
        b_src_snode_id = 4;
        b_tgt_snode_id = 6;
        b_type = FACTOR;
        expected_result= false;
        passed = passed && (compare(a_src_snode_id,a_tgt_snode_id,a_type,b_src_snode_id,b_tgt_snode_id,b_type) == expected_result);
        expected_result= true;
        passed = passed && (compare(b_src_snode_id,b_tgt_snode_id,b_type,a_src_snode_id,a_tgt_snode_id,a_type) == expected_result);

        //4th case : A2-4 vs F1-4 : F1-4 before A2-4
        a_src_snode_id = 2;
        a_tgt_snode_id = 4;
        a_type = AGGREGATE;
        b_src_snode_id = 1;
        b_tgt_snode_id = 4;
        b_type = FACTOR;
        expected_result= false;
        passed = passed && (compare(a_src_snode_id,a_tgt_snode_id,a_type,b_src_snode_id,b_tgt_snode_id,b_type) == expected_result);
        expected_result= true;
        passed = passed && (compare(b_src_snode_id,b_tgt_snode_id,b_type,a_src_snode_id,a_tgt_snode_id,a_type) == expected_result);

        //5th case : A4-4 vs F4-4 : A4-4 before F4-4
        a_src_snode_id = 4;
        a_tgt_snode_id = 4;
        a_type = AGGREGATE;
        b_src_snode_id = 4;
        b_tgt_snode_id = 4;
        b_type = FACTOR;
        expected_result= false;
        passed = passed && (compare(a_src_snode_id,a_tgt_snode_id,a_type,b_src_snode_id,b_tgt_snode_id,b_type) == expected_result);
        expected_result= true;
        passed = passed && (compare(b_src_snode_id,b_tgt_snode_id,b_type,a_src_snode_id,a_tgt_snode_id,a_type) == expected_result);





        return passed;
    }

    bool operator()(const FBDelayedComm<T> & a,const FBDelayedComm<T> & b) const
    {
      return compare(a.src_data->Id(),a.tgt_snode_id,a.type,b.src_data->Id(),b.tgt_snode_id,b.type);
//      return compare(a.src_snode_id,a.tgt_snode_id,a.type,b.src_snode_id,b.tgt_snode_id,b.type);
    }
  };
#else
  class FBDelayedComm;
  class FBDelayedCommCompare;
#ifdef _DEADLOCK_
  typedef std::queue<FBDelayedComm,deque<FBDelayedComm > > FBCommList;
#else
  typedef std::priority_queue<FBDelayedComm,deque<FBDelayedComm >,FBDelayedCommCompare > FBCommList;
#endif

  class FBDelayedComm{
    public:
    FBTaskType type;
    void * src_data;
    Int tgt_snode_id;
    Int src_snode_id;

    Int src_nzblk_idx;
    Int src_first_row;
    Int src_last_row;

    FBDelayedComm(FBTaskType a_type,void * a_src_data,Int a_src_snode_id, Int a_tgt_snode_id, Int a_src_nzblk_idx,
                   Int a_src_first_row, Int a_target, Int a_tag, Int a_count=0){
          type = a_type;
          src_data = a_src_data;
          tgt_snode_id = a_tgt_snode_id;
          src_snode_id = a_src_snode_id;
          src_first_row = a_src_first_row;
          src_nzblk_idx = a_src_nzblk_idx;
//          target = a_target;
//          tag = a_tag;
//          count = a_count;
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

    inline bool base_compare2(const Int & a_src_snode_id, const Int & a_tgt_snode_id, const FBTaskType & a_type,
        const Int & b_src_snode_id, const Int & b_tgt_snode_id, const FBTaskType & b_type) const{
        if(a_tgt_snode_id>b_tgt_snode_id){
          return true;
        }
        else if(a_tgt_snode_id==b_tgt_snode_id){
          if(a_src_snode_id==b_src_snode_id){
            //factor comes first
            //case Fx,y   vs   Ax,y    =>   F comes first
            return !(a_type == FACTOR);
          }
          else{
             return a_src_snode_id>b_src_snode_id;
          }
        }
        else{
          return false;
        }
    }


    //Logic is: does a go after b ?
    inline bool compare(const Int & a_src_snode_id, const Int & a_tgt_snode_id, const FBTaskType & a_type,
        const Int & b_src_snode_id, const Int & b_tgt_snode_id, const FBTaskType & b_type) const{


//      return base_compare2(a_src_snode_id, a_tgt_snode_id, a_type, b_src_snode_id, b_tgt_snode_id, b_type);


      //If they are the same type, sort by tgt id then by src_id
      if(a_type == b_type){
        return base_compare(a_src_snode_id, a_tgt_snode_id, b_src_snode_id, b_tgt_snode_id);
      }
      //case Fx,y   vs   Ax,y    =>   F comes first
      else if(a_type == FACTOR && a_src_snode_id == b_src_snode_id && a_tgt_snode_id == b_tgt_snode_id){
        return false;
      }
      //case Fy,*   vs   A*,y    =>   A comes first
      else if(a_type == FACTOR && a_src_snode_id == b_tgt_snode_id){
        return true;
      }
      //case Ax,*   vs   Fx,*    =>   F comes first
      else if(b_type == FACTOR && a_src_snode_id == b_src_snode_id && a_tgt_snode_id == b_tgt_snode_id){
        return true;
      }
      //case A*,y   vs   Fy,*    =>   A comes first
      else if(b_type == FACTOR && b_src_snode_id == a_tgt_snode_id){
        return false;
      }
      else{
        return base_compare(a_src_snode_id, a_tgt_snode_id, b_src_snode_id, b_tgt_snode_id);
      }

//      if(a_type == b_type){
//        return base_compare(a_src_snode_id, a_tgt_snode_id, b_src_snode_id, b_tgt_snode_id);
//      }
//      //If they are not of the same type
//      else{
//        if(a_type == FACTOR){
//          //case Fx,* vs Ax,*
//          if(a_src_snode_id == b_src_snode_id){
//            return true; //a is factor ?
//          }
//          else{
//            //case Fx,* vs A*,x
//            if(a_src_snode_id == b_tgt_snode_id){
//              return true;
//            }
//            else{
//              //              return a_src_snode_id<b_tgt_snode_id;
//              //              return b_src_snode_id<a_src_snode_id;
//              return base_compare(a_src_snode_id, a_tgt_snode_id, b_src_snode_id, b_tgt_snode_id);
//            } 
//          }
//        }
//        else{
//
//          //case Fx,* vs Ax,*
//          if(b_src_snode_id == a_src_snode_id){
//            return false; //a is factor ?
//          }
//          else{
//            //case Fx,* vs A*,x
//            if(b_src_snode_id == a_tgt_snode_id){
//              return false;
//            }
//            else{
//              return base_compare(a_src_snode_id, a_tgt_snode_id, b_src_snode_id, b_tgt_snode_id);
//            } 
//          }
//        }
//      }
    }

    bool operator()(const FBDelayedComm & a,const FBDelayedComm & b) const
    {
//      return compare(a.src_data->Id(),a.tgt_snode_id,a.type,b.src_data->Id(),b.tgt_snode_id,b.type);
      return compare(a.src_snode_id,a.tgt_snode_id,a.type,b.src_snode_id,b.tgt_snode_id,b.type);
    }
  };
#endif















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
    Int src_snode_id;
    Int tgt_snode_id;
  };


  struct SnodeUpdateFBCompare{
    bool operator()(const SnodeUpdateFB & a,const SnodeUpdateFB & b) const
    {

      //check whether it is an update or a factorization
      bool a_factor = a.tgt_snode_id == a.src_snode_id;
      bool b_factor = b.tgt_snode_id == b.src_snode_id;

      //if same type apply the usual priorities
//      if(a_factor == b_factor){
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
      //head = 0; 
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
