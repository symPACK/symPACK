#ifndef _COMM_PULL_DECL_HPP_
#define _COMM_PULL_DECL_HPP_


#include <list>
#include <vector>
#include <queue>

#include <upcxx.h>



#include "sympack/Environment.hpp"
//#include "sympack/CommTypes.hpp"





#ifdef NO_INTRA_PROFILE
#if defined (SPROFILE)
#define SYMPACK_TIMER_START(a)
#define SYMPACK_TIMER_STOP(a)
#endif
#endif

namespace symPACK{

struct SnodeUpdateFB;
//class SupernodalMatrixBase;

  struct MsgMetadata{
    //sender taskid
    int src;
    //receiver taskid
    int tgt;
    
    int GIndex;
    //type of message
//    TaskType type;
  };

class BackupBuffer{
    protected:
    bool inUse;
    char * local_ptr;
    size_t size;
    
    public:
    BackupBuffer();
    ~BackupBuffer();
    bool InUse(){return inUse;}
    bool Allocated(){return local_ptr!=NULL;}
    bool AllocLocal(size_t size);
    void DeallocLocal();
    char * GetPtr();
    void ReleasePtr();
};

class IncomingMessage{
  public:
    upcxx::global_ptr<char> remote_ptr;
    size_t msg_size;
    MsgMetadata meta;
    upcxx::event * event_ptr;
    char * local_ptr;
    SnodeUpdateFB * task_ptr;
    bool allocated;
    bool ownLocalStorage;
    bool isDone; 
    bool remoteDealloc; 
    bool isLocal;

    IncomingMessage();
    ~IncomingMessage();
    int Sender();
    bool Wait();
    bool IsDone();
    bool IsLocal();
    bool IsAsync();
    bool AllocLocal();
    upcxx::global_ptr<char> GetRemotePtr();
    char * GetLocalPtr();
    void SetLocalPtr(char * ptr,bool ownStorage = true);
    size_t Size(){return msg_size;}
    void AsyncGet();
    void DeallocRemote();
    void DeallocLocal();
};



struct MSGCompare{
    bool operator()(const IncomingMessage * & a,const IncomingMessage * & b) const
    {
      bool retval = false;


      //use the classic priorities otherwise
      if(a->meta.tgt>b->meta.tgt){
        retval = true;
      }
      else if(a->meta.tgt==b->meta.tgt){
        retval = a->meta.src>b->meta.src;
      }
      else{
        retval = false;
      }

      return retval;
    }
    bool operator()( IncomingMessage * a, IncomingMessage * b) const
    {
      bool retval = false;


      //use the classic priorities otherwise
      if(a->meta.tgt>b->meta.tgt){
        retval = true;
      }
      else if(a->meta.tgt==b->meta.tgt){
        retval = a->meta.src>b->meta.src;
      }
      else{
        retval = false;
      }

      return retval;
    }
  };


  extern size_t gVolComm;
  extern size_t gNumMsg;


  //extern std::priority_queue< IncomingMessage * ,  std::vector<IncomingMessage *>, MSGCompare > gIncomingRecv;
  extern std::list< IncomingMessage * > gIncomingRecv;

  extern std::list< IncomingMessage * > gIncomingRecvAsync;
  extern std::list< IncomingMessage * > gIncomingRecvLocal;
  extern int gMaxIrecv;
  //extern SupernodalMatrixBase * gSuperMatrixPtr;

  void signal_data(upcxx::global_ptr<char> local_ptr, size_t pMsg_size, int dest, MsgMetadata & meta);
  void rcv_async(upcxx::global_ptr<char> pRemote_ptr, size_t pMsg_size, MsgMetadata meta);

  inline void signal_data(upcxx::global_ptr<char> local_ptr, size_t pMsg_size, int dest, MsgMetadata & meta){
      SYMPACK_TIMER_START(SIGNAL_DATA);
      upcxx::async(dest)(rcv_async,local_ptr,pMsg_size,meta);
      SYMPACK_TIMER_STOP(SIGNAL_DATA);
  }


  inline void remote_delete(upcxx::global_ptr<char> pRemote_ptr){
      SYMPACK_TIMER_START(REMOTE_DELETE);
      if(upcxx::myrank()!=pRemote_ptr.where()){
        //logfileptr->OFS()<<"Performing remote delete on P"<<pRemote_ptr.where()<<std::endl;
        //upcxx::async(pRemote_ptr.where())(remote_delete,pRemote_ptr);
        //char * ptr = (char*)pRemote_ptr;
        //logfileptr->OFS()<<"Deallocating UPCXX "<<(uint64_t)ptr<<" "<<std::endl;

        upcxx::deallocate(pRemote_ptr);
      }
      else{

        //char * ptr = (char*)pRemote_ptr;
        //logfileptr->OFS()<<"Deallocating UPCXX "<<(uint64_t)ptr<<" "<<std::endl;

        upcxx::deallocate(pRemote_ptr);
      }
      SYMPACK_TIMER_STOP(REMOTE_DELETE);
  }

  inline void rcv_async(upcxx::global_ptr<char> pRemote_ptr, size_t pMsg_size, MsgMetadata meta){
    SYMPACK_TIMER_START(RCV_ASYNC);

#ifdef HANDLE_LOCAL_POINTER
      char * tryPtr = (char*)pRemote_ptr;
      if(tryPtr!=NULL){
//         logfileptr->OFS()<<"LOCAL POINTER"<<std::endl;
          gIncomingRecvLocal.push_back(new IncomingMessage() );
          IncomingMessage * msg_ptr = gIncomingRecvLocal.back();
          msg_ptr->meta = meta;
          msg_ptr->remote_ptr = pRemote_ptr;
          msg_ptr->local_ptr = tryPtr;
          msg_ptr->msg_size = pMsg_size;
          msg_ptr->isLocal = true;
      }
      else
#endif
      {

        //if we still have async buffers


        bool asyncComm = false;
        if(gIncomingRecvAsync.size() < gMaxIrecv || gMaxIrecv==-1){

            IncomingMessage * msg_ptr = new IncomingMessage();
//            IncomingMessage * msg_ptr = gIncomingRecvAsync.back();
            msg_ptr->meta = meta;
            msg_ptr->remote_ptr = pRemote_ptr;
            msg_ptr->msg_size = pMsg_size;


            //allocate receive buffer
            bool success = msg_ptr->AllocLocal();
            if(success){
              gIncomingRecvAsync.push_back( msg_ptr );

              msg_ptr->AsyncGet();
              //      logfileptr->OFS()<<gIncomingRecvAsync.size()<<" vs "<<gMaxIrecv<<std::endl;
              //add the function to the async queue
              //      upcxx::async(A.iam)(Aggregate_Compute_Async,Aptr,j,RemoteAggregate, async_copy_event,tstart);
              asyncComm = true;
            }
            else{
              delete msg_ptr;
            }
        }

        if(!asyncComm){
          IncomingMessage * msg_ptr = new IncomingMessage() ;
          //gIncomingRecv.push_back(new IncomingMessage() );
          //IncomingMessage * msg_ptr = gIncomingRecv.back();
          msg_ptr->remote_ptr = pRemote_ptr;
          msg_ptr->msg_size = pMsg_size;
          msg_ptr->meta = meta;
          //gIncomingRecv.push(msg_ptr);
          gIncomingRecv.push_back(msg_ptr);
        }

      }
    SYMPACK_TIMER_STOP(RCV_ASYNC);

  }

  inline std::list< IncomingMessage * >::iterator TestAsyncIncomingMessage(){
    //scope_timer(a,TEST_ASYNC);
    auto it = gIncomingRecvAsync.end();
    if(!gIncomingRecvAsync.empty()){
      //find if there is some finished async comm
      it = gIncomingRecvAsync.begin();
      for(; it!=gIncomingRecvAsync.end();++it){
        if( (*it)->IsDone() /*&& (*it)->IsAsync()*/ ){
//logfileptr->OFS()<<"ASYNC COMM DONE"<<std::endl;
          break;
        }
      }
    }
    return it;
  }


}



#ifdef NO_INTRA_PROFILE
#if defined (SPROFILE)
#define SYMPACK_TIMER_START(a) SYMPACK_FSTART(a);
#define SYMPACK_TIMER_STOP(a) SYMPACK_FSTOP(a);
#endif
#endif





#endif
