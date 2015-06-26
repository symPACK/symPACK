#ifndef _COMM_PULL_DECL_HPP_
#define _COMM_PULL_DECL_HPP_


#include <list>
#include <vector>

#include <upcxx.h>



#include "ngchol/Environment.hpp"
//#include "ngchol/CommTypes.hpp"


#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) TAU_FSTART(a);
#define TIMER_STOP(a) TAU_FSTOP(a);
#endif
#endif




namespace LIBCHOLESKY{

struct SnodeUpdateFB;
class SupernodalMatrixBase;

  struct MsgMetadata{
    //sender taskid
    int src;
    //receiver taskid
    int tgt;
    
    int GIndex;
    //type of message
//    TaskType type;
  };


class IncomingMessage{
  public:
    upcxx::global_ptr<char> remote_ptr;
    size_t msg_size;
    MsgMetadata meta;
    upcxx::event * event_ptr;
    char * local_ptr;
    SnodeUpdateFB * task_ptr;
    bool isDone; 
    IncomingMessage();
    ~IncomingMessage();
    int Sender();
    void Wait();
    bool IsDone();
    bool IsAsync();
    void AllocLocal();
    char * GetLocalPtr();
    size_t Size(){return msg_size;}

    upcxx::global_ptr<char> GetRemotePtr();
};

  extern std::list< IncomingMessage * > gIncomingRecv;
  extern std::list< IncomingMessage * > gIncomingRecvAsync;
  extern int gMaxIrecv;
  extern SupernodalMatrixBase * gSuperMatrixPtr;

  void signal_data(upcxx::global_ptr<char> local_ptr, size_t pMsg_size, int dest, MsgMetadata & meta);
  void rcv_async(upcxx::global_ptr<char> pRemote_ptr, size_t pMsg_size, MsgMetadata meta);

  inline void signal_data(upcxx::global_ptr<char> local_ptr, size_t pMsg_size, int dest, MsgMetadata & meta){
      scope_timer(SIGNAL_DATA);
      upcxx::async(dest)(rcv_async,local_ptr,pMsg_size,meta);
  }


  inline void remote_delete(upcxx::global_ptr<char> pRemote_ptr){
      scope_timer(REMOTE_DELETE);
      if(upcxx::myrank()!=pRemote_ptr.where()){
        //logfileptr->OFS()<<"Performing remote delete on P"<<pRemote_ptr.where()<<endl;
        //upcxx::async(pRemote_ptr.where())(remote_delete,pRemote_ptr);
        upcxx::deallocate(pRemote_ptr);
      }
      else{
        upcxx::deallocate(pRemote_ptr);
      }
  }

  inline void rcv_async(upcxx::global_ptr<char> pRemote_ptr, size_t pMsg_size, MsgMetadata meta){
      scope_timer(RCV_ASYNC);

    //if we still have async buffers
    if(gIncomingRecvAsync.size() < gMaxIrecv){
      gIncomingRecvAsync.push_back(new IncomingMessage() );
      IncomingMessage * msg_ptr = gIncomingRecvAsync.back();
      msg_ptr->meta = meta;
      msg_ptr->event_ptr = new upcxx::event;
      msg_ptr->remote_ptr = pRemote_ptr;
      msg_ptr->msg_size = pMsg_size;
      //allocate receive buffer
      msg_ptr->AllocLocal();

      upcxx::async_copy(pRemote_ptr,upcxx::global_ptr<char>(msg_ptr->GetLocalPtr()),pMsg_size,msg_ptr->event_ptr);

      //      logfileptr->OFS()<<gIncomingRecvAsync.size()<<" vs "<<gMaxIrecv<<endl;
      //add the function to the async queue
      //      upcxx::async(A.iam)(Aggregate_Compute_Async,Aptr,j,RemoteAggregate, async_copy_event,tstart);
    }
    else{
      gIncomingRecv.push_back(new IncomingMessage() );
      IncomingMessage * msg_ptr = gIncomingRecv.back();
      msg_ptr->remote_ptr = pRemote_ptr;
      msg_ptr->msg_size = pMsg_size;
      msg_ptr->meta = meta;
//      //allocate receive buffer
//      msg_ptr->AllocLocal();
//
//      upcxx::copy(pRemote_ptr,upcxx::global_ptr<T>(msg_ptr->GetLocalPtr()),pMsgSize);
      //call the function inline
      //      Aggregate_Compute_Async(Aptr, j, RemoteAggregate, NULL,tstart);
    }

  }

  inline std::list< IncomingMessage * >::iterator TestAsyncIncomingMessage(){
    auto it = gIncomingRecvAsync.end();
    if(!gIncomingRecvAsync.empty()){
      //find if there is some finished async comm
      it = gIncomingRecvAsync.begin();
      for(; it!=gIncomingRecvAsync.end();++it){
        if( (*it)->IsDone() /*&& (*it)->IsAsync()*/ ){
          break;
        }
      }
    }
    return it;
  }


}



#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a)
#define TIMER_STOP(a)
#endif
#endif




#endif
