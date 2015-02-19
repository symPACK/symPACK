#ifndef _COMM_PULL_DECL_HPP_
#define _COMM_PULL_DECL_HPP_


#include <list>
#include <vector>

#include <upcxx.h>



#include "ngchol/Environment.hpp"
//#include "ngchol/CommTypes.hpp"

namespace LIBCHOLESKY{

struct SnodeUpdateFB;

class IncomingMessage{
  public:
    upcxx::global_ptr<char> remote_ptr;
    size_t msg_size;
    upcxx::event * event_ptr;
    char * local_ptr;
    SnodeUpdateFB * task_ptr;
    bool isDone; 
    IncomingMessage(){
      event_ptr=NULL;
      task_ptr =NULL;
      local_ptr=NULL;
      isDone = false;
    }

    ~IncomingMessage(){
      //assert(IsDone());
      if(event_ptr!=NULL){
        delete event_ptr;
      }
      if(task_ptr!=NULL){
        delete task_ptr;
      }
      if(local_ptr!=NULL){
        delete local_ptr;
      }
    }

    int Sender(){
      return remote_ptr.where();
    }

    void Wait(){
      if(event_ptr!=NULL){
        event_ptr->wait();
        assert(event_ptr->isdone());
        delete event_ptr;
        event_ptr = NULL;
        isDone = true;
      }
      else{
        //allocate receive buffer
        AllocLocal();

        upcxx::copy(remote_ptr,upcxx::global_ptr<char>(GetLocalPtr()),msg_size);

        isDone = true;
      }
    }

    bool IsDone(){
      if(event_ptr!=NULL){
        return event_ptr->isdone();
      }
      else{
        return isDone;
      } 
    }

    bool IsAsync(){
      return (event_ptr==NULL);
    }

    void AllocLocal(){
      local_ptr = (char *)malloc(msg_size);
    }


    char * GetLocalPtr(){
      return (char*)local_ptr;
    }
    
    upcxx::global_ptr<char> GetRemotePtr(){
      return remote_ptr;
    }
    

};

  extern std::list< IncomingMessage * > gIncomingRecv;
  extern std::list< IncomingMessage * > gIncomingRecvAsync;
  extern int gMaxIrecv;

  void signal_data(char * local_ptr, size_t pMsg_size, int dest);
  void rcv_async(upcxx::global_ptr<char> pRemote_ptr, size_t pMsg_size);

  inline void signal_data(char * local_ptr, size_t pMsg_size, int dest){
      upcxx::async(dest)(rcv_async,upcxx::global_ptr<char>(local_ptr),pMsg_size);
  }


  inline void remote_delete(upcxx::global_ptr<char> pRemote_ptr){
      if(upcxx::myrank()!=pRemote_ptr.where()){
//        upcxx::async(pRemote_ptr.where())(remote_delete,pRemote_ptr);
cout<<"Performing remote delete on P"<<pRemote_ptr.where()<<endl;
        upcxx::deallocate(pRemote_ptr);
      }
      else{
        upcxx::deallocate(pRemote_ptr);
      }
  }

  inline void rcv_async(upcxx::global_ptr<char> pRemote_ptr, size_t pMsg_size){

    //if we still have async buffers
    if(gIncomingRecvAsync.size() < gMaxIrecv){
      gIncomingRecvAsync.push_back(new IncomingMessage() );
      IncomingMessage * msg_ptr = gIncomingRecvAsync.back();
      msg_ptr->event_ptr = new upcxx::event;
      msg_ptr->remote_ptr = pRemote_ptr;
      msg_ptr->msg_size = pMsg_size;
      //allocate receive buffer
      msg_ptr->AllocLocal();

      upcxx::async_copy(pRemote_ptr,upcxx::global_ptr<char>(msg_ptr->GetLocalPtr()),pMsg_size,msg_ptr->event_ptr);

      logfileptr->OFS()<<gIncomingRecvAsync.size()<<" vs "<<gMaxIrecv<<endl;
      //add the function to the async queue
      //      upcxx::async(A.iam)(Aggregate_Compute_Async,Aptr,j,RemoteAggregate, async_copy_event,tstart);
    }
    else{
      gIncomingRecv.push_back(new IncomingMessage() );
      IncomingMessage * msg_ptr = gIncomingRecv.back();
      msg_ptr->remote_ptr = pRemote_ptr;
      msg_ptr->msg_size = pMsg_size;
//      //allocate receive buffer
//      msg_ptr->AllocLocal();
//
//      upcxx::copy(pRemote_ptr,upcxx::global_ptr<T>(msg_ptr->GetLocalPtr()),pMsgSize);
      //call the function inline
      //      Aggregate_Compute_Async(Aptr, j, RemoteAggregate, NULL,tstart);
    }
  }
}

#endif
