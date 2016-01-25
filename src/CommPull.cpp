#include "sympack/CommPull.hpp"
#include "sympack/CommTypes.hpp"
#include "sympack/SupernodalMatrixBase.hpp"
#include "sympack/timer.hpp"

#define USE_LOCAL_ALLOCATE


namespace SYMPACK{
  std::list< IncomingMessage * > gIncomingRecv;
  //std::priority_queue< IncomingMessage *, vector<IncomingMessage *>, MSGCompare > gIncomingRecv;
  std::list< IncomingMessage * > gIncomingRecvAsync;
  std::list< IncomingMessage * > gIncomingRecvLocal;
  SupernodalMatrixBase * gSuperMatrixPtr = NULL;

  int gMaxIrecv = 0;

  size_t gVolComm=0;
  size_t gNumMsg=0;

 
    IncomingMessage::IncomingMessage(){
      event_ptr=NULL;
      task_ptr =NULL;
      local_ptr=NULL;
      isDone = false;
      isLocal = false;
      remoteDealloc = false;
    }

    IncomingMessage::~IncomingMessage(){
      //assert(IsDone());
      if(event_ptr!=NULL){
        delete event_ptr;
      }
      if(task_ptr!=NULL){
        delete task_ptr;
      }
      if(!isLocal && local_ptr!=NULL){
#ifndef USE_LOCAL_ALLOCATE
        delete local_ptr;
#else
        //TODO use upcxx::deallocate
        upcxx::global_ptr<char> tmp(local_ptr);
        upcxx::deallocate(tmp);
#endif
      }
    }

    void IncomingMessage::AsyncGet(){
      assert(event_ptr==NULL);
      event_ptr = new upcxx::event;
      upcxx::async_copy(remote_ptr,upcxx::global_ptr<char>(GetLocalPtr()),msg_size,event_ptr);
    }

    int IncomingMessage::Sender(){
      return remote_ptr.where();
    }

    void IncomingMessage::Wait(){
      scope_timer(a,IN_MSG_WAIT);
      if(isLocal){
        isDone = true;
      }
      else if(event_ptr!=NULL){
        //TODO wait is not necessary if calling async_try/isdone
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
#ifdef PROFILE_COMM
      gVolComm+= msg_size;
      gNumMsg++;
#endif
    }

    bool IncomingMessage::IsLocal(){
      return isLocal;
    }

    void IncomingMessage::DeallocRemote(){
        if(!remoteDealloc){
          remote_delete(GetRemotePtr());
          remoteDealloc=true;
        }
    }

    bool IncomingMessage::IsDone(){
      scope_timer(a,IN_MSG_ISDONE);
      if(event_ptr!=NULL){
        //return event_ptr->isdone();
        return event_ptr->async_try();
        //TODO also look at event_ptr async_try because it calls "progress"
      }
      else{
        return isDone;
      } 
    }

    bool IncomingMessage::IsAsync(){
      return (event_ptr==NULL);
    }

    void IncomingMessage::AllocLocal(){
#ifndef USE_LOCAL_ALLOCATE
      local_ptr = (char *)malloc(msg_size);
#else
      //TODO replace this by a upcxx::allocate
      upcxx::global_ptr<char> tmp = upcxx::allocate<char>(iam,msg_size);
      local_ptr=(char*)tmp; 
#endif
    }


    char * IncomingMessage::GetLocalPtr(){
      return (char*)local_ptr;
    }
    
    upcxx::global_ptr<char> IncomingMessage::GetRemotePtr(){
      return remote_ptr;
    }
    


#ifdef USE_LOCAL_ALLOCATE
#undef USE_LOCAL_ALLOCATE
#endif






}


