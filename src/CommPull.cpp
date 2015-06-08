#include "ngchol/CommPull.hpp"
#include "ngchol/CommTypes.hpp"
#include "ngchol/SupernodalMatrixBase.hpp"

namespace LIBCHOLESKY{
  std::list< IncomingMessage * > gIncomingRecv;
  std::list< IncomingMessage * > gIncomingRecvAsync;
  SupernodalMatrixBase * gSuperMatrixPtr = NULL;

  int gMaxIrecv = 0;


 
    IncomingMessage::IncomingMessage(){
      event_ptr=NULL;
      task_ptr =NULL;
      local_ptr=NULL;
      isDone = false;
    }

    IncomingMessage::~IncomingMessage(){
      //assert(IsDone());
      if(event_ptr!=NULL){
        delete event_ptr;
      }
      if(task_ptr!=NULL){
        delete task_ptr;
      }
      if(local_ptr!=NULL){
        //TODO use upcxx::deallocate
        //delete local_ptr;
        upcxx::global_ptr<char> tmp(local_ptr);
        upcxx::deallocate(tmp);
      }
    }

    int IncomingMessage::Sender(){
      return remote_ptr.where();
    }

    void IncomingMessage::Wait(){
      if(event_ptr!=NULL){
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
    }

    bool IncomingMessage::IsDone(){
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
      //local_ptr = (char *)malloc(msg_size);
      //TODO replace this by a upcxx::allocate
      upcxx::global_ptr<char> tmp = upcxx::allocate<char>(iam,msg_size);
      local_ptr=(char*)tmp; 
    }


    char * IncomingMessage::GetLocalPtr(){
      return (char*)local_ptr;
    }
    
    upcxx::global_ptr<char> IncomingMessage::GetRemotePtr(){
      return remote_ptr;
    }
    








}
