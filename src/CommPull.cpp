#include "ngchol/CommPull.hpp"
#include "ngchol/CommTypes.hpp"

namespace LIBCHOLESKY{
  std::list< IncomingMessage * > gIncomingRecv;
  std::list< IncomingMessage * > gIncomingRecvAsync;
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
        delete local_ptr;
      }
    }

    int IncomingMessage::Sender(){
      return remote_ptr.where();
    }

    void IncomingMessage::Wait(){
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

    bool IncomingMessage::IsDone(){
      if(event_ptr!=NULL){
        return event_ptr->isdone();
      }
      else{
        return isDone;
      } 
    }

    bool IncomingMessage::IsAsync(){
      return (event_ptr==NULL);
    }

    void IncomingMessage::AllocLocal(){
      local_ptr = (char *)malloc(msg_size);
    }


    char * IncomingMessage::GetLocalPtr(){
      return (char*)local_ptr;
    }
    
    upcxx::global_ptr<char> IncomingMessage::GetRemotePtr(){
      return remote_ptr;
    }
    








}
