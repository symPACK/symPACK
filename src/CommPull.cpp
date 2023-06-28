#include "sympack/Environment.hpp"
#include "sympack/CommPull.hpp"
#include "sympack/CommTypes.hpp"
#include "sympack/timer.hpp"
#include <map>

#define USE_LOCAL_ALLOCATE


namespace symPACK{



   
    
    #ifdef SP_THREADS
      upcxx_mutex_type upcxx_mutex;
#endif

  std::list< IncomingMessage * > gIncomingRecv;
  std::list< IncomingMessage * > gIncomingRecvAsync;
  std::list< IncomingMessage * > gIncomingRecvLocal;

  int gMaxIrecv = 0;

  size_t gVolComm=0;
  size_t gNumMsg=0;


  BackupBuffer::BackupBuffer(){
    inUse = false;
    local_ptr = nullptr;
    size = 0;
  }
  BackupBuffer::~BackupBuffer(){
    DeallocLocal();
  }

  bool BackupBuffer::AllocLocal(size_t psize){
    if(!Allocated()){
      local_ptr=nullptr;
#ifndef USE_LOCAL_ALLOCATE
      local_ptr = (char *)malloc(psize);
#else
      //TODO replace this by a upcxx::allocate
      upcxx::global_ptr<char> tmp = upcxx::allocate<char>(psize);
      local_ptr=(char*)tmp.local(); 
#endif
      if(local_ptr!=nullptr){
        size = psize;
      }
      else{
        throw MemoryAllocationException(psize);
      }
      return local_ptr!=nullptr;
    }
    else{
      return false;
    }
  }

  void BackupBuffer::DeallocLocal(){
    if(Allocated()){
#ifndef USE_LOCAL_ALLOCATE
      delete local_ptr;
#else
      upcxx::global_ptr<char> tmp = upcxx::to_global_ptr(local_ptr);
      upcxx::deallocate(tmp);
#endif
    }
  }

  char * BackupBuffer::GetPtr(){
    if(!inUse){
      inUse = true;
      return local_ptr; 
    }
    else{
      return nullptr;
    }
  }

  void BackupBuffer::ReleasePtr(){
    inUse = false;
  }




  IncomingMessage::IncomingMessage(){
    async_get = false;
    task_ptr =nullptr;
    local_ptr=nullptr;
    isDone = false;
    isLocal = false;
    remoteDealloc = false;

    allocated = false;
    ownLocalStorage = false;
  }

  IncomingMessage::~IncomingMessage(){
    if(async_get){
      {
        f_get.wait();
      }
    }
    if(task_ptr!=nullptr){
      delete task_ptr;
    }

    DeallocLocal();
  }

  void IncomingMessage::AsyncGet(){
    {
      f_get = upcxx::rget(remote_ptr,GetLocalPtr().get(),msg_size);
      async_get = true;
    }
  }

  int IncomingMessage::Sender(){
    return remote_ptr.where();
  }

  bool IncomingMessage::Wait(){
    scope_timer(a,IN_MSG_WAIT);
    bool success = false;
    if(isLocal){
      isDone = true;
      success = true;
    }
    else if(async_get){
      {
        f_get.wait();
        async_get = false;
        isDone = true;
        success = true;
      }
    }
    else{
      //allocate receive buffer
      success = AllocLocal();
      if(success){
        {
          upcxx::rget(remote_ptr,GetLocalPtr().get(),msg_size).wait();
          isDone = true;
        }
      }
    }
#ifdef PROFILE_COMM
    if(success){
      gVolComm+= msg_size;
      gNumMsg++;
    }
#endif
    return success;
  }

  bool IncomingMessage::IsLocal(){
    return isLocal;
  }

  void IncomingMessage::DeallocRemote(upcxx::dist_object<int> & remDealloc) {
    if(!remoteDealloc){
      auto ptr = GetRemotePtr();
      auto pdest = ptr.where();
      upcxx::rpc_ff(pdest,[ptr](upcxx::dist_object<int> & dealloc_cnt){ 
          static int init_cnt=*dealloc_cnt;
          static int call_cnt=0;
          call_cnt++;
          bassert(*dealloc_cnt > 0 );
          upcxx::deallocate(ptr); (*dealloc_cnt)--;},remDealloc);
    }
  }

  void IncomingMessage::DeallocLocal(){
    if(allocated && ownLocalStorage){
      if(!isLocal && local_ptr!=nullptr){
      local_ptr = nullptr;
      }
    }
  }


  bool IncomingMessage::IsDone(){
    scope_timer(a,IN_MSG_ISDONE);
    if(async_get){
      {
#if UPCXX_VERSION >= 20230305
	isDone = f_get.is_ready();
#else
        isDone = f_get.ready();
#endif
        if(!isDone){
          upcxx::progress();
#if UPCXX_VERSION >= 20230305
          isDone = f_get.is_ready();
#else
          isDone = f_get.ready();
#endif
        }
        return isDone; 
      }
    }
    else{
      return isDone;
    }
  }

  bool IncomingMessage::IsAsync(){
    return (async_get);
  }

  bool IncomingMessage::AllocLocal(){
    if(!allocated){
      local_ptr=nullptr;
      local_ptr=std::shared_ptr<char>( new char[msg_size], [=](char * ptr){ delete [] ptr; });

      allocated = local_ptr!=nullptr;
      ownLocalStorage = allocated;

      if(local_ptr==nullptr){
        throw MemoryAllocationException(msg_size);
      }
      return local_ptr!=nullptr;
    }
    else{
      return true;
    }
  }

  std::shared_ptr<char> & IncomingMessage::GetLocalPtr(){
    return std::ref(local_ptr);
  }

  void IncomingMessage::SetLocalPtr(std::shared_ptr<char> & ptr,bool ownStorage){
    local_ptr = ptr;
    allocated = true;
    ownLocalStorage = ownStorage;
  }

  upcxx::global_ptr<char> IncomingMessage::GetRemotePtr(){
    return remote_ptr;
  }



#ifdef USE_LOCAL_ALLOCATE
#undef USE_LOCAL_ALLOCATE
#endif
}

