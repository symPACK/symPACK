/*
   Copyright (c) 2016 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.  

Author: Mathias Jacquelin

This file is part of symPACK. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to Lawrence Berkeley National
Laboratory, without imposing a separate written license agreement for such
Enhancements, then you hereby grant the following license: a non-exclusive,
royalty-free perpetual license to install, use, modify, prepare derivative
works, incorporate into other computer software, distribute, and sublicense
such enhancements or derivative works thereof, in binary and source code form.
*/
#include "sympack/Environment.hpp"
#include "sympack/CommPull.hpp"
#include "sympack/CommTypes.hpp"
//#include "sympack/SupernodalMatrixBase.hpp"
#include "sympack/timer.hpp"

#define USE_LOCAL_ALLOCATE


namespace symPACK{

#ifdef SP_THREADS
  upcxx_mutex_type upcxx_mutex;
#endif

  std::list< IncomingMessage * > gIncomingRecv;
  //std::priority_queue< IncomingMessage *, std::vector<IncomingMessage *>, MSGCompare > gIncomingRecv;
  std::list< IncomingMessage * > gIncomingRecvAsync;
  std::list< IncomingMessage * > gIncomingRecvLocal;
  //SupernodalMatrixBase * gSuperMatrixPtr = NULL;

  int gMaxIrecv = 0;

  size_t gVolComm=0;
  size_t gNumMsg=0;


  BackupBuffer::BackupBuffer(){
    inUse = false;
    local_ptr = NULL;
    size = 0;
  }
  BackupBuffer::~BackupBuffer(){
    DeallocLocal();
  }

  bool BackupBuffer::AllocLocal(size_t psize){
    if(!Allocated()){
      local_ptr=NULL;
#ifndef USE_LOCAL_ALLOCATE
      local_ptr = (char *)malloc(psize);
#else
      //TODO replace this by a upcxx::allocate
      upcxx::global_ptr<char> tmp = upcxx::allocate<char>(upcxx::myrank(),psize);
      local_ptr=(char*)tmp; 
#endif
      if(local_ptr!=NULL){
        size = psize;
      }
      else{
        throw MemoryAllocationException(psize);
      }
      return local_ptr!=NULL;
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
      //TODO use upcxx::deallocate
      upcxx::global_ptr<char> tmp(local_ptr);
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
      return NULL;
    }
  }

  void BackupBuffer::ReleasePtr(){
    inUse = false;
  }




  IncomingMessage::IncomingMessage(){
    event_ptr=NULL;
    task_ptr =NULL;
    local_ptr=nullptr;
    isDone = false;
    isLocal = false;
    remoteDealloc = false;

    allocated = false;
    ownLocalStorage = false;

#ifdef SP_THREADS
//    references=0;
#endif
  }

  IncomingMessage::~IncomingMessage(){
    //assert(IsDone());
    if(event_ptr!=NULL){
#ifdef SP_THREADS
      if(Multithreading::NumThread>1){
        std::lock_guard<upcxx_mutex_type> lock(upcxx_mutex);
        delete event_ptr;
      }
      else
#endif
        delete event_ptr;
    }
    if(task_ptr!=NULL){
      delete task_ptr;
    }

    DeallocLocal();
  }

  void IncomingMessage::AsyncGet(){
#ifdef SP_THREADS
    if(Multithreading::NumThread>1){
      std::lock_guard<upcxx_mutex_type> lock(upcxx_mutex);
      assert(event_ptr==NULL);
      event_ptr = new upcxx::event;
      upcxx::async_copy(remote_ptr,upcxx::global_ptr<char>(GetLocalPtr().get()),msg_size,event_ptr);
    }
    else
#endif
    {
      assert(event_ptr==NULL);
      event_ptr = new upcxx::event;
      upcxx::async_copy(remote_ptr,upcxx::global_ptr<char>(GetLocalPtr().get()),msg_size,event_ptr);
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
    else if(event_ptr!=NULL){
#ifdef SP_THREADS
      if(Multithreading::NumThread>1){
        std::lock_guard<upcxx_mutex_type> lock(upcxx_mutex);
        //TODO wait is not necessary if calling async_try/isdone
        event_ptr->wait();
        assert(event_ptr->isdone());
        delete event_ptr;
        event_ptr = NULL;
        isDone = true;
        success = true;
      }
      else
#endif
      {
        //TODO wait is not necessary if calling async_try/isdone
        event_ptr->wait();
        assert(event_ptr->isdone());
        delete event_ptr;
        event_ptr = NULL;
        isDone = true;
        success = true;
      }
    }
    else{
      //allocate receive buffer
      success = AllocLocal();
      if(success){
#ifdef SP_THREADS
        if(Multithreading::NumThread>1){
          std::lock_guard<upcxx_mutex_type> lock(upcxx_mutex);
          upcxx::copy(remote_ptr,upcxx::global_ptr<char>(GetLocalPtr().get()),msg_size);
          isDone = true;
        }
        else
#endif
        {
          upcxx::copy(remote_ptr,upcxx::global_ptr<char>(GetLocalPtr().get()),msg_size);
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

  void IncomingMessage::DeallocRemote(){
    if(!remoteDealloc){
//#ifndef NDEBUG
//      logfileptr->OFS()<<"Deleting message from "<<meta.src<<" to "<<meta.tgt<<std::endl;
//#endif
      remote_delete(GetRemotePtr());
      remoteDealloc=true;
    }
  }
  void IncomingMessage::DeallocLocal(){
    if(allocated && ownLocalStorage){
      if(!isLocal && local_ptr!=nullptr){
#ifndef USE_LOCAL_ALLOCATE
//        delete local_ptr;
#else

#if 1

//#ifndef NDEBUG
//      logfileptr->OFS()<<"Deleting message from "<<meta.src<<" to "<<meta.tgt<<std::endl;
//#endif
      //upcxx::global_ptr<char> tmp(local_ptr);
//      upcxx::global_ptr<char> tmp(local_ptr.get());
//      remote_delete(tmp);
      local_ptr = nullptr;

#else
#ifdef SP_THREADS
        if(Multithreading::NumThread>1){
          std::lock_guard<upcxx_mutex_type> lock(upcxx_mutex);
          //TODO use upcxx::deallocate
          upcxx::global_ptr<char> tmp(local_ptr);
          upcxx::deallocate(tmp);
        }
        else
#endif
        {
          //TODO use upcxx::deallocate
          upcxx::global_ptr<char> tmp(local_ptr);
          upcxx::deallocate(tmp);
        }
#endif
#endif
      }
    }
  }


  bool IncomingMessage::IsDone(){
    scope_timer(a,IN_MSG_ISDONE);
    if(event_ptr!=NULL){
#ifdef SP_THREADS
      if(Multithreading::NumThread>1){
        std::lock_guard<upcxx_mutex_type> lock(upcxx_mutex);
        isDone = event_ptr->async_try();
        return isDone; 
      }
      else
#endif
      {
        isDone = event_ptr->async_try();
        return isDone; 
      }
    }
    else{
      return isDone;
    } 
  }

  bool IncomingMessage::IsAsync(){
    return (event_ptr==NULL);
  }

  bool IncomingMessage::AllocLocal(){
    if(!allocated){
      local_ptr=nullptr;
#ifndef USE_LOCAL_ALLOCATE
      //local_ptr = (char *)malloc(msg_size);
      local_ptr=std::shared_ptr<char>((char *)malloc(msg_size)); 
#else
      upcxx::global_ptr<char> tmp;
#ifdef SP_THREADS
      if(Multithreading::NumThread>1){
        std::lock_guard<upcxx_mutex_type> lock(upcxx_mutex);
        //tmp = upcxx::allocate<char>(upcxx::myrank(),msg_size);
        //local_ptr=(char*)tmp; 
        local_ptr=std::shared_ptr<char>( (char *)upcxx::allocate<char>(upcxx::myrank(),msg_size), [=](char * ptr){
      upcxx::global_ptr<char> tmp(ptr);
      remote_delete(tmp);
            }); 
      }
      else
#endif
      {
//        tmp = upcxx::allocate<char>(upcxx::myrank(),msg_size);
//        local_ptr=(char*)tmp; 
        local_ptr=std::shared_ptr<char>((char*)upcxx::allocate<char>(upcxx::myrank(),msg_size), [=](char * ptr){
      upcxx::global_ptr<char> tmp(ptr);
      remote_delete(tmp);
            }); 
      }
#endif

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


  //char * IncomingMessage::GetLocalPtr(){
  //  return (char*)local_ptr;
  //}

  std::shared_ptr<char> & IncomingMessage::GetLocalPtr(){
    return std::ref(local_ptr);
  }

//  void IncomingMessage::SetLocalPtr(char * ptr,bool ownStorage){
//    local_ptr = ptr;
//    allocated = true;
//    ownLocalStorage = ownStorage;
//  }


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


