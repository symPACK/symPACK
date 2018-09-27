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
#include <map>

#define USE_LOCAL_ALLOCATE


namespace symPACK{

  double maxWaitT = 0.0;
  double maxAWaitT = 0.0;

#ifdef NEW_UPCXX
  std::list< upcxx::future<> > gFutures;
#endif

  int last_key = 0;
  std::map<int,int> async_barriers;

   
    
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
#ifdef NEW_UPCXX
      upcxx::global_ptr<char> tmp = upcxx::allocate<char>(psize);
      local_ptr=(char*)tmp.local(); 
#else
      upcxx::global_ptr<char> tmp = upcxx::allocate<char>(upcxx::myrank(),psize);
      local_ptr=(char*)tmp; 
#endif
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
#if UPCXX_VERSION >= 20180305
      upcxx::global_ptr<char> tmp = upcxx::to_global_ptr(local_ptr);
#else
      upcxx::global_ptr<char> tmp(local_ptr);
#endif
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

#ifdef NEW_UPCXX
    async_get = false;
#else
    event_ptr=NULL;
#endif
    task_ptr =NULL;
    local_ptr=nullptr;
    isDone = false;
    isLocal = false;
    remoteDealloc = false;

    allocated = false;
    ownLocalStorage = false;
  }

  IncomingMessage::~IncomingMessage(){
    //assert(IsDone());
#ifdef NEW_UPCXX
    if(async_get){
#ifdef SP_THREADS
      if(Multithreading::NumThread>1){
        throw std::runtime_error("Multithreading is not yet supported in symPACK with the new version of UPCXX");
      }
      else
#endif
      {
        f_get.wait();
      }
    }
#else
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
#endif
    if(task_ptr!=NULL){
      delete task_ptr;
    }

    DeallocLocal();
  }

  void IncomingMessage::AsyncGet(){
#ifdef NEW_UPCXX
#ifdef SP_THREADS
    if(Multithreading::NumThread>1){
      throw std::runtime_error("Multithreading is not yet supported in symPACK with the new version of UPCXX");
    }
    else
#endif
    {
      f_get = upcxx::rget(remote_ptr,GetLocalPtr().get(),msg_size);
      async_get = true;
    }

#else
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
#endif
  }

  int IncomingMessage::Sender(){
    return remote_ptr.where();
  }

  bool IncomingMessage::Wait(){
    scope_timer(a,IN_MSG_WAIT);
#ifdef NEW_UPCXX
    bool success = false;
    if(isLocal){
      isDone = true;
      success = true;
    }
    else if(async_get){
#ifdef SP_THREADS
      if(Multithreading::NumThread>1){
        throw std::runtime_error("Multithreading is not yet supported in symPACK with the new version of UPCXX");
      }
      else
#endif
      {
double tstart = get_time();
        f_get.wait();
double tstop = get_time();
maxAWaitT = std::max(maxAWaitT, tstop-tstart);
        async_get = false;
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
          throw std::runtime_error("Multithreading is not yet supported in symPACK with the new version of UPCXX");
        }
        else
#endif
        {

double tstart = get_time();
          upcxx::rget(remote_ptr,GetLocalPtr().get(),msg_size).wait();
double tstop = get_time();
maxWaitT = std::max(maxWaitT, tstop-tstart);

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
#else
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
double tstart = get_time();
        event_ptr->wait();
double tstop = get_time();
maxAWaitT = std::max(maxAWaitT, tstop-tstart);
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
double tstart = get_time();
        event_ptr->wait();
double tstop = get_time();
maxAWaitT = std::max(maxAWaitT, tstop-tstart);
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

double tstart = get_time();
          upcxx::copy(remote_ptr,upcxx::global_ptr<char>(GetLocalPtr().get()),msg_size);
double tstop = get_time();
maxWaitT = std::max(maxWaitT, tstop-tstart);
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
#endif
    return success;
  }

  bool IncomingMessage::IsLocal(){
    return isLocal;
  }

#ifdef NEW_UPCXX
  void IncomingMessage::DeallocRemote( std::list< upcxx::future<> > & pFutures){
    if(!remoteDealloc){
      auto ptr = GetRemotePtr();
      auto pdest = ptr.where();
      auto fut = upcxx::rpc(pdest,[ptr](){upcxx::deallocate(ptr);});
      pFutures.push_back(fut);
    }
  }
#else
  void IncomingMessage::DeallocRemote(){
    if(!remoteDealloc){
      upcxx::deallocate(GetRemotePtr());
      remoteDealloc=true;
    }
  }
#endif
  void IncomingMessage::DeallocLocal(){
    if(allocated && ownLocalStorage){
      if(!isLocal && local_ptr!=nullptr){
      local_ptr = nullptr;
      }
    }
  }


  bool IncomingMessage::IsDone(){
    scope_timer(a,IN_MSG_ISDONE);
#ifdef NEW_UPCXX
    if(async_get){
#ifdef SP_THREADS
      if(Multithreading::NumThread>1){
        throw std::runtime_error("Multithreading is not yet supported in symPACK with the new version of UPCXX");
        return isDone; 
      }
      else
#endif
      {
        isDone = f_get.ready();
        if(!isDone){
          upcxx::progress();
          isDone = f_get.ready();
        }
        return isDone; 
      }
    }
    else{
      return isDone;
    }
#else
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
#endif
  }

  bool IncomingMessage::IsAsync(){
#ifdef NEW_UPCXX
    return (async_get);
#else
    return (event_ptr==NULL);
#endif
  }

  bool IncomingMessage::AllocLocal(){
#ifdef NEW_UPCXX
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
#else
    if(!allocated){
      local_ptr=nullptr;
#ifndef USE_LOCAL_ALLOCATE
      local_ptr=std::shared_ptr<char>( new char[msg_size], [=](char * ptr){ delete [] ptr; });
#else
      upcxx::global_ptr<char> tmp;
#ifdef SP_THREADS
      if(Multithreading::NumThread>1){
        std::lock_guard<upcxx_mutex_type> lock(upcxx_mutex);
        local_ptr=std::shared_ptr<char>( (char *)upcxx::allocate<char>(upcxx::myrank(),msg_size), [=](char * ptr){
            upcxx::global_ptr<char> tmp(ptr);
            upcxx::deallocate(tmp);
            }); 
      }
      else
#endif
      {
        local_ptr=std::shared_ptr<char>((char*)upcxx::allocate<char>(upcxx::myrank(),msg_size), [=](char * ptr){
            upcxx::global_ptr<char> tmp(ptr);
            upcxx::deallocate(tmp);
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
#endif
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


