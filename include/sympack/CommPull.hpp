
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

#ifndef _COMM_PULL_DECL_HPP_
#define _COMM_PULL_DECL_HPP_
#include <list>
#include <vector>
#include <queue>


#include <mutex>
#include <thread>
#include <atomic>
#include <numeric>

#include "sympack/Environment.hpp"


#ifdef NO_INTRA_PROFILE
#if defined (SPROFILE)
#define SYMPACK_TIMER_START(a)
#define SYMPACK_TIMER_STOP(a)
#endif
#endif

namespace symPACK{

  class RankGroup{
    public:
    RankGroup(MPI_Comm & comm){
      int size = 0;
      MPI_Comm_size(comm,&size);

      //resize our translation arrays
      l2g.resize(size);

      //get the MPI_Group from the communicator amd MPI_COMM_WORLD
      MPI_Group group, Wgroup;
      MPI_Comm_group(comm, &group);
      MPI_Comm_group(symPACK::world_comm, &Wgroup);
      std::vector<int> tmp(size);
      std::iota(tmp.begin(),tmp.end(),0);
      //Get the corresponding ranks in MPI_COMM_WORLD
      MPI_Group_translate_ranks(group, size, tmp.data(), Wgroup, l2g.data());

      for(int i = 0; i < size; i++ ){
        g2l[ l2g[i] ] = i;
      } 
    };

    int L2G(const int rank) const { return l2g[rank];}
    int G2L(const int rank) const { return g2l.at(rank);}

    int size() const { return l2g.size(); }
    std::vector<int> l2g;
    std::map<int,int> g2l;
  };


  struct SnodeUpdateFB;

  struct MsgMetadata{
    //sender taskid
    int src;
    //receiver taskid
    int tgt;

    int GIndex;
    size_t id;
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
      bool Allocated(){return local_ptr!=nullptr;}
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
      upcxx::future<> f_get; 
      bool async_get;
      std::shared_ptr<char> local_ptr;

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
      virtual bool IsDone();
      bool IsLocal();
      bool IsAsync();
      bool AllocLocal();
      upcxx::global_ptr<char> GetRemotePtr();
      std::shared_ptr<char> & GetLocalPtr();

      void SetLocalPtr(std::shared_ptr<char> & ptr,bool ownStorage = true);
      virtual size_t Size(){return msg_size;}
      void AsyncGet();
      void DeallocRemote(upcxx::dist_object<int> & remDealloc);
      void DeallocLocal();
  };

  template< typename C >
  class ChainedMessage: public IncomingMessage{
    public:
      std::shared_ptr<C> data;
      std::shared_ptr<IncomingMessage> chainedMsg;
      virtual size_t Size(){return msg_size;}

      virtual bool IsDone(){
        isDone = IncomingMessage::IsDone();
        return isDone;
      }

      ChainedMessage(std::shared_ptr<IncomingMessage> amsg):IncomingMessage(){
        data = nullptr;
        chainedMsg = amsg;

        isDone = amsg->IsDone(); 
        ownLocalStorage = false;
      }

      ChainedMessage(std::shared_ptr<C> adata, std::shared_ptr<IncomingMessage> amsg):IncomingMessage(){
        data = adata;
        chainedMsg = amsg;

        isDone = amsg->IsDone(); 
        ownLocalStorage = false;
      }
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


  extern std::list< IncomingMessage * > gIncomingRecv;

  extern std::list< IncomingMessage * > gIncomingRecvAsync;
  extern std::list< IncomingMessage * > gIncomingRecvLocal;
  extern int gMaxIrecv;

  void signal_data(upcxx::global_ptr<char> local_ptr, size_t pMsg_size, int dest, MsgMetadata & meta);

  void rcv_async(upcxx::global_ptr<char> pRemote_ptr, size_t pMsg_size, MsgMetadata meta);
  void dealloc_async(upcxx::global_ptr<char> ptr);

  inline void rcv_async(upcxx::global_ptr<char> pRemote_ptr, size_t pMsg_size, MsgMetadata meta){
    scope_timer(a,RCV_ASYNC);

    {
      //if we still have async buffers
      bool asyncComm = false;
      if(gIncomingRecvAsync.size() < gMaxIrecv || gMaxIrecv==-1){

        IncomingMessage * msg_ptr = new IncomingMessage();
        msg_ptr->meta = meta;
        msg_ptr->remote_ptr = pRemote_ptr;
        msg_ptr->msg_size = pMsg_size;


        //allocate receive buffer
        bool success = msg_ptr->AllocLocal();
        if(success){
          gIncomingRecvAsync.push_back( msg_ptr );
          msg_ptr->AsyncGet();
          asyncComm = true;
        }
        else{
          delete msg_ptr;
        }
      }

      if(!asyncComm){
        IncomingMessage * msg_ptr = new IncomingMessage() ;
        msg_ptr->remote_ptr = pRemote_ptr;
        msg_ptr->msg_size = pMsg_size;
        msg_ptr->meta = meta;
        gIncomingRecv.push_back(msg_ptr);
      }
    }
  }


  inline void signal_data(upcxx::global_ptr<char> local_ptr, size_t pMsg_size, int dest, MsgMetadata & meta){
    scope_timer(a,SIGNAL_DATA);
#ifdef SP_THREADS
#endif
    {

      upcxx::rpc_ff(dest,   
          [](upcxx::global_ptr<char> pRemote_ptr, size_t pMsg_size, MsgMetadata meta){
          scope_timer(a,RCV_ASYNC);
          //if we still have async buffers
          bool asyncComm = false;
          if(gIncomingRecvAsync.size() < gMaxIrecv || gMaxIrecv==-1){
          IncomingMessage * msg_ptr = new IncomingMessage();
          msg_ptr->meta = meta;
          msg_ptr->remote_ptr = pRemote_ptr;
          msg_ptr->msg_size = pMsg_size;


          //allocate receive buffer
          bool success = msg_ptr->AllocLocal();
          if(success){
          gIncomingRecvAsync.push_back( msg_ptr );
          msg_ptr->AsyncGet();
          asyncComm = true;
          }
          else{
            delete msg_ptr;
          }
          }

          if(!asyncComm){
            IncomingMessage * msg_ptr = new IncomingMessage() ;
            msg_ptr->remote_ptr = pRemote_ptr;
            msg_ptr->msg_size = pMsg_size;
            msg_ptr->meta = meta;
            gIncomingRecv.push_back(msg_ptr);
          }
          }
      ,local_ptr,pMsg_size,meta);
    }
  }

  inline void dealloc_async(upcxx::global_ptr<char> ptr){
        upcxx::deallocate(ptr);
  }


  inline std::list< IncomingMessage * >::iterator TestAsyncIncomingMessage(){
    auto it = gIncomingRecvAsync.end();
    if(!gIncomingRecvAsync.empty()){
      //find if there is some finished async comm
      it = gIncomingRecvAsync.begin();
      for(; it!=gIncomingRecvAsync.end();++it){
        if( (*it)->IsDone()  ){
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





#endif //_COMM_PULL_DECL_HPP_

