
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
//#include "sympack/CommTypes.hpp"





#ifdef NO_INTRA_PROFILE
#if defined (SPROFILE)
#define SYMPACK_TIMER_START(a)
#define SYMPACK_TIMER_STOP(a)
#endif
#endif

namespace symPACK{
  //extern double maxWaitT;
  //extern double maxAWaitT;
#ifdef NEW_UPCXX
  extern std::list< upcxx::future<> > gFutures;
#endif


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
//      int toto = upcxx::rank_me();
//      int titi = 0;
//      MPI_Comm_rank(symPACK::world_comm,&titi);
//gdb_lock();
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


  extern int last_key;
  extern std::map<int,int> async_barriers;

  bool barrier_done(int id);
  int get_barrier_id(int np);
  void signal_exit(int barrier_id, int np);
#ifdef NEW_UPCXX
   void barrier_wait(int barrier_id, const RankGroup & group);
#else
  void barrier_wait(int barrier_id);
#endif

  inline bool barrier_done(int id){
    return async_barriers[id]==0;
  }

  inline int get_barrier_id(int np){
    int id = last_key++;
    auto it = async_barriers.find(id);
    if(it ==async_barriers.end()){
      async_barriers[id] = np;
    }
    return id;
  }

  inline void signal_exit_am(int barrier_id,int np)
  {
    auto it = async_barriers.find(barrier_id);
    if(it ==async_barriers.end()){
      async_barriers[barrier_id] = np;
    }
    async_barriers[barrier_id]--;
  }

  inline void signal_exit(int barrier_id, int np)
  {
#ifdef NEW_UPCXX
    upcxx::future<> f;
    for (int i = 0; i < np; i++) {
      f = upcxx::when_all(f, upcxx::rpc(i,[](int barrier_id,int np)
            {
            auto it = async_barriers.find(barrier_id);
            if(it ==async_barriers.end()){
            async_barriers[barrier_id] = np;
            }
            async_barriers[barrier_id]--;
            },barrier_id,np) );
    }
    f.wait();
#else
    for (int i = 0; i < np; i++) {
      upcxx::async(i)(signal_exit_am,barrier_id,np);
    }

    //make sure we don't have anything outgoing in flight anymore 
    upcxx::async_wait();
#endif
  }

  inline void signal_exit(int barrier_id, const RankGroup & group)
  {
#ifdef NEW_UPCXX
//////    upcxx::future<> f;
//////    for (int i = 0; i < group.size(); i++) {
//////      int dest = group.L2G(i);
//////      f = upcxx::when_all(f, upcxx::rpc(dest, rpc_signal,barrier_id,group.size()) );
//////    }
//////    f.wait();
////
////    auto iam = upcxx::rank_me();
////    auto giam = group.L2G(iam);
////    auto rpc_signal = [](int barrier_id,int np) {
////      auto it = async_barriers.find(barrier_id);
////      if(it ==async_barriers.end()){
////        async_barriers[barrier_id] = np;
////      }
////      async_barriers[barrier_id]--;
////    };
////
////double tstart = get_time();
////    rpc_signal(barrier_id, group.size());
////
////    std::list< upcxx::future<> > fut;
////    for (int i = 0; i < group.size(); i++) {
////      int dest = group.L2G(i);
////      if(dest!=giam){
//////      fut.push_back(upcxx::rpc(dest, rpc_signal,barrier_id,group.size()));
////
////      fut.push_back(upcxx::rpc(dest, [](int barrier_id,int np) {
////      auto it = async_barriers.find(barrier_id);
////      if(it ==async_barriers.end()){
////        async_barriers[barrier_id] = np;
////      }
////      async_barriers[barrier_id]--;
////    },barrier_id,group.size()));
////
////
////      }
////    }
////double tstop = get_time();
////logfileptr->OFS()<<"launching rpcs time: "<<tstop-tstart<<std::endl;
////
////tstart = get_time();
////    while(!fut.empty()){
////      auto it = fut.begin();
////      while(it!=fut.end()){
////        if(!it->ready()){ break; }
////        it++;
////      }
////      fut.erase(fut.begin(),it);
////      
////      if(!fut.empty()){ upcxx::progress();}
////    }
////tstop = get_time();
////logfileptr->OFS()<<"barrier_wait progress time: "<<tstop-tstart<<std::endl;
#else
double tstart = get_time();
    for (int i = 0; i < group.size(); i++) {
      int dest = group.L2G(i);
      upcxx::async(dest)(signal_exit_am,barrier_id,group.size());
    }
double tstop = get_time();
logfileptr->OFS()<<"launching rpcs time: "<<tstop-tstart<<std::endl;
    //make sure we don't have anything outgoing in flight anymore
tstart = get_time();
    upcxx::async_wait();
tstop = get_time();
logfileptr->OFS()<<"barrier_wait progress time: "<<tstop-tstart<<std::endl;
#endif
  }

#ifdef NEW_UPCXX
    struct prom_handle{
      upcxx::promise<> prom;
      prom_handle(int np){
        prom.require_anonymous(np);
      }
    };
   inline void barrier_wait(int barrier_id, const RankGroup & group){
#if 0
    //upcxx::discharge();
    int np = async_barriers[barrier_id];
    upcxx::promise<> count;
    count.require_anonymous(np);
    upcxx::dist_object<int> dbarrier_id(barrier_id);
    int me = group.G2L(upcxx::rank_me());
    for ( int p = 0; p < np; p++ ){
      int dest = (me+1+p)%np;
      int pdest= group.L2G(dest);
      dbarrier_id.fetch(pdest).then( [&count,barrier_id](int val){bassert(val==barrier_id); count.fulfill_anonymous(1);});
    }
    count.finalize().wait();
#else
    int np = async_barriers[barrier_id];
    upcxx::promise<> * prom_ptr = new upcxx::promise<>();
    prom_ptr->require_anonymous(np);
    upcxx::dist_object<upcxx::promise<> * > dprom( prom_ptr );

    int me = group.G2L(upcxx::rank_me());
    upcxx::future<> fut = upcxx::make_future();
    for ( int p = 0; p < np; p++ ){
      int dest = (me+1+p)%np;
      int pdest= group.L2G(dest);
      fut = upcxx::when_all(fut,upcxx::rpc(pdest,[](upcxx::dist_object<upcxx::promise<> * > & pprom){(*pprom)->fulfill_anonymous(1);},dprom));
    }
    fut = upcxx::when_all(fut, prom_ptr->finalize() );
    fut.wait();
    delete ( prom_ptr );
#endif
   }
#else
   inline void barrier_wait(int barrier_id){
    while( !barrier_done(barrier_id) ){
      upcxx::advance(); 
    }
  }
#endif





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
#ifdef NEW_UPCXX
      upcxx::future<> f_get; 
      bool async_get;
#else
      upcxx::event * event_ptr;
#endif
      //char * local_ptr;
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
      //char * GetLocalPtr();
      std::shared_ptr<char> & GetLocalPtr();

      //void SetLocalPtr(char * ptr,bool ownStorage = true);
      void SetLocalPtr(std::shared_ptr<char> & ptr,bool ownStorage = true);
      virtual size_t Size(){return msg_size;}
      void AsyncGet();
#ifdef NEW_UPCXX
//      void DeallocRemote(std::list< upcxx::future<> > & pFutures);
      void DeallocRemote(upcxx::dist_object<int> & remDealloc);
#else
      void DeallocRemote();
#endif
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
        //local_ptr = amsg->GetLocalPtr();

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


  //extern std::priority_queue< IncomingMessage * ,  std::vector<IncomingMessage *>, MSGCompare > gIncomingRecv;
  extern std::list< IncomingMessage * > gIncomingRecv;

  extern std::list< IncomingMessage * > gIncomingRecvAsync;
  extern std::list< IncomingMessage * > gIncomingRecvLocal;
  extern int gMaxIrecv;
  //extern SupernodalMatrixBase * gSuperMatrixPtr;

#ifdef NEW_UPCXX
  //upcxx::future<> signal_data(upcxx::global_ptr<char> local_ptr, size_t pMsg_size, int dest, MsgMetadata & meta);
  void signal_data(upcxx::global_ptr<char> local_ptr, size_t pMsg_size, int dest, MsgMetadata & meta);
  upcxx::future<> remote_delete(upcxx::global_ptr<char> pRemote_ptr);
#else
  void signal_data(upcxx::global_ptr<char> local_ptr, size_t pMsg_size, int dest, MsgMetadata & meta);
  void remote_delete(upcxx::global_ptr<char> pRemote_ptr);
#endif

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


#ifdef NEW_UPCXX
  inline void signal_data(upcxx::global_ptr<char> local_ptr, size_t pMsg_size, int dest, MsgMetadata & meta){
    scope_timer(a,SIGNAL_DATA);
    upcxx::future<> f_signal;
#ifdef SP_THREADS
    if(Multithreading::NumThread>1){
      throw std::runtime_error("Multithreading is not yet supported in symPACK with the new version of UPCXX");
      //std::lock_guard<upcxx_mutex_type> lock(upcxx_mutex);
      //upcxx::async(dest)(rcv_async,local_ptr,pMsg_size,meta);
    }
    else
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


  inline upcxx::future<> remote_delete(upcxx::global_ptr<char> pRemote_ptr){
    scope_timer(a,REMOTE_DELETE);
    upcxx::future<> f_delete;
    auto dest = pRemote_ptr.where();
#ifdef SP_THREADS
    if(Multithreading::NumThread>1){
      throw std::runtime_error("Multithreading is not yet supported in symPACK with the new version of UPCXX");
      //throw an exception as multithreading is not supported yet in upcxx
      //std::lock_guard<upcxx_mutex_type> lock(upcxx_mutex);
      //upcxx::async(dest)(dealloc_async,pRemote_ptr);
    }
    else
#endif
    {
      f_delete = upcxx::rpc(dest,
          [](upcxx::global_ptr<char> ptr){
          upcxx::deallocate(ptr);
          }
          ,pRemote_ptr);
    }
    return f_delete;
  }

#else
  inline void signal_data(upcxx::global_ptr<char> local_ptr, size_t pMsg_size, int dest, MsgMetadata & meta){
    scope_timer(a,SIGNAL_DATA);
#ifdef SP_THREADS
    if(Multithreading::NumThread>1){
      std::lock_guard<upcxx_mutex_type> lock(upcxx_mutex);
      upcxx::async(dest)(rcv_async,local_ptr,pMsg_size,meta);
    }
    else
#endif
    {
      upcxx::async(dest)(rcv_async,local_ptr,pMsg_size,meta);
    }
  }


  inline void remote_delete(upcxx::global_ptr<char> pRemote_ptr){
    scope_timer(a,REMOTE_DELETE);
    int dest = pRemote_ptr.where();
#ifdef SP_THREADS
    if(Multithreading::NumThread>1){
      std::lock_guard<upcxx_mutex_type> lock(upcxx_mutex);
      upcxx::async(dest)(dealloc_async,pRemote_ptr);
    }
    else
#endif
    {
      upcxx::async(dest)(dealloc_async,pRemote_ptr);
    }
  }
#endif

  inline void dealloc_async(upcxx::global_ptr<char> ptr){
        upcxx::deallocate(ptr);
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





#endif //_COMM_PULL_DECL_HPP_

