/// @file run_sparse_fanboth.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2014-01-23
//#define _DEBUG_


#include <mpi.h>

#include <time.h>
#include <omp.h>

#include  "sympack.hpp"

#include  "sympack/sp_blas.hpp"
#include  "sympack/CommTypes.hpp"
#include  "sympack/CommPull.hpp"
#include  "sympack/Ordering.hpp"

extern "C" {
#include "bebop/util/config.h"
#include "bebop/smc/sparse_matrix.h"
#include "bebop/smc/csr_matrix.h"
#include "bebop/smc/csc_matrix.h"
#include "bebop/smc/sparse_matrix_ops.h"

#include "bebop/util/get_options.h"
#include "bebop/util/init.h"
#include "bebop/util/malloc.h"
#include "bebop/util/timer.h"
#include "bebop/util/util.h"
}

#include <upcxx.h>

#define MYSCALAR double


using namespace SYMPACK;


int main(int argc, char **argv) 
{
  upcxx::init(&argc, &argv);
  //MPI_Init(&argc,&argv);

  iam = upcxx::myrank();
  np = upcxx::ranks();

  logfileptr = new LogFile(iam);
  logfileptr->OFS()<<"********* LOGFILE OF P"<<iam<<" *********"<<endl;
  logfileptr->OFS()<<"**********************************"<<endl;

//  //Do a test for the upcxx pull thing
//  SYMPACK::vector<Icomm *> myComms(np,NULL);
//  for(int p = 0; p<np; ++p){
//    int dest = (p+iam)%np;
//    if(dest!=iam){
//      myComms[p] = new Icomm(2*sizeof(Int));
//      Serialize(*myComms.back(),&iam,1);
//      Serialize(*myComms.back(),&dest,1);
//    }
//  }
//
//  gMaxIrecv = 0;
// 
//  //this is a all to all
//  for(int p = 0; p<np; ++p){
//    int dest = (p+iam)%np;
//    if(dest!=iam){
//      signal_data((char*)myComms[p]->front(),myComms[p]->size(),dest);
//    logfileptr->OFS()<<"Sending to "<<dest<<endl;
//    }
//  }
//
//    upcxx::wait();
//    logfileptr->OFS()<<"recv queue size "<<gIncomingRecv.size()<<endl;
//
//  //check
//  while(!gIncomingRecv.empty()){
//    upcxx::advance();
//    IncomingMessage * msg = gIncomingRecv.front();
//    msg->Wait(); 
//    logfileptr->OFS()<<"Received msg from "<<msg->Sender()<<endl;
//    delete msg;
//    gIncomingRecv.pop_front();
//  }
//
//  upcxx::barrier();
//
//  for(auto it = myComms.begin();it!=myComms.end(); ++it){
//    if(*it!=NULL){
//      delete *it;
//    }
//  }

  int numint = 2;
  if(argc>3){
    numint = atoi(argv[3]);
  }

  int nummsg = 1;
  if(argc>2){
    nummsg = atoi(argv[2]);
  }

  SYMPACK::vector<upcxx::global_ptr<int>> myComms(nummsg*np,upcxx::global_ptr<int>(NULL));
  for(int p = 0; p<np; ++p){
    int dest = (p+iam)%np;
    for(int msg = 0; msg<nummsg; ++msg){
      myComms[p*nummsg + msg] = upcxx::allocate<int>(iam,numint);
      int* ptr = (int*)myComms[p*nummsg];

      for(int curi = 0; curi<numint; curi++){
        ptr[curi]=dest;
      }
      ptr[numint-1] = iam;
    }
  }

  gMaxIrecv = atoi(argv[1]);
 
  upcxx::barrier();
  double tstart = get_time();

  //this is a all to all
  for(int msg = 0; msg<nummsg; ++msg){
    for(int p = 0; p<np; ++p){
      int dest = (p+iam)%np;
      if(dest!=iam){
        int * ptr = (int*)myComms[p*nummsg + msg];
        MsgMetadata meta;
        meta.src=iam;
        meta.tgt=dest;
        signal_data((char*)ptr,numint*sizeof(int),dest,meta);
        logfileptr->OFS()<<"Sending {"<<ptr[0]<<"..."<<ptr[numint-1]<<"} to "<<dest<<endl;
      }
    }
  }


  //check
  int numRecv = nummsg*(np-1);
  while(numRecv>0){
    //upcxx::advance(1,1);
    upcxx::advance();
    bool comm_found = false;
    if(!gIncomingRecvAsync.empty()){

      //find if there is some finished async comm
      auto it = TestAsyncIncomingMessage();
      if(it!=gIncomingRecvAsync.end()){
        comm_found = true;

        IncomingMessage * msg = *it;
        msg->Wait(); 
        logfileptr->OFS()<<"Received async msg from "<<msg->Sender()<<endl;
        int * ptr = (int*)msg->GetLocalPtr();
        logfileptr->OFS()<<"Message is {"<<ptr[0]<<"..."<<ptr[numint-1]<<"} "<<endl;

        remote_delete(msg->GetRemotePtr());

        delete msg;

        gIncomingRecvAsync.erase(it);
        numRecv--;
      }
    }

    if(!comm_found && !gIncomingRecv.empty()){
      auto it = gIncomingRecv.begin();
      IncomingMessage * msg = *it;
      msg->Wait(); 
      logfileptr->OFS()<<"Received msg from "<<msg->Sender()<<endl;
      int * ptr = (int*)msg->GetLocalPtr();
      logfileptr->OFS()<<"Message is {"<<ptr[0]<<"..."<<ptr[numint-1]<<"} "<<endl;

      remote_delete(msg->GetRemotePtr());

      delete msg;

      gIncomingRecv.erase(it);
      numRecv--;
    }
  }

  upcxx::async_wait();

  upcxx::barrier();

  double tstop = get_time();

  for(int p = 0; p<np; ++p){
    int dest = (p+iam)%np;
    for(int msg = 0; msg<nummsg; ++msg){
      if( dest == iam){
        upcxx::deallocate(myComms[p*nummsg+msg]);
      }
    }
  }


  delete logfileptr;

  if(iam==0){
  cout<<"Time: "<<tstop-tstart<<endl;
  }
  
  upcxx::finalize();
  return 0;
}


