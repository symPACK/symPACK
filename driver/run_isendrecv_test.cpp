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


#define MYSCALAR double


using namespace SYMPACK;


int main(int argc, char **argv) 
{
  MPI_Init(&argc,&argv);

  MPI_Comm_rank(MPI_COMM_WORLD,&iam);
  MPI_Comm_size(MPI_COMM_WORLD,&np);

  logfileptr = new LogFile(iam);
  logfileptr->OFS()<<"********* LOGFILE OF P"<<iam<<" *********"<<endl;
  logfileptr->OFS()<<"**********************************"<<endl;


  int numint = 2;
  if(argc>3){
    numint = atoi(argv[3]);
  }

  int nummsg = 1;
  if(argc>2){
    nummsg = atoi(argv[2]);
  }

  std::vector<int *> myComms(nummsg*np,NULL);
  for(int p = 0; p<np; ++p){
    int dest = (p+iam)%np;
    for(int msg = 0; msg<nummsg; ++msg){
      myComms[p*nummsg + msg] = new int[numint];
      int* ptr = (int*)myComms[p];
      ptr[0]=iam;
      ptr[1]=dest;
    }
  }

  int gMaxIrecv = atoi(argv[1]);
  std::vector<MPI_Request> sendReqs;
  sendReqs.reserve(nummsg*np); 
  std::vector<int *> recvComms(gMaxIrecv,NULL);
  std::vector<MPI_Request> recvReqs;
  recvReqs.reserve(gMaxIrecv);


  double tstart = get_time();
  MPI_Barrier(MPI_COMM_WORLD);

  //this is a all to all
  for(int msg = 0; msg<nummsg; ++msg){
    for(int p = 0; p<np; ++p){
      int dest = (p+iam)%np;
      sendReqs.push_back(MPI_REQUEST_NULL);
      if(dest!=iam){
        MPI_Isend(myComms[p*nummsg + msg],numint*sizeof(int),MPI_BYTE,dest,dest,MPI_COMM_WORLD,&sendReqs.back());
        logfileptr->OFS()<<"Sending to "<<dest<<endl;
      }
    }
  }

  int recv_idx = 0;
  //launch the recv
  for(int msg = 0; msg<nummsg; ++msg){
    for(int p = 0; p<np; ++p){
      if(p!=iam){
        if(recvReqs.size()<gMaxIrecv){
          recvComms[recv_idx] = new int[numint];
          recvReqs.push_back(MPI_REQUEST_NULL);
          MPI_Irecv(recvComms[recv_idx],numint*sizeof(int),MPI_BYTE,MPI_ANY_SOURCE,iam,MPI_COMM_WORLD,&recvReqs.back());
          logfileptr->OFS()<<"Irecv from any"<<endl;
          recv_idx++;
        }
          
        if(recvReqs.size()>=gMaxIrecv){
          break;
        }
      }
    }

    if(recvReqs.size()>=gMaxIrecv){
      break;
    }
  }


  //check
  int numRecv = nummsg*(np-1);
  while(numRecv>0){
    bool comm_found = false;
      MPI_Status stat;

    if(gMaxIrecv>0){
      int recvIdx = 0;
      MPI_Waitany(gMaxIrecv,&recvReqs[0],&recvIdx,&stat);
      if(recvIdx!=MPI_UNDEFINED){
        comm_found = true;
        int sender = stat.MPI_SOURCE;
          
          logfileptr->OFS()<<"Received async msg from "<<sender<<endl;
          int * ptr = recvComms[recvIdx];
          logfileptr->OFS()<<" message is "<<ptr[0]<<" "<<ptr[1]<<endl;
          numRecv--;

      }
    }

    if(!comm_found){
      int * ptr = new int[numint];

      MPI_Recv(ptr,numint*sizeof(int),MPI_BYTE,MPI_ANY_SOURCE,iam,MPI_COMM_WORLD,&stat);
      logfileptr->OFS()<<"Received msg from "<<stat.MPI_SOURCE<<endl;
      logfileptr->OFS()<<" message is "<<ptr[0]<<" "<<ptr[1]<<endl;
      delete [] ptr;

      numRecv--;
    }
  }

  //wait

  MPI_Barrier(MPI_COMM_WORLD);
  
  double tstop = get_time();

  for(int p = 0; p<np; ++p){
    int dest = (p+iam)%np;
    for(int msg = 0; msg<nummsg; ++msg){
      delete [] myComms[p*nummsg + msg];
    }
  }

    for(int msg = 0; msg<recvComms.size(); ++msg){
      delete [] recvComms[msg];
    }

  delete logfileptr;
  if(iam==0){
    cout<<"Time: "<<tstop-tstart<<endl;
  }

  MPI_Finalize();
  return 0;
}


