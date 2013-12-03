/// @file run_cholesky.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2013-08-31
#define LAUNCH_ASYNC
#define DO_COMP

#define DBG_COL 0

//#define ADJUSTED_BUFFERS
//#define ASYNC_COPY

#include <upcxx.h>

#include <time.h>
#include <random>
#include <omp.h>

#include  "Environment.hpp"
//#include  "DistSparseMatrix.hpp"
#include  "NumVec.hpp"
#include  "NumMat.hpp"
#include  "utility.hpp"
#include  "ETree.hpp"
#include  "blas.hpp"
#include  "lapack.hpp"

#include  "LogFile.hpp"
#include  "upcxx_additions.hpp"

extern "C" {
#include <bebop/util/config.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/sparse_matrix_ops.h>

#include <bebop/util/get_options.h>
#include <bebop/util/init.h>
#include <bebop/util/malloc.h>
#include <bebop/util/timer.h>
#include <bebop/util/util.h>
}




#ifdef USE_TAU
#include "TAU.h"
#elif defined (PROFILE) || defined(PMPI)
#define TAU
#include "timer.hpp"
#endif

//#define MMAP(i,j,arr) arr(((i) % arr.m()),((j) % arr.n())) 
//#define MAP(i,j) MMAP(i,j,cmapp)
//#define r(i) ((i)%THREADS)
//#define c(i) 0
//#define MAP(i,j,bs) r(i/bs) + c(j/bs)*prow
//#define c(i) ((i)%THREADS)
//#define r(i) 0
//#define MAP(i,j,bs) r(i/bs)*pcol + c(j/bs)

using namespace LIBCHOLESKY;
using namespace std;

shared_lock fact_done;
namespace upcxx{
  template <typename T> inline void ldacopy(Int m, Int n, global_ptr<T> A, Int lda, global_ptr<T> B, Int ldb){
    for(Int j=0;j<n;j++){
      copy<T>(A + j*lda,B+j*ldb,m);
    }
  }


  template <typename T> inline void ldacopy_async(Int m, Int n, global_ptr<T> A, Int lda, global_ptr<T> B, Int ldb, event * e = NULL){
    for(Int j=0;j<n;j++){
      async_copy<T>(A + j*lda,B+j*ldb,m,&e[j]);
    }
  }
}


namespace LIBCHOLESKY{


class FBMatrix;

    LogFile * logfileptr;
void Update_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, upcxx::global_ptr<double> remoteFactorPtr);
void Aggregate_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, global_ptr<double> remoteAggregatePtr);
void Factor_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j);
void Gathercopy(upcxx::global_ptr<FBMatrix> Objptr, global_ptr<double> Adest, Int j);






class FBMatrix{
public:
  //TODO delete this
  DblNumMat W2;
  DblNumMat W;
  DblNumMat Achunk;
  upcxx::shared_array<upcxx::global_ptr<FBMatrix> > * RemoteObjPtr;

  //lock should be initialized with the number of contribution a block column is receiving
  IntNumVec AggLock;

  Int n;
  Int blksize;
  Int pcol, prow;
  Int np, iam;

  FBMatrix():blksize(1),n(0){
    np=THREADS;
    iam=MYTHREAD;
  }


  inline Int r(Int i){ return ((i)%THREADS);}
  inline Int c(Int i){return 0;}
  inline Int MAP(Int i, Int j) {return r(i/blksize) + c(j/blksize)*prow;}

  inline Int global_col_to_local(Int j){ return ((j)/(np*blksize))*blksize; }

  void Allocate(Int pn, Int pblksize){
    n=pn;
    blksize=pblksize;
      
    Int totBlk = n/blksize;
    Int remaining = n-totBlk*blksize;
    if(remaining>0) totBlk++;
  
    Int localBlk = totBlk/np;
    Int additionalBlock = totBlk - np*localBlk;
    if(iam<additionalBlock){ localBlk++; }
   
//    for(Int i=0;i<localBlk;i++){
//      Int prevBlk = (iam) +np*(i>0?i-1:0); 
//      logfileptr->OFS()<<"block "<<i<<" of size "<<min(n-(prevBlk)*blksize,blksize)<<endl;
//    } 
    

    Int prevBlk = (iam) +np*(localBlk-1); 
    Int chksize = (localBlk-1)*blksize + min(n-(prevBlk)*blksize,blksize);
//    logfileptr->OFS()<<"chksize="<<chksize<<endl;
    Achunk.Resize(n,chksize);
    SetValue(Achunk,(double)iam);

    //allocate global workspace
    W.Resize(n,n);
    SetValue(W,0.0);

    W2.Resize(n,n);
    SetValue(W2,0.0);
  }

  void Distribute( DblNumMat & Aorig){
    SetValue(Achunk,(double)iam);
    for(Int j = 0; j<n;j+=blksize){ 
      Int local_j = global_col_to_local(j);
      Int jb = min(blksize, n-j+1);
      if(iam==MAP(j,j)){
//        logfileptr->OFS()<<"Copying column "<<j<<endl;

        global_ptr<double> Aorigptr(&Aorig(0,j));
        global_ptr<double> Adestptr(&Achunk(0,local_j));
//        global_ptr<double> Aorigptr(&Aorig(j,j));
//        global_ptr<double> Adestptr(&Achunk(j,local_j));
//        upcxx::copy(Aorigptr,Adestptr,(n-j)*jb);
//        upcxx::ldacopy(n-j,jb,Aorigptr,n,Adestptr,n);
        upcxx::ldacopy(n,jb,Aorigptr,n,Adestptr,n);
      }
    }

   // SetValue(Achunk,(double)iam);


    Int numStep = ceil((double)n/(double)blksize);

    if(iam==0){
      cout<<"Factorization will require "<<numStep<<" steps"<<endl;
    }


    //Allocate the lock vector
    AggLock.Resize(numStep);
    for(int js=0; js<numStep;js++){
        int j = min(js*blksize,n-1);
        Int jb = min(blksize, n-j+1);

        int numProc = 0;
        if(j>=np*blksize){
          //col j must receive from jstep -1 proc of the previous steps
          numProc = np;
        }
        else{
          numProc = max(1,js -1);
        }
        AggLock[js]=1;//numProc;
    }
 
//    logfileptr->OFS()<<AggLock<<endl;

  }

  void Gather( DblNumMat & Adest){
    Adest.Resize(n,n);
    SetValue(Adest,0.0);

    for(Int j = 0; j<n;j+=blksize){ 
      Int local_j = global_col_to_local(j);
      Int jb = min(blksize, n-j+1);
      Int target = MAP(j,j);
//      logfileptr->OFS()<<"Copying column "<<j<<endl;
//      global_ptr<double> Adestptr(&Adest(0,j));
      global_ptr<double> Adestptr(&Adest(j,j));
      if(iam!=target){
        upcxx::async(target)(Gathercopy,(*RemoteObjPtr)[target],Adestptr,j);
      }
      else{

//        global_ptr<double> Asrcptr(&Achunk(0,local_j));
//        upcxx::copy(Asrcptr,Adestptr,(n)*jb);
        global_ptr<double> Asrcptr(&Achunk(j,local_j));
        upcxx::ldacopy(n-j,jb,Asrcptr,n,Adestptr,n);
      }
//      logfileptr->OFS()<<"Done Copying column "<<j<<endl;
    }
    upcxx::wait();
  }



void Aggregate_once(Int j, DblNumMat &DistW){
  //j is local
  Int local_j = global_col_to_local(j);

  Int jb = min(blksize, n-j+1);
  //aggregate previous updates
  TAU_FSTART(Aggregate);

    //proc holding update is
      for(int jj=0;jj<jb && n-j-jj>0;jj++){
        logfileptr->OFS()<<"Aggregating subcol "<<jj<< " to col "<<j<<" "<<DistW.m()<<" "<<DistW.n()<< endl;
#ifdef DO_COMP
#ifdef ADJUSTED_BUFFERS
if(j==DBG_COL*blksize){
        logfileptr->OFS()<<"FB Original second panel"<<endl;
//        for(Int i = j; i<n;i++){
        for(Int i = j; i<5;i++){
          for(Int jj = 0; jj<blksize;jj++){
            logfileptr->OFS()<<Achunk(i,local_j+jj)<<" ";
          }
          logfileptr->OFS()<<endl;
        }


        logfileptr->OFS()<<"FB Aggregating "<<n-j-jj<<" elements"<<endl;
        logfileptr->OFS()<<"FB DistW("<<0+jj<<","<<0+jj<<") = "<<DistW(0+jj,0+jj)<<endl;
        logfileptr->OFS()<<"FB Achunk("<<j+jj<<","<<local_j+jj<<") = "<<Achunk(j+jj,local_j+jj)<<endl;
}






        blas::Axpy((n-j-jj), -1.0, &DistW(0+jj,0+jj),1,&Achunk(j+jj,local_j+jj),1);
#else
if(j==DBG_COL*blksize){
        logfileptr->OFS()<<"FB Original second panel"<<endl;
//        for(Int i = 0; i<n;i++){
        for(Int i = 0; i<j+5;i++){
          for(Int jj = 0; jj<blksize;jj++){
            logfileptr->OFS()<<Achunk(i,local_j+jj)<<" ";
          }
          logfileptr->OFS()<<endl;
        }


        logfileptr->OFS()<<"FB Aggregating "<<n-j-jj<<" elements"<<endl;
        logfileptr->OFS()<<"FB DistW("<<j+jj<<","<<0+jj<<") = "<<DistW(j+jj,0+jj)<<endl;
        logfileptr->OFS()<<"FB Achunk("<<j+jj<<","<<local_j+jj<<") = "<<Achunk(j+jj,local_j+jj)<<endl;

        //col j is not supposed to be used so accumulate inside it
        blas::Axpy((n-j-jj), 1.0, &DistW(j+jj,0+jj),1,&W2(j+jj,j+jj),1);
}


        blas::Axpy((n-j-jj), -1.0, &DistW(j+jj,0+jj),1,&Achunk(j+jj,local_j+jj),1);
#endif
#endif
      }
  TAU_FSTOP(Aggregate);
}




void Factor(Int j){
  Int jb = min(blksize, n-j+1);
  //j is local
  Int local_j = global_col_to_local(j);

  TAU_FSTART(Factor);
#ifdef DO_COMP
  lapack::Potrf( 'L', jb, &Achunk(j,local_j ), n);
  if(n-j-jb>0){
    blas::Trsm('R','L','T','N',n-j-jb,jb, 1.0, &Achunk(j,local_j), n, &Achunk(j+jb,local_j), n);
  }
#endif
  TAU_FSTOP(Factor);
}

void Update(Int j, Int i, DblNumMat & Factor){
  //i is local, j is global
  Int jb = min(blksize, n-j+1);

  TAU_FSTART(Update);
  //call dgemm
  Int jbi = min(blksize, n-i+1);

#ifdef DO_COMP
#ifdef ADJUSTED_BUFFERS
  blas::Syrk('L','N', jbi ,jb,1.0,&Factor(i-j,0),n,1.0,&W(i,i),n); 

  if(n-i-jb>0){
    blas::Gemm('N','T',n-i-jb,jb,jb,1.0,&Factor(i-j+jb,0),n,&Factor(i-j,0),n,1.0,&W(i+jb,i),n);
  }
#else
  //blas::Syrk('L','N', jbi ,jb,1.0,&Factor(i,0),n,1.0,&W(i,i),n); 

  //if(n-i-jbi>0){
  //  blas::Gemm('N','T',n-i-jbi,jbi,jb,1.0,&Factor(i+jbi,0),n,&Factor(i,0),n,1.0,&W(i+jbi,i),n);
  //}
    blas::Gemm('N','T',n-i,jbi,jb,1.0,&Factor(i,0),n,&Factor(i,0),n,1.0,&W(i,i),n);
#endif
#endif

  TAU_FSTOP(Update);
}

};

  void Gathercopy(upcxx::global_ptr<FBMatrix> Objptr, global_ptr<double> Adest, Int j){
      assert(Objptr.tid()==MYTHREAD);
      FBMatrix & A = *Objptr;
//      logfileptr->OFS()<<"Performing the copy of column "<<j<<endl;
      Int local_j = A.global_col_to_local(j);
      Int jb = min(A.blksize, A.n-j+1);
      global_ptr<double> Asrc(&A.Achunk(j,local_j));
      upcxx::ldacopy(A.n-j,jb,Asrc,A.n,Adest,A.n);
  }



void Factor_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j){
  assert(Aptr.tid()==MYTHREAD);
  FBMatrix & A = *Aptr;


  upcxx::drain();


    logfileptr->OFS()<<"********************************************"<<endl;
    //Factor the column
    logfileptr->OFS()<<"Factoring column "<<j<<endl;
    A.Factor(j);
  
    logfileptr->OFS()<<"Factoring column done "<<j<<endl;
    //Launch the updates
    Int local_j = A.global_col_to_local(j);
    //for(Int i = j+A.blksize; i< min(A.n,j+(A.np)*A.blksize);i+=A.blksize){
    for(Int i = j+A.blksize; i< min(A.n,j+(A.np+1)*A.blksize);i+=A.blksize){
       Int target = A.MAP(i,j);
#ifdef ADJUSTED_BUFFERS
       global_ptr<double> AchkPtr(&A.Achunk(j,local_j));
#else
       global_ptr<double> AchkPtr(&A.Achunk(0,local_j));
#endif
       logfileptr->OFS()<<"Launching update from "<<j<<" on P"<<target<<endl;
#ifdef LAUNCH_ASYNC
        upcxx::async(target)(Update_Async,(*A.RemoteObjPtr)[target],j,AchkPtr);
#endif
    }

    if(j+A.blksize>=A.n){
     logfileptr->OFS()<<"Unlocking quit"<<endl;
     fact_done.unlock(); 
    }

    logfileptr->OFS()<<"Quitting Factor_Async"<<endl;
    logfileptr->OFS()<<"********************************************"<<endl;
}


void Aggregate_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, global_ptr<double> remoteAggregatePtr){
  assert(Aptr.tid()==MYTHREAD);
  FBMatrix & A = *Aptr;
  //fetch data
  logfileptr->OFS()<<"Aggregating Fetching data from P"<<remoteAggregatePtr.tid()<<endl;
  Int jb = min(A.blksize, A.n-j+1);

#ifdef ADJUSTED_BUFFERS
   DblNumMat RemoteAggregate(A.n-j,jb); 
#ifdef ASYNC_COPY
  TAU_FSTART(ldacopy_async);
   upcxx::event e[RemoteAggregate.n()]; 
   upcxx::ldacopy_async(RemoteAggregate.m(),RemoteAggregate.n(),remoteAggregatePtr,A.n-j,RemoteAggregate.GData(),A.n-j,e);
   for(int i=0;i<RemoteAggregate.n();i++){e[i].wait();}
  TAU_FSTOP(ldacopy_async);
#else
  TAU_FSTART(ldacopy);
  upcxx::ldacopy(RemoteAggregate.m(),RemoteAggregate.n(),remoteAggregatePtr,A.n-j,RemoteAggregate.GData(),A.n-j);
  logfileptr->OFS()<<"Aggregating update to "<<j<<endl;
  TAU_FSTOP(ldacopy);
#endif
#else
   DblNumMat RemoteAggregate(A.n,jb); 
  TAU_FSTART(copy);
    upcxx::copy(remoteAggregatePtr,RemoteAggregate.GData(),RemoteAggregate.m()*RemoteAggregate.n());
  TAU_FSTOP(copy);

#endif

if(j==DBG_COL*A.blksize){
//  logfileptr->OFS()<<"FB RemoteAggregate:"<<RemoteAggregate<<endl;
  logfileptr->OFS()<<"FB RemoteAggregate:"<<endl;
        for(Int i = 0; i<j+5;i++){
          for(Int jj = 0; jj<A.blksize;jj++){
            logfileptr->OFS()<<RemoteAggregate(i,jj)<<" ";
          }
          logfileptr->OFS()<<endl;
        }

}


   A.Aggregate_once(j, RemoteAggregate);
  logfileptr->OFS()<<"Done Aggregating update to "<<j<<endl;

if(j==DBG_COL*A.blksize){
        logfileptr->OFS()<<"FB accumulated aggregate"<<endl;
//        for(Int i = 0; i<A.n;i++){
        for(Int i = 0; i<j+5;i++){
          for(Int jj = 0; jj<A.blksize;jj++){
            logfileptr->OFS()<<A.W2(i,j+jj)<<" ";
          }
          logfileptr->OFS()<<endl;
        }
}






   A.AggLock[j/A.blksize]--;
   //logfileptr->OFS()<<A.AggLock<<" "<<j/A.blksize<<endl;
  //do the sum
   logfileptr->OFS()<<A.AggLock<<endl;
  if(A.AggLock[j/A.blksize]==0){
    //launch the factorization of column j
   logfileptr->OFS()<<"Col "<<j<<" updated. Can be factored "<<endl;
    


#ifdef ADJUSTED_BUFFERS
if(j==DBG_COL*A.blksize){
      Int local_j = A.global_col_to_local(j);
        logfileptr->OFS()<<endl;
        logfileptr->OFS()<<endl;
        logfileptr->OFS()<<"FB Second panel after aggregating the update"<<endl;
        for(Int i = j; i<j+5;i++){
//        for(Int i = j; i<A.n;i++){
          for(Int jj = 0; jj<A.blksize;jj++){
            logfileptr->OFS()<<A.Achunk(i,local_j+jj)<<" ";
          }
          logfileptr->OFS()<<endl;
        }

}
#else
if(j==DBG_COL*A.blksize){
      Int local_j = A.global_col_to_local(j);
        logfileptr->OFS()<<endl;
        logfileptr->OFS()<<endl;
        logfileptr->OFS()<<"FB Second panel after aggregating the update"<<endl;
//        for(Int i = 0; i<A.n;i++){
        for(Int i = 0; i<j+5;i++){
          for(Int jj = 0; jj<A.blksize;jj++){
            logfileptr->OFS()<<A.Achunk(i,local_j+jj)<<" ";
//            logfileptr->OFS()<<A.W(i,j+jj)<<" ";
          }
          logfileptr->OFS()<<endl;
        }

}
#endif




#ifdef LAUNCH_ASYNC
   upcxx::async(A.iam)(Factor_Async,Aptr,j);
#endif
  } 
    logfileptr->OFS()<<"Quitting Aggregate_Async"<<endl<<endl;
}


void Update_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, upcxx::global_ptr<double> remoteFactorPtr){
  assert(Aptr.tid()==MYTHREAD);
  FBMatrix & A = *Aptr;


  logfileptr->OFS()<<"Fetching factor from P"<<remoteFactorPtr.tid()<<endl;

#ifdef ADJUSTED_BUFFERS
  DblNumMat RemoteFactor(A.n-j,A.blksize);
#ifdef ASYNC_COPY
  TAU_FSTART(ldacopy_async);
   upcxx::event e[RemoteFactor.n()]; 
   upcxx::ldacopy_async(RemoteFactor.m(),RemoteFactor.n(),remoteFactorPtr,A.n,RemoteFactor.GData(),A.n-j,e);
   for(int i=0;i<RemoteFactor.n();i++){e[i].wait();}
  TAU_FSTOP(ldacopy_async);
#else
  TAU_FSTART(ldacopy);
  upcxx::ldacopy(RemoteFactor.m(),RemoteFactor.n(),remoteFactorPtr,A.n,RemoteFactor.GData(),A.n-j);
  TAU_FSTOP(ldacopy);
#endif
#else
  DblNumMat RemoteFactor(A.n,A.blksize);
  TAU_FSTART(copy);
  upcxx::copy(remoteFactorPtr,RemoteFactor.GData(),RemoteFactor.m()*RemoteFactor.n());
  TAU_FSTOP(copy);
#endif

  //do all my updates with i>j
  for(Int i = j+A.blksize; i<A.Achunk.m();i+=A.blksize){
    if(A.iam==A.MAP(i,j)){
//      logfileptr->OFS()<<"Updating column "<<i<<" with columns from P"<<remoteFactorPtr.tid()<<endl;
      A.Update(j, i, RemoteFactor);

    
    if(i==DBG_COL*A.blksize){
       logfileptr->OFS()<<"FB RemoteFactor is "<<endl;//<<RemoteFactor<<endl;
        for(Int ii = 0; ii<i+5;ii++){
          for(Int jj = 0; jj<A.blksize;jj++){
            logfileptr->OFS()<<RemoteFactor(ii,jj)<<" ";
          }
          logfileptr->OFS()<<endl;
        }
        
        logfileptr->OFS()<<endl;
        logfileptr->OFS()<<"FB Second panel aggregate"<<endl;
#ifdef ADJUSTED_BUFFERS
        for(Int ii = i; ii<A.n;ii++){
#else
        for(Int ii = 0; ii<i+5;ii++){
#endif
          for(Int jj = 0; jj<A.blksize;jj++){
            logfileptr->OFS()<<A.W(ii,i+jj)<<" ";
          }
          logfileptr->OFS()<<endl;
        }
    }


    }
  }

  IntNumVec myLastUpd(A.np);
  SetValue(myLastUpd,-1);
  
  for(Int i = j+A.blksize; i< min(A.n,j+(A.np+1)*A.blksize);i+=A.blksize){
    for(Int j2 = j; j2< i;j2+=A.blksize){
      if(A.iam==A.MAP(i,j2)){
        Int target = A.MAP(i,i);
        logfileptr->OFS()<<"j2="<<j2<<" vs "<<myLastUpd[target]<<" on P"<<target<<endl;
        myLastUpd[target] = max(j2,myLastUpd[target]);
      }
    }
  }

  logfileptr->OFS()<<myLastUpd<<endl;


  for(Int i = j+A.blksize; i< min(A.n,j+(A.np+1)*A.blksize);i+=A.blksize){
      if(A.iam==A.MAP(i,j)){
        Int target = A.MAP(i,i);
//  for(Int target = 0; target< A.np;target++){
        if(myLastUpd[target]==j){
            logfileptr->OFS()<<"Updating from "<<j<<", launching aggregate on P"<<target<<" for column "<<i<<endl;
myLastUpd[target]=-1;
#ifdef LAUNCH_ASYNC


#ifdef ADJUSTED_BUFFERS
        global_ptr<double> WbufPtr(&A.W(i,i));
#else
        global_ptr<double> WbufPtr(&A.W(0,i));
#endif


        upcxx::async(target)(Aggregate_Async,(*A.RemoteObjPtr)[target],i,WbufPtr);
#endif
      }
    }
  }
  
//  for(Int i = j+A.blksize; i< min(A.n,j+(A.np+1)*A.blksize);i+=A.blksize){
//  //for(Int i = j+A.blksize; i<A.Achunk.m();i+=A.blksize){
//    if(A.iam==A.MAP(i,j)){
//      int jstep = j/(A.np*A.blksize); 
//      int istep = i/(A.np*A.blksize);
//      logfileptr->OFS()<<" jstep="<<jstep<<" istep="<<istep<<endl;
//      if(jstep+1==istep || jstep==istep)
//      {
//        //If it's my last update for P(i,i)
//        //Launch Aggregate_async on processor P(i,i)
//        Int target = A.MAP(i,i);
//        logfileptr->ofs()<<"launching aggregate on P"<<target<<endl;
//#ifdef launch_async
//
//
//#ifdef adjusted_buffers
//        global_ptr<double> wbufptr(&a.w(i,i));
//#else
//        global_ptr<double> wbufptr(&a.w(0,i));
//#endif
//
//
//        upcxx::async(target)(aggregate_async,(*a.remoteobjptr)[target],i,wbufptr);
//#endif
//      } 
//    }
//  }
    logfileptr->OFS()<<"Quitting Update_Async"<<endl<<endl;
}








}










int main(int argc, char **argv) 
{
  upcxx::init(&argc,&argv);

#if defined(PROFILE) || defined(PMPI)
  TAU_PROFILE_INIT(argc, argv);
#endif


  int np = THREADS;
  int iam = MYTHREAD;

  try{


    logfileptr = new LogFile(iam);
    LogFile & logfile = *logfileptr;
  
    logfile<<"********* LOGFILE OF P"<<iam<<" *********"<<endl;
    logfile<<"**********************************"<<endl;







    // *********************************************************************
    // Input parameter
    // *********************************************************************
    std::map<std::string,std::string> options;

    OptionsCreate(argc, argv, options);

    std::string filename;
    if( options.find("-in") != options.end() ){
      filename= options["-in"];
    }
    std::string informatstr;
    if( options.find("-inf") != options.end() ){
      informatstr= options["-inf"];
    }
    Int blksize = 1;
    if( options.find("-b") != options.end() ){
      blksize= atoi(options["-b"].c_str());
    }


    upcxx::shared_array<upcxx::global_ptr<FBMatrix> > Aobjects;
    Aobjects.init(np);

    upcxx::global_ptr<FBMatrix> Afactptr = upcxx::Create<FBMatrix>();
    Aobjects[iam] = Afactptr;
    upcxx::barrier();

    
    //upcxx::global_ptr<FBMatrix> Afactptr = upcxx::Create<FBMatrix>();

    FBMatrix & Afact = *Afactptr;
    Afact.RemoteObjPtr = &Aobjects;

    //Read the input matrix
    sparse_matrix_file_format_t informat;
    informat = sparse_matrix_file_format_string_to_enum (informatstr.c_str());
    sparse_matrix_t* Atmp = load_sparse_matrix (informat, filename.c_str());
    sparse_matrix_expand_symmetric_storage (Atmp);
    sparse_matrix_convert (Atmp, CSR);
    Afact.n = ((csr_matrix_t *)Atmp->repr)->n;
    
    if(MYTHREAD==0){
      cout<<"Matrix order is "<<Afact.n<<endl;
    }

    DblNumMat A(Afact.n,Afact.n);
    csr_matrix_expand_to_dense (A.Data(), 0, Afact.n, Atmp->repr);
    destroy_sparse_matrix (Atmp);


    DblNumMat RefA = A;
    DblNumMat D = A;
    DblNumMat SecondD(A.m(),blksize);
    DblNumMat FSecondD(A.m(),blksize);
    DblNumMat FFSecondD(A.m(),blksize);
    DblNumMat AggD(A.m(),blksize);
    SetValue(AggD,0.0);


    Real timeSta, timeEnd;

    if(MYTHREAD==0)
    {
      timeSta =  omp_get_wtime( );
      TAU_FSTART(POTRF);
      lapack::Potrf( 'L', D.n(), D.Data(), D.n());
      TAU_FSTOP(POTRF);
      timeEnd =  omp_get_wtime( );
      cout<<"REF POTRF: "<<timeEnd-timeSta<<endl;


      Int i =DBG_COL*blksize;
      Int jbi = min(blksize, A.n()-i+1);
      std::copy(&A(0,i),&A(0,i+jbi),&SecondD(0,0)); 
      std::copy(&A(0,i),&A(0,i+jbi),&FSecondD(0,0)); 



      //compute AggD
      blas::Syrk('L','N', jbi ,i,1.0,&D(i,0),A.n(),1.0,&AggD(i,0),A.n()); 

//logfileptr->OFS()<<"AFTER SYRK"<<AggD<<endl;


      if(A.n()-i-jbi>0){
        blas::Gemm('N','T',A.n()-i-jbi,jbi,i,1.0,&D(i+jbi,0),A.n(),&D(i,0),A.n(),1.0,&AggD(i+jbi,0),A.n());
      }




      //Accumulate AggD in FSecondD
      for(int jj=0;jj<jbi && A.n()-i-jj>0;jj++){
        blas::Axpy((A.n()-i-jj), -1.0, &AggD(i+jj,0+jj),1,&FSecondD(i+jj,0+jj),1);
      }

      FFSecondD = FSecondD;
//      std::copy(&FSecondD(0,0),FSecondD.Data()+FSecondD.m()*blksize,&FFSecondD(0,0)); 
      lapack::Potrf( 'L', D.n()-blksize, &FFSecondD(i,0), D.n());
    }

    upcxx::barrier();

    //TODO This needs to be fixed
    Afact.prow = sqrt(np);
    Afact.pcol = Afact.prow;
    np = Afact.prow*Afact.pcol;
    if(MYTHREAD==0){
      cout<<"Number of cores to be used: "<<np<<endl;
    }



    //Allocate chunks of the matrix on each processor
    Afact.Allocate(A.n(),blksize);
    Afact.Distribute(A);


      TAU_FSTART(FANBOTH);
    if(MYTHREAD==Afact.MAP(Afact.n-1,Afact.n-1)){
      logfileptr->OFS()<<"LOCK IS TAKEN BY"<<MYTHREAD<<endl;
      fact_done.lock();
    }

   logfileptr->OFS()<<"initializing"<<endl;

    upcxx::barrier();

      if(MYTHREAD==Afact.MAP(0,0)){
        upcxx::async(Afact.iam)(Factor_Async,Afactptr,0);
      }


    while(fact_done.islocked()){upcxx::progress();};
//    while(fact_done.islocked()){upcxx::drain();};
//    while(fact_done.islocked()){upcxx::peek();};
    fact_done.lock();
      TAU_FSTOP(FANBOTH);
    fact_done.unlock();

    upcxx::barrier();
   logfileptr->OFS()<<"quitting"<<endl;

Int ncol = 5;
Int scol = DBG_COL*blksize;


    //gather all data on P0
    DblNumMat Afinished;
    Afact.Gather(Afinished);
    upcxx::barrier();
    for(Int i = scol+0; i<scol+ncol;i++){
      for(Int j = scol; j<scol+ncol;j++){
//    for(Int i = 0; i<Afinished.m();i++){
//      for(Int j = 0; j<Afinished.n();j++){
        logfileptr->OFS()<<Afinished(i,j)<<" ";
      }
      logfileptr->OFS()<<endl;
    }
      if(MYTHREAD==0)
      {
        logfileptr->OFS()<<endl;
        logfileptr->OFS()<<endl;
        for(Int i = scol; i<scol+ncol;i++){
          for(Int j = scol; j<scol+ncol;j++){
            logfileptr->OFS()<<D(i,j)<<" ";
          }
          logfileptr->OFS()<<endl;
        }



        logfileptr->OFS()<<endl;
        logfileptr->OFS()<<"Original second panel"<<endl;
        //for(Int i = 0; i<Afact.n;i++){
        for(Int i = 0; i<ncol;i++){
          for(Int j = 0; j<blksize;j++){
            logfileptr->OFS()<<SecondD(i,j)<<" ";
          }
          logfileptr->OFS()<<endl;
        }


        logfileptr->OFS()<<endl;
        logfileptr->OFS()<<"Second panel aggregate"<<endl;
        //for(Int i = 0; i<Afact.n;i++){
        for(Int i = 0; i<ncol;i++){
          for(Int j = 0; j<blksize;j++){
            logfileptr->OFS()<<AggD(i,j)<<" ";
          }
          logfileptr->OFS()<<endl;
        }


        logfileptr->OFS()<<endl;
        logfileptr->OFS()<<endl;
        logfileptr->OFS()<<"Second panel after aggregating the update"<<endl;
        //for(Int i = 0; i<Afact.n;i++){
        for(Int i = 0; i<ncol;i++){
          for(Int j = 0; j<blksize;j++){
            logfileptr->OFS()<<FSecondD(i,j)<<" ";
          }
          logfileptr->OFS()<<endl;
        }

        logfileptr->OFS()<<endl;
        logfileptr->OFS()<<"Second panel factored"<<endl;
        //for(Int i = 0; i<Afact.n;i++){
        for(Int i = 0; i<ncol;i++){
          for(Int j = 0; j<blksize;j++){
            logfileptr->OFS()<<FFSecondD(i,j)<<" ";
          }
          logfileptr->OFS()<<endl;
        }



      }

    upcxx::barrier();

    upcxx::barrier();
    upcxx::barrier();

    upcxx::barrier();
    //destroy
//    upcxx::Destroy<FBMatrix>(Afactptr);



      if(MYTHREAD==0)
      {
    //check the result
    double norm = 0;

for(Int j=0;j<Afinished.n();j++){
  for(Int i=0;i<j;i++){
    D(i,j)=0.0;
  }
}
    //do a solve
    Int n = Afinished.n();
    Int nrhs = n;
    DblNumMat RHS(n,nrhs);
    DblNumMat XTrue(n,nrhs);
    UniformRandom(XTrue);
    blas::Gemm('N','N',n,nrhs,n,1.0,&RefA(0,0),n,&XTrue(0,0),n,0.0,&RHS(0,0),n);

    DblNumMat X = RHS;
    lapack::Potrs('L',n,nrhs,&D(0,0),n,&X(0,0),n);
    blas::Axpy(n*nrhs,-1.0,&XTrue(0,0),1,&X(0,0),1);
    norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
    printf("Norm of residual after SOLVE for REF POTF2 is %10.2e\n",norm);

    X = RHS;
    lapack::Potrs('L',n,nrhs,&Afinished(0,0),n,&X(0,0),n);
    blas::Axpy(n*nrhs,-1.0,&XTrue(0,0),1,&X(0,0),1);
    norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
    printf("Norm of residual after SOLVE for FAN-BOTH is %10.2e\n",norm);
      

Int maxI = 0; 
Int maxJ = 0;
double maxDiff = 0.0;
for(Int j=0;j<Afinished.n();j++){
  for(Int i=j;i<Afinished.m();i++){
    double diff = abs(Afinished(i,j)-D(i,j));
    if(diff>=maxDiff){
      maxDiff=diff;
      maxI = i;
      maxJ = j;
    }
  }
}
    printf("Norm of max diff FAN-BOTH is %10.2e (%d,%d)\n",maxDiff,maxI,maxJ);

    blas::Axpy(n*n,-1.0,&D(0,0),1,&Afinished(0,0),1);
    norm = lapack::Lange('F',n,n,&Afinished(0,0),n);
    printf("Norm of residual POTRF - FAN-BOTH is %10.2e\n",norm);
}
////
////
////
////
////    X = RHS;
////    lapack::Posv('L',n,nrhs,&RefA(0,0),n,&X(0,0),n);
////    blas::Axpy(n*nrhs,-1.0,&XTrue(0,0),1,&X(0,0),1);
////    norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
////    printf("Norm of residual after POSV is %10.2e\n",norm);
    upcxx::finalize();
delete logfileptr;
return;


////    upcxx::wait();
  }
  catch( std::exception& e )
  {
    std::cerr << "Exception with message: "
      << e.what() << std::endl;
#ifndef _RELEASE_
    DumpCallStack();
#endif
  }


  return 0;
}





#if 0
//error check
if(0){
  for(int jj=0;jj<jb;jj++){
    DblNumVec res(n-j);
    memcpy(res.Data(),&A(j,j+jj),(n-j)*sizeof(double));
    blas::Axpy((n-j), -1.0, &D(j,j+jj),1,&res(0),1);
    double norm = lapack::Lange('F',n-j,1,&res(0),n)/lapack::Lange('F',n-j,1,&D(j,j+jj),n);
    statusOFS<<"Norm of res of col "<<j+jj<<": "<<norm<<std::endl;
    norm = lapack::Lange('M',n-j,1,&res(0),n)/lapack::Lange('M',n-j,1,&D(j,j+jj),n);
    statusOFS<<"Max Norm of res of col "<<j+jj<<": "<<norm<<std::endl;

    //DblNumVec res(n);
    //memcpy(res.Data(),&A(0,j+jj),(n)*sizeof(double));
    //blas::Axpy((n), -1.0, &D(0,j+jj),1,&res(0),1);

    //double norm = lapack::Lange('F',n,1,&res(0),n)/lapack::Lange('F',n,1,&D(j,j+jj),n);
    //statusOFS<<"Norm of res of col "<<j+jj<<": "<<norm<<std::endl;
    //norm = lapack::Lange('M',n,1,&res(0),n)/lapack::Lange('M',n,1,&D(j,j+jj),n);
    //statusOFS<<"Max Norm of res of col "<<j+jj<<": "<<norm<<std::endl;
  }
}

if(0)
{
  Int nrhs = j+jb;
  DblNumMat RHS(j+jb,nrhs);
  DblNumMat XTrue(j+jb,nrhs);
  UniformRandom(XTrue);
  blas::Gemm('N','N',j+jb,nrhs,j+jb,1.0,&RefA(0,0),n,&XTrue(0,0),j+jb,0.0,&RHS(0,0),j+jb);

  DblNumMat X = RHS;
  lapack::Potrs('L',j+jb,nrhs,&A(0,0),n,&X(0,0),j+jb);
  blas::Axpy((j+jb)*nrhs,-1.0,&XTrue(0,0),1,&X(0,0),1);
  double norm = lapack::Lange('F',j+jb,nrhs,&X(0,0),j+jb);
  statusOFS<< "Norm of solve at step "<<j<<" of size "<<j+jb<<"-by-"<<nrhs<<" is "<<norm<<endl;


}
#endif

