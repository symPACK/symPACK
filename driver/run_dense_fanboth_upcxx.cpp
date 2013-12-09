/// @file run_cholesky.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2013-08-31
//#define _DEBUG_
#define LAUNCH_ASYNC

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



#ifdef USE_TAU 
#define TIMER_START(a) TAU_START(TOSTRING(a));
#define TIMER_STOP(a) TAU_STOP(TOSTRING(a));
#elif defined (PROFILE)
#define TIMER_START(a) TAU_FSTART(a);
#define TIMER_STOP(a) TAU_FSTOP(a);
#else
#define TIMER_START(a)
#define TIMER_STOP(a)
#endif

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

void Update_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, upcxx::global_ptr<double> remoteFactorPtr);
void Aggregate_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, global_ptr<double> remoteAggregatePtr);
void Factor_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j);
void Gathercopy(upcxx::global_ptr<FBMatrix> Objptr, global_ptr<double> Adest, Int j);






class FBMatrix{
public:
//  DblNumMat W;
//  DblNumMat Achunk;
  vector<DblNumMat> AchunkLower;
  vector<DblNumMat> WLower;
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
  ~FBMatrix(){
  }


  inline Int r(Int i){ return ((i)%THREADS);}
  inline Int c(Int i){return 0;}
  inline Int MAP(Int i, Int j) {return r(i/blksize) + c(j/blksize)*prow;}

//  inline Int c(Int i){ return ((i)%THREADS);}
//  inline Int r(Int i){return 0;}
//  inline Int MAP(Int i, Int j) {return r(i/blksize)*pcol + c(j/blksize);}

//  inline Int r(Int i){ return ((i)%prow);}
//  inline Int c(Int i){ return ((i)%pcol);}
//  inline Int MAP(Int i, Int j) {
//logfileptr->OFS()<<"MAP("<<i<<","<<j<<") = "<< r(i/blksize) <<" + "<<c(j/blksize)<<"*"<<prow<<" = "<<r(i/blksize) + c(j/blksize)*prow<<endl;
//return r(i/blksize) + c(j/blksize)*prow;
//
//}



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
   

    Int prevBlk = (iam) +np*(localBlk-1); 
    Int chksize = (localBlk-1)*blksize + min(n-(prevBlk)*blksize,blksize);

    Int lowerChkSize = (localBlk-1)*(n-iam*blksize)*blksize - blksize*blksize*np*(localBlk-1)*localBlk/2 + (n-(prevBlk)*blksize)*min(n-(prevBlk)*blksize,blksize); 
    logfileptr->OFS()<<"localBlk="<<localBlk<<endl;
    logfileptr->OFS()<<"prevBlk="<<prevBlk<<endl;
    logfileptr->OFS()<<"chksize*n="<<chksize*n<<endl;
    logfileptr->OFS()<<"lowerChkSize="<<lowerChkSize<<endl;

    Int numStep = ceil((double)n/(double)blksize);
    WLower.resize(numStep, DblNumMat(0,0) );
    AchunkLower.resize(localBlk, DblNumMat(0,0) );    
    for(Int j = 0; j<n;j+=blksize){ 
      Int local_j = global_col_to_local(j);
      Int jb = min(blksize, n-j);
      if(iam==MAP(j,j)){
        AchunkLower[local_j/blksize].Resize(n-j,jb);
      }
      WLower[j/blksize].Resize(n-j,jb);
      SetValue(WLower[j/blksize],0.0);
    }



//    Achunk.Resize(n,chksize);
//    SetValue(Achunk,(double)iam);

    //allocate global workspace
//    W.Resize(n,n);
//    SetValue(W,0.0);

  }

  void Distribute( DblNumMat & Aorig){
//    SetValue(Achunk,(double)iam);
    for(Int j = 0; j<n;j+=blksize){ 
      Int local_j = global_col_to_local(j);
      Int jb = min(blksize, n-j);
      if(iam==MAP(j,j)){
//        global_ptr<double> Aorigptr(&Aorig(0,j));
//        global_ptr<double> Adestptr(&Achunk(0,local_j));
//        upcxx::copy(Aorigptr,Adestptr,n*jb);

//logfileptr->OFS()<<"j="<<j<<" jb="<<jb<<" size is "<<Aorig.n()<<endl;
{

        global_ptr<double> Aorigptr(&Aorig(j,j));
        global_ptr<double> Adestptr(&AchunkLower[local_j/blksize](0,0));
        upcxx::ldacopy(n-j,jb,Aorigptr,n,Adestptr,n-j);
}

      }
    }

    Int numStep = ceil((double)n/(double)blksize);

    if(iam==0){
      cout<<"Factorization will require "<<numStep<<" steps"<<endl;
    }


    //Allocate the lock vector
    AggLock.Resize(numStep);
    for(int js=0; js<numStep;js++){
        int j = min(js*blksize,n-1);
        Int jb = min(blksize, n-j);

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
      Int jb = min(blksize, n-j);
      Int target = MAP(j,j);
//      logfileptr->OFS()<<"Copying column "<<j<<endl;
//      global_ptr<double> Adestptr(&Adest(0,j));
      global_ptr<double> Adestptr(&Adest(j,j));
//      if(iam!=target){
        upcxx::async(target)(Gathercopy,(*RemoteObjPtr)[target],Adestptr,j);
//      }
//      else{

//        global_ptr<double> Asrcptr(&Achunk(0,local_j));
//        upcxx::copy(Asrcptr,Adestptr,(n)*jb);
//        global_ptr<double> Asrcptr(&Achunk(j,local_j));
//        upcxx::ldacopy(n-j,jb,Asrcptr,n,Adestptr,n);
//      }
//      logfileptr->OFS()<<"Done Copying column "<<j<<endl;
    }
    upcxx::wait();
  }



void Aggregate_once(Int j, DblNumMat &DistW){
  //j is local
  Int local_j = global_col_to_local(j);
        DblNumMat & LocalChunk = AchunkLower[local_j/blksize];

  Int jb = min(blksize, n-j);
  //aggregate previous updates
  TIMER_START(Aggregate);

    //proc holding update is
      for(int jj=0;jj<jb && n-j-jj>0;jj++){
#ifdef _DEBUG_
        logfileptr->OFS()<<"Aggregating subcol "<<jj<< " to col "<<j<<" "<<DistW.m()<<" "<<DistW.n()<< endl;
#endif

//        blas::Axpy((n-j-jj), -1.0, &DistW(j+jj,0+jj),1,&Achunk(j+jj,local_j+jj),1);

        blas::Axpy((n-j-jj), -1.0, &DistW(0+jj,0+jj),1,&LocalChunk(jj,jj),1);
      }

  TIMER_STOP(Aggregate);
}




void Factor(Int j){
  Int jb = min(blksize, n-j);
  //j is local
  Int local_j = global_col_to_local(j);

  TIMER_START(Factor);
//  lapack::Potrf( 'L', jb, &Achunk(j,local_j ), n);
//  if(n-j-jb>0){
//    blas::Trsm('R','L','T','N',n-j-jb,jb, 1.0, &Achunk(j,local_j), n, &Achunk(j+jb,local_j), n);
//  }



  DblNumMat & LocalChunk = AchunkLower[local_j/blksize];
  lapack::Potrf( 'L', jb, &LocalChunk(0,0 ), n-j);
  if(n-j-jb>0){
    blas::Trsm('R','L','T','N',n-j-jb,jb, 1.0, &LocalChunk(0,0), n-j, &LocalChunk(jb,0), n-j);
  }


  TIMER_STOP(Factor);
}

void Update(Int j, Int i, DblNumMat & Factor){
  //i is local, j is global
  Int jb = min(blksize, n-j);

  TIMER_START(Update);
  //call dgemm
  Int jbi = min(blksize, n-i);

    
//    blas::Gemm('N','T',n-i,jbi,jb,1.0,&Factor(i,0),n,&Factor(i,0),n,1.0,&W(i,i),n);

    DblNumMat & WChunk = WLower[i/blksize];
    blas::Gemm('N','T',n-i,jbi,jb,1.0,&Factor(i-j,0),n-j,&Factor(i-j,0),n-j,1.0,&WChunk(0,0),n-i);

//for(int ii = 0; ii< WChunk.m();ii++){
//  for(int jj = 0; jj< jbi;jj++){
//    logfileptr->OFS()<<W(i+ii,i+jj)<< " ";
//  }
//  logfileptr->OFS()<<"  |  ";
//
//  for(int jj = 0; jj< jbi;jj++){
//    logfileptr->OFS()<<WChunk(ii,jj)<< " ";
//  }
//  
//  logfileptr->OFS()<<endl;
//}




  TIMER_STOP(Update);
}

};

  void Gathercopy(upcxx::global_ptr<FBMatrix> Objptr, global_ptr<double> Adest, Int j){
      assert(Objptr.tid()==MYTHREAD);
      FBMatrix & A = *Objptr;
      Int local_j = A.global_col_to_local(j);
      Int jb = min(A.blksize, A.n-j);
//      global_ptr<double> Asrc(&A.Achunk(j,local_j));
//      upcxx::ldacopy(A.n-j,jb,Asrc,A.n,Adest,A.n);


      DblNumMat & LocalChunk = A.AchunkLower[local_j/A.blksize];
      global_ptr<double> Asrc = &LocalChunk(0,0);
      upcxx::ldacopy(A.n-j,jb,Asrc,A.n-j,Adest,A.n);

  }



void Factor_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j){
  assert(Aptr.tid()==MYTHREAD);
  FBMatrix & A = *Aptr;

#ifdef _DEBUG_
    logfileptr->OFS()<<"********************************************"<<endl;
    //Factor the column
    logfileptr->OFS()<<"Factoring column "<<j<<endl;
#endif

    A.Factor(j);
  
#ifdef _DEBUG_
    logfileptr->OFS()<<"Factoring column done "<<j<<endl;
#endif
    //Launch the updates
    Int local_j = A.global_col_to_local(j);
    for(Int i = j+A.blksize; i< min(A.n,j+(A.np+1)*A.blksize);i+=A.blksize){
       Int target = A.MAP(i,j);

       //global_ptr<double> AchkPtr(&A.Achunk(0,local_j));
       global_ptr<double> AchkPtr(&A.AchunkLower[local_j/A.blksize](0,0));

#ifdef _DEBUG_
       logfileptr->OFS()<<"Launching update from "<<j<<" on P"<<target<<endl;
#endif
#ifdef LAUNCH_ASYNC
        upcxx::async(target)(Update_Async,(*A.RemoteObjPtr)[target],j,AchkPtr);
#endif
    }

    if(j+A.blksize>=A.n){
#ifdef _DEBUG_
     logfileptr->OFS()<<"Unlocking quit"<<endl;
#endif
     fact_done.unlock(); 
    }

#ifdef _DEBUG_
    logfileptr->OFS()<<"Quitting Factor_Async"<<endl;
    logfileptr->OFS()<<"********************************************"<<endl;
#endif
}


void Aggregate_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, global_ptr<double> remoteAggregatePtr){
  assert(Aptr.tid()==MYTHREAD);
  FBMatrix & A = *Aptr;
  //fetch data
#ifdef _DEBUG_    
  logfileptr->OFS()<<"Aggregating Fetching data from P"<<remoteAggregatePtr.tid()<<endl;
#endif
  Int jb = min(A.blksize, A.n-j);

///#ifdef ADJUSTED_BUFFERS
///   DblNumMat RemoteAggregate(A.n-j,jb); 
///#ifdef ASYNC_COPY
///  TIMER_START(ldacopy_async);
///   upcxx::event e[RemoteAggregate.n()]; 
///   upcxx::ldacopy_async(RemoteAggregate.m(),RemoteAggregate.n(),remoteAggregatePtr,A.n-j,RemoteAggregate.GData(),A.n-j,e);
///   for(int i=0;i<RemoteAggregate.n();i++){e[i].wait();}
///  TIMER_STOP(ldacopy_async);
///#else
///  TIMER_START(ldacopy);
///  upcxx::ldacopy(RemoteAggregate.m(),RemoteAggregate.n(),remoteAggregatePtr,A.n-j,RemoteAggregate.GData(),A.n-j);
///  TIMER_STOP(ldacopy);
///
///#endif
///#else
///   DblNumMat RemoteAggregate(A.n,jb); 
///  TIMER_START(copy);
///    upcxx::copy(remoteAggregatePtr,RemoteAggregate.GData(),RemoteAggregate.m()*RemoteAggregate.n());
///  TIMER_STOP(copy);
///
///#endif



   DblNumMat RemoteAggregate(A.n-j,jb); 
  TIMER_START(copy);
    upcxx::copy(remoteAggregatePtr,RemoteAggregate.GData(),RemoteAggregate.m()*RemoteAggregate.n());
  TIMER_STOP(copy);



#ifdef _DEBUG_    
  logfileptr->OFS()<<"Aggregating update to "<<j<<endl;
#endif
   A.Aggregate_once(j, RemoteAggregate);

#ifdef _DEBUG_    
  logfileptr->OFS()<<"Done Aggregating update to "<<j<<endl;
#endif




   A.AggLock[j/A.blksize]--;
#ifdef _DEBUG_    
   logfileptr->OFS()<<A.AggLock<<endl;
#endif
  if(A.AggLock[j/A.blksize]==0){
    //launch the factorization of column j
#ifdef _DEBUG_    
   logfileptr->OFS()<<"Col "<<j<<" updated. Can be factored "<<endl;
#endif
    
#ifdef LAUNCH_ASYNC
   upcxx::async(A.iam)(Factor_Async,Aptr,j);
#endif
  } 
#ifdef _DEBUG_
    logfileptr->OFS()<<"Quitting Aggregate_Async"<<endl<<endl;
#endif
}


void Update_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, upcxx::global_ptr<double> remoteFactorPtr){
  assert(Aptr.tid()==MYTHREAD);
  FBMatrix & A = *Aptr;


#ifdef _DEBUG_    
  logfileptr->OFS()<<"Fetching factor from P"<<remoteFactorPtr.tid()<<endl;
#endif




//#ifdef ADJUSTED_BUFFERS
//  DblNumMat RemoteFactor(A.n-j,A.blksize);
//#ifdef ASYNC_COPY
//  TIMER_START(ldacopy_async);
//   upcxx::event e[RemoteFactor.n()]; 
//   upcxx::ldacopy_async(RemoteFactor.m(),RemoteFactor.n(),remoteFactorPtr,A.n,RemoteFactor.GData(),A.n-j,e);
//   for(int i=0;i<RemoteFactor.n();i++){e[i].wait();}
//  TIMER_STOP(ldacopy_async);
//#else
//  TIMER_START(ldacopy);
//  upcxx::ldacopy(RemoteFactor.m(),RemoteFactor.n(),remoteFactorPtr,A.n,RemoteFactor.GData(),A.n-j);
//  TIMER_STOP(ldacopy);
//#endif
//
//#else
//  DblNumMat RemoteFactor(A.n,A.blksize);
//  TIMER_START(copy);
//  upcxx::copy(remoteFactorPtr,RemoteFactor.GData(),RemoteFactor.m()*RemoteFactor.n());
//  TIMER_STOP(copy);
//#endif

  DblNumMat RemoteFactor(A.n-j,A.blksize);
  TIMER_START(copy);
  upcxx::copy(remoteFactorPtr,RemoteFactor.GData(),RemoteFactor.m()*RemoteFactor.n());
  TIMER_STOP(copy);

#ifdef _DEBUG_    
  logfileptr->OFS()<<"Fetched factor from P"<<remoteFactorPtr.tid()<<endl;
#endif




  //do all my updates with i>j
  for(Int i = j+A.blksize; i<A.n;i+=A.blksize){
    if(A.iam==A.MAP(i,j)){
#ifdef _DEBUG_
      logfileptr->OFS()<<"Updating column "<<i<<" with columns from P"<<remoteFactorPtr.tid()<<endl;
#endif
      A.Update(j, i, RemoteFactor);

    }
  }

  IntNumVec myLastUpd(A.np);
  SetValue(myLastUpd,-1);
  
  for(Int i = j+A.blksize; i< min(A.n,j+(A.np+1)*A.blksize);i+=A.blksize){
    for(Int j2 = j; j2< i;j2+=A.blksize){
      if(A.iam==A.MAP(i,j2)){
        Int target = A.MAP(i,i);
#ifdef _DEBUG_
        logfileptr->OFS()<<"j2="<<j2<<" vs "<<myLastUpd[target]<<" on P"<<target<<endl;
#endif
        myLastUpd[target] = max(j2,myLastUpd[target]);
      }
    }
  }

#ifdef _DEBUG_
  logfileptr->OFS()<<myLastUpd<<endl;
#endif


  for(Int i = j+A.blksize; i< min(A.n,j+(A.np+1)*A.blksize);i+=A.blksize){
      if(A.iam==A.MAP(i,j)){
        Int target = A.MAP(i,i);
        if(myLastUpd[target]==j){
#ifdef _DEBUG_
            logfileptr->OFS()<<"Updating from "<<j<<", launching aggregate on P"<<target<<" for column "<<i<<endl;
#endif

myLastUpd[target]=-1;
#ifdef LAUNCH_ASYNC
//#ifdef ADJUSTED_BUFFERS
//        global_ptr<double> WbufPtr(&A.W(i,i));
//#else
//        global_ptr<double> WbufPtr(&A.W(0,i));
//#endif

        global_ptr<double> WbufPtr(&A.WLower[i/A.blksize](0,0));

        upcxx::async(target)(Aggregate_Async,(*A.RemoteObjPtr)[target],i,WbufPtr);
#endif
      }
    }
  }
  
#ifdef _DEBUG_
    logfileptr->OFS()<<"Quitting Update_Async"<<endl<<endl;
#endif
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

    logfileptr->OFS()<<"********* LOGFILE OF P"<<iam<<" *********"<<endl;
    logfileptr->OFS()<<"**********************************"<<endl;

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
//    upcxx::global_ptr<FBMatrix> Afactptr;
//    Afactptr = upcxx::allocate<FBMatrix>(iam,1);// upcxx::Create2<FBMatrix>(Afactptr);
    FBMatrix * toto = Afactptr;
      logfileptr->OFS()<<"ptr = "<<(long)toto<<endl;
 


    Aobjects[iam] = Afactptr;
    upcxx::barrier();


    //upcxx::global_ptr<FBMatrix> Afactptr = upcxx::Create<FBMatrix>();

    FBMatrix & Afact = *Afactptr;
    Afact.RemoteObjPtr = &Aobjects;
    Real timeSta, timeEnd;

      DblNumMat RHS;
      DblNumMat XTrue;

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
    csr_matrix_expand_to_dense (A.Data(), 0, Afact.n, (const csr_matrix_t *) Atmp->repr);
    destroy_sparse_matrix (Atmp);



      Int n = Afact.n;
      Int nrhs = 5;

    if(MYTHREAD==0){

      RHS.Resize(n,nrhs);
      XTrue.Resize(n,nrhs);
      UniformRandom(XTrue);
      blas::Gemm('N','N',n,nrhs,n,1.0,&A(0,0),n,&XTrue(0,0),n,0.0,&RHS(0,0),n);
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
    Afact.Allocate(Afact.n,blksize);
    Afact.Distribute(A);

    upcxx::barrier();
    upcxx::wait();

if(MYTHREAD==0){
//      timeSta =  omp_get_wtime( );
//      TIMER_START(POTRF);
//      lapack::Potrf( 'L', A.n(), A.Data(), A.n());
//      TIMER_STOP(POTRF);
//      timeEnd =  omp_get_wtime( );
//      cout<<"REF POTRF: "<<timeEnd-timeSta<<endl;
//
//
//      //check the result
//      double norm = 0;
//
//      //do a solve
//      DblNumMat X = RHS;
//      lapack::Potrs('L',n,nrhs,&A(0,0),n,&X(0,0),n);
//
//
//      blas::Axpy(n*nrhs,-1.0,&XTrue(0,0),1,&X(0,0),1);
//      norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
//      printf("Norm of residual after SOLVE for POTRF is %10.2e\n",norm);
}




A.Clear();



    logfileptr->OFS()<<"distributed"<<endl;

      timeSta =  omp_get_wtime( );
    TIMER_START(FANBOTH);
    if(MYTHREAD==Afact.MAP(Afact.n-1,Afact.n-1)){
#ifdef _DEBUG_
      logfileptr->OFS()<<"LOCK IS TAKEN BY"<<MYTHREAD<<endl;
#endif
      fact_done.lock();
    }

#ifdef _DEBUG_
    logfileptr->OFS()<<"initializing"<<endl;
#endif

    upcxx::barrier();

    if(MYTHREAD==Afact.MAP(0,0)){
      upcxx::async(Afact.iam)(Factor_Async,Afactptr,0);
    }


    while(fact_done.islocked()){upcxx::progress();};
    fact_done.lock();
    TIMER_STOP(FANBOTH);
    timeEnd =  omp_get_wtime( );
    fact_done.unlock();

    upcxx::barrier();

#ifdef _DEBUG_
    logfileptr->OFS()<<"gathering"<<endl;
#endif

    //gather all data on P0
    DblNumMat Afinished;
    if(MYTHREAD==0)
    {
      cout<<"FAN-BOTH: "<<timeEnd-timeSta<<endl;
      Afact.Gather(Afinished);
    }
    upcxx::barrier();
#ifdef _DEBUG_
    logfileptr->OFS()<<"quitting"<<endl;
#endif


    upcxx::barrier();




    upcxx::Destroy<FBMatrix>(Afactptr);



    if(MYTHREAD==0)
    {
      //check the result
      double norm = 0;
//      DblNumMat X = RHS;
//      lapack::Potrs('L',n,nrhs,&A(0,0),n,&X(0,0),n);
//
//
//      blas::Axpy(n*nrhs,-1.0,&XTrue(0,0),1,&X(0,0),1);
//      norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
//      printf("Norm of residual after SOLVE for POTRF is %10.2e\n",norm);
//
      //do a solve
      DblNumMat X = RHS;
      // X = RHS;
      lapack::Potrs('L',n,nrhs,&Afinished(0,0),n,&X(0,0),n);

      blas::Axpy(n*nrhs,-1.0,&XTrue(0,0),1,&X(0,0),1);
      norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
      printf("Norm of residual after SOLVE for FAN-BOTH is %10.2e\n",norm);
    }

    upcxx::finalize();
    delete logfileptr;

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


