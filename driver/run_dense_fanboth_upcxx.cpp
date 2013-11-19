/// @file run_cholesky.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2013-08-31


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

#define MMAP(i,j,arr) arr(((i) % arr.m()),((j) % arr.n())) 
#define MAP(i,j) MMAP(i,j,cmapp)
//#define r(i) ((i)%numCore)
//#define c(i) 0
//#define MAP(i,j) r(i) + c(j)*prow
#define c(i) ((i)%THREADS)
#define r(i) 0
#define MAP(i,j,bs) r(i/bs)*pcol + c(j/bs)

using namespace LIBCHOLESKY;
using namespace std;


namespace LIBCHOLESKY{
class LogFile{
protected:
  std::ofstream myOFS_;

  void init_(const char * prefix, const char * suffix){
    stringstream  ss;
    ss << prefix << suffix;
    myOFS_.open( ss.str().c_str() );
  }
public:

  LogFile(const char * prefix, const char * suffix){
    init_(prefix,suffix);
  }

  LogFile(Int iam){
    std::stringstream ss;
    ss<<iam;
    init_("logTest",ss.str().c_str());
  }


  LogFile(){
  }




  ~LogFile(){
    if(myOFS_.is_open()){
      myOFS_.close();
    }
  }

//  template <typename T>
//  const LogFile& operator<<(const T& obj) 
//  {
//    // write obj to stream
//    mySS_<<obj;
//    myOFS_<<mySS_.str();
//    return *this;
//  }

  template <typename T>
  const std::ofstream& operator<<(const T& obj) 
  {
    // write obj to stream
    std::stringstream ss;
    ss<<obj;
    myOFS_<<ss.str();
    return myOFS_;
  }

  std::ofstream & OFS(){
    return myOFS_;
  }

};

    LogFile * logfileptr;

class FBMatrix{
public:
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

  void Allocate(Int pn, Int pblksize){
    n=pn;
    blksize=pblksize;
      
    Int totBlk = n/blksize;
    Int remaining = n-totBlk*blksize;
    if(remaining>0) totBlk++;
  
    Int localBlk = totBlk/np;
    Int additionalBlock = totBlk - np*localBlk;
    if(iam<additionalBlock){ localBlk++; }
   
    for(Int i=0;i<localBlk;i++){
      Int prevBlk = (iam) +np*(i>0?i-1:0); 
      logfileptr->OFS()<<"block "<<i<<" of size "<<min(n-(prevBlk)*blksize,blksize)<<endl;
    } 
    

    Int prevBlk = (iam) +np*(localBlk-1); 
    Int chksize = (localBlk-1)*blksize + min(n-(prevBlk)*blksize,blksize);
    logfileptr->OFS()<<"chksize="<<chksize<<endl;
    Achunk.Resize(n,chksize);
    SetValue(Achunk,(double)iam);

    //allocate global workspace
    W.Resize(n,n);
    SetValue(W,0.0);

  }


void Aggregate(Int j,Int blksize,Int pcol, Int numCore, DblNumMat * W){
  Int jb = min(blksize, n-j+1);
  //aggregate previous updates
  IntNumVec update_proc(numCore);
  SetValue(update_proc,0);
  for(int i=0;i<j;i+=jb)
  {

    //proc holding update is
    int p = MAP(j,i,blksize);
    if(update_proc(p)==0){
      update_proc(p)=1;

      Aggregate_once(j, W[p]);
    }
  }
}



void Aggregate_once(Int j, DblNumMat &DistW){
  //j is local
  Int local_j = j%(np*blksize);

  Int jb = min(blksize, n-j+1);
  //aggregate previous updates
  TAU_FSTART(Aggregate);

    //proc holding update is
      for(int jj=0;jj<jb && n-j-jj>0;jj++){
        blas::Axpy((n-j-jj), -1.0, &DistW(j+jj,0+jj),1,&Achunk(j+jj,local_j+jj),1);
      }
  TAU_FSTOP(Aggregate);
}




void Factor(Int j, Int blksize, DblNumMat &A){
  Int n= A.n();
  Int jb = min(blksize, n-j+1);

  TAU_FSTART(Factor);
  lapack::Potrf( 'L', jb, &A(j,j ), n);
  if(n-j-jb>0){
    blas::Trsm('R','L','T','N',n-j-jb,jb, 1.0, &A(j,j), n, &A(j+jb,j), n);
  }
  TAU_FSTOP(Factor);
}

void Update(Int j, Int i, DblNumMat & Factor){
  //i is local, j is global
  Int jb = min(blksize, n-j+1);

  TAU_FSTART(Update);
  //call dgemm
  Int jbi = min(blksize, n-i+1);

  blas::Syrk('L','N', jbi ,jb,1.0,&Factor(i,0),n,1.0,&W(i,i),n); 

  if(n-i-jb>0){
    blas::Gemm('N','T',n-i-jb,jb,jb,1.0,&Factor(i+jb,0),n,&Factor(i,0),n,1.0,&W(i+jb,i),n);
  }

  TAU_FSTOP(Update);

}

};











void Aggregate_Async(upcxx::global_ptr<FBMatrix> Aptr, Int i, Int j, Int rc, Int cc, global_ptr<double> remoteAggregatePtr){
    assert(Aptr.tid()==MYTHREAD);
  FBMatrix & A = *Aptr;
  //fetch data
  logfileptr->OFS()<<"Fetching data from col"<<j<<endl;
  Int local_j = j%(A.np*A.blksize);
  logfileptr->OFS()<<"which is "<<local_j<<" locally"<<endl;

//   DblNumMat RemoteAggregate(rc,cc);
//   upcxx::copy(remoteAggregatePtr,RemoteAggregate.GData(),rc*cc);

//   DblNumMat localAggregate(rc,cc);
//   SetValue(localAggregate,0.0);
   //A.Aggregate_once(j, localAggregate,RemoteAggregate);

//   logfileptr->OFS()<<"Remote aggregate "<<RemoteAggregate<<endl;
//   logfileptr->OFS()<<"Local aggregate "<< localAggregate<<endl;

////   A.lock[i]--;
////  //do the sum
////  if(!A.lock[i]){
////    //launch the factorization of column i
    logfileptr->OFS()<<"Col "<<i<<" updated. Can be factored "<<A.AggLock<<endl;
////  } 
}


void Update_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, Int i, Int rc, Int cc, upcxx::global_ptr<double> remoteFactorPtr){
    assert(Aptr.tid()==MYTHREAD);
    FBMatrix & A = *Aptr;

   logfileptr->OFS()<<"Fetching remote factor "<<endl;
   DblNumMat RemoteFactor(rc,cc);
   upcxx::copy(remoteFactorPtr,RemoteFactor.GData(),rc*cc);

   logfileptr->OFS()<<"Remote factor "<<RemoteFactor<<endl;


   A.Update(j, i, RemoteFactor);


   //If it's my last update for P(i,i)
   //Launch Aggregate_async on processor P(i,i)
   Int pcol = A.pcol;
   Int prow = A.prow;
   Int target = (MYTHREAD+1)%THREADS;//MAP(i,i,A.blksize);

    
   upcxx::async(target)(Aggregate_Async,(*A.RemoteObjPtr)[target],i,j,A.n,A.blksize,A.W.GData());

}








}







namespace upcxx{

    template <typename T> void new_fctn(global_ptr<T> gptr){
      new ((T*)gptr) T();
    }
    template <typename T> void delete_fctn(global_ptr<T> gptr){
      ((T*)gptr)->~T();
    }

    template <typename T> global_ptr<T>& Create(int where=MYTHREAD){
      //currently it only works with local assignments
      global_ptr<T> gptr = allocate<T>(where,1);
      async(where)(new_fctn<T>,gptr);
      upcxx::wait();
      return gptr;
    }
    template <typename T> void Destroy(global_ptr<T>& gptr){
      //currently it only works with local assignments
      async(gptr.tid())(delete_fctn<T>,gptr);
      upcxx::wait();
      deallocate<T>(gptr);
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

    //try to create a new object    
    

    int remote = (iam+1)%np;

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


    Real timeSta, timeEnd;

    if(MYTHREAD==0)
    {
      timeSta =  omp_get_wtime( );
      TAU_FSTART(POTRF);
      lapack::Potrf( 'L', D.n(), D.Data(), D.n());
      TAU_FSTOP(POTRF);
      timeEnd =  omp_get_wtime( );
      cout<<"REF POTRF: "<<timeEnd-timeSta<<endl;
    }

    upcxx::barrier();

    //TODO This needs to be fixed
    Afact.prow = sqrt(np);
    Afact.pcol = Afact.prow;
    np = Afact.prow*Afact.pcol;
    if(MYTHREAD==0){
      cout<<"Number of cores to be used: "<<np<<endl;
    }


    Int numStep = ceil((double)Afact.n/(double)Afact.blksize);

    if(MYTHREAD==0){
      cout<<"Factorization will require "<<numStep<<" steps"<<endl;
    }


    //Allocate chunks of the matrix on each processor
    Afact.Allocate(20,blksize);

    //Allocate the lock vector
    Afact.AggLock.Resize(numStep);
    for(int js=0; js<numStep;js++){
        int j = min(js*Afact.blksize,Afact.n-1);
        Int jb = min(Afact.blksize, Afact.n-j+1);

        Int prevBlock = j/Afact.blksize;
        Int prevProc = min(prevBlock?prevBlock/Afact.pcol+1:0,Afact.pcol);


        //test if j is my first block column or if I'm already in prevProc
        Int pcol = Afact.pcol;
      

        logfile<<"j="<<j<<" jb="<<jb<<" prevBlock="<<prevBlock<<" prevProc="<<prevProc<<endl;


        Afact.AggLock[j/Afact.blksize]=prevProc;

    }
 
    logfile<<Afact.AggLock<<endl;




    //setup distributed structure

  
//    if(iam == 0) chksize += n - np*chksize;


//    upcxx::shared_array<upcxx::global_ptr<double> > DistAptr;
//    DistAptr.init(n*n,n*Afact.blksize);


//    DblNumMat AChunk(n,Afact.blksize,false, (upcxx::global_ptr<double>)DistAptr[iam]);
//    SetValue(AChunk,(double)iam);


   logfileptr->OFS()<<"initializing"<<endl;

    upcxx::barrier();


   upcxx::async(remote)(Update_Async, Aobjects[remote],0,9,Afact.n,Afact.blksize,Afact.Achunk.GData());

   upcxx::wait();
   upcxx::barrier();
//   logfileptr->OFS()<<"Launching aggregate on P"<<remote<<endl;


//   upcxx::async(remote)(Aggregate_Async,(upcxx::global_ptr<FBMatrix>)DistAfact[remote],9,0,n,Afact.blksize,Afact.W.GData());

   
//   upcxx::barrier();

//    statusOFS<<"Reading data from P"<<remote<<endl;
//
//    int remoteChksize = n/np;
//    if( remote == 0) remoteChksize += n - np*remoteChksize;
//   DblNumMat RemoteAChunk(n,Afact.blksize);
//   upcxx::copy<double>(DistAptr[remote],RemoteAChunk.GData(),n);
//
//
//
////    logfile<<AChunk<<endl;
//    logfile<<RemoteAChunk<<endl;
//
//    upcxx::wait();
//   upcxx::barrier();

    //remote destroy
    upcxx::Destroy<FBMatrix>(Afactptr);

delete logfileptr;
return;




////
////    //start threads
////    timeSta =  omp_get_wtime( );
////    TAU_FSTART(Fan-Both);
////#pragma omp parallel shared(A,W) num_threads(numCore)
////    {
////      int iam = omp_get_thread_num();
////      //      ofstream statusOFS;
////      //      stringstream  ss;
////      //      ss << "logTest" << iam;
////      //      statusOFS.open( ss.str().c_str() );
////      //      statusOFS<<omp_get_num_threads()<<" P"<<iam<<endl;
////
////      for(int js=0; js<numStep;js++){
////        int j = min(js*blksize,n-1);
////        Int jb = min(blksize, n-j+1);
////
////        //do the aggregation and factorization
////        if(iam==MAP(j,j,blksize) ){
////
////          //aggregate previous updates
////          //TODO UPCXX : Aggregate might be launched asynchrouneously on MAP(j,j,blksize).
////          //             Aggregate_cpy should launch async copy from the source  
////          //             We now the col count of row j, i.e the number of senders.
////          //             When everything is received the sum is computed and Factor is launched.
////          Aggregate(j,blksize,pcol,numCore, A,W);
////
////          //Factor current column 
////          Factor(j,blksize,A);
////        }
////
////        #pragma omp barrier
////
////        //do the updates
////        for(int i=j+jb;i<n ;i+=jb){
////          //compute the update
////          if(iam==MAP(i,j,blksize)){
////            //TODO UPCXX async launch Update on MAP(i,j,blksize)
////            Update(j, i, n,A, W);
////          }
////        }
////
////        #pragma omp barrier
////
////      }  
////
////      //      statusOFS.close();
////
////    }
////
////    //    for(int j=0; j<numStep;j++){ OMP_DESTROY_LOCK(isLocked[j]); }
////    TAU_FSTOP(Fan-Both);
////
////    timeEnd =  omp_get_wtime( );
////    cout<<"Fan-Both: "<<timeEnd-timeSta<<endl;
////
////
////    //check the result
////    double norm = 0;
////
////    //do a solve
////    Int nrhs = n;
////    DblNumMat RHS(n,nrhs);
////    DblNumMat XTrue(n,nrhs);
////    UniformRandom(XTrue);
////    blas::Gemm('N','N',n,nrhs,n,1.0,&RefA(0,0),n,&XTrue(0,0),n,0.0,&RHS(0,0),n);
////
////    DblNumMat X = RHS;
////    lapack::Potrs('L',n,nrhs,&D(0,0),n,&X(0,0),n);
////    blas::Axpy(n*nrhs,-1.0,&XTrue(0,0),1,&X(0,0),1);
////    norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
////    printf("Norm of residual after SOLVE for REF POTF2 is %10.2e\n",norm);
////
////    X = RHS;
////    lapack::Potrs('L',n,nrhs,&A(0,0),n,&X(0,0),n);
////    blas::Axpy(n*nrhs,-1.0,&XTrue(0,0),1,&X(0,0),1);
////    norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
////    printf("Norm of residual after SOLVE for FAN-BOTH is %10.2e\n",norm);
////
////
////
////
////    X = RHS;
////    lapack::Posv('L',n,nrhs,&RefA(0,0),n,&X(0,0),n);
////    blas::Axpy(n*nrhs,-1.0,&XTrue(0,0),1,&X(0,0),1);
////    norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
////    printf("Norm of residual after POSV is %10.2e\n",norm);
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

