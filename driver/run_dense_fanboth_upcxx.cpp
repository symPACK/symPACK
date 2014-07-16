/// @file run_dense_fanboth_upcxx.cpp
/// @brief Test for the dense version of fanboth using upcxx.
/// @author Mathias Jacquelin
/// @date 2013-08-31
//#define _DEBUG_
#define LAUNCH_ASYNC
#define LAUNCH_FACTOR

//#define ADJUSTED_BUFFERS
//#define ASYNC_COPY

#define UPCXX

#include <time.h>
//#include <random>

#include  "Environment.hpp"
#include  "NumVec.hpp"
#include  "NumMat.hpp"
#include  "NumMat_upcxx.hpp"
#include  "utility.hpp"
#include  "ETree.hpp"
#include  "blas.hpp"
#include  "lapack.hpp"
#include  "FBMatrix.hpp"
#include  "FBMatrix_upcxx.hpp"
#include  "LogFile.hpp"

#include <upcxx.h>
#include <upcxx_runtime.h>
#include <async.h>
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


#include "timer.hpp"

using namespace LIBCHOLESKY;
using namespace std;



int main(int argc, char **argv) 
{

  upcxx::init(&argc,&argv);

  int np = THREADS;
  int iam = MYTHREAD;

#if defined(PROFILE) || defined(PMPI)
  TAU_PROFILE_INIT(argc, argv);
#endif

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
  Int maxPrefetch = 0;
  if( options.find("-pref") != options.end() ){
    maxPrefetch= atoi(options["-pref"].c_str());
  }
  bool loop= false;
  if( options.find("-loop") != options.end() ){
    loop = (atoi(options["-loop"].c_str())==1);
  }



  FBMatrix_upcxx * Afactptr = new FBMatrix_upcxx();

  upcxx::barrier();

  Afactptr->Initialize();

  Real timeSta, timeEnd;



  DblNumMat RHS;
  DblNumMat XTrue;

  //Read the input matrix in sparse format
  sparse_matrix_file_format_t informat;
  informat = sparse_matrix_file_format_string_to_enum (informatstr.c_str());
  sparse_matrix_t* Atmp = load_sparse_matrix (informat, filename.c_str());
  sparse_matrix_expand_symmetric_storage (Atmp);
  sparse_matrix_convert (Atmp, CSR);
  Afactptr->n = ((csr_matrix_t *)Atmp->repr)->n;


  if(iam==0){
    cout<<"Matrix order is "<<Afactptr->n<<endl;
  }

  //Convert the sparse matrix in dense format
  DblNumMat_upcxx A;
  A.Resize(Afactptr->n,Afactptr->n);
  csr_matrix_expand_to_dense (A.Data(), 0, Afactptr->n, (const csr_matrix_t *) Atmp->repr);
  destroy_sparse_matrix (Atmp);

  if(iam==0){
    cout<<"Matrix expanded to dense"<<endl;
  }


  //Create an exact solution and the corresponding RHS
  Int n = Afactptr->n;
  Int nrhs = 5;
  if(iam==0){
    RHS.Resize(n,nrhs);
    XTrue.Resize(n,nrhs);
    UniformRandom(XTrue);
    blas::Gemm('N','N',n,nrhs,n,1.0,&A(0,0),n,&XTrue(0,0),n,0.0,&RHS(0,0),n);
  }

  upcxx::barrier();


  //Setup the process grid
  Afactptr->prefetch = maxPrefetch;
//  Afactptr->prow = 1;
//  Afactptr->pcol = np;
//  np = Afactptr->prow*Afactptr->pcol;


  //Allocate chunks of the matrix on each processor

//gdb_lock();

  Afactptr->Allocate(np,Afactptr->n,blksize);


  if(iam==0){
    cout<<"Number of cores to be used: "<<np<<endl;
  }


  Afactptr->Distribute(A);
  if(iam==0){
    cout<<"FBMatrix distributed"<<endl;
  }

  upcxx::barrier();
  upcxx::wait();
  //Clear the dense matrix A because now Afactptr holds the data.
  A.Clear();


  timeSta =  get_time( );
  TIMER_START(FANBOTH);
  if(iam==Afactptr->MAP(Afactptr->n-1,Afactptr->n-1)){
#ifdef _DEBUG_
    cout<<"Last block column is processed by P"<<iam<<endl;
#endif
  }
  upcxx::barrier();

  //Start the factorization
  TIMER_START(FANBOTH_NUMFACT);
  if(loop){
    Afactptr->NumericalFactorizationLoop();
  }
  else{
    Afactptr->NumericalFactorization();
  }
  TIMER_STOP(FANBOTH_NUMFACT);

  //Wait until the end of the factorization
  TIMER_START(FANBOTH_WAITFACT);
  Afactptr->WaitFactorization();
  TIMER_STOP(FANBOTH_WAITFACT);
  timeEnd =  get_time( );

  TIMER_STOP(FANBOTH);


#ifdef TIMER_QUEUE
  upcxx::print_timer_queue();
#endif


  Afactptr->ClearTmp();

  //gather all data on P0
  DblNumMat Afinished;

  upcxx::barrier();
  logfileptr->OFS()<<"aggregate_comm_time: "<<Afactptr->aggregate_comm_time<<endl;
  logfileptr->OFS()<<"factor_comm_time: "<<Afactptr->factor_comm_time<<endl;
  if(iam==0) {
    cout<<"FAN-BOTH: "<<timeEnd-timeSta<<endl;
  }

  //Gather the distributed factored matrix on P0
  if(iam==0){
    Afactptr->Gather(Afinished);
  }
  upcxx::barrier();

  //Now check the result
  if(iam==0)
  {
    //check the result
    double norm = 0;

    //do a solve
    DblNumMat X = RHS;
    // X = RHS;
    lapack::Potrs('L',n,nrhs,&Afinished(0,0),n,&X(0,0),n);
    blas::Axpy(n*nrhs,-1.0,&XTrue(0,0),1,&X(0,0),1);
    norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
    printf("Norm of residual after SOLVE for FAN-BOTH is %10.2e\n",norm);
  }

  delete logfileptr;
  delete Afactptr;
  upcxx::finalize();

  return 0;
}


