/// @file run_sparse_fanboth_upcxx.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2014-01-23
//#define _DEBUG_


#include <upcxx.h>

#include <time.h>
#include <random>
#include <omp.h>

#include  "Environment.hpp"
#include  "DistSparseMatrix.hpp"
#include  "NumVec.hpp"
#include  "NumMat.hpp"
#include  "utility.hpp"
#include  "ETree.hpp"
#include  "blas.hpp"
#include  "lapack.hpp"
#include  "FBMatrix.hpp"

#include <async.h>
#include  "LogFile.hpp"
#include  "upcxx_additions.hpp"


extern "C" {
#include <bebop/util/config.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/csc_matrix.h>
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


    Aobjects[iam] = Afactptr;
    upcxx::barrier();
    upcxx::wait();
  
//    for(int i =0;i<Aobjects.size();i++){
//      logfileptr->OFS()<<(upcxx::global_ptr<FBMatrix> )Aobjects[i]<<endl;
//    }


    //upcxx::global_ptr<FBMatrix> Afactptr = upcxx::Create<FBMatrix>();

    FBMatrix & Afact = *Afactptr;
    Afact.Initialize(&Aobjects);

    Real timeSta, timeEnd;

      DblNumMat RHS;
      DblNumMat XTrue;

    //Read the input matrix
    sparse_matrix_file_format_t informat;
    informat = sparse_matrix_file_format_string_to_enum (informatstr.c_str());
    sparse_matrix_t* Atmp = load_sparse_matrix (informat, filename.c_str());
    sparse_matrix_expand_symmetric_storage (Atmp);
    sparse_matrix_convert (Atmp, CSC);
    Afact.n = ((csc_matrix_t *)Atmp->repr)->n;

    if(MYTHREAD==0){
      cout<<"Matrix order is "<<Afact.n<<endl;
    }


    DistSparseMatrix<Real> HMat;

    const csc_matrix_t * cscptr = (const csc_matrix_t *) Atmp->repr;
    //fill global structure info as we have it directly
    HMat.size = cscptr->n; 
    HMat.nnz = cscptr->nnz; 
    HMat.Global_.colptr.Resize(cscptr->n);
    std::copy(cscptr->colptr,cscptr->colptr+cscptr->n,HMat.Global_.colptr.Data());
    HMat.Global_.rowind.Resize(cscptr->nnz);
    std::copy(cscptr->rowidx,cscptr->rowidx+cscptr->nnz,HMat.Global_.rowind.Data());
    //Compute local structure info

logfileptr->OFS()<<HMat.Global_.rowind<<endl;
return;

	  // Compute the number of columns on each processor
	  IntNumVec numColLocalVec(np);
	  Int numColLocal, numColFirst;
	  numColFirst = HMat.size / np;
    SetValue( numColLocalVec, numColFirst );
    numColLocalVec[np-1] = HMat.size - numColFirst * (np-1);  // Modify the last entry	
  	numColLocal = numColLocalVec[np];

	  HMat.Local_.colptr.Resize( numColLocal + 1 );
  	for( Int i = 0; i < numColLocal + 1; i++ ){
	  	HMat.Local_.colptr[i] = HMat.Global_.colptr[np * numColFirst+i] - HMat.Global_.colptr[np * numColFirst] + 1;
	  }

	  // Calculate nnz_loc on each processor
	  HMat.Local_.nnz = HMat.Local_.colptr[numColLocal] - HMat.Local_.colptr[0];
    // Resize rowind and nzval appropriately 
    HMat.Local_.rowind.Resize( HMat.Local_.nnz );
	  HMat.nzvalLocal.Resize ( HMat.Local_.nnz );



    //Read my row indices
    Int prevRead = 0;
    Int numRead = 0;
		for( Int ip = 0; ip <iam; ip++ ){	
      prevRead += HMat.Global_.colptr[ip*numColFirst + numColLocalVec[ip]]
                     - HMat.Global_.colptr[ip*numColFirst];
    }

		numRead = HMat.Global_.colptr[iam*numColFirst + numColLocalVec[iam]] - 
			HMat.Global_.colptr[iam*numColFirst];
    std::copy(&HMat.Global_.rowind[prevRead],&HMat.Global_.rowind[prevRead+numRead],HMat.Local_.rowind.Data());


    return;
//
//    //copy appropriate nnz values
//    if(cscptr->value_type == REAL){
//      std::copy(&((const double*)cscptr->values)[prevRead],&((const double*)cscptr->values)[prevRead+numRead],HMat.nzvalLocal.Data());
//    }
//    else if(cscptr->value_type == COMPLEX){
//      std::copy(&((const double*)cscptr->values)[prevRead],&((const double*)cscptr->values)[prevRead+numRead],HMat.nzvalLocal.Data());
//    }
  
    return;

    DblNumMat A(Afact.n,Afact.n);
    //csc_matrix_expand_to_dense (A.Data(), 0, Afact.n, (const csc_matrix_t *) Atmp->repr);
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
    Afact.prow = 1;//sqrt(np);
    Afact.pcol = np;//Afact.prow;
    np = Afact.prow*Afact.pcol;
    if(MYTHREAD==0){
      cout<<"Number of cores to be used: "<<np<<endl;
    }



    //Allocate chunks of the matrix on each processor
    Afact.Allocate(Afact.n,blksize);
    Afact.Distribute(A);

    upcxx::barrier();
    upcxx::wait();
    A.Clear();

    logfileptr->OFS()<<"distributed"<<endl;

      timeSta =  omp_get_wtime( );
    TIMER_START(FANBOTH);
    if(MYTHREAD==Afact.MAP(Afact.n-1,Afact.n-1)){
#ifdef _DEBUG_
      cout<<"LOCK IS TAKEN BY"<<MYTHREAD<<endl;
#endif
    }

#ifdef _DEBUG_
    logfileptr->OFS()<<"initializing"<<endl;
#endif

    upcxx::barrier();

    Afact.NumericalFactorization();
    Afact.WaitFactorization();
    timeEnd =  omp_get_wtime( );

    TIMER_STOP(FANBOTH);

#ifdef _DEBUG_
    logfileptr->OFS()<<"gathering"<<endl;
#endif

    //gather all data on P0
    DblNumMat Afinished;
    upcxx::barrier();
    if(MYTHREAD==0)
    {
      cout<<"FAN-BOTH: "<<timeEnd-timeSta<<endl;
      Afact.Gather(Afinished);
    }

    upcxx::wait();

//
//    logfileptr->OFS()<<"Factor_Async is "<<(void *) Factor_Async<<endl;
//    logfileptr->OFS()<<"Update_Async is "<<(void *) Update_Async<<endl;
//    logfileptr->OFS()<<"Aggregate_Async is "<<(void *) Aggregate_Async<<endl;
//
//    logfileptr->OFS()<<"Waiting time in async queues:"<<endl;
//    for(int i = 0; i<queue_time.size(); i++){
//      logfileptr->OFS()<<queue_time[i]<<" ("<<queue_fptr[i]<<") "; 
//    }
//    logfileptr->OFS()<<endl;


#ifdef _DEBUG_
    logfileptr->OFS()<<"quitting "<<iam<<endl;
#endif

    if(iam==0)
    {
    logfileptr->OFS()<<"testing"<<endl;
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

    upcxx::finalize();
    upcxx::Destroy<FBMatrix>(Afactptr);
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


