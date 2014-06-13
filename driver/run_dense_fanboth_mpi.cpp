/// @file run_cholesky.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2013-08-31
//#define _DEBUG_


#include <time.h>

#ifndef __PGI
#include <random>
#endif

#include <omp.h>

#include  "FBMatrix.hpp"
#include  "FBMatrix_mpi.hpp"
#include  "Environment.hpp"
#include  "NumVec.hpp"
#include  "NumMat.hpp"
#include  "utility.hpp"
#include  "ETree.hpp"
#include  "blas.hpp"
#include  "lapack.hpp"
#include  "LogFile.hpp"

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


using namespace LIBCHOLESKY;
using namespace std;





int main(int argc, char **argv) 
{

  MPI_Init(&argc,&argv);
  int np, iam;
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Comm_rank(MPI_COMM_WORLD,&iam);

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


  FBMatrix_mpi * Afactptr = new FBMatrix_mpi();

  MPIGrid * grid = new MPIGrid(MPI_COMM_WORLD,1,np);
  Afactptr->Initialize(*grid);

  Real timeSta, timeEnd;

  DblNumMat RHS;
  DblNumMat XTrue;

  //Read the input matrix
  sparse_matrix_file_format_t informat;
  informat = sparse_matrix_file_format_string_to_enum (informatstr.c_str());
  sparse_matrix_t* Atmp = load_sparse_matrix (informat, filename.c_str());
  sparse_matrix_expand_symmetric_storage (Atmp);
  sparse_matrix_convert (Atmp, CSR);
  Afactptr->n = ((csr_matrix_t *)Atmp->repr)->n;

  if(iam==0){
    cout<<"Matrix order is "<<Afactptr->n<<endl;
  }


  DblNumMat A(Afactptr->n,Afactptr->n);
  csr_matrix_expand_to_dense (A.Data(), 0, Afactptr->n, (const csr_matrix_t *) Atmp->repr);
  destroy_sparse_matrix (Atmp);



  Int n = Afactptr->n;
  Int nrhs = 5;
  if(iam==0){
    RHS.Resize(n,nrhs);
    XTrue.Resize(n,nrhs);
    UniformRandom(XTrue);
    blas::Gemm('N','N',n,nrhs,n,1.0,&A(0,0),n,&XTrue(0,0),n,0.0,&RHS(0,0),n);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if(iam==0){
    cout<<"Number of cores to be used: "<<np<<endl;
  }



  //Allocate chunks of the matrix on each processor
  Afactptr->Allocate(np,Afactptr->n,blksize);
  Afactptr->Distribute(A);

#ifdef DRAW_GRAPH
  Afactptr->Draw_Graph();
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  A.Clear();

  logfileptr->OFS()<<"distributed"<<endl;

  timeSta =  omp_get_wtime( );
  TIMER_START(FANBOTH);
  //MPI_Barrier(MPI_COMM_WORLD);
  Afactptr->NumericalFactorization();
  Afactptr->WaitFactorization();
  TIMER_STOP(FANBOTH);
  timeEnd =  omp_get_wtime( );




  //gather all data on P0
  DblNumMat Afinished;

  MPI_Barrier(MPI_COMM_WORLD);
  if(iam==0)
  {
    cout<<"FAN-BOTH: "<<timeEnd-timeSta<<endl;
  }


  Afactptr->Gather(Afinished);

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
  delete Afactptr;
  delete grid;

  MPI_Finalize();

  return 0;
}


