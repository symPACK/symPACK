/// @file run_cholesky.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2013-08-31

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

extern "C"{

void Cblacs_pinfo( int* mypnum, int* nprocs);
void Cblacs_get( int context, int request, int* value);
int  Cblacs_gridinit( int* context, char* order, int np_row, int np_col);
void Cblacs_gridinfo( int context, int* np_row, int* np_col, int* my_row, int* my_col);
void Cblacs_gridexit( int context);
void Cblacs_exit( int error_code);
void Cblacs_barrier(int context, char* scope);
void Cblacs_pcoord(int ctxt, int myid, int * myrow, int * mycol);

void Cdgesd2d(int ctxt,int nr,int nc,double * A, int n, int sendr,int  sendc);
void Cdgerv2d(int ctxt,int nr,int nc,double * A_loc, int nrows, int dum, int my);
int numroc_(int *n, int *blksize, int *myrow, int *iZERO, int *procrows);

void descinit_(int * descA, int * m, int * n, int * mb, int *nb, int * irsrc, int * icsrc, int * ictxt, int * lld, int * info );
void pdgeadd_( char * trans, int *m, int *n, double * alpha, double * A, int * ia, int * ja, int * desca, double * beta, double * C, int * ic, int * jc, int * descc );
void pdpotrf_(char * uplo, int * n, double * A, int * ia, int * ja, int * desca, int * info);
}


using namespace LIBCHOLESKY;
using namespace std;

LogFile * logfileptr;

int main(int argc, char **argv) 
{

  int mpirank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);


  try{


    logfileptr = new LogFile(mpirank);
    LogFile & logfile = *logfileptr;

    logfile<<"********* LOGFILE OF P"<<mpirank<<" *********"<<endl;
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


    double * A=NULL; 
    double * Afact=NULL; 
    Int n;  
    if(mpirank==0){
    //Read the input matrix
    sparse_matrix_file_format_t informat;
    informat = sparse_matrix_file_format_string_to_enum (informatstr.c_str());
    sparse_matrix_t* Atmp = load_sparse_matrix (informat, filename.c_str());
    sparse_matrix_expand_symmetric_storage (Atmp);
    sparse_matrix_convert (Atmp, CSR);

    n = ((csr_matrix_t *)Atmp->repr)->n;

    A=new double[n*n];
      cout<<"Matrix order is "<<n<<endl;

    csr_matrix_expand_to_dense (A, 0, n, Atmp->repr);
    destroy_sparse_matrix (Atmp);

    Afact=new double[n*n];

    }

int ctxt, myid, myrow, mycol, numproc;
Cblacs_pinfo(&myid, &numproc);
int procrows = sqrt(numproc), proccols = sqrt(numproc);
Cblacs_get(0, 0, &ctxt);
Cblacs_gridinit(&ctxt, "Row-major", procrows, proccols);
Cblacs_pcoord(ctxt, myid, &myrow, &mycol);





/* Broadcast of the matrix dimensions */
MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
 
/* Reserve space for local matrices */
// Number of rows and cols owned by the current process
int iZERO = 0;
double one = 1.0;
double zero = 0.0;
int i_one = 1;
int Nb = ceil((double)n/(double)blksize);
int nrows = numroc_(&n, &blksize, &myrow, &iZERO, &procrows);
int ncols = numroc_(&n, &blksize, &mycol, &iZERO, &proccols);
for (int id = 0; id < numproc; ++id) {
    Cblacs_barrier(ctxt, "All");
}


double * A_loc = new double[nrows*ncols];
for (int i = 0; i < nrows*ncols; ++i) *(A_loc+i)=0.;

int info =0;
int descA_distr[9],descA[9];
// Initialize discriptors (local matrix A is considered as distributed with blocking parameters
// m, n,i.e. there is only one block - whole matrix A - which is located on process (0,0) )
int lld = MAX( numroc_( &n, &n, &myrow, &iZERO, &procrows ), 1 );
descinit_( descA, &n, &n, &n, &n, &iZERO, &iZERO, &ctxt, &lld, &info );
int lld_distr = MAX( nrows, 1 );
descinit_( descA_distr, &n, &n, &blksize, &blksize, &iZERO, &iZERO, &ctxt, &lld_distr, &info );

// Call pdgeadd_ to distribute matrix (i.e. copy A into A_loc)
pdgeadd_( "N", &n, &n, &one, A, &i_one, &i_one, descA, &zero, A_loc, &i_one, &i_one, descA_distr );



double timeSta, timeEnd;
      timeSta =  MPI_Wtime( );
pdpotrf_("L",&n,A_loc,&i_one,&i_one,descA_distr,&info);
      timeEnd =  MPI_Wtime( );

if(mpirank==0){
      cout<<"REF DPOTRF: "<<timeEnd-timeSta<<endl;
}


// Copy result into local matrix
pdgeadd_( "N", &n, &n, &one, A_loc, &i_one, &i_one, descA_distr, &zero, Afact, &i_one, &i_one, descA );








//
////Scattering the matrix
//
///* Scatter matrix */
//int sendr = 0, sendc = 0, recvr = 0, recvc = 0;
//for (int r = 0; r < n; r += blksize, sendr=(sendr+1)%procrows) {
//    sendc = 0;
//    // Number of rows to be sent
//    // Is this the last row block?
//    int nr = min(n-r,blksize);
// 
//    for (int c = 0; c < n; c += blksize, sendc=(sendc+1)%proccols) {
//        // Number of cols to be sent
//        // Is this the last col block?
//        int nc = min(n-c,blksize);
// 
//        if (mpirank==0) {
//            // Send a nr-by-nc submatrix to process (sendr, sendc)
//            Cdgesd2d(ctxt, nr, nc, A+n*c+r, n, sendr, sendc);
//        }
// 
//        if (myrow == sendr && mycol == sendc) {
//            // Receive the same data
//            // The leading dimension of the local matrix is nrows!
//            Cdgerv2d(ctxt, nr, nc, A_loc+nrows*recvc+recvr, nrows, 0, 0);
//            recvc = (recvc+nc)%ncols;
//        }
// 
//    }
// 
//    if (myrow == sendr)
//        recvr = (recvr+nr)%nrows;
//}
//
//
//
////call pdpotrf
//
//char UPLO = 'L';
//pdpotrf_(&UPLO,&n,A_loc,,,,&info);
//
//
//
///* Gather matrix */
//sendr = 0;
//for (int r = 0; r < n; r += blksize, sendr=(sendr+1)%procrows) {
//    sendc = 0;
//    // Number of rows to be sent
//    // Is this the last row block?
//    int nr = min(n-r,blksize);
// 
//    for (int c = 0; c < n; c += blksize, sendc=(sendc+1)%proccols) {
//        // Number of cols to be sent
//        // Is this the last col block?
//        int nc = min(n-c,blksize);
// 
//        if (myrow == sendr && mycol == sendc) {
//            // Send a nr-by-nc submatrix to process (sendr, sendc)
//            Cdgesd2d(ctxt, nr, nc, A_loc+nrows*recvc+recvr, nrows, 0, 0);
//            recvc = (recvc+nc)%ncols;
//        }
// 
//        if (mpirank==0) {
//            // Receive the same data
//            // The leading dimension of the local matrix is nrows!
//            Cdgerv2d(ctxt, nr, nc, Afact+n*c+r, n, sendr, sendc);
//        }
// 
//    }
// 
//    if (myrow == sendr)
//        recvr = (recvr+nr)%nrows;
//}

if(mpirank==0){
      //check the result
      double norm = 0;

      //do a solve
      Int nrhs = 5;
      double RHS[n*nrhs];
      double XTrue[n*nrhs];


      UniformRandom(XTrue,n*nrhs);

      blas::Gemm('N','N',n,nrhs,n,1.0,A,n,XTrue,n,0.0,RHS,n);

      double X[n*nrhs];
      std::copy(RHS,RHS+n*nrhs,X);

      lapack::Potrs('L',n,nrhs,Afact,n,X,n);
      blas::Axpy(n*nrhs,-1.0,XTrue,1,X,1);
      norm = lapack::Lange('F',n,nrhs,X,n);
      printf("Norm of residual after SOLVE for REF DPOTRF is %10.2e\n",norm);

    }


/* Release resources */
delete[] A;
delete[] Afact;
delete[] A_loc;

delete logfileptr;

Cblacs_gridexit(ctxt);
MPI_Finalize();


  }
  catch( std::exception& e )
  {
    std::cerr << "Exception with message: "
      << e.what() << std::endl;
  }


  return 0;
}


