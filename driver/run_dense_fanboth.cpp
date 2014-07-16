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
#define c(i) ((i)%numCore)
#define r(i) 0
#define MAP(i,j,bs) r(i/bs)*pcol + c(j/bs)

using namespace LIBCHOLESKY;
using namespace std;



void Usage(){
  std::cout << "Usage" << std::endl << "run_pselinv -T [isText] -F [doFacto -E [doTriSolve] -Sinv [doSelInv]]  -H <Hfile> -S [Sfile] -colperm [colperm] -npsymbfact [npsymbfact] -P [maxpipelinedepth] -SinvBcast [doSelInvBcast] -SinvPipeline [doSelInvPipeline] -SinvOriginal [doSelInvOriginal] -Shift [imaginary shift] " << std::endl;
}

int main(int argc, char **argv) 
{

  if( argc < 3 ) {
    Usage();
    return 0;
  }

#if defined(PROFILE) || defined(PMPI)
  TAU_PROFILE_INIT(argc, argv);
#endif


  try{

    // *********************************************************************
    // Input parameter
    // *********************************************************************
    std::map<std::string,std::string> options;

    OptionsCreate(argc, argv, options);

    Int n;

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

    //Read the input matrix
    sparse_matrix_file_format_t informat;
    informat = sparse_matrix_file_format_string_to_enum (informatstr.c_str());
    sparse_matrix_t* Atmp = load_sparse_matrix (informat, filename.c_str());
    sparse_matrix_expand_symmetric_storage (Atmp);
    sparse_matrix_convert (Atmp, CSR);
    n = ((csr_matrix_t *)Atmp->repr)->n;

    cout<<"Matrix order is "<<n<<endl;

    DblNumMat A(n,n);
    csr_matrix_expand_to_dense (A.Data(), 0, n, Atmp->repr);
    destroy_sparse_matrix (Atmp);


    DblNumMat RefA = A;
    DblNumMat D = A;

    Real timeSta, timeEnd;

    timeSta =  omp_get_wtime( );
//    TAU_FSTART(POTF2);
//    lapack::Potf2( 'L', n, D.Data(), n);
//    TAU_FSTOP(POTF2);
//    timeEnd =  omp_get_wtime( );
//    cout<<"REF POTF2: "<<timeEnd-timeSta<<endl;
    TAU_FSTART(POTRF);
    lapack::Potrf( 'L', n, D.Data(), n);
    TAU_FSTOP(POTRF);
    timeEnd =  omp_get_wtime( );
    cout<<"REF POTRF: "<<timeEnd-timeSta<<endl;


    //generate a 2d mapping of computations
    Int numCore = 0;
    #pragma omp parallel
    {
      #pragma omp master
      {
        numCore=omp_get_num_threads();
      } 
    }



    Int prow = sqrt(numCore);
    Int pcol = prow;
    numCore = pcol*prow;

    cout<<"Number of cores to be used: "<<numCore<<endl;

    //allocate global workspace
    DblNumMat W(n,numCore*n);
    SetValue(W,0.0);

    Int numStep = ceil((double)n/(double)blksize);

    cout<<"Factorization will require "<<numStep<<" steps "<< n <<"/"<<blksize<<endl;

    //start threads
    timeSta =  omp_get_wtime( );
    TAU_FSTART(Fan-Both);
    #pragma omp parallel shared(A,W) num_threads(numCore)
    {
      int iam = omp_get_thread_num();
//      ofstream statusOFS;
//      stringstream  ss;
//      ss << "logTest" << iam;
//      statusOFS.open( ss.str().c_str() );
//      statusOFS<<omp_get_num_threads()<<" P"<<iam<<endl;

      for(int js=0; js<numStep;js++){
        int j = min(js*blksize,n-1);
        Int jb = min(blksize, n-j+1);

        //do the aggregation and factorization
        if(iam==MAP(j,j,blksize) ){

          //aggregate previous updates
          TAU_FSTART(Aggregate);
          IntNumVec update_proc(numCore);
          SetValue(update_proc,0);
          for(int i=0;i<j;i+=jb)
          {

            //proc holding update is
            int p = MAP(j,i,blksize);
            if(update_proc(p)==0){
              update_proc(p)=1;

              for(int jj=0;jj<jb && n-j-jj>0;jj++){
                  blas::Axpy((n-j-jj), -1.0, &W(j+jj,p*n+j+jj),1,&A(j+jj,j+jj),1);
              }
            }
          }
          TAU_FSTOP(Aggregate);

          //Factor current column 
          TAU_FSTART(Factor);
          lapack::Potrf( 'L', jb, &A(j,j ), n);
          if(n-j-jb>0){
            blas::Trsm('R','L','T','N',n-j-jb,jb, 1.0, &A(j,j), n, &A(j+jb,j), n);
          }




          TAU_FSTOP(Factor);
        }

        #pragma omp barrier

        //do the updates

        for(int i=j+jb;i<n ;i+=jb){
          //compute the update
          if(iam==MAP(i,j,blksize)){
            TAU_FSTART(Update);
            //call dgemm


//dgemm( 'N', 'T', n-j-jb+1, jb,j-1, -1.0, a( j+jb, 1 ), n, a( j, 1 ), n, 1.0, a( j+jb, j ), n );
            //DblNumMat factors(jb,jb);
            //SetValue(factors,0.0);
            //for(int k=0;k<jb;k++){ factors(k,k) = A(i,j+k); } 
            //blas::Gemm('N','N',n-i,jb,jb,1.0,&A(i,j),n,&factors(0,0),jb,1.0,&W(i,iam*n+i),n);

            Int jbi = min(blksize, n-i+1);
            blas::Syrk('L','N', jbi ,jb,1.0,&A(i,j),n,1.0,&W(i,iam*n+i),n); 

            if(n-i-jb>0){
              blas::Gemm('N','T',n-i-jb,jb,jb,1.0,&A(i+jb,j),n,&A(i,j),n,1.0,&W(i+jb,iam*n+i),n);
            }

            TAU_FSTOP(Update);
          }
        }

        #pragma omp barrier

      }  

//      statusOFS.close();

    }

//    for(int j=0; j<numStep;j++){ OMP_DESTROY_LOCK(isLocked[j]); }
    TAU_FSTOP(Fan-Both);

    timeEnd =  omp_get_wtime( );
    cout<<"Fan-Both: "<<timeEnd-timeSta<<endl;

//    {
//      stringstream  ss;
//      ss << "logTest" << 0;
//      ofstream statusOFS;
//      statusOFS.open( ss.str().c_str(), std::ofstream::app);
//      //A = L L^T
//      statusOFS<<"Global POTRF:"<<D<<endl;
//      statusOFS<<"FanBOTH:"<<A<<endl;
//      statusOFS.close();
//    }



    //check the result
    double norm = 0;

    //do a solve
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
    lapack::Potrs('L',n,nrhs,&A(0,0),n,&X(0,0),n);
    blas::Axpy(n*nrhs,-1.0,&XTrue(0,0),1,&X(0,0),1);
    norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
    printf("Norm of residual after SOLVE for FAN-BOTH is %10.2e\n",norm);




    X = RHS;
    lapack::Posv('L',n,nrhs,&RefA(0,0),n,&X(0,0),n);
    blas::Axpy(n*nrhs,-1.0,&XTrue(0,0),1,&X(0,0),1);
    norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
    printf("Norm of residual after POSV is %10.2e\n",norm);
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

