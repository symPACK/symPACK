/// @file run_sparse_fanboth.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2014-01-23
//#define _DEBUG_


#include <mpi.h>

#include <time.h>
#include <random>
#include <omp.h>

#include  "Environment.hpp"
#include  "DistSparseMatrix.hpp"
#include  "NumVec.hpp"
#include  "NumMat.hpp"
#include  "utility.hpp"
#include  "ETree.hpp"
#include  "Mapping.hpp"

#include  "blas.hpp"
#include  "lapack.hpp"
#include  "NZBlock.hpp"
#include  "SuperNode.hpp"
#include  "SupernodalMatrix.hpp"

#include <async.h>
#include  "LogFile.hpp"


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
  MPI_Init(&argc,&argv);

    MPI_Comm worldcomm;
    MPI_Comm_dup(MPI_COMM_WORLD,&worldcomm);
  int np, iam;
  MPI_Comm_size(worldcomm,&np);
  MPI_Comm_rank(worldcomm,&iam);

    MAPCLASS mapping(np, sqrt(np), np, 1);



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


    Real timeSta, timeEnd;


    //Read the input matrix
    sparse_matrix_file_format_t informat;
    informat = sparse_matrix_file_format_string_to_enum (informatstr.c_str());


    Int nrhs = 5;

    DblNumMat STEPS;
    DblNumMat RHS;
    DblNumMat XTrue;
    DblNumMat A;
    {
      if(iam==0){
        sparse_matrix_t* Atmp = load_sparse_matrix (informat, filename.c_str());
        sparse_matrix_expand_symmetric_storage (Atmp);
        sparse_matrix_convert (Atmp, CSR);
        const csr_matrix_t * csrptr = (const csr_matrix_t *) Atmp->repr;
        A.Resize(csrptr->n,csrptr->n);
        csr_matrix_expand_to_dense (A.Data(), 0, A.m(), csrptr);

        Int n = A.m();

        STEPS.Resize(n,n);
        SetValue(STEPS,0.0);

        RHS.Resize(n,nrhs);
        XTrue.Resize(n,nrhs);
        SetValue(XTrue,1.0);
        //      UniformRandom(XTrue);
        blas::Gemm('N','N',n,nrhs,n,1.0,&A(0,0),n,&XTrue(0,0),n,0.0,&RHS(0,0),n);

        //cal dposv
        double norm = 0;

        //do a solve
        DblNumMat X = RHS;
        //      DblNumMat A2 = A;
        //      lapack::Posv('L',n,nrhs,&A2(0,0),n,&X(0,0),n);
        //      
        //      blas::Axpy(n*nrhs,-1.0,&XTrue(0,0),1,&X(0,0),1);
        //      norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
        //      logfileptr->OFS()<<"Norm of residual after POSV is "<<norm<<std::endl;
        //
        //      X = RHS;
        lapack::Potrf('L',n,&A(0,0),n);

        for(Int i = 0; i<A.m();++i){
          for(Int j = 0; j<A.n();++j){
            if(j>i){
              A(i,j)=LIBCHOLESKY::ZERO<double>();
            }
          }
        }


        //      logfileptr->OFS()<<"Cholesky factor of POTRF:"<<A<<std::endl;


        //simulate a potrs in order to get the intermediate values
        for(Int j = 0; j<nrhs;++j){
          for(Int k = 0; k<A.m();++k){
            if(X(k,j)!=0){
              X(k,j) = X(k,j) / A(k,k);
              if(j==0){
                STEPS(k,k) = X(k,j);
              }
              for(Int i = k+1;i<A.m();++i){
                X(i,j) -= X(k,j)*A(i,k);
                if(j==0){
                  STEPS(i,k) -= X(k,j)*A(i,k);

                  for(Int kk = 0;kk<k;++kk){
                    //    STEPS(i,k)+=STEPS(i,kk);
                  }
                }
              }
            }
          }
        }

        //blas::Trsm('L','L','N','N',n,nrhs,1.0,&A(0,0),n,&X(0,0),n);
        logfileptr->OFS()<<"Solution after forward substitution:"<<X<<std::endl;
        blas::Trsm('L','L','T','N',n,nrhs,1.0,&A(0,0),n,&X(0,0),n);
        logfileptr->OFS()<<"Solution after back substitution:"<<X<<std::endl;

        blas::Axpy(n*nrhs,-1.0,&XTrue(0,0),1,&X(0,0),1);
        norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
        logfileptr->OFS()<<"Norm of residual after MYPOSV is "<<norm<<std::endl;

        destroy_sparse_matrix (Atmp);

      }

    }





    sparse_matrix_t* Atmp = load_sparse_matrix (informat, filename.c_str());
    sparse_matrix_convert (Atmp, CSC);
    const csc_matrix_t * cscptr = (const csc_matrix_t *) Atmp->repr;
    DistSparseMatrix<Real> HMat(cscptr,worldcomm);

    if(iam==0){ cout<<"Matrix order is "<<HMat.size<<endl; }


    destroy_sparse_matrix (Atmp);

    //SparseMatrixStructure Global = HMat.GetGlobalStructure();
    //logfileptr->OFS()<<Global.colptr<<std::endl;
    //logfileptr->OFS()<<Global.rowind<<std::endl;

    //ETree etree(Global);
    //etree.PostOrderTree();


    //logfileptr->OFS()<<"etree is "<<etree<<endl;

    //do the symbolic factorization and build supernodal matrix
    SupernodalMatrix<double> SMat(HMat,mapping,worldcomm);

    SMat.Factorize(worldcomm);

    RHS.Resize(SMat.Size(),nrhs);
    MPI_Bcast(RHS.Data(),RHS.ByteSize(),MPI_BYTE,0,worldcomm);

    logfileptr->OFS()<<"RHS:"<<RHS<<endl;

    DblNumMat X = RHS;

    NumMat<Real> fullMatrix;
    SMat.GetFullFactors(fullMatrix,worldcomm);

    if(iam==0){
      blas::Axpy(fullMatrix.m()*fullMatrix.n(),-1.0,&A(0,0),1,&fullMatrix(0,0),1);

      logfileptr->OFS()<<fullMatrix<<endl;
    }


    SMat.Solve(&X,worldcomm);



    XTrue.Resize(SMat.Size(),nrhs);
    MPI_Bcast(XTrue.Data(),XTrue.ByteSize(),MPI_BYTE,0,worldcomm);


      blas::Axpy(X.m()*X.n(),-1.0,&XTrue(0,0),1,&X(0,0),1);
      double norm = lapack::Lange('F',X.m(),X.n(),&X(0,0),X.m());
      logfileptr->OFS()<<"Norm of residual after SPCHOL is "<<norm<<std::endl;



    MPI_Barrier(worldcomm);
    MPI_Comm_free(&worldcomm);
    delete logfileptr;
    MPI_Finalize();
  return 0;
}


