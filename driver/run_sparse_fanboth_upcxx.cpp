/// @file run_sparse_fanboth.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2014-01-23
//#define _DEBUG_


#include <mpi.h>

#include <time.h>
#include <random>
#include <omp.h>

#include  "ngchol.hpp"
//#include  "Environment.hpp"
//#include  "DistSparseMatrix.hpp"
//#include  "NumVec.hpp"
//#include  "NumMat.hpp"
//#include  "utility.hpp"
//#include  "ETree.hpp"
//#include  "Mapping.hpp"
//
//#include  "blas.hpp"
#include  "ngchol/sp_blas.hpp"
//#include  "lapack.hpp"
////#include  "NZBlock.hpp"
//#include  "SuperNode.hpp"
//#include  "SupernodalMatrix.hpp"
//
////#include <async.h>
//#include  "LogFile.hpp"


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




//#ifdef USE_TAU
//#include "TAU.h"
//#elif defined (PROFILE) || defined(PMPI)
//#define TAU
//#include "timer.hpp"
//#endif
//
//
//
//#ifdef USE_TAU 
//#define TIMER_START(a) TAU_START(TOSTRING(a));
//#define TIMER_STOP(a) TAU_STOP(TOSTRING(a));
//#elif defined (PROFILE)
//#define TIMER_START(a) TAU_FSTART(a);
//#define TIMER_STOP(a) TAU_FSTOP(a);
//#else
//#define TIMER_START(a)
//#define TIMER_STOP(a)
//#endif

using namespace LIBCHOLESKY;
using namespace std;



#define MYSCALAR double





int main(int argc, char **argv) 
{
  MPI_Init(&argc,&argv);

  MPI_Comm worldcomm;
  MPI_Comm_dup(MPI_COMM_WORLD,&worldcomm);
  //  int np, iam;
  MPI_Comm_size(worldcomm,&np);
  MPI_Comm_rank(worldcomm,&iam);

  MAPCLASS * mapping;
  if( np % (Int)sqrt(np) == 0){
    mapping = new MAPCLASS(np, sqrt(np), np, 1);
  }
  else{
    mapping = new MAPCLASS(np, np, np, 1);
  }


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
  Int maxSnode = -1;
  if( options.find("-b") != options.end() ){
    maxSnode= atoi(options["-b"].c_str());
  }

  Int maxIsend = 0;
  if( options.find("-is") != options.end() ){
    maxIsend= atoi(options["-is"].c_str());
  }

  Int maxIrecv = 0;
  if( options.find("-ir") != options.end() ){
    maxIrecv= atoi(options["-ir"].c_str());
  }



  Real timeSta, timeEnd;

  sparse_matrix_file_format_t informat;
  TIMER_START(READING_MATRIX);
DistSparseMatrix<Real> HMat(worldcomm);
  //Read the input matrix
  if(informatstr == "CSC"){
       ParaReadDistSparseMatrix( filename.c_str(), HMat, worldcomm ); 
  }
  else{

  informat = sparse_matrix_file_format_string_to_enum (informatstr.c_str());



  sparse_matrix_t* Atmp = load_sparse_matrix (informat, filename.c_str());
  sparse_matrix_convert (Atmp, CSC);
  const csc_matrix_t * cscptr = (const csc_matrix_t *) Atmp->repr;
  HMat.CopyData(cscptr->n,cscptr->nnz,cscptr->colptr,cscptr->rowidx,cscptr->values);


  destroy_sparse_matrix (Atmp);
  }
  TIMER_STOP(READING_MATRIX);
  if(iam==0){ cout<<"Matrix order is "<<HMat.size<<endl; }

#ifdef _CHECK_RESULT_

  Int nrhs = 1;
  DblNumMat RHS;
  DblNumMat XTrue;

  Int n = HMat.size;



#ifdef _CHECK_RESULT_SEQ_
  DblNumMat fwdSol;
  DblNumMat RHS2;
  DblNumMat XTrue2;
  DblNumMat A;
  {
    if(iam==0){
      sparse_matrix_t* Atmp = load_sparse_matrix (informat, filename.c_str());
      sparse_matrix_expand_symmetric_storage (Atmp);
      sparse_matrix_convert (Atmp, CSR);
      const csr_matrix_t * csrptr = (const csr_matrix_t *) Atmp->repr;
      A.Resize(csrptr->n,csrptr->n);
      csr_matrix_expand_to_dense (A.Data(), 0, A.m(), csrptr);

      RHS2.Resize(n,nrhs);
      XTrue2.Resize(n,nrhs);
      SetValue(XTrue2,1.0);
      //      UniformRandom(XTrue);
      blas::Gemm('N','N',n,nrhs,n,1.0,&A(0,0),n,&XTrue2(0,0),n,0.0,&RHS2(0,0),n);


      //cal dposv
      double norm = 0;

      //do a solve
      DblNumMat X = RHS2;
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
            for(Int i = k+1;i<A.m();++i){
              X(i,j) -= X(k,j)*A(i,k);
            }
          }
        }
      }

      //blas::Trsm('L','L','N','N',n,nrhs,1.0,&A(0,0),n,&X(0,0),n);
      fwdSol = X;
#ifdef _DEBUG_
      logfileptr->OFS()<<"Solution after forward substitution:"<<X<<std::endl;
#endif
      blas::Trsm('L','L','T','N',n,nrhs,1.0,&A(0,0),n,&X(0,0),n);
#ifdef _DEBUG_
      logfileptr->OFS()<<"Solution after back substitution:"<<X<<std::endl;
#endif

      blas::Axpy(n*nrhs,-1.0,&XTrue2(0,0),1,&X(0,0),1);
      norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
      logfileptr->OFS()<<"Norm of residual after MYPOSV is "<<norm<<std::endl;

      destroy_sparse_matrix (Atmp);

    }
    //
  }
  A.Clear();
#endif



  RHS.Resize(n,nrhs);
  XTrue.Resize(n,nrhs);
  SetValue(XTrue,1.0);
  //      UniformRandom(XTrue);

  sp_dgemm_dist("N","N", n, XTrue.n(), n, 
      LIBCHOLESKY::ONE<MYSCALAR>(), HMat, XTrue.Data(), XTrue.m(), 
      LIBCHOLESKY::ZERO<MYSCALAR>(), RHS.Data(), RHS.m());


#ifdef _CHECK_RESULT_SEQ_
  if(iam==0){
    DblNumMat RHSDiff = RHS;

    blas::Axpy(RHS.m()*RHS.n(),-1.0,&RHS2(0,0),1,&RHSDiff(0,0),1);
    double norm = lapack::Lange('F',RHS.m(),RHS.n(),&RHSDiff(0,0),RHS.m());
    logfileptr->OFS()<<"Norm of residual between RHS is "<<norm<<std::endl;

    if(abs(norm)>=1e-5){
      for(Int i = 0;i<RHS.m();++i){
        logfileptr->OFS()<<i+1<<"   "<<RHS(i,0)-RHS2(i,0)<<endl;
      }
      abort();
    }

  }
#endif




  //   logfileptr->OFS()<<RHS<<endl;
#endif



  //SparseMatrixStructure Global = HMat.GetGlobalStructure();
  //logfileptr->OFS()<<Global.colptr<<std::endl;
  //logfileptr->OFS()<<Global.rowind<<std::endl;

  //ETree etree(Global);
  //etree.PostOrderTree();


  //logfileptr->OFS()<<"etree is "<<etree<<endl;

  //do the symbolic factorization and build supernodal matrix
  SupernodalMatrix<double> SMat(HMat,maxSnode,*mapping,maxIsend,maxIrecv,worldcomm);

  timeSta = get_time();
  TIMER_START(SPARSE_FAN_OUT);
//  SMat.Factorize(worldcomm);
  SMat.FanOut(worldcomm);
//  SMat.FanBoth(worldcomm);
  TIMER_STOP(SPARSE_FAN_OUT);
  timeEnd = get_time();

  if(iam==0){
    cout<<"Factorization time: "<<timeEnd-timeSta<<endl;
  }

#ifdef _DEBUG_
  //    NumMat<Real> fullMatrix;
  //    SMat.GetFullFactors(fullMatrix,worldcomm);
  //    logfileptr->OFS()<<fullMatrix<<std::endl;
#endif

#ifdef _CHECK_RESULT_


#ifdef _CHECK_RESULT_SEQ_
  RHS2.Resize(SMat.Size(),nrhs);
  MPI_Bcast(RHS2.Data(),RHS2.ByteSize(),MPI_BYTE,0,worldcomm);

  XTrue2.Resize(SMat.Size(),nrhs);
  MPI_Bcast(XTrue2.Data(),XTrue2.ByteSize(),MPI_BYTE,0,worldcomm);
  //    logfileptr->OFS()<<"RHS:"<<RHS<<endl;
#endif

  //sort X the same way (permute rows)
  DblNumMat X(RHS.m(),RHS.n());
  for(Int row = 1; row<= RHS.m(); ++row){
    for(Int col = 1; col<= RHS.n(); ++col){
#ifdef _CHECK_RESULT_SEQ_
      X(row-1,col-1) = RHS2(SMat.perm_[row-1]-1,col-1);
#else
//      X(row-1,col-1) = RHS(tree.FromPostOrder(row)-1,col-1);
      X(row-1,col-1) = RHS(SMat.perm_[row-1]-1,col-1);
#endif
    }
  }

//      logfileptr->OFS()<<RHS<<endl;
//      logfileptr->OFS()<<X<<endl;

  //    if(iam==0){
  //      logfileptr->OFS()<<fullMatrix<<endl;
  //
  //
  //
  //      blas::Axpy(fullMatrix.m()*fullMatrix.n(),-1.0,&A(0,0),1,&fullMatrix(0,0),1);
  //      double norm = lapack::Lange('F',fullMatrix.m(),fullMatrix.n(),&fullMatrix(0,0),fullMatrix.m());
  //      for(int j = 0; j<fullMatrix.n();++j){
  //        int maxi = 0;
  //        double maxelem = 0.0;
  //        for(int i = 0; i<fullMatrix.m();++i){
  //          if(std::abs(maxelem)<= std::abs(fullMatrix(i,j))){
  //            maxelem = fullMatrix(i,j);
  //            maxi = i;
  //          }
  //        }
  //        logfileptr->OFS()<<"Max of col "<<j<<" is "<<maxelem<<" at line "<<maxi<<std::endl;
  //      }
  //
  //      logfileptr->OFS()<<"Norm of residual between full matrices is "<<norm<<std::endl;
  //    }


#ifdef _CHECK_RESULT_SEQ_
  fwdSol.Resize(SMat.Size(),nrhs);
  MPI_Bcast(fwdSol.Data(),fwdSol.ByteSize(),MPI_BYTE,0,worldcomm);

  DblNumMat poFwdSol(fwdSol.m(),fwdSol.n());
  for(Int row = 1; row<= X.m(); ++row){
    for(Int col = 1; col<= X.n(); ++col){
      poFwdSol(row-1,col-1) = fwdSol(tree.FromPostOrder(row)-1,col-1);
    }
  }
  SMat.Solve(&X,worldcomm,poFwdSol);
#else
  SMat.Solve(&X,worldcomm);
#endif
  SMat.GetSolution(X,worldcomm);



  //Sort back X
  DblNumMat X2(X.m(),X.n());
  for(Int row = 1; row<= X.m(); ++row){
    for(Int col = 1; col<= X.n(); ++col){
//      X2(row-1,col-1) = X(tree.ToPostOrder(row)-1,col-1);
//      X2(tree.FromPostOrder(row)-1,col-1) = X(row-1,col-1);
      X2(SMat.perm_[row-1]-1,col-1) = X(row-1,col-1);
    }
  }


//      logfileptr->OFS()<<X2<<endl;

#ifdef _CHECK_RESULT_SEQ_
  blas::Axpy(X2.m()*X2.n(),-1.0,&XTrue2(0,0),1,&X2(0,0),1);
#else
  blas::Axpy(X2.m()*X2.n(),-1.0,&XTrue(0,0),1,&X2(0,0),1);
#endif
  double norm = lapack::Lange('F',X2.m(),X2.n(),&X2(0,0),X2.m());
  logfileptr->OFS()<<"Norm of residual after SPCHOL is "<<norm<<std::endl;
#endif

  delete mapping;


  MPI_Barrier(worldcomm);
  MPI_Comm_free(&worldcomm);
  delete logfileptr;
  MPI_Finalize();
  return 0;
}


