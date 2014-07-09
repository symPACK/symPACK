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

#include  "ngchol/sp_blas.hpp"

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


#define MYSCALAR double


using namespace LIBCHOLESKY;


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
  HMat.CopyData(cscptr->n,cscptr->nnz,cscptr->colptr,cscptr->rowidx,(MYSCALAR *)cscptr->values);


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


  //do the symbolic factorization and build supernodal matrix
  SupernodalMatrix<double> SMat(HMat,maxSnode,*mapping,maxIsend,maxIrecv,worldcomm);

  timeSta = get_time();
  TIMER_START(SPARSE_FAN_OUT2);
  SMat.FanOut2();
  TIMER_STOP(SPARSE_FAN_OUT2);
  timeEnd = get_time();

  if(iam==0){
    cout<<"Factorization time: "<<timeEnd-timeSta<<endl;
  }
    logfileptr->OFS()<<"Factorization time: "<<timeEnd-timeSta<<endl;



#ifdef _DEBUG_
  //    NumMat<Real> fullMatrix;
  //    SMat.GetFullFactors(fullMatrix);
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
  const IntNumVec & perm = SMat.GetColPerm(); 
  for(Int row = 1; row<= RHS.m(); ++row){
    for(Int col = 1; col<= RHS.n(); ++col){
#ifdef _CHECK_RESULT_SEQ_
      X(row-1,col-1) = RHS2(perm[row-1]-1,col-1);
#else
      X(row-1,col-1) = RHS(perm[row-1]-1,col-1);
#endif
    }
  }

{
  DblNumMat X3 = X;
  SupernodalMatrix<double> SMat2(HMat,maxSnode,*mapping,maxIsend,maxIrecv,worldcomm);

  timeSta = get_time();
  TIMER_START(SPARSE_FAN_OUT);
  SMat2.FanOut();
  TIMER_STOP(SPARSE_FAN_OUT);
  timeEnd = get_time();

  if(iam==0){
    cout<<"Factorization time: "<<timeEnd-timeSta<<endl;
  }
    logfileptr->OFS()<<"Factorization time: "<<timeEnd-timeSta<<endl;

  SMat.Solve(&X3);
  SMat.GetSolution(X3);



  //Sort back X
  DblNumMat X4(X.m(),X.n());
  for(Int row = 1; row<= X.m(); ++row){
    for(Int col = 1; col<= X.n(); ++col){
      X4(perm[row-1]-1,col-1) = X3(row-1,col-1);
    }
  }


//      logfileptr->OFS()<<X2<<endl;

  blas::Axpy(X4.m()*X4.n(),-1.0,&XTrue(0,0),1,&X4(0,0),1);
  double norm = lapack::Lange('F',X4.m(),X4.n(),&X4(0,0),X4.m());
  logfileptr->OFS()<<"Norm of residual after SPCHOL is "<<norm<<std::endl;


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
  SMat.Solve(&X,poFwdSol);
#else
  SMat.Solve(&X);
#endif
  SMat.GetSolution(X);



  //Sort back X
  DblNumMat X2(X.m(),X.n());
  for(Int row = 1; row<= X.m(); ++row){
    for(Int col = 1; col<= X.n(); ++col){
      X2(perm[row-1]-1,col-1) = X(row-1,col-1);
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


