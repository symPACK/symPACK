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
#include  "ngchol/CommTypes.hpp"
#include  "ngchol/Ordering.hpp"

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


  RHS.Resize(n,nrhs);
  XTrue.Resize(n,nrhs);

  //SetValue(XTrue,1.0);
  Int val = 1.0;
  for(Int i = 0; i<n;++i){ 
    for(Int j=0;j<nrhs;++j){
      XTrue(i,j) = val;
      val = -val;
    }
  }
//        UniformRandom(XTrue);

  if(iam==0){
    cout<<"Starting spGEMM"<<endl;
  }

  timeSta = get_time();

#if 1
{
  SparseMatrixStructure Local = HMat.GetLocalStructure();
  SparseMatrixStructure Global;
  Local.ToGlobal(Global);
  Global.ExpandSymmetric();

  Int numColFirst = n / np;


  SetValue(RHS,0.0);
  for(Int j = 1; j<=n; ++j){
    Int iOwner = std::min((j-1)/numColFirst,np-1);
    if(iam == iOwner){
      Int iLocal = (j-(numColFirst)*iOwner);
      //do I own the column ?
      double t = XTrue(j-1,0);
      //do a dense mat mat mul ?
      for(Int ii = Local.colptr[iLocal-1]; ii< Local.colptr[iLocal];++ii){
        Int row = Local.rowind[ii-1];
        RHS(row-1,0) += t*HMat.nzvalLocal[ii-1];
        if(row>j){
          RHS(j-1,0) += XTrue(row-1,0)*HMat.nzvalLocal[ii-1];
        }
      }
    }
  }
  //Do a reduce of RHS
  MPI_Allreduce(MPI_IN_PLACE,&RHS(0,0),RHS.Size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}
#else








#if 0
  sp_dgemm_dist("N","N", n, XTrue.n(), n, 
      LIBCHOLESKY::ONE<MYSCALAR>(), HMat, XTrue.Data(), XTrue.m(), 
      LIBCHOLESKY::ZERO<MYSCALAR>(), RHS.Data(), RHS.m());
#else
  for(Int i = 0; i<n;++i){ 
    for(Int j=0;j<nrhs;++j){
      RHS(i,j) = i+1;
    }
  }
#endif

#endif

  timeEnd = get_time();
  if(iam==0){
    cout<<"spGEMM time: "<<timeEnd-timeSta<<endl;
  }


  //   logfileptr->OFS()<<RHS<<endl;
#endif


  if(iam==0){
    cout<<"Starting allocation"<<endl;
  }
  timeSta = get_time();
  //do the symbolic factorization and build supernodal matrix
  SupernodalMatrix<double> SMat(HMat,maxSnode,*mapping,maxIsend,maxIrecv,worldcomm);

  timeEnd = get_time();
  if(iam==0){
    cout<<"Allocation time: "<<timeEnd-timeSta<<endl;
  }


#ifdef _CHECK_RESULT_
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

      RHS2 = RHS;//.Resize(n,nrhs);
      XTrue2 = XTrue;//.Resize(n,nrhs);
      //SetValue(XTrue2,1.0);
      //      UniformRandom(XTrue);
      blas::Gemm('N','N',n,nrhs,n,1.0,&A(0,0),n,&XTrue2(0,0),n,0.0,&RHS2(0,0),n);

      double norm = 0;
      DblNumMat tmp = RHS2;
      blas::Axpy(n*nrhs,-1.0,&RHS(0,0),1,&tmp(0,0),1);
      norm = lapack::Lange('F',n,nrhs,&tmp(0,0),n);
      cout<<"Norm between RHS is "<<norm<<std::endl;


      //Order the matrix
      Ordering order = SMat.GetOrdering();

      DblNumMat Aperm(A.m(),A.n());
      for(Int i = 0; i<A.m(); ++i){
        for(Int j = 0; j<A.n(); ++j){
          Aperm(i,j) = A(order.perm[i]-1,order.perm[j]-1);
        }
      }


      //cal dposv
      norm = 0;

      //do a solve
      DblNumMat X(RHS2.m(),RHS2.n());
      DblNumMat XTrue3(XTrue2.m(),XTrue2.n());
      for(Int i = 0; i<n;++i){ 
        for(Int j=0;j<nrhs;++j){
          X(i,j) = RHS2(order.perm[i]-1,j);
          XTrue3(i,j) = XTrue2(order.perm[i]-1,j);
        }
      }



//logfileptr->OFS()<<"RHS2 "<<X<<endl;
//logfileptr->OFS()<<"Aperm "<<Aperm<<endl;



{
      DblNumMat B = Aperm;
      DblNumMat XTrue4 = XTrue3;
      DblNumMat X3 = X;

      lapack::Potrf('L',n,&B(0,0),n);
      lapack::Potrs('L',n,nrhs,&B(0,0),n,&X3(0,0),n);

  blas::Axpy(X3.m()*X3.n(),-1.0,&XTrue4(0,0),1,&X3(0,0),1);
  double norm2 = lapack::Lange('F',X3.m(),X3.n(),&X3(0,0),X3.m());
    cout<<"Norm of residual after potrf (potrs) is "<<norm2<<std::endl;

}






      lapack::Potrf('L',n,&Aperm(0,0),n);


      for(Int i = 0; i<A.m();++i){
        for(Int j = 0; j<A.n();++j){
          if(j>i){
            Aperm(i,j)=LIBCHOLESKY::ZERO<double>();
          }
        }
      }
//logfileptr->OFS()<<"Lperm "<<Aperm<<endl;

      //simulate a potrs in order to get the intermediate values
      for(Int j = 0; j<nrhs;++j){
        for(Int k = 0; k<A.m();++k){
          if(X(k,j)!=0){
            X(k,j) = X(k,j) / Aperm(k,k);
            for(Int i = k+1;i<A.m();++i){
              X(i,j) -= X(k,j)*Aperm(i,k);
            }
          }
        }
      }

      //blas::Trsm('L','L','N','N',n,nrhs,1.0,&A(0,0),n,&X(0,0),n);
      fwdSol.Resize(X.m(),X.n());
      for(Int i = 0; i<n;++i){ 
        for(Int j=0;j<nrhs;++j){
          fwdSol(i,j) = X(order.invp[i]-1,j);
        }
      }


#ifdef _DEBUG_
      logfileptr->OFS()<<"Solution after forward substitution:"<<X<<std::endl;
#endif
      blas::Trsm('L','L','T','N',n,nrhs,1.0,&Aperm(0,0),n,&X(0,0),n);
#ifdef _DEBUG_
      logfileptr->OFS()<<"Solution after back substitution:"<<X<<std::endl;
#endif

      blas::Axpy(n*nrhs,-1.0,&XTrue3(0,0),1,&X(0,0),1);
      norm = lapack::Lange('F',n,nrhs,&X(0,0),n);
      cout<<"Norm of residual after MYPOSV is "<<norm<<std::endl;

      destroy_sparse_matrix (Atmp);

    }
    //
  }
  A.Clear();
#endif


#ifdef _CHECK_RESULT_SEQ_
//  if(iam==0){
//    DblNumMat RHSDiff = RHS;
//
//    blas::Axpy(RHS.m()*RHS.n(),-1.0,&RHS2(0,0),1,&RHSDiff(0,0),1);
//    double norm = lapack::Lange('F',RHS.m(),RHS.n(),&RHSDiff(0,0),RHS.m());
//    logfileptr->OFS()<<"Norm of residual between RHS is "<<norm<<std::endl;
//
//    if(abs(norm)>=1e-5){
//      for(Int i = 0;i<RHS.m();++i){
//        logfileptr->OFS()<<i+1<<"   "<<RHS(i,0)-RHS2(i,0)<<endl;
//      }
//      abort();
//    }
//
//  }
#endif




#endif

















  if(iam==0){
    cout<<"Starting Factorization"<<endl;
  }


  timeSta = get_time();
  TIMER_START(SPARSE_FAN_OUT);
  SMat.FanOut();
  TIMER_STOP(SPARSE_FAN_OUT);
  timeEnd = get_time();

  if(iam==0){
    cout<<"Factorization time: "<<timeEnd-timeSta<<endl;
  }
    logfileptr->OFS()<<"Factorization time: "<<timeEnd-timeSta<<endl;

#ifdef _CHECK_RESULT_
#ifdef _CHECK_RESULT_POTRS_SEQ_
      NumMat<Real> fullMatrix;
      SMat.GetFullFactors(fullMatrix);
  if(iam==0){
//      logfileptr->OFS()<<"FullRes "<<fullMatrix<<std::endl;
      Ordering order = SMat.GetOrdering();
      DblNumMat X(RHS.m(),RHS.n());
      DblNumMat XTrue3(XTrue.m(),XTrue.n());
      for(Int i = 0; i<n;++i){ 
        for(Int j=0;j<nrhs;++j){
          X(i,j) = RHS(order.perm[i]-1,j);
          XTrue3(i,j) = XTrue(order.perm[i]-1,j);
        }
      }


      lapack::Potrs('L',n,nrhs,&fullMatrix(0,0),n,&X(0,0),n);

  blas::Axpy(X.m()*X.n(),-1.0,&XTrue3(0,0),1,&X(0,0),1);
  double norm2 = lapack::Lange('F',X.m(),X.n(),&X(0,0),X.m());
    cout<<"Norm of residual after SPCHOL (potrs) is "<<norm2<<std::endl;

  }
#endif



#ifdef _CHECK_RESULT_SEQ_
  RHS2.Resize(SMat.Size(),nrhs);
  MPI_Bcast(RHS2.Data(),RHS2.ByteSize(),MPI_BYTE,0,worldcomm);

  XTrue2.Resize(SMat.Size(),nrhs);
  MPI_Bcast(XTrue2.Data(),XTrue2.ByteSize(),MPI_BYTE,0,worldcomm);
  //    logfileptr->OFS()<<"RHS:"<<RHS<<endl;
#endif

  //sort X the same way (permute rows)
#ifdef _CHECK_RESULT_SEQ_
  DblNumMat X = RHS2;
#else
  DblNumMat X = RHS;
#endif

  if(iam==0){
    cout<<"Starting solve"<<endl;
  }

  timeSta = get_time();
#ifdef _CHECK_RESULT_SEQ_
  fwdSol.Resize(SMat.Size(),nrhs);
  MPI_Bcast(fwdSol.Data(),fwdSol.ByteSize(),MPI_BYTE,0,worldcomm);

  DblNumMat poFwdSol = fwdSol;
  SMat.Solve(&X,poFwdSol);
#else
  SMat.Solve(&X);
#endif
  timeEnd = get_time();

  if(iam==0){
    cout<<"Solve time: "<<timeEnd-timeSta<<endl;
  }


  SMat.GetSolution(X);


  if(iam==0){
//  blas::Axpy(X.m()*X.n(),-1.0,&XTrue(0,0),1,&X(0,0),1);
//  double norm = lapack::Lange('F',X.m(),X.n(),&X(0,0),X.m());
//    cout<<"Norm of residual after SPCHOL is "<<norm/normB<<std::endl;
}

{
  SparseMatrixStructure Local = HMat.GetLocalStructure();

  Int numColFirst = n / np;

  DblNumMat AX(n,nrhs);
  SetValue(AX,0.0);
  
  for(Int j = 1; j<=n; ++j){
    Int iOwner = std::min((j-1)/numColFirst,np-1);
    if(iam == iOwner){
      Int iLocal = (j-(numColFirst)*iOwner);
      //do I own the column ?
      double t = X(j-1,0);
      //do a dense mat mat mul ?
      for(Int ii = Local.colptr[iLocal-1]; ii< Local.colptr[iLocal];++ii){
        Int row = Local.rowind[ii-1];
        AX(row-1,0) += t*HMat.nzvalLocal[ii-1];
        if(row>j){
          AX(j-1,0) += X(row-1,0)*HMat.nzvalLocal[ii-1];
        }
      }
    }
  }
  //Do a reduce of RHS
  MPI_Allreduce(MPI_IN_PLACE,&AX(0,0),AX.Size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  if(iam==0){
  blas::Axpy(AX.m()*AX.n(),-1.0,&RHS(0,0),1,&AX(0,0),1);
  double normAX = lapack::Lange('F',AX.m(),AX.n(),&AX(0,0),AX.m());
  double normRHS = lapack::Lange('F',RHS.m(),RHS.n(),&RHS(0,0),RHS.m());
    cout<<"Norm of residual after SPCHOL is "<<normAX/normRHS<<std::endl;
  }
}





#endif

  delete mapping;


  MPI_Barrier(worldcomm);
  MPI_Comm_free(&worldcomm);
  delete logfileptr;
  MPI_Finalize();
  return 0;
}


