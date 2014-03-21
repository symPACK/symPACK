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
#include  "Mapping.hpp"

#include  "blas.hpp"
#include  "lapack.hpp"
#include  "FBMatrix_upcxx.hpp"
#include  "NZBlock.hpp"
#include  "SuperNode.hpp"
#include  "SupernodalMatrix.hpp"

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

//Return the row structure in the permuted matrix
void getRowStruct(const SparseMatrixStructure & global, const ETree & etree, const Int iPORow, std::vector<Int> & rowStruct){
  for(Int iPOCurCol = 1; iPOCurCol<=iPORow;++iPOCurCol){
    Int iCurCol = etree.FromPostOrder(iPOCurCol);

    Int iFirstRowPtr = global.colptr(iCurCol-1);
    Int iLastRowPtr = global.colptr(iCurCol)-1;
    for(Int iCurRowPtr=iFirstRowPtr ;iCurRowPtr<=iLastRowPtr;++iCurRowPtr){
      Int iPOCurRow = etree.ToPostOrder(global.rowind(iCurRowPtr-1));
      if(iPOCurRow == iPORow){
        rowStruct.push_back(iPOCurCol);
      }

      if(iPOCurRow >= iPORow){
        break;
      }
    }
  }
}

void getLRowStruct(const SparseMatrixStructure & global, const ETree & etree, const Int iPORow, const std::vector<Int> & ARowStruct, std::set<Int> & LRowStruct){
  LRowStruct.clear();
  for(Int i = 0; i<ARowStruct.size();++i){
    Int iCurNode = ARowStruct[i];
    //tracing from iCurRow to iRow;
    logfileptr->OFS()<<"Starting from node "<<iCurNode<<std::endl;
    if(iCurNode==iPORow){
      logfileptr->OFS()<<"Adding node "<<iCurNode<<std::endl;
      LRowStruct.insert(iCurNode);
    }
    else{
      while(iCurNode != iPORow && etree.PostParent(iCurNode-1) != 0){
        logfileptr->OFS()<<"Adding node "<<iCurNode<<std::endl;
        LRowStruct.insert(iCurNode);
        iCurNode = etree.PostParent(iCurNode-1);
        logfileptr->OFS()<<"Parent is "<<iCurNode<<std::endl;
      }
    }
  } 
}

int main(int argc, char **argv) 
{
  upcxx::init(&argc,&argv);
//  MPI_Init(&argc,&argv);


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


    FBMatrix_upcxx * Afactptr = new FBMatrix_upcxx();
  

    //upcxx::global_ptr<FBMatrix> Afactptr = upcxx::Create<FBMatrix>();

    Afactptr->Initialize();
    Afactptr->pcol = np;
    Afactptr->np = np;
    Afactptr->prow = sqrt(np);
    Afactptr->blksize = 1;

    Real timeSta, timeEnd;

      DblNumMat RHS;
      DblNumMat XTrue;

    //Read the input matrix
    sparse_matrix_file_format_t informat;
    informat = sparse_matrix_file_format_string_to_enum (informatstr.c_str());
    sparse_matrix_t* Atmp = load_sparse_matrix (informat, filename.c_str());
    //sparse_matrix_expand_symmetric_storage (Atmp);
    sparse_matrix_convert (Atmp, CSC);

    const csc_matrix_t * cscptr = (const csc_matrix_t *) Atmp->repr;
    DistSparseMatrix<Real> HMat(cscptr);

    Afactptr->n = HMat.size;
    if(MYTHREAD==0){ cout<<"Matrix order is "<<Afactptr->n<<endl; }


    destroy_sparse_matrix (Atmp);

    //SparseMatrixStructure Global = HMat.GetGlobalStructure();
    //logfileptr->OFS()<<Global.colptr<<std::endl;
    //logfileptr->OFS()<<Global.rowind<<std::endl;

    //ETree etree(Global);
    //etree.PostOrderTree();


    //logfileptr->OFS()<<"etree is "<<etree<<endl;

    //do the symbolic factorization and build supernodal matrix
    MAPCLASS mapping(np, sqrt(np), np, 1);
    MPI_Comm worldcomm;
    MPI_Comm_dup(MPI_COMM_WORLD,&worldcomm);
    SupernodalMatrix<double> SMat(HMat,mapping,worldcomm);

    SMat.Factorize(worldcomm);

//
//
//
//    //Try to determine nth row structure in A
//    std::vector< std::vector<Int> > fullRowStruct(HMat.size);
//    Int sizerowstruct = 0;
//    ETree & etree = SMat.GetETree();
//    SparseMatrixStructure global = SMat.GetGlobalStructure();
//    for(Int i=0;i<SMat.GetLocalSupernodes().size();++i){
//      SuperNode & snode = SMat.GetLocalSupernode(i);
//      
//
////      logfileptr->OFS()<<"Supernode "<<snode.id<<" owns columns "<<snode.firstCol<<" to "<<snode.lastCol<<std::endl;
//      logfileptr->OFS()<<"Supernode "<<snode.Id()<<" owns blocks:"<<std::endl;
//
//
//      for(int blkidx=0;blkidx<snode.NZBlockCnt();++blkidx){
//        NZBlock<double> & nzblk = snode.GetNZBlock(blkidx);
//        Int lastRow = nzblk.GIndex() + nzblk.NRows() -1;
//
//
//        for(Int iRow = nzblk.GIndex(); iRow<=lastRow; ++iRow){
//          std::vector<Int> & rowStruct = fullRowStruct[iRow-1];
//          if(rowStruct.size()==0){
//            global.GetARowStruct(etree, iRow,  rowStruct);
//            sizerowstruct+=rowStruct.size()*sizeof(Int);
//
//            
//            logfileptr->OFS()<<"Row "<<iRow<<" structure is "<<rowStruct<<std::endl;
//            std::set<Int>  LRowStruct;
//            global.GetLRowStruct(etree, iRow,  rowStruct, LRowStruct);
//            logfileptr->OFS()<<"Row "<<iRow<<" of L structure is "<<LRowStruct.size()<<std::endl;
//            for (std::set<Int>::iterator it=LRowStruct.begin(); it!=LRowStruct.end(); ++it){
//              logfileptr->OFS() << ' ' << *it;
//            }
//            
//            logfileptr->OFS() << std::endl;
//          }
//        }
//
////        logfileptr->OFS()<<nzblk<<std::endl;
//        logfileptr->OFS()<<"L("<<nzblk.GIndex()<<".."<<lastRow<<" , "<<snode.FirstCol()<<".."<<snode.LastCol()<<")"<<std::endl;
//      }
//
//
//      logfileptr->OFS()<<"--------------------------------------------------"<<std::endl;
//    }
//
//
//      logfileptr->OFS()<<"total size of row structure is "<<sizerowstruct<<" bytes"<<std::endl;









    upcxx::barrier();
    MPI_Comm_free(&worldcomm);
    delete logfileptr;
    upcxx::finalize();
    return;



    //Alloc RHS
    //Call PDGEMM to compute XTrue

    upcxx::barrier();



    //TODO This needs to be fixed
    Afactptr->prow = 1;//sqrt(np);
    Afactptr->pcol = np;//Afact.prow;
    np = Afactptr->prow*Afactptr->pcol;
    if(MYTHREAD==0){
      cout<<"Number of cores to be used: "<<np<<endl;
//      cout<<"Matrix order is "<<Afact.n<<endl;
    }



    //Allocate chunks of the matrix on each processor
    Afactptr->Allocate(np,Afactptr->n,blksize);
    //TODO: recode this
    //Afactptr->Distribute(HMat);

    upcxx::barrier();
    upcxx::wait();
    //deallocate HMat
    //A.Clear();

    logfileptr->OFS()<<"distributed"<<endl;

      timeSta =  omp_get_wtime( );
    TIMER_START(FANBOTH);
    if(MYTHREAD==Afactptr->MAP(Afactptr->n-1,Afactptr->n-1)){
#ifdef _DEBUG_
      cout<<"LOCK IS TAKEN BY"<<MYTHREAD<<endl;
#endif
    }

#ifdef _DEBUG_
    logfileptr->OFS()<<"initializing"<<endl;
#endif

    upcxx::barrier();

    Afactptr->NumericalFactorization();
    Afactptr->WaitFactorization();
    timeEnd =  omp_get_wtime( );

    TIMER_STOP(FANBOTH);

#ifdef _DEBUG_
    logfileptr->OFS()<<"gathering"<<endl;
#endif

    //gather all data on P0
 //   DblNumMat_upcxx Afinished;
 //   upcxx::barrier();
 //   if(MYTHREAD==0)
 //   {
 //     cout<<"FAN-BOTH: "<<timeEnd-timeSta<<endl;
 //     Afactptr->Gather(Afinished);
 //   }

 //   upcxx::wait();




delete logfileptr;
    upcxx::finalize();
    return;

    DblNumMat A(Afactptr->n,Afactptr->n);
    //csc_matrix_expand_to_dense (A.Data(), 0, Afact.n, (const csc_matrix_t *) Atmp->repr);
    destroy_sparse_matrix (Atmp);



      Int n = Afactptr->n;
      Int nrhs = 5;

    if(MYTHREAD==0){

      RHS.Resize(n,nrhs);
      XTrue.Resize(n,nrhs);
      UniformRandom(XTrue);
      blas::Gemm('N','N',n,nrhs,n,1.0,&A(0,0),n,&XTrue(0,0),n,0.0,&RHS(0,0),n);
    }

    upcxx::barrier();

    //TODO This needs to be fixed
    Afactptr->prow = 1;//sqrt(np);
    Afactptr->pcol = np;//Afact.prow;
    np = Afactptr->prow*Afactptr->pcol;
    if(MYTHREAD==0){
      cout<<"Number of cores to be used: "<<np<<endl;
    }



    //Allocate chunks of the matrix on each processor
    Afactptr->Allocate(np,Afactptr->n,blksize);
    Afactptr->Distribute(A);

    upcxx::barrier();
    upcxx::wait();
    A.Clear();

    logfileptr->OFS()<<"distributed"<<endl;

      timeSta =  omp_get_wtime( );
    TIMER_START(FANBOTH);
    if(MYTHREAD==Afactptr->MAP(Afactptr->n-1,Afactptr->n-1)){
#ifdef _DEBUG_
      cout<<"LOCK IS TAKEN BY"<<MYTHREAD<<endl;
#endif
    }

#ifdef _DEBUG_
    logfileptr->OFS()<<"initializing"<<endl;
#endif

    upcxx::barrier();

    Afactptr->NumericalFactorization();
    Afactptr->WaitFactorization();
    timeEnd =  omp_get_wtime( );

    TIMER_STOP(FANBOTH);

#ifdef _DEBUG_
    logfileptr->OFS()<<"gathering"<<endl;
#endif

    //gather all data on P0
    DblNumMat_upcxx Afinished;
    upcxx::barrier();
    if(MYTHREAD==0)
    {
      cout<<"FAN-BOTH: "<<timeEnd-timeSta<<endl;
      Afactptr->Gather(Afinished);
    }

    upcxx::wait();


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
    delete Afactptr;

    upcxx::finalize();

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


