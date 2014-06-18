/// @file run_sparse_fanboth_upcxx.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2014-01-23
//#define _DEBUG_


#include <upcxx.h>

#include <time.h>
#ifndef __PGI
#include <random>
#endif
#include <omp.h>

#include  "Environment.hpp"
#include  "DistSparseMatrix.hpp"
#include  "NumVec.hpp"
#include  "NumMat.hpp"
#include  "utility.hpp"
#include  "ETree.hpp"
#include  "blas.hpp"
#include  "lapack.hpp"
#include  "FBMatrix_upcxx.hpp"
#include  "NZBlock.hpp"

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



    //Unit test of NZBlock class
    {
      Int iGIndex = 10;
      Int iLIndex = 42;
      Int iRows = 5;
      Int iCols = 2;
      Int iCnt = iRows*iCols;
      double pdInput[iCnt];
      for(Int i = 0; i<iCnt;++i){pdInput[i]=static_cast<double>(i);}

      logfileptr->OFS()<<"----NZBlock unit test----"<<std::endl;
      logfileptr->OFS()<<"creating the NZBlock with "<<iCnt<<" double"<<std::endl;
      NZBlock<double> test (iRows,iCols,iGIndex,iLIndex,pdInput);
      logfileptr->OFS()<<"data pointer is "<<test.Data()<<std::endl;
      logfileptr->OFS()<<"NRows is "<<test.NRows()<<std::endl;
      logfileptr->OFS()<<"NCols is "<<test.NCols()<<std::endl;
      logfileptr->OFS()<<"Nzcnt is "<<test.Nzcnt()<<std::endl;
      logfileptr->OFS()<<"Nzval pointer is "<<test.Nzval()<<std::endl;
      logfileptr->OFS()<<"global index is "<<test.GIndex()<<std::endl;
      logfileptr->OFS()<<"local index is "<<test.LIndex()<<std::endl;
      logfileptr->OFS()<<"parsing the nzvalues (1D index)"<<std::endl;
      for(Int i = 0; i< test.Nzcnt(); ++i){ logfileptr->OFS()<<test.Nzval(i)<<" "; }
      logfileptr->OFS()<<std::endl;

      logfileptr->OFS()<<"parsing the nzvalues (2D indices) "<<std::endl;
      for(Int i = 0; i< test.NRows(); ++i){
        for(Int j = 0; j< test.NCols(); ++j){
         logfileptr->OFS()<<test.Nzval(i,j)<<" ";
        }
        logfileptr->OFS()<<std::endl;
      }
      logfileptr->OFS()<<std::endl;
      logfileptr->OFS()<<std::endl;
      logfileptr->OFS()<<test<<std::endl;


      logfileptr->OFS()<<"----End of NZBlock unit test----"<<std::endl;
    }




    delete logfileptr;
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


