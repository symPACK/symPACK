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


    MPI_Barrier(worldcomm);
    MPI_Comm_free(&worldcomm);
    delete logfileptr;
    MPI_Finalize();
  return 0;
}


