/// @file run_cholesky.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2013-08-31
#include  <mpi.h>

#include  "Environment.hpp"
#include  "DistSparseMatrix.hpp"
#include  "NumVec.hpp"
#include  "NumMat.hpp"
#include  "utility.hpp"
#include  "ETree.hpp"


#ifdef USE_TAU
#include "TAU.h"
#elif defined (PROFILE) || defined(PMPI)
#define TAU
#include "timer.hpp"
#endif



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

#if defined(PROFILE) || defined(USE_TAU) || defined(PMPI)
  TAU_PROFILE_INIT(argc, argv);
#endif

  MPI_Init( &argc, &argv );
  int mpirank, mpisize;
  MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpisize );


  try{
    MPI_Comm world_comm;

    // *********************************************************************
    // Input parameter
    // *********************************************************************
    std::map<std::string,std::string> options;

    OptionsCreate(argc, argv, options);

    // Default processor number
    Int nprow = 1;
    Int npcol = mpisize;

    if( options.find("-r") != options.end() ){
      if( options.find("-c") != options.end() ){
        nprow= atoi(options["-r"].c_str());
        npcol= atoi(options["-c"].c_str());
        if(nprow*npcol > mpisize){
          throw std::runtime_error("The number of used processors cannot be higher than the total number of available processors." );
        } 
      }
      else{
        throw std::runtime_error( "When using -r option, -c also needs to be provided." );
      }
    }
    else if( options.find("-c") != options.end() ){
      if( options.find("-r") != options.end() ){
        nprow= atoi(options["-r"].c_str());
        npcol= atoi(options["-c"].c_str());
        if(nprow*npcol > mpisize){
          throw std::runtime_error("The number of used processors cannot be higher than the total number of available processors." );
        } 
      }
      else{
        throw std::runtime_error( "When using -c option, -r also needs to be provided." );
      }
    }

    //Create a communicator with npcol*nprow processors
    MPI_Comm_split(MPI_COMM_WORLD, mpirank<nprow*npcol, mpirank, &world_comm);

cout<<"ALIVE"<<endl;

    if (mpirank<nprow*npcol){

      MPI_Comm_rank(world_comm, &mpirank );
      MPI_Comm_size(world_comm, &mpisize );


#if defined (PROFILE) || defined(PMPI) 
		TAU_PROFILE_SET_CONTEXT(world_comm);
#endif

    logfileptr = new LogFile(mpirank);

      //if( mpisize != nprow * npcol || nprow != npcol ){
      //  throw std::runtime_error( "nprow == npcol is assumed in this test routine." );
      //}

      if( mpirank == 0 )
        cout << "nprow = " << nprow << ", npcol = " << npcol << endl;

      std::string Hfile, Sfile;
      int isCSC = true;
      if( options.find("-T") != options.end() ){ 
        isCSC= ! atoi(options["-T"].c_str());
      }

//
//      int checkAccuracy = false;
//      if( options.find("-E") != options.end() ){ 
//        checkAccuracy= atoi(options["-E"].c_str());
//      }
//
//      int doFacto = true;
//      if( options.find("-F") != options.end() ){ 
//        doFacto= atoi(options["-F"].c_str());
//      }
//
//      int doSelInv = 1;
//      if( options.find("-Sinv") != options.end() ){ 
//        doSelInv= atoi(options["-Sinv"].c_str());
//      }
//
//      int doSymbfact = true;
//      if( options.find("-Symb") != options.end() ){ 
//        doSymbfact= atoi(options["-Symb"].c_str());
//      }
//
//      int doToDist = true;
//      if( options.find("-ToDist") != options.end() ){ 
//        doToDist= atoi(options["-ToDist"].c_str());
//      }
//
//
//      doFacto = doFacto && doSymbfact;
//
//
      if( options.find("-H") != options.end() ){ 
        Hfile = options["-H"];
      }
      else{
        throw std::logic_error("Hfile must be provided.");
      }

//
//      if( options.find("-S") != options.end() ){ 
//        Sfile = options["-S"];
//      }
//      else{
//        logfileptr->OFS() << "-S option is not given. " 
//          << "Treat the overlap matrix as an identity matrix." 
//          << std::endl << std::endl;
//      }
//
//      Int maxPipelineDepth = -1;
//      if( options.find("-P") != options.end() ){ 
//        maxPipelineDepth = atoi(options["-P"].c_str());
//      }
//      else{
//        logfileptr->OFS() << "-P option is not given. " 
//          << "Do not limit SelInv pipelining depth." 
//          << std::endl << std::endl;
//      }
//
//
//      Int numProcSymbFact;
//      if( options.find("-npsymbfact") != options.end() ){ 
//        numProcSymbFact = atoi( options["-npsymbfact"].c_str() );
//      }
//      else{
//        logfileptr->OFS() << "-npsymbfact option is not given. " 
//          << "Use default value (maximum number of procs)." 
//          << std::endl << std::endl;
//        numProcSymbFact = 0;
//      }
//
//
//
//
//
//
//
//      Int doSinv_Original = 0;
//      if( options.find("-SinvOriginal") != options.end() ){ 
//        doSinv_Original = atoi(options["-SinvOriginal"].c_str());
//      }
//      Int doSinv_Bcast = 1;
//      if( options.find("-SinvBcast") != options.end() ){ 
//        doSinv_Bcast = atoi(options["-SinvBcast"].c_str());
//      }
//
//      Int doSinvPipeline = 0;
//      if( options.find("-SinvPipeline") != options.end() ){ 
//        doSinvPipeline = atoi(options["-SinvPipeline"].c_str());
//      }
//
//
//
//      Int doConstructPattern = 1;
//      if( options.find("-Pattern") != options.end() ){ 
//        doConstructPattern = atoi(options["-Pattern"].c_str());
//      }
//
//      Int doPreSelinv = 1;
//      if( options.find("-PreSelinv") != options.end() ){ 
//        doPreSelinv = atoi(options["-PreSelinv"].c_str());
//      }
//
//      Int doSelinv = 1;
//      if( options.find("-Selinv") != options.end() ){ 
//        doSelinv = atoi(options["-Selinv"].c_str());
//      }
//
//
//      Real cshift = 0;
//      if( options.find("-Shift") != options.end() ){ 
//        cshift = atof(options["-Shift"].c_str());
//      }
//
//
//      std::string ColPerm;
//      if( options.find("-colperm") != options.end() ){ 
//        ColPerm = options["-colperm"];
//      }
//      else{
//        logfileptr->OFS() << "-colperm option is not given. " 
//          << "Use MMD_AT_PLUS_A." 
//          << std::endl << std::endl;
//        ColPerm = "MMD_AT_PLUS_A";
//      }
//
      // *********************************************************************
      // Read input matrix
      // *********************************************************************


      int      m, n;
      DistSparseMatrix<Real> HMat;
      Real timeSta, timeEnd;
      GetTime( timeSta );
      if(isCSC)
        ParaReadDistSparseMatrix( Hfile.c_str(), HMat, world_comm ); 
      else{
        ReadDistSparseMatrixFormatted( Hfile.c_str(), HMat, world_comm ); 
        HMat.comm = world_comm;
        ParaWriteDistSparseMatrix( "H.csc", HMat, world_comm ); 
      }

      GetTime( timeEnd );
      if( mpirank == 0 ){
        cout << "Time for reading H is " << timeEnd - timeSta << endl;
        cout << "H.size = " << HMat.size << endl;
        cout << "H.nnz  = " << HMat.nnz  << endl;
      }


     ETree elimTree;
     HMat.ToGlobalStruct();

     HMat.ConstructETree(elimTree);
//     logfileptr->OFS()<<elimTree<<std::endl;
     elimTree.PostOrderTree();
     logfileptr->OFS()<<elimTree<<std::endl;

     IntNumVec cc,rc;
     HMat.GetLColRowCount(elimTree,cc,rc);
     logfileptr->OFS()<<"Col count "<<cc<<std::endl;
     logfileptr->OFS()<<"Row count "<<rc<<std::endl;



    IntNumVec super; 
    HMat.FindSupernodes(elimTree,cc,super);
    logfileptr->OFS()<<"supernode indexes "<<super<<std::endl;

    HMat.SymbolicFactorization(elimTree,cc,super);
delete logfileptr;
    }
  }
  catch( std::exception& e )
  {
    std::cerr << "Processor " << mpirank << " caught exception with message: "
      << e.what() << std::endl;
#ifndef _RELEASE_
    DumpCallStack();
#endif
  }

  MPI_Finalize();

  return 0;
}
