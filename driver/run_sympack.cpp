#include <mpi.h>

#include  "sympack.hpp"
#include <sympack/symPACKMatrix.hpp>
#include  "sympack/CommTypes.hpp"
#include  "sympack/Ordering.hpp"

#include "utils.hpp"

/******* TYPE used in the computations ********/
#define SCALAR double
//#define SCALAR std::complex<double>

/******* TYPE in the input matrix ********/
#define RSCALAR double
#define CSCALAR std::complex<double>


using namespace symPACK;


int main(int argc, char **argv) 
{
  int success = symPACK_Init(&argc,&argv);
  if (success==-1)
    return 1;
  {
    int iam = 0;
    int np = 1;
    MPI_Comm worldcomm;
    MPI_Comm_size(MPI_COMM_WORLD,&np);
    assert(np==upcxx::rank_n());
    symPACK_Rank(&iam);
    MPI_Comm_split(MPI_COMM_WORLD, 0, upcxx::rank_me(), &worldcomm);

    MPI_Comm_size(worldcomm,&np);
    symPACK_Rank(&iam);

    //Initialize a logfile per rank
    logfileptr = new LogFile(iam);
    logfileptr->OFS()<<"********* LOGFILE OF P"<<iam<<" *********"<<std::endl;
    logfileptr->OFS()<<"**********************************"<<std::endl;

    // *********************************************************************
    // Input parameter
    // *********************************************************************
    symPACKOptions optionsFact;
    std::string filename;
    std::string informatstr;
    bool complextype=false;
    int nrhs = 0;
    process_options(argc, argv, optionsFact, filename, informatstr, complextype, nrhs);
    //-----------------------------------------------------------------

    Real timeSta, timeEnd;

    Int all_np = np;
    np = optionsFact.used_procs(np);
    optionsFact.MPIcomm = worldcomm;

    DistSparseMatrix<SCALAR> HMat(worldcomm);
    if(complextype){
      ReadMatrix<SCALAR,CSCALAR>(filename , informatstr,  HMat);
    }
    else{
      ReadMatrix<SCALAR,RSCALAR>(filename , informatstr,  HMat);
    }

    Int n = HMat.size;
    std::vector<SCALAR> RHS,XTrue;
    generate_rhs(HMat,RHS,XTrue,nrhs);


    std::vector<SCALAR> XFinal;
    auto SMat = std::make_shared<symPACKMatrix<SCALAR> >();
    try
    {
#ifdef CUDA_MODE
    logfileptr->OFS()<< "CUDA Mode enabled" << std::endl;
    symPACK_cuda_setup(optionsFact);
    upcxx::barrier();
#endif
      //do the symbolic factorization and build supernodal matrix
      /************* ALLOCATION AND SYMBOLIC FACTORIZATION PHASE ***********/
      timeSta = get_time();
      SMat->Init(optionsFact);
      SMat->SymbolicFactorization(HMat);
      SMat->DistributeMatrix(HMat);
      timeEnd = get_time();

      if(iam==0){
        std::cout<<"Initialization time: "<<timeEnd-timeSta<<" seconds"<<std::endl;
      }

      /************* NUMERICAL FACTORIZATION PHASE ***********/
      timeSta = get_time();
      SYMPACK_TIMER_START(FACTORIZATION);
      SMat->Factorize();
      SYMPACK_TIMER_STOP(FACTORIZATION);
      timeEnd = get_time();

      if(iam==0){
        std::cout<<"Factorization time: "<<timeEnd-timeSta<<" seconds"<<std::endl;
      }
      logfileptr->OFS()<<"Factorization time: "<<timeEnd-timeSta<<" seconds"<<std::endl;


    }
    catch(const std::bad_alloc& e){
      std::cout << "Allocation failed: " << e.what() << '\n';
      SMat = nullptr;
      abort();
    }
    catch(const std::runtime_error& e){
      std::cerr << "Runtime error: " << e.what() << '\n';
    }


    /**************** SOLVE PHASE ***********/
    if(nrhs>0){
      if(iam==0){
        std::cout<<"Starting solve"<<std::endl;
      }
      XFinal = RHS;

      timeSta = get_time();
      SMat->Solve(&XFinal[0],nrhs, XFinal.size());
      timeEnd = get_time();

      if(iam==0){
        std::cout<<"Solve time: "<<timeEnd-timeSta<<" seconds"<<std::endl;
      }

      SMat->GetSolution(&XFinal[0],nrhs);

      check_solution(HMat,RHS,XFinal);
    }

    MPI_Barrier(worldcomm);
    MPI_Comm_free(&worldcomm);
    delete logfileptr;
  }
  //This will also finalize MPI
  symPACK_Finalize();
  return 0;
}


