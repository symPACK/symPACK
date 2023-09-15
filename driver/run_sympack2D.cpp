#include <mpi.h>

#include "sympack.hpp"
#include <sympack/symPACKMatrix.hpp>
#include <sympack/symPACKMatrix2D.hpp>
#include "sympack/CommTypes.hpp"
#include "sympack/Ordering.hpp"

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
    logfileptr->OFS()<<"Matrix dimension: " << n << std::endl;
    std::vector<SCALAR> RHS,XTrue;
    generate_rhs(HMat,RHS,XTrue,nrhs);

    if (upcxx::rank_me()==0)
        std::cout<<"Ordering: "<<optionsFact.orderingStr<<std::endl;
    
    std::vector<SCALAR> XFinal;
    auto SMat2D = std::make_shared<symPACKMatrix2D<Ptr,Idx,SCALAR> >();
    try{
#ifdef CUDA_MODE
      logfileptr->OFS()<< "CUDA Mode enabled" << std::endl;
      symPACK_cuda_setup(optionsFact);
      upcxx::barrier();
#endif
      //do the symbolic factorization and build supernodal matrix
      /************* ALLOCATION AND SYMBOLIC FACTORIZATION PHASE ***********/
      SMat2D->Init(optionsFact);
      timeSta = get_time();
      SMat2D->SymbolicFactorization(HMat);
      logfileptr->OFS()<<"Distributing Matrix"<<std::endl;
      SMat2D->DistributeMatrix(HMat);
      timeEnd = get_time();
      if(upcxx::rank_me()==0){
        std::cout<<"Initialization time: "<<timeEnd-timeSta<<" seconds"<<std::endl;
      }
      
      timeSta = get_time();
      SMat2D->Factorize();
      timeEnd = get_time();
      if(iam==0){
        std::cout<<"Factorization time: "<<timeEnd-timeSta<<" seconds"<<std::endl;
      }
      logfileptr->OFS()<<"Factorization time: "<<timeEnd-timeSta<<" seconds"<<std::endl;
    }
    catch(const std::bad_alloc& e){
      std::cout << "Allocation failed: " << e.what() << '\n';
      SMat2D = nullptr;
      abort();
    }
    catch(const std::runtime_error& e){
      std::cerr << "Runtime error: " << e.what() << '\n';
    }


    /**************** SOLVE PHASE ***********/
    if (nrhs>0){
      if(iam==0){
        std::cout<<"Starting solve 2D"<<std::endl;
      }
      XFinal = RHS;

      timeSta = get_time();
      SMat2D->Solve(&XFinal[0],nrhs, XFinal.size());
      timeEnd = get_time();

      if(iam==0){
        std::cout<<"Solve 2D time: "<<timeEnd-timeSta<<" seconds"<<std::endl;
      }

      SMat2D->GetSolution(&XFinal[0],nrhs);

      check_solution(HMat,RHS,XFinal);
    }

    upcxx::barrier();
    MPI_Barrier(worldcomm);
    MPI_Comm_free(&worldcomm);
    delete logfileptr;
    delete statfileptr;
  }
  //This will also finalize MPI
  symPACK_Finalize();
  return 0;
}


