#include <mpi.h>

#include <time.h>
#include <omp.h>
#include <sympack/symPACKMatrix.hpp>

#include  "sympack.hpp"
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
  int handle = -1;
  symPACK_Init(&argc,&argv);
  {
    int iam = 0;
    int np = 1;
    MPI_Comm worldcomm;
    MPI_Comm_size(MPI_COMM_WORLD,&np);
    
    //std::cout<<np<<" "<<upcxx::rank_n()<<std::endl;
    //assert(np==upcxx::rank_n());

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
    //this is a dummy object
    symPACKOptions optionsFact;
    std::string filename;
    std::string informatstr;
    bool complextype=false;
    int nrhs = 0;
    process_options(argc, argv, optionsFact, filename, informatstr, complextype, nrhs);
    bool is2D = optionsFact.distribution == SYMPACK_DATA_2D;
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

  if(iam==0){
    std::cout<<"Starting allocation"<<std::endl;
  }

  std::vector<SCALAR> XFinal;
  {
    //Gather enverything on P0
    auto graph = HMat.GetLocalGraph();
    std::vector<int> colptr(graph.colptr.size());
    for(Idx col = 0; col< colptr.size(); col++){ colptr[col] = graph.colptr[col]; }
    std::vector<int> rowind(graph.rowind.size());
    for(Ptr ptr = 0; ptr< rowind.size(); ptr++){ rowind[ptr] = graph.rowind[ptr]; }
    std::vector<SCALAR> nzval(HMat.nzvalLocal.size());
    for(Ptr ptr = 0; ptr< nzval.size(); ptr++){ nzval[ptr] = HMat.nzvalLocal[ptr]; }

    std::vector<int> sizes(np);
    std::vector<int> displs(np+1);

    //auto vertexDistbak = graph.vertexDist;

    {
      int cnt = (colptr.size())*sizeof(int);
      MPI_Gather(&cnt,sizeof(int),MPI_BYTE,sizes.data(),sizeof(int),MPI_BYTE,0,worldcomm);
      displs[0] = 0;
      std::partial_sum(sizes.begin(),sizes.end(),displs.begin()+1);
      std::vector<int> tcolptr(iam==0?displs.back()/sizeof(int):2);

      MPI_Gatherv(colptr.data(), cnt, MPI_BYTE,
          tcolptr.data(), sizes.data(), displs.data(),MPI_BYTE,0,worldcomm);

      if(iam==0){
        //make colptr great again
        for(int p = 1; p < np; p++){
          int offset = displs[p]/sizeof(int);

          //for(int i = 0; i < sizes[p]/sizeof(int); i++){ logfileptr->OFS()<<tcolptr[offset+i]<<" "; } logfileptr->OFS()<<std::endl;

          int prev_entry_offset = displs[p]/sizeof(int)-p; 
          int increment = tcolptr[prev_entry_offset] - graph.GetBaseval();
          for(int i = 0; i < sizes[p]/sizeof(int); i++){
            tcolptr[prev_entry_offset+i] = tcolptr[offset+i] + increment; 
          }

          //for(int i = 0; i < sizes[p]/sizeof(int); i++){ logfileptr->OFS()<<tcolptr[prev_entry_offset+i]<<" "; } logfileptr->OFS()<<std::endl;
        }
        tcolptr.resize(HMat.size+1);
        tcolptr.back() = HMat.nnz+graph.GetBaseval();
      }
      else{
        tcolptr[0] = graph.GetBaseval();
        tcolptr[1] = tcolptr[0];
      }
      colptr.swap(tcolptr);
    }
    {
      int cnt = rowind.size()*sizeof(int);
      MPI_Gather(&cnt,sizeof(int),MPI_BYTE,sizes.data(),sizeof(int),MPI_BYTE,0,worldcomm);
      displs[0] = 0;
      std::partial_sum(sizes.begin(),sizes.end(),displs.begin()+1);
      std::vector<int> trowind(iam==0?displs.back()/sizeof(int):0);
      MPI_Gatherv(rowind.data(), rowind.size()*sizeof(int), MPI_BYTE,
          trowind.data(), sizes.data(), displs.data(),MPI_BYTE,0,worldcomm);
      rowind.swap(trowind);
    }

    {
      int cnt = nzval.size()*sizeof(SCALAR);
      MPI_Gather(&cnt,sizeof(int),MPI_BYTE,sizes.data(),sizeof(int),MPI_BYTE,0,worldcomm);
      displs[0] = 0;
      std::partial_sum(sizes.begin(),sizes.end(),displs.begin()+1);
      std::vector<SCALAR> tnzval(iam==0?displs.back()/sizeof(SCALAR):0);
      MPI_Gatherv(nzval.data(), nzval.size()*sizeof(SCALAR), MPI_BYTE,
          tnzval.data(), sizes.data(), displs.data(),MPI_BYTE,0,worldcomm);
      nzval.swap(tnzval);
    }
    for(int p = 1; p<=np;p++){ graph.vertexDist[p] = HMat.size+graph.GetBaseval(); }


    timeSta = get_time();

    handle = symPACK_C_InitInstanceDouble(worldcomm/*,is2D*/);
    symPACK_SymbolicFactorize(&handle, &HMat.size, colptr.data() , rowind.data() );
    symPACK_DistributeDouble(&handle, nzval.data() );
    timeEnd = get_time();

    if(iam==0){
      std::cout<<"Initialization time: "<<timeEnd-timeSta<<std::endl;
    }


    /************* NUMERICAL FACTORIZATION PHASE ***********/
    if(iam==0){
      std::cout<<"Starting Factorization"<<std::endl;
    }
    timeSta = get_time();
    SYMPACK_TIMER_START(FACTORIZATION);
    symPACK_NumericalFactorize(&handle);
    SYMPACK_TIMER_STOP(FACTORIZATION);
    timeEnd = get_time();

    if(iam==0){
      std::cout<<"Factorization time: "<<timeEnd-timeSta<<std::endl;
    }
    logfileptr->OFS()<<"Factorization time: "<<timeEnd-timeSta<<std::endl;

    if(nrhs>0){
      /**************** SOLVE PHASE ***********/
      if(iam==0){
        std::cout<<"Starting solve"<<std::endl;
      }
      XFinal = RHS;

      timeSta = get_time();
      symPACK_NumericalSolveDouble(&handle, &nrhs, XFinal.data());

      timeEnd = get_time();

      if(iam==0){
        std::cout<<"Solve time: "<<timeEnd-timeSta<<std::endl;
      }

      logfileptr->OFS()<<XFinal<<std::endl;
      check_solution( HMat, RHS, XFinal );
      //check_solution( HMat.size, graph.vertexDist.data(),colptr.data(),rowind.data(), nzval.data(), RHS , XFinal, worldcomm );
    }
  }

  MPI_Barrier(worldcomm);
  MPI_Comm_free(&worldcomm);
  delete logfileptr;
  }


  symPACK_FinalizeInstance(&handle);
  return 0;
}


