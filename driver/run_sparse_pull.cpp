/// @file run_sparse_fanboth.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2014-01-23
//#define _DEBUG_


#include <mpi.h>

#include <time.h>
#include <omp.h>
#include <sympack/symPACKMatrix.hpp>

#include  "sympack.hpp"
#include  "sympack/CommTypes.hpp"
#include  "sympack/Ordering.hpp"


//#define DUMP_MATLAB
//#define DUMP_MATLAB_SOL

/******* TYPE used in the computations ********/
#define SCALAR double
#define SCALAR std::complex<double>

/******* TYPE in the input matrix ********/
#define RSCALAR double
#define CSCALAR std::complex<double>


using namespace symPACK;


int main(int argc, char **argv) 
{
  if(0){
    std::vector<SCALAR,AlignedAllocator<SCALAR> > tmp(5);
    tmp[0]=1.0;
    tmp[1]=2.0;
    tmp[2]=3.0;
    tmp[3]=4.0;
    tmp[4]=5.0;

    for(auto && val: tmp){
      std::cout<<val<<" ";
    }
    std::cout<<std::endl;

    tmp.resize(0);
  }


  symPACK_Init(&argc,&argv);

  symPACKOptions optionsFact;

  int iam = 0;
  int np = 1;
  MPI_Comm worldcomm;
  MPI_Comm_dup(MPI_COMM_WORLD,&worldcomm);
  MPI_Comm_size(worldcomm,&np);
  MPI_Comm_rank(worldcomm,&iam);


#if defined(SPROFILE) || defined(PMPI)
  SYMPACK_SPROFILE_INIT(argc, argv);
#endif

  //Initialize a logfile per rank
  logfileptr = new LogFile(iam);
  logfileptr->OFS()<<"********* LOGFILE OF P"<<iam<<" *********"<<std::endl;
  logfileptr->OFS()<<"**********************************"<<std::endl;

  // *********************************************************************
  // Input parameter
  // *********************************************************************
  std::map<std::string,std::vector<std::string> > options;

  OptionsCreate(argc, argv, options);

  std::string filename;
  if( options.find("-in") != options.end() ){
    filename= options["-in"].front();
  }
  std::string informatstr;
  if( options.find("-inf") != options.end() ){
    informatstr= options["-inf"].front();
  }

  Int maxSnode = -1;
  if( options.find("-b") != options.end() ){
    maxSnode= atoi(options["-b"].front().c_str());
  }

  optionsFact.print_stats=false;
  if (options.find("-ps") != options.end()){
    optionsFact.print_stats = atoi(options["-ps"].front().c_str()) == 1;
  }

  bool complextype=false;
  if (options.find("-z") != options.end()){
    complextype = true;
  }
  
  if( options.find("-npord") != options.end() ){
    optionsFact.NpOrdering= atoi(options["-npord"].front().c_str());
  }

  optionsFact.relax.SetMaxSize(maxSnode);
  if( options.find("-relax") != options.end() ){
    if(options["-relax"].size()==3){
      optionsFact.relax.SetNrelax0(atoi(options["-relax"][0].c_str()));
      optionsFact.relax.SetNrelax1(atoi(options["-relax"][1].c_str()));
      optionsFact.relax.SetNrelax2(atoi(options["-relax"][2].c_str()));
    }
    else{
      //disable relaxed supernodes
      optionsFact.relax.SetNrelax0(0);
    }
  }



  Int numThreads = 1;
  if( options.find("-t") != options.end() ){
    numThreads = atoi(options["-t"].front().c_str());
  }


  Int maxIsend = -1;
  if( options.find("-is") != options.end() ){
    if(options["-is"].front() != "inf"){
      maxIsend= atoi(options["-is"].front().c_str());
    }
  }

  Int maxIrecv = -1;
  if( options.find("-ir") != options.end() ){
    if(options["-ir"].front() != "inf"){
      maxIrecv= atoi(options["-ir"].front().c_str());
    }
  }


  Int nrhs = 1;
  if( options.find("-nrhs") != options.end() ){
    nrhs= atoi(options["-nrhs"].front().c_str());
  }

  if( options.find("-dumpPerm") != options.end() ){
    optionsFact.dumpPerm = atoi(options["-dumpPerm"].front().c_str());
  }

  optionsFact.factorization = FANBOTH;
  if( options.find("-alg") != options.end() ){
    if(options["-alg"].front()=="FANBOTH"){
      optionsFact.factorization = FANBOTH;
    }
    else if(options["-alg"].front()=="FANOUT"){
      optionsFact.factorization = FANOUT;
    }
  }


  if( options.find("-refine") != options.end() ){
    optionsFact.order_refinement_str = options["-refine"].front();
  }

  if( options.find("-fact") != options.end() ){
    if(options["-fact"].front()=="LL"){
      optionsFact.decomposition = DecompositionType::LL;
    }
    else if(options["-fact"].front()=="LDL"){
      optionsFact.decomposition = DecompositionType::LDL;
    }
  }

  if( options.find("-lb") != options.end() ){
    optionsFact.load_balance_str = options["-lb"].front();
  }

  optionsFact.orderingStr = "MMD" ;
  if( options.find("-ordering") != options.end() ){
    optionsFact.orderingStr = options["-ordering"].front();
  }

  if( options.find("-scheduler") != options.end() ){
    if(options["-scheduler"].front()=="MCT"){
      optionsFact.scheduler = MCT;
    }
    if(options["-scheduler"].front()=="FIFO"){
      optionsFact.scheduler = FIFO;
    }
    else if(options["-scheduler"].front()=="PR"){
      optionsFact.scheduler = PR;
    }
    else{
      optionsFact.scheduler = DL;
    }
  }



  Real timeSta, timeEnd;


  if( options.find("-map") != options.end() ){
    optionsFact.mappingTypeStr = options["-map"].front();
  }
  else{
    optionsFact.mappingTypeStr = "ROW2D";
  }


  Int all_np = np;
  np = optionsFact.used_procs(np);

  //should be worldcomm at some point
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
    if(nrhs>0){
      RHS.resize(n*nrhs);
      XTrue.resize(n*nrhs);

      //Initialize XTrue;
      Int val = 1.0;
      for(Int i = 0; i<n;++i){ 
        for(Int j=0;j<nrhs;++j){
          XTrue[i+j*n] = val;
          val = -val;
        }
      }



      if(iam==0){
        std::cout<<"Starting spGEMM"<<std::endl;
      }

      timeSta = get_time();

      {
        const DistSparseMatrixGraph & Local = HMat.GetLocalGraph();
        Int baseval = Local.GetBaseval();
        Idx firstCol = Local.LocalFirstVertex()+(1-Local.GetBaseval());//1-based
        Idx lastCol = Local.LocalFirstVertex()+Local.LocalVertexCount()+(1-Local.GetBaseval());//1-based

        RHS.assign(n*nrhs,0.0);
        std::vector<Idx> rrowind;
        std::vector<SCALAR> rnzval;
        std::vector<SCALAR> ts(nrhs);
        for(Idx j = 1; j<=n; ++j){
          Int iOwner = 0; for(iOwner = 0; iOwner<np;iOwner++){ if(Local.vertexDist[iOwner]+(1-baseval)<=j && j<Local.vertexDist[iOwner+1]+(1-baseval)){ break; } }


          Ptr nnz = 0;
          Idx * rowind = NULL;
          SCALAR * nzval = NULL;

          //do I own the column ?

          if(iam==iOwner){
            Idx iLocal = j-firstCol;//0-based
            Ptr colbeg = Local.colptr[iLocal]-(1-Local.GetBaseval());//1-based
            Ptr colend = Local.colptr[iLocal+1]-(1-Local.GetBaseval());//1-based
            nnz = colend-colbeg;
            nzval = &HMat.nzvalLocal[colbeg-1];
            rowind = const_cast<Idx*>(&Local.rowind[colbeg-1]);
          }

          if(iam==iOwner){
            for(Int k = 0; k<nrhs; ++k){
              ts[k] = XTrue[j-1+k*n];
            }
            //do a dense mat mat mul ?
            for(Ptr ii = 1; ii<= nnz;++ii){
              Idx row = rowind[ii-1]-(1-Local.GetBaseval());//1-based
              for(Int k = 0; k<nrhs; ++k){
                RHS[row-1+k*n] += ts[k]*nzval[ii-1];
                if(row>j){
                  RHS[j-1+k*n] += XTrue[row-1+k*n]*nzval[ii-1];
                }
              }
            }
          }
        }
      }

      mpi::Allreduce( (SCALAR*)MPI_IN_PLACE, &RHS[0],  n*nrhs, MPI_SUM, worldcomm);

      timeEnd = get_time();
      if(iam==0){
        std::cout<<"spGEMM time: "<<timeEnd-timeSta<<std::endl;
      }

#ifdef DUMP_MATLAB
          {
            std::size_t N = RHS.size();
            logfileptr->OFS()<<"RHS = [ ";
            for(std::size_t i = 0; i<N;i++){
              logfileptr->OFS()<<ToMatlabScalar(RHS[i])<<" ";
            }
            logfileptr->OFS()<<"];"<<std::endl;
          }
#endif

#ifdef DUMP_MATLAB
          {
            std::size_t N = XTrue.size();
            logfileptr->OFS()<<"XTrue = [ ";
            for(std::size_t i = 0; i<N;i++){
              logfileptr->OFS()<<ToMatlabScalar(XTrue[i])<<" ";
            }
            logfileptr->OFS()<<"];"<<std::endl;
          }
#endif

    }





    if(iam==0){
      std::cout<<"Starting allocation"<<std::endl;
    }

#ifdef EXPLICIT_PERMUTE
    std::vector<Int> perm;
#endif
    std::vector<SCALAR> XFinal;
    {
      //do the symbolic factorization and build supernodal matrix
      optionsFact.maxIsend = maxIsend;
      optionsFact.maxIrecv = maxIrecv;
      optionsFact.numThreads = numThreads;

      //Set up the threading environment
      Multithreading::NumThread = numThreads;


      symPACKMatrix<SCALAR>*  SMat;

      /************* ALLOCATION AND SYMBOLIC FACTORIZATION PHASE ***********/
#ifndef NOTRY
      try
#endif
      {
        timeSta = get_time();
        SMat = new symPACKMatrix<SCALAR>();
        SMat->Init(optionsFact);
        SMat->SymbolicFactorization(HMat);
        SMat->DistributeMatrix(HMat);
        timeEnd = get_time();
#ifdef EXPLICIT_PERMUTE
        perm = SMat->GetOrdering().perm;
#endif
      }
#ifndef NOTRY
      catch(const std::bad_alloc& e){
        std::cout << "Allocation failed: " << e.what() << '\n';
        SMat = NULL;
        abort();
      }
#endif

      if(iam==0){
        std::cout<<"Initialization time: "<<timeEnd-timeSta<<std::endl;
      }

#ifdef DUMP_MATLAB
      if(iam==0){
        logfileptr->OFS()<<"A= ";
      }
      SMat->DumpMatlab();
#endif
#if 0
      {
        logfileptr->OFS()<<"Supernode partition"<<SMat->GetSupernodalPartition()<<std::endl;
      }
#endif


      /************* NUMERICAL FACTORIZATION PHASE ***********/
      if(iam==0){
        std::cout<<"Starting Factorization"<<std::endl;
      }
      timeSta = get_time();
      SYMPACK_TIMER_START(FACTORIZATION);
      SMat->Factorize();
      SYMPACK_TIMER_STOP(FACTORIZATION);
      timeEnd = get_time();

      if(iam==0){
        std::cout<<"Factorization time: "<<timeEnd-timeSta<<std::endl;
      }
      logfileptr->OFS()<<"Factorization time: "<<timeEnd-timeSta<<std::endl;

#ifdef DUMP_MATLAB
      if(iam==0){
        logfileptr->OFS()<<"L= ";
      }
      SMat->DumpMatlab();
#endif

      if(nrhs>0){
        /**************** SOLVE PHASE ***********/
        if(iam==0){
          std::cout<<"Starting solve"<<std::endl;
        }
        XFinal = RHS;

        timeSta = get_time();
        SMat->Solve(&XFinal[0],nrhs);
        timeEnd = get_time();

        if(iam==0){
          std::cout<<"Solve time: "<<timeEnd-timeSta<<std::endl;
        }

        SMat->GetSolution(&XFinal[0],nrhs);

#if defined(DUMP_MATLAB) || defined(DUMP_MATLAB_SOL)
        if(nrhs>0 && XFinal.size()>0) {
          {
            std::size_t N = XFinal.size();
            logfileptr->OFS()<<"XFinal = [ ";
            for(std::size_t i = 0; i<N;i++){
              logfileptr->OFS()<<ToMatlabScalar(XFinal[i])<<" ";
            }
            logfileptr->OFS()<<"];"<<std::endl;
          }
        }
#endif


      }
      delete SMat;

    }

  if(nrhs>0 && XFinal.size()>0) {
      const DistSparseMatrixGraph & Local = HMat.GetLocalGraph();
      Idx firstCol = Local.LocalFirstVertex()+(1-Local.GetBaseval());//1-based
      Idx lastCol = Local.LocalFirstVertex()+Local.LocalVertexCount()+(1-Local.GetBaseval());//1-based

      std::vector<SCALAR> AX(n*nrhs,SCALAR(0.0));

      for(Int k = 0; k<nrhs; ++k){
        for(Int j = 1; j<=n; ++j){
          //do I own the column ?
          if(j>=firstCol && j<lastCol){
            Int iLocal = j-firstCol;//0-based
#ifdef EXPLICIT_PERMUTE
            Int tgtCol = perm[j-1];
#else
            Int tgtCol = j;
#endif
            SCALAR t = XFinal[tgtCol-1+k*n];
            Ptr colbeg = Local.colptr[iLocal]-(1-Local.GetBaseval());//1-based
            Ptr colend = Local.colptr[iLocal+1]-(1-Local.GetBaseval());//1-based
            //do a dense mat mat mul ?
            for(Ptr ii = colbeg; ii< colend;++ii){
              Int row = Local.rowind[ii-1]-(1-Local.GetBaseval());
#ifdef EXPLICIT_PERMUTE
              Int tgtRow = perm[row-1];
#else
              Int tgtRow = row;
#endif
              AX[tgtRow-1+k*n] += t*HMat.nzvalLocal[ii-1];
              if(row>j){
                AX[tgtCol-1+k*n] += XFinal[tgtRow-1+k*n]*HMat.nzvalLocal[ii-1];
              }
            }
          }
        }
      }

      //Do a reduce of RHS
      mpi::Allreduce((SCALAR*)MPI_IN_PLACE,&AX[0],AX.size(),MPI_SUM,worldcomm);

      if(iam==0){
        blas::Axpy(AX.size(),-1.0,&RHS[0],1,&AX[0],1);
        double normAX = lapack::Lange('F',n,nrhs,&AX[0],n);
        double normRHS = lapack::Lange('F',n,nrhs,&RHS[0],n);
        std::cout<<"Norm of residual after SPCHOL is "<<normAX/normRHS<<std::endl;
      }
    }



#ifdef _TRACK_MEMORY_
    MemoryAllocator::printStats();
#endif



  MPI_Barrier(worldcomm);
  MPI_Comm_free(&worldcomm);


  delete logfileptr;

  //This will also finalize MPI
  symPACK_Finalize();
  return 0;
}


