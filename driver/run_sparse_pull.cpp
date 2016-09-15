/// @file run_sparse_fanboth.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2014-01-23
//#define _DEBUG_


#include <mpi.h>

#include <time.h>
#include <omp.h>

#include  "sympack.hpp"
#include  "sympack/SupernodalMatrix.hpp"

#include  "sympack/CommTypes.hpp"
#include  "sympack/Ordering.hpp"

#include <upcxx.h>

/******* TYPE used in the computations ********/
#define SCALAR double
//#define SCALAR std::complex<double>

/******* TYPE in the input matrix ********/
#define INSCALAR double


using namespace SYMPACK;


int main(int argc, char **argv) 
{
  MPI_Init(&argc,&argv);
  upcxx::init(&argc, &argv);

  symPACKOptions optionsFact;

  MPI_Comm worldcomm;
  MPI_Comm_dup(MPI_COMM_WORLD,&worldcomm);
  MPI_Comm_size(worldcomm,&np);
  MPI_Comm_rank(worldcomm,&iam);


#if defined(SPROFILE) || defined(PMPI)
  SYMPACK_SPROFILE_INIT(argc, argv);
#endif


  //Initialize a logfile per rank
  logfileptr = new LogFile(iam);
  logfileptr->OFS()<<"********* LOGFILE OF P"<<iam<<" *********"<<endl;
  logfileptr->OFS()<<"**********************************"<<endl;

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

  Int maxIsend = -1;
  if( options.find("-is") != options.end() ){
    maxIsend= atoi(options["-is"].front().c_str());
  }

  Int maxIrecv = -1;
  if( options.find("-ir") != options.end() ){
    maxIrecv= atoi(options["-ir"].front().c_str());
  }


  Int nrhs = 1;
  if( options.find("-nrhs") != options.end() ){
    nrhs= atoi(options["-nrhs"].front().c_str());
  }

  optionsFact.factorization = FANBOTH;
  if( options.find("-fb") != options.end() ){
    if(options["-fb"].front()=="static"){
      optionsFact.factorization = FANBOTH_STATIC;
    }
  }

  if( options.find("-refine") != options.end() ){
    optionsFact.order_refinement_str = options["-refine"].front();
  }

  if( options.find("-fact") != options.end() ){
    if(options["-fact"].front()=="LL"){
      optionsFact.decomposition = LL;
    }
    else if(options["-fact"].front()=="LDL"){
      optionsFact.decomposition = LDL;
    }
  }




  if( options.find("-lb") != options.end() ){
    optionsFact.load_balance_str = options["-lb"].front();
  }

  if( options.find("-ordering") != options.end() ){
    if(options["-ordering"].front()=="AMD"){
      optionsFact.ordering = AMD;
    }
    else if(options["-ordering"].front()=="MMD"){
      optionsFact.ordering = MMD;
    }
    else if(options["-ordering"].front()=="NDBOX"){
      optionsFact.ordering = NDBOX;
    }
    else if(options["-ordering"].front()=="NDGRID"){
      optionsFact.ordering = NDGRID;
    }
#ifdef USE_SCOTCH
    else if(options["-ordering"].front()=="SCOTCH"){
      optionsFact.ordering = SCOTCH;
    }
#endif
#ifdef USE_METIS
    else if(options["-ordering"].front()=="METIS"){
      optionsFact.ordering = METIS;
    }
#endif
#ifdef USE_PARMETIS
    else if(options["-ordering"].front()=="PARMETIS"){
      optionsFact.ordering = PARMETIS;
    }
#endif
#ifdef USE_PTSCOTCH
    else if(options["-ordering"].front()=="PTSCOTCH"){
      optionsFact.ordering = PTSCOTCH;
    }
#endif
    else{
      optionsFact.ordering = NATURAL;
    }
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

  MPI_Comm workcomm;
  MPI_Comm_split(worldcomm,iam<np,iam,&workcomm);
  Int new_rank = (iam<np)?iam:iam-np;
  upcxx::team_all.split(iam<np,new_rank, workteam);

  if(iam<np){
    MPI_Comm_size(workcomm,&np);
    MPI_Comm_rank(workcomm,&iam);


    DistSparseMatrix<SCALAR> HMat(workcomm);
    ReadMatrix<SCALAR,INSCALAR>(filename , informatstr,  HMat);


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
        cout<<"Starting spGEMM"<<endl;
      }

      timeSta = get_time();

      {
        //TODO HANDLE MULTIPLE RHS

      const DistSparseMatrixGraph & Local = HMat.GetLocalGraph();
      Idx firstCol = Local.LocalFirstVertex()+(1-Local.GetBaseval());//1-based
      Idx lastCol = Local.LocalFirstVertex()+Local.LocalVertexCount()+(1-Local.GetBaseval());//1-based

        RHS.assign(n*nrhs,0.0);
        for(Int k = 0; k<nrhs; ++k){
          for(Idx j = 1; j<=n; ++j){
            if(j>=firstCol && j<lastCol){
              Idx iLocal = j-firstCol;//0-based

              //do I own the column ?
              SCALAR t = XTrue[j-1+k*n];
              //do a dense mat mat mul ?
              Ptr colbeg = Local.colptr[iLocal]-(1-Local.GetBaseval());//1-based
              Ptr colend = Local.colptr[iLocal+1]-(1-Local.GetBaseval());//1-based
              for(Ptr ii = colbeg; ii< colend;++ii){
                Idx row = Local.rowind[ii-1]-(1-Local.GetBaseval());//1-based
                RHS[row-1+k*n] += t*HMat.nzvalLocal[ii-1];
                if(row>j){
                  RHS[j-1+k*n] += XTrue[row-1+k*n]*HMat.nzvalLocal[ii-1];
                }
              }
            }
          }
        }
        //Do a reduce of RHS
        mpi::Allreduce((SCALAR*)MPI_IN_PLACE,&RHS[0],RHS.size(),MPI_SUM,workcomm);
      }

      timeEnd = get_time();
      if(iam==0){
        cout<<"spGEMM time: "<<timeEnd-timeSta<<endl;
      }
    }

    if(iam==0){
      cout<<"Starting allocation"<<endl;
    }

#ifdef EXPLICIT_PERMUTE
    std::vector<Int> perm;
#endif
    std::vector<SCALAR> XFinal;
    {
      //do the symbolic factorization and build supernodal matrix
      optionsFact.maxIsend = maxIsend;
      optionsFact.maxIrecv = maxIrecv;


      optionsFact.commEnv = new CommEnvironment(workcomm);
      SupernodalMatrix<SCALAR>*  SMat;

      /************* ALLOCATION AND SYMBOLIC FACTORIZATION PHASE ***********/
#ifndef NOTRY
      try
#endif
      {
        timeSta = get_time();
        SMat = new SupernodalMatrix<SCALAR>();
        SMat->team_ = workteam;
        //SMat->Init2(HMat,optionsFact);
        SMat->Init(HMat,optionsFact);
        SMat->SymbolicFactorization(HMat);
        SMat->DistributeA(HMat);
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
        cout<<"Initialization time: "<<timeEnd-timeSta<<endl;
      }


      /************* NUMERICAL FACTORIZATION PHASE ***********/
      if(iam==0){
        cout<<"Starting Factorization"<<endl;
      }
      timeSta = get_time();
      SYMPACK_TIMER_START(FACTORIZATION);
      SMat->Factorize();
      SYMPACK_TIMER_STOP(FACTORIZATION);
      timeEnd = get_time();

      if(iam==0){
        cout<<"Factorization time: "<<timeEnd-timeSta<<endl;
      }
      logfileptr->OFS()<<"Factorization time: "<<timeEnd-timeSta<<endl;

      SMat->DumpMatlab();

//      SMat->DumpMatlab();

//only for debug purpose
#if 0
      SMat->Dump();
      for(Int k = 0; k<nrhs; ++k){
        for(Int j = 0; j<n; ++j){
          RHS[j*nrhs+k] = 1.0;
        }
      }
      logfileptr->OFS()<<"RHS: "<<RHS<<endl;
#endif


      if(nrhs>0){
        /**************** SOLVE PHASE ***********/
        if(iam==0){
          cout<<"Starting solve"<<endl;
        }
        XFinal = RHS;

        timeSta = get_time();
        SMat->Solve(&XFinal[0],nrhs);
        timeEnd = get_time();

        if(iam==0){
          cout<<"Solve time: "<<timeEnd-timeSta<<endl;
        }

        SMat->GetSolution(&XFinal[0],nrhs);
      }
      delete SMat;

    }

  if(nrhs>0 && XFinal.size()>0) {
      const DistSparseMatrixGraph & Local = HMat.GetLocalGraph();
      Idx firstCol = Local.LocalFirstVertex()+(1-Local.GetBaseval());//1-based
      Idx lastCol = Local.LocalFirstVertex()+Local.LocalVertexCount()+(1-Local.GetBaseval());//1-based

      std::vector<SCALAR> AX(n*nrhs,0.0);

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
      mpi::Allreduce((SCALAR*)MPI_IN_PLACE,&AX[0],AX.size(),MPI_SUM,workcomm);

      if(iam==0){
        blas::Axpy(AX.size(),-1.0,&RHS[0],1,&AX[0],1);
        double normAX = lapack::Lange('F',n,nrhs,&AX[0],n);
        double normRHS = lapack::Lange('F',n,nrhs,&RHS[0],n);
        cout<<"Norm of residual after SPCHOL is "<<normAX/normRHS<<std::endl;
      }
    }



#ifdef _TRACK_MEMORY_
    MemoryAllocator::printStats();
#endif


    delete optionsFact.commEnv;
  }

  MPI_Barrier(workcomm);
  MPI_Comm_free(&workcomm);

  MPI_Barrier(worldcomm);
  MPI_Comm_free(&worldcomm);




  delete logfileptr;


  //This will also finalize MPI
  upcxx::finalize();
  return 0;
}


