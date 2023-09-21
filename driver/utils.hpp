#ifndef _DRIVER_UTILS_HEADER_
#define _DRIVER_UTILS_HEADER_

#include <vector>
#include <map>
#include "sympack.hpp"

template<typename SCALAR>
void generate_rhs( symPACK::DistSparseMatrix<SCALAR> & HMat, std::vector<SCALAR> & RHS, std::vector<SCALAR> & XTrue, int nrhs){
  using namespace symPACK;
  if(nrhs>0){
    MPI_Comm & worldcomm = HMat.comm;
    int np = 1;
    MPI_Comm_size(worldcomm,&np);
    int iam = 0;
    MPI_Comm_rank(worldcomm,&iam);
    int n = HMat.size;
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

    double timeSta = get_time();
    {
      const DistSparseMatrixGraph & Local = HMat.GetLocalGraph();
      Int baseval = Local.GetBaseval();
      Idx firstCol = Local.LocalFirstVertex()+(1-Local.GetBaseval());//1-based
      Idx lastCol = Local.LocalFirstVertex()+Local.LocalVertexCount()+(1-Local.GetBaseval());//1-based

      RHS.assign(n*nrhs,0.0);
      std::vector<Idx> rrowind;
      std::vector<SCALAR> ts(nrhs);
      for(Idx j = 1; j<=n; ++j){
        Int iOwner = 0; for(iOwner = 0; iOwner<np;iOwner++){ if(Local.vertexDist[iOwner]+(1-baseval)<=j && j<Local.vertexDist[iOwner+1]+(1-baseval)){ break; } }


        Ptr nnz = 0;
        Idx * rowind = nullptr;
        SCALAR * nzval = nullptr;

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

    double timeEnd = get_time();
    if(iam==0){
      std::cout<<"spGEMM time: "<<timeEnd-timeSta<<" seconds"<<std::endl;
    }

  }
}


template<typename SCALAR>
void check_solution( symPACK::DistSparseMatrix<SCALAR> & HMat, std::vector<SCALAR> & RHS, std::vector<SCALAR> & XFinal ) {
  using namespace symPACK;
  if(XFinal.size()>0) {
    MPI_Comm & worldcomm = HMat.comm;
    int iam = 0;
    MPI_Comm_rank(worldcomm,&iam);
    int n = HMat.size;
    int nrhs = XFinal.size()/n;
    const DistSparseMatrixGraph & Local = HMat.GetLocalGraph();
    Idx firstCol = Local.LocalFirstVertex()+(1-Local.GetBaseval());//1-based
    Idx lastCol = Local.LocalFirstVertex()+Local.LocalVertexCount()+(1-Local.GetBaseval());//1-based

    std::vector<SCALAR> AX(n*nrhs,SCALAR(0.0));

    logfileptr->OFS()<<"XFinal: "<<XFinal<<std::endl;

    for(Int k = 0; k<nrhs; ++k){
      for(Int j = 1; j<=n; ++j){
        //do I own the column ?
        if(j>=firstCol && j<lastCol){
          Int iLocal = j-firstCol;//0-based
          Int tgtCol = j;
          SCALAR t = XFinal[tgtCol-1+k*n];
          Ptr colbeg = Local.colptr[iLocal]-(1-Local.GetBaseval());//1-based
          Ptr colend = Local.colptr[iLocal+1]-(1-Local.GetBaseval());//1-based
          //do a dense mat mat mul ?
          for(Ptr ii = colbeg; ii< colend;++ii){
            Int row = Local.rowind[ii-1]-(1-Local.GetBaseval());
            Int tgtRow = row;
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
      std::cout<<"Norm of residual after factorization/solve is "<<normAX/normRHS<<std::endl;
    }
  }
}

template<typename vdist_int, typename ptr_t, typename ind_t, typename SCALAR>
void check_solution( int n, vdist_int * vertexDist, ptr_t * colptr, ind_t * rowind, SCALAR * nzvalLocal, std::vector<SCALAR> & RHS, std::vector<SCALAR> & XFinal, MPI_Comm & worldcomm ) {
  using namespace symPACK;
  if(XFinal.size()>0) {
    int iam = 0;
    MPI_Comm_rank(worldcomm,&iam);
    int nrhs = XFinal.size()/n;
    int nLocal = int(vertexDist[iam+1]-vertexDist[iam]);
    int baseval = nLocal>0?colptr[0]:0;

    Idx firstCol = 1;//1-based
    Idx lastCol = baseval+(nLocal==0?0:nLocal-1)+(1-baseval);//1-based

    std::vector<SCALAR> AX(n*nrhs,SCALAR(0.0));

    for(Int k = 0; k<nrhs; ++k){
      for(Int j = 1; j<=n; ++j){
        //do I own the column ?
        if(j>=firstCol && j<lastCol){
          Int iLocal = j-firstCol;//0-based
          Int tgtCol = j;
          SCALAR t = XFinal[tgtCol-1+k*n];
          Ptr colbeg = colptr[iLocal]-(1-baseval);//1-based
          Ptr colend = colptr[iLocal+1]-(1-baseval);//1-based
          //do a dense mat mat mul ?
          for(Ptr ii = colbeg; ii< colend;++ii){
            Int row = rowind[ii-1]-(1-baseval);
            Int tgtRow = row;
            AX[tgtRow-1+k*n] += t*nzvalLocal[ii-1];
            if(row>j){
              AX[tgtCol-1+k*n] += XFinal[tgtRow-1+k*n]*nzvalLocal[ii-1];
            }
          }
        }
      }
    }

    //Do a reduce of RHS
    mpi::Allreduce((SCALAR*)MPI_IN_PLACE,&AX[0],AX.size(),MPI_SUM,worldcomm);

    logfileptr->OFS()<<"XFinal "<<XFinal<<std::endl;
    logfileptr->OFS()<<"AX "<<AX<<std::endl;
    logfileptr->OFS()<<"RHS "<<RHS<<std::endl;

    if(iam==0){
      blas::Axpy(AX.size(),-1.0,&RHS[0],1,&AX[0],1);
      double normAX = lapack::Lange('F',n,nrhs,&AX[0],n);
      double normRHS = lapack::Lange('F',n,nrhs,&RHS[0],n);
      std::cout<<"Norm of residual after facotrization/solve is "<<normAX/normRHS<<std::endl;
    }
  }
}

inline size_t parse_size(const char *desc, const std::vector<std::string> &args, size_t minsz=0) {
    // concatenate all the arguments
    auto argstring = std::accumulate(args.begin(), args.end(), std::string(), 
                        [](const std::string& a, const std::string& b) { return a + b; });
    auto errstr = std::string("error parsing option ") + desc + std::string(" ") + argstring;
    size_t pos = 0;
    double sz_d = std::stod(argstring, &pos);
    std::string unitstr = argstring.substr(pos);
    unitstr.erase(0, unitstr.find_first_not_of(" \t\n")); // strip

    size_t sz;
    if (unitstr.size() == 0) {
      sz = sz_d; // no units, assume bytes
    } else {
      switch (unitstr[0]) {
        case 'b': case 'B': sz = sz_d; break;
        case 'k': case 'K': sz = sz_d * (1<<10); break;
        case 'm': case 'M': sz = sz_d * (1<<20); break;
        case 'g': case 'G': sz = sz_d * (1<<30); break;
        case 't': case 'T': sz = sz_d * (1ULL<<40); break;
        default: 
          throw std::runtime_error(errstr);
      }
    }
    if (sz < minsz) throw std::runtime_error(errstr + ": value too small");
    
    return sz;
}

inline void process_options(int argc, char **argv, symPACK::symPACKOptions & optionsFact,std::string & filename, std::string & informatstr, bool & complextype, int & nrhs){
  using namespace symPACK;
  using option_t = std::map<std::string,std::vector<std::string> > ;
  // *********************************************************************
  // Input parameter
  // *********************************************************************
  option_t options;
  OptionsCreate(argc, argv, options);
  

  if( options.find("-in") != options.end() ){
    filename= options["-in"].front();
  }

  informatstr = "HARWELL_BOEING";
  if( options.find("-inf") != options.end() ){
    informatstr= options["-inf"].front();
  }

  complextype=false;
  if (options.find("-z") != options.end()){
    complextype = true;
  }

// Options related to GPU functionality
#ifdef CUDA_MODE
 optionsFact.gpu_solve = false;
 if (options.find("-gpu_solve") != options.end()) {
    optionsFact.gpu_solve = true;
 } 

 optionsFact.gpu_alloc_size = 0;
 if (options.find("-gpu_mem") != options.end()) {
    optionsFact.gpu_alloc_size = parse_size("-gpu_mem", options["-gpu_mem"], 1<<20);
 } 

 optionsFact.gpu_block_limit = 100000;
 if (options.find("-gpu_blk") != options.end()) {
    optionsFact.gpu_block_limit = parse_size("-gpu_blk", options["-gpu_blk"]);
 } 
 
 optionsFact.trsm_limit = 90000;
 if (options.find("-trsm_limit") != options.end()) {
    optionsFact.trsm_limit = parse_size("-trsm_limit", options["-trsm_limit"]);
 } 
 
 optionsFact.potrf_limit = 1500000;
 if (options.find("-potrf_limit") != options.end()) {
    optionsFact.potrf_limit = parse_size("-potrf_limit", options["-potrf_limit"]);
 } 
 
 optionsFact.gemm_limit = 500000;
 if (options.find("-gemm_limit") != options.end()) {
    optionsFact.gemm_limit = parse_size("-gemm_limit", options["-gemm_limit"]);
 } 

 optionsFact.syrk_limit = 200000;
 if (options.find("-syrk_limit") != options.end()) {
    optionsFact.syrk_limit = parse_size("-syrk_limit", options["-syrk_limit"]);
 } 
 
 optionsFact.gpu_verbose = false;
 if (options.find("-gpu_v") != options.end()) {
    optionsFact.gpu_verbose = true;
 }

 optionsFact.fallback_type = FallbackType::TERMINATE;
 if (options.find("-fallback") !=options.end()) {
    std::string fallback = options["-fallback"].front();
    if (fallback.compare("terminate")==0) {
        //do nothing, since this is the default option
    } else if (fallback.compare("cpu")==0) {
        optionsFact.fallback_type = FallbackType::CPU;
    } else {
        throw std::invalid_argument("Error: " + fallback + " is not a valid fallback option. Use 'terminate' or 'cpu' instad."); 
    }
 }
#endif
  
  //-----------------------------------------------------------------
  optionsFact.memory_limit=-1.0;
  if (options.find("-mem") != options.end()){
    optionsFact.memory_limit = atof(options["-mem"].front().c_str()) ;
  }

  optionsFact.print_stats=false;
  if (options.find("-ps") != options.end()){
    optionsFact.print_stats = atoi(options["-ps"].front().c_str()) == 1;
  }

  if( options.find("-dumpPerm") != options.end() ){
    optionsFact.dumpPerm = atoi(options["-dumpPerm"].front().c_str());
  }
  //-----------------------------------------------------------------

  optionsFact.orderingStr = "MMD" ;
  if( options.find("-ordering") != options.end() ){
    optionsFact.orderingStr = options["-ordering"].front();
  }

  if( options.find("-npord") != options.end() ){
    optionsFact.NpOrdering= atoi(options["-npord"].front().c_str());
  }

  if( options.find("-refine") != options.end() ){
    optionsFact.order_refinement_str = options["-refine"].front();
  }

  //-----------------------------------------------------------------
  Int verbose = 0;
  if( options.find("-v") != options.end() ){
    verbose= atoi(options["-v"].front().c_str());
  }
  optionsFact.verbose = verbose;

  Int maxSnode = 0;
  if( options.find("-b") != options.end() ){
    if(options["-b"].front() != "inf"){
      maxSnode= atoi(options["-b"].front().c_str());
    }
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

  //-----------------------------------------------------------------

  Int numThreads = 1;
  if( options.find("-t") != options.end() ){
    numThreads = atoi(options["-t"].front().c_str());
  }

  if( options.find("-dist") != options.end() ){
    if(options["-dist"].front() != "2D"){
      optionsFact.distribution = SYMPACK_DATA_1D;
    }
    else {
      optionsFact.distribution = SYMPACK_DATA_2D;
    }
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

  if( options.find("-lb") != options.end() ){
    optionsFact.load_balance_str = options["-lb"].front();
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

  optionsFact.mappingTypeStr = "ROW2D";
  if( options.find("-map") != options.end() ){
    optionsFact.mappingTypeStr = options["-map"].front();
  }

  //-----------------------------------------------------------------

  optionsFact.factorization = FANBOTH;
  if( options.find("-alg") != options.end() ){
    if(options["-alg"].front()=="FANBOTH"){
      optionsFact.factorization = FANBOTH;
    }
    else if(options["-alg"].front()=="FANOUT"){
      optionsFact.factorization = FANOUT;
    }
  }

  if( options.find("-fact") != options.end() ){
    if(options["-fact"].front()=="LL"){
      optionsFact.decomposition = DecompositionType::LL;
    }
    else if(options["-fact"].front()=="LDL"){
      optionsFact.decomposition = DecompositionType::LDL;
    }
  }

  //-----------------------------------------------------------------

  nrhs = 0;
  if( options.find("-nrhs") != options.end() ){
    nrhs= atoi(options["-nrhs"].front().c_str());
  }

  //-----------------------------------------------------------------
  optionsFact.maxIsend = maxIsend;
  optionsFact.maxIrecv = maxIrecv;
  optionsFact.numThreads = numThreads;

  //Set up the threading environment
  Multithreading::NumThread = numThreads;
}


#endif // _DRIVER_UTILS_HEADER_
