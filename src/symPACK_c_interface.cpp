#include  "sympack.hpp"

using namespace symPACK;


namespace symPACK{

  struct symPACK_handle{
    std::unique_ptr<DistSparseMatrixBase> inputMat;
    std::unique_ptr<symPACKMatrixBase> factorMat;
  };
  std::map<int, symPACK_handle  > symPACK_handles;
}

  //returns a integer corresponding to a symPACK handle
extern "C" 
  int symPACK_C_InitInstanceFloat(MPI_Comm ccomm, bool is2D) {
    symPACK_Init(nullptr,nullptr);

    symPACKOptions optionsFact;
    optionsFact.MPIcomm = ccomm;
    optionsFact.decomposition = DecompositionType::LDL;
    optionsFact.verbose = 0;
   
    symPACK_handle handle;
    if (!is2D) handle.factorMat = std::unique_ptr< symPACKMatrixBase > ( new symPACKMatrix<float>() );
    else handle.factorMat = std::unique_ptr< symPACKMatrixBase > ( new symPACKMatrix2D< Ptr, Idx, float>() );
    
    handle.inputMat = std::unique_ptr< DistSparseMatrixBase>( new DistSparseMatrix<float>(optionsFact.MPIcomm) );
    //store the handle id in the matrix as well
    int handle_id =handle.factorMat->sp_handle;

    symPACK_handles[handle_id] = std::move(handle);
    //TODO restore this
    symPACK_handles[handle_id].factorMat->Init(optionsFact);

    return handle_id;
  }

extern "C" 
  int symPACK_InitInstanceFloat(MPI_Fint * Fcomm, int is2D) { 
    MPI_Comm comm =  MPI_COMM_NULL;
    comm = MPI_Comm_f2c(*Fcomm);
    return symPACK_C_InitInstanceFloat(comm, is2D == 0);
  }

extern "C" 
  int symPACK_InitInstanceDouble(MPI_Fint * Fcomm, int is2D){
    MPI_Comm comm =  MPI_COMM_NULL;
    comm = MPI_Comm_f2c(*Fcomm);
    return symPACK_C_InitInstanceDouble(comm, is2D == 0);
  }

extern "C" 
  int symPACK_C_InitInstanceComplex(MPI_Comm ccomm, bool is2D){
    symPACK_Init(nullptr,nullptr);

    symPACKOptions optionsFact;
    optionsFact.MPIcomm = ccomm;
    optionsFact.decomposition = DecompositionType::LDL;
    optionsFact.verbose = 0;
   
    symPACK_handle handle;
    if (!is2D) handle.factorMat = std::unique_ptr< symPACKMatrixBase > ( new symPACKMatrix<std::complex<float>>() );
    else handle.factorMat = std::unique_ptr< symPACKMatrixBase > ( new symPACKMatrix2D< Ptr, Idx, std::complex<float>>() );
    
    handle.inputMat = std::unique_ptr< DistSparseMatrixBase>( new DistSparseMatrix<std::complex<float>>(optionsFact.MPIcomm) );
    //store the handle id in the matrix as well
    int handle_id =handle.factorMat->sp_handle;

    symPACK_handles[handle_id] = std::move(handle);
    //TODO restore this
    symPACK_handles[handle_id].factorMat->Init(optionsFact);

    return handle_id;
  }

extern "C" 
  int symPACK_InitInstanceComplex(MPI_Fint * Fcomm, int is2D){
    MPI_Comm comm =  MPI_COMM_NULL;
    comm = MPI_Comm_f2c(*Fcomm);
    return symPACK_C_InitInstanceComplex(comm, is2D == 0);
  }

extern "C" 
  int symPACK_C_InitInstanceDoubleComplex(MPI_Comm ccomm, bool is2D){
    symPACK_Init(nullptr,nullptr);

    symPACKOptions optionsFact;
    optionsFact.MPIcomm = ccomm;
    optionsFact.decomposition = DecompositionType::LDL;
    optionsFact.verbose = 0;
   
    symPACK_handle handle;
    if (!is2D) handle.factorMat = std::unique_ptr< symPACKMatrixBase > ( new symPACKMatrix<std::complex<double>>() );
    else handle.factorMat = std::unique_ptr< symPACKMatrixBase > ( new symPACKMatrix2D< Ptr, Idx, std::complex<double>>() );
    
    handle.inputMat = std::unique_ptr< DistSparseMatrixBase>( new DistSparseMatrix<std::complex<double>>(optionsFact.MPIcomm) );
    //store the handle id in the matrix as well
    int handle_id =handle.factorMat->sp_handle;

    symPACK_handles[handle_id] = std::move(handle);
    //TODO restore this
    symPACK_handles[handle_id].factorMat->Init(optionsFact);

    return handle_id;
  }

extern "C" 
  int symPACK_InitInstanceDoubleComplex(MPI_Fint * Fcomm, int is2D){
    MPI_Comm comm =  MPI_COMM_NULL;
    comm = MPI_Comm_f2c(*Fcomm);
    return symPACK_C_InitInstanceDoubleComplex(comm, is2D == 0);
  }


extern "C" 
  void symPACK_SymbolicFactorize(int * sp_handle, int * n, int * colptr , int * rowind){
    auto & handle = symPACK_handles[*sp_handle];

    //if(auto pSMat = dynamic_cast<symPACKMatrix<double> * >( handle.factorMat.get() )){
    //  auto pHMat = dynamic_cast<DistSparseMatrix<double> * >(handle.inputMat.get() );
    //  symPACK_LoadGraph(pHMat,n,colptr,rowind);
    //  pSMat->SymbolicFactorization(*pHMat);
    //}
    //else if(auto pSMat = dynamic_cast<symPACKMatrix2D<double> * >( handle.factorMat.get() )){
    //  auto pHMat = dynamic_cast<DistSparseMatrix<double> * >(handle.inputMat.get() );
    //  symPACK_LoadGraph(pHMat,n,colptr,rowind);
    //  pSMat->SymbolicFactorization(*pHMat);
    //}
    if (auto pSMat = dynamic_cast<symPACKMatrixMeta<float> * >( handle.factorMat.get() )) {
      auto pHMat = dynamic_cast<DistSparseMatrix<float> * >(handle.inputMat.get() );
      symPACK_LoadGraph(pHMat,n,colptr,rowind);
      pSMat->SymbolicFactorization(*pHMat);
    }
    else if (auto pSMat = dynamic_cast<symPACKMatrixMeta<double> * >( handle.factorMat.get() )) {
      auto pHMat = dynamic_cast<DistSparseMatrix<double> * >(handle.inputMat.get() );
      symPACK_LoadGraph(pHMat,n,colptr,rowind);
      pSMat->SymbolicFactorization(*pHMat);
    }
    else if (auto pSMat = dynamic_cast<symPACKMatrixMeta<std::complex<float>> * >( handle.factorMat.get() )) {
      auto pHMat = dynamic_cast<DistSparseMatrix<std::complex<float>> * >(handle.inputMat.get() );
      symPACK_LoadGraph(pHMat,n,colptr,rowind);
      pSMat->SymbolicFactorization(*pHMat);
    }
    else if (auto pSMat = dynamic_cast<symPACKMatrixMeta<std::complex<double>> * >( handle.factorMat.get() )) {
      auto pHMat = dynamic_cast<DistSparseMatrix<std::complex<double>> * >(handle.inputMat.get() );
      symPACK_LoadGraph(pHMat,n,colptr,rowind);
      pSMat->SymbolicFactorization(*pHMat);
    }
    else {
      abort();
    }
  }


extern "C" 
  void symPACK_DistributeFloat(int * sp_handle, float * nzvals){
    auto & handle = symPACK_handles[*sp_handle];
    auto pSMat = dynamic_cast<symPACKMatrix<float> * >( handle.factorMat.get() );
    auto pHMat = dynamic_cast<DistSparseMatrix<float> * >(handle.inputMat.get() );

    int mpirank = 0;
    MPI_Comm_rank(pHMat->comm,&mpirank);
    int mpisize = 1;
    MPI_Comm_size(pHMat->comm,&mpisize);
    auto & graph = pHMat->GetLocalGraph();
    auto nnzlocal = graph.LocalEdgeCount();
    pHMat->nzvalLocal.resize(nnzlocal);

    if(mpirank==0){
      for(Ptr ptr = 0; ptr < nnzlocal; ptr++){
        pHMat->nzvalLocal[ptr] = nzvals[ptr];
      }
    }

    pSMat->DistributeMatrix(*pHMat);
  }

extern "C" 
  void symPACK_DistributeDouble(int * sp_handle, double * nzvals){
    auto & handle = symPACK_handles[*sp_handle];
    auto pSMat = dynamic_cast<symPACKMatrix<double> * >( handle.factorMat.get() );
    auto pHMat = dynamic_cast<DistSparseMatrix<double> * >(handle.inputMat.get() );

    int mpirank = 0;
    MPI_Comm_rank(pHMat->comm,&mpirank);
    int mpisize = 1;
    MPI_Comm_size(pHMat->comm,&mpisize);
    auto & graph = pHMat->GetLocalGraph();
    auto nnzlocal = graph.LocalEdgeCount();
    pHMat->nzvalLocal.resize(nnzlocal);

    if(mpirank==0){
      for(Ptr ptr = 0; ptr < nnzlocal; ptr++){
        pHMat->nzvalLocal[ptr] = nzvals[ptr];
      }
    }

    pSMat->DistributeMatrix(*pHMat);
  }

extern "C" 
  void symPACK_DistributeComplex(int * sp_handle, float * nzvals){
    auto & handle = symPACK_handles[*sp_handle];
    auto pSMat = dynamic_cast<symPACKMatrix<std::complex<float>> * >( handle.factorMat.get() );
    auto pHMat = dynamic_cast<DistSparseMatrix<std::complex<float>> * >(handle.inputMat.get() );

    int mpirank = 0;
    MPI_Comm_rank(pHMat->comm,&mpirank);
    int mpisize = 1;
    MPI_Comm_size(pHMat->comm,&mpisize);
    auto & graph = pHMat->GetLocalGraph();
    auto nnzlocal = graph.LocalEdgeCount();
    pHMat->nzvalLocal.resize(nnzlocal);

    if(mpirank==0){
      for(Ptr ptr = 0; ptr < nnzlocal; ptr++){
        pHMat->nzvalLocal[ptr] = std::complex<float>(nzvals[2*ptr],nzvals[2*ptr+1]);
      }
    }

    pSMat->DistributeMatrix(*pHMat);
  }

extern "C" 
  void symPACK_DistributeDoubleComplex(int * sp_handle, double * nzvals){
    auto & handle = symPACK_handles[*sp_handle];
    auto pSMat = dynamic_cast<symPACKMatrix<std::complex<double>> * >( handle.factorMat.get() );
    auto pHMat = dynamic_cast<DistSparseMatrix<std::complex<double>> * >(handle.inputMat.get() );

    int mpirank = 0;
    MPI_Comm_rank(pHMat->comm,&mpirank);
    int mpisize = 1;
    MPI_Comm_size(pHMat->comm,&mpisize);
    auto & graph = pHMat->GetLocalGraph();
    auto nnzlocal = graph.LocalEdgeCount();
    pHMat->nzvalLocal.resize(nnzlocal);

    if(mpirank==0){
      for(Ptr ptr = 0; ptr < nnzlocal; ptr++){
        pHMat->nzvalLocal[ptr] = std::complex<double>(nzvals[2*ptr],nzvals[2*ptr+1]);
      }
    }

    pSMat->DistributeMatrix(*pHMat);
  }


extern "C" 
  void symPACK_NumericalFactorize(int * sp_handle){
    auto & handle = symPACK_handles[*sp_handle];
    auto pSMat = ( handle.factorMat.get() );
    pSMat->Factorize();
    //if(auto pSMat = dynamic_cast<symPACKMatrix<double> * >( handle.factorMat.get() )){
    //  pSMat->Factorize();
    //}
    //else if(auto pSMat = dynamic_cast<symPACKMatrix<std::complex<double> > * >( handle.factorMat.get() )){
    //  pSMat->Factorize();
    //}
    //else if(auto pSMat = dynamic_cast<symPACKMatrix<float> * >( handle.factorMat.get() )){
    //  pSMat->Factorize();
    //}
    //else if(auto pSMat = dynamic_cast<symPACKMatrix<std::complex<float> > * >( handle.factorMat.get() )){
    //  pSMat->Factorize();
    //}
    //else{
    //  abort();
    //}
  }

extern "C" 
  void symPACK_NumericalSolveFloat(int * sp_handle, int * nrhs, float * rhs){
    auto & handle = symPACK_handles[*sp_handle];
    auto pSMat = dynamic_cast<symPACKMatrixMeta<float> * >( handle.factorMat.get() );
    pSMat->Solve(rhs,*nrhs);
    pSMat->GetSolution(rhs,*nrhs);
  }

extern "C" 
  void symPACK_NumericalSolveDouble(int * sp_handle, int * nrhs, double * rhs){
    auto & handle = symPACK_handles[*sp_handle];
    auto pSMat = dynamic_cast<symPACKMatrixMeta<double> * >( handle.factorMat.get() );
    pSMat->Solve(rhs,*nrhs);
    pSMat->GetSolution(rhs,*nrhs);
  }

extern "C" 
  void symPACK_NumericalSolveComplex(int * sp_handle, int * nrhs, float * rhs){
    std::complex<float> * prhs = (std::complex<float> *)rhs;
    auto & handle = symPACK_handles[*sp_handle];
    auto pSMat = dynamic_cast<symPACKMatrixMeta<std::complex<float> > * >( handle.factorMat.get() );
    pSMat->Solve(prhs,*nrhs);
    pSMat->GetSolution(prhs,*nrhs);
  }

extern "C" 
  void symPACK_NumericalSolveDoubleComplex(int * sp_handle, int * nrhs, double * rhs){
    std::complex<double> * prhs = (std::complex<double> *)rhs;
    auto & handle = symPACK_handles[*sp_handle];
    auto pSMat = dynamic_cast<symPACKMatrixMeta<std::complex<double> > * >( handle.factorMat.get() );
    pSMat->Solve(prhs,*nrhs);
    pSMat->GetSolution(prhs,*nrhs);
  }

extern "C" 
  void symPACK_FinalizeInstance(int * sp_handle){
    symPACK_handles.erase(*sp_handle);

    if(symPACK_handles.empty()){
      symPACK_Finalize();
    }


  }

extern "C" 
  int symPACK_C_InitInstanceDouble(MPI_Comm ccomm, bool is2D){
    symPACK_Init(nullptr,nullptr);

    symPACKOptions optionsFact;
    optionsFact.MPIcomm = ccomm;
    optionsFact.decomposition = DecompositionType::LDL;
    optionsFact.verbose = 0;
   
    symPACK_handle handle;
    if (!is2D) handle.factorMat = std::unique_ptr< symPACKMatrixBase > ( new symPACKMatrix<double>() );
    else handle.factorMat = std::unique_ptr< symPACKMatrixBase > ( new symPACKMatrix2D< Ptr, Idx, double>() );
    
    handle.inputMat = std::unique_ptr< DistSparseMatrixBase>( new DistSparseMatrix<double>(optionsFact.MPIcomm) );
    //store the handle id in the matrix as well
    int handle_id =handle.factorMat->sp_handle;

    symPACK_handles[handle_id] = std::move(handle);
    //TODO restore this
    symPACK_handles[handle_id].factorMat->Init(optionsFact);

    return handle_id;
  }


