//#include  "sympack/supernodalTaskGraph.hpp"
//#include  "sympack/symPACKMatrix.hpp"
#include  "sympack.hpp"
//#include  "sympack/Environment.hpp"
//#include  "sympack/utility.hpp"
//#include  "sympack/SuperNode.hpp"
//

    using namespace symPACK;


namespace symPACK{

  struct symPACK_handle{
    std::unique_ptr<DistSparseMatrixBase> inputMat;
    std::unique_ptr<symPACKMatrixBase> factorMat;
  };
//
//
  std::map<int, symPACK_handle  > symPACK_handles;

  //int last_id = 0;
}

  //returns a integer corresponding to a symPACK handle
extern "C" 
  int symPACK_C_InitInstanceFloat(MPI_Comm ccomm) {}
extern "C" 
  int symPACK_InitInstanceFloat(MPI_Fint * Fcomm) {}


//extern "C" 
//  int symPACK_C_InitInstanceDouble(MPI_Comm ccomm){
//    symPACK_Init(nullptr,nullptr);
//
//    symPACK::symPACKOptions optionsFact;
//    optionsFact.MPIcomm = ccomm;
//    optionsFact.decomposition = symPACK::DecompositionType::LDL;
//    optionsFact.verbose = 0;
//   
//    symPACK::symPACK_handle handle;
//    handle.factorMat = std::unique_ptr< symPACK::symPACKMatrixBase > ( new symPACK::symPACKMatrix<double>() );//std::static_pointer_cast<symPACKMatrixBase >(SMat);
//    handle.inputMat = std::unique_ptr< symPACK::DistSparseMatrixBase>( new symPACK::DistSparseMatrix<double>(optionsFact.MPIcomm) );
//    symPACK::symPACK_handles[symPACK::last_id++] = std::move(handle);
//    //TODO restore this
//    static_cast<symPACK::symPACKMatrix<double> *>(handle.factorMat.get())->Init(optionsFact);
//
//    return symPACK::last_id-1;
//  }

extern "C" 
  int symPACK_InitInstanceDouble(MPI_Fint * Fcomm){
    MPI_Comm comm =  MPI_COMM_NULL;
    comm = MPI_Comm_f2c(*Fcomm);

    return symPACK_C_InitInstanceDouble(comm);
  }

extern "C" 
  int symPACK_C_InitInstanceComplex(MPI_Comm ccomm){
  
  }

extern "C" 
  int symPACK_InitInstanceComplex(MPI_Fint * Fcomm){

  }

extern "C" 
  int symPACK_C_InitInstanceDoubleComplex(MPI_Comm ccomm){

  }

extern "C" 
  int symPACK_InitInstanceDoubleComplex(MPI_Fint * Fcomm){

  }


extern "C" 
  void symPACK_SymbolicFactorize(int * sp_handle, int * n, int * colptr , int * rowind){
    auto & handle = symPACK_handles[*sp_handle];

    if(auto pSMat = dynamic_cast<symPACKMatrix<double> * >( handle.factorMat.get() )){
      auto pHMat = dynamic_cast<DistSparseMatrix<double> * >(handle.inputMat.get() );
      symPACK_LoadGraph(pHMat,n,colptr,rowind);
      pSMat->SymbolicFactorization(*pHMat);
    }
    //else if(auto pSMat = std::dynamic_pointer_cast<symPACKMatrix<std::complex<double> > >(handle.factorMat)){
    //  auto pHMat = std::dynamic_pointer_cast<DistSparseMatrix<std::complex<double> > >(handle.inputMat);
    //  symPACK_LoadGraph(pHMat,n,colptr,rowind);
    //  pSMat->SymbolicFactorization(*pHMat);
    //}
    //else if(auto pSMat = std::dynamic_pointer_cast<symPACKMatrix<float> >(handle.factorMat)){
    //  auto pHMat = std::dynamic_pointer_cast<DistSparseMatrix<float> >(handle.inputMat);
    //  symPACK_LoadGraph(pHMat,n,colptr,rowind);
    //  pSMat->SymbolicFactorization(*pHMat);
    //}
    //else if(auto pSMat = std::dynamic_pointer_cast<symPACKMatrix<std::complex<float> > >(handle.factorMat)){
    //  auto pHMat = std::dynamic_pointer_cast<DistSparseMatrix<std::complex<float> > >(handle.inputMat);
    //  symPACK_LoadGraph(pHMat,n,colptr,rowind);
    //  pSMat->SymbolicFactorization(*pHMat);
    //}

  }


extern "C" 
  void symPACK_DistributeFloat(int * sp_handle, float * nzvals){
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
  }

extern "C" 
  void symPACK_DistributeDoubleComplex(int * sp_handle, double * nzvals){
  }


extern "C" 
  void symPACK_NumericalFactorize(int * sp_handle){
    auto & handle = symPACK_handles[*sp_handle];
    if(auto pSMat = dynamic_cast<symPACKMatrix<double> * >( handle.factorMat.get() )){
      pSMat->Factorize();
    }
    else if(auto pSMat = dynamic_cast<symPACKMatrix<std::complex<double> > * >( handle.factorMat.get() )){
      pSMat->Factorize();
    }
    else if(auto pSMat = dynamic_cast<symPACKMatrix<float> * >( handle.factorMat.get() )){
      pSMat->Factorize();
    }
    else if(auto pSMat = dynamic_cast<symPACKMatrix<std::complex<float> > * >( handle.factorMat.get() )){
      pSMat->Factorize();
    }
    else{
      abort();
    }
  }

extern "C" 
  void symPACK_NumericalSolveFloat(int * sp_handle, int * nrhs, float * rhs){
    auto & handle = symPACK_handles[*sp_handle];
    auto pSMat = dynamic_cast<symPACKMatrix<float> * >( handle.factorMat.get() );
    pSMat->Solve(rhs,*nrhs);
    pSMat->GetSolution(rhs,*nrhs);
  }

extern "C" 
  void symPACK_NumericalSolveDouble(int * sp_handle, int * nrhs, double * rhs){
    auto & handle = symPACK_handles[*sp_handle];
    auto pSMat = dynamic_cast<symPACKMatrix<double> * >( handle.factorMat.get() );
    pSMat->Solve(rhs,*nrhs);
    pSMat->GetSolution(rhs,*nrhs);
  }

extern "C" 
  void symPACK_NumericalSolveComplex(int * sp_handle, int * nrhs, float * rhs){
    std::complex<float> * prhs = (std::complex<float> *)rhs;
    auto & handle = symPACK_handles[*sp_handle];
    auto pSMat = dynamic_cast<symPACKMatrix<std::complex<float> > * >( handle.factorMat.get() );
    pSMat->Solve(prhs,*nrhs);
    pSMat->GetSolution(prhs,*nrhs);
  }

extern "C" 
  void symPACK_NumericalSolveDoubleComplex(int * sp_handle, int * nrhs, double * rhs){
    std::complex<double> * prhs = (std::complex<double> *)rhs;
    auto & handle = symPACK_handles[*sp_handle];
    auto pSMat = dynamic_cast<symPACKMatrix<std::complex<double> > * >( handle.factorMat.get() );
    pSMat->Solve(prhs,*nrhs);
    pSMat->GetSolution(prhs,*nrhs);
  }

extern "C" 
  void symPACK_FinalizeInstance(int * sp_handle){
    //auto handle = symPACK_handles[*sp_handle];
    //handle.factorMat.reset(nullptr);
    //handle.inputMat.reset(nullptr);
    symPACK_handles.erase(*sp_handle);

    if(symPACK_handles.empty()){
      symPACK_Finalize();
    }


  }

extern "C" 
  int symPACK_C_InitInstanceDouble(MPI_Comm ccomm){
    //static int last_id = 0;
    symPACK_Init(nullptr,nullptr);

    symPACKOptions optionsFact;
    optionsFact.MPIcomm = ccomm;
    optionsFact.decomposition = DecompositionType::LDL;
    optionsFact.verbose = 0;
   
    symPACK_handle handle;
    handle.factorMat = std::unique_ptr< symPACKMatrixBase > ( new symPACKMatrix<double>() );//std::static_pointer_cast<symPACKMatrixBase >(SMat);
    handle.inputMat = std::unique_ptr< DistSparseMatrixBase>( new DistSparseMatrix<double>(optionsFact.MPIcomm) );
    //int handle_id = last_id++;
    //store the handle id in the matrix as well
    int handle_id =handle.factorMat->sp_handle;

    symPACK_handles[handle_id] = std::move(handle);
    //TODO restore this
    static_cast<symPACKMatrix<double> *>(symPACK_handles[handle_id].factorMat.get())->Init(optionsFact);

    return handle_id;
  }


