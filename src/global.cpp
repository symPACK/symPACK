#include  "sympack/Environment.hpp"
#include  "sympack/utility.hpp"
#include  "sympack/SuperNode.hpp"

#include  "sympack/symPACKMatrix2D.hpp"



bool libMPIInit = false;
bool libUPCXXInit = false;


namespace symPACK{
  int mpi_already_init = 0;
  MPI_Comm world_comm = MPI_COMM_NULL;
  upcxx::persona_scope * master_scope = nullptr;


  void liberate_master_scope(){
    delete master_scope;
    master_scope = nullptr;
  }

  void capture_master_scope() {
    if ( master_scope == nullptr ) {
      master_scope = new upcxx::persona_scope(upcxx::master_persona());
    }
  }
}

extern "C"
int symPACK_Rank(int * rank){
  int retval = 0;
  *rank = upcxx::rank_me();
    return 0;
}


#ifdef CUDA_MODE
void symPACK_cuda_setup() {
  int n_gpus = upcxx::gpu_default_device::device_n();
  symPACK::handlers.reserve(n_gpus); 
  symPACK::cusolver_handlers.reserve(n_gpus);
  cublasStatus_t status;
  cusolverStatus_t cusolver_status;
  int gpu_id = upcxx::rank_me() % n_gpus;

  cublasHandle_t handle;
  symPACK::handlers[gpu_id] = handle; 
  cudaSetDevice(gpu_id);
  CUBLAS_ERROR_CHECK(cublasCreate(&symPACK::handlers[gpu_id]));

  cusolverDnHandle_t cusolver_handle;
  symPACK::cusolver_handlers[gpu_id] = cusolver_handle;
  CUSOLVER_ERROR_CHECK(cusolverDnCreate(&symPACK::cusolver_handlers[gpu_id]));
}
#endif


extern "C"
int symPACK_Init(int *argc, char ***argv){
  int retval = 1;
  // init UPC++
  if ( ! libUPCXXInit ) {
    upcxx::init();
    upcxx::liberate_master_persona();
    symPACK::capture_master_scope();
    libUPCXXInit = true;
  }
#ifdef CUDA_MODE
  symPACK_cuda_setup();
#endif
  // init MPI, if necessary
  MPI_Initialized(&symPACK::mpi_already_init);
  if (!symPACK::mpi_already_init) MPI_Init(argc, argv);
  if ( symPACK::world_comm == MPI_COMM_NULL ) MPI_Comm_split(MPI_COMM_WORLD, 0, upcxx::rank_me(), &symPACK::world_comm);
  return retval;
}





extern "C"
int symPACK_Finalize(){
  int retval = 1;
  int rank = 0;
  symPACK_Rank(&rank);

  MPI_Comm_free(&symPACK::world_comm);
  symPACK::world_comm = MPI_COMM_NULL;

  if (!symPACK::mpi_already_init) MPI_Finalize();

  if(libUPCXXInit){
    symPACK::capture_master_scope();

    upcxx::finalize();
    symPACK::liberate_master_scope();

    libUPCXXInit = false;
    libMPIInit = false;
  }

#ifdef CUDA_MODE
  int n_gpus = upcxx::gpu_default_device::device_n();
  int gpu_id = upcxx::rank_me() % n_gpus;
  cudaSetDevice(gpu_id);
  cublasDestroy(symPACK::handlers[gpu_id]);
  cusolverDnDestroy(symPACK::cusolver_handlers[gpu_id]);
  //symPACK::gpu_allocator.destroy();
#endif

  return retval;
}






namespace symPACK{

  int symPACKMatrixBase::last_id = 0;

  std::map<int, symPACKMatrixBase *  > g_sp_handle_to_matrix;

} // namespace SYMPACK



namespace symPACK{
  namespace Multithreading{ 
    int NumThread = 1;
  } // namespace Multithreading
} // namespace SYMPACK


