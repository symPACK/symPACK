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
extern "C"
void compute_device_alloc_size() {
  int n_gpus = upcxx::gpu_default_device::device_n();
  if (n_gpus==0) {
    std::cerr<<"Error: found no CUDA devices"<<std::endl;
    abort();
  }

  MPI_Comm shmcomm;
  int tasks_per_node;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
  MPI_Comm_size(shmcomm, &tasks_per_node);
  
  size_t alloc_size, free, total;
  CUDA_ERROR_CHECK(cudaMemGetInfo(&free, &total));
  alloc_size = (free/2) / (std::max(tasks_per_node, n_gpus) / n_gpus);
  symPACK::gpu_alloc_size = alloc_size;
}
 
extern "C"
void symPACK_cuda_setup(symPACK::symPACKOptions optionsFact) {
  if (optionsFact.gpu_alloc_size != 0) {
    symPACK::gpu_alloc_size = optionsFact.gpu_alloc_size;  
  } 

  symPACK::gpu_allocator = upcxx::make_gpu_allocator<upcxx::gpu_default_device>(symPACK::gpu_alloc_size);

  int gpu_id = symPACK::gpu_allocator.device_id();
  cudaSetDevice(gpu_id);
  
  cublasHandle_t handle;
  symPACK::cublas_handler= handle; 
  CUBLAS_ERROR_CHECK(cublasCreate(&symPACK::cublas_handler));

  cusolverDnHandle_t cusolver_handle;
  symPACK::cusolver_handler = cusolver_handle;
  CUSOLVER_ERROR_CHECK(cusolverDnCreate(&symPACK::cusolver_handler));
  
  symPACK::gpu_block_limit = optionsFact.gpu_block_limit;
  symPACK::trsm_limit = optionsFact.trsm_limit;
  symPACK::potrf_limit = optionsFact.potrf_limit;
  symPACK::gemm_limit = optionsFact.gemm_limit;
  symPACK::syrk_limit = optionsFact.syrk_limit;
}
#endif


extern "C"
int symPACK_Init(int *argc, char ***argv){
  int retval = 1;
  // init MPI, if necessary
  MPI_Initialized(&symPACK::mpi_already_init);
  if (!symPACK::mpi_already_init) MPI_Init(argc, argv);
  int iam = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &iam);
  if ( symPACK::world_comm == MPI_COMM_NULL ) MPI_Comm_split(MPI_COMM_WORLD, 0, iam, &symPACK::world_comm);
  MPI_Barrier(MPI_COMM_WORLD);
   // init UPC++
  if ( ! libUPCXXInit ) {
    upcxx::init();
    upcxx::liberate_master_persona();
    symPACK::capture_master_scope();
    libUPCXXInit = true;
  }
  upcxx::barrier();
#ifdef CUDA_MODE  
  compute_device_alloc_size();
#endif
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

  
#ifdef CUDA_MODE
  cublasDestroy(symPACK::cublas_handler);
  cusolverDnDestroy(symPACK::cusolver_handler);
#endif
  if(libUPCXXInit){
    symPACK::capture_master_scope();

    upcxx::finalize();
    symPACK::liberate_master_scope();

    libUPCXXInit = false;
    libMPIInit = false;
  }

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


