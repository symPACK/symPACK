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

  int tasks_per_node = upcxx::local_team().rank_n();
 
  /* NOTE: This algorithm is known to underestimate the amount of space 
   * available on nodes that have more than one physical GPU, so using the -gem_mem argument is encouraged
   * in such cases. See README.md for details.
   */ 
  size_t alloc_size, free, total;
  CUDA_ERROR_CHECK(cudaMemGetInfo(&free, &total));
  alloc_size = (size_t)(free*0.8) / (std::max(tasks_per_node, n_gpus) / n_gpus);
  symPACK::gpu_alloc_size = alloc_size;
}
 
extern "C"
void symPACK_cuda_setup(symPACK::symPACKOptions optionsFact) {
  if (optionsFact.gpu_alloc_size != 0) {
    symPACK::gpu_alloc_size = optionsFact.gpu_alloc_size;  
  }

  symPACK::fallback_type = optionsFact.fallback_type; 

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
  symPACK::gpu_solve = optionsFact.gpu_solve;
  symPACK::gpu_verbose = optionsFact.gpu_verbose;
  
  if (optionsFact.gpu_verbose) {
    if (upcxx::rank_me()==0) {
        std::cout<<"========= GPU CONFIGURATION OPTIONS ========="<<std::endl;
        std::cout<<"-gpu_mem: "<<(symPACK::gpu_alloc_size/(double)(1<<20))<<" MiB"<<std::endl;
        std::cout<<"-gpu_blk: "<<optionsFact.gpu_block_limit<<" bytes"<<std::endl;
        std::cout<<"-trsm_limit: "<<optionsFact.trsm_limit<<" nonzeros"<<std::endl;
        std::cout<<"-potrf_limit: "<<optionsFact.potrf_limit<<" nonzeros"<<std::endl;
        std::cout<<"-gemm_limit: "<<optionsFact.gemm_limit<<" nonzeros"<<std::endl;
        std::cout<<"-syrk_limit: "<<optionsFact.syrk_limit<<" nonzeros"<<std::endl;
    }
    symPACK::cpu_ops.insert({"trsm", 0});
    symPACK::cpu_ops.insert({"potrf", 0});
    symPACK::cpu_ops.insert({"syrk", 0});
    symPACK::cpu_ops.insert({"gemm", 0});
    symPACK::gpu_ops.insert({"trsm", 0});
    symPACK::gpu_ops.insert({"potrf", 0});
    symPACK::gpu_ops.insert({"syrk", 0});
    symPACK::gpu_ops.insert({"gemm", 0});
  }

}
#endif


extern "C"
int symPACK_Init(int *argc, char ***argv){
  int retval = 1;

  // init MPI, if necessary
  // init MPI before GASNet to avoid bug 4638 on HPE Cray EX
  MPI_Initialized(&symPACK::mpi_already_init);
  if (!symPACK::mpi_already_init) MPI_Init(argc, argv);
  MPI_Barrier(MPI_COMM_WORLD);

   // init UPC++
  if ( ! libUPCXXInit ) {
    upcxx::init();
    upcxx::liberate_master_persona();
    symPACK::capture_master_scope();
    libUPCXXInit = true;
  }
#ifdef CUDA_MODE  
  compute_device_alloc_size();
#endif
  upcxx::barrier();

  // renumber MPI ranks to ensure they match UPC++ rank ordering:
  if ( symPACK::world_comm == MPI_COMM_NULL ) MPI_Comm_split(MPI_COMM_WORLD, 0, upcxx::rank_me(), &symPACK::world_comm);
  MPI_Barrier(symPACK::world_comm);

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
  if (symPACK::gpu_verbose) {
    std::cout <<"========= GPU Kernel Info ========="<<std::endl;
    std::cout<<"========= Rank: "<<upcxx::rank_me()<<" ========="<<std::endl;
    for (auto const& op : symPACK::gpu_ops) {
        std::cout<<op.second<<" "<<op.first<<" calls on the GPU"<<std::endl;
    }
    std::cout <<"========= CPU Kernel Info ========="<<std::endl;
    std::cout<<"========= Rank: "<<upcxx::rank_me()<<" ========="<<std::endl;
    for (auto const& op : symPACK::cpu_ops) {
        std::cout<<op.second<<" "<<op.first<<" calls on the CPU"<<std::endl;
    }
  }
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


