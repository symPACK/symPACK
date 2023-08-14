#ifndef _ENVIRONMENT_DECL_HPP_
#define _ENVIRONMENT_DECL_HPP_

//#define _MAP_DEBUG_
//#define SP_THREADS

#include "sympack_definitions.hpp"
#include "sympack_config.hpp"
#ifdef CUDA_MODE
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cublas_v2.h"
#include "cusolverDn.h"
#endif

#include <sys/time.h>
inline double get_time()
  {
    struct timeval tv;
    gettimeofday(&tv, 0);
    return tv.tv_sec + ((double) tv.tv_usec / 1000000);
  }


#include <upcxx/upcxx.hpp>

#if defined(CUDA_MODE) && !UPCXX_KIND_CUDA
#error "symPACK's CUDA mode requires the UPC++ library to be configured with CUDA support."
#endif

//debug
#include <sys/types.h>
#include <unistd.h>
#include <execinfo.h>
#include <exception>



// STL libraries
#include <iostream> 
#include <iomanip> 
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <cfloat>
#include <cstdint>
#include <complex>
#include <string>
#include <set>
#include <map>
#include <stack>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <mutex>

// MPI
#include <mpi.h>


namespace symPACK{
  void liberate_master_scope();
  void capture_master_scope();


  extern MPI_Comm world_comm;
#ifdef CUDA_MODE
  extern cublasHandle_t cublas_handler;
  extern cusolverDnHandle_t cusolver_handler;
  extern std::vector<cudaStream_t> streams;
  extern upcxx::device_allocator<upcxx::cuda_device> gpu_allocator;
  extern size_t gpu_alloc_size, gpu_block_limit, trsm_limit, potrf_limit, gemm_limit, syrk_limit;
  extern bool gpu_solve;
  enum class FallbackType {TERMINATE, CPU};
  extern FallbackType fallback_type;
  extern bool gpu_verbose;
#endif
  extern std::map<std::string, int> cpu_ops, gpu_ops;
}

namespace symPACK{
#ifdef SP_THREADS
  using upcxx_mutex_type = std::recursive_mutex;
  extern upcxx_mutex_type upcxx_mutex;
#endif
}

namespace symPACK{
  namespace Multithreading{ 
    extern int NumThread;
  } // namespace Multithreading
} // namespace SYMPACK



#include "sympack/datatypes.hpp"

/***********************************************************************
 *  Error handling
 **********************************************************************/

namespace symPACK{
  template <typename Map>
  inline void increment_counter(Map &map, std::string const &op) {
  #ifdef CUDA_MODE
    if (symPACK::gpu_verbose) {
        map[op]+=1;
    }
  #endif
  }

  inline void gdb_lock(){
      pid_t pid = getpid();
int iam = upcxx::rank_me();
      std::cerr<<"P"<<iam<<" is locked, pid is "<<pid<<std::endl;
    volatile int lock = 1;
    while (lock == 1){ }
      std::cerr<<"P"<<iam<<" is unlocked"<<std::endl;
  }

  inline void gdb_lock(Int proc){
int iam = upcxx::rank_me();
    if(iam==proc){
      pid_t pid = getpid();
      std::cerr<<"P"<<iam<<" is locked, pid is "<<pid<<std::endl;
      volatile int lock = 1;
      while (lock == 1){ }
      std::cerr<<"P"<<iam<<" is unlocked"<<std::endl;
    }
  }


#ifndef NDEBUG
  inline void bassert(bool cond){
    if(!(cond)){
      gdb_lock();
    }
  }
#else
  #define bassert(cond)
#endif
}


#ifndef _USE_NUMVEC_
namespace symPACK{
template <typename T> using vector = std::vector<T>;
}
#else
#include "NumVec.hpp"
namespace symPACK{
  template <typename T> using vector = symPACK::NumVec<T,size_t>;
}
#endif


// Always use complex data.
#define _USE_COMPLEX_
#include "sympack/Types.hpp"



#include "sympack/LogFile.hpp"


#define CUBLAS_ERROR_CHECK(s)                                                                       \
do {                                                                                                   \
    cublasStatus_t err;                                                                             \
    if ((err = (s)) != CUBLAS_STATUS_SUCCESS)                                                       \
    {                                                                                               \
        std::cout << "cuBLAS Error " << err << " at " << __FILE__ << ":" << __LINE__ << "\n";       \
        exit(1);                                                                                    \
    }                                                                                               \
} while (0)

#define CUDA_ERROR_CHECK(s)                                                                         \
do {                                                                                                   \
    cudaError_t error = s;                                                                          \
    if (error != cudaSuccess) {                                                                     \
        std::cout << "CUDA Error " << error << " at " << __FILE__ << ":" << __LINE__ << "\n";       \
        std::cout << cudaGetErrorString(error) << "\n";                                             \
        exit(1);                                                                                    \
    }                                                                                               \
} while (0)

#define CUSOLVER_ERROR_CHECK(s)                                                                          \
do {                                                                                                        \
    cusolverStatus_t error = s;                                                                          \
    if (error != CUSOLVER_STATUS_SUCCESS) {                                                              \
        std::cout << "cuSOLVER Error " << error << " at " << __FILE__ << ":" << __LINE__ << "\n";        \
        exit(1);                                                                                         \
    }                                                                                                    \
} while (0)







// TODO Remove environment_impl.hpp. Move things to utility.hpp and only
// keep environment.hpp
// Update numXXX_*.hpp and tinyvec*.hpp

// *********************************************************************
// Redefine the global macros
// *********************************************************************


// The verbose level of debugging information
#ifdef  DEBUG
#define _DEBUGlevel_ DEBUG
#endif

// Release mode. For speed up the calculation and reduce verbose level.
// Note that RELEASE overwrites DEBUG level.
#ifdef RELEASE
#define _RELEASE_
#define _DEBUGlevel -1
#endif


//extern LogFile * logfileptr;


/***********************************************************************
 *  Data types and constants
 **********************************************************************/

namespace symPACK{

// Basic data types

//#define SYMSOLVE_FC FC_GLOBAL(symsolve, SYMSOLVE)

#define FORTRAN(name) FC_GLOBAL(name, DUMMY)
#define BLAS(name) FC_GLOBAL(name, DUMMY)
#define LAPACK(name) FC_GLOBAL(name, DUMMY)

//#ifndef Add_
//
//#define FORTRAN(name) name
//#define BLAS(name) name
//#define LAPACK(name) name
//#else
//#define FORTRAN(name) name##_
//#define BLAS(name) name##_
//#define LAPACK(name) name##_
//#endif

// IO
extern  std::ofstream  statusOFS;

// *********************************************************************
// Define constants
// *********************************************************************
// Commonly used
const char UPPER = 'U';
const char LOWER = 'L';


} // namespace SYMPACK



namespace symPACK{

// We define an output stream that does nothing. This is done so that the 
// root process can be used to print data to a file's ostream while all other 
// processes use a null ostream. 
struct NullStream : std::ostream
{            
	struct NullStreamBuffer : std::streambuf
	{
		Int overflow( Int c ) { return traits_type::not_eof(c); }
	} nullStreamBuffer_;

	NullStream() 
		: std::ios(&nullStreamBuffer_), std::ostream(&nullStreamBuffer_)
		{ }
};  

/////////////////////////////////////////////

class ExceptionTracer
{
public:
	ExceptionTracer()
	{
		void * array[25];
		int nSize = backtrace(array, 25);
		char ** symbols = backtrace_symbols(array, nSize);

		for (int i = 0; i < nSize; i++)
		{
			symPACKOS << symbols[i] << std::endl;
		}

		free(symbols);
	}
};

// *********************************************************************
// Global utility functions 
// These utility functions do not depend on local definitions
// *********************************************************************
// Return the closest integer to a real number
Int iround( Real a );

// Read the options from command line
void OptionsCreate(Int argc, char** argv, 
		std::map<std::string,std::string>& options);

} // namespace SYMPACK


// *********************************************************************
// Profiling functions 
// *********************************************************************
#include "sympack/timer.hpp"

#define VAL(str) #str
#define TOSTRING(str) VAL(str)


#if defined (SPROFILE)
#define SYMPACK_TIMER_START(a) SYMPACK_FSTART(a);
#define SYMPACK_TIMER_STOP(a) SYMPACK_FSTOP(a);
#else
#define SYMPACK_TIMER_START(a)
#define SYMPACK_TIMER_STOP(a)
#endif



#include "sympack/impl/Environment_impl.hpp"



#endif // _ENVIRONMENT_DECL_HPP_
