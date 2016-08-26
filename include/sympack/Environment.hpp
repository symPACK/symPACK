#ifndef _ENVIRONMENT_DECL_HPP_
#define _ENVIRONMENT_DECL_HPP_

#include "sympack_definitions.hpp"
#include "sympack_config.hpp"

//debug
#include <sys/types.h>
#include <unistd.h>

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
#include <execinfo.h>
//#include <signal.h>
#include <exception>

// MPI
#include <mpi.h>
#include <upcxx.h>

#include "sympack/datatypes.hpp"

#ifndef _USE_NUMVEC_
namespace SYMPACK{
template <typename T> using vector = std::vector<T>;
}
#else
#include "NumVec.hpp"
namespace SYMPACK{
  template <typename T> using vector = SYMPACK::NumVec<T,size_t>;
}
#endif



// Always use complex data.
#define _USE_COMPLEX_
#include "sympack/Types.hpp"



#include "sympack/LogFile.hpp"










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

namespace SYMPACK{

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

/***********************************************************************
 *  Error handling
 **********************************************************************/

namespace SYMPACK{

  extern upcxx::team * workteam;
  extern Int iam;
  extern Int np;


  inline void gdb_lock(){
      pid_t pid = getpid();
      std::cout<<"P"<<iam<<" is locked, pid is "<<pid<<std::endl;
    volatile int lock = 1;
    while (lock == 1){ }
      std::cout<<"P"<<iam<<" is unlocked"<<std::endl;
  }

  inline void gdb_lock(Int proc){
    if(iam==proc){
      pid_t pid = getpid();
      std::cout<<"P"<<iam<<" is locked, pid is "<<pid<<std::endl;
      volatile int lock = 1;
      while (lock == 1){ }
      std::cout<<"P"<<iam<<" is unlocked"<<std::endl;
    }
  }


  inline void bassert(bool cond){if(!(cond)){gdb_lock();}}



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
			std::cout << symbols[i] << std::endl;
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
#include <sys/time.h>




inline double get_time()
  {
    struct timeval tv;
    gettimeofday(&tv, 0);
    return tv.tv_sec + ((double) tv.tv_usec / 1000000);
  }


#ifdef USE_TAU
#define PROFILING_ON
#include "TAU.h"
#elif defined (PROFILE) || defined(PMPI)
#define TAU
#endif
#include "sympack/timer.hpp"

#define VAL(str) #str
#define TOSTRING(str) VAL(str)


#ifdef USE_TAU 
#define TIMER_START(a) TAU_START(TOSTRING(a));
#define TIMER_STOP(a) TAU_STOP(TOSTRING(a));
#elif defined (PROFILE)
#define TIMER_START(a) TAU_FSTART(a);
#define TIMER_STOP(a) TAU_FSTOP(a);
#else
#define TIMER_START(a)
#define TIMER_STOP(a)
#endif


#include "sympack/Environment_impl.hpp"






#endif // _ENVIRONMENT_DECL_HPP_
