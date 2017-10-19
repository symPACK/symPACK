/*
	 Copyright (c) 2016 The Regents of the University of California,
	 through Lawrence Berkeley National Laboratory.  

   Author: Mathias Jacquelin
	 
   This file is part of symPACK. All rights reserved.

	 Redistribution and use in source and binary forms, with or without
	 modification, are permitted provided that the following conditions are met:

	 (1) Redistributions of source code must retain the above copyright notice, this
	 list of conditions and the following disclaimer.
	 (2) Redistributions in binary form must reproduce the above copyright notice,
	 this list of conditions and the following disclaimer in the documentation
	 and/or other materials provided with the distribution.
	 (3) Neither the name of the University of California, Lawrence Berkeley
	 National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
	 be used to endorse or promote products derived from this software without
	 specific prior written permission.

	 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
	 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
	 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
	 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
	 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

	 You are under no obligation whatsoever to provide any bug fixes, patches, or
	 upgrades to the features, functionality or performance of the source code
	 ("Enhancements") to anyone; however, if you choose to make your Enhancements
	 available either publicly, or directly to Lawrence Berkeley National
	 Laboratory, without imposing a separate written license agreement for such
	 Enhancements, then you hereby grant the following license: a non-exclusive,
	 royalty-free perpetual license to install, use, modify, prepare derivative
	 works, incorporate into other computer software, distribute, and sublicense
	 such enhancements or derivative works thereof, in binary and source code form.
*/
#ifndef _ENVIRONMENT_DECL_HPP_
#define _ENVIRONMENT_DECL_HPP_

//#define _MAP_DEBUG_
#define SP_THREADS

#include "sympack_definitions.hpp"
#include "sympack_config.hpp"

//#define EXPLICIT_PERMUTE


#ifdef NEW_UPCXX
#include <upcxx/upcxx.hpp>
#else
#include <upcxx.h>
#endif

//debug
#include <sys/types.h>
#include <unistd.h>
#include <execinfo.h>
//#include <signal.h>
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



  inline void gdb_lock(){
      pid_t pid = getpid();
#ifdef NEW_UPCXX
int iam = upcxx::rank_me();
#else
int iam = upcxx::myrank();
#endif
      std::cout<<"P"<<iam<<" is locked, pid is "<<pid<<std::endl;
    volatile int lock = 1;
    while (lock == 1){ }
      std::cout<<"P"<<iam<<" is unlocked"<<std::endl;
  }

  inline void gdb_lock(Int proc){
#ifdef NEW_UPCXX
int iam = upcxx::rank_me();
#else
int iam = upcxx::myrank();
#endif
    if(iam==proc){
      pid_t pid = getpid();
      std::cout<<"P"<<iam<<" is locked, pid is "<<pid<<std::endl;
      volatile int lock = 1;
      while (lock == 1){ }
      std::cout<<"P"<<iam<<" is unlocked"<<std::endl;
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
