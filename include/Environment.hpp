#ifndef _ENVIRONMENT_DECL_HPP_
#define _ENVIRONMENT_DECL_HPP_

// STL libraries
#include <iostream> 
#include <iomanip> 
#include <fstream>
#include <sstream>
#include <unistd.h>

#include <cfloat>
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



// Always use complex data.
#define _USE_COMPLEX_
#include "datatypes.hpp"













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

#include "LogFile.hpp"

extern LogFile * logfileptr;


/***********************************************************************
 *  Data types and constants
 **********************************************************************/

namespace LIBCHOLESKY{

// Basic data types

#ifndef Add_
#define FORTRAN(name) name
#define BLAS(name) name
#define LAPACK(name) name
#else
#define FORTRAN(name) name##_
#define BLAS(name) name##_
#define LAPACK(name) name##_
#endif

// IO
extern  std::ofstream  statusOFS;

// *********************************************************************
// Define constants
// *********************************************************************
// Commonly used
const char UPPER = 'U';
const char LOWER = 'L';


} // namespace LIBCHOLESKY

/***********************************************************************
 *  Error handling
 **********************************************************************/

namespace LIBCHOLESKY{





  extern Int iam,np;

  inline void gdb_lock(){
    int lock = 1;
    while (lock == 1){ }
  }










#ifndef _RELEASE_
void PushCallStack( std::string s );
void PopCallStack();
void DumpCallStack();
#endif // ifndef _RELEASE_

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

} // namespace LIBCHOLESKY


// *********************************************************************
// Profiling functions 
// *********************************************************************


#ifdef USE_TAU
#define PROFILING_ON
#include "TAU.h"
#elif defined (PROFILE) || defined(PMPI)
#define TAU
#include "timer.hpp"
#endif

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


#include "Environment_impl.hpp"


#endif // _ENVIRONMENT_DECL_HPP_
