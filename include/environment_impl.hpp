#ifndef _ENVIRONMENT_IMPL_HPP_
#define _ENVIRONMENT_IMPL_HPP_

#include "Environment.hpp"

// *********************************************************************
// Global utility functions 
// These utility functions do not depend on local definitions
// *********************************************************************
namespace LIBCHOLESKY{
inline Int 
	iround(Real a){ 
		Int b = 0;
		if(a>0) b = (a-Int(a)<0.5)?Int(a):(Int(a)+1);
		else b = (Int(a)-a<0.5)?Int(a):(Int(a)-1);
		return b; 
	}

inline void OptionsCreate(Int argc, char** argv, std::map<std::string,std::string>& options)
{
	options.clear();
	for(Int k=1; k<argc; k=k+2) {
		options[ std::string(argv[k]) ] = std::string(argv[k+1]);
	}
}

// *********************************************************************
// Stringstream
// *********************************************************************

// Size information.
// Like sstm.str().length() but without making the copy
inline Int Size( std::stringstream& sstm ){
	Int length;
	sstm.seekg (0, std::ios::end);
	length = sstm.tellg();
	sstm.seekg (0, std::ios::beg);
	return length;
}

} // namespace LIBCHOLESKY

#endif // _ENVIRONMENT_IMPL_HPP_
