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
#ifndef _ENVIRONMENT_IMPL_HPP_
#define _ENVIRONMENT_IMPL_HPP_


// *********************************************************************
// Global utility functions 
// These utility functions do not depend on local definitions
// *********************************************************************
namespace symPACK{
inline Int 
	iround(Real a){ 
		Int b = 0;
		if(a>0) b = (a-Int(a)<0.5)?Int(a):(Int(a)+1);
		else b = (Int(a)-a<0.5)?Int(a):(Int(a)-1);
		return b; 
	}
inline void OptionsCreate(Int argc, char** argv, std::map<std::string,std::string >& options)
{
	options.clear();
	for(Int k=1; k<argc; k=k+2) {
		options[ std::string(argv[k]) ] = std::string(argv[k+1]);
	}
}
inline void OptionsCreate(Int argc, char** argv, std::map<std::string,std::vector<std::string> >& options)
{
	options.clear();
  Int k =1;
  Int key = 1;
  while(k<argc){
    if(*argv[k]=='-'){
      options[std::string(argv[k])].clear();
      key = k;
    }
    else{
      options[std::string(argv[key])].push_back(std::string(argv[k]));
    }
    k++;
  }
//	for(Int k=1; k<argc; k=k+2) {
//		options[ std::string(argv[k]) ] = std::string(argv[k+1]);
//	}
}

// *********************************************************************
// std::stringstream
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

} // namespace SYMPACK

#endif // _ENVIRONMENT_IMPL_HPP_
