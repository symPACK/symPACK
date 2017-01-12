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
#include  "sympack/Environment.hpp"
#include  "sympack/utility.hpp"
#include  "sympack/SuperNode.hpp"

#include <upcxx.h>

  extern "C"
  int symPACK_Init(int *argc=NULL, char ***argv=NULL){
    int retval = 0;
    retval = MPI_Init(argc,argv);
    retval = retval && upcxx::init(argc, argv);
    int mpirank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD,&mpirank);
    std::cout<<"My rank is "<<mpirank<<std::endl;

//    retval = upcxx::init(argc, argv);
//    retval = retval && MPI_Init(argc,argv);
    return retval;
  }

  extern "C"
  int symPACK_Finalize(){
    return upcxx::finalize();
  }

namespace symPACK{



#ifdef _TRACK_MEMORY_
  std::map<char*,size_t> MemoryAllocator::cnt_;
  size_t MemoryAllocator::total_=0;
  size_t MemoryAllocator::hwm_=0;
#endif



// *********************************************************************
// IO
// *********************************************************************
  std::ofstream  statusOFS;

  std::vector<std::ofstream>  statusOFSs;

// *********************************************************************
// Error handling
// *********************************************************************
	// If we are not in RELEASE mode, then implement wrappers for a
	// CallStack
#ifndef _RELEASE_
	std::stack<std::string> callStack;	

	void PushCallStack( std::string s )
	{ callStack.push(s); }

	void PopCallStack()
	{ callStack.pop(); }

	void DumpCallStack()
	{
		std::ostringstream msg;
		while( ! callStack.empty() )
		{
			msg << "Stack[" << callStack.size() << "]: " << callStack.top() << "\n";
			callStack.pop();
		}
		std::cerr << msg.str() << std::endl;
	}

#endif // ifndef _RELEASE_
} // namespace SYMPACK
