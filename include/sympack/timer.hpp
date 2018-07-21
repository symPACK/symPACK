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
#ifndef __TIMER_H__
#define __TIMER_H__

#include <mpi.h>
#include <sys/time.h>
#include "sympack/Environment.hpp"
#include "sympack/LogFile.hpp"
#include <iostream>
#include <string>


namespace symPACK{



class symPACK_timer{
  public:
    char const * timer_name;
    int64_t index;
    int64_t exited;
    int64_t original;

    vector<int64_t> arr_index;
    vector<int64_t> arr_exited;
    vector<int64_t> arr_original;
  
  public:
    symPACK_timer(char const * name);
    ~symPACK_timer();
    void stop();
    void start();
    void exit();
    
};

void symPACK_set_main_args(int64_t argc, char * const * argv);
void symPACK_set_context(MPI_Comm ctxt);

//extern Int iam;

class symPACK_scope_timer{
  std::string name;
  public:

  symPACK_scope_timer(const char * pname){
    name = pname;
    symPACK_timer t(name.c_str()); 
    t.start();
  }
  ~symPACK_scope_timer(){
    symPACK_timer t(name.c_str());
    t.stop();
  }
};

#ifdef SPROFILE
#define SYMPACK_FSTART(ARG)                                           \
  do { symPACK::symPACK_timer t(#ARG);/* logfileptr->OFS()<<"DEBUG START "<<#ARG<<std::endl;*/ t.start(); } while (0)

#define SYMPACK_FSTOP(ARG)                                            \
  do { symPACK::symPACK_timer t(#ARG);/* logfileptr->OFS()<<"DEBUG STOP "<<#ARG<<std::endl;*/ t.stop(); } while (0)


#define SYMPACK_SPROFILE_INIT(argc, argv)                              \
  symPACK::symPACK_set_main_args(argc, argv);

#define SYMPACK_PROFILE_SET_CONTEXT(ARG)                              \
  if (ARG==0) symPACK::symPACK_set_context(symPACK::world_comm);                    \
  else symPACK::symPACK_set_context((MPI_Comm)ARG);

#define SYMPACK_TIMER_START(a) SYMPACK_FSTART(a);
#define SYMPACK_TIMER_STOP(a) SYMPACK_FSTOP(a);
#define SYMPACK_TIMER_SPECIAL_START(a) SYMPACK_FSTART(a);
#define SYMPACK_TIMER_SPECIAL_STOP(a) SYMPACK_FSTOP(a);
#define scope_timer(b,a) symPACK::symPACK_scope_timer b(#a)
#define scope_timer_special(b,a) symPACK::symPACK_scope_timer b(#a)

#else
#define SYMPACK_TIMER_START(a)
#define SYMPACK_TIMER_STOP(a)
#define SYMPACK_TIMER_SPECIAL_START(a)
#define SYMPACK_TIMER_SPECIAL_STOP(a)
#define scope_timer(b,a)
#define scope_timer_special(b,a)
#endif



}

#endif //__TIMER_H__

