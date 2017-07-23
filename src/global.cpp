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

//#include  "sympack/symPACKMatrix.hpp"

//namespace symPACK{
//template class symPACKMatrix<float>;
//template class symPACKMatrix<double>;
//template class symPACKMatrix<std::complex<float> >;
//template class symPACKMatrix<std::complex<double> >;
//template class supernodalTaskGraph<FBTask>;
//template class supernodalTaskGraph<CompTask>;
//}
//namespace symPACK{
//
//extern template class symPACKMatrix<float>;
//extern template class symPACKMatrix<double>;
//extern template class symPACKMatrix<std::complex<float> >;
//extern template class symPACKMatrix<std::complex<double> >;
//}

#include <upcxx.h>


bool libMPIInit = false;
bool libUPCXXInit = false;

extern "C"
int symPACK_Init(int *argc=NULL, char ***argv=NULL){
  int retval = 0;
  MPI_Initialized(&retval);
  if(retval==0){
    retval = MPI_Init(argc,argv)==MPI_SUCCESS;
    if(retval){
      libMPIInit = true;
    }
    else{
      abort();
    }
  }

  char *orig_pmi_gni_cookie = getenv("PMI_GNI_COOKIE");
  if (orig_pmi_gni_cookie) {
    char *new_pmi_gni_cookie = (char *)malloc(32);
    sprintf(new_pmi_gni_cookie, "PMI_GNI_COOKIE=%d",
        1+atoi(orig_pmi_gni_cookie));
    putenv(new_pmi_gni_cookie);
  }


  if(!upcxx::is_init()){
    retval = retval && upcxx::init(NULL,NULL);
    assert(upcxx::is_init());
    libUPCXXInit = true;
  }

  //    retval = upcxx::init(argc, argv);
  //
  //char *orig_pmi_gni_cookie = getenv("PMI_GNI_COOKIE");
  //if (orig_pmi_gni_cookie) {
  //     char *new_pmi_gni_cookie = (char *)malloc(32);
  //     sprintf(new_pmi_gni_cookie, "PMI_GNI_COOKIE=%d",
  //                 1+atoi(orig_pmi_gni_cookie));
  //     putenv(new_pmi_gni_cookie);
  //}
  //
  //
  //    retval = retval && MPI_Init(argc,argv);
  return retval;
}

extern "C"
int symPACK_Finalize(){
  if(libUPCXXInit){
    return upcxx::finalize();
  }
  else{
    return 0;
  }
}


extern "C"
int symPACK_Rank(int * rank){
  int retval = 0;
  MPI_Initialized(&retval);
  if(retval!=0 || libMPIInit){
    MPI_Comm_rank(MPI_COMM_WORLD,rank);
  }
  else{
    if(upcxx::is_init() || libUPCXXInit){
      *rank = upcxx::myrank();
    }
    else{
      *rank = -1;
      return -1;
    }
  }
  return 0;
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

} // namespace SYMPACK



namespace symPACK{
  namespace Multithreading{ 
    int NumThread = 1;
  } // namespace Multithreading





} // namespace SYMPACK


