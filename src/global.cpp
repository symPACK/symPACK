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

#ifdef NEW_UPCXX
   #include  "sympack/symPACKMatrix2D.hpp"
#endif

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



bool libMPIInit = false;
bool libUPCXXInit = false;


namespace symPACK{
  int mpi_already_init = 0;
  MPI_Comm world_comm = MPI_COMM_NULL;
}

extern "C"
int symPACK_Rank(int * rank){
  int retval = 0;
#ifdef NEW_UPCXX
  *rank = upcxx::rank_me();
#else
  MPI_Initialized(&retval);
  if(retval!=0 || libMPIInit){
    MPI_Comm_rank(MPI_COMM_WORLD,rank);
  }
  else{
#ifdef NEW_UPCXX
    if( libUPCXXInit){
      *rank = upcxx::rank_me();
    }
#else
    if(upcxx::is_init() || libUPCXXInit){
      *rank = upcxx::myrank();
    }
#endif
    else{
      *rank = -1;
      return -1;
    }
  }
#endif
  return 0;
}




extern "C"
int symPACK_Init(int *argc, char ***argv){
#ifdef NEW_UPCXX
  int retval = 1;
  // init UPC++
  if ( libUPCXXInit ) symPACK::gdb_lock();
  upcxx::init();
  libUPCXXInit = true;
  
  // init MPI, if necessary
  //int mpi_already_init;
  MPI_Initialized(&symPACK::mpi_already_init);
  if (!symPACK::mpi_already_init) MPI_Init(argc, argv);

  MPI_Comm_split(MPI_COMM_WORLD, 0, upcxx::rank_me(), &symPACK::world_comm);

  /* program goes here */
  //int mpi_rankme, mpi_rankn;
  //MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rankme);
  //MPI_Comm_size(MPI_COMM_WORLD, &mpi_rankn);
  
  //char hostname[MPI_MAX_PROCESSOR_NAME];
  //int junk;
  //MPI_Get_processor_name(hostname, &junk);
  //std::cout << "Hello world from UPC++ rank " << upcxx::rank_me() << "/" << upcxx::rank_n()
  //  << ", MPI rank " << mpi_rankme << "/" << mpi_rankn << " : "
  //  << hostname << std::endl;
 
    //
    //
    //

#else
  int retval = 1;

//    upcxx::init(argc,argv);
//
//    int mpiinit = 0;
//    MPI_Initialized(&mpiinit);
//    if(mpiinit==0){
//      if(MPI_Init(argc,argv)!=MPI_SUCCESS){
//        symPACK::gdb_lock();
//      }
//    }
//  return 1;

  if(!libUPCXXInit){
    //char *orig_pmi_gni_cookie = getenv("PMI_GNI_COOKIE");
    //if (orig_pmi_gni_cookie) {
    //  char *new_pmi_gni_cookie = (char *)malloc(32);
    //  sprintf(new_pmi_gni_cookie, "PMI_GNI_COOKIE=%d",
    //      1+atoi(orig_pmi_gni_cookie));
    //  putenv(new_pmi_gni_cookie);
    //}


#ifdef NEW_UPCXX
      upcxx::init();
      libUPCXXInit = true;

#else
    if(!upcxx::is_init()){
//std::cerr<<"upcxx init"<<std::endl;
      upcxx::init(argc,argv);
//std::cerr<<"past upcxx init"<<std::endl;
      libUPCXXInit = true;
    }
#endif
  }

  if(!libMPIInit){
    symPACK::mpi_already_init = 0;
    MPI_Initialized(&symPACK::mpi_already_init);
    if(symPACK::mpi_already_init==0){
      if(MPI_Init(argc,argv)==MPI_SUCCESS){
        libMPIInit = true;
        retval = retval && 1;
      }
//std::cerr<<"past MPI_Initialize"<<std::endl;
    }
  }

  MPI_Comm_split(MPI_COMM_WORLD, 0, upcxx::rank_me(), &symPACK::world_comm);

    int mpiinit = 0;
    MPI_Initialized(&mpiinit);
    assert(mpiinit==1);
//    assert(upcxx::is_init());

#endif
  return retval;
}


extern "C"
int symPACK_Finalize(){
  int retval = 1;
  int rank = 0;
  symPACK_Rank(&rank);

  MPI_Comm_free(&symPACK::world_comm);

  if (!symPACK::mpi_already_init) MPI_Finalize();

  if(libUPCXXInit){
    //if(rank==0){std::cerr<<"upcxx finalize"<<std::endl;}
    upcxx::finalize();
    //if(rank==0){std::cerr<<"past upcxx finalize"<<std::endl;}
    libUPCXXInit = false;
    libMPIInit = false;
  }

  
  ////if(libMPIInit)
  ////{
  ////  int mpiflag = 0;
  ////  if(rank==0){std::cerr<<"MPI_Finalized"<<std::endl;};
  ////  MPI_Finalized(&mpiflag);
  ////  if(mpiflag==0){
  ////    if(rank==0){std::cerr<<"MPI_Finalize"<<std::endl;};
  ////    retval = retval && MPI_SUCCESS == MPI_Finalize();
  ////    if(rank==0){std::cerr<<"past MPI_Finalize"<<std::endl;};
  ////    libMPIInit = false;
  ////  }
  ////}

  ////int mpiinit = 0;
  ////MPI_Finalized(&mpiinit);
  ////assert(mpiinit==1);
 // assert(!upcxx::is_init());
  return retval;
}






namespace symPACK{

  int symPACKMatrixBase::last_id = 0;

#ifdef NEW_UPCXX
  //std::map<int, std::deque<incoming_data_t>  > g_sp_handle_incoming;

  std::map<int, symPACKMatrixBase *  > g_sp_handle_to_matrix;
#endif

#ifdef _TRACK_MEMORY_
  std::map<char*,size_t> MemoryAllocator::cnt_;
  size_t MemoryAllocator::total_=0;
  size_t MemoryAllocator::hwm_=0;
#endif

} // namespace SYMPACK



namespace symPACK{
  namespace Multithreading{ 
    int NumThread = 1;
  } // namespace Multithreading
} // namespace SYMPACK


