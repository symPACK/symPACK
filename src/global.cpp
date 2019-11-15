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

   #include  "sympack/symPACKMatrix2D.hpp"



bool libMPIInit = false;
bool libUPCXXInit = false;


namespace symPACK{
  int mpi_already_init = 0;
  MPI_Comm world_comm = MPI_COMM_NULL;
  upcxx::persona_scope * master_scope = nullptr;


  void liberate_master_scope(){
    delete master_scope;
    master_scope = nullptr;
  }

  void capture_master_scope() {
    if ( master_scope == nullptr ) {
      master_scope = new upcxx::persona_scope(upcxx::master_persona());
    }
  }
}

extern "C"
int symPACK_Rank(int * rank){
  int retval = 0;
  *rank = upcxx::rank_me();
    return 0;
}




extern "C"
int symPACK_Init(int *argc, char ***argv){
  int retval = 1;
  // init UPC++
  if ( libUPCXXInit ) symPACK::gdb_lock();
  upcxx::init();
  upcxx::liberate_master_persona();
  symPACK::capture_master_scope();
  libUPCXXInit = true;
  
  // init MPI, if necessary
  MPI_Initialized(&symPACK::mpi_already_init);
  if (!symPACK::mpi_already_init) MPI_Init(argc, argv);

  MPI_Comm_split(MPI_COMM_WORLD, 0, upcxx::rank_me(), &symPACK::world_comm);
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
    symPACK::capture_master_scope();

    upcxx::finalize();
    symPACK::liberate_master_scope();

    libUPCXXInit = false;
    libMPIInit = false;
  }

  return retval;
}






namespace symPACK{

  int symPACKMatrixBase::last_id = 0;

  std::map<int, symPACKMatrixBase *  > g_sp_handle_to_matrix;

} // namespace SYMPACK



namespace symPACK{
  namespace Multithreading{ 
    int NumThread = 1;
  } // namespace Multithreading
} // namespace SYMPACK


