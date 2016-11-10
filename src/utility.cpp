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
/// @file utility.cpp
/// @brief Implementation of the non-templated utility subroutines.
/// @author Lin Lin
/// @date 2012-09-20
#include "sympack/Environment.hpp"
#include "sympack/utility.hpp"

#include <typeinfo>

using namespace std;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::cerr;

namespace symPACK{


//Int iam,np;

// *********************************************************************
// IO functions
// TODO Move this to utility.hpp and make them inline functions
// *********************************************************************
//---------------------------------------------------------
Int SeparateRead(std::string name, std::istringstream& is)
{
  MPI_Barrier(MPI_COMM_WORLD);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  char filename[100];
  sprintf(filename, "%s_%d_%d", name.c_str(), mpirank, mpisize);  //cerr<<filename<<std::endl;
  std::ifstream fin(filename);
	if( !fin.good() ){
		throw std::logic_error( "File cannot be openeded!" );
	}
 
 	is.str( std::string(std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>()) );
  fin.close();
  //
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}

//---------------------------------------------------------
Int SeparateWrite(std::string name, std::ostringstream& os)
{
   MPI_Barrier(MPI_COMM_WORLD);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  char filename[100];
  sprintf(filename, "%s_%d_%d", name.c_str(), mpirank, mpisize);
  ofstream fout(filename);
	if( !fout.good() ){
		throw std::logic_error( "File cannot be openeded!" );
	}
  fout<<os.str();
  fout.close();
  //
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}


//---------------------------------------------------------
Int SharedWrite(std::string name, std::ostringstream& os)
{
  MPI_Barrier(MPI_COMM_WORLD);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  if(mpirank==0) {
    ofstream fout(name.c_str());
		if( !fout.good() ){
			throw std::logic_error( "File cannot be openeded!" );
		}
    fout<<os.str();
    fout.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}


//---------------------------------------------------------
Int SeparateWriteAscii(std::string name, std::ostringstream& os)
{
  MPI_Barrier(MPI_COMM_WORLD);
/////  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
/////  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
/////  //
/////  char filename[100];
/////  sprintf(filename, "%s_%d_%d", name.c_str(), mpirank, mpisize);
/////  ofstream fout(filename, ios::trunc);
/////	if( !fout.good() ){
/////		throw std::logic_error( "File cannot be opened!" );
/////	}
/////  fout<<os;
/////  fout.close();
  //
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}


//---------------------------------------------------------



void SetValue( std::vector<char>& vec, bool val ){
  fill(vec.begin(),vec.end(),val);
}




  void PtrSum( void *in, void *inout, int *len, MPI_Datatype *dptr ) 
  { 
    int i; 

    Ptr * pinout = (Ptr*)inout;
    Ptr * pin = (Ptr*)in;
#pragma unroll
    for (i=0; i< *len; ++i) { 
      pinout[i] += pin[i];
    } 
  } 


  void PtrMax( void *in, void *inout, int *len, MPI_Datatype *dptr ) 
  { 
    int i; 

    Ptr * pinout = (Ptr*)inout;
    Ptr * pin = (Ptr*)in;
#pragma unroll
    for (i=0; i< *len; ++i) { 
      pinout[i] = std::max(pinout[i], pin[i]);
    } 
  } 




  bool is_complex_type<std::complex<float> >::value = true;
  bool is_complex_type<std::complex<double> >::value = true;





}  // namespace SYMPACK
