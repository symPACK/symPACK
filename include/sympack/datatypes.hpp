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
/// @file datatypes.hpp
/// @brief Abstract datatypes
/// @author Mathias Jacquelin
/// @date 2013-08-31
#ifndef _DATATYPES_HPP_
#define _DATATYPES_HPP_


#include  "sympack/Environment.hpp"
#include <vector>

namespace symPACK {

  //typedef    long long int              LongInt;
  typedef    long long int              LongInt;
  typedef    int                   Int;
  typedef    uint32_t              Idx32;
  typedef    uint64_t              Idx64;
  typedef    uint64_t              Idx64;
  typedef    double                Real;
  typedef    std::complex<double>  Complex; // Must use elemental form of complex
#ifdef _USE_COMPLEX_
  typedef    std::complex<double>  Scalar;  // Must use elemental form of complex
#else
  typedef    double                Scalar;
#endif

  // *********************************************************************
  // Define constants
  // *********************************************************************
  // Commonly used
  const Int I_ZERO = 0;
  const Int I_ONE  = 1;
  const Int I_MINUS_ONE  = -1;
  const Real D_ZERO = 0.0;
  const Real D_ONE  = 1.0;
  const Real D_MINUS_ONE  = -1.0;
  const Complex Z_ZERO = Complex(0.0, 0.0);
  const Complex Z_ONE  = Complex(1.0, 0.0);
  const Complex Z_MINUS_ONE  = Complex(-1.0, 0.0);
  const Complex Z_I    = Complex(0.0, 1.0);
  const Complex Z_MINUS_I    = Complex(0.0, -1.0);
  const Scalar SCALAR_ZERO    = static_cast<Scalar>(0.0);
  const Scalar SCALAR_ONE     = static_cast<Scalar>(1.0);
  const Scalar SCALAR_MINUS_ONE = static_cast<Scalar>(-1.0);

  template<typename T>
    const T ZERO(){ return static_cast<T>(0.0);};
  template<typename T>
    const T ONE(){ return static_cast<T>(1.0);};
  template<typename T>
    const T MINUS_ONE(){ return static_cast<T>(-1.0);};



  typedef Idx64 Ptr;
  typedef Idx32 Idx;



}

#endif // ifdef _DATATYPES_HPP_
