#ifndef _DATATYPES_HPP_
#define _DATATYPES_HPP_


#include  "sympack/Environment.hpp"
#include <vector>

namespace symPACK {

  //typedef    long long int              LongInt;
  typedef    long long int              LongInt;
  typedef    int                   Int;
  //typedef    int64_t               Int;
  typedef    uint32_t              Idx32;
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
    const T ZERO(){ return static_cast<T>(0.0);}
  template<typename T>
    const T ONE(){ return static_cast<T>(1.0);}
  template<typename T>
    const T MINUS_ONE(){ return static_cast<T>(-1.0);}



  typedef Idx64 Ptr;
  typedef Idx32 Idx;
  //typedef Idx64 Idx;



}

#endif // ifdef _DATATYPES_HPP_
