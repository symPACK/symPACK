#ifndef _UTILITY_HPP_ 
#define _UTILITY_HPP_





#include  <mpi.h>
#include  <stdlib.h>
#include  <numeric>

#include <exception>
#include <stdexcept>
#include <sstream>
#include <mutex>
#include <type_traits>

#include  "sympack/Environment.hpp"
#include  "sympack/Types.hpp"
#include  "sympack/DistSparseMatrix.hpp"
#include  "sympack/ETree.hpp"

#include <atomic>


#include <sys/types.h>  /* Primitive System Data Types */ 
#include <stdio.h>      /* Input/Output */
#include <stdlib.h>     /* General Utilities */
#include <time.h>	/* time and difftime */
#include <string.h>	/* string stuff */

#include <sys/resource.h> /* getrusage */




namespace symPACK{

  class SharedNode{
    public:
      MPI_Comm shmcomm;
      int shmrank;
      int shmsize;
      SharedNode(){
        shmrank=0;
        shmsize=1;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
        MPI_Comm_rank(shmcomm, &shmrank);
        MPI_Comm_size(shmcomm, &shmsize);
        if(shmsize>1){
          shmsize--;
        }
      }
      ~SharedNode(){
        MPI_Comm_free(&shmcomm);
      }
  };

  namespace Multithreading{

    class SpinLock
    {
      public:
        void lock()
        {
          while(lck.test_and_set(std::memory_order_acquire))
          {}
        }

        void unlock()
        {
          lck.clear(std::memory_order_release);
        }

      private:
        std::atomic_flag lck = ATOMIC_FLAG_INIT;
    };


  }
}


namespace symPACK{
  // Function to return gcd of a and b 
  template <typename F> 
   inline F gcd(F a, F b) { 
       if (a == 0) 
         return b; 
       return gcd(b % a, a); 
   } 
                     
  // Function to find gcd of array of numbers 
  template <typename F, typename S> 
   inline F findGCD(F * arr, S n) { 
         typename std::remove_cv<F>::type result = arr[0]; 
         for (S i = 1; i < n; i++) 
             result = gcd(arr[i], result); 
                                       
         return result; 
   } 


  //exact cost of a mxn panel
#define CHOLESKY_COST(m,n)  ((n)*pow((m),2.0) + 2*(n)*(m)-pow((n),3.0)/3.0 - 3.0*pow((n),2.0)/2.0 - (n)/6.0 -1.0)


#define FMULS_TRMM_2(__m, __n) (0.5 * (double)(__n) * (double)(__m) * ((double)(__m)+1.))
#define FADDS_TRMM_2(__m, __n) (0.5 * (double)(__n) * (double)(__m) * ((double)(__m)-1.))


#define FMULS_TRMM(__side, __m, __n) ( ( (__side) == 'L' ) ? FMULS_TRMM_2((__m), (__n)) : FMULS_TRMM_2((__n), (__m)) )
#define FADDS_TRMM(__side, __m, __n) ( ( (__side) == 'L' ) ? FADDS_TRMM_2((__m), (__n)) : FADDS_TRMM_2((__n), (__m)) )

#define FMULS_TRSM FMULS_TRMM
#define FADDS_TRSM FMULS_TRMM

#define FLOPS_ZTRSM(__side, __m, __n) (6. * FMULS_TRSM(__side, (__m), (__n)) + 2.0 * FADDS_TRSM(__side, (__m), (__n)) )
#define FLOPS_CTRSM(__side, __m, __n) (6. * FMULS_TRSM(__side, (__m), (__n)) + 2.0 * FADDS_TRSM(__side, (__m), (__n)) )
#define FLOPS_DTRSM(__side, __m, __n) (     FMULS_TRSM(__side, (__m), (__n)) +       FADDS_TRSM(__side, (__m), (__n)) )
#define FLOPS_STRSM(__side, __m, __n) (     FMULS_TRSM(__side, (__m), (__n)) +       FADDS_TRSM(__side, (__m), (__n)) )


#define FMULS_POTRF(__n) ((double)(__n) * (((1. / 6.) * (double)(__n) + 0.5) * (double)(__n) + (1. / 3.)))
#define FADDS_POTRF(__n) ((double)(__n) * (((1. / 6.) * (double)(__n)      ) * (double)(__n) - (1. / 6.)))


#define FLOPS_ZPOTRF(__n) (6. * FMULS_POTRF((__n)) + 2.0 * FADDS_POTRF((__n)) )
#define FLOPS_CPOTRF(__n) (6. * FMULS_POTRF((__n)) + 2.0 * FADDS_POTRF((__n)) )
#define FLOPS_DPOTRF(__n) (     FMULS_POTRF((__n)) +       FADDS_POTRF((__n)) )
#define FLOPS_SPOTRF(__n) (     FMULS_POTRF((__n)) +       FADDS_POTRF((__n)) )


#define FMULS_GEMM(__m, __n, __k) ((double)(__m) * (double)(__n) * (double)(__k))
#define FADDS_GEMM(__m, __n, __k) ((double)(__m) * (double)(__n) * (double)(__k))


#define FLOPS_ZGEMM(__m, __n, __k) (6. * FMULS_GEMM((__m), (__n), (__k)) + 2.0 * FADDS_GEMM((__m), (__n), (__k)) )
#define FLOPS_CGEMM(__m, __n, __k) (6. * FMULS_GEMM((__m), (__n), (__k)) + 2.0 * FADDS_GEMM((__m), (__n), (__k)) )
#define FLOPS_DGEMM(__m, __n, __k) (     FMULS_GEMM((__m), (__n), (__k)) +       FADDS_GEMM((__m), (__n), (__k)) )
#define FLOPS_SGEMM(__m, __n, __k) (     FMULS_GEMM((__m), (__n), (__k)) +       FADDS_GEMM((__m), (__n), (__k)) )





  template <typename F> void SetValue( std::vector<F>& vec, F val ){
    fill(vec.begin(),vec.end(),val);
  }

  void SetValue( std::vector<char>& vec, bool val );


  // *********************************************************************
  // Define constants
  // *********************************************************************

  // Write format control parameters 
  const int LENGTH_VAR_NAME = 8;
  const int LENGTH_DBL_DATA = 16;
  const int LENGTH_INT_DATA = 5;
  const int LENGTH_VAR_UNIT = 6;
  const int LENGTH_DBL_PREC = 8;
  const int LENGTH_VAR_DATA = 16;



  // *********************************************************************
  // Formatted output stream
  // *********************************************************************

  // Bool is NOT defined due to ambiguity with Int

  inline Int PrintBlock(std::ostream &os, const std::string name){

    os << std::endl<< "*********************************************************************" << std::endl;
    os << name << std::endl;
    os << "*********************************************************************" << std::endl << std::endl;
    return 0;
  }

  // std::string
  inline Int Print(std::ostream &os, const std::string name) {
    os << std::setiosflags(std::ios::left) << name << std::endl;
    return 0;
  }

  inline Int Print(std::ostream &os, const char* name) {
    os << std::setiosflags(std::ios::left) << std::string(name) << std::endl;
    return 0;
  }

  inline Int Print(std::ostream &os, const std::string name, std::string val) {
    os << std::setiosflags(std::ios::left) 
      << std::setw(LENGTH_VAR_NAME) << name
      << std::setw(LENGTH_VAR_DATA) << val
      << std::endl;
    return 0;
  }

  inline Int Print(std::ostream &os, const std::string name, const char* val) {
    os << std::setiosflags(std::ios::left) 
      << std::setw(LENGTH_VAR_NAME) << name
      << std::setw(LENGTH_VAR_DATA) << std::string(val)
      << std::endl;
    return 0;
  }


  // Real

  // one real number

  inline Int Print(std::ostream &os, const std::string name, Real val) {
    os << std::setiosflags(std::ios::left) 
      << std::setw(LENGTH_VAR_NAME) << name
      << std::setiosflags(std::ios::scientific)
      << std::setiosflags(std::ios::showpos)
      << std::setw(LENGTH_DBL_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val
      << std::resetiosflags(std::ios::scientific)
      << std::resetiosflags(std::ios::showpos)
      << std::endl;
    return 0;
  }

  inline Int Print(std::ostream &os, const char* name, Real val) {
    os << std::setiosflags(std::ios::left) 
      << std::setw(LENGTH_VAR_NAME) << std::string(name)
      << std::setiosflags(std::ios::scientific)
      << std::setiosflags(std::ios::showpos)
      << std::setw(LENGTH_DBL_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val
      << std::resetiosflags(std::ios::scientific)
      << std::resetiosflags(std::ios::showpos)
      << std::endl;
    return 0;
  }


  inline Int Print(std::ostream &os, const std::string name, Real val, const std::string unit) {
    os << std::setiosflags(std::ios::left) 
      << std::setw(LENGTH_VAR_NAME) << name
      << std::setiosflags(std::ios::scientific)
      << std::setiosflags(std::ios::showpos)
      << std::setw(LENGTH_DBL_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val
      << std::resetiosflags(std::ios::scientific)
      << std::resetiosflags(std::ios::showpos)
      << std::setw(LENGTH_VAR_UNIT) << unit 
      << std::endl;
    return 0;
  }

  inline Int Print(std::ostream &os, const char *name, Real val, const char *unit) {
    os << std::setiosflags(std::ios::left) 
      << std::setw(LENGTH_VAR_NAME) << std::string(name)
      << std::setiosflags(std::ios::scientific)
      << std::setiosflags(std::ios::showpos)
      << std::setw(LENGTH_DBL_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val
      << std::resetiosflags(std::ios::scientific)
      << std::resetiosflags(std::ios::showpos)
      << std::setw(LENGTH_VAR_UNIT) << std::string(unit) 
      << std::endl;
    return 0;
  }

  // Two real numbers
  inline Int Print(std::ostream &os, const std::string name1, Real val1, const std::string unit1,
      const std::string name2, Real val2, const std::string unit2) {
    os << std::setiosflags(std::ios::left)
      << std::setw(LENGTH_VAR_NAME) << name1
      << std::setiosflags(std::ios::scientific)
      << std::setiosflags(std::ios::showpos)
      << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val1
      << std::setw(LENGTH_VAR_UNIT) << unit1 
      << std::setw(LENGTH_VAR_NAME) << name2
      << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val2
      << std::resetiosflags(std::ios::scientific)
      << std::resetiosflags(std::ios::showpos)
      << std::setw(LENGTH_VAR_UNIT) << unit2 
      << std::endl;
    return 0;
  }

  inline Int Print(std::ostream &os, const char *name1, Real val1, const char *unit1,
      char *name2, Real val2, char *unit2) {
    os << std::setiosflags(std::ios::left)
      << std::setw(LENGTH_VAR_NAME) << std::string(name1)
      << std::setiosflags(std::ios::scientific)
      << std::setiosflags(std::ios::showpos)
      << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val1
      << std::setw(LENGTH_VAR_UNIT) << std::string(unit1) 
      << std::setw(LENGTH_VAR_NAME) << std::string(name2)
      << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val2
      << std::resetiosflags(std::ios::scientific)
      << std::resetiosflags(std::ios::showpos)
      << std::setw(LENGTH_VAR_UNIT) << std::string(unit2) 
      << std::endl;
    return 0;
  }

  // Int and Real
  inline Int Print(std::ostream &os, const std::string name1, Int val1, const std::string unit1,
      const std::string name2, Real val2, const std::string unit2) {
    os << std::setiosflags(std::ios::left)
      << std::setw(LENGTH_VAR_NAME) << name1
      << std::setw(LENGTH_INT_DATA) << val1
      << std::setw(LENGTH_VAR_UNIT) << unit1 
      << std::setw(LENGTH_VAR_NAME) << name2
      << std::setiosflags(std::ios::scientific)
      << std::setiosflags(std::ios::showpos)
      << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val2
      << std::resetiosflags(std::ios::scientific)
      << std::resetiosflags(std::ios::showpos)
      << std::setw(LENGTH_VAR_UNIT) << unit2 
      << std::endl;
    return 0;
  }

  inline Int Print(std::ostream &os, const char *name1, Int val1, const char *unit1,
      char *name2, Real val2, char *unit2) {
    os << std::setiosflags(std::ios::left)
      << std::setw(LENGTH_VAR_NAME) << std::string(name1)
      << std::setw(LENGTH_INT_DATA) << val1
      << std::setw(LENGTH_VAR_UNIT) << std::string(unit1) 
      << std::setw(LENGTH_VAR_NAME) << std::string(name2)
      << std::setiosflags(std::ios::scientific)
      << std::setiosflags(std::ios::showpos)
      << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val2
      << std::resetiosflags(std::ios::scientific)
      << std::resetiosflags(std::ios::showpos)
      << std::setw(LENGTH_VAR_UNIT) << std::string(unit2) 
      << std::endl;
    return 0;
  }

  inline Int Print(std::ostream &os, 
      const char *name1, Int val1, 
      const char *name2, Real val2, 
      char *name3, Real val3) {
    os << std::setiosflags(std::ios::left)
      << std::setw(LENGTH_VAR_NAME) << std::string(name1)
      << std::setw(LENGTH_INT_DATA) << val1 
      << std::setiosflags(std::ios::scientific)
      << std::setiosflags(std::ios::showpos)
      << std::setw(LENGTH_VAR_NAME) << std::string(name2)
      << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val2
      << std::setw(LENGTH_VAR_NAME) << std::string(name3) 
      << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val3
      << std::resetiosflags(std::ios::scientific)
      << std::resetiosflags(std::ios::showpos)
      << std::endl;
    return 0;
  }


  inline Int Print(std::ostream &os, 
      const char *name1, Int val1, 
      const char *name2, Real val2, 
      const char *name3, Real val3,
      const char *name4, Real val4 ) {
    os << std::setiosflags(std::ios::left)
      << std::setw(LENGTH_VAR_NAME) << std::string(name1)
      << std::setw(LENGTH_INT_DATA) << val1 
      << std::setiosflags(std::ios::scientific)
      << std::setiosflags(std::ios::showpos)
      << std::setw(LENGTH_VAR_NAME) << std::string(name2)
      << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val2
      << std::setw(LENGTH_VAR_NAME) << std::string(name3) 
      << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val3
      << std::setw(LENGTH_VAR_NAME) << std::string(name4) 
      << std::setw(LENGTH_VAR_DATA) << std::setprecision(LENGTH_DBL_PREC)<< val4
      << std::resetiosflags(std::ios::scientific)
      << std::resetiosflags(std::ios::showpos)
      << std::endl;
    return 0;
  }

  // Int

  // one integer number
  inline Int Print(std::ostream &os, std::string name, Int val) {
    os << std::setiosflags(std::ios::left)
      << std::setw(LENGTH_VAR_NAME) << name
      << std::setw(LENGTH_VAR_DATA) << val
      << std::endl;
    return 0;
  }


  inline Int Print(std::ostream &os, const char *name, Int val) {
    os << std::setiosflags(std::ios::left)
      << std::setw(LENGTH_VAR_NAME) << std::string(name)
      << std::setw(LENGTH_VAR_DATA) << val
      << std::endl;
    return 0;
  }


  inline Int Print(std::ostream &os, const std::string name, Int val, const std::string unit) {
    os << std::setiosflags(std::ios::left)
      << std::setw(LENGTH_VAR_NAME) << name
      << std::setw(LENGTH_VAR_DATA) << val
      << std::setw(LENGTH_VAR_UNIT) << unit 
      << std::endl;
    return 0;
  }


  inline Int Print(std::ostream &os, const char* name, Int val, const std::string unit) {
    os << std::setiosflags(std::ios::left)
      << std::setw(LENGTH_VAR_NAME) << std::string(name)
      << std::setw(LENGTH_VAR_DATA) << val
      << std::setw(LENGTH_VAR_UNIT) << unit 
      << std::endl;
    return 0;
  }



  // two integer numbers
  inline Int Print(std::ostream &os, const std::string name1, Int val1, const std::string unit1,
      const std::string name2, Int val2, const std::string unit2) {
    os << std::setiosflags(std::ios::left)
      << std::setw(LENGTH_VAR_NAME) << name1
      << std::setw(LENGTH_VAR_DATA) << val1
      << std::setw(LENGTH_VAR_UNIT) << unit1 
      << std::setw(LENGTH_VAR_NAME) << name2
      << std::setw(LENGTH_VAR_DATA) << val2
      << std::setw(LENGTH_VAR_UNIT) << unit2 
      << std::endl;
    return 0;
  }

  inline Int Print(std::ostream &os, const char *name1, Int val1, const char *unit1,
      char *name2, Int val2, char *unit2) {
    os << std::setiosflags(std::ios::left)
      << std::setw(LENGTH_VAR_NAME) << std::string(name1)
      << std::setw(LENGTH_VAR_DATA) << val1
      << std::setw(LENGTH_VAR_UNIT) << std::string(unit1) 
      << std::setw(LENGTH_VAR_NAME) << std::string(name2)
      << std::setw(LENGTH_VAR_DATA) << val2
      << std::setw(LENGTH_VAR_UNIT) << std::string(unit2) 
      << std::endl;
    return 0;
  }

  // Bool

  // one boolean number
  inline Int Print(std::ostream &os, const std::string name, bool val) {
    os << std::setiosflags(std::ios::left)
      << std::setw(LENGTH_VAR_NAME) << name;
    if( val == true )
      os << std::setw(LENGTH_VAR_NAME) << "true" << std::endl;
    else
      os << std::setw(LENGTH_VAR_NAME) << "false" << std::endl;
    return 0;
  }


  inline Int Print(std::ostream &os, const char* name, bool val) {
    os << std::setiosflags(std::ios::left)
      << std::setw(LENGTH_VAR_NAME) << std::string(name);
    if( val == true )
      os << std::setw(LENGTH_VAR_NAME) << "true" << std::endl;
    else
      os << std::setw(LENGTH_VAR_NAME) << "false" << std::endl;
    return 0;
  }

  // *********************************************************************
  // Overload << and >> operators for basic data types
  // *********************************************************************



  // std::set
  template <class F> inline std::ostream& operator<<( std::ostream& os, const std::set<F>& s)
  {
    os<<s.size()<<std::endl;
    os.setf(std::ios_base::scientific, std::ios_base::floatfield);
    for(typename std::set<F>::iterator it = s.begin();it!=s.end();it++){
      os<<" "<<*it;
    }
    os<<std::endl;
    return os;
  }







  // std::vector
  template <class F> inline std::ostream& operator<<( std::ostream& os, const std::vector<F>& vec)
  {
    os<<vec.size()<<"| [";
    os.setf(std::ios_base::scientific, std::ios_base::floatfield);
    for(Int i=0; i<vec.size(); i++)	 
      os<<" "<<vec[i];
    os<<"]";
    return os;
  }

  template <> inline std::ostream& operator<<( std::ostream& os, const std::vector<Int>& vec)
  {
    os<<vec.size()<<"| [";
    for(Int i=0; i<vec.size(); i++)	 
      os<<" "<<vec[i];
    os<<"]";
    return os;
  }
  // Etree
  inline std::ostream& operator<<( std::ostream& os, const ETree& tree) 
  {
    os<<tree.n()<<std::endl;
    for(Int i = 1; i<=tree.Size(); i++){
      os<<tree.Parent(i-1)<<" ";
    }
    os<<std::endl;

    if(tree.IsPostOrdered()){
      for(Int i = 1; i<=tree.Size(); i++){
        os<<tree.PostParent(i-1)<<" ";
      }
      os<<std::endl;
    }

    return os;
  }

  // *********************************************************************
  // serialize/deserialize for basic types
  // More specific serialize/deserialize will be defined in individual
  // class files
  // *********************************************************************

  // standard case for most serialization/deserialization process.
  const std::vector<Int> NO_MASK(1);

  //bool
  inline Int serialize(const bool& val, std::ostream& os, const std::vector<Int>& mask)
  {
    os.write((char*)&val, sizeof(bool));
    return 0;
  }

  inline Int deserialize(bool& val, std::istream& is, const std::vector<Int>& mask)
  {
    is.read((char*)&val, sizeof(bool));
    return 0;
  }

  //char
  inline Int serialize(const char& val, std::ostream& os, const std::vector<Int>& mask)
  {
    os.write((char*)&val, sizeof(char));
    return 0;
  }

  inline Int deserialize(char& val, std::istream& is, const std::vector<Int>& mask)
  {
    is.read((char*)&val, sizeof(char));
    return 0;
  }

  inline Int combine(char& val, char& ext)
  {
    throw  std::logic_error( "Combine operation not implemented." );
  }

  //-------------------
  //Int
  inline Int serialize(const Int& val, std::ostream& os, const std::vector<Int>& mask)
  {
    os.write((char*)&val, sizeof(Int));
    return 0;
  }

  inline Int deserialize(Int& val, std::istream& is, const std::vector<Int>& mask)
  {
    is.read((char*)&val, sizeof(Int));
    return 0;
  }

  inline Int combine(Int& val, Int& ext)
  {
    val += ext;
    return 0;
  }

  //-------------------
  //Real
  inline Int serialize(const Real& val, std::ostream& os, const std::vector<Int>& mask)
  {
    os.write((char*)&val, sizeof(Real));
    return 0;
  }

  inline Int deserialize(Real& val, std::istream& is, const std::vector<Int>& mask)
  {
    is.read((char*)&val, sizeof(Real));
    return 0;
  }

  inline Int combine(Real& val, Real& ext)
  {
    val += ext;
    return 0;
  }

  //-------------------
  //Complex
  inline Int serialize(const Complex& val, std::ostream& os, const std::vector<Int>& mask)
  {
    os.write((char*)&val, sizeof(Complex));
    return 0;
  }

  inline Int deserialize(Complex& val, std::istream& is, const std::vector<Int>& mask)
  {
    is.read((char*)&val, sizeof(Complex));
    return 0;
  }

  inline Int combine(Complex& val, Complex& ext)
  {
    val += ext;
    return 0;
  }

  //-------------------
  //std::vector
  template<class T>
    Int serialize(const std::vector<T>& val, std::ostream& os, const std::vector<Int>& mask)
    {
      Int sz = val.size();
      os.write((char*)&sz, sizeof(Int));
      for(Int k=0; k<sz; k++)
        serialize(val[k], os, mask);
      return 0;
    }

  template<class T>
    Int deserialize(std::vector<T>& val, std::istream& is, const std::vector<Int>& mask)
    {
      Int sz;
      is.read((char*)&sz, sizeof(Int));
      val.resize(sz);
      for(Int k=0; k<sz; k++)
        deserialize(val[k], is, mask);
      return 0;
    }

  template<class T>
    Int combine(std::vector<T>& val, std::vector<T>& ext)
    {
      throw  std::logic_error( "Combine operation not implemented." );
      return 0;
    }

  //-------------------
  //std::set
  template<class T>
    Int serialize(const std::set<T>& val, std::ostream& os, const std::vector<Int>& mask)
    {
      Int sz = val.size();
      os.write((char*)&sz, sizeof(Int));
      for(typename std::set<T>::const_iterator mi=val.begin(); mi!=val.end(); mi++) 
        serialize((*mi), os, mask);
      return 0;
    }

  template<class T>
    Int deserialize(std::set<T>& val, std::istream& is, const std::vector<Int>& mask)
    {
      val.clear();
      Int sz;
      is.read((char*)&sz, sizeof(Int));
      for(Int k=0; k<sz; k++) {
        T t; deserialize(t, is, mask);
        val.insert(t);
      }
      return 0;
    }

  template<class T>
    Int combine(std::set<T>& val, std::set<T>& ext)
    {
      throw  std::logic_error( "Combine operation not implemented." );
      return 0;
    }

  //-------------------
  //std::map
  template<class T, class S>
    Int serialize(const std::map<T,S>& val, std::ostream& os, const std::vector<Int>& mask)
    {
      Int sz = val.size();
      os.write((char*)&sz, sizeof(Int));
      for(typename std::map<T,S>::const_iterator mi=val.begin(); mi!=val.end(); mi++) {
        serialize((*mi).first, os, mask);
        serialize((*mi).second, os, mask);
      }
      return 0;
    }

  template<class T, class S>
    Int deserialize(std::map<T,S>& val, std::istream& is, const std::vector<Int>& mask)
    {
      val.clear();
      Int sz;
      is.read((char*)&sz, sizeof(Int));
      for(Int k=0; k<sz; k++) {
        T t;	deserialize(t, is, mask);
        S s;	deserialize(s, is, mask);
        val[t] = s;
      }
      return 0;
    }

  template<class T, class S>
    Int combine(std::map<T,S>& val, std::map<T,S>& ext)
    {
      throw  std::logic_error( "Combine operation not implemented." );
      return 0;
    }

  //-------------------
  //std::pair
  template<class T, class S>
    Int serialize(const std::pair<T,S>& val, std::ostream& os, const std::vector<Int>& mask)
    {
      serialize(val.first, os, mask);
      serialize(val.second, os, mask);
      return 0;
    }

  template<class T, class S>
    Int deserialize(std::pair<T,S>& val, std::istream& is, const std::vector<Int>& mask)
    {
      deserialize(val.first, is, mask);
      deserialize(val.second, is, mask);
      return 0;
    }

  template<class T, class S>
    Int combine(std::pair<T,S>& val, std::pair<T,S>& ext)
    {
      throw  std::logic_error( "Combine operation not implemented." );
      return 0;
    }



  //-------------------
  //DistSparseMatrix
  template<class T>
    Int inline serialize(const DistSparseMatrix<T>& val, std::ostream& os, const std::vector<Int>& mask)
    {
      serialize( val.size,        os, mask );
      serialize( val.nnz,         os, mask );
      serialize( val.Local_.nnz,    os, mask );
      serialize( val.Local_.colptr, os, mask );
      serialize( val.Local_.rowind, os, mask );
      serialize( val.nzvalLocal,  os, mask );
      // No need to serialize the communicator
      return 0;
    }

  template<class T>
    Int inline deserialize(DistSparseMatrix<T>& val, std::istream& is, const std::vector<Int>& mask)
    {
      deserialize( val.size,        is, mask );
      deserialize( val.nnz,         is, mask );
      deserialize( val.Local_.nnz,    is, mask );
      deserialize( val.Local_.colptr, is, mask );
      deserialize( val.Local_.rowind, is, mask );
      deserialize( val.nzvalLocal,  is, mask );
      // No need to deserialize the communicator
      return 0;
    }

  template<class T>
    Int inline combine(DistSparseMatrix<T>& val, DistSparseMatrix<T>& ext)
    {
      throw  std::logic_error( "Combine operation not implemented." );
      return 0;
    }


  // *********************************************************************
  // Parallel IO functions
  // *********************************************************************

  Int SeparateRead(std::string name, std::istringstream& is);

  Int SeparateWrite(std::string name, std::ostringstream& os);

  Int SeparateWriteAscii(std::string name, std::ostringstream& os);

  Int SharedRead(std::string name, std::istringstream& is);

  Int SharedWrite(std::string name, std::ostringstream& os);


  // *********************************************************************
  // Comparator
  // *********************************************************************

  // Real
  inline bool PairLtComparator( const std::pair<Real, Int>& l, 
      const std::pair<Real, Int>& r ){
    return l.first < r.first;
  }

  inline bool PairGtComparator( const std::pair<Real, Int>& l, 
      const std::pair<Real, Int>& r ){
    return l.first > r.first;
  }

  // For sorting with indices
  // Example usage:
  //   std::sort(val.begin(), val.end(), IndexComp<std::vector<int>&>(indices));
  template<class T> 
    struct IndexComp {
      private: 
        const T indices_;
      public:
        IndexComp (const T indices) : indices_(indices) {}
        bool operator()(const size_t a, const size_t b) const
        { return indices_[a] < indices_[b]; }
    };




  // *********************************************************************
  // Sparse Matrix
  // *********************************************************************


  // Functions for DistSparseMatrix

  template <class F1, class F2> 
    void
    CopyPattern	( const DistSparseMatrix<F1>& A, DistSparseMatrix<F2>& B )
    {
      B.size        = A.size;
      B.nnz         = A.nnz;
      B.Local_.nnz    = A.Local_.nnz;
      B.Local_.colptr = A.Local_.colptr;
      B.Local_.rowind = A.Local_.rowind;
      B.nzvalLocal.Resize( A.Local_.nnz );
      B.comm        = A.comm;
      return ;
    }		// -----  end of template function CopyPattern  ----- 

  // *********************************************************************
  // Other numerical routines
  // *********************************************************************

  void
    LinearInterpolation ( 
        const std::vector<Real>& x, 
        const std::vector<Real>& y,
        const std::vector<Real>& xx,
        std::vector<Real>& yy );

} // namespace SYMPACK




namespace symPACK{

  template<typename T> struct is_complex_type{
    static bool value;
  };

  template< > struct is_complex_type<std::complex<float> >{
    static bool value;
  };

  template< > struct is_complex_type<std::complex<double> >{
    static bool value;
  };

  template<typename T> bool is_complex_type<T>::value = false;





  template <typename T, typename Compare>
    std::vector<std::size_t> & sort_permutation(
        const std::vector<T>& vec,
        Compare compare,
        std::vector<std::size_t>& p)
    {
      p.resize(vec.size());
      std::iota(p.begin(), p.end(), 0);
      std::sort(p.begin(), p.end(),
          [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
      return p;
    }

  template <typename T, typename Compare>
    std::vector<std::size_t> sort_permutation(
        const std::vector<T>& vec,
        Compare compare)
    {
      std::vector<std::size_t> p(vec.size());
      std::iota(p.begin(), p.end(), 0);
      std::sort(p.begin(), p.end(),
          [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
      return p;
    }

  template <typename T, typename Compare>
    std::vector<std::size_t> & sort_permutation(
        const T*begin, const T* end,
        Compare compare,
        std::vector<std::size_t>& p)
    {
      p.resize(end-begin);
      std::iota(p.begin(), p.end(), 0);
      std::sort(p.begin(), p.end(),
          [&](std::size_t i, std::size_t j){ return compare(begin[i], begin[j]); });
      return p;
    }

  template <typename T, typename Compare>
    std::vector<std::size_t> sort_permutation(
        const T*begin, const T* end,
        Compare compare)
    {
      std::vector<std::size_t> p(end-begin);
      std::iota(p.begin(), p.end(), 0);
      std::sort(p.begin(), p.end(),
          [&](std::size_t i, std::size_t j){ return compare(begin[i], begin[j]); });
      return p;
    }




  template <typename T>
    std::vector<T> apply_permutation(
        const std::vector<T>& vec,
        const std::vector<std::size_t>& p)
    {
      std::vector<T> sorted_vec(p.size());
      std::transform(p.begin(), p.end(), sorted_vec.begin(),
          [&](std::size_t i){ return vec[i]; });
      return sorted_vec;
    }

  template <typename T>
    std::vector<T> apply_permutation(
        const T*vec,
        const std::vector<std::size_t>& p)
    {
      std::vector<T> sorted_vec(p.size());
      std::transform(p.begin(), p.end(), sorted_vec.begin(),
          [&](std::size_t i){ return vec[i]; });
      return sorted_vec;
    }

  template <typename T>
    void apply_permutation(
        T*begin, T* end,
        const std::vector<std::size_t>& p)
    {
      std::vector<T> sorted_vec(p.size());
      std::transform(p.begin(), p.end(), sorted_vec.begin(),
          [&](std::size_t i){ return begin[i]; });
      std::copy(sorted_vec.begin(),sorted_vec.end(),begin);
    }

  template <typename T>
    std::vector<T> & apply_permutation(
        const std::vector<T>& vec,
        const std::vector<std::size_t>& p,
        std::vector<T> &sorted_vec)
    {
      sorted_vec.resize(p.size());
      std::transform(p.begin(), p.end(), sorted_vec.begin(),
          [&](std::size_t i){ return vec[i]; });
      return sorted_vec;
    }

  template <typename T>
    std::vector<T> & apply_permutation(
        const T*vec,
        const std::vector<std::size_t>& p,
        std::vector<T> &sorted_vec)
    {
      sorted_vec.resize(p.size());
      std::transform(p.begin(), p.end(), sorted_vec.begin(),
          [&](std::size_t i){ return vec[i]; });
      return sorted_vec;
    }

  template <typename T>
    void apply_permutation(
        T*begin, T* end,
        const std::vector<std::size_t>& p,
        std::vector<T> &sorted_vec)
    {
      sorted_vec.resize(p.size());
      std::transform(p.begin(), p.end(), sorted_vec.begin(),
          [&](std::size_t i){ return begin[i]; });
      std::copy(sorted_vec.begin(),sorted_vec.end(),begin);
    }








  template <typename SCALAR, typename INSCALAR >
    int ReadHB_PARA_MPIIO(std::string & filename, DistSparseMatrix<SCALAR> & HMat){
      abort();
            return 0;
    }




  template <typename SCALAR,typename INSCALAR>
    void ParaReadDistSparseMatrix ( const char* filename, DistSparseMatrix<SCALAR>& pspmat, MPI_Comm comm );

  template <typename SCALAR,typename INSCALAR>
    void ParaReadDistSparseMatrix ( const char* filename, DistSparseMatrix<SCALAR>& pspmat, MPI_Comm comm )
    {
      MPI_Datatype typeInt;
      MPI_Type_contiguous( sizeof(Int), MPI_BYTE, &typeInt );
      MPI_Type_commit(&typeInt);

      MPI_Datatype typePtr;
      MPI_Type_contiguous( sizeof(Ptr), MPI_BYTE, &typePtr );
      MPI_Type_commit(&typePtr);

      MPI_Datatype typeVal;
      MPI_Type_contiguous( sizeof(INSCALAR), MPI_BYTE, &typeVal );
      MPI_Type_commit(&typeVal);

      // Get the processor information within the current communicator
      MPI_Barrier( comm );
      Int mpirank;  MPI_Comm_rank(comm, &mpirank);
      Int mpisize;  MPI_Comm_size(comm, &mpisize);
      MPI_Status mpistat;
      MPI_Datatype type;
      int lens[3];
      MPI_Aint disps[3];
      MPI_Datatype types[3];

      Int err = 0;

      int filemode = MPI_MODE_RDONLY;

      MPI_File fin;
      MPI_Status status;

      // FIXME Maybe change to MPI_Comm_Dup.
      pspmat.comm = comm;

      err = MPI_File_open(comm,(char*) filename, filemode, MPI_INFO_NULL,  &fin);

      if (err != MPI_SUCCESS) {
#ifdef USE_ABORT
        abort();
#endif
        throw std::logic_error( "File cannot be opened!" );
      }

      // FIXME Note that nnz uses the Int data type for consistency of writing / reading
      // Read header
      if( mpirank == 0 ){
        err = MPI_File_read_at(fin, 0,(char*)&pspmat.size, sizeof(Int), MPI_BYTE, &status);
        err = MPI_File_read_at(fin, sizeof(Int),(char*)&pspmat.nnz, sizeof(Int), MPI_BYTE, &status);
        Ptr tmp = (Ptr)*((Int*)&pspmat.nnz);
        pspmat.nnz = tmp;
      }


      /* define a struct that describes all our data */
      lens[0] = sizeof(Int);
      lens[1] = sizeof(Ptr);
      MPI_Get_address(&pspmat.size, &disps[0]);
      MPI_Get_address(&pspmat.nnz, &disps[1]);
      types[0] = MPI_BYTE;
      types[1] = MPI_BYTE;
      MPI_Type_create_struct(2, lens, disps, types, &type);
      MPI_Type_commit(&type);


      /* broadcast the header data to everyone */
      MPI_Bcast(MPI_BOTTOM, 1, type, 0, comm);

      MPI_Type_free(&type);

      Idx row,col;
      Ptr colbeg,colend,ptr;
      int p;

      //compute local number of columns
      pspmat.Localg_.size = pspmat.size;
      pspmat.Localg_.SetComm(comm);
      pspmat.Localg_.SetBaseval(1);
      pspmat.Localg_.vertexDist.resize(mpisize+1);
      int colPerProc = std::max(1,(int)(pspmat.size/mpisize));
      for(p = 1; p <mpisize;p++){
        pspmat.Localg_.vertexDist[p] = std::min(pspmat.size-1,p*colPerProc)+pspmat.Localg_.GetBaseval();
      }
      pspmat.Localg_.vertexDist.front()= pspmat.Localg_.GetBaseval();
      pspmat.Localg_.vertexDist.back()= pspmat.size+pspmat.Localg_.GetBaseval();

      //initialize an identity permutation
      pspmat.cinvp.resize(pspmat.Localg_.LocalVertexCount());
      std::iota(pspmat.cinvp.begin(),pspmat.cinvp.end(),pspmat.Localg_.LocalFirstVertex());

      Idx firstNode = pspmat.Localg_.LocalFirstVertex() - pspmat.Localg_.GetBaseval();//0 based
      Idx numColLocal = pspmat.Localg_.LocalVertexCount();

      //initialize local arrays
      pspmat.Localg_.colptr.resize(numColLocal+1);
      pspmat.Localg_.SetKeepDiag(1);
      pspmat.Localg_.SetSorted(0);


      // Compute the number of columns on each processor
      MPI_Offset myColPtrOffset = (2 + ((mpirank==0)?0:1) )* sizeof(Int)  + (firstNode)*sizeof(Int) + (mpirank==0?sizeof(Int):0);

      Int np1 = 0;
      err= MPI_File_read_at_all(fin, myColPtrOffset, &pspmat.Localg_.colptr[0], numColLocal+1, typeInt, &status);

      if (err != MPI_SUCCESS) {
#ifdef USE_ABORT
        abort();
#endif
        throw std::logic_error( "error reading colptr" );
      }

      bassert(sizeof(Int)<=sizeof(Ptr));
      volatile Int * tmpColptr = (Int*)&pspmat.Localg_.colptr[0];
      for (col = pspmat.Localg_.colptr.size();col>0;col--){
        pspmat.Localg_.colptr[col-1] = (Ptr)tmpColptr[col-1];
      }

      // Calculate nnz_loc on each processor
      Ptr mynnz = pspmat.Localg_.colptr[numColLocal] - pspmat.Localg_.colptr[0];

      pspmat.Localg_.nnz = pspmat.nnz;
      pspmat.Localg_.rowind.resize( mynnz );

      //read rowIdx
      MPI_Offset myRowIdxOffset = (3+((mpirank==0)?0:1))*sizeof(Int) + (pspmat.size+1 + (pspmat.Localg_.colptr[0]-1))*sizeof(Int) + (mpirank==0?sizeof(Int):0);
      err= MPI_File_read_at_all(fin, myRowIdxOffset, &pspmat.Localg_.rowind[0], mynnz, typeInt,&status);

      if (err != MPI_SUCCESS) {
#ifdef USE_ABORT
        abort();
#endif
        throw std::logic_error( "error reading rowind" );
      }

      bassert(sizeof(Int)<=sizeof(Idx));
      volatile Int * tmpRowind = (Int*)&pspmat.Localg_.rowind[0];
      for (ptr = pspmat.Localg_.rowind.size(); ptr>0;ptr--){
        pspmat.Localg_.rowind[ptr-1] = (Idx)tmpRowind[ptr-1];
      }

      //read nzval
      MPI_Offset myNzValOffset = (4+ ((mpirank==0)?0:1)  )*sizeof(Int) + (pspmat.size+1 + pspmat.nnz)*sizeof(Int) + (pspmat.Localg_.colptr[0]-1)*sizeof(INSCALAR) + (mpirank==0?sizeof(Int):0);

      std::vector<INSCALAR, Mallocator<INSCALAR> > tmpBuf;
      if(std::is_same<SCALAR,INSCALAR>() ){
        pspmat.nzvalLocal.resize ( mynnz );
        err = MPI_File_read_at_all(fin, myNzValOffset, &pspmat.nzvalLocal[0], mynnz, typeVal,&status);
      }
      else{
        tmpBuf.resize(mynnz);
        err = MPI_File_read_at_all(fin, myNzValOffset, &tmpBuf[0], mynnz, typeVal,&status);
        pspmat.nzvalLocal.resize ( mynnz );
        if( is_complex_type<INSCALAR>::value ){
          for(ptr = 0; ptr<mynnz;ptr++){
            pspmat.nzvalLocal[ptr] = SCALAR( std::abs(tmpBuf[ptr]) );
          }
        }
        else{
          for(ptr = 0; ptr<mynnz;ptr++){
            pspmat.nzvalLocal[ptr] = SCALAR(std::real(tmpBuf[ptr]));
          }
        }
      }



      int error_code = err;
      if(error_code!=MPI_SUCCESS){
        char error_string[BUFSIZ];
        int length_of_error_string, error_class;

        MPI_Error_class(error_code, &error_class);
        MPI_Error_string(error_class, error_string, &length_of_error_string);
        logfileptr->OFS()<<error_string<<std::endl;
        MPI_Error_string(error_code, error_string, &length_of_error_string);
        logfileptr->OFS()<<error_string<<std::endl;

        //now check the status
        error_code = status.MPI_ERROR;
        if(error_code != MPI_SUCCESS){
          MPI_Error_class(error_code, &error_class);
          MPI_Error_string(error_class, error_string, &length_of_error_string);
          logfileptr->OFS()<<error_string<<std::endl;
          MPI_Error_string(error_code, error_string, &length_of_error_string);
          logfileptr->OFS()<<error_string<<std::endl;
        }
        gdb_lock();
      }

      if (err != MPI_SUCCESS) {
#ifdef USE_ABORT
        abort();
#endif
        throw std::logic_error( "error reading nzval" );
      }



      //convert to local references
      for( col = 1; col < numColLocal + 1; col++ ){
        pspmat.Localg_.colptr[col] = pspmat.Localg_.colptr[col] -  pspmat.Localg_.colptr[0] + 1;
      }
      pspmat.Localg_.colptr[0]=1;

      MPI_Barrier( comm );

      MPI_File_close(&fin);


      //Sort rowind + nzval accordingly
      for(col = 1; col <= numColLocal; col++){
        colbeg = pspmat.Localg_.colptr[col-1];
        colend = pspmat.Localg_.colptr[col]-1;
        std::vector<size_t> lperm = sort_permutation(&pspmat.Localg_.rowind[colbeg-1],
            &pspmat.Localg_.rowind[colend-1]+1,std::less<Idx>());
        apply_permutation(&pspmat.Localg_.rowind[colbeg-1],&pspmat.Localg_.rowind[colend-1]+1,lperm);
        apply_permutation(&pspmat.nzvalLocal[colbeg-1],&pspmat.nzvalLocal[colend-1]+1,lperm);
      }
      //force sorted boolean
      pspmat.Localg_.sorted = 1;

      pspmat.Localg_.SetSorted(1);

      MPI_Type_free(&typeVal);
      MPI_Type_free(&typePtr);
      MPI_Type_free(&typeInt);
      return ;

    }		// -----  end of function ParaReadDistSparseMatrix  ----- 




  template <typename SCALAR, typename INSCALAR >
    int ReadHB_PARA(std::string & filename, DistSparseMatrix<SCALAR> & HMat);


  template <typename SCALAR, typename INSCALAR >
    inline int ReadHB_PARA(std::string & filename, DistSparseMatrix<SCALAR> & HMat){
      MPI_Comm & workcomm = HMat.comm;

      int mpirank;
      MPI_Comm_rank(workcomm,&mpirank);
      int mpisize;
      MPI_Comm_size(workcomm,&mpisize);

      std::ifstream infile;
      infile.open(filename.c_str());

      std::string line;
      //read xadj on the first line of the input file
      std::stringstream iss;
      //skip 1st line
      if(std::getline(infile, line)){}
      Idx colptrCnt;
      Ptr rowindCnt;
      Ptr nzvalCnt;
      if(std::getline(infile, line)){
        iss.str("");
        iss.clear();
        iss<<line;
        Ptr dummy;
        iss>>dummy;
        iss>>colptrCnt>>rowindCnt>>nzvalCnt;
      }
      //read from third line

      auto & n = HMat.size;
      auto & nnz = HMat.nnz;


      auto m = n;
      if(std::getline(infile, line))
      {
        iss.str("");
        iss.clear();
        iss<<line;
        std::string type;
        iss>>type;
        iss>>m>>n>>nnz;
      }

      //compute local number of columns
      HMat.Localg_.SetComm(workcomm);
      HMat.Localg_.SetBaseval(1);
      HMat.Localg_.vertexDist.resize(mpisize+1);
      int colPerProc = std::max(1,(int)(n/mpisize));
      for(int p = 1; p <mpisize;p++){
        HMat.Localg_.vertexDist[p] = std::min(n-1,p*colPerProc)+HMat.Localg_.GetBaseval();
      }
      HMat.Localg_.vertexDist.front()= HMat.Localg_.GetBaseval();
      HMat.Localg_.vertexDist.back()= n+HMat.Localg_.GetBaseval();

      //initialize an identity permutation
      HMat.cinvp.resize(HMat.Localg_.LocalVertexCount());
      std::iota(HMat.cinvp.begin(),HMat.cinvp.end(),HMat.Localg_.LocalFirstVertex());

      Idx firstNode = HMat.Localg_.LocalFirstVertex();
      Idx nlocal = HMat.Localg_.LocalVertexCount();

      //initialize local arrays
      HMat.Localg_.colptr.resize(nlocal+1);
      HMat.Localg_.SetKeepDiag(1);
      HMat.Localg_.SetSorted(0);


      //read from 4th line
      int colptrWidth = 0;
      int rowindWidth = 0;
      int nzvalWidth = 0;
      int colptrCntPerRow = 0;
      int rowindCntPerRow = 0;
      int nzvalCntPerRow = 0;
      if(std::getline(infile, line))
      {
        iss.str("");
        iss.clear();
        iss<<line;
        std::string format;
        iss>>format;
        int dummy;
        sscanf(format.c_str(),"(%dI%d)",&colptrCntPerRow,&colptrWidth);
        iss>>format;
        sscanf(format.c_str(),"(%dI%d)",&rowindCntPerRow,&rowindWidth);
        iss>>format;
        sscanf(format.c_str(),"(%dE%d.%d)",&nzvalCntPerRow,&nzvalWidth,&dummy);
      }

      //now compute the actual number of rows
      //colptr
      size_t curPos = infile.tellg();
      size_t lineLastNode = std::ceil(( double(firstNode+nlocal) / double(colptrCntPerRow)));
      size_t lineFirstNode = std::ceil(( double(firstNode) / double(colptrCntPerRow)));
      size_t skip = (firstNode - 1)*colptrWidth + (lineFirstNode - 1);
      size_t readBytes = (nlocal+1)*colptrWidth + (lineLastNode - lineFirstNode);
      size_t skipAfter = (n+1 - (firstNode+nlocal))*colptrWidth + (colptrCnt - lineLastNode +1) ;

      infile.seekg(skip,std::ios_base::cur);

      std::string rdStr;
      rdStr.resize(readBytes);

      infile.read(&rdStr[0], readBytes);

      iss.str(rdStr);
      Ptr cptr;
      Idx locPos = 0;
      while(iss>> cptr){
        HMat.Localg_.colptr[locPos++]=cptr;
      }

      infile.seekg(skipAfter,std::ios_base::cur);
      size_t curEnd = infile.tellg();

      //convert to local colptr and compute local nnz
      Ptr first_idx = HMat.Localg_.colptr.front();
      Ptr last_idx = HMat.Localg_.colptr.back();
      for(Idx i=nlocal+1;i>=1;i--){
        HMat.Localg_.colptr[i-1] = HMat.Localg_.colptr[i-1] - HMat.Localg_.colptr[0] + 1;
      }
      Ptr nnzLocal = HMat.Localg_.colptr.back()-1;

      HMat.Localg_.rowind.resize(nnzLocal);

      Ptr elem_idx;
      //rowind
      curPos = infile.tellg();
      size_t lineLastEdge = std::ceil(( double(last_idx-1) / double(rowindCntPerRow)));
      size_t lineFirstEdge = std::ceil(( double(first_idx) / double(rowindCntPerRow)));
      skip = (first_idx - 1)*rowindWidth + (lineFirstEdge - 1);
      readBytes = (last_idx - first_idx)*rowindWidth + (lineLastEdge - lineFirstEdge);
      skipAfter = (nnz+1 - last_idx)*rowindWidth + (rowindCnt - lineLastEdge +1) ;

      infile.seekg(skip,std::ios_base::cur);

      rdStr.resize(readBytes);

      infile.read(&rdStr[0], readBytes);
      iss.str(rdStr);
      iss.clear(); // Clear state flags.
      Idx rind;
      locPos = 0;
      while(iss>> rind){
        HMat.Localg_.rowind[locPos++]=rind;
      }

      infile.seekg(skipAfter,std::ios_base::cur);
      curEnd = infile.tellg();


      HMat.nzvalLocal.resize(nnzLocal);

      //nzval
      int scalarPerNzval = is_complex_type<INSCALAR>::value?2:1;
      curPos = infile.tellg();
      lineLastEdge = std::ceil(( double(last_idx-1)*scalarPerNzval / double(nzvalCntPerRow)));
      lineFirstEdge = std::ceil(( double(first_idx*scalarPerNzval) / double(nzvalCntPerRow)));
      skip = (first_idx - 1)*scalarPerNzval*nzvalWidth + (lineFirstEdge - 1);
      readBytes = (last_idx - first_idx)*scalarPerNzval*nzvalWidth + (lineLastEdge - lineFirstEdge);

      infile.seekg(skip,std::ios_base::cur);
      rdStr.resize(readBytes);
      infile.read(&rdStr[0], readBytes);

      iss.str(rdStr);
      iss.clear(); // Clear state flags.

      SCALAR tmp(0.0);
      auto j = std::abs(tmp);
      locPos = 0;
      if(scalarPerNzval>1){
        decltype(j) * nzval_ptr = (decltype(j)*)(&HMat.nzvalLocal[0]);
        while(iss>> j){
          nzval_ptr[locPos++]=j;
        }
      }
      else{
        SCALAR * nzval_ptr = &HMat.nzvalLocal[0];
        while(iss>> j){
          nzval_ptr[locPos++]=SCALAR(j);
        }
      }

      infile.close();

      HMat.Localg_.size = HMat.size;
      HMat.Localg_.nnz = nnz;
      HMat.Localg_.expanded = 0;

      //Sort rowind + nzval accordingly
      Idx col;
      Ptr colbeg,colend;
      for(col = 1; col <= nlocal; col++){
        colbeg = HMat.Localg_.colptr[col-1];
        colend = HMat.Localg_.colptr[col]-1;
        std::vector<size_t> lperm = sort_permutation(&HMat.Localg_.rowind[colbeg-1],
            &HMat.Localg_.rowind[colend-1]+1,std::less<Idx>());
        apply_permutation(&HMat.Localg_.rowind[colbeg-1],&HMat.Localg_.rowind[colend-1]+1,lperm);
        apply_permutation(&HMat.nzvalLocal[colbeg-1],&HMat.nzvalLocal[colend-1]+1,lperm);
      }
      //force sorted boolean
      HMat.Localg_.sorted = 1;

      return 0;
    }

  template <typename SCALAR, typename INSCALAR >
    void ReadDistSparseMatrixFormatted( const char* filename, DistSparseMatrix<SCALAR>& pspmat, MPI_Comm comm );
  template <typename SCALAR, typename INSCALAR >
    void ReadDistSparseMatrixFormatted ( const char* filename, DistSparseMatrix<SCALAR>& pspmat, MPI_Comm comm )
    {
      // Get the processor information within the current communicator
      MPI_Barrier( comm );
      Int mpirank;  MPI_Comm_rank(comm, &mpirank);
      Int mpisize;  MPI_Comm_size(comm, &mpisize);
      MPI_Status mpistat;
      std::ifstream fin;

      // Read basic information
      if( mpirank == 0 ){
        fin.open(filename);
        if( !fin.good() ){
          throw std::logic_error( "File cannot be openeded!" );
        }
        Int dummy;
        fin >> pspmat.size >> dummy;
        fin >> pspmat.nnz >> dummy;
      }

      MPI_Bcast(&pspmat.size, 1, MPI_INT, 0, comm);
      MPI_Bcast(&pspmat.nnz,  1, MPI_INT, 0, comm);

      Int n = pspmat.size;
      pspmat.Localg_.size = n;
      pspmat.Localg_.nnz = pspmat.nnz;
      //compute local number of columns
      pspmat.Localg_.SetComm(comm);
      int baseval = 1;
      pspmat.Localg_.SetBaseval(baseval);
      pspmat.Localg_.vertexDist.resize(mpisize+1);
      int colPerProc = std::max(1,(int)(n/mpisize));
      for(int p = 1; p <mpisize;p++){
        pspmat.Localg_.vertexDist[p] = std::min(n-1,p*colPerProc)+baseval;
      }
      pspmat.Localg_.vertexDist.front()= baseval;
      pspmat.Localg_.vertexDist.back()= n+baseval;



      //initialize an identity permutation
      pspmat.cinvp.resize(pspmat.Localg_.LocalVertexCount());
      std::iota(pspmat.cinvp.begin(),pspmat.cinvp.end(),pspmat.Localg_.LocalFirstVertex());

      Idx firstNode = pspmat.Localg_.LocalFirstVertex();
      Idx nlocal = pspmat.Localg_.LocalVertexCount();

      pspmat.Localg_.SetKeepDiag(1);
      pspmat.Localg_.SetSorted(0);

      // Read colptr
      std::vector<Ptr>  gcolptr;
      if( mpirank == 0 ){
        gcolptr.resize(pspmat.size+1);
        Ptr* ptr = &gcolptr[0];
        for( Int i = 0; i < pspmat.size+1; i++ )
          fin >> *(ptr++);
      }

      std::vector<int> sizes(mpisize,0);
      for(int p = 0; p<mpisize; p++){
        sizes[p] = pspmat.Localg_.vertexDist[p+1] - pspmat.Localg_.vertexDist[p];
      }
      std::vector<int> disps(mpisize+1,0);
      std::partial_sum(sizes.begin(),sizes.end(),disps.begin()+1);
      //recomput the sizes to add one extra element
      for(int p = 0; p<mpisize; p++){
        sizes[p] +=1;
      }


      MPI_Datatype type;
      MPI_Type_contiguous( sizeof(Ptr), MPI_BYTE, &type );
      MPI_Type_commit(&type);
      //scatter colptr
      pspmat.Localg_.colptr.resize(nlocal+1);
      MPI_Scatterv(mpirank==0?&gcolptr[0]:nullptr,&sizes[0],&disps[0],type,&pspmat.Localg_.colptr[0],nlocal+1,type,0,comm);
      MPI_Type_free(&type);

      logfileptr->OFS()<<pspmat.Localg_.colptr<<std::endl;
      for( Int i = nlocal; i >=0; i-- ){
        pspmat.Localg_.colptr[i] -= (pspmat.Localg_.colptr[0] - baseval);
      }
      logfileptr->OFS()<<pspmat.Localg_.colptr<<std::endl;

      // Calculate nnz_loc on each processor
      Ptr nnzLocal = pspmat.Localg_.colptr.back()-1;

      pspmat.Localg_.rowind.resize( nnzLocal );
      pspmat.nzvalLocal.resize ( nnzLocal );

      // Read and distribute the row indices
      bool isLowerTri = true;
      MPI_Type_contiguous( sizeof(Idx), MPI_BYTE, &type );
      MPI_Type_commit(&type);
      if( mpirank == 0 ){
        Int tmp;
        std::vector<Idx> buf;
        Ptr numRead;
        for( Int ip = 0; ip < mpisize; ip++ ){
          numRead = gcolptr[ pspmat.Localg_.vertexDist[ip+1] - baseval ] 
            - gcolptr[ pspmat.Localg_.vertexDist[ip] - baseval ];
          Idx * ptr = nullptr;
          if(ip>0){
            buf.resize(numRead);
            ptr = &buf[0];
          }
          else{
            ptr = &pspmat.Localg_.rowind[0];
          }

          for(Idx col = pspmat.Localg_.vertexDist[ip];
              col < pspmat.Localg_.vertexDist[ip+1]; col++){
            for( Ptr rptr = gcolptr[ col  - baseval ];
                rptr < gcolptr[ col+1 - baseval ]; rptr++){
              Idx tmp = 0;
              fin >> tmp;
              if(tmp<col-baseval+1){
                isLowerTri = false;
              }
              *(ptr++) = tmp;
            }
          }
          logfileptr->OFS()<<"Sending "<<numRead<<" to P"<<ip<<std::endl;
          if( ip > 0 ){
            MPI_Send(&buf[0], numRead, type, ip, 1, comm);
          }
        }
      }
      else{
        logfileptr->OFS()<<"Expecting "<<nnzLocal<<std::endl;
        MPI_Recv( &pspmat.Localg_.rowind[0], nnzLocal, type, 0, 1, comm, &mpistat );
      }
      MPI_Type_free(&type);

      // Read and distribute the nonzero values
      MPI_Type_contiguous( sizeof(SCALAR), MPI_BYTE, &type );
      MPI_Type_commit(&type);
      if( mpirank == 0 ){
        Int tmp;
        std::vector<SCALAR> buf;
        Ptr numRead;
        for( Int ip = 0; ip < mpisize; ip++ ){
          numRead = gcolptr[ pspmat.Localg_.vertexDist[ip+1] - baseval ] 
            - gcolptr[ pspmat.Localg_.vertexDist[ip] - baseval ];
          SCALAR * ptr = nullptr;
          if(ip>0){
            buf.resize(numRead);
            ptr = &buf[0];
          }
          else{
            ptr = &pspmat.nzvalLocal[0];
          }

          if(is_complex_type<INSCALAR>::value && !is_complex_type<SCALAR>::value){
            throw std::logic_error( "Cannot convert from COMPLEX to REAL." );
          }
          else if(!is_complex_type<INSCALAR>::value && is_complex_type<SCALAR>::value){
            INSCALAR * tptr = (INSCALAR*)(ptr);
            for( Int i = 0; i < numRead; i+=2 ){
              INSCALAR val;
              fin >> val;
              *(tptr++) = val;
              fin >> val;
              *(tptr++) = val;
            }
          }
          else{
            for( Int i = 0; i < numRead; i++ ){
              SCALAR tmp;
              fin >> tmp;
              *(ptr++) = (SCALAR)tmp;
            }
          }
          if( ip > 0 ){
            MPI_Send(&buf[0], numRead, type, ip, 1, comm);
          }
        }
      }
      else{
        MPI_Recv( &pspmat.nzvalLocal[0], nnzLocal, type, 0, 1, comm, &mpistat );
      }
      MPI_Type_free(&type);

      // Close the file
      if( mpirank == 0 ){
        fin.close();
      }

      pspmat.Localg_.SetSorted(1);

      //Enforce lower triangular format
      MPI_Bcast(&isLowerTri,sizeof(bool),MPI_BYTE,0,comm);
      if(!isLowerTri){
        if(mpirank==0){
          symPACKOS<<"Input matrix is not in lower triangular format. symPACK is converting it."<<std::endl;
        }
        pspmat.Localg_.expanded = 1;
        pspmat.ToLowerTriangular();
      }

      return ;
    }		// -----  end of function ReadDistSparseMatrixFormatted  ----- 






  template <typename SCALAR, typename INSCALAR >
    void ReadMatrix(std::string & filename, std::string & informatstr,  DistSparseMatrix<SCALAR> & HMat){
      MPI_Comm & workcomm = HMat.comm;

      int mpirank;
      MPI_Comm_rank(workcomm,&mpirank);
      int mpisize;
      MPI_Comm_size(workcomm,&mpisize);
      if(mpirank==0){ symPACKOS<<"Start reading the matrix"<<std::endl; }
      SYMPACK_TIMER_START(READING_MATRIX);
      double tstart = get_time();
      //Read the input matrix
      if(informatstr == "CSC"){
        ParaReadDistSparseMatrix<SCALAR,INSCALAR>( filename.c_str(), HMat, workcomm ); 
      }
      else if (informatstr == "HB" || informatstr == "RB" || informatstr == "HARWELL_BOEING"){
        ReadHB_PARA<SCALAR,INSCALAR>(filename, HMat);
      }
      else if(informatstr == "matrix"){
        ReadDistSparseMatrixFormatted<SCALAR,INSCALAR>( filename.c_str(), HMat, workcomm );
      }
      else{
        throw std::logic_error( "Unknown matrix format." );
      }
      double tstop = get_time();
      SYMPACK_TIMER_STOP(READING_MATRIX);
      if(mpirank==0){ symPACKOS<<"Matrix read time: "<<tstop - tstart<<" seconds"<<std::endl; }
      if(mpirank==0){ symPACKOS<<"Matrix order is "<<HMat.size<<std::endl; }
    }



  void PtrSum( void *in, void *inout, int *len, MPI_Datatype *dptr ); 
  void PtrMax( void *in, void *inout, int *len, MPI_Datatype *dptr );
  template<typename T> void GenericMPIMax( void *in, void *inout, int *len, MPI_Datatype *dptr );

  template<typename T>
    void GenericMPIMax( void *in, void *inout, int *len, MPI_Datatype *dptr ) 
    { 
      int i; 

      T * pinout = (T*)inout;
      T * pin = (T*)in;
#pragma unroll
      for (i=0; i< *len; ++i) { 
        pinout[i] = std::max(pinout[i], pin[i]);
      } 
    }

  template<typename T>
    inline std::string ToMatlabScalar( std::complex<T> val){
      std::stringstream s;
      s.precision(std::numeric_limits< std::complex<T> >::max_digits10);
      s.precision(15);
      s<<"complex("<<std::scientific<<std::real(val)<<","<<std::imag(val)<<")";
      return s.str();
    }

  template<typename T>
    inline std::string ToMatlabScalar( T val){
      std::stringstream s;
      s.precision(std::numeric_limits< T >::max_digits10);
      s<<std::scientific<<val;
      return s.str();
    }




} // namespace SYMPACK

#endif // _UTILITY_HPP_
