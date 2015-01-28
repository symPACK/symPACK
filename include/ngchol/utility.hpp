/// @file utility.hpp
/// @brief Various utility subroutines.
/// @author Lin Lin & Mathias Jacquelin
/// @date 2012-09-27
#ifndef _UTILITY_HPP_ 
#define _UTILITY_HPP_

#include  <mpi.h>
#include  <stdlib.h>
#include  "ngchol/Environment.hpp"
#include  "ngchol/NumVec.hpp"
#include  "ngchol/NumMat.hpp"
#include  "ngchol/SparseMatrix.hpp"
#include  "ngchol/DistSparseMatrix.hpp"
#include  "ngchol/ETree.hpp"

namespace LIBCHOLESKY{


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

// String
inline Int Print(std::ostream &os, const std::string name) {
  os << std::setiosflags(std::ios::left) << name << std::endl;
  return 0;
};

inline Int Print(std::ostream &os, const char* name) {
  os << std::setiosflags(std::ios::left) << std::string(name) << std::endl;
  return 0;
};

inline Int Print(std::ostream &os, const std::string name, std::string val) {
  os << std::setiosflags(std::ios::left) 
     << std::setw(LENGTH_VAR_NAME) << name
     << std::setw(LENGTH_VAR_DATA) << val
     << std::endl;
  return 0;
};

inline Int Print(std::ostream &os, const std::string name, const char* val) {
  os << std::setiosflags(std::ios::left) 
     << std::setw(LENGTH_VAR_NAME) << name
     << std::setw(LENGTH_VAR_DATA) << std::string(val)
     << std::endl;
  return 0;
};


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
};

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
};


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
};

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
};

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
};

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
};

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
};

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
};

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
};


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
};

// Int

// one integer number
inline Int Print(std::ostream &os, std::string name, Int val) {
  os << std::setiosflags(std::ios::left)
     << std::setw(LENGTH_VAR_NAME) << name
     << std::setw(LENGTH_VAR_DATA) << val
     << std::endl;
  return 0;
};


inline Int Print(std::ostream &os, const char *name, Int val) {
  os << std::setiosflags(std::ios::left)
     << std::setw(LENGTH_VAR_NAME) << std::string(name)
     << std::setw(LENGTH_VAR_DATA) << val
     << std::endl;
  return 0;
};


inline Int Print(std::ostream &os, const std::string name, Int val, const std::string unit) {
  os << std::setiosflags(std::ios::left)
     << std::setw(LENGTH_VAR_NAME) << name
     << std::setw(LENGTH_VAR_DATA) << val
     << std::setw(LENGTH_VAR_UNIT) << unit 
     << std::endl;
  return 0;
};


inline Int Print(std::ostream &os, const char* name, Int val, const std::string unit) {
  os << std::setiosflags(std::ios::left)
     << std::setw(LENGTH_VAR_NAME) << std::string(name)
     << std::setw(LENGTH_VAR_DATA) << val
     << std::setw(LENGTH_VAR_UNIT) << unit 
     << std::endl;
  return 0;
};



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
};

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
};

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
};


inline Int Print(std::ostream &os, const char* name, bool val) {
  os << std::setiosflags(std::ios::left)
     << std::setw(LENGTH_VAR_NAME) << std::string(name);
	if( val == true )
     os << std::setw(LENGTH_VAR_NAME) << "true" << std::endl;
	else
     os << std::setw(LENGTH_VAR_NAME) << "false" << std::endl;
  return 0;
};

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
	os<<vec.size()<<std::endl;
	os.setf(std::ios_base::scientific, std::ios_base::floatfield);
	for(Int i=0; i<vec.size(); i++)	 
		os<<" "<<vec[i];
	os<<std::endl;
	return os;
}

// NumVec
template <typename F, typename TIdx> inline std::ostream& operator<<( std::ostream& os, const NumVec<F,TIdx>& vec)
{
	os<<vec.m()<<std::endl;
	os.setf(std::ios_base::scientific, std::ios_base::floatfield);
	for(TIdx i=0; i<vec.m(); i++){ 
		os<<" "<<vec(i);
    if(((i+1)%10)==0){
	    os<<std::endl;
    }
  }
	os<<std::endl;
	return os;
}

//template <class F> inline std::istream& operator>>( std::istream& is, NumVec<F>& vec)
//{
//	Int m;  is>>m;  vec.resize(m);
//	for(Int i=0; i<vec.m(); i++)	 
//		is >> vec(i);
//	return is;
//}

// NumMat
template <class F> inline std::ostream& operator<<( std::ostream& os, const NumMat<F>& mat)
{
  os<<mat.m()<<" "<<mat.n()<<std::endl;
  os.setf(std::ios_base::scientific, std::ios_base::floatfield);
  os << std::setprecision(16);
  for(Int i=0; i<mat.m(); i++) {
    for(Int j=0; j<mat.n(); j++)
      os<<" "<<mat(i,j);
    os<<std::endl;
  }
  return os;
}


// Etree
inline std::ostream& operator<<( std::ostream& os, const ETree& tree) 
{
  os<<tree.n()<<std::endl;
  for(Int i = 1; i<=tree.Size(); i++){
//    os<<i<<" ["<<tree.parent_(i-1)<<"] ";
    os<<tree.Parent(i-1)<<" ";
  }
  os<<std::endl;

  if(tree.IsPostOrdered()){
    //for(Int i = 1; i<=tree.parent_.m(); i++){
    for(Int i = 1; i<=tree.Size(); i++){
//      os<<i<<" ["<<tree.PostParent(i-1)<<"] ";
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

/*
//-------------------
//BolNumVec
inline Int serialize(const BolNumVec& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)(val.Data()), m*sizeof(bool));
  return 0;
}

inline Int deserialize(BolNumVec& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  is.read((char*)&m, sizeof(Int));
  val.Resize(m);
  is.read((char*)(val.Data()), m*sizeof(bool));
  return 0;
}

//-------------------
//BolNumMat
inline Int serialize(const BolNumMat& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  Int n = val.n();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)&n, sizeof(Int));
  os.write((char*)(val.Data()), m*n*sizeof(bool));
  return 0;
}

inline Int deserialize(BolNumMat& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  Int n;
  is.read((char*)&m, sizeof(Int));
  is.read((char*)&n, sizeof(Int));
  val.Resize(m,n);
  is.read((char*)(val.Data()), m*n*sizeof(bool));
  return 0;
}

*/

//-------------------
//IntNumVec
inline Int serialize(const IntNumVec& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)(val.Data()), m*sizeof(Int));
  return 0;
}

inline Int deserialize(IntNumVec& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  is.read((char*)&m, sizeof(Int));
  val.Resize(m);
  is.read((char*)(val.Data()), m*sizeof(Int));
  return 0;
}

inline Int combine(IntNumVec& val, IntNumVec& ext)
{
  //val.resize(ext.m());
  assert(val.m()==ext.m());
  for(Int i=0; i<val.m(); i++)    val(i) += ext(i);
  return 0;
}


//-------------------
//IntNumMat
inline Int serialize(const IntNumMat& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  Int n = val.n();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)&n, sizeof(Int));
  os.write((char*)(val.Data()), m*n*sizeof(Int));
  return 0;
}

inline Int deserialize(IntNumMat& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  Int n;
  is.read((char*)&m, sizeof(Int));
  is.read((char*)&n, sizeof(Int));
  val.Resize(m,n);
  is.read((char*)(val.Data()), m*n*sizeof(Int));
  return 0;
}

inline Int combine(IntNumMat& val, IntNumMat& ext)
{
  //val.resize(ext.m(),ext.n());
  assert(val.m()==ext.m() && val.n()==ext.n());
  for(Int i=0; i<val.m(); i++)
    for(Int j=0; j<val.n(); j++)
      val(i,j) += ext(i,j);
  return 0;
}

//-------------------
//DblNumVec
inline Int serialize(const DblNumVec& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)(val.Data()), m*sizeof(Real));
  return 0;
}

inline Int deserialize(DblNumVec& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  is.read((char*)&m, sizeof(Int));
  val.Resize(m);
  is.read((char*)(val.Data()), m*sizeof(Real));
  return 0;
}

inline Int combine(DblNumVec& val, DblNumVec& ext)
{
  //val.resize(ext.m());
  assert(val.m()==ext.m());
  for(Int i=0; i<val.m(); i++)    val(i) += ext(i);
  return 0;
}

//-------------------
//DblNumMat
inline Int serialize(const DblNumMat& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  Int n = val.n();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)&n, sizeof(Int));
  os.write((char*)(val.Data()), m*n*sizeof(Real));
  return 0;
}

inline Int deserialize(DblNumMat& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  Int n;
  is.read((char*)&m, sizeof(Int));
  is.read((char*)&n, sizeof(Int));
  val.Resize(m,n);
  is.read((char*)(val.Data()), m*n*sizeof(Real));
  return 0;
}

inline Int combine(DblNumMat& val, DblNumMat& ext)
{
  //val.resize(ext.m(),ext.n());
  assert(val.m()==ext.m() && val.n()==ext.n());
  for(Int i=0; i<val.m(); i++)
    for(Int j=0; j<val.n(); j++)
      val(i,j) += ext(i,j);
  return 0;
}




//-------------------
//CpxNumVec
inline Int serialize(const CpxNumVec& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)(val.Data()), m*sizeof(Complex));
  return 0;
}

inline Int deserialize(CpxNumVec& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  is.read((char*)&m, sizeof(Int));
  val.Resize(m);
  is.read((char*)(val.Data()), m*sizeof(Complex));
  return 0;
}

inline Int combine(CpxNumVec& val, CpxNumVec& ext)
{
  //val.resize(ext.m());
  assert(val.m()==ext.m());
  for(Int i=0; i<val.m(); i++)    val(i) += ext(i);
  return 0;
}

//-------------------
//CpxNumMat
inline Int serialize(const CpxNumMat& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  Int n = val.n();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)&n, sizeof(Int));
  os.write((char*)(val.Data()), m*n*sizeof(Complex));
  return 0;
}

inline Int deserialize(CpxNumMat& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  Int n;
  is.read((char*)&m, sizeof(Int));
  is.read((char*)&n, sizeof(Int));
  val.Resize(m,n);
  is.read((char*)(val.Data()), m*n*sizeof(Complex));
  return 0;
}

inline Int combine(CpxNumMat& val, CpxNumMat& ext)
{
  //val.resize(ext.m(),ext.n());
  assert(val.m()==ext.m() && val.n()==ext.n());
  for(Int i=0; i<val.m(); i++)
    for(Int j=0; j<val.n(); j++)
      val(i,j) += ext(i,j);
  return 0;
}



//-------------------
//NumVec
template<typename T, typename TIdx>
Int inline serialize(const NumVec<T,TIdx>& val, std::ostream& os, const std::vector<Int>& mask)
{
  TIdx m = val.m();
  os.write((char*)&m, sizeof(TIdx));
  for(TIdx i=0; i<m; i++)
	serialize(val(i), os, mask);
  return 0;
}

template<typename T, typename TIdx>
Int inline deserialize(NumVec<T,TIdx>& val, std::istream& is, const std::vector<Int>& mask)
{
  TIdx m;
  is.read((char*)&m, sizeof(TIdx));
  val.Resize(m);
  for(TIdx i=0; i<m; i++)
    deserialize(val(i), is, mask);
  return 0;
}

template<typename T, typename TIdx>
Int inline combine(NumVec<T,TIdx>& val, NumVec<T,TIdx>& ext)
{
	throw  std::logic_error( "Combine operation not implemented." );
  return 0;
}

//-------------------
//NumMat
template<class T>
Int inline serialize(const NumMat<T>& val, std::ostream& os, const std::vector<Int>& mask)
{
  Int m = val.m();
  Int n = val.n();
  os.write((char*)&m, sizeof(Int));
  os.write((char*)&n, sizeof(Int));
  for(Int j=0; j<n; j++)
	for(Int i=0; i<m; i++)
	  serialize(val(i,j), os, mask);
  return 0;
}
template<class T>
Int inline deserialize(NumMat<T>& val, std::istream& is, const std::vector<Int>& mask)
{
  Int m;
  Int n;
  is.read((char*)&m, sizeof(Int));
  is.read((char*)&n, sizeof(Int));
  val.Resize(m,n);
  for(Int j=0; j<n; j++)
	for(Int i=0; i<m; i++)
	  deserialize(val(i,j), is, mask);
  return 0;
}

template<class T>
Int inline combine(NumMat<T>& val, NumMat<T>& ext)
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
// Ej
// *********************************************************************
inline void IdentityCol( Int col, NumVec<Real,Int>& vec )
{
	for(Int i=0; i<std::min(col,vec.m()); i++)
		vec(i) = 0.0;
  
  if(col<vec.m())
    vec(col) = 1.0;


	for(Int i=col+1; i<vec.m(); i++)
		vec(i) = 0.0;
}

inline void IdentityCol( Int col, NumVec<Complex,Int>& vec )
{
	for(Int i=0; i<std::min(col,vec.m()); i++)
		vec(i) = Complex(0.0,0.0);
  
  if(col<vec.m())
    vec(col) = Complex(1.0,0.0);


	for(Int i=col+1; i<vec.m(); i++)
		vec(i) = Complex(0.0,0.0);
}



inline void IdentityCol( IntNumVec & cols, NumMat<Real>& mat )
{
  for(Int j=0;j<cols.m();j++){
    Int col= cols(j);
    for(Int i=0; i<std::min(col,mat.m()); i++)
      mat(i,j) = 0.0;

    if(col<mat.m())
      mat(col,j) = 1.0;


    for(Int i=col+1; i<mat.m(); i++)
      mat(i,j) = 0.0;
  }
}


inline void IdentityCol( IntNumVec & cols, NumMat<Complex>& mat )
{
  for(Int j=0;j<cols.m();j++){
    Int col= cols(j);
	for(Int i=0; i<std::min(col,mat.m()); i++)
		mat(i,j) = Complex(0.0,0.0);
  
  if(col<mat.m())
    mat(col,j) = Complex(1.0,0.0);


	for(Int i=col+1; i<mat.m(); i++)
		mat(i,j) = Complex(0.0,0.0);
}
}




// *********************************************************************
// Random numbers
// *********************************************************************
inline void SetRandomSeed(long int seed){
	srand48(seed);
}

inline Real UniformRandom(){
	return (Real)drand48();
}


inline void UniformRandom( Real * vec, Int size )
{
	for(Int i=0; i<size; i++)
		vec[i] = UniformRandom();
}




inline void UniformRandom( NumVec<Real,Int>& vec )
{
	for(Int i=0; i<vec.m(); i++)
		vec(i) = UniformRandom();
}

inline void UniformRandom( NumVec<Complex,Int>& vec )
{
	for(Int i=0; i<vec.m(); i++)
		vec(i) = Complex(UniformRandom(), UniformRandom());
}

inline void UniformRandom( NumMat<Real>& M )
{
	Real *ptr = M.Data();
  for(Int i=0; i < M.m() * M.n(); i++) 
		*(ptr++) = UniformRandom(); 
}

inline void UniformRandom( NumMat<Complex>& M )
{
	Complex *ptr = M.Data();
  for(Int i=0; i < M.m() * M.n(); i++) 
		*(ptr++) = Complex(UniformRandom(), UniformRandom()); 
}




// *********************************************************************
// Timing
// *********************************************************************
inline void GetTime(Real&  t){
	t = MPI_Wtime();
}

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

// TODO Complex format
void ReadSparseMatrix ( const char* filename, SparseMatrix<Real>& spmat );

void ReadDistSparseMatrix( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm );

void ParaReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm );

void ParaWriteDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm );

void ReadDistSparseMatrixFormatted( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm );

template <class F1, class F2> 
void
CopyPattern	( const SparseMatrix<F1>& A, SparseMatrix<F2>& B )
{
	B.size        = A.size;
	B.nnz         = A.nnz;
	B.colptr      = A.colptr;
	B.rowind      = A.rowind;
	B.nzval.Resize( A.nnz );
	return ;
}		// -----  end of template function CopyPattern  ----- 


// TODO Real format
void
GetDiagonal ( const DistSparseMatrix<Complex>& A, 
		NumVec<Complex,Int>& diag );


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





















} // namespace LIBCHOLESKY
#endif // _UTILITY_HPP_
