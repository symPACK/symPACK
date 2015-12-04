/// @file utility.hpp
/// @brief Various utility subroutines.
/// @author Lin Lin & Mathias Jacquelin
/// @date 2012-09-27
#ifndef _UTILITY_HPP_ 
#define _UTILITY_HPP_

#include  "ngchol/ETree.hpp"

#include <stdlib.h>
#include <iostream>
#include <set>
#include <vector>

namespace LIBCHOLESKY{


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
	for(int i=0; i<vec.size(); i++)	 
		os<<" "<<vec[i];
	os<<std::endl;
	return os;
}


// Etree
inline std::ostream& operator<<( std::ostream& os, const ETree& tree) 
{
  os<<tree.n()<<std::endl;
  for(int i = 1; i<=tree.Size(); i++){
//    os<<i<<" ["<<tree.parent_(i-1)<<"] ";
    os<<tree.Parent(i-1)<<" ";
  }
  os<<std::endl;

  if(tree.IsPostOrdered()){
    //for(int i = 1; i<=tree.parent_.m(); i++){
    for(int i = 1; i<=tree.Size(); i++){
//      os<<i<<" ["<<tree.PostParent(i-1)<<"] ";
      os<<tree.PostParent(i-1)<<" ";
    }
    os<<std::endl;
  }

  return os;
}




} // namespace LIBCHOLESKY
#endif // _UTILITY_HPP_
