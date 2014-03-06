/// @file SuperNode.hpp
/// @brief TODO
/// @author Mathias Jacquelin
/// @date 2013-09-10
#ifndef _SUPERNODE_FACT_HPP_
#define _SUPERNODE_FACT_HPP_

#include "Environment.hpp"
#include "NumVec.hpp"
#include "NZBlock.hpp"

namespace LIBCHOLESKY{

class SuperNode{
  public:
  Int id;
  Int size;
  Int firstCol;
  Int lastCol;
  SuperNode(){};      
  std::vector<NZBlock<double> > Lcol;
  std::vector<Int > LocalIndices;
};



} // namespace LIBCHOLESKY


#endif // _SUPERNODE_FACT_HPP_
