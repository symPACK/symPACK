/// @file Cholesky.hpp
/// @brief TODO
/// @author Mathias Jacquelin
/// @date 2013-09-09
#ifndef _CHOLESKY_FACT_HPP_
#define _CHOLESKY_FACT_HPP_

#include "Environment.hpp"
#include "NumVec.hpp"
#include "SuperNode.hpp"

namespace LIBCHOLESKY{

template <class F> void LL_SS_Fact(DistSparseMatrix<F> & Amat, IntNumVec & xsuper, IntNumVec & cc, IntNumVec & rc );




} // namespace LIBCHOLESKY


#endif // _CHOLESKY_FACT_HPP_
