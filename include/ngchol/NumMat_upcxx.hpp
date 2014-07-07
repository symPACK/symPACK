/// @file NumMat_upcxx.hpp
/// @brief Numerical matrix.
/// @author Mathias Jacquelin
/// @date 2014-02-05
#ifndef _NUMMAT_UPCXX_HPP_
#define _NUMMAT_UPCXX_HPP_

#include "ngchol/Environment.hpp"
#include "ngchol/NumMat.hpp"

#include <upcxx.h>

namespace LIBCHOLESKY{

/// @class NumMat_upcxx
///
/// @brief Numerical matrix.
///
/// NumMat is a portable encapsulation of a pointer to represent a 2D
/// matrix, which can either own (owndata == true) or view (owndata ==
/// false) a piece of data.  
  template <typename F> class NumMat_upcxx : public NumMat<F>{
    protected:
      virtual void error_message(std::stringstream & ss, Int i, Int j);
      virtual void error_message(std::stringstream & ss, Int j);
      virtual inline void alloc_data();
      virtual inline void delete_data();
    public:

      upcxx::global_ptr<F> gdata_;

      NumMat_upcxx(Int m=0, Int n=0);
      NumMat_upcxx(Int m, Int n, bool owndata, F* data);
      NumMat_upcxx(const NumMat_upcxx<F>& C);
      virtual ~NumMat_upcxx();
      virtual NumMat_upcxx& Copy(const NumMat_upcxx<F>& C);
      virtual void Resize(Int m, Int n);
      virtual void Clear();

      upcxx::global_ptr<F> GData() const;
      upcxx::global_ptr<F> GVecData(Int j) const;
  };

// Commonly used
typedef NumMat_upcxx<bool>     BolNumMat_upcxx;
typedef NumMat_upcxx<Int>      IntNumMat_upcxx;
typedef NumMat_upcxx<Real>     DblNumMat_upcxx;
typedef NumMat_upcxx<Complex>  CpxNumMat_upcxx;

} // namespace LIBCHOLESKY

#include "ngchol/NumMat_upcxx_impl.hpp"

#endif // _NUMMAT_UPCXX_HPP_
