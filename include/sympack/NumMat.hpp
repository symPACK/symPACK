/// @file nummat_decl.hpp
/// @brief Numerical matrix.
/// @author Lexing Ying and Lin Lin
/// @date 2010-09-27
#ifndef _NUMMAT_DECL_HPP_
#define _NUMMAT_DECL_HPP_

#include "sympack/Environment.hpp"
//#include "sympack/SuperNode.hpp"



namespace SYMPACK{

//template <typename T> class SuperNode;

/// @class NumMat
///
/// @brief Numerical matrix.
///
/// NumMat is a portable encapsulation of a pointer to represent a 2D
/// matrix, which can either own (owndata == true) or view (owndata ==
/// false) a piece of data.  
template <typename F>
class NumMat
{
  protected:
    virtual void error_message(std::stringstream & ss, Int i, Int j);
    virtual void error_message(std::stringstream & ss, Int j);
    virtual inline void alloc_data();
    virtual inline void delete_data();
  public:
    /// @brief The size of the first dimension.
    Int m_ ; 

    /// @brief The size of second dimension.
    Int n_ ;

    /// @brief Whether it owns the data.
    bool owndata_;

#ifdef _ASSERT_
    bool allocated_ ;
#endif

    /// @brief The pointer for the actual data.
    F* data_;

    NumMat(Int m=0, Int n=0);
    NumMat(Int m, Int n, bool owndata, F* data);
    NumMat(const NumMat& C);
//    NumMat(const SuperNode<F>& S);
    virtual ~NumMat();
    NumMat& Copy(const NumMat& C);
//    NumMat& SnodeToDense(const SuperNode<F>& S);
    virtual void Resize(Int m, Int n);
    virtual void Clear();

    NumMat& operator=(const NumMat& C);
//    NumMat& operator=(const SuperNode<F>& S);

    inline const F& operator()(Int i, Int j) const;
    inline F& operator()(Int i, Int j);

    inline const F& at(Int i, Int j) const;
#ifdef _DEBUG_
    F& at(Int i, Int j);
#else
    inline F& at(Int i, Int j);
#endif



    F* Data() const;
    F* VecData(Int j) const;
    Int m() const { return m_; }
    Int n() const { return n_; }
    Int Size() const { return m_*n_; }
    Int ByteSize() const { return m_*n_*sizeof(F); }
};

// Commonly used
typedef NumMat<bool>     BolNumMat;
typedef NumMat<Int>      IntNumMat;
typedef NumMat<Real>     DblNumMat;
typedef NumMat<Complex>  CpxNumMat;

// *********************************************************************
// Utility functions
// *********************************************************************
/// @brief SetValue sets a numerical matrix to a constant val.
template <typename F> inline void SetValue(NumMat<F>& M, F val);

/// @brief Energy computes the L2 norm of a matrix (treated as a SYMPACK::vector).
template <typename F> Real Energy(const NumMat<F>& M);
 
template <typename F> inline void Transpose ( const NumMat<F>& A, NumMat<F>& B );
template <typename F> inline void Symmetrize( NumMat<F>& A );
} // namespace SYMPACK

#include "sympack/NumMat_impl.hpp"

#endif // _NUMMAT_DECL_HPP_
