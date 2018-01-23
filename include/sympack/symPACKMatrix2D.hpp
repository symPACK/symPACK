/*
   "symPACK" Copyright (c) 2016, The Regents of the University of California,
   through Lawrence Berkeley National Laboratory (subject to receipt of any
   required approvals from the U.S. Dept. of Energy).  All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   (1) Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

   (2) Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

   (3) Neither the name of the University of California, Lawrence Berkeley
   National Laboratory, U.S. Dept. of Energy nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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
#ifndef _SYMPACK_MATRIX2D_DECL_HPP_
#define _SYMPACK_MATRIX2D_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/symPACKMatrix.hpp"




namespace symPACK{


  template <typename colptr_t, typename rowind_t, typename T> 
    class symPACKMatrix2D: public symPACKMatrixMeta<T>{
    public:

    class block_t {
      public:
        std::vector<T> nzval;
    };

    using blockCell_t = std::vector<block_t>;

    protected:
      std::vector<blockCell_t> cells_;

    public:

      symPACKMatrix2D();
      ~symPACKMatrix2D();



      void Init(symPACKOptions & options );
      void Init(DistSparseMatrix<T> & pMat, symPACKOptions & options );

      void SymbolicFactorization(DistSparseMatrix<T> & pMat);
      void DistributeMatrix(DistSparseMatrix<T> & pMat);




  };

  template <typename colptr_t, typename rowind_t, typename T>
  symPACKMatrix2D<colptr_t,rowind_t,T>::symPACKMatrix2D():
    symPACKMatrixMeta<T>()
  {
  }

  template <typename colptr_t, typename rowind_t, typename T>
  symPACKMatrix2D<colptr_t,rowind_t,T>::~symPACKMatrix2D(){
  }

  template <typename colptr_t, typename rowind_t, typename T>
   void symPACKMatrix2D<colptr_t,rowind_t,T>::Init(symPACKOptions & options ){
   } 

  template <typename colptr_t, typename rowind_t, typename T>
   void symPACKMatrix2D<colptr_t,rowind_t,T>::Init(DistSparseMatrix<T> & pMat,symPACKOptions & options ){
   } 

  template <typename colptr_t, typename rowind_t, typename T>
   void symPACKMatrix2D<colptr_t,rowind_t,T>::SumbolicFactorization(DistSparseMatrix<T> & pMat ){
   } 

  template <typename colptr_t, typename rowind_t, typename T>
   void symPACKMatrix2D<colptr_t,rowind_t,T>::DistributeMatrix(DistSparseMatrix<T> & pMat ){
   } 
}

//#include <sympack/impl/symPACKMatrix2D_impl.hpp>
#endif //_SYMPACK_MATRIX2D_DECL_HPP_

