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
/// @file sparse_matrix_decl.hpp
/// @brief Sparse matrix and Distributed sparse matrix in compressed
/// column format.
/// @author Lin Lin
/// @date 2012-11-10
#ifndef _DIST_SPARSE_MATRIX_DECL_HPP_
#define _DIST_SPARSE_MATRIX_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/ETree.hpp"

#include <mpi.h>


namespace symPACK{


class DistSparseMatrixGraph;
template <typename T> class symPACKMatrix;

template <typename F> class DistSparseMatrix{
  friend class symPACKMatrix<F>;
  friend class symPACKMatrix<F>;
  //friend functions
  friend void ReadDistSparseMatrixFormatted ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm );
  friend void ReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm );
  friend void ParaWriteDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm );
  friend void ParaReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm );

  friend void ReadDistSparseMatrixFormatted ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm );
  friend void ReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm );
  friend void ParaWriteDistSparseMatrix ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm );
  friend void ParaReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm );

  template <typename SCALAR, typename INSCALAR >
  friend int ReadHB_PARA(std::string & filename, DistSparseMatrix<SCALAR> & HMat);

  template <typename SCALAR, typename INSCALAR >
  friend int ReadHB_PARA_MPIIO(std::string & filename, DistSparseMatrix<SCALAR> & HMat);


  public:
  DistSparseMatrixGraph Localg_;
	/// @brief Matrix dimension.
	Int          size;

	/// @brief Total number of nonzeros elements.
	//Idx64          nnz;                             
	Ptr          nnz;                             

	/// @brief Dimension nnzLocal, storing the nonzero values.
	std::vector<F>    nzvalLocal;                      

	/// @brief MPI communicator
	MPI_Comm     comm;        

  /// @brief Permutation currently used in the matrix. Empty vector means unpermuted matrix
  std::vector<Int> cinvp;

  inline bool isPermuted(Int * invp, Idx * vertexDist = NULL);

  DistSparseMatrix();
  DistSparseMatrix(MPI_Comm aComm);
  void CopyData(const int n, const int nnz, const int * colptr, const int * rowidx, const F * nzval,bool onebased=false);
  DistSparseMatrix(const int n, const int nnz, const int * colptr, const int * rowidx, const F * nzval , MPI_Comm oComm);
  
  template <typename T> void ConvertData(const int n, const int nnz, const int * colptr, const int * rowidx, const T * nzval,bool onebased=false);

  const DistSparseMatrixGraph & GetLocalGraph() const;
  DistSparseMatrixGraph & GetLocalGraph();

  void SetLocalGraph(const DistSparseMatrixGraph & pgraph);

  void Dump() const;
  void DumpMatlab() const;

  void Permute(Int * invp, Idx * newVertexDist = NULL);
  void ExpandSymmetric();
  void ToLowerTriangular();
  void SortGraph();
};

// Commonly used
typedef DistSparseMatrix<Real>       DblDistSparseMatrix;
typedef DistSparseMatrix<Complex>    CpxDistSparseMatrix;



} // namespace SYMPACK

#include "sympack/DistSparseMatrix_impl.hpp"

#endif // _DIST_SPARSE_MATRIX_DECL_HPP_
