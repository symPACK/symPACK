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
#ifndef _DIST_SPARSE_MATRIX_GRAPH_HPP_
#define _DIST_SPARSE_MATRIX_GRAPH_HPP_

#include "sympack/Environment.hpp"
#include "sympack/Ordering.hpp"
#include "sympack/Types.hpp"

namespace symPACK{
class ETree;
class Ordering;
class DistSparseMatrixGraph;
class SparseMatrixGraph;
template <typename F> class DistSparseMatrix;

class SparseMatrixGraph{
  friend class Ordering;
  friend class DistSparseMatrixGraph;
  public:
  Idx          size;                            // Matrix dimension (global)
  Ptr          nnz;                             // Number of nonzeros (global)
  std::vector<Ptr>  colptr;                 // Column index pointer
  std::vector<Idx>  rowind;                 // Starting row index pointer
  public:
  int baseval;
  int keepDiag;
  int sorted;
  int expanded;

  public:
  void SetBaseval(int aBaseval);
  void SetKeepDiag(int aKeepDiag);
  void SetSorted(int aSorted);

  int IsExpanded() const {return expanded;}
  int GetBaseval() const {return baseval;}
  int GetKeepDiag()const {return keepDiag;}
  int GetSorted() const {return sorted;}

  SparseMatrixGraph();
  void SortEdges();
  void DistributeGraph(DistSparseMatrixGraph & dg);
  void ExpandSymmetric();
  void BroadcastGraph(MPI_Comm comm, Int root = 0);

  Idx VertexCount() const { return colptr.size()==0?0:colptr.size()-1;}
  Ptr EdgeCount() const{ return rowind.size();}
  Idx LocalFirstVertex() const {return baseval;}
  Idx LocalVertexCount() const { return VertexCount();}
  Ptr LocalEdgeCount() const{ return EdgeCount();}


};

class DistSparseMatrixGraph{
  friend class Ordering;
  template <typename F> friend class DistSparseMatrix;


  protected:
  public:

  public:
	Idx          size;                            // Matrix dimension (global)
  Ptr          nnz;                             // Number of nonzeros (global)
	std::vector<Ptr>  colptr;                 // Column index pointer
	std::vector<Idx>  rowind;                 // Starting row index pointer
	std::vector<Idx>  vertexDist;             // vertex distribution array (size np+1)
  MPI_Comm comm;

  public:
  int mpirank;
  int mpisize;
  int baseval;
  int keepDiag;
  int sorted;
  int expanded;


  public:
    void SetComm(const MPI_Comm & aComm);
    void SetBaseval(int aBaseval);
    void SetKeepDiag(int aKeepDiag);
    void SetSorted(int aSorted);
    void SetExpanded(int aExpanded);


    MPI_Comm const & GetComm() const {return comm;}
    int GetBaseval() const {return baseval;}
    int GetKeepDiag()const {return keepDiag;}
    int GetSorted() const {return sorted;}

  public:
  DistSparseMatrixGraph();
  ~DistSparseMatrixGraph();
  DistSparseMatrixGraph( const DistSparseMatrixGraph& g );
  //constructor from global to local

  DistSparseMatrixGraph& operator=( const DistSparseMatrixGraph& g );
  //accessors
  bool IsExpanded() const {return expanded;}
  Idx LocalFirstVertex() const {return vertexDist[mpirank];}
  Idx LocalVertexCount() const { Idx count = vertexDist[mpirank+1] - vertexDist[mpirank]; bassert(colptr.size()==count+1 || colptr.empty() ); return count; }
  Ptr LocalEdgeCount() const{ return rowind.size();}

  //utility

  void SortEdges();
  void ExpandSymmetric();
  void ToSymmetric();

  void Permute(Int * invp);
  void Permute(Int * invp, Idx * newVertexDist);
  void Permute(Int * invp, Int invpbaseval);
  void Permute(Int * invp, Idx * newVertexDist, Int invpbaseval);
  void Redistribute(Idx * newVertexDist);
 
  //redistribute the graph according to the supernodal partition
  void RedistributeSupernodal(Int nsuper, Int * xsuper, Int * xsuperdist, Int * supMembership );
  void AllGatherStructure(SparseMatrixGraph & g);
  void GatherStructure(SparseMatrixGraph & g, int proot);
protected:
  void permute_(Int * invp, Idx * newVertexDist=nullptr, Int invpbaseval=1);
};

}

#endif //_DIST_SPARSE_MATRIX_GRAPH_HPP_

