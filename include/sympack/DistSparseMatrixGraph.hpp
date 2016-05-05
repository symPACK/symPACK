#ifndef _DIST_SPARSE_MATRIX_GRAPH_HPP_
#define _DIST_SPARSE_MATRIX_GRAPH_HPP_

#include "sympack/Environment.hpp"
#include "sympack/Ordering.hpp"
#include "sympack/Types.hpp"

namespace SYMPACK{
class ETree;
class Ordering;
class SparseMatrixStructure;

class DistSparseMatrixGraph{
  friend class Ordering;


  protected:
  bool bIsExpanded;
  

  public:
  int baseval;
  int keepDiag;
  int sorted;
	Idx          size;                            // Matrix dimension (global)
	Ptr          nnz;                             // Number of nonzeros (local)
	SYMPACK::vector<Ptr>  colptr;                          // Column index pointer
	SYMPACK::vector<Idx>  rowind;                          // Starting row index pointer
  MPI_Comm comm;

  DistSparseMatrixGraph();
  //constructor from global to local
  DistSparseMatrixGraph(const SparseMatrixStructure & A);

  //accessors
  bool IsExpanded() const {return bIsExpanded;}
  Idx LocalVertexCount() const { return colptr.size()-1;}
  Ptr LocalEdgeCount() const{ return rowind.size();}

  //utility
  void FromStructure(const SparseMatrixStructure & A);
  void SortEdges();
  void ExpandSymmetric();
  void Permute(Int * invp);
  //redistribute the graph according to the supernodal partition
  void RedistributeSupernodal(Int nsuper, Int * xsuper, Int * supMembership );

//  void FindSupernodes(ETree& tree, Ordering & aOrder, SYMPACK::vector<Int> & cc,SYMPACK::vector<Int> & supMembership, SYMPACK::vector<Int> & xsuper, Int maxSize = -1);
//
//  void RelaxSupernodes(ETree& tree, SYMPACK::vector<Int> & cc,SYMPACK::vector<Int> & supMembership, SYMPACK::vector<Int> & xsuper, RelaxationParameters & params  );
//
  //Analysis related functions
//  void GetLColRowCount(ETree & tree, Ordering & aOrder, SYMPACK::vector<Int> & cc, SYMPACK::vector<Int> & rc);
//#ifdef REFINED_SNODE
//  void RefineSupernodes(ETree& tree, Ordering & aOrder, SYMPACK::vector<Int> & supMembership, SYMPACK::vector<Int> & xsuper, PtrVec & xlindx, IdxVec & lindx, SYMPACK::vector<Int> & perm);
//#endif
//  void RelaxSupernodes(ETree& tree, SYMPACK::vector<Int> & cc,SYMPACK::vector<Int> & supMembership, SYMPACK::vector<Int> & xsuper, RelaxationParameters & params  );
//  void SymbolicFactorizationRelaxed(ETree& tree, Ordering & aOrder,const SYMPACK::vector<Int> & cc,const SYMPACK::vector<Int> & xsuper,const SYMPACK::vector<Int> & SupMembership, PtrVec & xlindx, IdxVec & lindx);
//  void SymbolicFactorization(ETree& tree, Ordering & aOrder,const SYMPACK::vector<Int> & cc,const SYMPACK::vector<Int> & xsuper,const SYMPACK::vector<Int> & SupMembership, PtrVec & xlindx, IdxVec & lindx);

};

}

#endif
