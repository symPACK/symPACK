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
	vector<Ptr>  colptr;                          // Column index pointer
	vector<Idx>  rowind;                          // Starting row index pointer
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

  //Analysis related functions
//  void GetLColRowCount(ETree & tree, Ordering & aOrder, vector<Int> & cc, vector<Int> & rc);
//  void FindSupernodes(ETree& tree, Ordering & aOrder, vector<Int> & cc,vector<Int> & supMembership, vector<Int> & xsuper, Int maxSize = -1);
//#ifdef REFINED_SNODE
//  void RefineSupernodes(ETree& tree, Ordering & aOrder, vector<Int> & supMembership, vector<Int> & xsuper, PtrVec & xlindx, IdxVec & lindx, vector<Int> & perm);
//#endif
//  void RelaxSupernodes(ETree& tree, vector<Int> & cc,vector<Int> & supMembership, vector<Int> & xsuper, RelaxationParameters & params  );
//  void SymbolicFactorizationRelaxed(ETree& tree, Ordering & aOrder,const vector<Int> & cc,const vector<Int> & xsuper,const vector<Int> & SupMembership, PtrVec & xlindx, IdxVec & lindx);
//  void SymbolicFactorization(ETree& tree, Ordering & aOrder,const vector<Int> & cc,const vector<Int> & xsuper,const vector<Int> & SupMembership, PtrVec & xlindx, IdxVec & lindx);

};

}

#endif
