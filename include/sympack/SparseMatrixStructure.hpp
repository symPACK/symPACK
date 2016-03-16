#ifndef _DIST_SPARSE_MATRIX_STRUCTURE_HPP_
#define _DIST_SPARSE_MATRIX_STRUCTURE_HPP_

#include "sympack/Environment.hpp"
#include "sympack/Ordering.hpp"
#include "sympack/Types.hpp"

namespace SYMPACK{
class ETree;
class Ordering;
class DistSparseMatrixGraph;

class SparseMatrixStructure{
  friend class Ordering;
  friend class DistSparseMatrixGraph;


  protected:
  bool bIsGlobal;
  bool bIsExpanded;
  

  public:
  int baseval;
  int keepDiag;
  int sorted;
	Int          size;                            // Matrix dimension
	Int          nnz;                             // Number of nonzeros
	vector<Ptr>  colptr;                          // Column index pointer
	vector<Idx>  rowind;                          // Starting row index pointer

	vector<Ptr>  expColptr;                          // Column index pointer expanded
	vector<Idx>  expRowind;                          // Starting row index pointer expanded

  SparseMatrixStructure();

  //Distributed graph to global SparseMatrixStructure
  SparseMatrixStructure( const DistSparseMatrixGraph & G );


  bool IsExpanded() const { return bIsExpanded;}


  void ClearExpandedSymmetric();

  void ExpandSymmetric(MPI_Comm * workcomm = NULL);

  void ToGlobal(SparseMatrixStructure & pGlobal,MPI_Comm & comm);


  void GetLColRowCount(ETree & tree, Ordering & aOrder, vector<Int> & cc, vector<Int> & rc);
  void FindSupernodes(ETree& tree, Ordering & aOrder, vector<Int> & cc,vector<Int> & supMembership, vector<Int> & xsuper, Int maxSize = -1);

#ifdef REFINED_SNODE
  void RefineSupernodes(ETree& tree, Ordering & aOrder, vector<Int> & supMembership, vector<Int> & xsuper, PtrVec & xlindx, IdxVec & lindx, vector<Int> & perm);
#endif

  void RelaxSupernodes(ETree& tree, vector<Int> & cc,vector<Int> & supMembership, vector<Int> & xsuper, RelaxationParameters & params  );
  void SymbolicFactorizationRelaxed(ETree& tree, Ordering & aOrder,const vector<Int> & cc,const vector<Int> & xsuper,const vector<Int> & SupMembership, PtrVec & xlindx, IdxVec & lindx);

  void SymbolicFactorization(ETree& tree, Ordering & aOrder,const vector<Int> & cc,const vector<Int> & xsuper,const vector<Int> & SupMembership, PtrVec & xlindx, IdxVec & lindx);







//TRASH
//  void GetLColRowCountDEPRECATED(ETree & tree, vector<Int> & cc, vector<Int> & rc);
//  void SymbolicFactorizationDEPRECATED(ETree& tree,const vector<Int> & cc,const vector<Int> & xsuper, vector<Int> & xlindx, vector<Int> & lindx);
//  void GetARowStruct(const ETree & etree, const Int iPORow, std::vector<Int> & rowStruct);
//  void GetLRowStruct(const ETree & etree, const Int iPORow, const std::vector<Int> & ARowStruct, std::set<Int> & LRowStruct);
//
//  void GetSuperARowStruct(const ETree & etree, const vector<Int> & Xsuper, const vector<Int> & SupMembership, const Int iSupNo, std::vector<Int> & SuperRowStruct);
//  void GetSuperLRowStruct(const ETree & etree, const vector<Int> & Xsuper, const vector<Int> & SupMembership, const Int iSupNo, std::set<Int> & SuperLRowStruct);
//    


};


//class LocalSparseMatrixStructure: public SparseMatrixStructure{
//  public:
//  void ToGlobal(GlobalSparseMatrixStructure & pGlobal);
//};
//
//class GlobalSparseMatrixStructure: public SparseMatrixStructure{
//  public:
//  void GetLColRowCount(ETree & tree, vector<Int> & cc, vector<Int> & rc);
//  void FindSupernodes(ETree& tree, vector<Int> & cc, vector<Int> & xsuper);
//  void SymbolicFactorization(ETree& tree,const vector<Int> & cc,const vector<Int> & xsuper, vector<Int> & xlindx, vector<Int> & lindx);
//};




}

#endif
