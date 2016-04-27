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
	SYMPACK::vector<Ptr>  colptr;                          // Column index pointer
	SYMPACK::vector<Idx>  rowind;                          // Starting row index pointer

	SYMPACK::vector<Ptr>  expColptr;                          // Column index pointer expanded
	SYMPACK::vector<Idx>  expRowind;                          // Starting row index pointer expanded

  SparseMatrixStructure();

  //Distributed graph to global SparseMatrixStructure
  SparseMatrixStructure( const DistSparseMatrixGraph & G );


  bool IsExpanded() const { return bIsExpanded;}


  void ClearExpandedSymmetric();

  void ExpandSymmetric(MPI_Comm * workcomm = NULL);

  void ToGlobal(SparseMatrixStructure & pGlobal,MPI_Comm & comm);


  void GetLColRowCount(ETree & tree, Ordering & aOrder, SYMPACK::vector<Int> & cc, SYMPACK::vector<Int> & rc);
  void FindSupernodes(ETree& tree, Ordering & aOrder, SYMPACK::vector<Int> & cc,SYMPACK::vector<Int> & supMembership, SYMPACK::vector<Int> & xsuper, Int maxSize = -1);

#ifdef REFINED_SNODE
  void RefineSupernodes(ETree& tree, Ordering & aOrder, SYMPACK::vector<Int> & supMembership, SYMPACK::vector<Int> & xsuper, PtrVec & xlindx, IdxVec & lindx, SYMPACK::vector<Int> & perm);
#endif

  void RelaxSupernodes(ETree& tree, SYMPACK::vector<Int> & cc,SYMPACK::vector<Int> & supMembership, SYMPACK::vector<Int> & xsuper, RelaxationParameters & params  );
  void SymbolicFactorizationRelaxedDist(ETree& tree, Ordering & aOrder,const SYMPACK::vector<Int> & cc,const SYMPACK::vector<Int> & xsuper,const SYMPACK::vector<Int> & SupMembership, PtrVec & xlindx, IdxVec & lindx,MPI_Comm & comm);
  void SymbolicFactorizationRelaxed(ETree& tree, Ordering & aOrder,const SYMPACK::vector<Int> & cc,const SYMPACK::vector<Int> & xsuper,const SYMPACK::vector<Int> & SupMembership, PtrVec & xlindx, IdxVec & lindx);
  //void SymbolicFactorizationRelaxed_upcxx(ETree& tree,Ordering & aOrder, const SYMPACK::vector<Int> & cc,const SYMPACK::vector<Int> & xsuper,const SYMPACK::vector<Int> & SupMembership, upcxx::shared_array<Ptr> & xlindx, upcxx::shared_array<Idx> & lindx);

  void SymbolicFactorization(ETree& tree, Ordering & aOrder,const SYMPACK::vector<Int> & cc,const SYMPACK::vector<Int> & xsuper,const SYMPACK::vector<Int> & SupMembership, PtrVec & xlindx, IdxVec & lindx);







//TRASH
//  void GetLColRowCountDEPRECATED(ETree & tree, SYMPACK::vector<Int> & cc, SYMPACK::vector<Int> & rc);
//  void SymbolicFactorizationDEPRECATED(ETree& tree,const SYMPACK::vector<Int> & cc,const SYMPACK::vector<Int> & xsuper, SYMPACK::vector<Int> & xlindx, SYMPACK::vector<Int> & lindx);
//  void GetARowStruct(const ETree & etree, const Int iPORow, SYMPACK::vector<Int> & rowStruct);
//  void GetLRowStruct(const ETree & etree, const Int iPORow, const SYMPACK::vector<Int> & ARowStruct, std::set<Int> & LRowStruct);
//
//  void GetSuperARowStruct(const ETree & etree, const SYMPACK::vector<Int> & Xsuper, const SYMPACK::vector<Int> & SupMembership, const Int iSupNo, SYMPACK::vector<Int> & SuperRowStruct);
//  void GetSuperLRowStruct(const ETree & etree, const SYMPACK::vector<Int> & Xsuper, const SYMPACK::vector<Int> & SupMembership, const Int iSupNo, std::set<Int> & SuperLRowStruct);
//    


};


//class LocalSparseMatrixStructure: public SparseMatrixStructure{
//  public:
//  void ToGlobal(GlobalSparseMatrixStructure & pGlobal);
//};
//
//class GlobalSparseMatrixStructure: public SparseMatrixStructure{
//  public:
//  void GetLColRowCount(ETree & tree, SYMPACK::vector<Int> & cc, SYMPACK::vector<Int> & rc);
//  void FindSupernodes(ETree& tree, SYMPACK::vector<Int> & cc, SYMPACK::vector<Int> & xsuper);
//  void SymbolicFactorization(ETree& tree,const SYMPACK::vector<Int> & cc,const SYMPACK::vector<Int> & xsuper, SYMPACK::vector<Int> & xlindx, SYMPACK::vector<Int> & lindx);
//};




}

#endif
