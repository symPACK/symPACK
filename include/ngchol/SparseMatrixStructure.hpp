#ifndef _DIST_SPARSE_MATRIX_STRUCTURE_HPP_
#define _DIST_SPARSE_MATRIX_STRUCTURE_HPP_

#include "ngchol/Environment.hpp"
#include "ngchol/NumVec.hpp"
#include "ngchol/Ordering.hpp"

namespace LIBCHOLESKY{
class ETree;
class Ordering;

class SparseMatrixStructure{
  friend class Ordering;


  protected:
  bool bIsGlobal;
  bool bIsExpanded;

  public:
	Int          size;                            // Matrix dimension
	Int          nnz;                             // Number of nonzeros
	IntNumVec    colptr;                          // Column index pointer
	IntNumVec    rowind;                          // Starting row index pointer

	IntNumVec    expColptr;                          // Column index pointer expanded
	IntNumVec    expRowind;                          // Starting row index pointer expanded


  void ClearExpandedSymmetric();

  void ExpandSymmetric();

  void ToGlobal(SparseMatrixStructure & pGlobal,MPI_Comm & comm);


  void GetLColRowCount(ETree & tree, Ordering & aOrder, IntNumVec & cc, IntNumVec & rc);
  void FindSupernodes(ETree& tree, Ordering & aOrder, IntNumVec & cc,IntNumVec & supMembership, IntNumVec & xsuper, Int maxSize = -1);

#ifdef REFINED_SNODE
  void RefineSupernodes(ETree& tree, Ordering & aOrder, IntNumVec & supMembership, IntNumVec & xsuper, PtrVec & xlindx, IdxVec & lindx, IntNumVec & perm);
#endif

#ifdef RELAXED_SNODE
  void RelaxSupernodes(ETree& tree, IntNumVec & cc,IntNumVec & supMembership, IntNumVec & xsuper, RelaxationParameters & params  );
  void SymbolicFactorizationRelaxed(ETree& tree, Ordering & aOrder,const IntNumVec & cc,const IntNumVec & xsuper,const IntNumVec & SupMembership, PtrVec & xlindx, IdxVec & lindx);
#endif

  void SymbolicFactorization(ETree& tree, Ordering & aOrder,const IntNumVec & cc,const IntNumVec & xsuper,const IntNumVec & SupMembership, PtrVec & xlindx, IdxVec & lindx);





  SparseMatrixStructure();


//TRASH
//  void GetLColRowCountDEPRECATED(ETree & tree, IntNumVec & cc, IntNumVec & rc);
//  void SymbolicFactorizationDEPRECATED(ETree& tree,const IntNumVec & cc,const IntNumVec & xsuper, IntNumVec & xlindx, IntNumVec & lindx);
//  void GetARowStruct(const ETree & etree, const Int iPORow, std::vector<Int> & rowStruct);
//  void GetLRowStruct(const ETree & etree, const Int iPORow, const std::vector<Int> & ARowStruct, std::set<Int> & LRowStruct);
//
//  void GetSuperARowStruct(const ETree & etree, const IntNumVec & Xsuper, const IntNumVec & SupMembership, const Int iSupNo, std::vector<Int> & SuperRowStruct);
//  void GetSuperLRowStruct(const ETree & etree, const IntNumVec & Xsuper, const IntNumVec & SupMembership, const Int iSupNo, std::set<Int> & SuperLRowStruct);
//    


};


//class LocalSparseMatrixStructure: public SparseMatrixStructure{
//  public:
//  void ToGlobal(GlobalSparseMatrixStructure & pGlobal);
//};
//
//class GlobalSparseMatrixStructure: public SparseMatrixStructure{
//  public:
//  void GetLColRowCount(ETree & tree, IntNumVec & cc, IntNumVec & rc);
//  void FindSupernodes(ETree& tree, IntNumVec & cc, IntNumVec & xsuper);
//  void SymbolicFactorization(ETree& tree,const IntNumVec & cc,const IntNumVec & xsuper, IntNumVec & xlindx, IntNumVec & lindx);
//};




}

#endif
