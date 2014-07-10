#ifndef _DIST_SPARSE_MATRIX_STRUCTURE_HPP_
#define _DIST_SPARSE_MATRIX_STRUCTURE_HPP_

#include "ngchol/Environment.hpp"
#include "ngchol/NumVec.hpp"

namespace LIBCHOLESKY{

class ETree;

class SparseMatrixStructure{
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

  //Permutation to have the ordered matrix
  // Computed by the Ordering class
  IntNumVec    perm;
  IntNumVec    invp;


  void ClearExpandedSymmetric();

  void ExpandSymmetric();

  void ToGlobal(SparseMatrixStructure & pGlobal);


  void GetLColRowCount(ETree & tree, IntNumVec & cc, IntNumVec & rc);
  void FindSupernodes(ETree& tree, IntNumVec & cc,IntNumVec & supMembership, IntNumVec & xsuper, Int maxSize = -1);

#ifdef REFINED_SNODE
  void RefineSupernodes(ETree& tree, IntNumVec & supMembership, IntNumVec & xsuper, IntNumVec & xlindx, IntNumVec & lindx, IntNumVec & perm);
#endif

#ifdef RELAXED_SNODE
  void RelaxSupernodes(ETree& tree, IntNumVec & cc,IntNumVec & supMembership, IntNumVec & xsuper, Int maxSize );
  void SymbolicFactorizationRelaxed(ETree& tree,const IntNumVec & cc,const IntNumVec & xsuper,const IntNumVec & SupMembership, IntNumVec & xlindx, IntNumVec & lindx);
#endif

  void SymbolicFactorization(ETree& tree,const IntNumVec & cc,const IntNumVec & xsuper,const IntNumVec & SupMembership, IntNumVec & xlindx, IntNumVec & lindx);





  SparseMatrixStructure();




//TODO make an ordering class friend with this class
void MMD(IntNumVec & pperm, IntNumVec & pinvp);

  inline Int Perm(Int i){
    if (perm.m()!=size){
      perm.Resize(size);
      for(Int i = 0; i<size;++i){
        perm[i]=i+1;
      }
    }
    return perm[i];
  }
  inline Int Invp(Int i){
    if (invp.m()!=size){
      invp.Resize(size);
      for(Int i = 0; i<size;++i){
        invp[i]=i+1;
      }
    }
    return invp[i];
  }







void Permute(IntNumVec & perm);


//TRASH
  void GetLColRowCountDEPRECATED(ETree & tree, IntNumVec & cc, IntNumVec & rc);
  void SymbolicFactorizationDEPRECATED(ETree& tree,const IntNumVec & cc,const IntNumVec & xsuper, IntNumVec & xlindx, IntNumVec & lindx);
  void GetARowStruct(const ETree & etree, const Int iPORow, std::vector<Int> & rowStruct);
  void GetLRowStruct(const ETree & etree, const Int iPORow, const std::vector<Int> & ARowStruct, std::set<Int> & LRowStruct);

  void GetSuperARowStruct(const ETree & etree, const IntNumVec & Xsuper, const IntNumVec & SupMembership, const Int iSupNo, std::vector<Int> & SuperRowStruct);
  void GetSuperLRowStruct(const ETree & etree, const IntNumVec & Xsuper, const IntNumVec & SupMembership, const Int iSupNo, std::set<Int> & SuperLRowStruct);
    


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
