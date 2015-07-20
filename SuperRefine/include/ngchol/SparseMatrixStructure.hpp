#ifndef _DIST_SPARSE_MATRIX_STRUCTURE_HPP_
#define _DIST_SPARSE_MATRIX_STRUCTURE_HPP_

#include <vector>
//#include "ngchol/Environment.hpp"
//#include "ngchol/NumVec.hpp"
#include "ngchol/Ordering.hpp"
#include <mpi.h>

using namespace std;
namespace LIBCHOLESKY{

class ETree;
class Ordering;

class SparseMatrixStructure{
  friend class Ordering;


  protected:
  bool bIsGlobal;
  bool bIsExpanded;

  public:
	int          size;                            // Matrix dimension
	int          nnz;                             // Number of nonzeros
	vector<int>    colptr;                          // Column index pointer
	vector<int>    rowind;                          // Starting row index pointer

	vector<int>    expColptr;                          // Column index pointer expanded
	vector<int>    expRowind;                          // Starting row index pointer expanded


  void ClearExpandedSymmetric();

  void ExpandSymmetric();

  void ToGlobal(SparseMatrixStructure & pGlobal,MPI_Comm & comm);


  void GetLColRowCount(ETree & tree, Ordering & aOrder, vector<int> & cc, vector<int> & rc);
  void FindSupernodes(ETree& tree, Ordering & aOrder, vector<int> & cc,vector<int> & supMembership, vector<int> & xsuper, int maxSize = -1);
  void RefineSupernodes(ETree& tree, Ordering & aOrder, vector<int> & supMembership, vector<int> & xsuper, vector<int64_t> & xlindx, vector<int32_t> & lindx, vector<int> & perm, vector<int> & origPerm, vector<int> & newxsuper);

  void RelaxSupernodes(ETree& tree, vector<int> & cc,vector<int> & supMembership, vector<int> & xsuper, int maxSize );
  void SymbolicFactorizationRelaxed(ETree& tree, Ordering & aOrder,const vector<int> & cc,const vector<int> & xsuper,const vector<int> & SupMembership, vector<int64_t> & xlindx, vector<int32_t> & lindx);
  void SymbolicFactorization(ETree& tree, Ordering & aOrder,const vector<int> & cc,const vector<int> & xsuper,const vector<int> & SupMembership, vector<int64_t> & xlindx, vector<int32_t> & lindx);





  SparseMatrixStructure();


};



}

#endif
