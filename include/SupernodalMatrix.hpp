#ifndef _SUPERNODAL_MATRIX_DECL_HPP_
#define _SUPERNODAL_MATRIX_DECL_HPP_

#include "Environment.hpp"
#include "SuperNode.hpp"

#include "NumVec.hpp"
#include "DistSparseMatrix.hpp"
#include "ETree.hpp"
#include "FBMatrix.hpp"



#include <vector>
#include "Mapping.hpp"



namespace LIBCHOLESKY{

template <typename T> class SupernodalMatrix{
  protected:
  bool globalAllocated = false;

  IntNumVec Xsuper_;
  IntNumVec SupMembership_;
  MAPCLASS Mapping_;
 
  SparseMatrixStructure Local_;
  SparseMatrixStructure Global_;

  IntNumVec UpdateCount_;

  ETree ETree_;
  ETree SupETree_;
  Int iSize_;
  std::vector<SuperNode<T> * > LocalSupernodes_;

  void GetUpdatingSupernodeCount(const IntNumVec & Xsuper, const IntNumVec & xlindx, const IntNumVec & lindx, const IntNumVec & SupMembership, IntNumVec & sc);
  inline bool FindNextUpdate(SuperNode<T> & src_snode, Int & src_nzblk_idx, Int & src_first_row,  Int & src_last_row, Int & tgt_snode_id);
  inline bool FindPivot(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode,Int & pivot_idx, Int & pivot_fr, Int & pivot_lr);

  void UpdateSuperNode(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode,Int & pivot_idx, Int  pivot_fr = I_ZERO);
  public:

	/// @brief MPI communicator
	//MPI_Comm     comm = MPI_COMM_NULL;        


  SupernodalMatrix();
  SupernodalMatrix(const DistSparseMatrix<T> & pMat, MAPCLASS & pMapping, MPI_Comm & pComm );

  //Accessors
  Int Size(){return iSize_;}
  IntNumVec & GetSupernodes(){ return Xsuper_;}
  std::vector<SuperNode<T> *  > & GetLocalSupernodes(){ return LocalSupernodes_; } 
  SuperNode<T> & GetLocalSupernode(Int i){ return *LocalSupernodes_[i]; } 
  ETree & GetETree(){return ETree_;}

  SparseMatrixStructure GetGlobalStructure();
  SparseMatrixStructure GetLocalStructure() const;

  void Factorize(MPI_Comm & pComm);

};



} // namespace LIBCHOLESKY

#include "SupernodalMatrix_impl.hpp"

#endif 
