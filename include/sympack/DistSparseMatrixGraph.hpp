#ifndef _DIST_SPARSE_MATRIX_GRAPH_HPP_
#define _DIST_SPARSE_MATRIX_GRAPH_HPP_

#include "sympack/Environment.hpp"
#include "sympack/Ordering.hpp"
#include "sympack/Types.hpp"

namespace SYMPACK{
class ETree;
class Ordering;
class SparseMatrixStructure;
class DistSparseMatrixGraph;
class SparseMatrixGraph;

class SparseMatrixGraph{
  friend class Ordering;
  friend class DistSparseMatrixGraph;
  public:
  Idx          size;                            // Matrix dimension (global)
  Ptr          nnz;                             // Number of nonzeros (global)
  SYMPACK::vector<Ptr>  colptr;                 // Column index pointer
  SYMPACK::vector<Idx>  rowind;                 // Starting row index pointer
  protected:
  int baseval;
  int keepDiag;
  int sorted;
  bool bIsExpanded;

  public:
  void SetBaseval(int aBaseval);
  void SetKeepDiag(int aKeepDiag);
  void SetSorted(int aSorted);

  bool IsExpanded() const {return bIsExpanded;}
  int GetBaseval() const {return baseval;}
  int GetKeepDiag()const {return keepDiag;}
  int GetSorted() const {return sorted;}

  SparseMatrixGraph();
  void SortEdges();
  void DistributeGraph(DistSparseMatrixGraph & dg);

  Idx VertexCount() const { return colptr.size()==0?0:colptr.size()-1;}
  Ptr EdgeCount() const{ return rowind.size();}
};

class DistSparseMatrixGraph{
  friend class Ordering;


  protected:
  bool bIsExpanded;
  

  public:
	Idx          size;                            // Matrix dimension (global)
	Ptr          nnz;                             // Number of nonzeros (global)
	SYMPACK::vector<Ptr>  colptr;                 // Column index pointer
	SYMPACK::vector<Idx>  rowind;                 // Starting row index pointer
	SYMPACK::vector<Idx>  vertexDist;             // vertex distribution array (size np+1)

  protected:
  MPI_Comm comm;
  int mpirank;
  int mpisize;
  int baseval;
  int keepDiag;
  int sorted;


  public:
    void SetComm(const MPI_Comm & aComm);
    void SetBaseval(int aBaseval);
    void SetKeepDiag(int aKeepDiag);
    void SetSorted(int aSorted);

    MPI_Comm const & GetComm() const {return comm;}
    int GetBaseval() const {return baseval;}
    int GetKeepDiag()const {return keepDiag;}
    int GetSorted() const {return sorted;}

  public:
  DistSparseMatrixGraph();
  ~DistSparseMatrixGraph();
  DistSparseMatrixGraph( const DistSparseMatrixGraph& g );
  //constructor from global to local
  DistSparseMatrixGraph(const SparseMatrixStructure & A);

  DistSparseMatrixGraph& operator=( const DistSparseMatrixGraph& g );
  //accessors
  bool IsExpanded() const {return bIsExpanded;}
  Idx LocalFirstVertex() const {return vertexDist[mpirank];}
  Idx LocalVertexCount() const { return vertexDist[mpirank+1] - vertexDist[mpirank];}
  Ptr LocalEdgeCount() const{ return rowind.size();}

  //utility
  void FromStructure(const SparseMatrixStructure & A);
  void SortEdges();
  void ExpandSymmetric();
  void Permute(Int * invp);
  void Permute(Int * invp, Idx * newVertexDist);
  void Permute(Int * invp, Int invpbaseval);
  void Permute(Int * invp, Idx * newVertexDist, Int invpbaseval);
 
  //redistribute the graph according to the supernodal partition
  void RedistributeSupernodal(Int nsuper, Int * xsuper, Int * xsuperdist, Int * supMembership );
  void AllGatherStructure(SparseMatrixGraph & g);
  void GatherStructure(SparseMatrixGraph & g, int proot);
protected:
  void permute_(Int * invp, Idx * newVertexDist=NULL, Int invpbaseval=1);
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
