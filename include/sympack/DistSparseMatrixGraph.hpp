#ifndef _DIST_SPARSE_MATRIX_GRAPH_HPP_
#define _DIST_SPARSE_MATRIX_GRAPH_HPP_

#include "sympack/Environment.hpp"
#include "sympack/Ordering.hpp"
#include "sympack/Types.hpp"

namespace symPACK{
class ETree;
class Ordering;
class SparseMatrixStructure;
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
  void ExpandSymmetric();
  void BroadcastGraph(MPI_Comm comm, Int root = 0);

  Idx VertexCount() const { return colptr.size()==0?0:colptr.size()-1;}
  Ptr EdgeCount() const{ return rowind.size();}
};

class DistSparseMatrixGraph{
  friend class Ordering;
  template <typename F> friend class DistSparseMatrix;


  protected:
  public:
  bool bIsExpanded;
  

  public:
	Idx          size;                            // Matrix dimension (global)
	Ptr          nnz;                             // Number of nonzeros (global)
	std::vector<Ptr>  colptr;                 // Column index pointer
	std::vector<Idx>  rowind;                 // Starting row index pointer
	std::vector<Idx>  vertexDist;             // vertex distribution array (size np+1)
  MPI_Comm comm;

  protected:
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
  Idx LocalVertexCount() const { Idx count = vertexDist[mpirank+1] - vertexDist[mpirank]; bassert(colptr.size()==count+1 || colptr.empty() ); return count; }
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
};

}

#endif
