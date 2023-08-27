#ifndef _SYMPACK_MATRIX_BASE_DECL_HPP_
#define _SYMPACK_MATRIX_BASE_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/DistSparseMatrix.hpp"
#include "sympack/ETree.hpp"
#include "sympack/Mapping.hpp"
#include "sympack/CommTypes.hpp"
#include "sympack/Ordering.hpp"

#include "sympack/supernodalTaskGraph.hpp"
#include "sympack/CommPull.hpp"
#include "sympack/Types.hpp"
#include "sympack/Task.hpp"
#include "sympack/Scheduler.hpp"

#include <list>
#include <tuple>
#include <deque>
#include <queue>
#include <vector>
#include <limits>
#include <numeric>
#include <tuple>

#include "SuperNode.hpp"
#include "SuperNodeInd.hpp"

#ifdef SP_THREADS
#if !UPCXX_BACKEND_GASNET_PAR && !UPCXX_THREADMODE
  // UPCXX_BACKEND_GASNET_PAR is an unspecified identifer in versions <= 2021.3.0
  // UPCXX_THREADMODE is the supported way to test this starting in 2020.10.0
  #error "UPCXX_THREADMODE=par required when using threads in symPACK."
#endif
#endif

namespace symPACK{

  //Forward declarations
  template<typename Task> class supernodalTaskGraph;

  class symPACKMatrixBase{
    protected:
      static int last_id;
    public:
      int sp_handle; 
      virtual ~symPACKMatrixBase(){
        sp_handle = last_id++;
      }
      //core functionalities
      virtual void Factorize() = 0;

      virtual void Init(symPACKOptions & options ) = 0;
  };

  template <typename T> class symPACKMatrixMeta: public symPACKMatrixBase{
    public:
      virtual void SymbolicFactorization(DistSparseMatrix<T> & pMat) = 0;
      virtual void DistributeMatrix(DistSparseMatrix<T> & pMat) = 0;
      //Solve routines
      //note: RHS & B are stored in column major format
      virtual void Solve(T * RHS, int nrhs, int rhs_size,  T * Xptr=nullptr) = 0;
      virtual void GetSolution(T * B, int nrhs) = 0;

      virtual void DumpMatlab() {};

      symPACKMatrixMeta():symPACKMatrixBase(),workteam_(nullptr){
#ifndef NO_MPI
        fullcomm_ = MPI_COMM_NULL;
        non_workcomm_ = MPI_COMM_NULL;
        workcomm_ = MPI_COMM_NULL;
#endif

      }

      ~symPACKMatrixMeta(){
#ifndef NO_MPI
        if(this->non_workcomm_ != MPI_COMM_NULL){
          MPI_Comm_free(&non_workcomm_);
        }

        if(this->workcomm_ != MPI_COMM_NULL){
          MPI_Comm_free(&this->workcomm_);
        }

        if(this->fullcomm_ != MPI_COMM_NULL){
          MPI_Comm_free(&this->fullcomm_);
        }
#endif
        if (workteam_!=nullptr) {
          workteam_->destroy();
        }
      }

      //Accessors
      Int Size(){return iSize_;}
      const ETree & GetETree(){return ETree_;}
      const Ordering & GetOrdering(){return Order_;}
      symPACKOptions GetOptions(){ return options_;}
      std::vector<Int> & GetSupernodalPartition(){ return Xsuper_;}
      const std::vector<Int> & GetSupMembership(){return SupMembership_;}

      Idx TotalSupernodeCnt() { return Xsuper_.empty()?0:Xsuper_.size()-1;}
    protected:
      //MPI/UPCXX ranks and sizes
      int iam, np,all_np;
      std::shared_ptr<RankGroup> group_;

      symPACKOptions options_;

      //Order of the matrix
      Int iSize_;
      //Column-based elimination tree
      ETree ETree_;
      //Column permutation
      Ordering Order_;

      //Local and Global structure of the matrix (CSC format)
      DistSparseMatrixGraph graph_;

      //Supernodal partition array: supernode I ranges from column Xsuper_[I-1] to Xsuper_[I]-1
      std::vector<Int> Xsuper_;
      std::vector<Int> XsuperDist_;
      //Supernode membership array: column i belongs to supernode SupMembership_[i-1]
      std::vector<Int> SupMembership_;

#ifndef NO_MPI
      MPI_Comm fullcomm_;
      MPI_Comm non_workcomm_;
      MPI_Comm workcomm_;
#endif
      std::unique_ptr<upcxx::team> workteam_;

      //CSC structure of L factor
      PtrVec locXlindx_;
      IdxVec locLindx_;

    protected:
      void getLColRowCount(SparseMatrixGraph & sgraph, std::vector<Int> & cc, std::vector<Int> & rc);
      void getLColRowCount(DistSparseMatrixGraph & dgraph, std::vector<Int> & cc, std::vector<Int> & rc);
      void findSupernodes(ETree& tree, Ordering & aOrder, std::vector<Int> & cc,std::vector<Int> & supMembership, std::vector<Int> & xsuper, Int maxSize = -1);
      void relaxSupernodes(ETree& tree, std::vector<Int> & cc,std::vector<Int> & supMembership, std::vector<Int> & xsuper, RelaxationParameters & params  );

      void symbolicFactorizationRelaxedDist(std::vector<Int> & cc);

      void refineSupernodes(int ordflag,int altflag,DistSparseMatrix<T>* pMat = nullptr);



      void gatherLStructure(std::vector<Ptr>& xlindx, std::vector<Idx> & lindx);


  };

} // namespace SYMPACK

#include <sympack/impl/symPACKMatrixBase_impl.hpp>


#endif //_SYMPACK_MATRIX_BASE_DECL_HPP_

