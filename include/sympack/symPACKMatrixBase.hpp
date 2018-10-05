/*
   "symPACK" Copyright (c) 2016, The Regents of the University of California,
   through Lawrence Berkeley National Laboratory (subject to receipt of any
   required approvals from the U.S. Dept. of Energy).  All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   (1) Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

   (2) Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

   (3) Neither the name of the University of California, Lawrence Berkeley
   National Laboratory, U.S. Dept. of Energy nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
   Copyright (c) 2016 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.  

Author: Mathias Jacquelin

This file is part of symPACK. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to Lawrence Berkeley National
Laboratory, without imposing a separate written license agreement for such
Enhancements, then you hereby grant the following license: a non-exclusive,
royalty-free perpetual license to install, use, modify, prepare derivative
works, incorporate into other computer software, distribute, and sublicense
such enhancements or derivative works thereof, in binary and source code form.
 */
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
  };

  template <typename T> class symPACKMatrixMeta: public symPACKMatrixBase{
    public:
      symPACKMatrixMeta():symPACKMatrixBase(){
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

      //CSC structure of L factor
      //      PtrVec xlindx_;
      //      IdxVec lindx_;
      PtrVec locXlindx_;
      IdxVec locLindx_;

    protected:
      void getLColRowCount(SparseMatrixGraph & sgraph, std::vector<Int> & cc, std::vector<Int> & rc);
      void getLColRowCount(DistSparseMatrixGraph & dgraph, std::vector<Int> & cc, std::vector<Int> & rc);
      void findSupernodes(ETree& tree, Ordering & aOrder, std::vector<Int> & cc,std::vector<Int> & supMembership, std::vector<Int> & xsuper, Int maxSize = -1);
      void relaxSupernodes(ETree& tree, std::vector<Int> & cc,std::vector<Int> & supMembership, std::vector<Int> & xsuper, RelaxationParameters & params  );

      void symbolicFactorizationRelaxedDist(std::vector<Int> & cc);

      void refineSupernodes(int ordflag,int altflag,DistSparseMatrix<T>* pMat = NULL);



      void gatherLStructure(std::vector<Ptr>& xlindx, std::vector<Idx> & lindx);


  };

  namespace TSP{
    typedef struct SymbolBlok_ {
      int frownum;  /*< First row index            */
      int lrownum;  /*< Last row index (inclusive) */
      int lcblknm;  /*< Local column block         */
      int fcblknm;  /*< Facing column block        */ //>> CBLK que tu update (ou cblk couran pr le bloc diagonal)
    } SymbolBlok;

    typedef struct SymbolCblk_ {
      int fcolnum;  /*< First column index               */
      int lcolnum;  /*< Last column index (inclusive)    */
      int bloknum;  /*< First block in column (diagonal) */
      int brownum;  /*< First block in row facing the diagonal block in browtab, 0-based */
    } SymbolCblk;



    typedef struct SymbolMatrix_ {
      int            baseval;  /*< Base value for numberings               */
      int            dof;      /*< Degrees of freedom per node
                                 (constant if > 0, unconstant if 0 (not implemented)) */
      int            cblknbr;  /*< Number of column blocks                 */
      int            bloknbr;  /*< Number of blocks                        */
      int            nodenbr;  /*< Number of node in the compressed symbol */
      SymbolCblk   * cblktab;  /*< Array of column blocks [+1,based]       */
      SymbolBlok   * bloktab;  /*< Array of blocks [based]                 */
      int * browtab;  /*< Global index of a given block number in bloktab 0-based    */
    } SymbolMatrix;

    typedef struct Order_ {
      int  baseval;   /*< base value used for numbering       */
      int  vertnbr;   /*< Number of vertices                  */
      int  cblknbr;   /*< Number of column blocks             */
      int *permtab;   /*< Permutation array [based]           */ //scotch
      int *peritab;   /*< Inverse permutation array [based]   */
      int *rangtab;   /*< Column block range array [based,+1] */
      int *treetab;   /*< Partitioning tree [based]           */
    } Order;


    void countBlock(const vector<Int> & Xsuper, const vector<Ptr> & Xlindx, const vector<Idx> & Lindx, const vector<Int> & supMembership, vector<int> & blkCount);
    void symbolBuildRowtab(SymbolMatrix *symbptr);
    Order * GetPastixOrder(SymbolMatrix * symbptr,const vector<Int> & xsuper, const ETree & etree, const Int * perm, const Int * iperm);
    SymbolMatrix *  GetPastixSymbolMatrix(const vector<Int> & xsuper,const vector<Int> & supMembership, vector<Ptr> & xlindx, vector<Idx> & lindx);


    void symbolReordering( const SymbolMatrix *symbptr, Order *order, int split_level, int stop_criteria, int stop_when_fitting );
  }
} // namespace SYMPACK

#include <sympack/impl/symPACKMatrixBase_impl.hpp>


#endif //_SYMPACK_MATRIX_BASE_DECL_HPP_
