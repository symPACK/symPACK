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
#ifndef _SYMPACK_MATRIX_DECL_HPP_
#define _SYMPACK_MATRIX_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/DistSparseMatrix.hpp"
#include "sympack/ETree.hpp"
#include "sympack/Mapping.hpp"
#include "sympack/CommTypes.hpp"
#include "sympack/Ordering.hpp"

#include "sympack/CommPull.hpp"
#include "sympack/Types.hpp"
#include "sympack/Task.hpp"
#include "sympack/Scheduler.hpp"

#include <upcxx.h>

#include <list>
#include <tuple>
#include <deque>
#include <queue>
#include <vector>
#include <limits>
#include <numeric>
#include <tuple>
#include <memory>

#include "SuperNode.hpp"
#include "SuperNodeInd.hpp"



//#define PREFETCH_STRUCTURE

namespace symPACK{
  template<typename T> class symPACKMatrix;

  class supernodalTaskGraph{
    template<typename T> friend class symPACKMatrix;
    protected:
    std::vector<std::list<FBTask> * > taskLists_;
    Int localTaskCount_;

    public:

    supernodalTaskGraph( );
    supernodalTaskGraph( const supernodalTaskGraph& g );
    supernodalTaskGraph& operator=( const supernodalTaskGraph& g );
    ~supernodalTaskGraph();

    void removeTask(std::list<FBTask>::iterator & taskit); 
    std::list<FBTask>::iterator addTask(FBTask & task);
    Int getTaskCount();
    Int setTaskCount(Int value);
    Int increaseTaskCount();
    Int decreaseTaskCount();
    std::list<FBTask>::iterator find_task(Int src, Int tgt, TaskType type );
  };



  template <typename T> class symPACKMatrix{


    public:

      //Constructors
      symPACKMatrix();
      symPACKMatrix(DistSparseMatrix<T> & pMat, symPACKOptions & options );
      //TODO
      symPACKMatrix( symPACKMatrix & M){};
      //Destructor
      ~symPACKMatrix();

      //operators
      //TODO
      symPACKMatrix & operator=( symPACKMatrix & M){return M;};


      void Init(symPACKOptions & options );
      void Init(DistSparseMatrix<T> & pMat, symPACKOptions & options );

      void SymbolicFactorization(DistSparseMatrix<T> & pMat);
      void DistributeMatrix(DistSparseMatrix<T> & pMat);

      //Accessors
      Int Size(){return iSize_;}
      Int LocalSupernodeCnt(){ return LocalSupernodes_.size(); } 
      std::vector<Int> & GetSupernodalPartition(){ return Xsuper_;}
      const ETree & GetETree(){return ETree_;}
      const Ordering & GetOrdering(){return Order_;}
      const Mapping * GetMapping(){return Mapping_;}
      symPACKOptions GetOptions(){ return options_;}
      const std::vector<Int> & GetSupMembership(){return SupMembership_;}
      std::vector<SuperNode<T> *  > & GetLocalSupernodes(){ return LocalSupernodes_; }
      //TODO Check if that's useful
      SuperNode<T> & GetLocalSupernode(Int i){ return *LocalSupernodes_[i]; } 


      //core functionalities
      void Factorize();

      //Solve routines
      //note: RHS & B are stored in column major format
      void Solve(T * RHS, int nrhs,  T * Xptr=NULL);
      void GetSolution(T * B, int nrhs);

      void FanBoth( );


      //debug routines
      void DumpMatlab();
      void Dump();
      void DumpContrib();

      Idx TotalSupernodeCnt() { return Xsuper_.empty()?0:Xsuper_.size()-1;}








      /******************* Global to Local Indexes utility routines ******************/
      //returns the 1-based index of supernode id global in the local supernode array
      Int snodeLocalIndex(Int global);
      //returns a reference to  a local supernode with id global
      SuperNode<T> * snodeLocal(Int global);
      template< class Alloc>
        SuperNode<T,Alloc> * snodeLocal(Int global, std::vector<SuperNode<T,Alloc> *> & snodeColl);

      template <class Allocator = UpcxxAllocator>
        SuperNode<T,Allocator> * CreateSuperNode(DecompositionType type);
      template <class Allocator = UpcxxAllocator>
        SuperNode<T,Allocator> * CreateSuperNode(DecompositionType type,Int aiId, Int aiFc, Int aiLc, Int ai_num_rows, Int aiN, Int aiNZBlkCnt=-1);
      template <class Allocator = UpcxxAllocator>
        SuperNode<T,Allocator> * CreateSuperNode(DecompositionType type,char * dataPtr,size_t size, Int firstRow = -1);
      template <class Allocator = UpcxxAllocator>
        SuperNode<T,Allocator> * CreateSuperNode(DecompositionType type,Int aiId, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices);


    protected:
      //MPI/UPCXX ranks and sizes
      int iam, np,all_np;
      std::shared_ptr<RankGroup> group_;

      symPACKOptions options_;
      CommEnvironment * CommEnv_;
      MPI_Comm non_workcomm_;
      MPI_Comm fullcomm_;

      //Order of the matrix
      Int iSize_;
      //Column-based elimination tree
      ETree ETree_;
      //Column permutation
      Ordering Order_;
      //MAPCLASS describing the Mapping of the computations
      Mapping * Mapping_;
      LoadBalancer * Balancer_;

      //Local and Global structure of the matrix (CSC format)
      DistSparseMatrixGraph graph_;

      //CSC structure of L factor
      //      PtrVec xlindx_;
      //      IdxVec lindx_;
      PtrVec locXlindx_;
      IdxVec locLindx_;

      std::vector<Int> numBlk_;
      std::vector<Int> cc_,rc_;


      BackupBuffer backupBuffer_;

      //Supernodal partition array: supernode I ranges from column Xsuper_[I-1] to Xsuper_[I]-1
      std::vector<Int> Xsuper_;
      std::vector<Int> XsuperDist_;
      //Supernode membership array: column i belongs to supernode SupMembership_[i-1]
      std::vector<Int> SupMembership_;

      //TODO Task lists
      Scheduler<std::list<FBTask>::iterator> * scheduler_;
      Scheduler<FBTask> * scheduler2_;

      //backup for factorization
      std::vector<Int> UpdatesToDo_;
      std::vector<std::list<FBTask> * > origTaskLists_;
      Int localTaskCount_;
      void generateTaskGraph(Int & localTaskCount, std::vector<std::list<FBTask> * > & taskLists);
      std::list<FBTask>::iterator find_task(std::vector<std::list<FBTask> * > & taskLists, Int src, Int tgt, TaskType type );

      supernodalTaskGraph taskGraph_;
      void generateTaskGraph(supernodalTaskGraph & taskGraph,std::vector<Int> & AggregatesToRecv,  std::vector<Int>& LocalAggregates);

      std::vector<std::list<Int> > chSupTree_;
      void dfs_traversal(std::vector<std::list<Int> > & tree,int node,std::list<Int> & frontier);


      std::list<std::list<FBTask>::iterator > readyTasks_;

      //Array storing the supernodal update count to a target supernode
      std::vector<Int> UpdateCount_;
      //Array storing the width of the widest supernode updating a target supernode
      std::vector<Int> UpdateWidth_;
      std::vector<Int> UpdateHeight_;


      //This has to be moved to an option structure
      Int maxIsend_;
      Int maxIrecv_;
      Int incomingRecvCnt_;


      //Vector holding pointers to local SuperNode2 objects (L factor)
      std::vector<SuperNode<T> * > LocalSupernodes_;
      //std::vector<upcxx::global_ptr<SuperNodeDesc > > remoteFactors_;
      std::vector< std::tuple< upcxx::global_ptr<SuperNodeDesc >,Int> > remoteFactors_;



      //Vector holding pointers to local contributions
      //This has to be renamed because it contains the distributed solution
      std::vector<SuperNode<T,MallocAllocator> *> Contributions_;


#ifndef ITREE2
      std::vector<Int> globToLocSnodes_;
#else
      ITree globToLocSnodes_;
#endif


      TempUpdateBuffers<T> tmpBufs;

      void findSupernodes(ETree& tree, Ordering & aOrder, std::vector<Int> & cc,std::vector<Int> & supMembership, std::vector<Int> & xsuper, Int maxSize = -1);
      void relaxSupernodes(ETree& tree, std::vector<Int> & cc,std::vector<Int> & supMembership, std::vector<Int> & xsuper, RelaxationParameters & params  );

      void symbolicFactorizationRelaxedDist(std::vector<Int> & cc);

      void refineSupernodes(int ordflag,int altflag,DistSparseMatrix<T>* pMat = NULL);

      void getLColRowCount(std::vector<Int> & cc, std::vector<Int> & rc);
      void getLColRowCount(SparseMatrixGraph & sgraph, std::vector<Int> & cc, std::vector<Int> & rc);


    protected:


      //FanBoth related routines
      Int FBUpdate(Int I,Int prevJ=-1);
      void FBGetUpdateCount(std::vector<Int> & UpdatesToDo, std::vector<Int> & AggregatesToRecv, std::vector<Int> & LocalAggregates);
      void GetUpdatingSupernodeCount( std::vector<Int> & sc,std::vector<Int> & mw, std::vector<Int> & mh, std::vector<Int> & numBlk);

      void FBFactorizationTask(supernodalTaskGraph & taskGraph, FBTask & curTask, Int iLocalI, std::vector< SuperNode<T> * > & aggVectors, bool is_static = false);
      void FBAggregationTask(supernodalTaskGraph & taskGraph, FBTask & curTask, Int iLocalI, bool is_static = false);
      void FBUpdateTask(supernodalTaskGraph & taskGraph, FBTask & curTask, std::vector<Int> & UpdatesToDo, std::vector< SuperNode<T> * > & aggVectors, bool is_static = false);


      //Communication related routines
      void AddOutgoingComm(AsyncComms & outgoingSend, Icomm * send_buffer);
      void AdvanceOutgoing(AsyncComms & outgoingSend);
      template< class Alloc>
        void SendDelayedMessagesUp(Int cur_snode_id, CommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode<T,Alloc> *> & snodeColl);
      template< class Alloc>
        void SendDelayedMessagesDown(Int iLocalI, DownCommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode<T,Alloc> *> & snodeColl);
      template< class Alloc>
        void SendMessage(const DelayedComm & comm, AsyncComms & OutgoingSend, std::vector<SuperNode<T,Alloc> *> & snodeColl);

      void CheckIncomingMessages(supernodalTaskGraph & taskGraph, std::vector< SuperNode<T> * > & aggVectors, bool is_static = false);


  };

} // namespace SYMPACK

#include <sympack/symPACKMatrix_impl.hpp>




#endif //_SYMPACK_MATRIX_DECL_HPP_

