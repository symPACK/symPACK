#ifndef _SUPERNODAL2_MATRIX_DECL_HPP_
#define _SUPERNODAL2_MATRIX_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/SupernodalMatrixBase.hpp"
#include "sympack/SuperNode2.hpp"

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
#include <deque>
#include <queue>
#include <vector>

#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) 
#define TIMER_STOP(a) 
#define scope_timer(b,a)
#endif
#endif



namespace SYMPACK{
  template <typename T> class SupernodalMatrix2: public SupernodalMatrixBase{


    public:
      upcxx::team * team_;

      //Constructors
      SupernodalMatrix2();
      SupernodalMatrix2(const DistSparseMatrix<T> & pMat, NGCholOptions & options );
      //TODO
      SupernodalMatrix2( SupernodalMatrix2 & M){};
      //Destructor
      ~SupernodalMatrix2();

      //operators
      //TODO
      SupernodalMatrix2 & operator=( SupernodalMatrix2 & M){return M;};


      void Init(const DistSparseMatrix<T> & pMat, NGCholOptions & options );

      //Accessors
      Int Size(){return iSize_;}
      Int SupernodeCnt(){ return LocalSupernodes_.size(); } 
      SYMPACK::vector<Int> & GetSupernodalPartition(){ return Xsuper_;}
      const ETree & GetETree(){return ETree_;}
      const Ordering & GetOrdering(){return Order_;}
      const SYMPACK::vector<Int> & GetSupMembership(){return SupMembership_;}
      SYMPACK::vector<SuperNode2<T> *  > & GetLocalSupernodes(){ return LocalSupernodes_; } 
      //TODO Check if that's useful
      SuperNode2<T> & GetLocalSupernode(Int i){ return *LocalSupernodes_[i]; } 

      SparseMatrixStructure GetGlobalStructure();
      SparseMatrixStructure GetLocalStructure() const;

      //core functionalities
      void Factorize();
      void Solve(NumMat<T> * RHS, NumMat<T> * Xptr=NULL);
      void GetSolution(NumMat<T> & B);

      void FanBoth( );
      void FanBoth_Static( );


      void Dump();


    protected:
      NGCholOptions options_;
      CommEnvironment * CommEnv_;

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
      SparseMatrixStructure * Local_;
      SparseMatrixStructure * Global_;
      //Is the global structure of the matrix allocated
      bool isGlobStructAllocated_;

      //CSC structure of L factor
      PtrVec xlindx_;
      IdxVec lindx_;
      bool isXlindxAllocated_;
      bool isLindxAllocated_;

      BackupBuffer backupBuffer_;

      //Supernodal partition array: supernode I ranges from column Xsuper_[I-1] to Xsuper_[I]-1
      SYMPACK::vector<Int> Xsuper_;
      //Supernode membership array: column i belongs to supernode SupMembership_[i-1]
      SYMPACK::vector<Int> SupMembership_;

      //TODO Task lists
      Scheduler<std::list<FBTask>::iterator> * scheduler_;
      SYMPACK::vector<std::list<FBTask> * > taskLists_;
      std::list<std::list<FBTask>::iterator > readyTasks_;
      std::list<FBTask>::iterator find_task(Int src, Int tgt, TaskType type );

      //Array storing the supernodal update count to a target supernode
      SYMPACK::vector<Int> UpdateCount_;
      //Array storing the width of the widest supernode updating a target supernode
      SYMPACK::vector<Int> UpdateWidth_;
      SYMPACK::vector<Int> UpdateHeight_;


      //This has to be moved to an option structure
      Int maxIsend_;
      Int maxIrecv_;
      Int incomingRecvCnt_;


      //Vector holding pointers to local SuperNode2 objects (L factor)
      SYMPACK::vector<SuperNode2<T> * > LocalSupernodes_;



      //Vector holding pointers to local contributions
      //This has to be renamed because it contains the distributed solution
      SYMPACK::vector<SuperNode2<T> *> Contributions_;


#ifndef ITREE2
      SYMPACK::vector<Int> globToLocSnodes_;
#else
      ITree globToLocSnodes_;
#endif


      TempUpdateBuffers<T> tmpBufs;


      Int localTaskCount_;


      protected:


      /******************* Global to Local Indexes utility routines ******************/
      //returns the 1-based index of supernode id global in the local supernode array
      Int snodeLocalIndex(Int global);
      //returns a reference to  a local supernode with id global
      SuperNode2<T> * snodeLocal(Int global);
      SuperNode2<T> * snodeLocal(Int global, SYMPACK::vector<SuperNode2<T> *> & snodeColl);




      void GetUpdatingSupernodeCount( SYMPACK::vector<Int> & sc,SYMPACK::vector<Int> & mw, SYMPACK::vector<Int> & mh);
      inline bool FindNextUpdate(SuperNode2<T> & src_snode, Int & tgt_snode_id, Int & f_ur, Int & f_ub, Int & n_ur, Int & n_ub);

      //FanBoth related routines
      Int FBUpdate(Int I,Int prevJ=-1);
      void FBGetUpdateCount(SYMPACK::vector<Int> & UpdatesToDo, SYMPACK::vector<Int> & AggregatesToRecv, SYMPACK::vector<Int> & LocalAggregates);

      void FBFactorizationTask(FBTask & curTask, Int iLocalI, bool is_static = false);
      void FBAggregationTask(FBTask & curTask, Int iLocalI, bool is_static = false);
      void FBUpdateTask(FBTask & curTask, SYMPACK::vector<Int> & UpdatesToDo, SYMPACK::vector<Int> & AggregatesDone,SYMPACK::vector< SuperNode2<T> * > & aggVectors,  SYMPACK::vector<Int> & FactorsToRecv, SYMPACK::vector<Int> & AggregatesToRecv,Int & localTaskCount, bool is_static = false);

      //Solve related routines
      void forward_update(SuperNode2<T> * src_contrib,SuperNode2<T> * tgt_contrib);
      void back_update(SuperNode2<T> * src_contrib,SuperNode2<T> * tgt_contrib);


      //Communication related routines
      void AddOutgoingComm(AsyncComms & outgoingSend, Icomm * send_buffer);
      void AdvanceOutgoing(AsyncComms & outgoingSend);
      void SendDelayedMessagesUp(Int cur_snode_id, CommList & MsgToSend, AsyncComms & OutgoingSend, SYMPACK::vector<SuperNode2<T> *> & snodeColl);
      void SendDelayedMessagesDown(Int iLocalI, DownCommList & MsgToSend, AsyncComms & OutgoingSend, SYMPACK::vector<SuperNode2<T> *> & snodeColl);
      void SendMessage(const DelayedComm & comm, AsyncComms & OutgoingSend, SYMPACK::vector<SuperNode2<T> *> & snodeColl);

      void CheckIncomingMessages(bool is_static = false);


  };
} // namespace SYMPACK

#include "sympack/SupernodalMatrix2_impl.hpp"


#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) TAU_FSTART(a);
#define TIMER_STOP(a) TAU_FSTOP(a);
#define scope_timer(b,a) CTF_scope_timer b(#a)
#endif
#endif



#endif 
