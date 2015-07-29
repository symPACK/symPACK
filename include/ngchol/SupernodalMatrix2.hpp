#ifndef _SUPERNODAL2_MATRIX_DECL_HPP_
#define _SUPERNODAL2_MATRIX_DECL_HPP_

#include "ngchol/Environment.hpp"
#include "ngchol/SupernodalMatrixBase.hpp"
#include "ngchol/SuperNode2.hpp"

#include "ngchol/NumVec.hpp"
#include "ngchol/DistSparseMatrix.hpp"
#include "ngchol/ETree.hpp"
#include "ngchol/Mapping.hpp"
#include "ngchol/CommTypes.hpp"
#include "ngchol/Ordering.hpp"

#include "ngchol/CommPull.hpp"
#include "ngchol/Types.hpp"
#include "ngchol/Task.hpp"
#include "ngchol/Scheduler.hpp"

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



namespace LIBCHOLESKY{
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



      //Accessors
      Int Size(){return iSize_;}
      Int SupernodeCnt(){ return LocalSupernodes_.size(); } 
      IntNumVec & GetSupernodalPartition(){ return Xsuper_;}
      const ETree & GetETree(){return ETree_;}
      const Ordering & GetOrdering(){return Order_;}
      const IntNumVec & GetSupMembership(){return SupMembership_;}
      std::vector<SuperNode2<T> *  > & GetLocalSupernodes(){ return LocalSupernodes_; } 
      //TODO Check if that's useful
      SuperNode2<T> & GetLocalSupernode(Int i){ return *LocalSupernodes_[i]; } 

      SparseMatrixStructure GetGlobalStructure();
      SparseMatrixStructure GetLocalStructure() const;

      //core functionalities
      void Factorize();
      void Solve(NumMat<T> * RHS, NumMat<T> * Xptr=NULL);
      void GetSolution(NumMat<T> & B);

      void FanBoth( );


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


      //Supernodal partition array: supernode I ranges from column Xsuper_[I-1] to Xsuper_[I]-1
      IntNumVec Xsuper_;
      //Supernode membership array: column i belongs to supernode SupMembership_[i-1]
      IntNumVec SupMembership_;

      //TODO Task lists
      Scheduler<std::list<FBTask>::iterator> * scheduler_;
      std::vector<std::list<FBTask> * > taskLists_;
      std::list<std::list<FBTask>::iterator > readyTasks_;
      std::list<FBTask>::iterator find_task(Int src, Int tgt, TaskType type );

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
      std::vector<SuperNode2<T> * > LocalSupernodes_;



      //Vector holding pointers to local contributions
      //This has to be renamed because it contains the distributed solution
      std::vector<SuperNode2<T> *> Contributions_;


#ifndef ITREE2
      std::vector<Int> globToLocSnodes_;
#else
      ITree globToLocSnodes_;
#endif


      TempUpdateBuffers<T> tmpBufs;


      Int localTaskCount_;


      protected:

      void Init(const DistSparseMatrix<T> & pMat, NGCholOptions & options );

      /******************* Global to Local Indexes utility routines ******************/
      //returns the 1-based index of supernode id global in the local supernode array
      Int snodeLocalIndex(Int global);
      //returns a reference to  a local supernode with id global
      SuperNode2<T> * snodeLocal(Int global);
      SuperNode2<T> * snodeLocal(Int global, std::vector<SuperNode2<T> *> & snodeColl);




      void GetUpdatingSupernodeCount( std::vector<Int> & sc,std::vector<Int> & mw, std::vector<Int> & mh);
      inline bool FindNextUpdate(SuperNode2<T> & src_snode, Int & tgt_snode_id, Int & f_ur, Int & f_ub, Int & n_ur, Int & n_ub);

      //FanBoth related routines
      Int FBUpdate(Int I,Int prevJ=-1);
      void FBGetUpdateCount(std::vector<Int> & UpdatesToDo, std::vector<Int> & AggregatesToRecv, std::vector<Int> & LocalAggregates);

      void FBFactorizationTask(FBTask & curTask, Int iLocalI);
      void FBAggregationTask(FBTask & curTask, Int iLocalI);
      void FBUpdateTask(FBTask & curTask, std::vector<Int> & UpdatesToDo, std::vector<Int> & AggregatesDone,std::vector< SuperNode2<T> * > & aggVectors,  std::vector<Int> & FactorsToRecv, std::vector<Int> & AggregatesToRecv,Int & localTaskCount);

      //Solve related routines
      void forward_update(SuperNode2<T> * src_contrib,SuperNode2<T> * tgt_contrib);
      void back_update(SuperNode2<T> * src_contrib,SuperNode2<T> * tgt_contrib);


      //Communication related routines
      void AddOutgoingComm(AsyncComms & outgoingSend, Icomm * send_buffer);
      void AdvanceOutgoing(AsyncComms & outgoingSend);
      void SendDelayedMessagesUp(Int cur_snode_id, CommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode2<T> *> & snodeColl);
      void SendDelayedMessagesDown(Int iLocalI, DownCommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode2<T> *> & snodeColl);
      void SendMessage(const DelayedComm & comm, AsyncComms & OutgoingSend, std::vector<SuperNode2<T> *> & snodeColl);

      void CheckIncomingMessages();


  };
} // namespace LIBCHOLESKY

#include "ngchol/SupernodalMatrix2_impl.hpp"


#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) TAU_FSTART(a);
#define TIMER_STOP(a) TAU_FSTOP(a);
#define scope_timer(b,a) CTF_scope_timer b(#a)
#endif
#endif



#endif 
