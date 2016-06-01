#ifndef _SUPERNODAL_MATRIX_DECL_HPP_
#define _SUPERNODAL_MATRIX_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/SupernodalMatrixBase.hpp"
#include "sympack/SuperNode.hpp"

#include "sympack/DistSparseMatrix.hpp"
#include "sympack/ETree.hpp"
#include "sympack/Mapping.hpp"
#include "sympack/CommTypes.hpp"
#include "sympack/Ordering.hpp"

#include "sympack/Types.hpp"
#include "sympack/Task.hpp"

#include <list>
#include <deque>
#include <queue>
#include <vector>

#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) 
#define TIMER_STOP(a) 
#endif
#endif



namespace SYMPACK{
  inline Int DEFAULT_TARGET(Mapping * map,Int src, Int tgt){ return map->Map(tgt-1,tgt-1);}
  inline Int AGG_TARGET(Mapping * map,Int src, Int tgt){ return map->Map(tgt-1,tgt-1);}
  inline Int FACT_TARGET(Mapping * map,Int src, Int tgt){ return map->Map(tgt-1,src-1);}

  inline Int DEFAULT_TAG(Int src, Int tgt){ return (2*tgt);}
  inline Int AGG_TAG(Int src, Int tgt){ return (2*tgt+1);}
  inline Int FACT_TAG(Int src, Int tgt){ return (2*tgt);}

  inline Int AGG_TAG_TO_ID(Int tag){ return ((Int)(tag-1)/2);}
  inline Int FACT_TAG_TO_ID(Int tag){ return ((Int)(tag)/2);}


  template<typename T> 
    Int getAggBufSize(const SnodeUpdateFB & curTask, const SYMPACK::vector<Int> & Xsuper, const SYMPACK::vector<Int> & UpdateHeight){
      //Int max_bytes = 6*sizeof(Int); 
      Int max_bytes = sizeof(SuperNodeDesc);//6*sizeof(Int); 
      //The upper bound must be of the width of the "largest" child
      Int nrows = UpdateHeight[curTask.tgt_snode_id-1];
      Int ncols = Xsuper[curTask.tgt_snode_id] - Xsuper[curTask.tgt_snode_id-1];
      Int nz_cnt = nrows * ncols;
      Int nblocks = nrows;
      max_bytes += (nblocks)*sizeof(NZBlockDesc);
      max_bytes += nz_cnt*sizeof(T);
      //extra int to store the number of updates within the aggregate
      //max_bytes += sizeof(Int); 

      return max_bytes;
    }

  template<typename T> 
    Int getFactBufSize(const SnodeUpdateFB & curTask, const SYMPACK::vector<Int> & UpdateWidth, const SYMPACK::vector<Int> & UpdateHeight){
      Int max_bytes = sizeof(SuperNodeDesc);//6*sizeof(Int); 
      //The upper bound must be of the width of the "largest" child
      Int nrows = UpdateHeight[curTask.tgt_snode_id-1];
      Int ncols = UpdateWidth[curTask.tgt_snode_id-1];
      Int nz_cnt = nrows * ncols;
      Int nblocks = nrows;
      max_bytes += (nblocks)*sizeof(NZBlockDesc);
      max_bytes += nz_cnt*sizeof(T);

      return max_bytes;
    } 

  template<typename T> 
    Int getMaxBufSize(const SYMPACK::vector<Int> & UpdateWidth, const SYMPACK::vector<Int> & UpdateHeight){
      Int max_bytes = sizeof(SuperNodeDesc);//6*sizeof(Int); 
      //Int max_bytes = 6*sizeof(Int); 
      //The upper bound must be of the width of the "largest" child
      Int nrows = 0;
      for(Int i = 0; i<UpdateHeight.size(); ++i ){ nrows = max(nrows, UpdateHeight[i]); }
      Int ncols = 0;
      for(Int i = 0; i<UpdateWidth.size(); ++i ){ ncols = max(ncols, UpdateWidth[i]); }

      Int nz_cnt = nrows * ncols;
      Int nblocks = nrows;
      max_bytes += (nblocks)*sizeof(NZBlockDesc);
      max_bytes += nz_cnt*sizeof(T);

      return max_bytes;
    }

  template <typename T> class SupernodalMatrix: public SupernodalMatrixBase{


    public:

      //Constructors
      SupernodalMatrix();
      SupernodalMatrix(const DistSparseMatrix<T> & pMat, NGCholOptions & options );
      //Destructor
      ~SupernodalMatrix();

      //Accessors
      Int Size(){return iSize_;}
      SYMPACK::vector<Int> & GetSupernodalPartition(){ return Xsuper_;}

      ETree & GetETree(){return ETree_;}
      const Ordering & GetOrdering(){return Order_;}
      const SYMPACK::vector<Int> & GetSupMembership(){return SupMembership_;}

      Int SupernodeCnt(){ return LocalSupernodes_.size(); } 
      SYMPACK::vector<SuperNode<T, MallocAllocator> *  > & GetLocalSupernodes(){ return LocalSupernodes_; } 
      SuperNode<T, MallocAllocator> & GetLocalSupernode(Int i){ return *LocalSupernodes_[i]; } 

      SparseMatrixStructure GetGlobalStructure();
      SparseMatrixStructure GetLocalStructure() const;

      void Factorize();
      void FanOut( );
      void FanOutTask( );
      void FanBoth( );

      void FanOut2( );




      void FanBothPull( );
      void FBPullAsyncRecv(Int iLocalI, SYMPACK::vector<AsyncComms> & incomingRecvAggArr, SYMPACK::vector<AsyncComms * > & incomingRecvFactArr, SYMPACK::vector<Int> & AggregatesToRecv, SYMPACK::vector<Int> & FactorsToRecv);
      void FBPullFactorizationTask(SnodeUpdateFB & curTask, Int iLocalI, SYMPACK::vector<Int> & AggregatesDone, SYMPACK::vector<Int> & AggregatesToRecv, SYMPACK::vector<char> & src_blocks,SYMPACK::vector<AsyncComms> & incomingRecvAggArr);
      void FBPullUpdateTask(SnodeUpdateFB & curTask, SYMPACK::vector<Int> & UpdatesToDo, SYMPACK::vector<Int> & AggregatesDone, SYMPACK::vector< SuperNode<T, MallocAllocator> * > & aggVectors, SYMPACK::vector<char> & src_blocks,SYMPACK::vector<AsyncComms * > & incomingRecvFactArr, SYMPACK::vector<Int> & FactorsToRecv);

#ifdef _SEPARATE_COMM_
      SuperNode<T, MallocAllocator> * FBPullRecvFactor(const SnodeUpdateFB & curTask, SYMPACK::vector<char> & src_blocks,AsyncComms * cur_incomingRecv,AsyncComms::iterator & it, SYMPACK::vector<Int> & FactorsToRecv, Int & recv_tgt_id);
#else
      SuperNode<T, MallocAllocator> * FBPullRecvFactor(const SnodeUpdateFB & curTask, SYMPACK::vector<char> & src_blocks,AsyncComms * cur_incomingRecv,AsyncComms::iterator & it, SYMPACK::vector<Int> & FactorsToRecv);
#endif












#ifdef _CHECK_RESULT_SEQ_
      void Solve(NumMat<T> * RHS, NumMat<T> & forwardSol, NumMat<T> * Xptr=NULL);
#else
      void Solve(NumMat<T> * RHS, NumMat<T> * Xptr=NULL);
#endif

      void GetFullFactors( NumMat<T> & fullMatrix);
      void GetSolution(NumMat<T> & B);

      void Dump();


    protected:
      NGCholOptions options_;
      CommEnvironment * CommEnv_;


      //Is the global structure of the matrix allocated
      bool globalAllocated;

      //Local and Global structure of the matrix (CSC format)
      SparseMatrixStructure Local_;
      SparseMatrixStructure * Global_;
      //CSC structure of L factor
      PtrVec xlindx_;
      IdxVec lindx_;


      //Supernodal partition array: supernode I ranges from column Xsuper_[I-1] to Xsuper_[I]-1
      SYMPACK::vector<Int> Xsuper_;
      //Supernode membership array: column i belongs to supernode SupMembership_[i-1]
      SYMPACK::vector<Int> SupMembership_;

      //TODO Task lists
      SYMPACK::vector<std::list<SnodeUpdateFB> * > taskLists_;
      std::list<SnodeUpdateFB*> readyTasks_;
      
      //MAPCLASS describing the Mapping of the computations
      Mapping * Mapping_;
      LoadBalancer * Balancer_;



      //Array storing the supernodal update count to a target supernode
      SYMPACK::vector<Int> UpdateCount_;
      //Array storing the width of the widest supernode updating a target supernode
      SYMPACK::vector<Int> UpdateWidth_;
      SYMPACK::vector<Int> UpdateHeight_;

      //Column-based elimination tree
      ETree ETree_;

      //Column permutation
      Ordering Order_;

      //Order of the matrix
      Int iSize_;

      //This has to be moved to an option structure
      Int maxIsend_;
      Int maxIrecv_;
      Int incomingRecvCnt_;


      //Vector holding pointers to local SuperNode objects (L factor)
      SYMPACK::vector<SuperNode<T, MallocAllocator> * > LocalSupernodes_;








      //Vector holding pointers to local contributions
      //This has to be renamed because it contains the distributed solution
      SYMPACK::vector<SuperNode<T, MallocAllocator> *> Contributions_;



#ifndef ITREE2
      SYMPACK::vector<Int> globToLocSnodes_;
#else
      ITree globToLocSnodes_;
#endif

      /******************* Global to Local Indexes utility routines ******************/
      //returns the 1-based index of supernode id global in the local supernode array
      Int snodeLocalIndex(Int global);
      //returns a reference to  a local supernode with id global
      SuperNode<T, MallocAllocator> * snodeLocal(Int global);
      SuperNode<T, MallocAllocator> * snodeLocal(Int global, SYMPACK::vector<SuperNode<T, MallocAllocator> *> & snodeColl);



#ifdef _SEPARATE_COMM_
      CommEnvironment * FBAggCommEnv_;
#endif

      AsyncComms outgoingSend;
      FBCommList MsgToSend;
      FBTasks LocalTasks;
      TempUpdateBuffers<T> tmpBufs;

#ifdef _STAT_COMM_
      size_t maxAggreg_ ;
      size_t sizesAggreg_;
      Int countAggreg_;
      size_t maxFactors_ ;
      size_t sizesFactors_;
      Int countFactors_;

      size_t maxAggregRecv_ ;
      size_t sizesAggregRecv_;
      Int countAggregRecv_;
      size_t maxFactorsRecv_ ;
      size_t sizesFactorsRecv_;
      Int countFactorsRecv_;
#endif

      Int FBTaskAsyncRecv(Int iLocalI, SnodeUpdateFB & curTask, SYMPACK::vector<AsyncComms> & incomingRecvAggArr, SYMPACK::vector<AsyncComms * > & incomingRecvFactArr, SYMPACK::vector<Int> & AggregatesToRecv, SYMPACK::vector<Int> & FactorsToRecv);
      void FBAsyncRecv(Int iLocalI, SYMPACK::vector<AsyncComms> & incomingRecvAggArr, SYMPACK::vector<AsyncComms * > & incomingRecvFactArr, SYMPACK::vector<Int> & AggregatesToRecv, SYMPACK::vector<Int> & FactorsToRecv);
      void FBFactorizationTask(SnodeUpdateFB & curTask, Int iLocalI, SYMPACK::vector<Int> & AggregatesDone, SYMPACK::vector<Int> & FactorsToRecv, SYMPACK::vector<Int> & AggregatesToRecv, SYMPACK::vector<char> & src_blocks,SYMPACK::vector<AsyncComms> & incomingRecvAggArr, SYMPACK::vector<AsyncComms * > & incomingRecvFactArr);
      void FBUpdateTask(SnodeUpdateFB & curTask, SYMPACK::vector<Int> & UpdatesToDo, SYMPACK::vector<Int> & AggregatesDone, SYMPACK::vector< SuperNode<T, MallocAllocator> * > & aggVectors, SYMPACK::vector<char> & src_blocks,SYMPACK::vector<AsyncComms> & incomingRecvAggArr, SYMPACK::vector<AsyncComms * > & incomingRecvFactArr, SYMPACK::vector<Int> & FactorsToRecv, SYMPACK::vector<Int> & AggregatesToRecv);






      void Init(const DistSparseMatrix<T> & pMat, NGCholOptions & options );



      //Do the packing of a FO update and enqueue it to outgoingSend
      void AddOutgoingComm(AsyncComms & outgoingSend, Int src_snode_id, Int src_snode_size ,Int src_first_row, NZBlockDesc & pivot_desc, Int nzblk_cnt, T * nzval_ptr, Int nz_cnt);
      void AddOutgoingComm(AsyncComms & outgoingSend, Icomm * send_buffer);
      //Wait for completion of some outgoing communication in outgoingSend
      void AdvanceOutgoing(AsyncComms & outgoingSend);


      //FanOut communication routines
      //AsyncRecvFactors tries to post some MPI_Irecv on the next updates (targetting current to next supernodes)
      void AsyncRecvFactors(Int iLocalI, SYMPACK::vector<AsyncComms> & incomingRecvArr,SYMPACK::vector<Int> & FactorsToRecv,SYMPACK::vector<Int> & UpdatesToDo);
      void AsyncRecv(Int iLocalI, SYMPACK::vector<AsyncComms> * incomingRecvFactArr, Int * FactorsToRecv, Int * UpdatesToDo);





      //WaitIncomingFactors returns an iterator to the first completed asynchronous factor receive. It must be called in a while loop.
      //Returns cur_incomingRecv.end() if everything has been received.
      inline AsyncComms::iterator WaitIncomingFactors(AsyncComms & cur_incomingRecv, MPI_Status & recv_status, AsyncComms & outgoingSend);




      void GetUpdatingSupernodeCount( SYMPACK::vector<Int> & sc,SYMPACK::vector<Int> & mw, SYMPACK::vector<Int> & mh);



      inline bool FindNextUpdate(SuperNode<T, MallocAllocator> & src_snode, Int & tgt_snode_id, Int & f_ur, Int & f_ub, Int & n_ur, Int & n_ub);

      //FanOut related routines



      void SendMessage(const FBDelayedComm & comm, AsyncComms & OutgoingSend);

      void SendMessage(const DelayedComm & comm, AsyncComms & OutgoingSend, SYMPACK::vector<SuperNode<T, MallocAllocator> *> & snodeColl);


      void SendDelayedMessages(Int cur_snode_id, CommList & MsgToSend, AsyncComms & OutgoingSend, SYMPACK::vector<SuperNode<T, MallocAllocator> *> & snodeColl, bool reverse=false);

      void SendDelayedMessagesUp(Int cur_snode_id, CommList & MsgToSend, AsyncComms & OutgoingSend, SYMPACK::vector<SuperNode<T, MallocAllocator> *> & snodeColl);
      void SendDelayedMessagesDown(Int iLocalI, DownCommList & MsgToSend, AsyncComms & OutgoingSend, SYMPACK::vector<SuperNode<T, MallocAllocator> *> & snodeColl);
#ifdef SINGLE_BLAS
      inline void UpdateSuperNode(SuperNode<T, MallocAllocator> & src_snode, SuperNode<T, MallocAllocator> & tgt_snode,Int & pivot_idx, NumMat<T> & tmpBuf,SYMPACK::vector<Int> & src_colindx, SYMPACK::vector<Int> & src_rowindx, SYMPACK::vector<Int> & src_to_tgt_offset
          , Int  pivot_fr = I_ZERO);
#else
      inline void UpdateSuperNode(SuperNode<T, MallocAllocator> & src_snode, SuperNode<T, MallocAllocator> & tgt_snode,Int & pivot_idx, Int  pivot_fr = I_ZERO);
#endif

      void SendDelayedMessagesUp(FBCommList & MsgToSend, AsyncComms & OutgoingSend, const SnodeUpdateFB * nextTask);
      void SendDelayedMessagesUp(FBCommList & MsgToSend, AsyncComms & OutgoingSend, FBTasks & taskList);




      //FanBoth related routines
      Int FBUpdate(Int I,Int prevJ=-1);
      void FBGetUpdateCount(SYMPACK::vector<Int> & UpdatesToDo, SYMPACK::vector<Int> & AggregatesToRecv);
#ifdef _SEPARATE_COMM_
      SuperNode<T, MallocAllocator> * FBRecvFactor(const SnodeUpdateFB & curTask, SYMPACK::vector<char> & src_blocks,AsyncComms * cur_incomingRecv,AsyncComms::iterator & it, SYMPACK::vector<Int> & FactorsToRecv, Int & recv_tgt_id);
#else
      SuperNode<T, MallocAllocator> * FBRecvFactor(const SnodeUpdateFB & curTask, SYMPACK::vector<char> & src_blocks,AsyncComms * cur_incomingRecv,AsyncComms::iterator & it, SYMPACK::vector<Int> & FactorsToRecv);
#endif

      //Solve related routines
      void forward_update(SuperNode<T, MallocAllocator> * src_contrib,SuperNode<T, MallocAllocator> * tgt_contrib);
      void back_update(SuperNode<T, MallocAllocator> * src_contrib,SuperNode<T, MallocAllocator> * tgt_contrib);


    protected:


#ifdef PREALLOC_IRECV
      //IRecv buffers
      list<Icomm *> availRecvBuffers_;
      list<Icomm *> usedRecvBuffers_;
#endif      



#if 0
      //Supernodal elimination tree //deprecated
      ETree SupETree_;
      //deprecated
#ifdef UPDATE_LIST
      inline void FindUpdates(SuperNode<T> & src_snode, std::list<SnodeUpdateOld> & updates  );
#endif
      inline bool FindNextUpdateOld(SuperNode<T> & src_snode, Int & src_nzblk_idx, Int & src_first_row,  Int & src_last_row, Int & tgt_snode_id);
#endif







  };

} // namespace SYMPACK

#include "sympack/SupernodalMatrix_impl.hpp"


#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) TAU_FSTART(a);
#define TIMER_STOP(a) TAU_FSTOP(a);
#endif
#endif



#endif 
