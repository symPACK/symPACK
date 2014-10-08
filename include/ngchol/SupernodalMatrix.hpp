#ifndef _SUPERNODAL_MATRIX_DECL_HPP_
#define _SUPERNODAL_MATRIX_DECL_HPP_

#include "ngchol/Environment.hpp"
#include "ngchol/SuperNode.hpp"

#include "ngchol/NumVec.hpp"
#include "ngchol/DistSparseMatrix.hpp"
#include "ngchol/ETree.hpp"
#include "ngchol/Mapping.hpp"
#include "ngchol/CommTypes.hpp"
#include "ngchol/Ordering.hpp"



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



namespace LIBCHOLESKY{


  Int DEFAULT_TARGET(Mapping * map,Int src, Int tgt){ return map->Map(tgt-1,tgt-1);}
  Int AGG_TARGET(Mapping * map,Int src, Int tgt){ return map->Map(tgt-1,tgt-1);}
  Int FACT_TARGET(Mapping * map,Int src, Int tgt){ return map->Map(tgt-1,src-1);}

  Int DEFAULT_TAG(Int src, Int tgt){ return (2*tgt);}
  Int AGG_TAG(Int src, Int tgt){ return (2*tgt+1);}
  Int FACT_TAG(Int src, Int tgt){ return (2*tgt);}

  Int AGG_TAG_TO_ID(Int tag){ return ((Int)(tag-1)/2);}
  Int FACT_TAG_TO_ID(Int tag){ return ((Int)(tag)/2);}

enum MappingType {ROW2D,COL2D,MODWRAP2D,MODWRAP2DNZ};
enum FactorizationType {FANOUT,FANBOTH};
struct NGCholOptions{
  MappingType mapping;
  FactorizationType factorization;
};

template<typename T> 
Int getAggBufSize(const SnodeUpdateFB & curTask, const IntNumVec & Xsuper, const IntNumVec & UpdateHeight){
  Int max_bytes = 5*sizeof(Int); 
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
Int getFactBufSize(const SnodeUpdateFB & curTask, const IntNumVec & UpdateWidth, const IntNumVec & UpdateHeight){
  Int max_bytes = 5*sizeof(Int); 
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
Int getMaxBufSize(const IntNumVec & UpdateWidth, const IntNumVec & UpdateHeight){
  Int max_bytes = 5*sizeof(Int); 
  //The upper bound must be of the width of the "largest" child
  Int nrows = 0;
  for(Int i = 0; i<UpdateHeight.m(); ++i ){ nrows = max(nrows, UpdateHeight[i]); }
  Int ncols = 0;
  for(Int i = 0; i<UpdateWidth.m(); ++i ){ ncols = max(ncols, UpdateWidth[i]); }

  Int nz_cnt = nrows * ncols;
  Int nblocks = nrows;
  max_bytes += (nblocks)*sizeof(NZBlockDesc);
  max_bytes += nz_cnt*sizeof(T);

  return max_bytes;
}

template <typename T> class SupernodalMatrix{


  public:

  //Constructors
  SupernodalMatrix();
  SupernodalMatrix(const DistSparseMatrix<T> & pMat, Int maxSnode,Mapping * pMapping, Int maxIsend, Int maxIrecv, MPI_Comm & pComm );
  //Destructor
  ~SupernodalMatrix();

  //Accessors
  Int Size(){return iSize_;}
  IntNumVec & GetSupernodalPartition(){ return Xsuper_;}

  ETree & GetETree(){return ETree_;}
  const Ordering & GetOrdering(){return Order_;}
  const IntNumVec & GetSupMembership(){return SupMembership_;}

  Int SupernodeCnt(){ return LocalSupernodes_.size(); } 
  std::vector<SuperNode<T> *  > & GetLocalSupernodes(){ return LocalSupernodes_; } 
  SuperNode<T> & GetLocalSupernode(Int i){ return *LocalSupernodes_[i]; } 

  SparseMatrixStructure GetGlobalStructure();
  SparseMatrixStructure GetLocalStructure() const;

  void Factorize();
  void FanOut( );
  void FanOutTask( );
  void FanBoth( );

  void FanOut2( );


#ifdef _CHECK_RESULT_SEQ_
  void Solve(NumMat<T> * RHS, NumMat<T> & forwardSol, NumMat<T> * Xptr=NULL);
#else
  void Solve(NumMat<T> * RHS, NumMat<T> * Xptr=NULL);
#endif

  void GetFullFactors( NumMat<T> & fullMatrix);
  void GetSolution(NumMat<T> & B);

  void Dump();


  protected:
  CommEnvironment * CommEnv_;


  //Is the global structure of the matrix allocated
  bool globalAllocated;

  //Local and Global structure of the matrix (CSC format)
  SparseMatrixStructure Local_;
  SparseMatrixStructure Global_;
  //CSC structure of L factor
  IntNumVec xlindx_;
  IntNumVec lindx_;


  //Supernodal partition array: supernode I ranges from column Xsuper_[I-1] to Xsuper_[I]-1
  IntNumVec Xsuper_;
  //Supernode membership array: column i belongs to supernode SupMembership_[i-1]
  IntNumVec SupMembership_;


  //MAPCLASS describing the Mapping of the computations
  Mapping * Mapping_;
 
  //Array storing the supernodal update count to a target supernode
  IntNumVec UpdateCount_;
  //Array storing the width of the widest supernode updating a target supernode
  IntNumVec UpdateWidth_;
  IntNumVec UpdateHeight_;

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
  std::vector<SuperNode<T> * > LocalSupernodes_;








    //Vector holding pointers to local contributions
    //This has to be renamed because it contains the distributed solution
    std::vector<SuperNode<T> *> Contributions_;


    
    ITree globToLocSnodes_;

#ifdef _SEPARATE_COMM_
    CommEnvironment * FBAggCommEnv_;
#endif

    AsyncComms outgoingSend;
    FBCommList MsgToSend;
    FBTasks LocalTasks;
    TempUpdateBuffers<T> tmpBufs;


  void FBAsyncRecv(Int iLocalI, std::vector<AsyncComms> & incomingRecvAggArr, std::vector<AsyncComms * > & incomingRecvFactArr, IntNumVec & AggregatesToRecv, IntNumVec & FactorsToRecv);
  void FBFactorizationTask(SnodeUpdateFB & curTask, Int iLocalI, IntNumVec & AggregatesDone, IntNumVec & AggregatesToRecv, std::vector<char> & src_blocks,std::vector<AsyncComms> & incomingRecvAggArr);
  void FBUpdateTask(SnodeUpdateFB & curTask, IntNumVec & UpdatesToDo, IntNumVec & AggregatesDone, std::vector< SuperNode<T> * > & aggVectors, std::vector<char> & src_blocks,std::vector<AsyncComms * > & incomingRecvFactArr, IntNumVec & FactorsToRecv);






  void Init(const DistSparseMatrix<T> & pMat, Int maxSnode,Mapping * pMapping, Int maxIsend, Int maxIrecv, MPI_Comm & pComm );
  


  //Do the packing of a FO update and enqueue it to outgoingSend
  void AddOutgoingComm(AsyncComms & outgoingSend, Int src_snode_id, Int src_snode_size ,Int src_first_row, NZBlockDesc & pivot_desc, Int nzblk_cnt, T * nzval_ptr, Int nz_cnt);
  void AddOutgoingComm(AsyncComms & outgoingSend, Icomm * send_buffer);
  //Wait for completion of some outgoing communication in outgoingSend
  void AdvanceOutgoing(AsyncComms & outgoingSend);


  //FanOut communication routines
  //AsyncRecvFactors tries to post some MPI_Irecv on the next updates (targetting current to next supernodes)
  void AsyncRecvFactors(Int iLocalI, std::vector<AsyncComms> & incomingRecvArr,IntNumVec & FactorsToRecv,IntNumVec & UpdatesToDo);
  void AsyncRecv(Int iLocalI, std::vector<AsyncComms> * incomingRecvFactArr, Int * FactorsToRecv, Int * UpdatesToDo);





  //WaitIncomingFactors returns an iterator to the first completed asynchronous factor receive. It must be called in a while loop.
  //Returns cur_incomingRecv.end() if everything has been received.
  inline AsyncComms::iterator WaitIncomingFactors(AsyncComms & cur_incomingRecv, MPI_Status & recv_status, AsyncComms & outgoingSend);




  void GetUpdatingSupernodeCount( IntNumVec & sc,IntNumVec & mw, IntNumVec & mh);



  inline bool FindNextUpdate(SuperNode<T> & src_snode, Int & tgt_snode_id, Int & f_ur, Int & f_ub, Int & n_ur, Int & n_ub);

  //FanOut related routines



  void SendMessage(const FBDelayedComm & comm, AsyncComms & OutgoingSend);

  void SendMessage(const DelayedComm & comm, AsyncComms & OutgoingSend, std::vector<SuperNode<T> *> & snodeColl);


  void SendDelayedMessages(Int cur_snode_id, CommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode<T> *> & snodeColl, bool reverse=false);

  void SendDelayedMessagesUp(Int cur_snode_id, CommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode<T> *> & snodeColl);
  void SendDelayedMessagesDown(Int iLocalI, DownCommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode<T> *> & snodeColl);
#ifdef SINGLE_BLAS
  inline void UpdateSuperNode(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode,Int & pivot_idx, NumMat<T> & tmpBuf,IntNumVec & src_colindx, IntNumVec & src_rowindx, IntNumVec & src_to_tgt_offset
, Int  pivot_fr = I_ZERO);
#else
  inline void UpdateSuperNode(SuperNode<T> & src_snode, SuperNode<T> & tgt_snode,Int & pivot_idx, Int  pivot_fr = I_ZERO);
#endif

  void SendDelayedMessagesUp(FBCommList & MsgToSend, AsyncComms & OutgoingSend, const SnodeUpdateFB * nextTask);
  void SendDelayedMessagesUp(FBCommList & MsgToSend, AsyncComms & OutgoingSend, FBTasks & taskList);




  //FanBoth related routines
  Int FBUpdate(Int I,Int prevJ=-1);
  void FBGetUpdateCount(IntNumVec & UpdatesToDo, IntNumVec & AggregatesToRecv);
#ifdef _SEPARATE_COMM_
  SuperNode<T> * FBRecvFactor(const SnodeUpdateFB & curTask, std::vector<char> & src_blocks,AsyncComms * cur_incomingRecv,AsyncComms::iterator & it, IntNumVec & FactorsToRecv, Int & recv_tgt_id);
#else
  SuperNode<T> * FBRecvFactor(const SnodeUpdateFB & curTask, std::vector<char> & src_blocks,AsyncComms * cur_incomingRecv,AsyncComms::iterator & it, IntNumVec & FactorsToRecv);
#endif

  //Solve related routines
  void forward_update(SuperNode<T> * src_contrib,SuperNode<T> * tgt_contrib);
  void back_update(SuperNode<T> * src_contrib,SuperNode<T> * tgt_contrib);


protected:
  //Supernodal elimination tree //deprecated
  ETree SupETree_;

#if 0
  //deprecated
  #ifdef UPDATE_LIST
  inline void FindUpdates(SuperNode<T> & src_snode, std::list<SnodeUpdateOld> & updates  );
  #endif
  inline bool FindNextUpdateOld(SuperNode<T> & src_snode, Int & src_nzblk_idx, Int & src_first_row,  Int & src_last_row, Int & tgt_snode_id);
#endif







};



} // namespace LIBCHOLESKY

#include "ngchol/SupernodalMatrix_impl.hpp"


#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) TAU_FSTART(a);
#define TIMER_STOP(a) TAU_FSTOP(a);
#endif
#endif



#endif 
