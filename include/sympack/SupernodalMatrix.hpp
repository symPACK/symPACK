#ifndef _SUPERNODAL2_MATRIX_DECL_HPP_
#define _SUPERNODAL2_MATRIX_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/SupernodalMatrixBase.hpp"
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
#include <limits>
#include <numeric>
#include "SuperNode.hpp"
#include "SuperNodeInd.hpp"

#ifdef MULTITHREADING
#include <omp.h>
#endif

#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) 
#define TIMER_STOP(a) 
#define scope_timer(b,a)
#endif
#endif



namespace SYMPACK{
      template<typename T> class SupernodalMatrix;

      class supernodalTaskGraph{
        template<typename T> friend class SupernodalMatrix;
        protected:
        SYMPACK::vector<std::list<FBTask> * > taskLists_;
        Int localTaskCount_;
        
#ifdef MULTITHREADING
        SYMPACK::vector<omp_lock_t> superLocks_;
#endif
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



  template <typename T> class SupernodalMatrix: public SupernodalMatrixBase{


    public:
      upcxx::team * team_;

      //Constructors
      SupernodalMatrix();
      SupernodalMatrix(DistSparseMatrix<T> & pMat, NGCholOptions & options );
      //TODO
      SupernodalMatrix( SupernodalMatrix & M){};
      //Destructor
      ~SupernodalMatrix();

      //operators
      //TODO
      SupernodalMatrix & operator=( SupernodalMatrix & M){return M;};


      void Init(DistSparseMatrix<T> & pMat, NGCholOptions & options );

      //Accessors
      Int Size(){return iSize_;}
      Int LocalSupernodeCnt(){ return LocalSupernodes_.size(); } 
      SYMPACK::vector<Int> & GetSupernodalPartition(){ return Xsuper_;}
      const ETree & GetETree(){return ETree_;}
      const Ordering & GetOrdering(){return Order_;}
      const SYMPACK::vector<Int> & GetSupMembership(){return SupMembership_;}
      SYMPACK::vector<SuperNode<T> *  > & GetLocalSupernodes(){ return LocalSupernodes_; } 
      //TODO Check if that's useful
      SuperNode<T> & GetLocalSupernode(Int i){ return *LocalSupernodes_[i]; } 

      SparseMatrixStructure GetGlobalStructure();
      SparseMatrixStructure GetLocalStructure() const;

      //core functionalities
      void Factorize();

      //Solve routines
      //note: RHS & B are stored in column major format
      void Solve(T * RHS, int nrhs,  T * Xptr=NULL);
      void GetSolution(T * B, int nrhs);

      void FanBoth( );
      void FanBoth_Static( );


      void Dump();
      void DumpContrib();

      Idx TotalSupernodeCnt() { return Xsuper_.empty()?0:Xsuper_.size()-1;}

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
      DistSparseMatrixGraph graph_;
      SparseMatrixStructure * Local_;
      SparseMatrixStructure * Global_;
      //Is the global structure of the matrix allocated
      bool isGlobStructAllocated_;

      //CSC structure of L factor
      PtrVec xlindx_;
      IdxVec lindx_;

      PtrVec locXlindx_;
      IdxVec locLindx_;

#ifdef MULTITHREADING
      SYMPACK::vector<omp_lock_t> superLocks_;
#endif

      BackupBuffer backupBuffer_;

      //Supernodal partition array: supernode I ranges from column Xsuper_[I-1] to Xsuper_[I]-1
      SYMPACK::vector<Int> Xsuper_;
      SYMPACK::vector<Int> XsuperDist_;
      //Supernode membership array: column i belongs to supernode SupMembership_[i-1]
      SYMPACK::vector<Int> SupMembership_;

      //TODO Task lists
      Scheduler<std::list<FBTask>::iterator> * scheduler_;
      Scheduler<FBTask> * scheduler2_;
      //SYMPACK::vector<std::list<FBTask> * > taskLists_;

      //backup for factorization
      SYMPACK::vector<Int> UpdatesToDo_;
      SYMPACK::vector<std::list<FBTask> * > origTaskLists_;
      Int localTaskCount_;
      void generateTaskGraph(Int & localTaskCount, SYMPACK::vector<std::list<FBTask> * > & taskLists);
      std::list<FBTask>::iterator find_task(SYMPACK::vector<std::list<FBTask> * > & taskLists, Int src, Int tgt, TaskType type );

      supernodalTaskGraph taskGraph_;
      void generateTaskGraph(supernodalTaskGraph & taskGraph,SYMPACK::vector<Int> & AggregatesToRecv,  SYMPACK::vector<Int>& LocalAggregates);

      SYMPACK::vector<std::list<Int> > chSupTree_;
      void dfs_traversal(SYMPACK::vector<std::list<Int> > & tree,int node,std::list<Int> & frontier);

template <class Allocator = UpcxxAllocator>
      SuperNode<T,Allocator> * CreateSuperNode(DecompositionType type,Int aiId, Int aiFc, Int aiLc, Int ai_num_rows, Int aiN, Int aiNZBlkCnt=-1);
template <class Allocator = UpcxxAllocator>
      SuperNode<T,Allocator> * CreateSuperNode(DecompositionType type,char * dataPtr,size_t size, Idx firstRow = -1);
template <class Allocator = UpcxxAllocator>
      SuperNode<T,Allocator> * CreateSuperNode(DecompositionType type,Int aiId, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices);

      std::list<std::list<FBTask>::iterator > readyTasks_;

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
      SYMPACK::vector<SuperNode<T> * > LocalSupernodes_;



      //Vector holding pointers to local contributions
      //This has to be renamed because it contains the distributed solution
      SYMPACK::vector<SuperNode<T,MallocAllocator> *> Contributions_;


#ifndef ITREE2
      SYMPACK::vector<Int> globToLocSnodes_;
#else
      ITree globToLocSnodes_;
#endif


      TempUpdateBuffers<T> tmpBufs;

  void findSupernodes(ETree& tree, Ordering & aOrder, SYMPACK::vector<Int> & cc,SYMPACK::vector<Int> & supMembership, SYMPACK::vector<Int> & xsuper, Int maxSize = -1);
  void relaxSupernodes(ETree& tree, SYMPACK::vector<Int> & cc,SYMPACK::vector<Int> & supMembership, SYMPACK::vector<Int> & xsuper, RelaxationParameters & params  );

  void symbolicFactorizationRelaxedDist(SYMPACK::vector<Int> & cc);

  void refineSupernodes(int ordflag,int altflag);

  void getLColRowCount(SYMPACK::vector<Int> & cc, SYMPACK::vector<Int> & rc);
  void getLColRowCount(SparseMatrixGraph & sgraph, SYMPACK::vector<Int> & cc, SYMPACK::vector<Int> & rc);


      protected:


      /******************* Global to Local Indexes utility routines ******************/
      //returns the 1-based index of supernode id global in the local supernode array
      Int snodeLocalIndex(Int global);
      //returns a reference to  a local supernode with id global
      SuperNode<T> * snodeLocal(Int global);
      template< class Alloc>
      SuperNode<T,Alloc> * snodeLocal(Int global, SYMPACK::vector<SuperNode<T,Alloc> *> & snodeColl);



      //FanBoth related routines
      Int FBUpdate(Int I,Int prevJ=-1);
      void FBGetUpdateCount(SYMPACK::vector<Int> & UpdatesToDo, SYMPACK::vector<Int> & AggregatesToRecv, SYMPACK::vector<Int> & LocalAggregates);
      void GetUpdatingSupernodeCount( SYMPACK::vector<Int> & sc,SYMPACK::vector<Int> & mw, SYMPACK::vector<Int> & mh, SYMPACK::vector<Int> & numBlk);

      void FBFactorizationTask(supernodalTaskGraph & taskGraph, FBTask & curTask, Int iLocalI, bool is_static = false);
      void FBAggregationTask(supernodalTaskGraph & taskGraph, FBTask & curTask, Int iLocalI, bool is_static = false);
      void FBUpdateTask(supernodalTaskGraph & taskGraph, FBTask & curTask, SYMPACK::vector<Int> & UpdatesToDo, SYMPACK::vector< SuperNode<T> * > & aggVectors, bool is_static = false);


      //Communication related routines
      void AddOutgoingComm(AsyncComms & outgoingSend, Icomm * send_buffer);
      void AdvanceOutgoing(AsyncComms & outgoingSend);
      template< class Alloc>
      void SendDelayedMessagesUp(Int cur_snode_id, CommList & MsgToSend, AsyncComms & OutgoingSend, SYMPACK::vector<SuperNode<T,Alloc> *> & snodeColl);
      template< class Alloc>
      void SendDelayedMessagesDown(Int iLocalI, DownCommList & MsgToSend, AsyncComms & OutgoingSend, SYMPACK::vector<SuperNode<T,Alloc> *> & snodeColl);
      template< class Alloc>
      void SendMessage(const DelayedComm & comm, AsyncComms & OutgoingSend, SYMPACK::vector<SuperNode<T,Alloc> *> & snodeColl);

      void CheckIncomingMessages(supernodalTaskGraph & taskGraph, bool is_static = false);


  };

} // namespace SYMPACK

#include "sympack/SupernodalMatrix_impl.hpp"


#ifdef NO_INTRA_PROFILE
#if defined (PROFILE)
#define TIMER_START(a) TAU_FSTART(a);
#define TIMER_STOP(a) TAU_FSTOP(a);
#define scope_timer(b,a) CTF_scope_timer b(#a)
#endif
#endif



#endif 
