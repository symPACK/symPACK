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
#ifndef _SYMPACK_MATRIX_DECL_HPP_
#define _SYMPACK_MATRIX_DECL_HPP_

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
#include "symPACKMatrixBase.hpp"

#include "symPACKMatrix2D.hpp"

#define FANIN_OPTIMIZATION
#define _SEQ_SPECIAL_CASE_
//#define _LAMBDAS_


//#define _SEQ_SPECIAL_CASE_
//#define PREFETCH_STRUCTURE

namespace symPACK{

  //Forward declarations
  template<typename Task> class supernodalTaskGraph;


  template <typename colptr_t, typename rowind_t, typename T> 
    class symPACKMatrix2D;




  template <typename T> class symPACKMatrix: public symPACKMatrixMeta<T>{
    public:

      //Constructors
      symPACKMatrix();
      symPACKMatrix(DistSparseMatrix<T> & pMat, symPACKOptions & options );
      //TODO
      symPACKMatrix( symPACKMatrix & M){};




      symPACKMatrix( symPACKMatrix2D<Ptr,Idx,T> & M){

    Mapping_ = nullptr;
    Balancer_ = nullptr;
    CommEnv_ = nullptr;
    scheduler_ = nullptr;
    scheduler2_ = nullptr;

    if(logfileptr==NULL){
      logfileptr = new LogFile(upcxx::rank_me(),false);
      logfileptr->OFS()<<"********* LOGFILE OF P"<<upcxx::rank_me()<<" *********"<<std::endl;
      logfileptr->OFS()<<"**********************************"<<std::endl;
    }
    logfileptr->OFS()<<"Shared node size "<<shmNode_.shmsize<<", rank "<<shmNode_.shmrank<<std::endl;


this->fullcomm_ = MPI_COMM_NULL;
        this->iSize_ = M.iSize_; 
        this->Init(M.options_);        

        this->graph_ = M.graph_;
        this->Order_ = M.Order_;
        //this->Order_.NpOrdering = this->options_.NpOrdering;
        //this->Order_.invp = M.Order_.invp;
        //this->Order_.perm = M.Order_.perm;
        this->ETree_ = M.ETree_;
        this->Xsuper_ = M.Xsuper_;
        this->SupMembership_ = M.SupMembership_;
        this->iam = M.iam;
        this->np = M.np;

        delete this->CommEnv_;
        this->CommEnv_ = new CommEnvironment(this->workcomm_);
        this->group_.reset( new RankGroup( this->workcomm_ ) ); 
        this->XsuperDist_ = M.XsuperDist_;
        this->locXlindx_ = M.locXlindx_;
        this->locLindx_ = M.locLindx_;

        //no balancer, balancing is done in 2D matrix
        this->Balancer_ = NULL;


      Int pmapping = this->np;
      Int pr = (Int)sqrt((double)pmapping);
      if(this->options_.mappingTypeStr ==  "MODWRAP2DTREE"){
        this->Mapping_ = new ModWrap2DTreeMapping(pmapping, pr, pr, 1);
      }
      else if(this->options_.mappingTypeStr ==  "WRAP2D"){
        this->Mapping_ = new Wrap2D(pmapping, pr, pr, 1);
      }
      else if(this->options_.mappingTypeStr ==  "WRAP2DFORCED"){
        this->Mapping_ = new Wrap2DForced(pmapping, pr, pr, 1);
      }
      else if(this->options_.mappingTypeStr ==  "MODWRAP2DNS"){
        this->Mapping_ = new Modwrap2DNS(pmapping, pr, pr, 1);
      }
      else if(this->options_.mappingTypeStr ==  "ROW2D"){
        this->Mapping_ = new Row2D(pmapping, pmapping, pmapping, 1);
      }
      else if(this->options_.mappingTypeStr ==  "COL2D"){
        this->Mapping_ = new Col2D(pmapping, pmapping, pmapping, 1);
      }
      else{
        this->Mapping_ = new Modwrap2D(pmapping, pr, pr, 1);
      }



        this->scheduler2_ = NULL;


        //Matrix is distributed in 2D fashion, need to pick one of the process holding blocks of a supernode and gather all blocks on it
      std::vector<int> load(M.np,0);
      std::vector<int> sup_mapp(M.nsuper,-1);

//      for ( auto & it : M.cells_ ) {
//        auto & cell = it.second;
//        if ( cell == nullptr ) continue;
//        if ( cell->owner==M.iam) {
//          auto blockptr = (typename symPACKMatrix2D<Ptr,Idx,T>::snodeBlock_t *)cell.get();
//          logfileptr->OFS()<<"orig "<<blockptr->nnz()<<std::endl;
//        }
//      }

      for (int supid=1;supid<=M.nsuper;supid++) {
        auto cell = M.pQueryCELL(supid-1,supid-1);
        if ( cell == nullptr ) continue;

        assert( supid == cell->j );
        int fc = M.Xsuper_[supid-1];
        auto & proot = sup_mapp[supid-1];
        if ( proot == -1 ) {
          logfileptr->OFS()<<"Working on supernode "<<cell->j<<std::endl;
          //find the least loaded processor among those owning blocks of supernode supid
          int cell_count = 0;
          std::set<int> proc_group;
          std::vector< std::shared_ptr<typename symPACKMatrix2D<Ptr,Idx,T>::snodeBlock_t> > localBlocks;

          int minLoadP = -1;
          int minLoad = -1;
          for (int I=supid;I<=M.nsuper;I++) {
            auto cellptr = M.pQueryCELL(I-1,supid-1);
            if ( cellptr != nullptr ) {
              cell_count++;
              proc_group.insert(cellptr->owner);
              if ( cellptr->owner==M.iam) {
                auto blockptr = M.pCELL(I-1,supid-1);
                localBlocks.push_back(blockptr);
//                logfileptr->OFS()<<"orig "<<blockptr->nnz()<<std::endl;
              }

              if ( minLoadP == -1 || minLoad > load[cellptr->owner]) {
                minLoadP = cellptr->owner;
                minLoad = load[cellptr->owner];
              }
            }
          }
          proot = minLoadP;
          load[proot] += 1;
//          logfileptr->OFS()<<"Proot is "<<proot<<std::endl;

          //Gather all cells on proot
          //create a dist_object only to access the vector from within the RPC
          if ( M.iam != proot )
            cell_count = 0;
          

          using cell_desc_t = std::tuple</*std::unique_ptr<typename symPACKMatrix2D<Ptr,Idx,T>::snodeBlock_t>,*/ std::vector<char>, size_t, int,int, size_t, size_t, Idx, int >;
          
          struct cell_list_t {
            std::vector< cell_desc_t > cells;
//            std::set< cell_desc_t, cell_desc_comp > cells;
            upcxx::promise<> prom;
            cell_list_t(int count){
              this->prom.require_anonymous(count);
              cells.reserve(count);
            }
          };

          struct cell_desc_comp {
              bool operator () (const cell_desc_t & a, const cell_desc_t & b){return std::get<2>(a) < std::get<2>(b);}
          } cell_desc_comp;

          auto cell_list = new cell_list_t(cell_count);

          upcxx::dist_object< cell_list_t * > remote_ptrs( cell_list );
upcxx::barrier();
          upcxx::promise<> prom;//(localBlocks.size());
          for ( auto & cellptr : localBlocks ) {
            assert( cellptr->j == supid );
          }
  
          for ( auto & cellptr : localBlocks ) {
            upcxx::rpc_ff(proot, upcxx::source_cx::as_promise(prom), 
                            [&remote_ptrs,supid,fc,proot]( upcxx::global_ptr< char > ptr, std::size_t size, size_t nnz, size_t nblocks, Idx width, int sender, int i ) {
                              //allocate a landing_zone
                              (*remote_ptrs)->cells.push_back( std::make_tuple( std::vector<char>(size), size,i, supid, nnz, nblocks, width, sender) );
                              //auto res = (*remote_ptrs)->cells.insert( std::move(std::make_tuple( std::vector<char>(size), size,i, supid, nnz, nblocks, width, sender)) );
                              cell_desc_t * descptr;
                              descptr = &(*remote_ptrs)->cells.back();
                              //descptr = (typename cell_list_t::cell_desc_t *)&*res.first;
                              
                              //issue a rget
                              upcxx::rget(ptr,std::get<0>(*descptr).data(),size).then([&remote_ptrs]() {(*remote_ptrs)->prom.fulfill_anonymous(1); });
//upcxx::operation_cx::as_lpc( upcxx::master_persona(), [&remote_ptrs,descptr,fc,width,nnz,nblocks,supid,size,i]() {
//                                auto & vect = std::get<1>(*descptr);
//                assert (vect.size() == size);
//                //logfileptr->OFS()<<(*remote_ptrs)<<" "<<descptr<<" Executing lpc post rget "<<vect.size()<<" vs "<<size<<std::endl;
//                                auto lptr = std::get<1>(*descptr).data();
//                                //std::get<0>(*descptr) = std::unique_ptr<typename symPACKMatrix2D<Ptr,Idx,T>::snodeBlock_t>( new typename symPACKMatrix2D<Ptr,Idx,T>::snodeBlock_t(i,supid,lptr,fc,width,nnz,nblocks)); 
//                                (*remote_ptrs)->prom.fulfill_anonymous(1); }));
                            }, cellptr->_gstorage, cellptr->_storage_size,cellptr->nnz(), cellptr->nblocks(), std::get<0>(cellptr->_dims), upcxx::rank_me() , cellptr->i );
          }
       
          //make sure everything has been sent 
          prom.finalize().wait(); 

          //wait until we get all the data 
          cell_list->prom.finalize().wait();
          if ( M.iam == proot ) {
            assert ( cell_list->cells.size() == cell_count );
          } 
        
          //everything has been created from this point 
          if ( M.iam == proot ) {
            //do the conversion
            Int I = supid;
            int fc = this->Xsuper_[I-1];
            int lc = this->Xsuper_[I]-1;

#ifndef ITREE2
            globToLocSnodes_.push_back(I-1);
#else 
            ITree::Interval snode_inter = { I, I, LocalSupernodes_.size() };
            globToLocSnodes_.Insert(snode_inter);
#endif
            SuperNode<T> * newSnode = NULL;

                        try{
            size_t height = 0;
            size_t nzBlockCnt = 0;
            
            std::sort(cell_list->cells.begin(), cell_list->cells.end(), cell_desc_comp );

            for ( auto & cell: cell_list->cells ) { 
              std::vector<char> storage;
              size_t size;
              int i,j;
              size_t nnz;
              size_t nblocks;
              Idx width;
              int sender;
              std::tie( storage, size, i,j, nnz, nblocks, width, sender) = cell;
              assert(I==j);
              height += nnz/width;
              nzBlockCnt += nblocks;
            }
            newSnode = CreateSuperNode(this->options_.decomposition,I,fc,fc,lc,height,this->iSize_,nzBlockCnt,this->options_.panel);
            auto nzval = newSnode->GetNZval(0);

            //now add blocks
            for ( auto & cell: cell_list->cells ) { 
              std::vector<char> storage;
              size_t size;
              int i,j;
              size_t nnz;
              size_t nblocks;
              Idx width;
              int sender;
              std::tie( storage, size, i,j, nnz, nblocks, width, sender) = cell;
              std::unique_ptr<typename symPACKMatrix2D<Ptr,Idx,T>::snodeBlock_t> snode2D( new typename symPACKMatrix2D<Ptr,Idx,T>::snodeBlock_t(i,j,storage.data(),fc,width,nnz,nblocks));
//              logfileptr->OFS()<<snode2D->nnz()<<std::endl;
              auto & blocks = snode2D->blocks();
              for ( auto & block: blocks ) {
                newSnode->AddNZBlock(snode2D->block_nrows(block) , block.first_row);
              }
              //now copy numerical values
              std::copy(snode2D->_nzval,snode2D->_nzval+nnz,nzval);
              nzval+=nnz;
            }

                        }
                        catch(const MemoryAllocationException & e){
                          std::stringstream sstr;
                          sstr<<"There is not enough memory on the system to allocate the factors. Try to increase GASNET_MAX_SEGSIZE value."<<std::endl;
                          logfileptr->OFS()<<sstr.str();
                          std::cerr<<sstr.str();
                          throw(e);
                        }

            LocalSupernodes_.push_back(newSnode);
          }

//          upcxx::barrier(); 
          delete cell_list;
        }

 
      }
      upcxx::barrier(); 

//      logfileptr->OFS()<<sup_mapp<<std::endl;
//      for( auto p: sup_mapp){ assert(p!=-1); }

        //create a Mapping
        //this->Mapping_ has to be updated to reflect the distribution 
            this->Mapping_->Update(sup_mapp);

        
            //starting from now, this only concerns the working processors
    SYMPACK_TIMER_START(Get_UpdateCount);
    GetUpdatingSupernodeCount(UpdateCount_,UpdateWidth_,UpdateHeight_,numBlk_);
    SYMPACK_TIMER_STOP(Get_UpdateCount);

    return;
    if(this->iam<this->np){

      {
        double timeSta = get_time();
        std::vector<Int> AggregatesToRecv;
        std::vector<Int> LocalAggregates;


        FBGetUpdateCount(UpdatesToDo_,AggregatesToRecv,LocalAggregates);
        generateTaskGraph(taskGraph_, AggregatesToRecv, LocalAggregates);
          {
#ifdef _MEM_PROFILER_
          utility::scope_memprofiler m("symPACKMatrix_task_graph");
#endif
        generateTaskGraph_New(taskGraph_New_, AggregatesToRecv, LocalAggregates,UpdateWidth_,UpdateHeight_);
          }


        //#define _OUTPUT_TASK_GRAPH_
#ifdef _OUTPUT_TASK_GRAPH_
        logfileptr->OFS()<<"tasks: ";
        for(auto binit = taskGraph_.taskLists_.begin(); binit != taskGraph_.taskLists_.end(); binit++){
          if((*binit)!=NULL){
            for(auto taskit = (*binit)->begin(); taskit!= (*binit)->end(); taskit++){
              logfileptr->OFS()<<"t_"<<taskit->src_snode_id<<"_"<<taskit->tgt_snode_id<<" ";
            }
          }

        }
        logfileptr->OFS()<<std::endl; 
        logfileptr->OFS()<<"taskmap: ";
        for(auto binit = taskGraph_.taskLists_.begin(); binit != taskGraph_.taskLists_.end(); binit++){
          if((*binit)!=NULL){
            for(auto taskit = (*binit)->begin(); taskit!= (*binit)->end(); taskit++){
              logfileptr->OFS()<<this->iam<<" ";
            }
          }
        }
        logfileptr->OFS()<<std::endl; 
#endif

        double timeStop = get_time();
        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Task graph generation time: "<<timeStop - timeSta<<std::endl;
        }
      }
      MPI_Barrier(CommEnv_->MPI_GetComm());
    }


                logfileptr->OFS()<<"DONE"<<std::endl;

      };



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
      Int LocalSupernodeCnt(){ return LocalSupernodes_.size(); } 
      const Mapping * GetMapping(){return Mapping_;}
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
      void FanBoth_New( );


      //debug routines
      void DumpMatlab();
      void Dump();
      void DumpContrib();









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
        SuperNode<T,Allocator> * CreateSuperNode(DecompositionType type,Int aiId, Int aiFr, Int aiFc, Int aiLc, Int ai_num_rows, Int aiN, Int aiNZBlkCnt=-1, Int panel=-1);


      template <class Allocator = UpcxxAllocator>
        SuperNode<T,Allocator> * CreateSuperNode(DecompositionType type,char * dataPtr,size_t size, Int firstRow = -1);

      template <class Allocator = UpcxxAllocator>
        SuperNode<T,Allocator> * CreateSuperNode(DecompositionType type,Int aiId, Int aiFr, Int aiFc, Int aiLc, Int aiN, std::set<Idx> & rowIndices, Int panel=-1);


    protected:
      CommEnvironment * CommEnv_;
      //upcxx::team * team_;

      SharedNode shmNode_;

      //MAPCLASS describing the Mapping of the computations
      Mapping * Mapping_;
      LoadBalancer * Balancer_;



      std::vector<Int> numBlk_;
      std::vector<Int> cc_,rc_;


      BackupBuffer backupBuffer_;

#ifdef NEW_UPCXX
      std::list< upcxx::future<> > gFutures;
#endif

      //TODO Task lists
      Scheduler<std::list<FBTask>::iterator> * scheduler_;
      Scheduler<FBTask> * scheduler2_;
      
      std::shared_ptr< Scheduler< std::shared_ptr<GenericTask> > > scheduler_new_;

      //backup for factorization
      std::vector<Int> UpdatesToDo_;
      Int localTaskCount_;

      supernodalTaskGraph<FBTask> taskGraph_;
      void generateTaskGraph(supernodalTaskGraph<FBTask> & taskGraph,std::vector<Int> & AggregatesToRecv,  std::vector<Int>& LocalAggregates);

      taskGraph taskGraph_New_;
      void generateTaskGraph_New(taskGraph & taskGraph, std::vector<Int> & AggregatesToRecv,  std::vector<Int>& LocalAggregates, std::vector<Int> & mw, std::vector<Int> & mh);

      std::vector<std::list<Int> > chSupTree_;
      void dfs_traversal(std::vector<std::list<Int> > & tree,int node,std::list<Int> & frontier);
      void solve_(T * RHS, int nrhs,  T * Xptr=NULL);
      void solveNew_(T * RHS, int nrhs,  T * Xptr=NULL);
      void solveNew2_(T * RHS, int nrhs,  T * Xptr=NULL);

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
#ifdef NEW_UPCXX
      std::vector< std::tuple< upcxx::global_ptr<char>,Int> > remoteFactors_;
#else
      std::vector< std::tuple< upcxx::global_ptr<SuperNodeDesc >,Int> > remoteFactors_;
#endif


      //Vector holding pointers to local contributions
      //This has to be renamed because it contains the distributed solution
      std::vector<SuperNode<T,MallocAllocator> *> Contributions_;
      std::vector< std::shared_ptr<SuperNodeBase<T> > > Contributions2_;


#ifndef ITREE2
      std::vector<Int> globToLocSnodes_;
#else
      ITree globToLocSnodes_;
#endif


      TempUpdateBuffers<T> tmpBufs;
#ifdef SP_THREADS
      std::map<std::thread::id,TempUpdateBuffers<T> > tmpBufs_th;
#endif


    protected:


      //FanBoth related routines
      Int FBUpdate(Int I,Int prevJ=-1);
      void FBGetUpdateCount(std::vector<Int> & UpdatesToDo, std::vector<Int> & AggregatesToRecv, std::vector<Int> & LocalAggregates);
      void GetUpdatingSupernodeCount( std::vector<Int> & sc,std::vector<Int> & mw, std::vector<Int> & mh, std::vector<Int> & numBlk);

      void FBFactorizationTask(supernodalTaskGraph<FBTask> & taskGraph, FBTask & curTask, Int iLocalI, std::vector< SuperNode<T> * > & aggVectors, bool is_static = false);
      void FBAggregationTask(supernodalTaskGraph<FBTask> & taskGraph, FBTask & curTask, Int iLocalI, bool is_static = false);
      void FBUpdateTask(supernodalTaskGraph<FBTask> & taskGraph, FBTask & curTask, std::vector<Int> & UpdatesToDo, std::vector< SuperNode<T> * > & aggVectors, bool is_static = false);


      //Communication related routines
      void AddOutgoingComm(AsyncComms & outgoingSend, Icomm * send_buffer);
      void AdvanceOutgoing(AsyncComms & outgoingSend);
      template< class Alloc>
        void SendDelayedMessagesUp(Int cur_snode_id, CommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode<T,Alloc> *> & snodeColl);
      template< class Alloc>
        void SendDelayedMessagesDown(Int iLocalI, DownCommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode<T,Alloc> *> & snodeColl);
      template< class Alloc>
        void SendMessage(const DelayedComm & comm, AsyncComms & OutgoingSend, std::vector<SuperNode<T,Alloc> *> & snodeColl);

      void CheckIncomingMessages(supernodalTaskGraph<FBTask> & taskGraph, std::vector< SuperNode<T> * > & aggVectors, bool is_static = false);



      template<typename Task>
      void CheckIncomingMessages_Solve(supernodalTaskGraph<Task> & taskGraph, std::shared_ptr<Scheduler<Task> > scheduler);





      //functions replacing lambdas
      //inline void _factorTask1D(Int src, Int tgt, const std::shared_ptr<GenericTask> & pTask, const Factorization::op_type & type, taskGraph & graph);
      //inline void _updateTask1D(Int src, Int tgt, const std::shared_ptr<GenericTask> & pTask, const Factorization::op_type & type, std::vector<Int> & UpdatesToDo,std::vector< SuperNode<T>* >& aggVectors, taskGraph & graph);
      //inline void _dec_ref_task(taskGraph & graph, taskGraph::task_iterator & taskit, Int loc, Int rem);






  };

typedef symPACKMatrix<float>       symPACKMatrixFloat;
typedef symPACKMatrix<double>       symPACKMatrixDouble;
typedef symPACKMatrix<std::complex<float> >       symPACKMatrixComplex;
typedef symPACKMatrix<std::complex<double> >       symPACKMatrixDoubleComplex;

} // namespace SYMPACK

#include <sympack/impl/symPACKMatrix_impl.hpp>

//namespace symPACK{
//struct symPACK_handle;
////class symPACKMatrixBase;
////class DistSparseMatrixBase;
//
//
//  extern std::map<int, symPACK_handle  > symPACK_handles;
//
//
//
//}




#endif //_SYMPACK_MATRIX_DECL_HPP_

