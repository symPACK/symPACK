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
#include "sympack/mpi_interf.hpp"

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


//#define _SEQ_SPECIAL_CASE_
//#define PREFETCH_STRUCTURE

namespace symPACK{

  //Forward declarations
  template<typename Task> class supernodalTaskGraph;


  template <typename colptr_t, typename rowind_t, typename T, typename int_t> 
    class symPACKMatrix2D;




  template <typename T> class symPACKMatrix: public symPACKMatrixMeta<T>{
    public:

      //Constructors
      symPACKMatrix();
      symPACKMatrix(DistSparseMatrix<T> & pMat, symPACKOptions & options );
      //TODO
      symPACKMatrix( symPACKMatrix & M){};

      symPACKMatrix( symPACKMatrix2D<Ptr,Idx,T,int> & M){

        using cell_desc_t = std::tuple< std::vector<char>, size_t, int,int, size_t, size_t, Idx, int >;

        struct cell_list_t {
          std::vector< cell_desc_t > cells;
          upcxx::promise<> prom;
          cell_list_t(int count){
            this->prom.require_anonymous(count);
            cells.reserve(count);
          }
        };

        struct cell_desc_comp {
          bool operator () (const cell_desc_t & a, const cell_desc_t & b){return std::get<2>(a) < std::get<2>(b);}
        } cell_desc_comp;


        upcxx::dist_object< std::vector<cell_list_t *> > remote_ptrs( std::vector<cell_list_t *>(M.nsuper,nullptr) );


        Mapping_ = nullptr;
        Balancer_ = nullptr;
        CommEnv_ = nullptr;
        scheduler_ = nullptr;
        scheduler2_ = nullptr;

        if(logfileptr==nullptr){
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
        this->Balancer_ = nullptr;


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
        else if(this->options_.mappingTypeStr ==  "FANOUT"){
          this->Mapping_ = new Row2D(pmapping, pmapping, pmapping, 1);
        }
        else if(this->options_.mappingTypeStr ==  "FANIN"){
          this->Mapping_ = new Col2D(pmapping, pmapping, pmapping, 1);
        }
        else{
          this->Mapping_ = new Modwrap2D(pmapping, pr, pr, 1);
        }



        this->scheduler2_ = nullptr;

        double timeSta = get_time();
        //Matrix is distributed in 2D fashion, need to pick one of the process holding blocks of a supernode and gather all blocks on it
        std::vector<int> sup_mapp(M.nsuper,-1);

        std::vector<Int> cell_counts(this->Xsuper_.size(),0);
        //Communicate the number of blocks in each supernode.
        assert( this->locXlindx_.size()-1 == this->XsuperDist_[this->iam+1] - this->XsuperDist_[this->iam] );
        for (int locSupid = 1; locSupid <= this->locXlindx_.size()-1; locSupid++ ) {

          int supid = this->XsuperDist_[this->iam]+locSupid-1;
          auto colbeg = this->locXlindx_[locSupid-1];
          auto colend = this->locXlindx_[locSupid]-1;
          int prevSnode = -1;
          int cell_count = 0;
          for ( auto ptr = colbeg; ptr <= colend; ptr++ ) {
            auto row = this->locLindx_[ptr-1];
            int J = this->SupMembership_[row-1];
            if(J!=prevSnode){
              cell_count++;
            }
            prevSnode = J;
          }
          cell_counts[supid-1] = cell_count;
        }

        {
          MPI_Datatype type;
          MPI_Type_contiguous( sizeof(int), MPI_BYTE, &type );
          MPI_Type_commit(&type);

          int supidBeg = this->XsuperDist_[this->iam];
          int supidEnd = this->XsuperDist_[this->iam+1];
          //allgather receive sizes
          size_t ssize = supidEnd - supidBeg;
          vector<int> rsizes(this->np,0);
          MPI_Allgather(&ssize,sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,this->workcomm_);

          //compute receive displacements
          vector<int> rdispls(this->np+1,0);
          rdispls[0] = 0;
          std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);


          //Now do the allgatherv
          vector<int> cell_counts2(rdispls.back());
          MPI_Allgatherv(&cell_counts[supidBeg-1],ssize,type,cell_counts2.data(),rsizes.data(),rdispls.data(),type,this->workcomm_);
          MPI_Type_free(&type);
          cell_counts.swap(cell_counts2);
        }

        using cell_meta_t = std::tuple< std::size_t , size_t, size_t , Idx , int, int > ;
        std::vector< char > sendbuf;
        std::vector<size_t> ssizes(this->np,0);

        auto local_it = M.localBlocks_.begin();
        for (int supid=1;supid<=M.nsuper;supid++) {
          auto cell = M.pQueryCELL(supid-1,supid-1);
          assert( cell != nullptr );
          assert( supid == cell->j );
          int fc = M.Xsuper_[supid-1];
          auto proot = cell->owner;
          sup_mapp[supid-1] = proot;

          while ( local_it != M.localBlocks_.end() && (*local_it)->j<supid  ) { local_it++; }
          while ( local_it != M.localBlocks_.end() && (*local_it)->j == supid ) {
            auto cellptr = *local_it;
            ssizes[proot]+=sizeof(cell_meta_t)+cellptr->_storage_size;
            local_it++;
          }
        }

        vector<size_t> sdispls(this->np+1,0);
        sdispls[0] = 0;
        std::partial_sum(ssizes.begin(),ssizes.end(),&sdispls[1]);

        sendbuf.resize(sdispls.back());
        local_it = M.localBlocks_.begin();
        for (int supid=1;supid<=M.nsuper;supid++) {
          auto cell = M.pQueryCELL(supid-1,supid-1);
          assert( cell != nullptr );
          assert( supid == cell->j );
          int fc = M.Xsuper_[supid-1];
          auto proot = cell->owner;
          sup_mapp[supid-1] = proot;

          //find the least loaded processor among those owning blocks of supernode supid
          while ( local_it != M.localBlocks_.end() && (*local_it)->j<supid  ) { local_it++; }
          while ( local_it != M.localBlocks_.end() && (*local_it)->j == supid ) {
            auto cellptr = *local_it;
            cell_meta_t tuple = std::make_tuple(cellptr->_storage_size,cellptr->nnz(), cellptr->nblocks(), std::get<0>(cellptr->_dims), cellptr->i, cellptr->j);
            char * tail = &sendbuf[sdispls[proot]];
            *(cell_meta_t*)tail = tuple;
            tail+=sizeof(tuple);
            std::copy(cellptr->_gstorage.local(),cellptr->_gstorage.local()+cellptr->_storage_size,tail);
            tail+=cellptr->_storage_size;
            sdispls[proot]+=sizeof(tuple)+cellptr->_storage_size;
            local_it++;
          }
        }

        sdispls[0] = 0;
        std::partial_sum(ssizes.begin(),ssizes.end(),&sdispls[1]);

        
 
        vector<char> recvbuf;
        vector<size_t> rsizes(this->np,0);
        vector<size_t> rdispls(this->np+1,0);

        std::function<void(std::vector<char> &,size_t)> resize_lambda =
            [](std::vector<char> & container, size_t sz){
            container.resize(sz); 
            };

        mpi::Alltoallv(sendbuf, ssizes.data(),sdispls.data(), MPI_BYTE,
            recvbuf,rsizes.data(),rdispls.data(),this->fullcomm_,resize_lambda );
        
        //clear the send buffer
        { vector<char> tmp; sendbuf.swap( tmp ); }

        Int curSnode = 1;
      

        auto heads = rdispls; 
        while ( std::any_of(heads.begin(),heads.end()-1,[&heads,&rdispls](size_t & a){ size_t idx = &a-heads.data(); return a<rdispls[idx+1];}) ){
          bassert(curSnode<=M.nsuper);
          SuperNode<T> * newSnode = nullptr;
          std::vector< std::tuple<cell_meta_t,char*> > blocks;

          size_t height = 0;
          size_t nzBlockCnt = 0;


          for ( int psend = 0; psend < this->np; psend++ ) {
            Int localPrevSnode = -1;
            cell_meta_t tuple; 

            while ( heads[psend] < rdispls[psend+1] ) {
              tuple = *((cell_meta_t*)&recvbuf[heads[psend]]);
              bassert(std::get<5>(tuple)<=M.nsuper && std::get<5>(tuple)>0);
              if ( !((std::get<5>(tuple) == localPrevSnode  && localPrevSnode == curSnode) || localPrevSnode == -1) ) {
                break;
              }
              if ( std::get<5>(tuple) == curSnode ) {
                size_t size, nnz, nblocks;
                int i,j;
                Idx width;
                std::tie( size, nnz, nblocks, width, i,j) = tuple;         
                char * ptr = &recvbuf[heads[psend]+sizeof(tuple)];
                blocks.push_back(std::make_tuple(tuple,ptr));

                height += nnz/width;
                nzBlockCnt += nblocks;
                heads[psend] += sizeof(tuple) + size; 
              }
              localPrevSnode = std::get<5>(tuple);
            }
          }

          if ( blocks.size()>0 ){
#ifndef ITREE2
              this->globToLocSnodes_.push_back(curSnode-1);
#else 
              ITree::Interval snode_inter = { curSnode, curSnode, this->LocalSupernodes_.size() };
              this->globToLocSnodes_.Insert(snode_inter);
#endif

            int fc = this->Xsuper_[curSnode-1];
            int lc = this->Xsuper_[curSnode]-1;
            newSnode = CreateSuperNode(this->options_.decomposition,curSnode,fc,fc,lc,height,this->iSize_,nzBlockCnt,this->options_.panel);

            auto nzval = newSnode->GetNZval(0);


            std::sort(blocks.begin(),blocks.end(),[](const std::tuple<cell_meta_t,char*> & a, const std::tuple<cell_meta_t,char*> & b){ 
                auto ai = std::get<4>(std::get<0>(a));
                auto aj = std::get<5>(std::get<0>(a));
                auto bi = std::get<4>(std::get<0>(b));
                auto bj = std::get<5>(std::get<0>(b));
              return aj < bj || (ai < bi &&  aj == bj ) ;
              }); 


            for(auto & tpl: blocks){
              cell_meta_t & tuple = std::get<0>(tpl);
              size_t size, nnz, nblocks;
              int i,j;
              Idx width;
              std::tie( size, nnz, nblocks, width, i,j) = tuple;         
              char * ptr = std::get<1>(tpl);

              std::unique_ptr<typename symPACKMatrix2D<Ptr,Idx,T,int>::snodeBlock_t> snode2D;
              
              if(this->options_.decomposition == DecompositionType::LDL){
                snode2D.reset( new typename symPACKMatrix2D<Ptr,Idx,T,int>::snodeBlockLDL_t(i,j,ptr,fc,width,nnz,nblocks));
              }
              else {
                snode2D.reset( new typename symPACKMatrix2D<Ptr,Idx,T,int>::snodeBlock_t(i,j,ptr,fc,width,nnz,nblocks));
              }
              
              auto & blocks2D = snode2D->blocks();
              for ( auto & block: blocks2D ) {
                newSnode->AddNZBlock(snode2D->block_nrows(block) , block.first_row);
              }
              //now copy numerical values
              std::copy(snode2D->_nzval,snode2D->_nzval+nnz,nzval);
              nzval+=nnz;

              if(this->options_.decomposition == DecompositionType::LDL){
                if ( snode2D->i == snode2D->j ){
                  auto snode2DLDL = (typename symPACKMatrix2D<Ptr,Idx,T,int>::snodeBlockLDL_t*)snode2D.get();
                  //now copy diagonal
                  auto srcdiag = (T*)snode2DLDL->GetDiag();
                  auto diag = ((SuperNodeInd<T>*)newSnode)->GetDiag();
                  std::copy(srcdiag,srcdiag+width,diag);
                }
              }
            }

            this->LocalSupernodes_.push_back(newSnode);
          }
          curSnode++;
        }

        double timeEnd = get_time();
        if(this->iam==0){
          symPACKOS<<"Conversion to 1D time: "<<timeEnd-timeSta<<std::endl;
        }

        //create a Mapping
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
            {
              generateTaskGraph(taskGraph_, AggregatesToRecv, LocalAggregates,UpdateWidth_,UpdateHeight_);
            }
            double timeStop = get_time();
            if(this->iam==0 && this->options_.verbose){
              symPACKOS<<"Task graph generation time: "<<timeStop - timeSta<<std::endl;
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


      virtual void Init(symPACKOptions & options ) override;
      void Init(DistSparseMatrix<T> & pMat, symPACKOptions & options );

      virtual void SymbolicFactorization(DistSparseMatrix<T> & pMat) override;
      virtual void DistributeMatrix(DistSparseMatrix<T> & pMat) override;

      //Accessors
      Int LocalSupernodeCnt(){ return LocalSupernodes_.size(); } 
      const Mapping * GetMapping(){return Mapping_;}
      std::vector<SuperNode<T> *  > & GetLocalSupernodes(){ return LocalSupernodes_; }
      //TODO Check if that's useful
      SuperNode<T> & GetLocalSupernode(Int i){ return *LocalSupernodes_[i]; } 


      //core functionalities
      virtual void Factorize() override;
      //Solve routines
      //note: RHS & B are stored in column major format
      virtual void Solve(T * RHS, int nrhs, int dev_size,  T * Xptr=nullptr) override;
      virtual void GetSolution(T * B, int nrhs) override;

      void FanBoth( );


      //debug routines
      virtual void DumpMatlab() override;
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

      SharedNode shmNode_;

      //MAPCLASS describing the Mapping of the computations
      Mapping * Mapping_;
      LoadBalancer * Balancer_;



      std::vector<Int> numBlk_;
      std::vector<Int> cc_,rc_;


      BackupBuffer backupBuffer_;

      upcxx::dist_object<int> * remDealloc;
#ifdef SP_THREADS
  // declare an agreed upon persona for the progress thread
  upcxx::persona progress_persona_;
  std::thread progress_thread_;
#endif

      //TODO Task lists
      Scheduler<std::list<FBTask>::iterator> * scheduler_;
      Scheduler<FBTask> * scheduler2_;

      std::shared_ptr< Scheduler< std::shared_ptr<GenericTask> > > scheduler_new_;

      //backup for factorization
      std::vector<Int> UpdatesToDo_;
      Int localTaskCount_;


      taskGraph taskGraph_;
      void generateTaskGraph(taskGraph & taskGraph, std::vector<Int> & AggregatesToRecv,  std::vector<Int>& LocalAggregates, std::vector<Int> & mw, std::vector<Int> & mh);

      std::vector<std::list<Int> > chSupTree_;
      void dfs_traversal(std::vector<std::list<Int> > & tree,int node,std::list<Int> & frontier);
      void solve_(T * RHS, int nrhs,  T * Xptr=nullptr);

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
      std::vector< std::tuple< upcxx::global_ptr<char>,Int> > remoteFactors_;

      //Vector holding pointers to local contributions
      //This has to be renamed because it contains the distributed solution
      std::vector<SuperNode<T,MallocAllocator> *> Contributions_;
    public:
      std::vector< std::shared_ptr<SuperNodeBase<T> > > Contributions2_;
    protected:

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



      //Communication related routines
      void AddOutgoingComm(AsyncComms & outgoingSend, Icomm * send_buffer);
      void AdvanceOutgoing(AsyncComms & outgoingSend);
      template< class Alloc>
        void SendDelayedMessagesUp(Int cur_snode_id, CommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode<T,Alloc> *> & snodeColl);
      template< class Alloc>
        void SendDelayedMessagesDown(Int iLocalI, DownCommList & MsgToSend, AsyncComms & OutgoingSend, std::vector<SuperNode<T,Alloc> *> & snodeColl);
      template< class Alloc>
        void SendMessage(const DelayedComm & comm, AsyncComms & OutgoingSend, std::vector<SuperNode<T,Alloc> *> & snodeColl);
  };

typedef symPACKMatrix<float>       symPACKMatrixFloat;
typedef symPACKMatrix<double>       symPACKMatrixDouble;
typedef symPACKMatrix<std::complex<float> >       symPACKMatrixComplex;
typedef symPACKMatrix<std::complex<double> >       symPACKMatrixDoubleComplex;

} // namespace SYMPACK

#include <sympack/impl/symPACKMatrix_impl.hpp>

#endif //_SYMPACK_MATRIX_DECL_HPP_

