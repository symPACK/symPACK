#ifndef _TYPES_DECL_HPP_
#define _TYPES_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/timer.hpp"
#include "sympack/CommTypes.hpp"
#ifdef CUDA_MODE
#include "cuda_runtime.h"
#include "cublas_v2.h"
#include "cusolverDn.h"
#endif
#include <string>
#include <type_traits>
#include <cstdlib>

namespace symPACK{

  typedef vector<Ptr> PtrVec;
  typedef vector<Idx> IdxVec;


  struct RelaxationParameters{
    Int nrelax0;
    Int nrelax1;
    Int nrelax2;
    Int maxSize;

    double zrelax0 = 0.8;
    double zrelax1 = 0.1;
    double zrelax2 = 0.05;

    RelaxationParameters(){
      nrelax0 = 8 ;
      nrelax1 = 32;
      nrelax2 = 64;
    }
    RelaxationParameters(Int pmaxSize):RelaxationParameters(){
      maxSize = pmaxSize;
      if(maxSize>0){
        nrelax0 = std::min(nrelax0,maxSize);
        nrelax1 = std::min(nrelax1,maxSize);
        nrelax2 = std::min(nrelax2,maxSize);
      }

    }

    RelaxationParameters(Int pnrelax0, Int pnrelax1, Int pnrelax2, Int pmaxSize){
      maxSize = pmaxSize;
      if(maxSize>0){
        nrelax0 = std::min(pnrelax0,maxSize);
        nrelax1 = std::min(pnrelax1,maxSize);
        nrelax2 = std::min(pnrelax2,maxSize);
      }
      else{
        nrelax0 = pnrelax0;
        nrelax1 = pnrelax1;
        nrelax2 = pnrelax2;
      }

    }

    void SetMaxSize(Int pmaxSize){
      maxSize = pmaxSize;
      if(maxSize>0){
        nrelax0 = std::min(nrelax0,pmaxSize);
        nrelax1 = std::min(nrelax1,pmaxSize);
        nrelax2 = std::min(nrelax2,pmaxSize);
      }
    }
    void SetNrelax0(Int pnrelax0){
      if(maxSize>0){
        nrelax0 = std::min(pnrelax0,maxSize);
      }
      else{
        nrelax0 = pnrelax0;
      }
    }
    void SetNrelax1(Int pnrelax1){
      if(maxSize>0){
        nrelax1 = std::min(pnrelax1,maxSize);
      }
      else{
        nrelax1 = pnrelax1;
      }
    }
    void SetNrelax2(Int pnrelax2){
      if(maxSize>0){
        nrelax2 = std::min(pnrelax2,maxSize);
      }
      else{
        nrelax2 = pnrelax2;
      }
    }

  };

  enum class DecompositionType {LL,LDL};

  enum MappingType {ROW2D,COL2D,MODWRAP2D,MODWRAP2DNS,WRAP2D,WRAP2DFORCED};
  enum FactorizationType {FANOUT,FANBOTH,FANBOTH_STATIC};
  enum LoadBalanceType {NOLB,NNZ,NCOLS,WORK,SUBCUBE,SUBCUBE_NNZ};
  enum OrderingType {NATURAL,RCM,MMD,AMD,NDBOX,NDGRID,SCOTCH,PTSCOTCH,METIS,PARMETIS,USER};
  enum SchedulerType {DL,MCT,PR,FIFO};
  enum DataDistribution {SYMPACK_DATA_1D,SYMPACK_DATA_2D};
  //enum OrderRefinementType {BarryDL,MCT,PR,FIFO};
  class symPACKOptions{
    public:
      int panel;
      int numThreads;
      bool iterRefinement;
      int NpOrdering;
      bool print_stats;
#ifdef CUDA_MODE
      size_t gpu_alloc_size, gpu_block_limit, trsm_limit, potrf_limit, gemm_limit, syrk_limit;
      bool gpu_solve;
      bool gpu_verbose;
      FallbackType fallback_type;
#endif
      DecompositionType decomposition;
      MappingType mappingType;
      std::string mappingTypeStr;
      FactorizationType factorization;
      LoadBalanceType load_balance;
      std::string load_balance_str;
      std::string order_refinement_str;
      OrderingType ordering;
      std::string orderingStr;
      SchedulerType scheduler;
      Int maxIsend;
      Int maxIrecv;
      CommEnvironment * commEnv;
      RelaxationParameters relax;
      DataDistribution distribution;
      MPI_Comm MPIcomm;
      int verbose;
      int dumpPerm;
      int * perm;
      double memory_limit;
    protected:
      bool isSqrtP(){
        bool val = false;
        if(mappingTypeStr == "MODWRAP2D" || mappingTypeStr == "MODWRAP2DNS" || mappingTypeStr == "WRAP2D" || mappingTypeStr == "WRAP2DFORCED"){
          val = true;
        }

        //        switch(mappingType){
        //          case MODWRAP2D: case MODWRAP2DNS: case WRAP2D: case WRAP2DFORCED:
        //            val = true;
        //            break;
        //          default:
        //            val = false;
        //            break;
        //        }
        return val;
      }

    public:
      symPACKOptions(){
        panel = -1;
        numThreads = 1;
        verbose = 0;
        NpOrdering = 0;
        decomposition = DecompositionType::LL; 
        factorization = FANBOTH;
        mappingTypeStr = "ROW2D";
        load_balance_str = "SUBCUBE-FO";
        orderingStr = "MMD";
        order_refinement_str = "NO";
        scheduler = DL;
        maxIsend = 0;
        maxIrecv=0;
        distribution = SYMPACK_DATA_1D;
        relax = RelaxationParameters(150);
        relax.SetMaxSize(150);

        commEnv = nullptr;
        MPIcomm = MPI_COMM_NULL;

        perm = nullptr;
        dumpPerm = 0;
        //        mappingType = ROW2D;
        ///        ordering = MMD;
        //        load_balance = SUBCUBE;
        print_stats=false;
        iterRefinement=false;

        memory_limit = -1.0;
      }

      Int used_procs(Int np){

        if(isSqrtP()){
          Int nr = (Int)sqrt((double)np);
          Int nc = np/nr;
          //Int nc = nr;
          np= nr*nc;
        }
        return np;
      }

  };







  struct LocalUpdate{
    Int src_snode_id;
    Int src_nzblk_idx;
    Int src_first_row;
    LocalUpdate(Int snode_id,Int nzblk_idx,Int first_row):src_snode_id(snode_id),src_nzblk_idx(nzblk_idx),src_first_row(first_row){};
  };

  struct SnodeUpdate{
    Int src_snode_id;
    Int tgt_snode_id;
    Int src_first_row;
    Int src_next_row;
    Int blkidx;
    Int next_blkidx;
    //std::vector<bool> is_factor_sent;

    SnodeUpdate(){
      src_snode_id = 0;
      tgt_snode_id = 0;
      src_first_row = 0;
      src_next_row = 0;
      blkidx = 0;
      next_blkidx = 0;
      //is_factor_sent.assign(np,false);
    }
  };






  template <class T>
    class AlignedAllocator {
      public:
        // type definitions
        //typename std::aligned_storage<sizeof(T),alignof(T)>::type T_pod;
        using T_pod = typename std::aligned_storage<sizeof(T),alignof(T)>::type;

        typedef T        value_type;
        typedef T*       pointer;
        typedef const T* const_pointer;
        typedef T&       reference;
        typedef const T& const_reference;
        typedef std::size_t    size_type;
        typedef std::ptrdiff_t difference_type;

      protected:
        T_pod * storage_pod_;
      public:

        // rebind allocator to type U
        template <class U>
          struct rebind {
            typedef AlignedAllocator<U> other;
          };

        // return address of values
        pointer address (reference value) const {
          return &value;
        }
        const_pointer address (const_reference value) const {
          return &value;
        }

        /* constructors and destructor
         * - nothing to do because the allocator has no state
         */
        AlignedAllocator() throw() {
          storage_pod_ = nullptr;
        }

        AlignedAllocator(const AlignedAllocator&) throw() {
          storage_pod_ = nullptr;
        }

        template <class U>
          AlignedAllocator (const AlignedAllocator<U>&) throw() {
            storage_pod_ = nullptr;
          }

        ~AlignedAllocator() throw() {
        }

        // return maximum number of elements that can be allocated
        size_type max_size () const throw() {
          return std::numeric_limits<std::size_t>::max() / sizeof(T);
        }

        // allocate but don't initialize num elements of type T
        pointer allocate (size_type num, const void* = 0) {
          // print message and allocate memory with global new
          //std::cerr << "allocate " << num << " element(s)"
          //          << " of size " << sizeof(T) << std::endl;

          //pointer ret = (pointer)(::operator new(num*sizeof(T)));

          //storage_pod_ = reinterpret_cast<T_pod *>(::operator new( T_pod[num] ));
          storage_pod_ = new T_pod[num];
          pointer ret = reinterpret_cast<pointer>(storage_pod_);

          //std::cerr << " allocated at: " << (void*)ret << std::endl;

          return ret;
        }

        // initialize elements of allocated storage p with value value
        void construct (pointer p, const T& value) {
          // initialize memory with placement new
          //new((void*)p)T(value);
          new((void*)p)T(value);
        }

        // destroy elements of initialized storage p
        void destroy (pointer p) {
          // destroy objects by calling their destructor
          p->~T();
        }

        // deallocate storage p of deleted elements
        void deallocate (pointer p, size_type num) {
          // print message and deallocate memory with global delete
          //std::cerr << "deallocate " << num << " element(s)"
          //          << " of size " << sizeof(T)
          //          << " at: " << (void*)p << std::endl;

          //::operator delete((void*)p);
          ::operator delete[](this->storage_pod_);
          storage_pod_ = nullptr;
        }
    };

  // return that all specializations of this allocator are interchangeable
  template <class T1, class T2>
    bool operator== (const AlignedAllocator<T1>&,
        const AlignedAllocator<T2>&) throw() {
      return true;
    }
  template <class T1, class T2>
    bool operator!= (const AlignedAllocator<T1>&,
        const AlignedAllocator<T2>&) throw() {
      return false;
    }






  template<typename T>
    class TempUpdateBuffers{
      public:
        //std::vector<T  ,AlignedAllocator<T> > tmpBuf;
        //std::vector<Int,AlignedAllocator<Int> > src_colindx;
        //std::vector<Int,AlignedAllocator<Int> > src_to_tgt_offset;
        std::vector<T   > tmpBuf;
       
        std::vector<Int > src_colindx;
        std::vector<Int > src_to_tgt_offset;

        void Resize(Int size, Int mw){
          if(size*mw > tmpBuf.size()){
            tmpBuf.resize(size*mw);
          }
          if(mw > src_colindx.size()){
            src_colindx.resize(mw);
          }
          if(size > src_to_tgt_offset.size()){
            src_to_tgt_offset.resize(size);
          }
        }

        void Clear(){
          {
            std::vector<T> tmp;
            tmpBuf.swap(tmp);            
          }
          {
            std::vector<Int> tmp;
            src_colindx.swap(tmp);
          }
          {
            std::vector<Int> tmp;
            src_to_tgt_offset.swap(tmp);
          }
          //tmpBuf.clear();
          //src_colindx.clear();
          //src_to_tgt_offset.clear();
        }

        TempUpdateBuffers(){}
        TempUpdateBuffers(Int size, Int mw){
          Resize(size,mw);
        }
    };

  struct duet{
    Idx row;
    Idx col;
  };

  struct sortDuet {
    bool operator() (const duet & a, const duet & b){
      bool retval = a.row<b.row;
      if(a.row==b.row){
        retval = a.col<b.col;
      }
      return retval;
    }
  };

  struct sortDuetInv {
    bool operator() (const duet & a, const duet & b){
      bool retval = a.row>b.row;
      if(a.row==b.row){
        retval = a.col>b.col;
      }
      return retval;
    }
  };





  template<typename T>
    struct triplet{
      Idx row;
      Idx col;
      T val;
    };

  template<typename T>
    struct sortTriplet {
      bool operator() (const triplet<T> & a, const triplet<T> & b){
        bool retval = a.row<b.row;
        if(a.row==b.row){
          retval = a.col<b.col;
        }
        return retval;
      }
    };

  template<typename T>
    struct sortTripletInv {
      bool operator() (const triplet<T> & a, const triplet<T> & b){
        bool retval = a.row>b.row;
        if(a.row==b.row){
          retval = a.col>b.col;
        }
        return retval;
      }
    };




  struct OrderStats{
    int64_t totalSnodeBlocks = 0;
    int64_t totalBlocks = 0;
    double blocksPerSnode = 0;
    double blocksPerCol = 0;
    double avgSnodeBlockSize = 0;
    double avgBlockSize = 0;
    void reset(){
      totalSnodeBlocks = 0;
      totalBlocks = 0;
      blocksPerSnode = 0;
      blocksPerCol = 0;
      avgSnodeBlockSize = 0;
      avgBlockSize = 0;
    }
    void print(){
      symPACKOS<<"totalSnodeBlocks: "<<totalSnodeBlocks<<std::endl;
      symPACKOS<<"totalBlocks: "<<totalBlocks<<std::endl;
      symPACKOS<<"blocksPerSnode: "<<blocksPerSnode<<std::endl;
      symPACKOS<<"blocksPerCol: "<<blocksPerCol<<std::endl;
      symPACKOS<<"avgSnodeBlockSize (lines): "<<avgSnodeBlockSize<<std::endl;
      symPACKOS<<"avgBlockSize (nnz): "<<avgBlockSize<<std::endl;
    }


    void get(const std::vector<Int> & Xsuper, const std::vector<Int> & XsuperDist, const std::vector<Ptr> & locXlindx, const std::vector<Idx> & locLindx, const MPI_Comm & comm)
    {
      this->reset();

      int iam = 0;
      int mpisize = 1;
      MPI_Comm_size(comm,&mpisize);
      MPI_Comm_rank(comm,&iam);
      int64_t supNrows = 0;

      Int numLocSnode = XsuperDist[iam+1]-XsuperDist[iam];
      Int firstSnode = XsuperDist[iam];

      for(Int locsupno = 1; locsupno<locXlindx.size(); ++locsupno){
        Int I = locsupno + firstSnode-1;

        //count number of contiguous blocks (at least one diagonal block)
        Idx fc = Xsuper[I-1];
        Idx lc = Xsuper[I]-1;
        Ptr fi = locXlindx[locsupno-1];
        Ptr li = locXlindx[locsupno]-1;
        Idx iPrevRow = locLindx[fi-1]-1;
        Idx iFirstRow = locLindx[fi-1];

        Int width = lc - fc + 1; 

        for(Idx col = fc; col<=lc;col++){ 
          //1 to count the diagonal block, 0 to skip it
          int32_t nzBlockCnt = 0;//1;
          int32_t height = li - fi + 1;
          for(Ptr idx = fi; idx<=li;idx++){
            Idx iRow = locLindx[idx-1];
            if(iRow<col){
              --height; 
            }
            //enforce the first block to be a square diagonal block
            if(nzBlockCnt==1 && iRow>col){
              nzBlockCnt++;
              this->avgBlockSize+=col-iFirstRow+1;
              if(col==fc){
                this->avgSnodeBlockSize+=width;
              }
              iFirstRow=iRow;
            }
            else if(iRow!=iPrevRow+1){
              nzBlockCnt++;
              this->avgBlockSize+=iPrevRow-iFirstRow+1;
              if(col==fc){
                this->avgSnodeBlockSize+=iPrevRow-iFirstRow+1;
              }
              iFirstRow=iRow;
            }

            iPrevRow=iRow;
          }


          this->totalBlocks+=nzBlockCnt;
          if(col==lc){
            this->totalSnodeBlocks+=nzBlockCnt;
            supNrows+=height;
          }
        }
      }

      //Now reduce everything and then compute averages

      auto locavgSnodeBlockSize = this->avgSnodeBlockSize;
      auto locavgBlockSize =      this->avgBlockSize;
      auto loctotalBlocks =       this->totalBlocks; 
      auto loctotalSnodeBlocks =  this->totalSnodeBlocks;

      // logfileptr->OFS()<<"totalSnodeBlocks: "<<totalSnodeBlocks<<std::endl;
      // logfileptr->OFS()<<"totalBlocks: "<<totalBlocks<<std::endl;
      // logfileptr->OFS()<<"blocksPerSnode: "<<blocksPerSnode<<std::endl;
      // logfileptr->OFS()<<"blocksPerCol: "<<blocksPerCol<<std::endl;
      // logfileptr->OFS()<<"avgSnodeBlockSize: "<<avgSnodeBlockSize<<std::endl;
      // logfileptr->OFS()<<"avgBlockSize: "<<avgBlockSize<<std::endl;

      this->avgSnodeBlockSize = 0;
      this->avgBlockSize = 0;
      this->totalBlocks = 0; 
      this->totalSnodeBlocks = 0;


      MPI_Allreduce(&locavgSnodeBlockSize, &this->avgSnodeBlockSize, 1,MPI_DOUBLE,MPI_SUM,comm);
      MPI_Allreduce(&locavgBlockSize,      &this->avgBlockSize,      1,MPI_DOUBLE,MPI_SUM,comm);
      MPI_Allreduce(&loctotalBlocks,       &this->totalBlocks,       1,MPI_LONG_LONG,MPI_SUM,comm);
      MPI_Allreduce(&loctotalSnodeBlocks,  &this->totalSnodeBlocks,  1,MPI_LONG_LONG,MPI_SUM,comm);


      this->avgBlockSize/=this->totalBlocks;
      this->avgSnodeBlockSize/=this->totalSnodeBlocks;

      this->blocksPerCol=this->totalBlocks/(Xsuper.back()-1);
      this->blocksPerSnode=this->totalSnodeBlocks/(Xsuper.size()-1);
    }



  };


}


namespace symPACK{
  class MemoryAllocationException: public std::runtime_error {
    public:
      MemoryAllocationException(size_t psz)
        : std::runtime_error("Memory allocation error"), sz(psz)
      {
        //gdb_lock();
        std::stringstream err_sstr;
        err_sstr<<std::runtime_error::what()<< " of size "<<sz<<std::endl;
        err_str = err_sstr.str();
      }

      virtual char const * what() const throw() { 
        return err_str.c_str();
      }
    protected:
      size_t sz;
      std::string err_str;  
  };



  class MemoryAllocator{
    public:
      static char * allocate(size_t count){return nullptr;};

      static void deallocate(char* ptr) {};
  };


  class MallocAllocator: public MemoryAllocator{
    public:
      static char * allocate(size_t count){
        char * locTmpPtr = nullptr; 
        try{
          locTmpPtr = new char[count];
        }
        catch (const std::bad_alloc& e) {
          throw MemoryAllocationException(count);
        }

        return locTmpPtr;
      }

      static void deallocate(char* ptr){
        if(ptr!=nullptr){
          delete [] ptr;
        }
      }
  };





  class UpcxxAllocator: public MemoryAllocator{
    public:
      static char * allocate(size_t count){
        upcxx::global_ptr<char> tmpPtr;
        tmpPtr = upcxx::allocate<char>(count);
        char * locTmpPtr = (char*)tmpPtr.local();


        if(locTmpPtr == nullptr){
          throw MemoryAllocationException(count);
        }

        return locTmpPtr;
      }

      static void deallocate(char* ptr){
        if(ptr!=nullptr){

          upcxx::global_ptr<char> tmpPtr = upcxx::to_global_ptr((char*)ptr);
          upcxx::deallocate(tmpPtr);
        }
      }
  };



  template <class T>
    struct Mallocator {
      typedef T value_type;
      Mallocator() = default;
      template <class U> Mallocator(const Mallocator<U>&) {}
      T* allocate(std::size_t n) { return static_cast<T*>(std::malloc(n*sizeof(T))); }
      void deallocate(T* p, std::size_t) { std::free(p); }
    };
  template <class T, class U>
    bool operator==(const Mallocator<T>&, const Mallocator<U>&) { return true; }
  template <class T, class U>
    bool operator!=(const Mallocator<T>&, const Mallocator<U>&) { return false; }



}




#endif //_TYPES_DECL_HPP_

