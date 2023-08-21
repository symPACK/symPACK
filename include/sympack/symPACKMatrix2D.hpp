#ifndef _SYMPACK_MATRIX2D_DECL_HPP_
#define _SYMPACK_MATRIX2D_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/symPACKMatrixBase.hpp"
#include "sympack/symPACKMatrix.hpp"


#include "sympack/mpi_interf.hpp"

#ifdef CUDA_MODE
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cublas_v2.h"
#include "cusolverDn.h"
#endif

//#define _NO_COMPUTATION_
//#define LOCK_SRC 2

#define _SUBCUBE2D_
#define _EAGER_FETCH_
//#define _SUBCUBE_SQRT_LIMIT_
//#define _NO_SUBCUBE2D_

#define _POINTER_EXCHANGE_
#define _COMPACT_DEPENDENCIES_

#define _PRIORITY_QUEUE_AVAIL_
#define _PRIORITY_QUEUE_RDY_


#define _USE_PROM_AVAIL_
#define _USE_PROM_RDY_

#include <functional>
#include <chrono>


#ifdef _PRIORITY_QUEUE_RDY_
#define push_ready(sched,ptr) sched->ready_tasks.push(ptr);
#define top_ready() ready_tasks.top();
#define pop_ready() ready_tasks.pop();
#else
#define push_ready(sched,ptr) sched->ready_tasks.push_back(ptr);
#define top_ready() ready_tasks.front();
#define pop_ready() ready_tasks.pop_front();
#endif


#ifdef _PRIORITY_QUEUE_AVAIL_
#define push_avail(sched,ptr) sched->avail_tasks.push(ptr);
#define top_avail() avail_tasks.top();
#define pop_avail() avail_tasks.pop();
#else
#define push_avail(sched,ptr) sched->avail_tasks.push_back(ptr);
#define top_avail() avail_tasks.front();
#define pop_avail() avail_tasks.pop_front();
#endif


namespace symPACK{
  template <std::size_t alignment>
    inline bool is_aligned(void * ptr) noexcept {
      std::size_t max = 1u;
      return std::align(alignment, 1u, ptr, max)==ptr;
    }


  namespace scheduling {
    struct key_t { 
      unsigned int cell_J; 
      unsigned int cell_I;
      unsigned int src;
      Factorization::op_type type;
      int owner;
      key_t(unsigned int j, unsigned int i, unsigned int s, Factorization::op_type t, int o):
        cell_J(j),cell_I(i),src(s),type(t), owner(o) {}
      key_t():cell_J(0),cell_I(0),src(0),type(Factorization::op_type::FACTOR), owner(0) {}

      bool operator==( symPACK::scheduling::key_t& rhs)
      {
        return cell_J==rhs.cell_J && cell_I == rhs.cell_I && src == rhs.src && type == rhs.type && owner==rhs.owner;
      }

      bool operator<( symPACK::scheduling::key_t& rhs)
      {
        if (src != rhs.src) return src<rhs.src;
        else if (cell_I != rhs.cell_I) return cell_I < rhs.cell_I;
        else if (cell_J != rhs.cell_J) return cell_J < rhs.cell_J;
        else if (owner != rhs.owner) return owner < rhs.owner;
        return (int)type < (int)rhs.type;
      }

      bool operator==( const symPACK::scheduling::key_t& rhs) const
      {
        return cell_J==rhs.cell_J && cell_I == rhs.cell_I && src == rhs.src && type == rhs.type && owner == rhs.owner;
      }

      bool operator<( const symPACK::scheduling::key_t& rhs) const
      {
        if (src != rhs.src) return src<rhs.src;
        else if (cell_I != rhs.cell_I) return cell_I < rhs.cell_I;
        else if (cell_J != rhs.cell_J) return cell_J < rhs.cell_J;
        else if (owner != rhs.owner) return owner < rhs.owner;
        return (int)type < (int)rhs.type;
      }
    };
  }
}

// custom specialization of std::hash can be injected in namespace std
namespace std
{
  template<> struct hash<symPACK::scheduling::key_t>
  {
    typedef symPACK::scheduling::key_t argument_type;
    typedef std::size_t result_type;
    result_type operator()(argument_type const& s) const noexcept
    {
      //this is inspired by https://gist.github.com/nsuke/5990643
      result_type seed = 0;
      result_type const h1 ( std::hash<unsigned int>{}(s.cell_J) );
      seed ^= h1 + 0x9e3779b9 + (seed<<6) + (seed>>2);
      result_type const h2 ( std::hash<unsigned int>{}(s.cell_I) );
      seed ^= h2 + 0x9e3779b9 + (seed<<6) + (seed>>2);
      result_type const h3 ( std::hash<unsigned int>{}(s.src) );
      seed ^= h3 + 0x9e3779b9 + (seed<<6) + (seed>>2);
      result_type const h5 ( std::hash<int>{}((int)s.owner) );
      seed ^= h5 + 0x9e3779b9 + (seed<<6) + (seed>>2);
      result_type const h4 ( std::hash<int>{}((int)s.type) );
      seed ^= h4 + 0x9e3779b9 + (seed<<6) + (seed>>2);
      return seed; // or use boost::hash_combine (see Discussion)
    }
  };
}


namespace symPACK{


  class blockCellBase_t {
    using intrank_t = upcxx::intrank_t;
    public:
    intrank_t owner;
    int i;
    int j;
    blockCellBase_t():owner(0),i(0),j(0)
    {}
    virtual ~blockCellBase_t() {
    }

    virtual int I() {return i;}
  };


#ifdef SP_THREADS
  template< typename cellptr >
    class cell_lock {
      public:
        cell_lock(cellptr cell) {
          cell_ = cell;
        }

        ~cell_lock() {
          if (Multithreading::NumThread>2) cell_->in_use = false;
        }
      protected:
        cellptr cell_;
    };

  template<>
    class cell_lock< std::atomic<bool> > {
      public:
        cell_lock(std::atomic<bool> & celllock):lock_(celllock) {
        }

        ~cell_lock() {
          if (Multithreading::NumThread>2) lock_ =  false;
        }
      protected:
        std::atomic<bool> & lock_;
    };
#endif







  namespace scheduling {

    template <typename ttask_t, typename tmeta_t>
      class incoming_data_t;

    //(from I, to J,Factorization::op_type, updated first col, facing row));
    using meta_t = std::tuple<Idx,Idx,Factorization::op_type,Idx,Idx>;

    inline std::ostream& operator<<( std::ostream& os, const meta_t& s)
    {
      os<<"T"<<" "<<std::get<0>(s)<<" "<<std::get<1>(s)<<" ";

      switch(std::get<2>(s)) {
        case Factorization::op_type::TRSM:
          {
            os<<"TRSM";
          }
          break;
        case Factorization::op_type::FACTOR:
          {
            os<<"FACTOR";
          }
          break;
        case Factorization::op_type::UPDATE2D_COMP:
          {
            os<<"UPDATE"<< " "<<std::get<4>(s);
          }
          break;
        default:
	  os<<"UNSUPPORTED OP TYPE";
          break;
      }
      os<<std::endl;
      return os;
    }


    using depend_t = std::tuple<Int,Int>; //(I,J)

    template
      < typename tmeta_t = meta_t,
      typename tdepend_t = depend_t >
        class task_t{
          public:
            //cell description
            using depend_task_t = tdepend_t;
            using meta_t = tmeta_t;
            using data_t = incoming_data_t<task_t,tmeta_t>;
            using lock_t = std::atomic<bool>;
            lock_t * _lock_ptr;

            meta_t _meta;

            std::function< void() > execute;

            //byte storage for tasks 
            std::unordered_map<int,std::list<std::size_t> > out_dependencies;
            //need counter here as same REMOTE cell can be input to many tasks.
            std::vector< std::shared_ptr< data_t > > input_msg;

            int in_remote_dependencies_cnt;
            int in_local_dependencies_cnt;

#ifdef _USE_PROM_AVAIL_
            //promise to wait for all incoming_data_t to be created
            upcxx::promise<> in_avail_prom;
            int in_avail_counter;
#else
            std::atomic<int> in_avail_counter;
#endif

#ifdef _USE_PROM_RDY_
            //promise to sync on all incoming_data_t fetch
            upcxx::promise<> in_prom;
            int in_counter;
#else
            std::atomic<int> in_counter;
#endif

            template<typename T> 
              void satisfy_dep(int cnt,T&sched) {
#ifdef _USE_PROM_RDY_
                upcxx::master_persona().lpc_ff( 
                    [this,cnt] () {
                    bassert(this->in_counter-cnt >= 0);
                    this->in_counter-=cnt;
                    this->in_prom.fulfill_anonymous(cnt);
                    });
#else
                T * sched_ptr = &sched;
                upcxx::master_persona().lpc_ff( 
                    [this,sched_ptr,cnt] () {
                    this->in_counter-=cnt;
                    if ( this->in_counter == 0 ) {
#ifdef SP_THREADS
                    if (sched_ptr->extraTaskHandle_!=nullptr) {
                    bool delay = sched_ptr->extraTaskHandle_(this);
                    if (delay) {
                    sched_ptr->delayedTasks_[this->_lock_ptr].push_back(this);
                    return;
                    }
                    }
#endif
                    push_ready(sched_ptr,this)
                    }
                    });
#endif
              } 

            template<typename T> 
              void avail_dep(int cnt,T&sched) {
#ifdef _USE_PROM_AVAIL_
                upcxx::master_persona().lpc_ff( 
                    [this,cnt] () {
                    bassert(this->in_avail_counter-cnt >= 0);
                    this->in_avail_counter-=cnt;
                    this->in_avail_prom.fulfill_anonymous(cnt);
                    });
#else
                T * sched_ptr = &sched;
                upcxx::master_persona().lpc_ff( 
                    [this,sched_ptr,cnt] () {
                    bassert(this->in_avail_counter-cnt >= 0);
                    this->in_avail_counter-=cnt;
                    if ( this->in_avail_counter == 0 ) {
                    push_avail(sched_ptr,this)
                    }
                    });
#endif
              } 



            bool executed;

            task_t( ):executed(false),
            _lock_ptr(nullptr),
            in_remote_dependencies_cnt(0),
            in_local_dependencies_cnt(0) { }

            virtual ~task_t() = default;

            void reset() {
              input_msg.clear();
#ifdef _USE_PROM_RDY_
              in_prom = upcxx::promise<>();
#endif
#ifdef _USE_PROM_AVAIL_
              in_avail_prom = upcxx::promise<>();
#endif
              executed = false;
              _lock_ptr = nullptr;
            }

        };


    template <typename ttask_t = task_t<meta_t, depend_t> , typename tmeta_t = task_t<meta_t, depend_t>::meta_t>
      class incoming_data_t {
        public:
          using task_t = ttask_t;
          using meta_t = tmeta_t;
          upcxx::promise<incoming_data_t *> on_fetch;
          upcxx::future<incoming_data_t *> on_fetch_future;

          //this should be made generic, not specific to sparse matrices
          std::vector< task_t *> target_tasks;

          meta_t in_meta;
          //a pointer to be used by the user if he wants to attach some data
          std::shared_ptr<blockCellBase_t> extra_data;
	      incoming_data_t() : transfered(false),size(0),extra_data(nullptr),landing_zone(nullptr)
#ifdef CUDA_MODE
          ,d_landing_zone(nullptr), is_gpu_block(false), d_size(0)
#endif
          {
            on_fetch_future = on_fetch.get_future();
          };
#ifdef CUDA_MODE
          using dev_ptr = upcxx::global_ptr<char, upcxx::memory_kind::cuda_device>;
          dev_ptr d_landing_zone;
          bool is_gpu_block;
          size_t d_size;
#endif
          

          bool transfered;
          upcxx::global_ptr<char> remote_gptr;
          size_t size;
          char * landing_zone;

          void allocate() {
            if (landing_zone == nullptr) landing_zone = new char[size];
#ifdef CUDA_MODE
            if (d_landing_zone==nullptr) d_landing_zone = symPACK::gpu_allocator.allocate<char>(size);
            if (d_landing_zone==nullptr) is_gpu_block=false;
#endif
          }

          upcxx::future<incoming_data_t *> fetch() {
            if (!transfered) {
              transfered=true;
#ifdef CUDA_MODE
            if (is_gpu_block) {
              upcxx::when_all(
                upcxx::rget(remote_gptr, landing_zone, size),
                upcxx::copy(remote_gptr, d_landing_zone, size)                    
              ).then([this]() {
                  on_fetch.fulfill_result(this);
                  return;});
            } else {
                upcxx::rget(remote_gptr,landing_zone,size).then([this]() {
                    on_fetch.fulfill_result(this);
                    return;});
            }
#else
              upcxx::rget(remote_gptr,landing_zone,size).then([this]() {
                  on_fetch.fulfill_result(this);
                  return;});
#endif
            }
	    return on_fetch_future;
          }

          ~incoming_data_t() {
            if (landing_zone) delete [] landing_zone;
#ifdef CUDA_MODE
            if (d_landing_zone) symPACK::gpu_allocator.deallocate(d_landing_zone);
#endif
          }

      };

    template <typename bin_label_t = Int, typename ttask_t = task_t<meta_t, depend_t> >
      class task_graph_t: public std::unordered_map< bin_label_t,  std::deque< ttask_t * > > {
      };

    template <typename label_t = key_t, typename ttask_t = task_t<meta_t, depend_t> >
      class task_graph_t3: public std::vector< std::unique_ptr<ttask_t> > {
      };


    template <typename ttask_t = task_t<meta_t, depend_t> , typename ttaskgraph_t = task_graph_t<Int,task_t<meta_t, depend_t> > >
      class Scheduler2D{
        public:
          int sp_handle;


#ifdef SP_THREADS

          std::list< int > free_workers_;
          std::vector< int > worker_loads_;
          std::vector<worker_ptrs<ttask_t>> work_pointers_;
          bool assignWork(ttask_t * ptask);

          std::map< typename ttask_t::lock_t *, std::list<ttask_t*> > delayedTasks_;
          std::recursive_mutex scheduler_mutex_;
          std::function<void()> threadInitHandle_;
          std::function<void()> quiesceHandle_;
          std::function<bool(ttask_t *)> extraTaskHandle_;
          Scheduler2D():threadInitHandle_(nullptr),extraTaskHandle_(nullptr),quiesceHandle_(nullptr) {}
          virtual ~Scheduler2D() {
          }
#endif
#ifdef _PRIORITY_QUEUE_AVAIL_
          struct avail_comp{
            bool operator()(ttask_t *& a,ttask_t*& b) { 
              auto tgt_a = std::get<1>(a->_meta);
              auto tgt_b = std::get<1>(b->_meta);
              if ( tgt_a != tgt_b )
                return tgt_a<tgt_b;
              else if ( std::get<2>(a->_meta) != std::get<2>(b->_meta) )
                return std::get<2>(a->_meta) < std::get<2>(b->_meta);
              else
                return std::get<4>(a->_meta)<std::get<4>(b->_meta); 
            }
          };
          std::priority_queue<ttask_t*,std::vector<ttask_t*>,avail_comp> avail_tasks;
#else
          std::list<ttask_t*> avail_tasks;
#endif

#ifdef _PRIORITY_QUEUE_RDY_
          struct rdy_comp{
            bool operator()(ttask_t *& a,ttask_t*& b) { 
              auto tgt_a = std::get<1>(a->_meta);
              auto tgt_b = std::get<1>(b->_meta);
              if ( tgt_a != tgt_b )
                return tgt_a<tgt_b;
              else if ( std::get<2>(a->_meta) != std::get<2>(b->_meta) )
                return std::get<2>(a->_meta) < std::get<2>(b->_meta);
              else
                return std::get<4>(a->_meta)<std::get<4>(b->_meta); 
            }
          };
          std::priority_queue<ttask_t*,std::vector<ttask_t*>,rdy_comp> ready_tasks;
#else
          std::list<ttask_t*> ready_tasks;
#endif

          void execute(ttaskgraph_t & graph, double & mem_budget );
      };

  }

  extern std::map<int, symPACKMatrixBase *  > g_sp_handle_to_matrix;

  template <typename colptr_t, typename rowind_t, typename T, typename int_t = int> 
    class blockCell_t: public blockCellBase_t {
      protected:
        using intrank_t = upcxx::intrank_t;

      public:
        class block_t {
          public:
            //this is just an offset from the begining of the cell;
            size_t offset;
            rowind_t first_row;
        };
        //virtual int I() {return i;}

#ifdef SP_THREADS
        std::atomic<bool> in_use;
#endif
#ifdef CUDA_MODE
        using dev_ptr = upcxx::global_ptr<T, upcxx::memory_kind::cuda_device>;
        dev_ptr _d_nzval;
        bool is_gpu_block;
#endif
        rowind_t first_col;
        std::tuple<rowind_t> _dims;
        upcxx::global_ptr< char > _gstorage;
        char * _storage;

        T* _nzval;

        size_t _cblocks; //capacity
        size_t _cnz; //capacity
        size_t _nnz;
        size_t _storage_size;
        
        rowind_t _total_rows;

        //relevant only for LDL but more convenient here

        class block_container_t {
          public:
            block_t * _blocks;
            size_t _nblocks;

            block_container_t():_blocks(nullptr),_nblocks(0) {}

            size_t size() const{
              return _nblocks;
            }

            block_t * data() {
              return _blocks;
            }

            class iterator
            {
              public:
                typedef iterator self_type;
                typedef block_t value_type;
                typedef block_t& reference;
                typedef block_t* pointer;
                typedef std::forward_iterator_tag iterator_category;
                typedef int difference_type;
                iterator(pointer ptr) : ptr_(ptr) { }
                self_type operator++() { self_type i = *this; ptr_++; return i; }
                self_type operator++(int junk) { ptr_++; return *this; }
                reference operator*() { return *ptr_; }
                pointer operator->() { return ptr_; }
                bool operator==(const self_type& rhs) { return ptr_ == rhs.ptr_; }
                bool operator!=(const self_type& rhs) { return ptr_ != rhs.ptr_; }
              private:
                pointer ptr_;
            };

            class const_iterator
            {
              public:
                typedef const_iterator self_type;
                typedef block_t value_type;
                typedef block_t const& reference;
                typedef block_t const* pointer;
                typedef int difference_type;
                typedef std::forward_iterator_tag iterator_category;
                const_iterator(pointer ptr) : ptr_(ptr) { }
                self_type operator++() { self_type i = *this; ptr_++; return i; }
                self_type operator++(int junk) { ptr_++; return *this; }
                reference operator*() { return *ptr_; }
                pointer operator->() { return ptr_; }
                bool operator==(const self_type& rhs) { return ptr_ == rhs.ptr_; }
                bool operator!=(const self_type& rhs) { return ptr_ != rhs.ptr_; }
              private:
                pointer ptr_;
            };


            iterator begin()
            {
              return iterator(_blocks);
            }

            iterator end()
            {
              return iterator(_blocks+_nblocks);
            }

            const_iterator begin() const
            {
              return const_iterator(_blocks);
            }

            const_iterator end() const
            {
              return const_iterator(_blocks+_nblocks);
            }

            block_t& operator[](size_t index)
            {
              bassert(index < _nblocks);
              return _blocks[index];
            }

            const block_t& operator[](size_t index) const
            {
              bassert(index < _nblocks);
              return _blocks[index];
            }



        };

        block_container_t _block_container;


        
        blockCell_t(): 
          first_col(0),_dims(std::make_tuple(0)),_total_rows(0),
#ifdef SP_THREADS
          in_use(false),
#endif
          _storage(nullptr),_nzval(nullptr),
#ifdef CUDA_MODE
          _d_nzval(nullptr), is_gpu_block(false),
#endif
	      _gstorage(nullptr),_nnz(0),  _storage_size(0), _own_storage(true) {}
        bool _own_storage;

        virtual ~blockCell_t() {
          if (_own_storage) {
#ifdef CUDA_MODE
            if (is_gpu_block && !_d_nzval.is_null())
              symPACK::gpu_allocator.deallocate(_d_nzval);
#endif
            if ( !_gstorage.is_null() ) {
              bassert(_storage == _gstorage.local());
              upcxx::deallocate( _gstorage );
            }
            else {
              delete [] _storage;
            }
          }
        }

        rowind_t total_rows() {
          return _nnz / width(); 
        }

        //This is only called during FUC and BUC
        blockCell_t ( int_t i, int_t j, rowind_t firstcol, rowind_t width, size_t nzval_cnt, size_t block_cnt, bool shared_segment = true ): blockCell_t() {
            this->i = i;
            this->j = j;
            _dims = std::make_tuple(width);
            first_col = firstcol;
            allocate(nzval_cnt,block_cnt, shared_segment);
	}

	blockCell_t (  int_t i, int_t j,char * ext_storage, rowind_t firstcol, rowind_t width, size_t nzval_cnt, size_t block_cnt ): blockCell_t() {
          this->i = i;
          this->j = j;
          _own_storage = false;
          _gstorage = nullptr;
          _storage = ext_storage;   

          _dims = std::make_tuple(width);
          first_col = firstcol;
          initialize(nzval_cnt,block_cnt);
          _nnz = nzval_cnt;
          _block_container._nblocks = block_cnt;
        }
#ifdef CUDA_MODE
    // GPU block constructor
    blockCell_t (  int_t i, int_t j,char * ext_storage, 
                   upcxx::global_ptr<char, upcxx::memory_kind::cuda_device> d_storage, 
                   rowind_t firstcol, rowind_t width, size_t nzval_cnt, 
                   size_t block_cnt ): blockCell_t() {
          this->i = i;
          this->j = j;
          _own_storage = false;
          _gstorage = nullptr;
          _storage = ext_storage;   
          
          upcxx::global_ptr<block_t, upcxx::memory_kind::cuda_device> d_blocks;
          d_blocks = upcxx::reinterpret_pointer_cast<block_t>(d_storage);
          _d_nzval = upcxx::reinterpret_pointer_cast<T>(d_blocks + block_cnt);
          is_gpu_block = true;
          _dims = std::make_tuple(width);
          first_col = firstcol;
          initialize(nzval_cnt,block_cnt);
          _nnz = nzval_cnt;
          _block_container._nblocks = block_cnt;
        }

#endif

        //THIS CONSTRUCTOR IS NEVER CALLED
        blockCell_t (  int_t i, int_t j,upcxx::global_ptr<char> ext_gstorage, rowind_t firstcol, rowind_t width, size_t nzval_cnt, size_t block_cnt ): blockCell_t() {
          this->i = i;
          this->j = j;
          _own_storage = false;
          _gstorage = ext_gstorage;
          _storage = _gstorage.local();
          _dims = std::make_tuple(width);
          first_col = firstcol;
          initialize(nzval_cnt,block_cnt);
          _nnz = nzval_cnt;
          _block_container._nblocks = block_cnt;
        }

        // Copy constructor.  
        blockCell_t ( const blockCell_t & other ): blockCell_t() {
#ifdef CUDA_MODE
          UPCXX_ASSERT(!other._d_nzval && !other.is_gpu_block, "This function is unimplemented for blocks using host-bypass communication");
#endif
          i           = other.i;
          j           = other.j;
          owner       = other.owner;
          first_col   = other.first_col;
          _dims       = other._dims;
          _total_rows = other._total_rows;
          _own_storage = other._own_storage;
          allocate( other.nz_capacity(), other.block_capacity() , !other._gstorage.is_null() );
          //now copy the data
          std::copy( other._storage, other._storage + other._storage_size, _storage );

          _block_container._nblocks = other._block_container._nblocks;
          _nnz = other._nnz;        
        }

        // Move constructor.  
        blockCell_t ( blockCell_t && other ): blockCell_t() {
#ifdef CUDA_MODE
          UPCXX_ASSERT(!other._d_nzval && !other.is_gpu_block, "This function is unimplemented for blocks using host-bypass communication");
#endif
          i           = other.i;
          j           = other.j;
          owner       = other.owner;
          first_col   = other.first_col;
          _dims       = other._dims;
          _total_rows = other._total_rows;
          _own_storage = other._own_storage;

          _gstorage = other._gstorage;
          _storage = other._storage;

          initialize( other._cnz , other._cblocks );

          //invalidate other
          other._gstorage = nullptr;
          other._storage = nullptr;
          other._storage_size = 0;
          other._cblocks = 0;
          other._cnz = 0;
          other._nnz = 0;
          other._nzval = nullptr;
          other._block_container._nblocks = 0;
          other._block_container._blocks = nullptr;
          other._own_storage = false;
        }




        // Copy assignment operator.  
        blockCell_t& operator=(const blockCell_t& other)  {  
#ifdef CUDA_MODE
          UPCXX_ASSERT(!other._d_nzval && !other.is_gpu_block, "This function is unimplemented for blocks using host-bypass communication");
#endif
          if (this != &other)  
          {  
            // Free the existing resource.  
            if ( ! _gstorage.is_null() ) {
              upcxx::deallocate( _gstorage );
            }
            else {
              delete [] _storage;
            }

            i           = other.i;
            j           = other.j;
            owner       = other.owner;
            first_col   = other.first_col;
            _dims       = other._dims;
            _total_rows = other._total_rows;
            _own_storage = other._own_storage;

            allocate( other._cnz, other._cblocks , !other._gstorage.is_null() );
            //now copy the data
            std::copy( other._storage, other._storage + other._storage_size, _storage );

            _block_container._nblocks = other._block_container._nblocks;
            _nnz = other._nnz;

          } 
          return *this;  
        }  

        // Move assignment operator.  
        blockCell_t& operator=(blockCell_t&& other)  {  
#ifdef CUDA_MODE
          UPCXX_ASSERT(!other._d_nzval && !other.is_gpu_block, "This function is unimplemented for blocks using host-bypass communication");
#endif
          if (this != &other)  
          {  
            // Free the existing resource.  
            if ( !_gstorage.is_null() ) {
              upcxx::deallocate( _gstorage );
            }
            else {
              delete [] _storage;
            }

            i           = other.i;
            j           = other.j;
            owner       = other.owner;
            first_col   = other.first_col;
            _dims       = other._dims;
            _total_rows = other._total_rows;
            _own_storage = other._own_storage;

            _gstorage = other._gstorage;
            _storage = other._storage;

            initialize( other._cnz , other._cblocks );

            //invalidate other
            other._gstorage = nullptr;
            other._storage = nullptr;
            other._storage_size = 0;
            other._cblocks = 0;
            other._cnz = 0;
            other._nnz = 0;
            other._nzval = nullptr; 
            other._block_container._nblocks = 0;
            other._block_container._blocks = nullptr;
            other._own_storage = false;
          } 
          return *this;  
        }  




        inline int_t block_capacity() const { return _cblocks; }
        inline int_t nblocks() const { return _block_container.size(); }
        inline const block_container_t & blocks() const { return _block_container; }
        inline int_t nz_capacity() const { return _cnz; }
        inline int_t nnz() const { return _nnz; }
        inline int_t width() const { return std::get<0>(_dims); }


        //THIS IS VERY DANGEROUS: block HAS to be IN _block_container, not a copy
        inline rowind_t block_nrows(const block_t & block) const{
          auto blkidx = &block - _block_container._blocks;
          bassert(blkidx<nblocks());
          size_t end = (blkidx<nblocks()-1)?_block_container[blkidx+1].offset:_nnz;
          rowind_t nRows = (end-block.offset)/width();
          return nRows;
        }

        inline rowind_t block_nrows(const int_t blkidx) const{
          bassert(blkidx<nblocks());
          size_t end = (blkidx<nblocks()-1)?_block_container[blkidx+1].offset:_nnz;
          rowind_t nRows = (end-_block_container[blkidx].offset)/width();
          return nRows;
        }

        void add_block( rowind_t first_row, rowind_t nrows) {
          bassert( this->_block_container.data()!=nullptr );
          bassert( this->_block_container.size() + 1 <= block_capacity() );
          bassert( nnz() + nrows*std::get<0>(_dims) <= nz_capacity() );
          auto & new_block = _block_container[_block_container._nblocks++];
          new_block.first_row = first_row;
          new_block.offset = _nnz;
          _nnz += nrows*std::get<0>(_dims);
         
        }

        

        void initialize ( size_t nzval_cnt, size_t block_cnt ) {
          _storage_size = nzval_cnt*sizeof(T) + block_cnt*sizeof(block_t);
          _cnz = nzval_cnt;
          _cblocks = block_cnt;
          _nnz = 0;
          _block_container._nblocks = 0;

#ifdef _ALIGNED_
          _nzval = reinterpret_cast<T*>( _storage );
          _block_container._blocks = reinterpret_cast<block_t*>( _nzval + _cnz );
#else
          _block_container._blocks = reinterpret_cast<block_t*>( _storage );
          _nzval = reinterpret_cast<T*>( _block_container._blocks + _cblocks );
#endif


        }

        void allocate ( size_t nzval_cnt, size_t block_cnt, bool shared_segment ) {
          bassert(nzval_cnt!=0 && block_cnt!=0);
 
#ifdef _ALIGNED_
          size_t aligned_size = std::ceil((nzval_cnt*sizeof(T) + block_cnt*sizeof(block_t))/alignof(T))*alignof(T);
#endif

          if ( shared_segment ) {

#ifdef _ALIGNED_
            _gstorage = upcxx::allocate<char, alignof(T) >( aligned_size ); 
#else
            _gstorage = upcxx::allocate<char>( nzval_cnt*sizeof(T) + block_cnt*sizeof(block_t) );
#endif
            _storage = _gstorage.local();
            if ( this->_storage==nullptr ) symPACKOS<<"Trying to allocate "<<nzval_cnt<<" for cell ("<<i<<","<<j<<")"<<std::endl;
            assert( this->_storage!=nullptr );
          }
          else {
            _gstorage = nullptr;

#ifdef _ALIGNED_
            _storage = (char *)new T[aligned_size];
#else
            _storage = new char[nzval_cnt*sizeof(T) + block_cnt*sizeof(block_t)];
#endif
          }
          initialize( nzval_cnt, block_cnt );
        }
	
	//Debugging utility functions

        void print_block(blockCell_t & block, std::string prefix) const{  
          
          logfileptr->OFS()<<prefix<<" ("<<block.i<<","<<block.j<<")"<<std::endl;
          logfileptr->OFS()<<prefix<<" nzval:"<<std::endl;
          for (auto & nzblock: block.blocks() ) {
            // iterate through block_t structs. block.width() is the number of columns in the block
            logfileptr->OFS()<<nzblock.first_row<<"-|"<<block.first_col<<"---------------------"<<block.first_col+block.width()-1<<std::endl;
            // block_nrows(block_t) is used to get the number of rows in the block denoted by block_t
            for (int vrow = 0; vrow< block.block_nrows(nzblock); vrow++) {
              std::streamsize p = logfileptr->OFS().precision();
              logfileptr->OFS().precision(std::numeric_limits< T >::max_digits10);
              // nzval contains the actual data, stored in row major order, each block_t contains a subset of the data starting at block_t.offset
              for (int vcol = 0; vcol< block.width() ; vcol++) {
                logfileptr->OFS()<<std::scientific<<ToMatlabScalar(block._nzval[nzblock.offset+vrow*block.width()+vcol])<<" ";
              }
              logfileptr->OFS().precision(p);
              logfileptr->OFS()<<std::endl;
            }
            logfileptr->OFS()<<nzblock.first_row+block.block_nrows(nzblock)-1<<"-|"<<block.first_col<<"---------------------"<<block.first_col+block.width()-1<<std::endl;
          }
        }  

        template <typename U=T>
        bool check_close(U * correct, U * actual, int n) {
          double epsilon = 1.0;
          logfileptr->OFS()<<"BUFFER SIZE: "<<n<<std::endl;
          bool result = true;
          for (int i=0; i<n; i++) {
            if (fabs(correct[i] - actual[i]) > epsilon && correct[i]!=actual[i]) {
              logfileptr->OFS()<<"Buffers are not equal at index "<<i<<std::endl;
              logfileptr->OFS()<<"Expected "<<correct[i]<<", got "<<actual[i]<<std::endl;
              logfileptr->OFS()<<"Difference: "<<correct[i]-actual[i]<<std::endl;
              result = false;
            }
          }
          if (result)
            logfileptr->OFS()<<"Buffers are equal"<<std::endl;
          return result;
        }

#ifdef CUDA_MODE
        template <typename U=T>
        void check_close(U * correct, upcxx::global_ptr<U, upcxx::memory_kind::cuda_device> actual, int n) {
          double epsilon = 1.0;
          bool iscorrect = true;
          logfileptr->OFS()<<"BUFFER SIZE: "<<n<<std::endl;
          U * h_actual = new U[n];
          CUDA_ERROR_CHECK(cudaMemcpy(h_actual, symPACK::gpu_allocator.local(actual), sizeof(U)*n, cudaMemcpyDeviceToHost));
          cudaDeviceSynchronize();
          for (int i=0; i<n; i++) {
            if (fabs(correct[i] - h_actual[i]) > epsilon && correct[i]!=h_actual[i]) {
              logfileptr->OFS()<<"Buffers are not equal at index "<<i<<std::endl;
              logfileptr->OFS()<<"Expected "<<correct[i]<<", got "<<h_actual[i]<<std::endl;
              logfileptr->OFS()<<"Difference: "<<correct[i]-h_actual[i]<<std::endl;
              iscorrect = false;
            }
          }
          if (iscorrect)
            logfileptr->OFS()<<"Buffers are equal"<<std::endl;
          delete [] h_actual;
        }

        //actual must be a device pointer
        void check_close_dev(T * h_correct, T * actual, int n) {
          double epsilon = 1.0;
          bool iscorrect = true;
          logfileptr->OFS()<<"BUFFER SIZE: "<<n<<std::endl;
          T * h_actual = new T[n];
          CUDA_ERROR_CHECK(cudaMemcpy(h_actual, actual, sizeof(T)*n, cudaMemcpyDeviceToHost));
          for (int i=0; i<n; i++) {
            if (fabs(h_correct[i] - h_actual[i]) > epsilon && h_correct[i]!=h_actual[i]) {
              logfileptr->OFS()<<"Buffers are not equal at index "<<i<<std::endl;
              logfileptr->OFS()<<"Expected "<<h_correct[i]<<", got "<<h_actual[i]<<std::endl;
              logfileptr->OFS()<<"Difference: "<<h_correct[i]-h_actual[i]<<std::endl;
              iscorrect = false;
            }
          }
          if (iscorrect)
            logfileptr->OFS()<<"Buffers are equal"<<std::endl;
          delete [] h_actual;
        }
#endif
	
        void dump_nnz(std::string prelude) {
            progressptr->OFS()<<prelude<<std::endl;
            for (int i=0; i<_nnz; i++) {
                progressptr->OFS()<<_nzval[i]<<",";
            }
            progressptr->OFS()<<std::endl;
        }

#ifdef CUDA_MODE
        void terminate_fallback_handler() {
#if UPCXX_SPEC_VERSION >= 20230300
            std::string err_msg("Error: GPU device segment of " + 
                        std::to_string(symPACK::gpu_allocator.segment_size()) + 
                        " bytes exhausted with " + std::to_string(symPACK::gpu_allocator.segment_used()) +                        " bytes in use. Try increasing `-gpu_mem` (see README.md)");
#else
            std::string err_msg("Error: GPU device segment exhausted, try increasing `-gpu_mem` (see README.md)");
#endif
            throw std::runtime_error(err_msg);
        }
    
#endif

        virtual int factorize( TempUpdateBuffers<T> & tmpBuffers) {
          scope_timer(a,blockCell_t::factorize);
#if defined(_NO_COMPUTATION_)
          return 0;
#endif

          auto snode_size = std::get<0>(_dims);
          auto diag_nzval = _nzval;
          try {
            auto potrf_cpu = [snode_size, diag_nzval](){
                symPACK::increment_counter(symPACK::cpu_ops, "potrf");
                lapack::Potrf( 'U', snode_size, diag_nzval, snode_size);
            };
#ifdef CUDA_MODE
            if (snode_size*snode_size > symPACK::potrf_limit) {
                
            dev_ptr d_diag_nzval = symPACK::gpu_allocator.allocate<T>(snode_size*snode_size);
            
            /* Couldn't allocate enough space on GPU, so use CPU */
            if (d_diag_nzval==nullptr) {
                switch(symPACK::fallback_type) {
                    case FallbackType::TERMINATE :
                        terminate_fallback_handler();
                        break;
                    case FallbackType::CPU :
                        potrf_cpu();
                        break; 
                }
            } else {
                upcxx::copy(diag_nzval, d_diag_nzval, snode_size*snode_size).wait();	
                lapack::cusolver_potrf(symPACK::cusolver_handler, 'U', snode_size, 
                               symPACK::gpu_allocator.local(d_diag_nzval), snode_size);
                symPACK::increment_counter(symPACK::gpu_ops, "potrf");
                CUDA_ERROR_CHECK(cudaDeviceSynchronize());
                upcxx::copy(d_diag_nzval, diag_nzval, snode_size*snode_size).wait();
                symPACK::gpu_allocator.deallocate(d_diag_nzval);
            }

            } else {
                potrf_cpu();
            }
#else
                potrf_cpu();
#endif
	  }
          catch(const std::runtime_error& e) {
            std::cerr << "Runtime error: " << e.what() << '\n';
            gdb_lock();
          }
          return 0;
        }

        virtual int trsm( const blockCellBase_t * pdiag, TempUpdateBuffers<T> & tmpBuffers) {
#if defined(_NO_COMPUTATION_)
          return 0;
#endif

          bassert(dynamic_cast<const blockCell_t*>(pdiag));
          auto diag = dynamic_cast<const blockCell_t*>(pdiag);

          bassert(diag->nblocks()>0);

          auto snode_size = std::get<0>(_dims);
          auto diag_nzval = diag->_nzval;
          auto nzblk_nzval = _nzval;
          int n_rows = total_rows();
          auto trsm_cpu = [snode_size, n_rows, diag_nzval, nzblk_nzval](){
            symPACK::increment_counter(symPACK::cpu_ops, "trsm");
            blas::Trsm('L','U','T','N',
                        snode_size, n_rows, 
                        T(1.0),  diag_nzval, snode_size, 
                        nzblk_nzval, snode_size);
          };
#ifdef CUDA_MODE
          if (diag->is_gpu_block || snode_size*total_rows() > symPACK::trsm_limit) {
            
            dev_ptr diag_nzval_inter, nzblk_nzval_inter;
            nzblk_nzval_inter = symPACK::gpu_allocator.allocate<T>(snode_size*total_rows());                

            if (diag->is_gpu_block) {
                diag_nzval_inter = diag->_d_nzval;
            } else {
                diag_nzval_inter = symPACK::gpu_allocator.allocate<T>(snode_size*snode_size);
            }
            
            if (diag_nzval_inter==nullptr || nzblk_nzval_inter==nullptr) {
                switch(symPACK::fallback_type) {
                    case FallbackType::TERMINATE :
                        terminate_fallback_handler();
                        break;
                    case FallbackType::CPU :
                        trsm_cpu();
                        break; 
                }
            } else {
                
                if (diag->is_gpu_block) {
                    upcxx::copy(nzblk_nzval, nzblk_nzval_inter, _nnz).wait();
                } else {
                    upcxx::when_all(
                        upcxx::copy(nzblk_nzval, nzblk_nzval_inter, _nnz),
                        upcxx::copy(diag_nzval, diag_nzval_inter, diag->_nnz)
                    ).wait();
                }
                
                cublas::cublas_trsm_wrapper2(CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT,
                          snode_size, total_rows(), T(1.0), symPACK::gpu_allocator.local(diag_nzval_inter), snode_size, 
                          symPACK::gpu_allocator.local(nzblk_nzval_inter), snode_size);
                symPACK::increment_counter(symPACK::gpu_ops, "trsm");
                CUDA_ERROR_CHECK(cudaDeviceSynchronize());
                
                upcxx::copy(nzblk_nzval_inter, nzblk_nzval, snode_size*total_rows()).wait();
            }

            if (!diag->is_gpu_block) {
                symPACK::gpu_allocator.deallocate(diag_nzval_inter);
            }
            symPACK::gpu_allocator.deallocate(nzblk_nzval_inter);

          } else {
            trsm_cpu();
          }
#else
            trsm_cpu();
#endif
          return 0;
        }

        virtual int update( blockCellBase_t * ppivot, blockCellBase_t * pfacing, TempUpdateBuffers<T> & tmpBuffers, T* diag = nullptr) {
#if defined(_NO_COMPUTATION_)
          return 0;
#endif
          //do the owner compute update first
          {
            bassert(dynamic_cast<blockCell_t*>(ppivot));
            bassert(dynamic_cast<blockCell_t*>(pfacing));
            blockCell_t & pivot = *dynamic_cast<blockCell_t*>(ppivot);
            blockCell_t & facing = *dynamic_cast<blockCell_t*>(pfacing);
            bassert(nblocks()>0);

            bassert(pivot.nblocks()>0);
            auto pivot_fr = pivot._block_container[0].first_row;

            bassert(facing.nblocks()>0);
            auto facing_fr = facing._block_container[0].first_row;
            auto facing_lr = facing._block_container[facing.nblocks()-1].first_row+facing.block_nrows(facing.nblocks()-1)-1;

            auto src_snode_size = std::get<0>(pivot._dims);
            auto tgt_snode_size = std::get<0>(_dims);


            //find the first row updated by src_snode
            auto tgt_fc = pivot_fr;
            auto tgt_lc = pivot._block_container[pivot._block_container._nblocks-1].first_row
              + pivot.block_nrows(pivot._block_container._nblocks-1) -1;

            int_t first_pivot_idx = 0;
            int_t last_pivot_idx = pivot.nblocks()-1;

            //determine the first column that will be updated in the target supernode
            rowind_t tgt_local_fc =  tgt_fc - this->first_col;
            rowind_t tgt_local_lc =  tgt_lc - this->first_col;

            rowind_t pivot_nrows = pivot.total_rows();
            rowind_t tgt_nrows = this->total_rows();
            rowind_t src_nrows = facing.total_rows();

            //condensed update width is the number of rows in the pivot block
            int_t tgt_width = pivot_nrows;

            T * pivot_nzval = pivot._nzval;
            T * facing_nzval = facing._nzval;
            T * tgt = this->_nzval;

            //Pointer to the output buffer of the GEMM
            T * buf = nullptr;
            T beta = T(0);

            //If the target supernode has the same structure,
            //The GEMM is directly done in place
            size_t tgt_offset = 0;
            bool in_place = ( first_pivot_idx == last_pivot_idx );

            if (in_place) { 
              int tgt_first_upd_blk_idx  = 0;
              for ( ; tgt_first_upd_blk_idx < this->_block_container.size(); tgt_first_upd_blk_idx++ ) {
                auto & block = this->_block_container[tgt_first_upd_blk_idx];
                if (facing_fr >= block.first_row && facing_fr <= block.first_row + this->block_nrows(block) -1) break;
              }

              tgt_offset = this->_block_container[tgt_first_upd_blk_idx].offset
                + (facing_fr - this->_block_container[tgt_first_upd_blk_idx].first_row) * this->width() 
                + tgt_local_fc; 

              //find the last block updated
              int tgt_last_upd_blk_idx  = this->_block_container.size()-1;
              for ( ; tgt_last_upd_blk_idx > tgt_first_upd_blk_idx; tgt_last_upd_blk_idx-- ) {
                auto & block = this->_block_container[tgt_last_upd_blk_idx];
                if (facing_lr >= block.first_row && facing_lr <= block.first_row + this->block_nrows(block) -1) break;
              }
              //make sure that in between these two blocks, everything matches
              int upd_blk_cnt = tgt_last_upd_blk_idx - tgt_first_upd_blk_idx +1;
              if ( in_place && upd_blk_cnt == facing.nblocks() && upd_blk_cnt>=1) {
                for ( int blkidx = tgt_first_upd_blk_idx; blkidx <= tgt_last_upd_blk_idx; blkidx++) {
                  int facingidx = blkidx - tgt_first_upd_blk_idx;
                  int f_fr = facing._block_container[facingidx].first_row; 
                  int f_lr = facing.block_nrows(facingidx) + f_fr -1;
                  int t_fr = std::max(facing_fr, this->_block_container[blkidx].first_row); 
                  int t_lr = std::min(facing_lr, this->_block_container[blkidx].first_row+this->block_nrows(blkidx)-1); 
                  if (f_fr != t_fr || f_lr != t_lr) {
                    in_place = false;
                    break;
                  }
                }
              }
              else {
                in_place = false;
              }
            }
            int ldbuf = tgt_width;
            if (in_place) {
              buf = tgt + tgt_offset;
              beta = T(1);
              ldbuf = tgt_snode_size;
            }
            else {
              //Compute the update in a temporary buffer
              tmpBuffers.tmpBuf.resize(tgt_width*src_nrows + src_snode_size*tgt_width);
              buf = &tmpBuffers.tmpBuf[0];
            }

            if ( in_place && this->i == this->j ) {
              bassert(src_nrows==tgt_width);
              SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_SYRK);
              auto syrk_cpu = [tgt_width, pivot_nzval, src_snode_size, beta, buf, ldbuf]() { 
                    symPACK::increment_counter(symPACK::cpu_ops, "syrk");
                    blas::Syrk('U','T',tgt_width, src_snode_size, 
                                T(-1.0), pivot_nzval, src_snode_size, 
                                beta, buf, ldbuf);
              };
#ifdef CUDA_MODE
              if (tgt_width * src_snode_size > symPACK::syrk_limit) {
            
                dev_ptr d_pivot_nzval = symPACK::gpu_allocator.allocate<T>(tgt_width * src_snode_size);
                dev_ptr d_buf = symPACK::gpu_allocator.allocate<T>(tgt_width * ldbuf);
                
                if (d_pivot_nzval==nullptr || d_buf==nullptr) {
                    switch(symPACK::fallback_type) {
                        case FallbackType::TERMINATE :
                            terminate_fallback_handler();
                            break;
                        case FallbackType::CPU :
                            syrk_cpu(); 
                            break;
                    }
                } else {
                    upcxx::when_all(
                        upcxx::copy(pivot_nzval, d_pivot_nzval, tgt_width*src_snode_size),
                        upcxx::copy(buf, d_buf, tgt_width * ldbuf)
                    ).wait();

                    cublas::cublas_syrk_wrapper2(CUBLAS_FILL_MODE_UPPER,CUBLAS_OP_T,tgt_width, src_snode_size,
                            T(-1.0), symPACK::gpu_allocator.local(d_pivot_nzval), 
                            src_snode_size, beta, 
                            symPACK::gpu_allocator.local(d_buf), ldbuf);
                    symPACK::increment_counter(symPACK::gpu_ops, "syrk");
                    CUDA_ERROR_CHECK(cudaDeviceSynchronize());

                    upcxx::copy(d_buf, buf, tgt_width * ldbuf).wait();

                }
                symPACK::gpu_allocator.deallocate(d_pivot_nzval);
                symPACK::gpu_allocator.deallocate(d_buf);
              } else {
                syrk_cpu();
              }
#else
              syrk_cpu();
#endif                  
              SYMPACK_TIMER_SPECIAL_STOP(UPDATE_SNODE_SYRK);
            }
            else {
              //everything is in row-major
              SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_GEMM);
              auto gemm_cpu = [tgt_width, src_nrows, src_snode_size, 
                               pivot_nzval, facing_nzval, beta, buf, ldbuf]() {
                    symPACK::increment_counter(symPACK::cpu_ops, "gemm");
                    blas::Gemm('T','N',tgt_width, src_nrows,src_snode_size,
                                T(-1.0),pivot_nzval,src_snode_size,
                                facing_nzval,src_snode_size,beta,buf,ldbuf);
              };
#ifdef CUDA_MODE
              if (tgt_width * src_nrows > symPACK::gemm_limit ||
                  tgt_width * src_snode_size > symPACK::gemm_limit ||
                  src_nrows * src_snode_size > symPACK::gemm_limit) {
                dev_ptr d_pivot_nzval = symPACK::gpu_allocator.allocate<T>(tgt_width * src_snode_size);
                dev_ptr d_facing_nzval = symPACK::gpu_allocator.allocate<T>(src_nrows * src_snode_size);
                dev_ptr d_buf = symPACK::gpu_allocator.allocate<T>(ldbuf * src_nrows);
                if (d_pivot_nzval==nullptr || d_facing_nzval==nullptr || d_buf==nullptr) {
                    switch(symPACK::fallback_type) {
                        case FallbackType::TERMINATE :
                            terminate_fallback_handler();
                            break;
                        case FallbackType::CPU :
                            gemm_cpu();
                            break; 
                    }
                } else {
                    upcxx::when_all(
                        upcxx::copy(pivot_nzval, d_pivot_nzval, tgt_width*src_snode_size),
                        upcxx::copy(facing_nzval, d_facing_nzval, src_nrows*src_snode_size),
                        upcxx::copy(buf, d_buf, ldbuf*src_nrows)
                    ).wait();
                    
                    cublas::cublas_gemm_wrapper2(CUBLAS_OP_T, CUBLAS_OP_N,
                              tgt_width, src_nrows, src_snode_size,
                              T(-1.0), 
                              symPACK::gpu_allocator.local(d_pivot_nzval), src_snode_size,
                              symPACK::gpu_allocator.local(d_facing_nzval), src_snode_size, beta, 
                              symPACK::gpu_allocator.local(d_buf), ldbuf);
                    symPACK::increment_counter(symPACK::gpu_ops, "gemm");
                    CUDA_ERROR_CHECK(cudaDeviceSynchronize());
                    

                    upcxx::copy(d_buf, buf, ldbuf*src_nrows).wait();
                    
                }

                symPACK::gpu_allocator.deallocate(d_pivot_nzval);
                symPACK::gpu_allocator.deallocate(d_facing_nzval);
                symPACK::gpu_allocator.deallocate(d_buf);

              } else {
                gemm_cpu();
              }
#else
                gemm_cpu();
#endif                  
              SYMPACK_TIMER_SPECIAL_STOP(UPDATE_SNODE_GEMM);
            }

            //If the GEMM wasn't done in place we need to aggregate the update
            //This is the assembly phase
            if (!in_place) {
              SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_INDEX_MAP);
              tmpBuffers.src_colindx.resize(tgt_width);
              tmpBuffers.src_to_tgt_offset.resize(src_nrows);
              colptr_t colidx = 0;
              colptr_t rowidx = 0;
              size_t offset = 0;
              block_t * tgt_ptr = _block_container._blocks;   
              for ( auto & cur_block: facing.blocks() ) {
                rowind_t cur_src_nrows = facing.block_nrows(cur_block);
                rowind_t cur_src_lr = cur_block.first_row + cur_src_nrows -1;
                rowind_t cur_src_fr = cur_block.first_row;

                //The other one MUST reside into a single block in the target
                rowind_t row = cur_src_fr;
                while (row<=cur_src_lr) {
                  do {
                    if (tgt_ptr->first_row <= row 
                        && row< tgt_ptr->first_row + block_nrows(*tgt_ptr) ) {
                      break;
                    }                 
                  } while ( ++tgt_ptr<_block_container._blocks + _block_container._nblocks );                  

                  int_t lr = std::min(cur_src_lr,tgt_ptr->first_row + block_nrows(*tgt_ptr)-1);
                  int_t tgtOffset = tgt_ptr->offset + (row - tgt_ptr->first_row)*tgt_snode_size;
		  

                  for (int_t cr = row ;cr<=lr;++cr) {
                    offset+=tgt_width;
                    tmpBuffers.src_to_tgt_offset[rowidx] = tgtOffset + (cr - row)*tgt_snode_size;
                    rowidx++;
                  }
                  row += (lr-row+1);
                }
              }

              for ( auto & cur_block: pivot.blocks() ) {
                rowind_t cur_src_ncols = pivot.block_nrows(cur_block);
                rowind_t cur_src_lc = std::min(cur_block.first_row + cur_src_ncols -1, this->first_col+this->width()-1);
                rowind_t cur_src_fc = std::max(cur_block.first_row,this->first_col);
                for (rowind_t col = cur_src_fc ;col<=cur_src_lc;++col) {
                  bassert(this->first_col <= col && col< this->first_col+this->width() );
                  tmpBuffers.src_colindx[colidx++] = col;
                }
              }

              //Multiple cases to consider
              SYMPACK_TIMER_SPECIAL_STOP(UPDATE_SNODE_INDEX_MAP);
              SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_ADD);
              if (first_pivot_idx==last_pivot_idx) {
                // Updating contiguous columns
                rowind_t tgt_offset = (tgt_fc - this->first_col);
                for (rowind_t rowidx = 0; rowidx < src_nrows; ++rowidx) {
                 
                  T * A = &buf[rowidx*tgt_width];
                  T * B = &tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_offset];                 
    			  blas::Axpy(tgt_width, T(1.0), A, 1, B, 1);
                }
              }
              else {
                // full sparse case
                for (rowind_t rowidx = 0; rowidx < src_nrows; ++rowidx) {
                  for (colptr_t colidx = 0; colidx< tmpBuffers.src_colindx.size();++colidx) {
                    rowind_t col = tmpBuffers.src_colindx[colidx];
                    rowind_t tgt_colidx = col - this->first_col;           
                    tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_colidx] 
                      += buf[rowidx*tgt_width+colidx];                       
                  }
                }
              }
              SYMPACK_TIMER_SPECIAL_STOP(UPDATE_SNODE_ADD);
            }
          }
          return 0;
        }

        virtual void copy_row_structure( int_t width, blockCell_t * other ) { 
          this->i = other->i;
          this->j = other->j;
          //check if the block has been initialized yet
          int_t nnz = width * other->total_rows();
          this->_dims = std::make_tuple(width);
          this->first_col = other->first_col;
          this->allocate(nnz,other->nblocks(),true);
          this->initialize(nnz,other->nblocks());
     	  for ( auto & cur_block: other->blocks() ) {
            this->add_block(cur_block.first_row,other->block_nrows(cur_block));          
          }
        }


        virtual int forward_update( blockCellBase_t * psrc_cell) {
          //src_cell is a descendant of the supernode corresponding to *this
          blockCell_t & src_cell = *dynamic_cast<blockCell_t*>(psrc_cell);
          Int ldsol = this->width();
          bassert(ldsol==src_cell.width());

          int_t tgt_blk_idx  = 0;
          for ( auto & src_block: src_cell.blocks() ) {
            for ( ; tgt_blk_idx < this->nblocks(); tgt_blk_idx++ ) {
              auto & block = this->_block_container[tgt_blk_idx];
              if (src_block.first_row >= block.first_row && src_block.first_row <= block.first_row + this->block_nrows(block) -1)
                break;
            }

            auto & block = this->_block_container[tgt_blk_idx];
            T * src = src_cell._nzval + src_block.offset;
            T * tgt = this->_nzval + block.offset + (src_block.first_row - block.first_row)*ldsol; 
            for ( int_t row = 0; row < src_cell.block_nrows(src_block); row++) {
              for ( int_t col = 0; col < ldsol; col++) {
                tgt[row*ldsol+col] += src[row*ldsol+col];
              }
            }
          }
          return 0;
        }


        virtual int forward_update_contrib( blockCellBase_t * ptgt_contrib, blockCellBase_t * pdiag_contrib = nullptr) {
          bassert(dynamic_cast<blockCell_t*>(ptgt_contrib));
          blockCell_t & tgt_contrib = *dynamic_cast<blockCell_t*>(ptgt_contrib);

          Int ldsol = tgt_contrib.width();
          Int ldfact = this->width();
          if ( this->i == this->j ) {
            bassert(this->blocks().size()==1);
            auto diag_nzval = this->_nzval;
            auto tgt_nzval = tgt_contrib._nzval;
            auto trsm_cpu = [ldsol, ldfact, diag_nzval, tgt_nzval]() {
                symPACK::increment_counter(symPACK::cpu_ops, "trsm");
                blas::Trsm('R','U','N','N',ldsol,ldfact, T(1.0),  
                            diag_nzval, ldfact, tgt_nzval, ldsol);
            };
#ifdef CUDA_MODE
            if (ldsol*ldfact > symPACK::trsm_limit && symPACK::gpu_solve) {
            
                dev_ptr d_diag_nzval = symPACK::gpu_allocator.allocate<T>(ldfact*ldfact);
                dev_ptr d_nzval = symPACK::gpu_allocator.allocate<T>(ldsol*ldfact);
                if (d_diag_nzval == nullptr || d_nzval == nullptr) { 
                    switch(symPACK::fallback_type) {
                            case FallbackType::TERMINATE :
                                terminate_fallback_handler();
                                break;
                            case FallbackType::CPU :
                                trsm_cpu();
                                break; 
                    }
                } else {
                    upcxx::when_all(
                        upcxx::copy(diag_nzval, d_diag_nzval, ldfact*ldfact),
                        upcxx::copy(tgt_contrib._nzval, d_nzval, ldsol*ldfact)
                    ).wait();

                    cublas::cublas_trsm_wrapper2(CUBLAS_SIDE_RIGHT,CUBLAS_FILL_MODE_UPPER, 
                        CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT,
                        ldsol,ldfact, T(1.0),  
                        symPACK::gpu_allocator.local(d_diag_nzval), ldfact, 
                        symPACK::gpu_allocator.local(d_nzval), ldsol);
                    symPACK::increment_counter(symPACK::gpu_ops, "trsm");
                    CUDA_ERROR_CHECK(cudaDeviceSynchronize());

                    upcxx::copy(d_nzval, tgt_contrib._nzval, ldsol*ldfact).wait();
                    
                }

                symPACK::gpu_allocator.deallocate(d_diag_nzval);
                symPACK::gpu_allocator.deallocate(d_nzval);

            } else {
                trsm_cpu();
            }
#else
            trsm_cpu();
#endif
          }
          else {
            bassert(pdiag_contrib);
            bassert(dynamic_cast<blockCell_t*>(pdiag_contrib));
            blockCell_t & diag_contrib = *dynamic_cast<blockCell_t*>(pdiag_contrib);

            int_t tgt_blk_idx  = 0;
            for ( auto & src_block: this->blocks() ) {
              for ( ; tgt_blk_idx < tgt_contrib.nblocks(); tgt_blk_idx++ ) {
                auto & block = tgt_contrib._block_container[tgt_blk_idx];
                if (src_block.first_row >= block.first_row && 
                    src_block.first_row <= block.first_row + tgt_contrib.block_nrows(block) -1) break;
              }
              if (tgt_blk_idx>=tgt_contrib.nblocks()) tgt_blk_idx=tgt_contrib.nblocks()-1;
                  auto & block = tgt_contrib._block_container[tgt_blk_idx];

                  T * src = this->_nzval + src_block.offset;
                  T * tgt = tgt_contrib._nzval + block.offset + (src_block.first_row - block.first_row)*ldsol; 
                  //Do -L*Y (gemm)
                  int nrows = this->block_nrows(src_block);
                  auto gemm_cpu = [nrows, diag_contrib, ldsol, src, ldfact, tgt](){
                      symPACK::increment_counter(symPACK::cpu_ops, "gemm");
                      blas::Gemm('N','N',ldsol,nrows,ldfact,
                            T(-1.0),diag_contrib._nzval,ldsol,src,ldfact,T(1.0),tgt,ldsol);
                  };
#ifdef CUDA_MODE
                  if ((ldsol * ldfact > symPACK::gemm_limit ||
                       this->block_nrows(src_block)*ldfact > symPACK::gemm_limit ||
                       ldsol * this->block_nrows(src_block) > symPACK::gemm_limit) && symPACK::gpu_solve) {
                      dev_ptr d_nzval = symPACK::gpu_allocator.allocate<T>(ldsol * ldfact);
                      dev_ptr d_src = symPACK::gpu_allocator.allocate<T>(this->block_nrows(src_block)*ldfact);
                      dev_ptr d_tgt = symPACK::gpu_allocator.allocate<T>(this->block_nrows(src_block)*ldsol);

                      if (d_src == nullptr || d_nzval == nullptr || d_tgt == nullptr) { 
                        switch(symPACK::fallback_type) {
                          case FallbackType::TERMINATE :
                            terminate_fallback_handler();
                            break;
                          case FallbackType::CPU :
                            gemm_cpu();
                            break;
                        }
                      } else {

                          upcxx::when_all(
                            upcxx::copy(diag_contrib._nzval, d_nzval, ldsol*ldfact),
                            upcxx::copy(src, d_src, this->block_nrows(src_block)*ldfact),
                            upcxx::copy(tgt, d_tgt, this->block_nrows(src_block)*ldsol)
                          ).wait();

                          cublas::cublas_gemm_wrapper2(CUBLAS_OP_N,CUBLAS_OP_N,ldsol,this->block_nrows(src_block),ldfact,
                            T(-1.0),symPACK::gpu_allocator.local(d_nzval),ldsol,
                            symPACK::gpu_allocator.local(d_src),ldfact,T(1.0),
                            symPACK::gpu_allocator.local(d_tgt),ldsol);
                          symPACK::increment_counter(symPACK::gpu_ops, "gemm");
                          CUDA_ERROR_CHECK(cudaDeviceSynchronize());
                          
                          upcxx::copy(d_tgt, tgt, ldsol*this->block_nrows(src_block)).wait();

                    }

                    symPACK::gpu_allocator.deallocate(d_nzval);
                    symPACK::gpu_allocator.deallocate(d_src);
                    symPACK::gpu_allocator.deallocate(d_tgt);

                  } else {
                    gemm_cpu();
                  }
#else
                  gemm_cpu();
#endif                 
            }
          }
          return 0;
        }

        virtual int back_update( blockCellBase_t * psrc_cell) {
#if defined(_NO_COMPUTATION_)
          return 0;
#endif
          //src_cell is an ancestor of the supernode corresponding to *this
          blockCell_t & src_cell = *dynamic_cast<blockCell_t*>(psrc_cell);
          Int ldsol = this->width();
          bassert(ldsol==src_cell.width());
          bassert(this->total_rows()==src_cell.total_rows());
          blas::Axpy( src_cell.nnz(),1.0,src_cell._nzval, 1, this->_nzval,1);
          return 0;
        }
	
	//NOTE: This is only called with TRANSA=T
        virtual int _tri_solve(char TRANSA, int_t M,int_t N,T ALPHA,T* A,int_t LDA,T* B,int_t LDB) {
          auto trsm_cpu = [TRANSA, M, N, ALPHA, A, LDA, B, LDB] () {
            symPACK::increment_counter(symPACK::cpu_ops, "trsm");
            blas::Trsm('R','U',TRANSA,'N',M,N, ALPHA,  A, LDA, B, LDB);
          };
#ifdef CUDA_MODE
          if (M*N > symPACK::trsm_limit && symPACK::gpu_solve) {
            dev_ptr d_A = symPACK::gpu_allocator.allocate<T>(LDA*N);
            dev_ptr d_B = symPACK::gpu_allocator.allocate<T>(LDB*N);
            
            if (d_A == nullptr || d_B == nullptr) { 
                switch(symPACK::fallback_type) {
                  case FallbackType::TERMINATE :
                    terminate_fallback_handler();
                    break;
                  case FallbackType::CPU :
                    trsm_cpu();
                    break;
                }
                if (symPACK::gpu_verbose) {
                    symPACK::cpu_ops["trsm"]+=1;
                }
            } else {
                upcxx::when_all(
                    upcxx::copy(A, d_A, LDA*N),
                    upcxx::copy(B, d_B, LDB*N)
                ).wait();

                cublas::cublas_trsm_wrapper2(CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_UPPER, 
                                CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT,M,N, ALPHA, 
                                symPACK::gpu_allocator.local(d_A), LDA, 
                                symPACK::gpu_allocator.local(d_B), LDB);
                symPACK::increment_counter(symPACK::gpu_ops, "trsm");
                CUDA_ERROR_CHECK(cudaDeviceSynchronize());

                upcxx::copy(d_B, B, LDB*N).wait();
            }

            symPACK::gpu_allocator.deallocate(d_A);
            symPACK::gpu_allocator.deallocate(d_B);

          } else {
            trsm_cpu();
          }
#else
          trsm_cpu();
#endif
          return 0;
        }

        virtual int back_update_contrib( blockCellBase_t * ptgt_contrib, blockCellBase_t * pdiag_contrib = nullptr) {
#if defined(_NO_COMPUTATION_)
          return 0;
#endif
          bassert(dynamic_cast<blockCell_t*>(ptgt_contrib));
          blockCell_t & tgt_contrib = *dynamic_cast<blockCell_t*>(ptgt_contrib);
          Int ldsol = tgt_contrib.width();
          Int ldfact = this->width();
          if ( this->i == this->j ) {
            bassert(this->blocks().size()==1);
            auto diag_nzval = this->_nzval;
	        auto contrib_nzval = tgt_contrib._nzval;
            this->_tri_solve('T',ldsol,ldfact, T(1.0), diag_nzval, ldfact, contrib_nzval, ldsol);
          }
          else {
            bassert(pdiag_contrib);
            bassert(dynamic_cast<blockCell_t*>(pdiag_contrib));
            blockCell_t & diag_contrib = *dynamic_cast<blockCell_t*>(pdiag_contrib);
            int_t src_blk_idx  = 0;
            for ( auto & fact_block: this->blocks() ) {
              for ( ; src_blk_idx < diag_contrib.nblocks(); src_blk_idx++ ) {
                auto & block = diag_contrib._block_container[src_blk_idx];
                if (fact_block.first_row >= block.first_row && fact_block.first_row <= block.first_row + diag_contrib.block_nrows(block) -1)
                  break;
              }

              auto & block = diag_contrib._block_container[src_blk_idx];
              T * src = diag_contrib._nzval + block.offset + (fact_block.first_row - block.first_row)*ldsol;
              T * fact = this->_nzval+fact_block.offset;
     	      //Do -X*LT (gemm)
              int nrows = this->block_nrows(fact_block);
              auto tgt_nzval = tgt_contrib._nzval;
              auto gemm_cpu = [ldsol, ldfact, nrows, src, fact, tgt_nzval]() {
                  symPACK::increment_counter(symPACK::cpu_ops, "gemm");
                  blas::Gemm('N','T',ldsol,ldfact,nrows, 
                        T(-1.0),src,ldsol,fact,ldfact,T(1.0),tgt_nzval,ldsol);
              };
#ifdef CUDA_MODE
              T alpha = T(-1.0);
              T beta = T(1.0);
              if ((ldsol * ldfact > symPACK::gemm_limit ||
                    this->block_nrows(fact_block)*ldsol > symPACK::gemm_limit ||
                    ldfact * this->block_nrows(fact_block) > symPACK::gemm_limit) && symPACK::gpu_solve) {
		  
                  dev_ptr d_nzval = symPACK::gpu_allocator.allocate<T>(ldsol * ldfact);
                  dev_ptr d_src = symPACK::gpu_allocator.allocate<T>(this->block_nrows(fact_block)*ldsol);
                  dev_ptr d_fact = symPACK::gpu_allocator.allocate<T>(this->block_nrows(fact_block)*ldfact);

                  if (d_src == nullptr || d_nzval == nullptr || d_fact == nullptr) { 
                        switch(symPACK::fallback_type) {
                          case FallbackType::TERMINATE :
                            terminate_fallback_handler();
                            break;
                          case FallbackType::CPU :
                            gemm_cpu();
                            break;
                        }
                  } else {
                      upcxx::when_all(
                        upcxx::copy(tgt_contrib._nzval, d_nzval, ldsol*ldfact),
                        upcxx::copy(src, d_src, this->block_nrows(fact_block)*ldsol),
                        upcxx::copy(fact, d_fact, this->block_nrows(fact_block)*ldfact)
                      ).wait();

                      cublas::cublas_gemm_wrapper2(CUBLAS_OP_N,CUBLAS_OP_T,
                                ldsol,ldfact,this->block_nrows(fact_block),
                                T(-1.0),symPACK::gpu_allocator.local(d_src),ldsol,
                                symPACK::gpu_allocator.local(d_fact),ldfact,T(1.0),
                                symPACK::gpu_allocator.local(d_nzval),ldsol);
                      symPACK::increment_counter(symPACK::gpu_ops, "gemm");
                      CUDA_ERROR_CHECK(cudaDeviceSynchronize());
                      
                      upcxx::copy(d_nzval, tgt_contrib._nzval, ldfact*ldsol).wait();
                }

                symPACK::gpu_allocator.deallocate(d_nzval);
                symPACK::gpu_allocator.deallocate(d_src);
                symPACK::gpu_allocator.deallocate(d_fact);

              } else {
                gemm_cpu();
              }
#else              
            gemm_cpu();
#endif              
            }
          }
          return 0;
        }

    };


  template <typename colptr_t, typename rowind_t, typename T, typename int_t = int> 
    class blockCellLDL_t: public blockCell_t<colptr_t,rowind_t,T,int_t> {
      protected:
        T * _diag;
        rowind_t first_row;

      public:
        //TODO MAKE THAT WORK IN THE COPY CONSTRUCTORS
        T * _bufLDL;
        std::atomic<int> local_pivot;
        using block_t = typename blockCell_t<colptr_t,rowind_t, T>::block_t;

        blockCellLDL_t():blockCell_t<colptr_t,rowind_t, T>(),_bufLDL(nullptr),_diag(nullptr) {
        }

        virtual ~blockCellLDL_t() {
          if ( this->_bufLDL ) delete [] this->_bufLDL;
        }

        blockCellLDL_t (  int_t i, int_t j, rowind_t firstcol, rowind_t width, size_t nzval_cnt, size_t block_cnt, bool shared_segment = true  ): blockCellLDL_t<colptr_t,rowind_t, T>() {
          this->i = i;
          this->j = j;
          this->_dims = std::make_tuple(width);
          this->first_col = firstcol;
          this->allocate(nzval_cnt,block_cnt, shared_segment);
          this->initialize(nzval_cnt,block_cnt);
        }


        blockCellLDL_t (  int_t i, int_t j, rowind_t firstcol, rowind_t width, size_t nzval_cnt, size_t block_cnt, int_t diag_cnt, bool shared_segment = true  ): blockCellLDL_t<colptr_t,rowind_t, T>() {
          this->i = i;
          this->j = j;
          this->_dims = std::make_tuple(width);
          this->first_col = firstcol;
          this->allocate(nzval_cnt,block_cnt, shared_segment);
          this->initialize(nzval_cnt,block_cnt);
        }



        blockCellLDL_t (  int_t i, int_t j,char * ext_storage, rowind_t firstcol, rowind_t width, size_t nzval_cnt, size_t block_cnt, int_t diag_cnt = -1 ): blockCellLDL_t<colptr_t,rowind_t, T>() {
          this->i = i;
          this->j = j;
          this->_own_storage = false;
          this->_gstorage = nullptr;
          this->_storage = ext_storage;
          this->_dims = std::make_tuple(width);
          this->first_col = firstcol;
          this->initialize(nzval_cnt,block_cnt);
          this->_nnz = nzval_cnt;
          this->_block_container._nblocks = block_cnt;
        }

        blockCellLDL_t (  int_t i, int_t j,upcxx::global_ptr<char> ext_gstorage, rowind_t firstcol, rowind_t width, size_t nzval_cnt, size_t block_cnt, int_t diag_cnt = -1 ): blockCellLDL_t<colptr_t,rowind_t, T>() {
          this->i = i;
          this->j = j;
          this->_own_storage = false;
          this->_gstorage = ext_gstorage;
          this->_storage = this->_gstorage.local();
          this->_dims = std::make_tuple(width);
          this->first_col = firstcol;
          this->initialize(nzval_cnt,block_cnt);
          this->_nnz = nzval_cnt;
          this->_block_container._nblocks = block_cnt;
        }

        // Copy constructor.  
        blockCellLDL_t ( const blockCellLDL_t & other ): blockCellLDL_t() {
          this->i           = other.i;
          this->j           = other.j;
          this->owner       = other.owner;
          this->first_col   = other.first_col;
          this->_dims       = other._dims;
          this->_total_rows = other._total_rows;
          this->_own_storage = other._own_storage;

          allocate( other.nz_capacity(), other.block_capacity() , !other._gstorage.is_null() );
          //now copy the data
          std::copy( other._storage, other._storage + other._storage_size, this->_storage );
          this->_block_container._nblocks = other._block_container._nblocks;
          this->_nnz = other._nnz;

        }

        // Move constructor.  
        blockCellLDL_t ( const blockCellLDL_t && other ): blockCellLDL_t() {
          this->i           = other.i;
          this->j           = other.j;
          this->owner       = other.owner;
          this->first_col   = other.first_col;
          this->_dims       = other._dims;
          this->_total_rows = other._total_rows;
          this->_own_storage = other._own_storage;

          this->_gstorage = other._gstorage;
          this->_storage = other._storage;
          initialize( other._cnz , other._cblocks );

          //invalidate other
          other._gstorage = nullptr;
          other._storage = nullptr;
          other._storage_size = 0;
          other._cblocks = 0;
          other._cnz = 0;
          other._nnz = 0;
          other._nzval = nullptr;
          other._diag = nullptr;        
          other._block_container._nblocks = 0;
          other._block_container._blocks = nullptr;
          other._own_storage = false;
        }

        // Copy assignment operator.  
        blockCellLDL_t& operator=(const blockCellLDL_t& other)  {  
          if (this != &other)  
          {  
            // Free the existing resource.  
            if ( ! this->_gstorage.is_null() ) {
              upcxx::deallocate( this->_gstorage );
            }
            else {
              delete [] this->_storage;
            }

            this->i           = other.i;
            this->j           = other.j;
            this->owner       = other.owner;
            this->first_col   = other.first_col;
            this->_dims       = other._dims;
            this->_total_rows = other._total_rows;
            this->_own_storage = other._own_storage;

            allocate( other._cnz, other._cblocks , !other._gstorage.is_null() );
            //now copy the data
            std::copy( other._storage, other._storage + other._storage_size, this->_storage );
            this->_block_container._nblocks = other._block_container._nblocks;
            this->_nnz = other._nnz;
          } 
          return *this;  
        }  

        // Move assignment operator.  
        blockCellLDL_t& operator=(const blockCellLDL_t&& other)  {  
          if (this != &other)  
          {  
            // Free the existing resource.  
            if ( !this->_gstorage.is_null() ) {
              upcxx::deallocate( this->_gstorage );
            }
            else {
              delete [] this->_storage;
            }

            this->i           = other.i;
            this->j           = other.j;
            this->owner       = other.owner;
            this->first_col   = other.first_col;
            this->_dims       = other._dims;
            this->_total_rows = other._total_rows;
            this->_own_storage = other._own_storage;

            this->_gstorage = other._gstorage;
            this->_storage = other._storage;
            initialize( other._cnz , other._cblocks );

            //invalidate other
            other._gstorage = nullptr;
            other._storage = nullptr;
            other._storage_size = 0;
            other._cblocks = 0;
            other._cnz = 0;
            other._nnz = 0;
            other._nzval = nullptr;        
            other._diag = nullptr;        
            other._block_container._nblocks = 0;
            other._block_container._blocks = nullptr;
            other._own_storage = false;

          } 
          return *this;  
        }  

        T * GetDiag() { return this->_diag;}

        upcxx::global_ptr<char> Diag() {
          auto ptr = upcxx::try_global_ptr((char*)this->_diag);
#ifndef _NDEBUG_
          if ( ptr.is_null() )
            throw std::runtime_error( "blockCellLDL_t must be allocated on the shared segment to use this function.");
#endif
          return ptr;
        }

        void initialize ( size_t nzval_cnt, size_t block_cnt ) {
          ((blockCell_t<colptr_t,rowind_t, T>*)this)->initialize(nzval_cnt, block_cnt);

          if ( this->i == this->j ) {
            this->_storage_size += this->width()*sizeof(T);
            this->_diag = reinterpret_cast<T*>( this->_nzval + nzval_cnt );
#ifdef _ALIGNED_
            this->_block_container._blocks = reinterpret_cast<block_t*>( this->_diag + this->width() );
#endif
          }
          else {
            this->_diag = nullptr;
          }
        }

        void allocate ( size_t nzval_cnt, size_t block_cnt, bool shared_segment ) {
          bassert(nzval_cnt!=0 && block_cnt!=0);

          size_t aligned_size = 0;
          if ( shared_segment ) {
            if ( this->i == this->j ) {
#ifdef _ALIGNED_
              aligned_size = std::ceil(( (nzval_cnt+this->width())*sizeof(T) + block_cnt*sizeof(block_t))/alignof(T))*alignof(T);
#else
              this->_gstorage = upcxx::allocate<char>( nzval_cnt*sizeof(T) + block_cnt*sizeof(block_t) + this->width()*sizeof(T) );
#endif
            }
            else {
#ifdef _ALIGNED_
              aligned_size = std::ceil((nzval_cnt*sizeof(T) + block_cnt*sizeof(block_t))/alignof(T))*alignof(T);
#else
              this->_gstorage = upcxx::allocate<char>( nzval_cnt*sizeof(T) + block_cnt*sizeof(block_t) );
#endif
            }
#ifdef _ALIGNED_
            this->_gstorage = upcxx::allocate<char, alignof(T) >( aligned_size );
#endif
            this->_storage = this->_gstorage.local();
          }
          else {
            this->_gstorage = nullptr;
            if ( this->i == this->j ) {
#ifdef _ALIGNED_
              aligned_size = std::ceil(( (nzval_cnt+this->width())*sizeof(T) + block_cnt*sizeof(block_t))/alignof(T))*alignof(T);
#else
              this->_storage = new char[nzval_cnt*sizeof(T) + block_cnt*sizeof(block_t) + this->width()*sizeof(T)];
#endif
            }
            else {
#ifdef _ALIGNED_
              aligned_size = std::ceil((nzval_cnt*sizeof(T) + block_cnt*sizeof(block_t))/alignof(T))*alignof(T);
#else
              this->_storage = new char[nzval_cnt*sizeof(T) + block_cnt*sizeof(block_t)];
#endif
            }
#ifdef _ALIGNED_
            this->_storage = (char *)new T[aligned_size];
#endif
          }
          assert(this->_storage!= nullptr);
          initialize( nzval_cnt, block_cnt );
        }


        void print_block(const blockCellLDL_t & block, std::string prefix) const{  
          logfileptr->OFS()<<prefix<<" ("<<block.i<<","<<block.j<<")"<<std::endl;
          logfileptr->OFS()<<prefix<<" nzval:"<<std::endl;
          for (auto & nzblock: block.blocks() ) {
            logfileptr->OFS()<<nzblock.first_row<<"-|"<<block.first_col<<"---------------------"<<block.first_col+block.width()-1<<std::endl;
            for (int vrow = 0; vrow< block.block_nrows(nzblock); vrow++) {
              std::streamsize p = logfileptr->OFS().precision();
              logfileptr->OFS().precision(std::numeric_limits< T >::max_digits10);
              for (int vcol = 0; vcol< block.width() ; vcol++) {
                logfileptr->OFS()<<std::scientific<<ToMatlabScalar(block._nzval[nzblock.offset+vrow*block.width()+vcol])<<" ";
              }
              logfileptr->OFS().precision(p);
              logfileptr->OFS()<<std::endl;
            }
            logfileptr->OFS()<<nzblock.first_row+block.block_nrows(nzblock)-1<<"-|"<<block.first_col<<"---------------------"<<block.first_col+block.width()-1<<std::endl;
          }
        }  

        virtual void copy_row_structure( int_t width, blockCell_t<colptr_t,rowind_t, T> * other ) override { 
          this->i = other->i;
          this->j = other->j;

          //check if the block has been initialized yet
          int_t nnz = width * other->total_rows();
          this->_dims = std::make_tuple(width);
          this->first_col = other->first_col;

          this->allocate(nnz,other->nblocks(),true);
          this->initialize(nnz,other->nblocks());
          for ( auto & cur_block: other->blocks() ) {
            this->add_block(cur_block.first_row,other->block_nrows(cur_block));
          }
        }




        virtual int factorize( TempUpdateBuffers<T> & tmpBuffers) override{
          scope_timer(a,blockCellLDL_t::factorize);
#if defined(_NO_COMPUTATION_)
          return 0;
#endif
          auto snode_size = std::get<0>(this->_dims);
          auto diag_nzval = this->_nzval;

          int_t INFO = 0;
          try{
            lapack::Potrf_LDL( "U", snode_size, diag_nzval, snode_size, tmpBuffers, INFO);
            //copy the diagonal entries into the diag portion of the supernode
#pragma unroll
            for (int_t col = 0; col< snode_size; col++) { this->_diag[col] = diag_nzval[col+ (col)*snode_size]; }
          }
          catch(const std::runtime_error& e) {
            std::cerr << "Runtime error: " << e.what() << '\n';
            gdb_lock();
          }
          return 0;
        }

        virtual int trsm( const blockCellBase_t * pdiag, TempUpdateBuffers<T> & tmpBuffers) override {
          scope_timer(a,blockCellLDL_t::trsm);
#if defined(_NO_COMPUTATION_)
          return 0;
#endif
          bassert(dynamic_cast<const blockCellLDL_t*>(pdiag));
          auto diag = dynamic_cast<const blockCellLDL_t*>(pdiag);

          bassert(diag->nblocks()>0);
          auto diag_nzval = diag->_nzval;

          auto snode_size = std::get<0>(this->_dims);
          auto nzblk_nzval = this->_nzval;
          blas::Trsm('L','U','T','U',snode_size, this->total_rows(), T(1.0),  diag_nzval, snode_size, nzblk_nzval, snode_size);

          //scale column I
          for ( int_t I = 1; I<=snode_size; I++) {
            blas::Scal( this->total_rows(), T(1.0)/diag->_diag[I-1], &nzblk_nzval[I-1], snode_size );
          }

          return 0;
        }


        int computeDLT( T * diag ) {
#if defined(_NO_COMPUTATION_)
          return 0;
#endif
          auto snode_size = this->width();
          // W = DLT 
          int_t tgt_width = this->total_rows();
          bool computeW = this->_bufLDL == nullptr;
          bassert(computeW);
          this->_bufLDL = new T[snode_size*tgt_width];
          auto bufLDL = this->_bufLDL;


          for ( int_t row = 0; row < snode_size; row++ ) {
            for (int_t col = 0; col<tgt_width; col++) {
              bufLDL[row*tgt_width+col] = diag[row]*this->_nzval[row+col*snode_size];
            }
          }
          return 0;
        }



        virtual int update(  blockCellBase_t * ppivot, blockCellBase_t * pfacing, TempUpdateBuffers<T> & tmpBuffers, T* diag = nullptr) override {
          scope_timer(a,blockCellLDL_t::update);
#if defined(_NO_COMPUTATION_)
          return 0;
#endif
          //do the owner compute update first

          {
            bassert(dynamic_cast<blockCellLDL_t*>(ppivot));
            bassert(dynamic_cast<blockCellLDL_t*>(pfacing));

            blockCellLDL_t & pivot = *dynamic_cast<blockCellLDL_t*>(ppivot);
            blockCellLDL_t & facing = *dynamic_cast<blockCellLDL_t*>(pfacing);


            bassert(this->nblocks()>0);

            bassert(pivot.nblocks()>0);
            auto pivot_fr = pivot._block_container[0].first_row;

            bassert(facing.nblocks()>0);
            auto facing_fr = facing._block_container[0].first_row;
            auto facing_lr = facing._block_container[facing.nblocks()-1].first_row+facing.block_nrows(facing.nblocks()-1)-1;

            auto src_snode_size = std::get<0>(pivot._dims);
            auto tgt_snode_size = std::get<0>(this->_dims);


            //find the first row updated by src_snode
            auto tgt_fc = pivot_fr;
            auto tgt_lc = pivot._block_container[pivot._block_container._nblocks-1].first_row
              + pivot.block_nrows(pivot._block_container._nblocks-1) -1;

            int_t first_pivot_idx = 0;
            int_t last_pivot_idx = pivot.nblocks()-1;

            //determine the first column that will be updated in the target supernode
            rowind_t tgt_local_fc =  tgt_fc - this->first_col;
            rowind_t tgt_local_lc =  tgt_lc - this->first_col;

            rowind_t pivot_nrows = pivot.total_rows();
            rowind_t tgt_nrows = this->total_rows();
            rowind_t src_nrows = facing.total_rows();

            //condensed update width is the number of rows in the pivot block
            int_t tgt_width = pivot_nrows;

            T * pivot_nzval = pivot._nzval;
            T * facing_nzval = facing._nzval;
            T * tgt = this->_nzval;

            //Pointer to the output buffer of the GEMM
            T * buf = nullptr;
            T beta = T(0);
            T * bufLDL = nullptr;
            //If the target supernode has the same structure,
            //The GEMM is directly done in place
            size_t tgt_offset = 0;
            bool in_place = ( first_pivot_idx == last_pivot_idx );

            if (in_place) { 
              int tgt_first_upd_blk_idx  = 0;
              for ( ; tgt_first_upd_blk_idx < this->_block_container.size(); tgt_first_upd_blk_idx++ ) {
                auto & block = this->_block_container[tgt_first_upd_blk_idx];
                if (facing_fr >= block.first_row && facing_fr <= block.first_row + this->block_nrows(block) -1)
                  break;
              }

              tgt_offset = this->_block_container[tgt_first_upd_blk_idx].offset
                + (facing_fr - this->_block_container[tgt_first_upd_blk_idx].first_row) * this->width() 
                + tgt_local_fc; 

              //find the last block updated
              int tgt_last_upd_blk_idx  = this->_block_container.size()-1;
              for ( ; tgt_last_upd_blk_idx > tgt_first_upd_blk_idx; tgt_last_upd_blk_idx-- ) {
                auto & block = this->_block_container[tgt_last_upd_blk_idx];
                if (facing_lr >= block.first_row && facing_lr <= block.first_row + this->block_nrows(block) -1)
                  break;
              }
              //make sure that in between these two blocks, everything matches
              int upd_blk_cnt = tgt_last_upd_blk_idx - tgt_first_upd_blk_idx +1;
              if ( in_place && upd_blk_cnt == facing.nblocks() && upd_blk_cnt>=1) {
                for ( int blkidx = tgt_first_upd_blk_idx; blkidx <= tgt_last_upd_blk_idx; blkidx++) {
                  int facingidx = blkidx - tgt_first_upd_blk_idx;
                  int f_fr = facing._block_container[facingidx].first_row; 
                  int f_lr = facing.block_nrows(facingidx) + f_fr -1;
                  int t_fr = std::max(facing_fr, this->_block_container[blkidx].first_row); 
                  int t_lr = std::min(facing_lr, this->_block_container[blkidx].first_row+this->block_nrows(blkidx)-1); 
                  if (f_fr != t_fr || f_lr != t_lr) {
                    in_place = false;
                    break;
                  }
                }
              }
              else {
                in_place = false;
              }
            }

            int ldbuf = tgt_width;

            if (in_place)
            {
              //TODO
              buf = tgt + tgt_offset;
              beta = T(1);
              ldbuf = tgt_snode_size;
            }
            else {
              //Compute the update in a temporary buffer
//#ifdef SP_THREADS
              tmpBuffers.tmpBuf.resize(tgt_width*src_nrows);
//#endif
              buf = &tmpBuffers.tmpBuf[0];
            }

            //everything is in row-major
            SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_GEMM);

            //First do W = DLT 
            bool computeW = pivot._bufLDL == nullptr;
            if ( computeW ) {
              bassert(pivot._own_storage);
              pivot.computeDLT(diag);
            }

            bassert(pivot._bufLDL!=nullptr);
            bufLDL = pivot._bufLDL;

            //Then do -L*W (gemm)
            blas::Gemm('N','N',tgt_width,src_nrows,src_snode_size,
                T(-1.0),bufLDL,tgt_width,facing_nzval,src_snode_size,beta,buf,ldbuf);

            if ( pivot._own_storage ) {
              if ( pivot.local_pivot.fetch_sub(1)==1) {
                delete [] pivot._bufLDL;
                pivot._bufLDL = nullptr;
              }
            }

            SYMPACK_TIMER_SPECIAL_STOP(UPDATE_SNODE_GEMM);

            //If the GEMM wasn't done in place we need to aggregate the update
            //This is the assembly phase
            if (!in_place)
            {
              {
                SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_INDEX_MAP);
                tmpBuffers.src_colindx.resize(tgt_width);
                tmpBuffers.src_to_tgt_offset.resize(src_nrows);
                colptr_t colidx = 0;
                colptr_t rowidx = 0;
                size_t offset = 0;

                block_t * tgt_ptr = this->_block_container._blocks;

                for ( auto & cur_block: facing.blocks() ) {
                  rowind_t cur_src_nrows = facing.block_nrows(cur_block);
                  rowind_t cur_src_lr = cur_block.first_row + cur_src_nrows -1;
                  rowind_t cur_src_fr = cur_block.first_row;

                  //The other one MUST reside into a single block in the target
                  rowind_t row = cur_src_fr;
                  while (row<=cur_src_lr) {
                    do {
                      if (tgt_ptr->first_row <= row 
                          && row< tgt_ptr->first_row + this->block_nrows(*tgt_ptr) ) {
                        break;
                      }
                    } while ( ++tgt_ptr<this->_block_container._blocks + this->_block_container._nblocks ); 

                    int_t lr = std::min(cur_src_lr,tgt_ptr->first_row + this->block_nrows(*tgt_ptr)-1);
                    int_t tgtOffset = tgt_ptr->offset + (row - tgt_ptr->first_row)*tgt_snode_size;

                    for (int_t cr = row ;cr<=lr;++cr) {
                      offset+=tgt_width;
                      tmpBuffers.src_to_tgt_offset[rowidx] = tgtOffset + (cr - row)*tgt_snode_size;
                      rowidx++;
                    }
                    row += (lr-row+1);
                  }
                }

                for ( auto & cur_block: pivot.blocks() ) {
                  rowind_t cur_src_ncols = pivot.block_nrows(cur_block);
                  rowind_t cur_src_lc = std::min(cur_block.first_row + cur_src_ncols -1, this->first_col+this->width()-1);
                  rowind_t cur_src_fc = std::max(cur_block.first_row,this->first_col);

                  for (rowind_t col = cur_src_fc ;col<=cur_src_lc;++col) {
                    bassert(this->first_col <= col && col< this->first_col+this->width() );
                    tmpBuffers.src_colindx[colidx++] = col;
                  }
                }

                //Multiple cases to consider
                SYMPACK_TIMER_SPECIAL_STOP(UPDATE_SNODE_INDEX_MAP);
                SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_ADD);
                if (first_pivot_idx==last_pivot_idx) {
                  // Updating contiguous columns
                  rowind_t tgt_offset = (tgt_fc - this->first_col);
                  for (rowind_t rowidx = 0; rowidx < src_nrows; ++rowidx) {
                    T * A = &buf[rowidx*tgt_width];
                    T * B = &tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_offset];
#pragma unroll
                    for (rowind_t i = 0; i < tgt_width; ++i) { B[i] += A[i]; }
                  }
                }
                else {
                  // full sparse case
                  for (rowind_t rowidx = 0; rowidx < src_nrows; ++rowidx) {
                    for (colptr_t colidx = 0; colidx< tmpBuffers.src_colindx.size();++colidx) {
                      rowind_t col = tmpBuffers.src_colindx[colidx];
                      rowind_t tgt_colidx = col - this->first_col;
                      tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_colidx] 
                        += buf[rowidx*tgt_width+colidx]; 
                    }
                  }
                }

                SYMPACK_TIMER_SPECIAL_STOP(UPDATE_SNODE_ADD);
              }
            }
          }
          return 0;
        }

        bool scaled = false;
        int scale_contrib(blockCellLDL_t * ptgt_contrib) {
#if defined(_NO_COMPUTATION_)
          return 0;
#endif
          blockCellLDL_t & tgt_contrib = *(ptgt_contrib);
          int_t ldsol = tgt_contrib.width();
          auto & block = tgt_contrib.blocks()[0];
          int_t ldfact = tgt_contrib.block_nrows(block);

          bassert(tgt_contrib.i == tgt_contrib.j);
          bassert(!tgt_contrib.scaled);
          for(int_t kk = 0; kk<ldfact; ++kk){
            blas::Scal( ldsol, T(1.0)/this->_diag[kk], &tgt_contrib._nzval[kk*ldsol], 1 );
          }
          tgt_contrib.scaled = true;
          return 0;
        }

        virtual int _tri_solve(char TRANSA, int_t M,int_t N,T ALPHA,T* A,int_t LDA,T* B,int_t LDB) override {
#if defined(_NO_COMPUTATION_)
          return 0;
#endif
          blas::Trsm('R','U',TRANSA,'U',M,N, ALPHA,  A, LDA, B, LDB);
          return 0;
        }


        virtual int forward_update_contrib( blockCellBase_t * ptgt_contrib, blockCellBase_t * pdiag_contrib = nullptr) override {
#if defined(_NO_COMPUTATION_)
          return 0;
#endif
          static bool first = false;

          bassert(dynamic_cast<blockCellLDL_t*>(ptgt_contrib));
          blockCellLDL_t & tgt_contrib = *dynamic_cast<blockCellLDL_t*>(ptgt_contrib);

          bassert(!tgt_contrib.scaled);
          int_t ldsol = tgt_contrib.width();
          int_t ldfact = this->width();

          if ( this->i == this->j ) {
            bassert(this->blocks().size()==1);
            auto diag_nzval = this->_nzval;
            this->_tri_solve('N',ldsol,ldfact, T(1.0),  diag_nzval, ldfact, tgt_contrib._nzval, ldsol);
          }
          else {
            bassert(pdiag_contrib);
            bassert(dynamic_cast<blockCellLDL_t*>(pdiag_contrib));
            blockCellLDL_t & diag_contrib = *dynamic_cast<blockCellLDL_t*>(pdiag_contrib);

            int_t tgt_blk_idx  = 0;
            for ( auto & src_block: this->blocks() ) {
              for ( ; tgt_blk_idx < tgt_contrib.nblocks(); tgt_blk_idx++ ) {
                auto & block = tgt_contrib._block_container[tgt_blk_idx];
                if (src_block.first_row >= block.first_row && 
                    src_block.first_row <= block.first_row + tgt_contrib.block_nrows(block) -1) break;
              }

              auto & block = tgt_contrib._block_container[tgt_blk_idx];
              T * src = this->_nzval + src_block.offset;
              T * tgt = tgt_contrib._nzval + block.offset + (src_block.first_row - block.first_row)*ldsol; 
              auto cur_nrows = this->block_nrows(src_block);

              bassert(diag_contrib.j != tgt_contrib.j);
              //Do -L*Y (gemm)
              blas::Gemm('N','N',ldsol,this->block_nrows(src_block),ldfact, T(-1.0),diag_contrib._nzval,ldsol,src,ldfact,T(1.0),tgt,ldsol);
            }
          }
          return 0;
        }



    };

  template <typename T> class symPACKMatrix;

  template <typename colptr_t, typename rowind_t, typename T, typename int_t = int> 
    class symPACKMatrix2D: public symPACKMatrixMeta<T>{


      friend class symPACKMatrix<T>;

      int dist_np;
      using snodeBlock_t = blockCell_t<colptr_t,rowind_t, T>;
      using snodeBlockLDL_t = blockCellLDL_t<colptr_t,rowind_t, T>;

      using snodeBlockBase_sptr_t = std::shared_ptr<blockCellBase_t>;
      using snodeBlock_sptr_t = std::shared_ptr<snodeBlock_t>;
      using snodeBlockLDL_sptr_t = std::shared_ptr<snodeBlockLDL_t>;

      using SparseTask2D = scheduling::task_t<scheduling::meta_t, scheduling::depend_t>;


      public:
      struct solve_data_t {
        T * rhs;
        int nrhs;
        //This containes AT MOST nsuper contribs.
        std::vector< std::tuple<int_t,snodeBlock_sptr_t> > contribs;
        upcxx::dist_object<int> remoteDeallocCounter;

        void deallocRemote(int owner,int sp_handle, int J, int dep_cnt) {
          int sender = upcxx::rank_me();
          upcxx::rpc_ff(owner,[sender,sp_handle,J,dep_cnt](upcxx::dist_object<int> & dealloc_cnt){ 
              auto matptr = (symPACKMatrix2D<colptr_t,rowind_t,T> *) g_sp_handle_to_matrix[sp_handle];

              auto & counter = std::get<0>(matptr->solve_data.contribs[J]);
              counter-=dep_cnt;
              if ( counter==0 ) {
              //TODO THIS DOES NOT PLAY WELL WITH MULTITHREADING
              std::get<1>(matptr->solve_data.contribs[J]).reset();
              }
              bassert(*dealloc_cnt > 0 );
              (*dealloc_cnt)-=dep_cnt;
              },remoteDeallocCounter);
        }

        std::atomic<bool> * contribs_lock;

        std::vector<int> update_right_cnt;
        std::vector<int> update_up_cnt;

        solve_data_t():remoteDeallocCounter(0) {
          contribs_lock = nullptr;
        }

        ~solve_data_t() {
          delete [] contribs_lock;
        }
      };
      solve_data_t solve_data;

      std::vector< snodeBlock_sptr_t > localBlocks_;
      int nsuper;
      double mem_budget;
      using cell_key_t = std::pair<uint64_t,uint64_t>;
      cell_key_t coord2supidx(uint64_t i, uint64_t j) { 
        cell_key_t val = std::make_pair(i,j);
        return val;
      }
      std::map<cell_key_t, snodeBlockBase_sptr_t > cells_;

      inline snodeBlockBase_sptr_t pQueryCELL (int a, int b)  { 
        auto idx = coord2supidx((a),(b)); 
        auto it = this->cells_.find(idx);
        if (it != this->cells_.end() ) { 
          return it->second;
        }
        else
          return nullptr;
      }

      inline snodeBlock_sptr_t pQueryCELL2 (int a, int b)  { 
#ifdef _TIMING_
        gasneti_tick_t start = gasneti_ticks_now();
#endif
        auto idx = coord2supidx((a),(b)); 
        auto it = this->cells_.find(idx);
        snodeBlock_sptr_t ptr = nullptr;
        if (it != this->cells_.end() ) { 
          ptr = std::dynamic_pointer_cast<snodeBlock_t>(it->second);
          bassert ( ptr != nullptr );
        }

#ifdef _TIMING_
        CELL_ticks += gasneti_ticks_to_ns(gasneti_ticks_now() - start);
#endif

        return ptr;
      }
      inline snodeBlock_sptr_t pCELL (int a, int b) { return std::static_pointer_cast<snodeBlock_t>(this->cells_[coord2supidx((a),(b))]); }

#ifdef _TIMING_
      uint64_t CELL_ticks = 0;
      uint64_t rpc_fact_ticks = 0;
      uint64_t rpc_trsm_ticks = 0;
      uint64_t rpc_upd_ticks = 0;
      uint64_t deps_fact_ticks = 0;
      uint64_t deps_trsm_ticks = 0;
      uint64_t deps_upd_ticks = 0;

      double comp_fact_ticks = 0;
      double comp_trsm_ticks = 0;
      double comp_upd_ticks = 0;
#endif

      inline snodeBlock_t & CELL (int a, int b) {
        auto ptr = pQueryCELL2(a,b); assert(ptr!=nullptr && (ptr->i-1==a && ptr->j-1==b)); 
        return *ptr;  
      }

      protected:


      //TODO ideally, this should be DistSparseMatrixGraph<colptr_t,rowind_t>
#ifdef SP_THREADS
      std::map<std::thread::id,TempUpdateBuffers<T> > tmpBufs_th;
#else
      TempUpdateBuffers<T> tmpBufs;
#endif

      std::vector< std::tuple<Int,upcxx::global_ptr<char>> > diag_pointers_;
      const upcxx::global_ptr<char> & find_diag_pointer( int supno ) const {
        auto it = std::lower_bound(this->diag_pointers_.begin(),this->diag_pointers_.end(), supno, [](const std::tuple<Int,upcxx::global_ptr<char> > & a, const Int & key ) { return std::get<0>(a) < key; });
        bassert(it != this->diag_pointers_.end());
        bassert(supno == std::get<0>(*it));
        return std::get<1>(*it);
      }

      public:
      using TaskGraph2D = scheduling::task_graph_t3<scheduling::key_t,SparseTask2D >; 
      std::vector< std::tuple< scheduling::key_t, std::size_t > > task_idx;
      TaskGraph2D task_graph;

      std::vector< std::tuple< scheduling::key_t, std::size_t > > task_idx_solve;
      TaskGraph2D task_graph_solve;

      scheduling::Scheduler2D<SparseTask2D,TaskGraph2D> scheduler;

      symPACKMatrix2D();
      ~symPACKMatrix2D();


      virtual void DumpMatlab() override;

      virtual void Init(symPACKOptions & options ) override;
      void Init(DistSparseMatrix<T> & pMat, symPACKOptions & options );

      virtual void SymbolicFactorization(DistSparseMatrix<T> & pMat) override;
      virtual void DistributeMatrix(DistSparseMatrix<T> & pMat) override;


      virtual void Factorize() override;
      virtual void Solve( T * RHS, int nrhs, int rhs_size,  T * Xptr = nullptr ) override;
      virtual void GetSolution(T * B, int nrhs) override;



    };

  template <typename colptr_t, typename rowind_t, typename T, typename int_t>
    symPACKMatrix2D<colptr_t,rowind_t,T,int_t>::symPACKMatrix2D():
      symPACKMatrixMeta<T>()
  {
#ifndef NO_MPI
#endif
    g_sp_handle_to_matrix[this->sp_handle] = this;
    this->dist_np = 0;
  }

  template <typename colptr_t, typename rowind_t, typename T, typename int_t>
    symPACKMatrix2D<colptr_t,rowind_t,T,int_t>::~symPACKMatrix2D() {

    }

  template <typename colptr_t, typename rowind_t, typename T, typename int_t>
    void symPACKMatrix2D<colptr_t,rowind_t,T,int_t>::Init(symPACKOptions & options ) {
      scope_timer(a,symPACKMatrix2D::Init);
      this->options_ = options;
      logfileptr->verbose = this->options_.verbose>0;
      if (this->options_.verbose==0) {
        symPACKOS.rdbuf(nullptr);
      }

      this->mem_budget = this->options_.memory_limit;

#ifndef NO_MPI
      this->all_np = 0;
      MPI_Comm_size(this->options_.MPIcomm,&this->all_np);
      MPI_Comm_rank(this->options_.MPIcomm,&this->iam);
#endif

      this->iam = upcxx::rank_me();
      this->all_np = upcxx::rank_n();

      this->np = this->all_np;
      this->dist_np = this->all_np;


#ifndef NO_MPI

      if(this->workcomm_!=MPI_COMM_NULL){
        MPI_Comm_free(&this->workcomm_);
      }
      if(this->non_workcomm_!=MPI_COMM_NULL){
        MPI_Comm_free(&this->non_workcomm_);
      }
      MPI_Comm_split(this->options_.MPIcomm,this->iam<this->np,this->iam,&this->workcomm_);

      //do another split to contain P0 and all the non working processors

      if (this->all_np!=this->np) {
        this->non_workcomm_ = MPI_COMM_NULL;
        MPI_Comm_split(this->options_.MPIcomm,this->iam==0||this->iam>=this->np,this->iam==0?0:this->iam-this->np+1,&this->non_workcomm_);
      }
#else
      //    upcxx::team_all.split(this->iam<this->np,new_rank, this->team_);
#endif

    } 

  template <typename colptr_t, typename rowind_t, typename T, typename int_t>
    void symPACKMatrix2D<colptr_t,rowind_t,T,int_t>::Init(DistSparseMatrix<T> & pMat,symPACKOptions & options ) {
    } 

  template <typename colptr_t, typename rowind_t, typename T, typename int_t>
    void symPACKMatrix2D<colptr_t,rowind_t,T,int_t>::SymbolicFactorization(DistSparseMatrix<T> & pMat ) {
      scope_timer(a,symPACKMatrix2D::SymbolicFactorization);
      std::vector<Int> cc;
      {
        //This has to be declared here to be able to debug ...
        std::vector<int, Mallocator<int> > xadj;
        std::vector<int, Mallocator<int> > adj;
        Idx row;
        int fc,lc,colbeg,colend,col;
        Ptr supbeg,supend,rowidx;
        Int I;

#ifndef NO_MPI
        if (this->fullcomm_!=MPI_COMM_NULL) {
          MPI_Comm_free(&this->fullcomm_);
        }
        MPI_Comm_dup(pMat.comm,&this->fullcomm_);
#endif

        this->iSize_ = pMat.size;

        this->graph_ = pMat.GetLocalGraph();
        this->graph_.SetBaseval(1);
        this->graph_.SetSorted(1);
        this->graph_.ExpandSymmetric();
        logfileptr->OFS()<<"Matrix structure expanded"<<std::endl;

        SparseMatrixGraph * sgraph = nullptr;
        {
          {
            double timeSta = get_time();
            scope_timer(c,symPACKMatrix2D::ordering);

            if (this->options_.orderingStr=="MMD") {
              this->options_.ordering = symPACK::MMD;
            }
            else if (this->options_.orderingStr=="RCM") {
              this->options_.ordering = symPACK::RCM;
            }
            else if (this->options_.orderingStr=="AMD") {
              this->options_.ordering = symPACK::AMD;
            }
#ifdef USE_METIS
            else if (this->options_.orderingStr=="METIS") {
              this->options_.ordering = symPACK::METIS;
            }
#endif
#ifdef USE_SCOTCH
            else if (this->options_.orderingStr=="SCOTCH") {
              this->options_.ordering = symPACK::SCOTCH;
            }
#endif
#ifdef USE_PARMETIS
            else if (this->options_.orderingStr=="PARMETIS") {
              this->options_.ordering = symPACK::PARMETIS;
            }
#endif
#ifdef USE_PTSCOTCH
            else if (this->options_.orderingStr=="PTSCOTCH") {
              this->options_.ordering = symPACK::PTSCOTCH;
            }
#endif
            else if (this->options_.orderingStr=="NATURAL") {
              this->options_.ordering = symPACK::NATURAL;
            }
            else if (this->options_.orderingStr=="USER") {
              if (this->options_.perm ==nullptr) {
                throw std::logic_error( "When using USER, symPACKOptions.perm must be provided.\n" );
              }
              else {
                this->options_.ordering = symPACK::USER;
              }
            }
            else {
              std::stringstream sstr;
              sstr<<"This ordering method is not supported by symPACK. Valid options are:";
              sstr<<"NATURAL MMD AMD RCM NDBOX NDGRID USER ";
#ifdef USE_SCOTCH
              sstr<<"SCOTCH ";
#endif
#ifdef USE_PTSCOTCH
              sstr<<"PTSCOTCH ";
#endif
#ifdef USE_METIS
              sstr<<"METIS ";
#endif
#ifdef USE_PARMETIS
              sstr<<"PARMETIS ";
#endif
              sstr<<std::endl;
              throw std::logic_error( sstr.str()  );

            }

            this->Order_.NpOrdering = this->options_.NpOrdering;
            switch(this->options_.ordering) {
              case MMD:
                {
                  sgraph = new SparseMatrixGraph();
                  this->graph_.GatherStructure(*sgraph,0);
                  sgraph->SetBaseval(1);
                  sgraph->SetKeepDiag(0);
                  this->Order_.MMD(*sgraph, this->graph_.comm);
                }
                break;

              case RCM:
                {
                  sgraph = new SparseMatrixGraph();
                  this->graph_.GatherStructure(*sgraph,0);
                  sgraph->SetBaseval(1);
                  sgraph->SetKeepDiag(0);
                  this->Order_.RCM(*sgraph, this->graph_.comm);
                }
                break;

              case AMD:
                {
                  sgraph = new SparseMatrixGraph();
                  this->graph_.GatherStructure(*sgraph,0);
                  sgraph->SetBaseval(1);
                  sgraph->SetKeepDiag(0);
                  this->Order_.AMD(*sgraph, this->graph_.comm);
                }
                break;

              case USER:
                {
                  this->Order_.perm.resize(this->iSize_);
                  //find baseval
                  auto baseval = std::min_element(&this->options_.perm[0],&this->options_.perm[0]+this->iSize_);
                  //now compute everything in 1 based 
                  for (int i=0;i<this->Order_.perm.size();i++) {this->Order_.perm[i] = this->options_.perm[i] - *baseval + 1; /*1 based*/}
                  this->Order_.invp.resize(this->iSize_);
                  for (int i=0;i<this->Order_.perm.size();i++) {this->Order_.invp[this->Order_.perm[i]-1] = i;} 
                }
                break;

              case NDBOX:
                this->Order_.NDBOX(this->Size(), this->graph_.comm);
                break;

              case NDGRID:
                this->Order_.NDGRID(this->Size(), this->graph_.comm);
                break;

#ifdef USE_SCOTCH
              case SCOTCH:
                {

                  sgraph = new SparseMatrixGraph();
                  this->graph_.GatherStructure(*sgraph,0);
                  sgraph->SetKeepDiag(0);
                  this->Order_.SCOTCH(*sgraph, this->graph_.comm);
                }
                break;
#endif
#ifdef USE_METIS
              case METIS:
                {
                  sgraph = new SparseMatrixGraph();
                  this->graph_.GatherStructure(*sgraph,0);
                  sgraph->SetKeepDiag(0);
                  this->Order_.METIS(*sgraph, this->graph_.comm);
                }
                break;
#endif
#ifdef USE_PARMETIS
              case PARMETIS:
                {
                  this->graph_.SetKeepDiag(0);
                  this->Order_.PARMETIS(this->graph_);
                }
                break;
#endif
#ifdef USE_PTSCOTCH
              case PTSCOTCH:
                {
                  this->graph_.SetKeepDiag(0);
                  this->Order_.PTSCOTCH(this->graph_);
                }
                break;
#endif
              case NATURAL:
                this->Order_.perm.resize(this->iSize_);
                for (int i=0;i<this->Order_.perm.size();i++) {this->Order_.perm[i]=i+1;} 
                this->Order_.invp = this->Order_.perm;
                break;
              default:
                {
                  std::stringstream sstr;
                  sstr<<"This ordering method is not supported by symPACK. Valid options are:";
                  sstr<<"MMD AMD RCM NDBOX NDGRID USER ";
#ifdef USE_SCOTCH
                  sstr<<"SCOTCH ";
#endif
#ifdef USE_PTSCOTCH
                  sstr<<"PTSCOTCH ";
#endif
#ifdef USE_METIS
                  sstr<<"METIS ";
#endif
#ifdef USE_PARMETIS
                  sstr<<"PARMETIS ";
#endif
                  sstr<<std::endl;
                  throw std::logic_error( sstr.str()  );
                  //do nothing: either natural or user provided ordering
                }
                break;
            }

            //The ordering is available on every processor of the full communicator
            double timeStop = get_time();
            if (this->iam==0 && this->options_.verbose) {
              symPACKOS<<"Ordering time: "<<timeStop - timeSta<<std::endl;
            }
            logfileptr->OFS()<<"Ordering done"<<std::endl;
          }
        }

        if (this->options_.dumpPerm>0) {
          logfileptr->OFS()<<"perm = [";
          for (auto i : this->Order_.perm) { 
            logfileptr->OFS()<<i<<" ";
          }
          logfileptr->OFS()<<"]"<<std::endl;
        }

        std::vector<Int> rc;


        if (sgraph==nullptr) { 
          logfileptr->OFS()<<"copying graph"<<std::endl;
          auto graph = this->graph_;
          double timeSta = get_time();
          graph.SetKeepDiag(1);
          graph.SetBaseval(1);
          logfileptr->OFS()<<"graph copied"<<std::endl;
          //Expand to unsymmetric storage
          graph.ExpandSymmetric();
          logfileptr->OFS()<<"graph expanded"<<std::endl;

          auto backinvp = this->Order_.invp;
          graph.Permute(this->Order_.invp.data());
          logfileptr->OFS()<<"graph permuted 1/2"<<std::endl;
          this->ETree_.ConstructETree(graph,this->Order_);
          this->ETree_.PostOrderTree(this->Order_);
          logfileptr->OFS()<<"ETree computed and postordered"<<std::endl;

          std::vector<Int> relinvp;
          this->Order_.GetRelativeInvp(backinvp,relinvp);
          graph.Permute(relinvp.data());
          logfileptr->OFS()<<"graph permuted 2/2"<<std::endl;

          double timeSta_cc = get_time();
          this->getLColRowCount(graph,cc,rc);
          double timeStop_cc = get_time();
          if (this->iam==0 && this->options_.verbose) {
            symPACKOS<<"Column count (distributed) construction time: "<<timeStop_cc - timeSta_cc<<std::endl;
          }

          if (this->options_.ordering != NATURAL) {
            this->ETree_.SortChildren(cc,this->Order_);
            if (this->options_.dumpPerm>0) {
              logfileptr->OFS()<<"perm = "<<this->Order_.perm<<std::endl;
            }
          }

          double timeStop = get_time();
          if (this->iam==0 && this->options_.verbose) {
            symPACKOS<<"Elimination tree construction time: "<<timeStop - timeSta<<std::endl;
          }
        }
        else
        {
          //gather the graph if necessary to build the elimination tree
          if (sgraph==nullptr) { 
            sgraph = new SparseMatrixGraph();
            this->graph_.GatherStructure(*sgraph,0);
          }
          this->graph_.SetKeepDiag(1);
          sgraph->SetBaseval(1);
          sgraph->SetKeepDiag(1);

          double timeSta = get_time();
          this->ETree_.ConstructETree(*sgraph,this->Order_,this->fullcomm_);
          this->ETree_.PostOrderTree(this->Order_);

          double timeSta_cc = get_time();
          this->getLColRowCount(*sgraph,cc,rc);
          double timeStop_cc = get_time();
          if (this->iam==0 && this->options_.verbose) {
            symPACKOS<<"Column count (gather + serial + bcast) construction time: "<<timeStop_cc - timeSta_cc<<std::endl;
          }

          if (this->options_.ordering != NATURAL) {
            this->ETree_.SortChildren(cc,this->Order_);
            if (this->options_.dumpPerm>0) {
              logfileptr->OFS()<<"perm = "<<this->Order_.perm<<std::endl;
            }
          }

          double timeStop = get_time();
          if (this->iam==0 && this->options_.verbose) {
            symPACKOS<<"Elimination tree construction time: "<<timeStop - timeSta<<std::endl;
          }
        }

        //compute some statistics
        if (this->iam==0 && this->options_.verbose) {
          double flops = 0.0;
          int64_t NNZ = 0;
          for (Int i = 0; i<cc.size();++i) {
            flops+= (double)pow((double)cc[i],2.0);
            NNZ+=cc[i];
          }
          symPACKOS<<"Flops: "<<flops<<std::endl;
          symPACKOS<<"NNZ in L factor: "<<NNZ<<std::endl;
        }

        //get rid of the sequential graph
        if (sgraph!=nullptr && this->options_.order_refinement_str.substr(0,4) != "TSPB") {delete sgraph;}

        { 
          double timeSta = get_time();
          this->findSupernodes(this->ETree_,this->Order_,cc,this->SupMembership_,this->Xsuper_,this->options_.relax.maxSize);
          logfileptr->OFS()<<"Supernodes found"<<std::endl;

          if (this->options_.relax.nrelax0>0)
          {
            this->relaxSupernodes(this->ETree_, cc,this->SupMembership_, this->Xsuper_, this->options_.relax );
            logfileptr->OFS()<<"Relaxation done"<<std::endl;
          }

          //modify this->np since it cannot be greater than the number of supernodes
          this->dist_np = std::min(this->all_np,this->options_.used_procs(this->Xsuper_.size()-1)); 
          this->np = this->all_np;

          if (this->workcomm_!=MPI_COMM_NULL) {
            MPI_Comm_free(&this->workcomm_);
          }
          MPI_Comm_split(this->options_.MPIcomm,this->iam<this->np,this->iam,&this->workcomm_);
          this->group_.reset( new RankGroup( this->workcomm_ ) ); 

          //do another split to contain P0 and all the non working processors
          if (this->all_np!=this->np) {
            if (this->non_workcomm_!=MPI_COMM_NULL) {
              MPI_Comm_free(&this->non_workcomm_);
            }
            this->non_workcomm_ = MPI_COMM_NULL;
            MPI_Comm_split(this->options_.MPIcomm,this->iam==0||this->iam>=this->np,this->iam==0?0:this->iam-this->np+1,&this->non_workcomm_);
          }

          //Compute this->XsuperDist_
          std::vector<Idx> newVertexDist;
          {
            Idx supPerProc = std::max((size_t)1,(this->Xsuper_.size()-1) / this->dist_np);
            this->XsuperDist_.resize(this->all_np+1,0);
            this->XsuperDist_[this->dist_np] = this->Xsuper_.size();
            for (int p = 0; p<this->dist_np; p++) {
              this->XsuperDist_[p]= std::min(this->XsuperDist_[this->dist_np],(Int)(p*supPerProc+1));
            }
            for (int p = this->dist_np+1; p<this->all_np; p++) {
              this->XsuperDist_[p]= this->XsuperDist_[p-1];
            }
            this->XsuperDist_[this->all_np] = this->Xsuper_.size();


            newVertexDist.resize(this->all_np+1,0);
            newVertexDist[this->all_np] = this->iSize_+1;
            for (int p = 0; p < this->all_np; p++) {
              Int S = this->XsuperDist_[p];
              newVertexDist[p] = this->Xsuper_[S-1];
            }
          }

          //TODO TEMPORARY UGLY FIX
          this->np = this->dist_np;


          double timeStaSymb = get_time();
          this->symbolicFactorizationRelaxedDist(cc);

          double timeStopSymb = get_time();
          if (this->iam==0 && this->options_.verbose) {
            symPACKOS<<"Symbolic factorization time: "<<timeStopSymb - timeStaSymb<<std::endl;
          }
          logfileptr->OFS()<<"Symbfact done"<<std::endl;


          if (this->options_.order_refinement_str != "NO") {
            double timeSta = get_time();
            if (this->options_.order_refinement_str == "SET") { 
              this->refineSupernodes(3,1,&pMat);
            }
            else if (this->options_.order_refinement_str == "SET10") { 
              this->refineSupernodes(1,0,&pMat);
            }
            else if (this->options_.order_refinement_str == "SET11") { 
              this->refineSupernodes(1,1,&pMat);
            }
            else if (this->options_.order_refinement_str == "SET20") { 
              this->refineSupernodes(2,0,&pMat);
            }
            else if (this->options_.order_refinement_str == "SET21") { 
              this->refineSupernodes(2,1,&pMat);
            }
            else if (this->options_.order_refinement_str == "SET30") { 
              this->refineSupernodes(3,0,&pMat);
            }
            else if (this->options_.order_refinement_str == "SET31") { 
              this->refineSupernodes(3,1,&pMat);
            }
            else if (this->options_.order_refinement_str.substr(0,4) == "TSPB") {

              if (sgraph==nullptr) { 
                sgraph = new SparseMatrixGraph();
                this->graph_.GatherStructure(*sgraph,0);
              }

              ETree& tree = this->ETree_;
              Ordering & aOrder = this->Order_;
              std::vector<Int> & supMembership = this->SupMembership_; 
              std::vector<Int> & xsuper = this->Xsuper_; 

              std::vector<int> ixlindx;
              std::vector<int> ilindx;

              std::vector<int>  new_invp;

              //Gather this->locXlindx_ and this->locLindx_
              {
                std::vector<Ptr> xlindx;
                std::vector<Idx> lindx;
                this->gatherLStructure(xlindx, lindx);

                if (this->iam==0) {
                  ixlindx.resize(xlindx.size());
                  for (int i = 0;i<xlindx.size();i++) {
                    ixlindx[i] = xlindx[i];
                  }
                  ilindx.resize(lindx.size());
                  for (int i = 0;i<lindx.size();i++) {
                    ilindx[i] = lindx[i];
                  }
                }
              }

              if (this->iam==0) {
                bassert(sgraph!=nullptr);
                int neqns = this->iSize_;

                int nofsub =ilindx.size();
                nsuper = xsuper.size()-1;


                new_invp.assign(neqns,0);


                std::vector<int>  new_perm(neqns,0);

                for (size_t i =0; i<new_perm.size(); i++) { new_invp[i] = this->Order_.invp[i];}
                for (size_t i =0; i<new_perm.size(); i++) { new_perm[i] = this->Order_.perm[i];}

                int supsiz = 0;
                for (I = 1; I <= nsuper; I++) {
                  Int fc = this->Xsuper_[I-1];
                  Int lc = this->Xsuper_[I]-1;
                  supsiz = std::max(supsiz,lc-fc+1);
                }

                sgraph->SetKeepDiag(0);
                xadj.resize(sgraph->colptr.size());
                adj.resize(sgraph->rowind.size());
                int nadj = adj.size();
                for (size_t i = 0; i< sgraph->colptr.size(); i++) { xadj[i] = int(sgraph->colptr[i]); }
                for (size_t i = 0; i< sgraph->rowind.size(); i++) { adj[i] = int(sgraph->rowind[i]); }

                if (sgraph!=nullptr) {delete sgraph;}

                std::vector<int, Mallocator<int> > etpar(neqns);
                for (int i = 0; i<neqns; i++) { etpar[i] = tree.PostParent(i); }


                std::vector<int, Mallocator<int> > xskadj(neqns+1);
                std::vector<int, Mallocator<int> > sklenf(neqns);
                std::vector<int, Mallocator<int> > sklenb(neqns);
                std::vector<int, Mallocator<int> > skadj(nadj);
                std::vector<int, Mallocator<int> > invp2(neqns);
                std::vector<int, Mallocator<int> > link(neqns);
                std::vector<int, Mallocator<int> > fstloc(nsuper);
                std::vector<int, Mallocator<int> > sperm(nsuper);
                std::vector<int, Mallocator<int> > fstloc2(neqns);
                std::vector<int, Mallocator<int> > dist1((supsiz+1)*(supsiz+1));
                std::vector<int, Mallocator<int> > suppar(nsuper);
                int iwsiz = 8*neqns+3;
                std::vector<int, Mallocator<int> > iwork(iwsiz);

                int iflag = 0;
                double timeSta = get_time();
                if (this->options_.order_refinement_str == "TSPB") {
                  std::vector<int, Mallocator<int> > rep(neqns);
                  FORTRAN(ordsup_ind_tsp_paths2)
                    ( 
                     &nadj, &neqns , &nofsub, &nsuper, &supsiz, 
                     xsuper.data(), ixlindx.data(), ilindx.data(), supMembership.data(), 
                     xadj.data(), adj.data(), 
                     etpar.data(),
                     new_perm.data(), new_invp.data(), 
                     &iflag , 
                     xskadj.data(), sklenf.data(), sklenb.data(), 
                     skadj.data() , invp2.data() , link.data()  , 
                     fstloc.data(), sperm.data() , fstloc2.data(), dist1.data(), 
                     suppar.data(), &iwsiz , iwork.data() , rep.data()             );
                }
                else {
                  FORTRAN(ordsup_ind_tsp_paths)
                    ( 
                     &nadj, &neqns , &nofsub, &nsuper, &supsiz, 
                     xsuper.data(), ixlindx.data(), ilindx.data(), supMembership.data(), 
                     xadj.data(), adj.data(), 
                     etpar.data(),
                     new_perm.data(), new_invp.data(), 
                     &iflag , 
                     xskadj.data(), sklenf.data(), sklenb.data(), 
                     skadj.data() , invp2.data() , link.data()  , 
                     fstloc.data(), sperm.data() , fstloc2.data(), dist1.data(), 
                     suppar.data(), &iwsiz , iwork.data() );
                }

                double timeStop = get_time();
                if (this->iam==0 && this->options_.verbose) {
                  symPACKOS<<"TSPB reordering done in "<<timeStop-timeSta<<std::endl;
                }


                //this->Order_.Compose(new_invp);
                for (size_t i =0; i<new_perm.size(); i++) { this->Order_.invp[i] = new_invp[i];}
                for (size_t i =0; i<new_perm.size(); i++) { this->Order_.perm[i] = new_perm[i];}
              }

              // broadcast invp
              Int N = aOrder.invp.size();
              MPI_Bcast(&aOrder.invp[0],N*sizeof(Int),MPI_BYTE,0,this->fullcomm_);
              MPI_Bcast(&aOrder.perm[0],N*sizeof(Int),MPI_BYTE,0,this->fullcomm_);

            }

            double timeStop = get_time();

            if (this->iam==0 && this->options_.verbose) {
              symPACKOS<<"Supernode reordering done in "<<timeStop-timeSta<<std::endl;
            }

            {
              double timeSta = get_time();
              this->symbolicFactorizationRelaxedDist(cc);
              double timeStop = get_time();
              if (this->iam==0 && this->options_.verbose) {
                symPACKOS<<"Symbolic factorization time: "<<timeStop - timeSta<<std::endl;
              }
              logfileptr->OFS()<<"Symbfact done"<<std::endl;
            }
          }

          double timeStop = get_time();
          if (this->iam==0 && this->options_.verbose) {
            symPACKOS<<"Total symbolic factorization time: "<<timeStop - timeSta<<std::endl;
          }

          //Print statistics
          if (this->options_.print_stats) {
            OrderStats stats;
            stats.get(this->Xsuper_, this->XsuperDist_, this->locXlindx_, this->locLindx_, this->fullcomm_);
            if (this->iam==0) {
              stats.print();
            }
          }
        }

#ifdef _DEBUG_
        logfileptr->OFS()<<"Membership list is "<<this->SupMembership_<<std::endl;
        logfileptr->OFS()<<"xsuper "<<this->Xsuper_<<std::endl;
#endif


#ifdef _OUTPUT_ETREE_
        logfileptr->OFS()<<"ETree is "<<this->ETree_<<std::endl;
        {
          auto supETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);
          logfileptr->OFS()<<"Supernodal ETree is "<<supETree<<std::endl;
          logfileptr->OFS()<<"Xsuper is "<<this->Xsuper_<<std::endl;
        } 
#endif
      }

      //TODO TEMPORARY UGLY FIX
      this->np = this->all_np;

      std::vector< std::tuple<Int,upcxx::global_ptr<char>> > local_diag_pointers;
      {
        //compute a mapping
        // Do a 2D cell-cyclic distribution for now.
        nsuper = this->Xsuper_.size()-1;
        auto iam = this->iam;
        auto np = this->np;
        Int npcol = std::floor(std::sqrt(np));
        Int nprow = std::floor(np / npcol);

        {
          using cell_tuple_t = std::tuple<Idx,Idx,Int,Int>;
          std::list< cell_tuple_t> cellsToSend;

          Int numLocSnode = this->XsuperDist_[this->iam+1]-this->XsuperDist_[this->iam];
          Int firstSnode = this->XsuperDist_[this->iam];

          //now create cell structures
          for (Int locsupno = this->locXlindx_.size()-1; locsupno >= 1; --locsupno)
          {
            Idx I = locsupno + firstSnode-1;
            Int first_col = this->Xsuper_[I-1];
            Int last_col = this->Xsuper_[I]-1;

            Ptr lfi = this->locXlindx_[locsupno-1];
            Ptr lli = this->locXlindx_[locsupno]-1;

            Int iWidth = last_col - first_col + 1;

            Int * pNrows = nullptr;
            Int * pNblock = nullptr;
#ifdef _SUBCUBE2D_
            double * pCost = nullptr;
#endif

            Idx iStartRow = this->locLindx_[lfi-1];
            Idx iPrevRow = iStartRow;
            Int iContiguousRows = 0;


            Int K = -1;
            Int K_prevSnode = -1;
            for (Ptr K_sidx = lfi; K_sidx<=lli;K_sidx++) {
              Idx K_row = this->locLindx_[K_sidx-1];
              K = this->SupMembership_[K_row-1];

              //Split at boundary or after diagonal block
              if (K!=K_prevSnode) {
                if (K>=I) {
                  //add last block from previous cell
                  if (pNblock) {
                    bassert(K>I);
                    bassert(iContiguousRows>0);
                    //now we are at the end of a block
                    (*pNblock)++;
                    //add a fake cell to describe the block
                    cellsToSend.push_back(std::make_tuple(-1,-1,iStartRow,iContiguousRows));
                    iContiguousRows = 0;
                  }

                  iContiguousRows = 0;
                  iStartRow = K_row;
                  iPrevRow = iStartRow;
                  cellsToSend.push_back(std::make_tuple(K,I,0,0));
                  pNblock = &std::get<2>(cellsToSend.back());
                  pNrows = &std::get<3>(cellsToSend.back());
                }
              }
              K_prevSnode = K;


              //counting contiguous rows
              if (K_row==iPrevRow+1 || K_row==iStartRow ) {
                ++iContiguousRows;
              }
              //add previous block from current cell
              else {
                bassert(iContiguousRows>0);
                //now we are at the end of a block
                (*pNblock)++;
                //add a fake cell to describe the block
                cellsToSend.push_back(std::make_tuple(-1,-1,iStartRow,iContiguousRows));
                iContiguousRows = 1;
                iStartRow = K_row;
              }
              iPrevRow=K_row;
              (*pNrows)++;
            }

            if (pNblock) {
              //now we are at the end of a block
              (*pNblock)++;
              //add a fake cell to describe the block
              bassert(iContiguousRows>0);
              cellsToSend.push_back(std::make_tuple(-1,-1,iStartRow,iContiguousRows));
              iContiguousRows = 0;
            }

          }



          MPI_Datatype type;
          MPI_Type_contiguous( sizeof(cell_tuple_t), MPI_BYTE, &type );
          MPI_Type_commit(&type);

          //then do an allgatherv
          //compute send sizes
          int ssize = cellsToSend.size();

          //Build the contiguous array
          vector<cell_tuple_t> sendbuf;
          sendbuf.reserve(ssize);
          sendbuf.insert(sendbuf.end(),cellsToSend.begin(),cellsToSend.end());

          //allgather receive sizes
          vector<int> rsizes(this->np,0);
          MPI_Allgather(&ssize,sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,this->workcomm_);

          //compute receive displacements
          vector<int> rdispls(this->np+1,0);
          rdispls[0] = 0;
          std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);


          //Now do the allgatherv
          vector<cell_tuple_t> recvbuf(rdispls.back());
          MPI_Allgatherv(&sendbuf[0],ssize,type,&recvbuf[0],&rsizes[0],&rdispls[0],type,this->workcomm_);
          MPI_Type_free(&type);

          std::map<std::tuple<int,int>, int> mapping;
          using ProcGroup = TreeLoadBalancer::ProcGroup;
          std::vector< ProcGroup > levelGroups_;
          std::vector< Int > groupIdx_;
          std::vector<double> NodeLoad;
          std::vector< std::map<Int,double> > LocalCellLoad(this->locXlindx_.size());

          std::map< int_t, std::map<int_t,double> > CellLoad;
          {

            auto factor_cost = [](Int m, Int n)->double{
              return FLOPS_DPOTRF(m);
            };

            auto trsm_cost = [](Int m, Int n)->double{
              return FLOPS_DTRSM('L',m,n);
            };

            auto update_cost = [](Int m, Int n, Int k)->double{
              return FLOPS_DGEMM(m,n,k);
            };


            auto supETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);
            Int numLevels = 1;

            std::vector<double> SubTreeLoad(supETree.Size()+1,0.0);
            std::vector<int> children(supETree.Size()+1,0);
            Int numLocSnode = this->XsuperDist_[iam+1]-this->XsuperDist_[iam];
            Int firstSnode = this->XsuperDist_[iam];

            for (Int locsupno = 1; locsupno<this->locXlindx_.size(); ++locsupno)
            {
              Int I = locsupno + firstSnode-1;

              Idx fc = this->Xsuper_[I-1];
              Idx lc = this->Xsuper_[I]-1;

              Int width = lc - fc + 1;
              Int height = cc[fc-1];

              //cost of factoring curent panel
              LocalCellLoad[locsupno-1][I]+=factor_cost(width,width);
              //cost of updating ancestor and doing TRSMs)s
              {
                Ptr lfi = this->locXlindx_[locsupno-1];
                Ptr lli = this->locXlindx_[locsupno]-1;

                auto getPerCellCost = [this,&update_cost,&LocalCellLoad,&CellLoad,width,I,firstSnode,numLocSnode,&trsm_cost](Int tgt_snode_id, Ptr rowptr, Ptr lli) {
                  Int K = -1;
                  Int K_prevSnode = -1;
                  Int nrows = 0;
                  Int odRows = 0;
                  Int locsupno = I - firstSnode+1;
                  for (Ptr K_sidx = rowptr; K_sidx<=lli;K_sidx++) {
                    Idx K_row = this->locLindx_[K_sidx-1];
                    K = this->SupMembership_[K_row-1];
                    //Split at boundary or after diagonal block
                    if (K!=K_prevSnode) {
                      if (K_prevSnode==tgt_snode_id) {
                        odRows = nrows;
                      }

                      if (K_prevSnode>=tgt_snode_id) {
                        Int m = nrows;
                        Int n = odRows;
                        Int k = width;
                        double cost = update_cost(m,n,k);
                        if (tgt_snode_id >= firstSnode+numLocSnode) {
                          CellLoad[tgt_snode_id][K_prevSnode]+=cost;
                        }
                        else {
                          LocalCellLoad[tgt_snode_id-firstSnode][K_prevSnode]+=cost;
                        }
                      }

                      if (K_prevSnode>I) {
                        Int m = nrows;
                        Int n = odRows;
                        Int k = width;
                        //cost of a TRSM = one fourth of gemm
                        double cost = trsm_cost(nrows,width);
                        LocalCellLoad[locsupno-1][K_prevSnode]+=cost;
                        bassert(K_prevSnode>0);
                      }

                      if (K>=tgt_snode_id) {
                        nrows = 0;
                      }
                    }
                    nrows++;
                    K_prevSnode = K;
                  }

                  if (K_prevSnode>=tgt_snode_id) {
                    Int m = nrows;
                    Int n = odRows;
                    Int k = width;
                    double cost = update_cost(m,n,k);
                    if (tgt_snode_id >= firstSnode+numLocSnode) {
                      CellLoad[tgt_snode_id][K_prevSnode]+=cost;
                    }
                    else {
                      LocalCellLoad[tgt_snode_id-firstSnode][K_prevSnode]+=cost;
                    }
                    bassert(K_prevSnode>0);
                  }

                  if (K_prevSnode>I) {
                    Int m = nrows;
                    Int n = odRows;
                    Int k = width;
                    //cost of a TRSM = one fourth of gemm
                    double cost = trsm_cost(nrows,width);
                    LocalCellLoad[locsupno-1][K_prevSnode]+=cost;
                    bassert(K_prevSnode>0);
                  }
                };

                //parse rows
                Int tgt_snode_id = I;
                Idx tgt_fr = fc;
                Idx tgt_lr = fc;
                Ptr tgt_fr_ptr = 0;
                Ptr tgt_lr_ptr = 0;

                getPerCellCost(tgt_snode_id,lfi,lli);

                for (Ptr rowptr = lfi; rowptr<=lli;rowptr++) {
                  Idx row = this->locLindx_[rowptr-1]; 
                  if (this->SupMembership_[row-1]==tgt_snode_id) { continue;}

                  //we have a new tgt_snode_id
                  tgt_snode_id = this->SupMembership_[row-1];

                  getPerCellCost(tgt_snode_id,rowptr,lli);
                }
              }
            }

            //LocalCellLoad has to be Allgatherv'd
            {
              using cell_load_t = std::tuple<int_t,double>; 
              std::vector< cell_load_t > load_sendbuf,load_recvbuf;
              int load_send_size = 0;
              for (auto & m : LocalCellLoad) {
                load_send_size += 1 + m.size();
              }
              int remote_load_send_size = 0;
              for (auto & m : CellLoad) {
                remote_load_send_size += 1 + m.second.size();
              }

              load_send_size+=remote_load_send_size;
              load_sendbuf.reserve(load_send_size);
              //now serialize
              int supno = firstSnode;
              for (auto & m : LocalCellLoad) {
                load_sendbuf.push_back(std::make_tuple(-supno,0.0));
                for (auto & c: m) {
                  load_sendbuf.push_back(std::make_tuple(c.first,c.second));
                }
                supno++;
              }


              for (auto & m : CellLoad) {
                load_sendbuf.push_back(std::make_tuple(-m.first,0.0));
                for (auto & c: m.second) {
                  load_sendbuf.push_back(std::make_tuple(c.first,c.second));
                }
              }


              MPI_Datatype cell_load_type;
              MPI_Type_contiguous( sizeof(cell_load_t), MPI_BYTE, &cell_load_type );
              MPI_Type_commit(&cell_load_type);
              std::vector<int> load_rsizes(this->np,0);
              std::vector<int> remote_load_rsizes(this->np,0);
              std::vector<int> load_rdispls(this->np+1,0);

              MPI_Allgather(&load_send_size,sizeof(load_send_size),MPI_BYTE,load_rsizes.data(),sizeof(load_send_size),MPI_BYTE,this->workcomm_);
              MPI_Allgather(&remote_load_send_size,sizeof(remote_load_send_size),MPI_BYTE,remote_load_rsizes.data(),sizeof(remote_load_send_size),MPI_BYTE,this->workcomm_);

              load_rdispls[0] = 0;
              std::partial_sum(load_rsizes.begin(),load_rsizes.end(),&load_rdispls[1]);
              load_recvbuf.resize(load_rdispls.back());

              //now communicate the task_idx arrays: allgatherv
              MPI_Allgatherv(load_sendbuf.data(),load_send_size,cell_load_type,load_recvbuf.data(),load_rsizes.data(),load_rdispls.data(),cell_load_type,this->workcomm_);

              CellLoad.clear();
              //content of load_recvbuf can now be added to CellLoad: load_rsizes - remote_load_rsizes can directly be copied while the other need to be accumulated
              for ( int p =0; p<np; p++) {
                int cur_tgt = 0;
                for ( int idx = load_rdispls[p]; idx < load_rdispls[p+1]; idx++ ) {
                  auto & current_cell = load_recvbuf[idx];
                  if ( std::get<0>(current_cell) < 0) {
                    cur_tgt = -std::get<0>(current_cell);
                  }
                  else {
                    bassert(cur_tgt!=0);
                    CellLoad[cur_tgt][std::get<0>(current_cell)] += std::get<1>(current_cell);
                  }
                }
              } 



              MPI_Type_free(&cell_load_type);
            }
            NodeLoad.assign(supETree.Size()+1,0.0);
            //now we need to reduce this CellLoad. AlltoAllv then merge/reduce? 
            for (auto & m : CellLoad) {
              NodeLoad[m.first] = 0.0;
              for (auto & c: m.second) {
                NodeLoad[m.first] += c.second;
              }
            }

            SubTreeLoad = NodeLoad;


            for (Int I=1;I<=supETree.Size();I++) {
              Int parent = supETree.Parent(I-1);
              ++children[parent];
              SubTreeLoad[parent]+=SubTreeLoad[I];
            }

#ifndef _NDEBUG_
            auto check = [&NodeLoad,&SubTreeLoad,&supETree]() { double val1 = std::accumulate(NodeLoad.begin(),NodeLoad.end(),0.0); double val2 = 0.0; for (Int I=1;I<=supETree.Size();I++){ Int parent = supETree.Parent(I-1); if (parent == 0 ) { val2 += SubTreeLoad[I];} } ; return val1==val2;};
            bassert( check() );
#endif

            using node_t = struct _node_t { int id; std::vector<_node_t *> children; };
            std::vector<node_t> nodes(nsuper+1);
            node_t * root = &nodes[0];
            for (Int I=nsuper; I>= 1;I--) { 
              Int parent = supETree.Parent(I-1);
              nodes[I].id = I;
              nodes[parent].children.push_back(&nodes[I]);
            }

            int n_ = nsuper;
            std::vector<Int> levels(supETree.Size()+1,0);
            for (Int I=n_; I>= 1;I--) { 
              Int parent = supETree.Parent(I-1);
              if (parent==0) {levels[I]=0;}
              else { levels[I] = levels[parent]+1; numLevels = std::max(numLevels,levels[I]);}
            }
            numLevels++;

            std::vector< ProcGroup > procGroups_;
            groupIdx_.resize(n_+1,0);
            procGroups_.resize(n_+1);
            procGroups_[0].Ranks().reserve(np);
            for (Int p = 0;p<np;++p) { procGroups_[0].Ranks().push_back(p);}

            std::vector<Int> pstart(n_+1,0);
            levelGroups_.reserve(numLevels);
            levelGroups_.push_back(ProcGroup());
            levelGroups_[0].Ranks().reserve(np);
            for (Int p = 0;p<np;++p) {levelGroups_[0].Ranks().push_back(p);}


            {
              auto & pstart2 = pstart;
              auto & procGroups2 = procGroups_;
              auto & groupIdx2 = groupIdx_;
              auto & levelGroups2 = levelGroups_;

#ifndef _NO_SUBCUBE2D_
              std::function<void(node_t *, node_t*)> recSubtree = [&](node_t * node, node_t * parent) {
                Int npParent = levelGroups2[groupIdx2[parent->id]].Ranks().size();
#ifdef _SUBCUBE_SQRT_LIMIT_
                Int procLimit = std::ceil(std::sqrt(np));
#else
                Int procLimit = 1;
#endif
                if ( npParent < procLimit || parent->children.size() == 1) {
                  groupIdx2[node->id] = groupIdx2[parent->id];
                }
                else {
                  Int I = node->id;
                  double parent_load = 0.0;
                  if (parent != root) {
                    Int fc = this->Xsuper_[parent->id-1];
                    Int width = this->Xsuper_[parent->id] - this->Xsuper_[parent->id-1];
                    Int height = cc[fc-1];
                    parent_load = NodeLoad[parent->id];
                  }

                  double proportion = std::min(1.0,SubTreeLoad[I]/(SubTreeLoad[parent->id]-parent_load));
                  Int pFirstIdx = std::min(pstart2[parent->id],npParent-1);
                  Int npIdeal =(Int)std::round(npParent*proportion);
                  Int numProcs = std::max((Int)1,std::min(npParent-pFirstIdx,npIdeal));
                  Int pFirst = levelGroups2[groupIdx2[parent->id]].Ranks().at(pFirstIdx);

                  pstart2[parent->id]+= numProcs;
                  pstart2[parent->id] = pstart2[parent->id] % npParent;

                  if (npParent!=numProcs) {
                    levelGroups2.push_back(ProcGroup());
                    levelGroups2.back().Ranks().reserve(numProcs);

                    std::vector<Int> & parentRanks = levelGroups2[groupIdx2[parent->id]].Ranks();
                    levelGroups2.back().Ranks().insert(levelGroups2.back().Ranks().begin(),parentRanks.begin()+pFirstIdx,parentRanks.begin()+pFirstIdx+numProcs);

                    groupIdx2[I] = levelGroups2.size()-1;
                  }
                  else {
                    groupIdx2[I] = groupIdx2[parent->id];
                  }
                }
                for ( auto & c: node->children) {
                  recSubtree(c,node);
                }
              };

              double timeSta = get_time();
              recSubtree(root,root);


              double timeEnd = get_time();
              logfileptr->OFS()<<"time subtree "<<timeEnd-timeSta<<std::endl;
#endif
            }
          }

          std::vector<double> procLoad(np,0);
          std::vector<std::vector<double>> group_load(nsuper);
          auto find_min_proc = [this,&levelGroups_,&groupIdx_,&group_load,&CellLoad,&procLoad](int i,int j) {
#ifdef _NO_SUBCUBE2D_
            Int nprow = std::floor(std::sqrt(this->np));
            auto pcol = j % nprow;
            auto prow = i % nprow;
            Int minLoadP = pcol*nprow+prow;
            double local_load = CellLoad[j][i];
            procLoad[minLoadP]+=local_load;
            return minLoadP;
#else
            Int minLoadP= -1;
            double minLoad = -1;
            ProcGroup & group = levelGroups_[groupIdx_[j-1]];

            auto & ranks = group.Ranks();
            for (Int i = 0; i<ranks.size();++i) {
              if (procLoad[ranks[i]]<minLoad || minLoad==-1) {
                minLoad = procLoad[ranks[i]];
                minLoadP = i;
              }
            }

            double local_load = CellLoad[j][i];
            procLoad[ranks[minLoadP]]+=local_load;
            return ranks[minLoadP];
#endif
          };


          snodeBlockBase_sptr_t pLast_cell = nullptr;

          for (int psend = np-1; psend>=0; --psend) {
            for (int idxcell = rdispls[psend]; idxcell<rdispls[psend+1]; idxcell++) {
              auto & cur_cell = recvbuf[idxcell];

              Int i = std::get<0>(cur_cell);
              Int j = std::get<1>(cur_cell);

              if (i==-1 && j==-1) {
                bassert(pLast_cell!=nullptr);
                if ( pLast_cell->owner == iam ) {
                  auto ptr = std::static_pointer_cast<snodeBlock_t>(pLast_cell);
                  Int first_row = std::get<2>(cur_cell);
                  Int nrows = std::get<3>(cur_cell);
                  ptr->add_block(first_row,nrows);
                }
              }
              else {
                Int nBlock = std::get<2>(cur_cell);
                Int nRows = std::get<3>(cur_cell);
                auto idx = coord2supidx(i-1,j-1);
                //j is supernode index
                auto p = find_min_proc(i,j);
                Idx fc = this->Xsuper_[j-1];
                Idx lc = this->Xsuper_[j]-1;
                Int iWidth = lc-fc+1;
                size_t block_cnt = nBlock;
                size_t nnz = nRows * iWidth;

                snodeBlockBase_sptr_t sptr(nullptr);
                if ( p == iam ) {
                  if (this->options_.decomposition == DecompositionType::LDL) {
                    sptr = std::static_pointer_cast<blockCellBase_t>(std::make_shared<snodeBlockLDL_t>(i,j,fc,iWidth,nnz,block_cnt));
                    this->localBlocks_.push_back(std::static_pointer_cast<snodeBlock_t>(sptr));
                    if ( i == j ) {
                      auto diagcell_ptr = std::dynamic_pointer_cast<snodeBlockLDL_t>(sptr);
                      assert(diagcell_ptr);
                      local_diag_pointers.push_back(std::make_tuple(i,diagcell_ptr->Diag()));
                    }
                  }
                  else {
                    sptr = std::static_pointer_cast<blockCellBase_t>(std::make_shared<snodeBlock_t>(i,j,fc,iWidth,nnz,block_cnt));
                    this->localBlocks_.push_back(std::static_pointer_cast<snodeBlock_t>(sptr));
                  }
                }
                else {
                  sptr = std::make_shared<blockCellBase_t>();
                }
                cells_[idx] = sptr;
                sptr->i = i;
                sptr->j = j;
                sptr->owner = p;
                pLast_cell = sptr;

#ifndef _NDEBUG_
                if ( p ==iam) {
                  auto & test = CELL(i-1,j-1);
                }
                else { 
                  auto test = pQueryCELL(i-1,j-1);
                  assert(test);
                }
#endif
              }
            }
          }

          std::sort(this->localBlocks_.begin(),this->localBlocks_.end(),[](snodeBlock_sptr_t & a, snodeBlock_sptr_t & b) {
              return a->j < b->j || (a->i < b->i &&  a->j == b->j ) ; 
              });

        }
      }

      logfileptr->OFS()<<"#supernodes = "<<nsuper<<" #cells = "<<cells_.size()<< " which is "<<cells_.size()*sizeof(snodeBlock_t)<<" bytes"<<std::endl;

      if (this->iam==0) {
        symPACKOS<<"#supernodes = "<<nsuper<<" #cells = "<<cells_.size()<< " which is "<<cells_.size()*sizeof(snodeBlock_t)<<" bytes"<<std::endl;
      }

      //generate task graph for the factorization
      {
        MPI_Datatype type;
        MPI_Type_contiguous( sizeof(SparseTask2D::meta_t), MPI_BYTE, &type );
        MPI_Type_commit(&type);
        vector<int> ssizes(this->np,0);
        vector<int> sdispls(this->np+1,0);
        vector<SparseTask2D::meta_t> sendbuf;
        auto supETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);

        vector<int> ssizes_dep(this->np,0);
        vector<int> sdispls_dep(this->np+1,0);
        vector<SparseTask2D::meta_t> sendbuf_dep;
        {  
          //we will need to communicate if only partial xlindx_, lindx_
          //idea: build tasklist per processor and then exchange
          //tuple is: src_snode,tgt_snode,op_type,lt_first_row = updated fc, facing_first_row = updated fr,
          std::map<Idx, std::list< SparseTask2D::meta_t> > Updates;
          std::map<Idx, std::list< SparseTask2D::meta_t> > Messages;
          std::vector<int> marker(this->np,0);

          Int numLocSnode = this->XsuperDist_[this->iam+1]-this->XsuperDist_[this->iam];
          Int firstSnode = this->XsuperDist_[this->iam];
          for (Int locsupno = 1; locsupno<this->locXlindx_.size(); ++locsupno) {
            Idx I = locsupno + firstSnode-1;
            Int first_col = this->Xsuper_[I-1];
            Int last_col = this->Xsuper_[I]-1;

            bassert( cells_.find(coord2supidx(I-1,I-1)) != cells_.end() );

            auto ptr_fcell = pQueryCELL(I-1,I-1);
            Int iOwner = ptr_fcell->owner;

            //Create the factor task on the owner
            Updates[iOwner].push_back(std::make_tuple(I,I,Factorization::op_type::FACTOR,first_col,first_col));

            Ptr lfi = this->locXlindx_[locsupno-1];
            Ptr lli = this->locXlindx_[locsupno]-1;
            std::list<Int> ancestor_rows;

            Int K = -1; 
            Idx K_prevSnode = (Idx)-1;
            for (Ptr K_sidx = lfi; K_sidx<=lli;K_sidx++) {
              Idx K_row = this->locLindx_[K_sidx-1];
              K = this->SupMembership_[K_row-1];

              //Split at boundary or after diagonal block
              if (K!=K_prevSnode) {
                if (K>I) {
                  ancestor_rows.push_back(K_row);
                }
              }
              K_prevSnode = K;
            }

            for (auto J_row : ancestor_rows) {
              Int J = this->SupMembership_[J_row-1];
              auto ptr_fodcell = pQueryCELL(J-1,I-1);
              Int iFODOwner = ptr_fodcell->owner;

              Messages[iOwner].push_back(std::make_tuple(I,I,Factorization::op_type::TRSM_SEND,J_row,J_row));
              Messages[iFODOwner].push_back(std::make_tuple(I,I,Factorization::op_type::TRSM_RECV,J_row,J_row));
              Updates[iFODOwner].push_back(std::make_tuple(I,I,Factorization::op_type::TRSM,J_row,J_row));

              //add the UPDATE_SEND from FODOwner tasks
              for (auto K_row : ancestor_rows) {
                K = this->SupMembership_[K_row-1];

                if (K>=J) {
                  auto ptr_tgtupdcell = pQueryCELL(K-1,J-1);
                  Int iTgtOwner = ptr_tgtupdcell->owner;
                  //TODO update this for non fan-out mapping
                  Int iUpdOwner = ptr_tgtupdcell->owner;

                  //cell(J,I) to cell(K,J)
                  Messages[iFODOwner].push_back(std::make_tuple(I,J,Factorization::op_type::UPDATE2D_SEND_OD,K_row,J_row));
                  Messages[iUpdOwner].push_back(std::make_tuple(I,J,Factorization::op_type::UPDATE2D_RECV_OD,K_row,J_row));
                }

                if (K<=J) {
                  auto ptr_tgtupdcell = pQueryCELL(J-1,K-1);
                  Int iTgtOwner = ptr_tgtupdcell->owner;
                  //TODO update this for non fan-out mapping
                  Int iUpdOwner = ptr_tgtupdcell->owner;
                  auto ptr_facingcell = pQueryCELL(J-1,I-1);
                  Int iFacingOwner = ptr_facingcell->owner;

                  //sender point of view
                  //cell(J,I) to cell(J,K)
                  Messages[iFacingOwner].push_back(std::make_tuple(I,K,Factorization::op_type::UPDATE2D_SEND,K_row,J_row));
                  Messages[iUpdOwner].push_back(std::make_tuple(I,K,Factorization::op_type::UPDATE2D_RECV,K_row,J_row));

                  if (this->options_.decomposition == DecompositionType::LDL) {
                    Messages[iOwner].push_back(std::make_tuple(I,K,Factorization::op_type::UPDATE2D_DIAG_SEND,K_row,J_row));
                    Messages[iUpdOwner].push_back(std::make_tuple(I,K,Factorization::op_type::UPDATE2D_DIAG_RECV,K_row,J_row));
                  }
                }
              }

              std::set<Int> dlt_sent;
              for (auto K_row : ancestor_rows) {
                K = this->SupMembership_[K_row-1];
                if (K>=J) {
                  auto ptr_tgtupdcell = pQueryCELL(K-1,J-1);
                  Int iTgtOwner = ptr_tgtupdcell->owner;
                  //TODO update this for non fan-out mapping
                  Int iUpdOwner = ptr_tgtupdcell->owner;

                  auto ptr_facingcell = pQueryCELL(J-1,I-1);
                  Int iFacingOwner = ptr_facingcell->owner;

                  //update on cell(K,J)
                  Updates[iUpdOwner].push_back(std::make_tuple(I,J,Factorization::op_type::UPDATE2D_COMP,J_row,K_row));

                  if (this->options_.decomposition == DecompositionType::LDL) {
                    if ( dlt_sent.count(iUpdOwner) == 0 ) {
                      dlt_sent.insert(iUpdOwner);
                      Updates[iUpdOwner].push_back(std::make_tuple(I,J,Factorization::op_type::DLT2D_COMP,J_row,J_row));
                    }
                  }

                  //cell(K,J) to cell (K,J)
                  Messages[iUpdOwner].push_back(std::make_tuple(I,J,Factorization::op_type::AGGREGATE2D_SEND,J_row,K_row));
                  Messages[iTgtOwner].push_back(std::make_tuple(I,J,Factorization::op_type::AGGREGATE2D_RECV,J_row,K_row));
                }
              }
            }
          }

          //then do an alltoallv
          //compute send sizes
          for (auto itp = Updates.begin();itp!=Updates.end();itp++) {
            ssizes[itp->first] = itp->second.size();
          }

          //compute send displacements
          sdispls[0] = 0;
          std::partial_sum(ssizes.begin(),ssizes.end(),&sdispls[1]);

          //Build the contiguous array of pairs
          sendbuf.reserve(sdispls.back());

          for (auto itp = Updates.begin();itp!=Updates.end();itp++) {
            sendbuf.insert(sendbuf.end(),itp->second.begin(),itp->second.end());
          }

          //then do an alltoallv
          //compute send sizes
          for (auto itp = Messages.begin();itp!=Messages.end();itp++) {
            ssizes_dep[itp->first] = itp->second.size();
          }

          //compute send displacements
          sdispls_dep[0] = 0;
          std::partial_sum(ssizes_dep.begin(),ssizes_dep.end(),&sdispls_dep[1]);

          //Build the contiguous array of pairs
          sendbuf_dep.reserve(sdispls_dep.back());

          for (auto itp = Messages.begin();itp!=Messages.end();itp++) {
            sendbuf_dep.insert(sendbuf_dep.end(),itp->second.begin(),itp->second.end());
          }

        }

        //gather receive sizes
        vector<int> rsizes(this->np,0);
        MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,this->workcomm_);

        //compute receive displacements
        vector<int> rdispls(this->np+1,0);
        rdispls[0] = 0;
        std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);

        //Now do the alltoallv
        vector<SparseTask2D::meta_t> recvbuf(rdispls.back());
        MPI_Alltoallv(&sendbuf[0],&ssizes[0],&sdispls[0],type,&recvbuf[0],&rsizes[0],&rdispls[0],type,this->workcomm_);

        //clear the send buffer
        {
          vector<SparseTask2D::meta_t> tmp;
          sendbuf.swap( tmp );
        }


        //gather receive sizes
        vector<int> rsizes_dep(this->np,0);
        MPI_Alltoall(&ssizes_dep[0],sizeof(int),MPI_BYTE,&rsizes_dep[0],sizeof(int),MPI_BYTE,this->workcomm_);

        //compute receive displacements
        vector<int> rdispls_dep(this->np+1,0);
        rdispls_dep[0] = 0;
        std::partial_sum(rsizes_dep.begin(),rsizes_dep.end(),&rdispls_dep[1]);

        //Now do the alltoallv
        vector<SparseTask2D::meta_t> recvbuf_dep(rdispls_dep.back());
        MPI_Alltoallv(&sendbuf_dep[0],&ssizes_dep[0],&sdispls_dep[0],type,&recvbuf_dep[0],&rsizes_dep[0],&rdispls_dep[0],type,this->workcomm_);
        MPI_Type_free(&type);

        //clear the send buffer
        {
          vector<SparseTask2D::meta_t> tmp;
          sendbuf_dep.swap( tmp );
        }

        //do a top down traversal of my local task list and create the computational tasks
        task_graph.clear();
        task_graph.reserve(recvbuf.size());
        task_idx.clear();
        task_idx.reserve(recvbuf.size());

        for (auto it = recvbuf.begin();it!=recvbuf.end();it++) {
          auto & cur_op = (*it);
          auto & src_snode = std::get<0>(cur_op);
          auto & tgt_snode = std::get<1>(cur_op);
          auto & type = std::get<2>(cur_op);
          auto & lt_first_row = std::get<3>(cur_op);
          auto & facing_first_row = std::get<4>(cur_op);

          SparseTask2D * ptask = nullptr;
          switch(type) {
            case Factorization::op_type::TRSM:
            case Factorization::op_type::FACTOR:
            case Factorization::op_type::UPDATE2D_COMP:
            case Factorization::op_type::DLT2D_COMP:
              {
                ptask = new SparseTask2D;
                size_t tuple_size = sizeof(SparseTask2D::meta_t);
                auto meta = &ptask->_meta;
                meta[0] = *it;
              }
              break;
            default:
              break;
          }

          auto I = src_snode;
          auto J = tgt_snode;

          if (ptask!=nullptr) {
            scheduling::key_t key(this->SupMembership_[std::get<4>(cur_op)-1], std::get<1>(cur_op), std::get<0>(cur_op), std::get<2>(cur_op),this->iam);
            task_graph.push_back( std::unique_ptr<SparseTask2D>(ptask) );
            task_idx.push_back( std::make_tuple(key,task_idx.size()));
          }
#ifdef _VERBOSE_
          switch(type) {
            case Factorization::op_type::TRSM:
              {
                auto J = this->SupMembership_[facing_first_row-1];
                logfileptr->OFS()<<"TRSM"<<" from "<<I<<" to ("<<I<<") cell ("<<J<<","<<I<<")"<<std::endl;
              }
              break;
            case Factorization::op_type::FACTOR:
              {
                logfileptr->OFS()<<"FACTOR"<<" cell ("<<J<<","<<I<<")"<<std::endl;
              }
              break;
            case Factorization::op_type::UPDATE2D_COMP:
              {
                auto K = this->SupMembership_[facing_first_row-1];
                logfileptr->OFS()<<"UPDATE"<<" from "<<I<<" to "<<J<<" facing cell ("<<K<<","<<I<<") and cell ("<<J<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;
              }
              break;
            default:
              break;
          }
#endif
        }

        //sort task_idx by keys
        std::sort(task_idx.begin(),task_idx.end(),[](std::tuple<scheduling::key_t, std::size_t > &a, std::tuple<scheduling::key_t, std::size_t >&b) { return std::get<0>(a) < std::get<0>(b);});

        MPI_Datatype task_idx_type;
        MPI_Type_contiguous( sizeof(std::tuple<scheduling::key_t, std::size_t >), MPI_BYTE, &task_idx_type );
        MPI_Type_commit(&task_idx_type);
        std::vector< std::tuple<scheduling::key_t, std::size_t > > task_idx_recvbuf;
        std::vector<int> task_idx_rsizes(this->np,0);
        std::vector<int> task_idx_rdispls(this->np+1,0);

        int task_idx_size = task_idx.size();
        MPI_Allgather(&task_idx_size,sizeof(task_idx_size),MPI_BYTE,task_idx_rsizes.data(),sizeof(task_idx_size),MPI_BYTE,this->workcomm_);

        task_idx_rdispls[0] = 0;
        std::partial_sum(task_idx_rsizes.begin(),task_idx_rsizes.end(),&task_idx_rdispls[1]);
        task_idx_recvbuf.resize(task_idx_rdispls.back());

        //now communicate the task_idx arrays: allgatherv
        MPI_Allgatherv(task_idx.data(),task_idx_size,task_idx_type,task_idx_recvbuf.data(),task_idx_rsizes.data(),task_idx_rdispls.data(),task_idx_type,this->workcomm_);
        MPI_Type_free(&task_idx_type);

        auto key_lb_comp = [](const std::tuple<scheduling::key_t, std::size_t > & a, const scheduling::key_t & key) { return std::get<0>(a) < key;};
        //now we can process the dependency tasks

        for (auto it = recvbuf_dep.begin();it!=recvbuf_dep.end();it++) {
          auto & cur_op = (*it);
          auto & src_snode = std::get<0>(cur_op);
          auto & tgt_snode = std::get<1>(cur_op);
          auto & type = std::get<2>(cur_op);
          auto & lt_first_row = std::get<3>(cur_op);
          auto & facing_first_row = std::get<4>(cur_op);

          auto I = src_snode;
          auto J = tgt_snode;

          switch(type) {
            case Factorization::op_type::UPDATE2D_DIAG_SEND:
              {
                Idx K = this->SupMembership_[facing_first_row-1];
                std::swap(K,J);
                auto k1 = scheduling::key_t(I,I,I,Factorization::op_type::FACTOR,this->iam);
                auto outtask_idx_it = std::lower_bound(task_idx.begin(), task_idx.end(), k1, key_lb_comp);
                auto taskptr = task_graph[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr);
                auto k2 = scheduling::key_t(this->SupMembership_[std::get<4>(taskptr->_meta)-1], std::get<1>(taskptr->_meta), std::get<0>(taskptr->_meta), std::get<2>(taskptr->_meta),this->iam);
                bassert(k2 == k1 );

                auto owner = pQueryCELL(J-1,K-1)->owner;
                auto task_idx_beg = &task_idx_recvbuf[task_idx_rdispls[owner]];
                auto task_idx_end = task_idx_beg + task_idx_rsizes[owner];

                std::size_t remote_task_idx = 0;
                scheduling::key_t key(J,K,I,Factorization::op_type::UPDATE2D_COMP,owner);

                if (this->options_.decomposition == DecompositionType::LDL) {
                  key = scheduling::key_t(K,K,I,Factorization::op_type::DLT2D_COMP,owner);
                }

                auto task_idx_it = std::lower_bound(task_idx_beg, task_idx_end, (key), key_lb_comp);
                bassert(task_idx_it != task_idx_end); 
                remote_task_idx = std::get<1>(*task_idx_it);

                taskptr->out_dependencies[owner].push_back(remote_task_idx);
              }
              break;
            case Factorization::op_type::UPDATE2D_DIAG_RECV:
              {
                Idx K = this->SupMembership_[facing_first_row-1];
                std::swap(K,J);

                if ( this->options_.decomposition == DecompositionType::LDL ) {
                  auto k1 = scheduling::key_t(K,K,I,Factorization::op_type::DLT2D_COMP,this->iam);
                  auto outtask_idx_it = std::lower_bound(task_idx.begin(), task_idx.end(), k1, key_lb_comp);
                  auto taskptr = task_graph[std::get<1>(*outtask_idx_it)].get();
                  bassert(taskptr!=nullptr);
                  bassert(std::get<2>(taskptr->_meta)==Factorization::op_type::DLT2D_COMP);
                  if (  pQueryCELL(I-1,I-1)->owner != pQueryCELL(J-1,K-1)->owner )
                    taskptr->in_remote_dependencies_cnt++;
                  else
                    taskptr->in_local_dependencies_cnt++;

                  //now add a local dependency for the UPDATE2D_COMP task
                  std::size_t remote_task_idx = 0;
                  auto k2 = scheduling::key_t(J,K,I,Factorization::op_type::UPDATE2D_COMP,this->iam);
                  auto intask_idx_it = std::lower_bound(task_idx.begin(), task_idx.end(), k2, key_lb_comp);
                  auto intaskptr = task_graph[std::get<1>(*intask_idx_it)].get();
                  bassert(intaskptr!=nullptr);
                  intaskptr->in_local_dependencies_cnt++;
                  auto owner = pQueryCELL(J-1,K-1)->owner;
                  remote_task_idx = std::get<1>(*intask_idx_it);
                  taskptr->out_dependencies[owner].push_back(remote_task_idx);
                }
                else {
                  //find the UPDATE2D_COMP and add cell J,I as incoming dependency
                  auto k1 = scheduling::key_t(J,K,I,Factorization::op_type::UPDATE2D_COMP,this->iam);
                  auto outtask_idx_it = std::lower_bound(task_idx.begin(), task_idx.end(), k1, key_lb_comp);
                  auto taskptr = task_graph[std::get<1>(*outtask_idx_it)].get();
                  bassert(taskptr!=nullptr);

                  if (  pQueryCELL(I-1,I-1)->owner != pQueryCELL(J-1,K-1)->owner )
                    taskptr->in_remote_dependencies_cnt++;
                  else
                    taskptr->in_local_dependencies_cnt++;
                }

              }
              break;

            case Factorization::op_type::TRSM_SEND:
              {
                auto J = this->SupMembership_[facing_first_row-1];

                auto k1 = scheduling::key_t(I,I,I,Factorization::op_type::FACTOR,this->iam);
                auto outtask_idx_it = std::lower_bound(task_idx.begin(), task_idx.end(), k1, key_lb_comp);
                auto taskptr = task_graph[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr);
                auto k2 = scheduling::key_t(this->SupMembership_[std::get<4>(taskptr->_meta)-1], std::get<1>(taskptr->_meta), std::get<0>(taskptr->_meta), std::get<2>(taskptr->_meta),this->iam);
                bassert(k2 == k1 );
#ifdef _VERBOSE_
                logfileptr->OFS()<<"TRSM_SEND"<<" from "<<I<<" to "<<I<<" cell ("<<J<<","<<I<<")"<<std::endl;
#endif
                //get index of task TRSM(J,I)
                auto owner = pQueryCELL(J-1,I-1)->owner;
                std::size_t remote_task_idx = 0;
                scheduling::key_t key(J,I,I,Factorization::op_type::TRSM,owner);
                auto task_idx_beg = &task_idx_recvbuf[task_idx_rdispls[owner]];
                auto task_idx_end = task_idx_beg + task_idx_rsizes[owner];
                auto task_idx_it = std::lower_bound(task_idx_beg, task_idx_end, key, key_lb_comp);
                bassert(task_idx_it != task_idx_end); 
                remote_task_idx = std::get<1>(*task_idx_it);
                taskptr->out_dependencies[owner].push_back(remote_task_idx);
#ifdef _VERBOSE_
                logfileptr->OFS()<<"        out dep added to FACTOR"<<" from "<<I<<" to "<<I<<std::endl;
#endif
              }
              break;
            case Factorization::op_type::TRSM_RECV:
              {
                auto J = this->SupMembership_[facing_first_row-1];

                //find the TRSM and add cell I,I as incoming dependency
                auto k1 = scheduling::key_t(J,I,I,Factorization::op_type::TRSM,this->iam);
                auto outtask_idx_it = std::lower_bound(task_idx.begin(), task_idx.end(), k1, key_lb_comp);
                auto taskptr = task_graph[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr); 
                auto k2 = scheduling::key_t(this->SupMembership_[std::get<4>(taskptr->_meta)-1], std::get<1>(taskptr->_meta), std::get<0>(taskptr->_meta), std::get<2>(taskptr->_meta),this->iam);
                bassert(k2 == k1 );
#ifdef _VERBOSE_
                logfileptr->OFS()<<"TRSM_RECV"<<" from "<<I<<" to "<<I<<" cell ("<<J<<","<<I<<")"<<std::endl;
#endif
                if (  pQueryCELL(J-1,I-1)->owner != pQueryCELL(I-1,I-1)->owner )
                  taskptr->in_remote_dependencies_cnt++;
                else
                  taskptr->in_local_dependencies_cnt++;
#ifdef _VERBOSE_
                logfileptr->OFS()<<"        in dep added to TRSM"<<" from "<<I<<" to "<<I<<std::endl;
#endif
              }
              break;
            case Factorization::op_type::AGGREGATE2D_RECV:
              {
                auto K = this->SupMembership_[facing_first_row-1];

#ifdef _VERBOSE_
                logfileptr->OFS()<<"AGGREGATE2D_RECV"<<" from "<<I<<" to "<<J<<" cell ("<<K<<","<<J<<") to cell("<<K<<","<<J<<")"<<std::endl;
#endif
                //find the FACTOR or TRSM and add cell K,J as incoming dependency
                if (K==J) {
                  auto k1 = scheduling::key_t(J,J,J,Factorization::op_type::FACTOR,this->iam);
                  auto outtask_idx_it = std::lower_bound(task_idx.begin(), task_idx.end(), k1, key_lb_comp);
                  auto taskptr = task_graph[std::get<1>(*outtask_idx_it)].get();
                  bassert(taskptr!=nullptr); 
                  auto k2 = scheduling::key_t(this->SupMembership_[std::get<4>(taskptr->_meta)-1], std::get<1>(taskptr->_meta), std::get<0>(taskptr->_meta), std::get<2>(taskptr->_meta),this->iam);
                  bassert(k2 == k1 );

                  taskptr->in_local_dependencies_cnt++;
#ifdef _VERBOSE_
                  logfileptr->OFS()<<"        in dep added to FACTOR"<<" from "<<J<<" to "<<J<<std::endl;
#endif
                }
                else {
                  auto k1 = scheduling::key_t(K,J,J,Factorization::op_type::TRSM,this->iam);
                  auto outtask_idx_it = std::lower_bound(task_idx.begin(), task_idx.end(), k1, key_lb_comp);
                  auto taskptr = task_graph[std::get<1>(*outtask_idx_it)].get();
                  bassert(taskptr!=nullptr); 
                  auto k2 = scheduling::key_t(this->SupMembership_[std::get<4>(taskptr->_meta)-1], std::get<1>(taskptr->_meta), std::get<0>(taskptr->_meta), std::get<2>(taskptr->_meta),this->iam);
                  bassert(k2 == k1 );
                  taskptr->in_local_dependencies_cnt++;
#ifdef _VERBOSE_
                  logfileptr->OFS()<<"        in dep added to TRSM"<<" from "<<J<<" to "<<J<<std::endl;
#endif
                }
              }
              break;
            case Factorization::op_type::AGGREGATE2D_SEND:
              {
                //TODO this might be useless
                auto K = this->SupMembership_[facing_first_row-1];

#ifdef _VERBOSE_
                logfileptr->OFS()<<"AGGREGATE2D_SEND"<<" from "<<I<<" to "<<J<<" cell ("<<K<<","<<J<<") to cell("<<K<<","<<J<<")"<<std::endl;
#endif
                //find the UPDATE2D_COMP and add cell K,J as outgoing dependency

                auto k1 = scheduling::key_t(K,J,I,Factorization::op_type::UPDATE2D_COMP,this->iam);
                auto outtask_idx_it = std::lower_bound(task_idx.begin(), task_idx.end(), k1, key_lb_comp);
                auto taskptr = task_graph[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr); 
                auto k2 = scheduling::key_t(this->SupMembership_[std::get<4>(taskptr->_meta)-1], std::get<1>(taskptr->_meta), std::get<0>(taskptr->_meta), std::get<2>(taskptr->_meta),this->iam);
                bassert(k2 == k1 );

                //find the FACTOR or TRSM target task
                auto owner = pQueryCELL(K-1,J-1)->owner;
                auto task_idx_beg = &task_idx_recvbuf[task_idx_rdispls[owner]];
                auto task_idx_end = task_idx_beg + task_idx_rsizes[owner];

                std::size_t remote_task_idx = 0;
                if (K==J) {
                  scheduling::key_t key(J,J,J,Factorization::op_type::FACTOR,owner);
                  auto task_idx_it = std::lower_bound(task_idx_beg, task_idx_end, (key), key_lb_comp);
                  bassert(task_idx_it != task_idx_end); 
                  remote_task_idx = std::get<1>(*task_idx_it);
                }
                else {
                  scheduling::key_t key(K,J,J,Factorization::op_type::TRSM,owner);
                  auto task_idx_it = std::lower_bound(task_idx_beg, task_idx_end, (key), key_lb_comp);
                  bassert(task_idx_it != task_idx_end); 
                  remote_task_idx = std::get<1>(*task_idx_it);
                }

                taskptr->out_dependencies[owner].push_back(remote_task_idx);
#ifdef _VERBOSE_
                logfileptr->OFS()<<"        out dep added to UPDATE_COMP"<<" from "<<I<<" to "<<J<<std::endl;
#endif
              }
              break;

            case Factorization::op_type::UPDATE2D_RECV:
              {
                //TODO this might be useless
                Idx K = this->SupMembership_[facing_first_row-1];
                std::swap(K,J);
#ifdef _VERBOSE_
                logfileptr->OFS()<<"UPDATE_RECV"<<" from "<<I<<" to "<<K<<" cell ("<<J<<","<<I<<") to cell ("<<J<<","<<K<<")"<<std::endl;
#endif
                //find the UPDATE2D_COMP and add cell J,I as incoming dependency
                auto k1 = scheduling::key_t(J,K,I,Factorization::op_type::UPDATE2D_COMP,this->iam);
                auto outtask_idx_it = std::lower_bound(task_idx.begin(), task_idx.end(), k1, key_lb_comp);
                auto taskptr = task_graph[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr);
                auto k2 = scheduling::key_t(this->SupMembership_[std::get<4>(taskptr->_meta)-1], std::get<1>(taskptr->_meta), std::get<0>(taskptr->_meta), std::get<2>(taskptr->_meta),this->iam);
                bassert(k2 == k1 );

                if (  pQueryCELL(J-1,I-1)->owner != pQueryCELL(J-1,K-1)->owner )
                  taskptr->in_remote_dependencies_cnt++;
                else
                  taskptr->in_local_dependencies_cnt++;

#ifdef _VERBOSE_
                logfileptr->OFS()<<"        in dep added to UPDATE_COMP"<<" from "<<I<<" to "<<K<<std::endl;
#endif
              }
              break;
            case Factorization::op_type::UPDATE2D_SEND:
              {
                Idx K = this->SupMembership_[facing_first_row-1];
                std::swap(K,J);
#ifdef _VERBOSE_
                logfileptr->OFS()<<"UPDATE_SEND"<<" from "<<I<<" to "<<K<<" cell ("<<J<<","<<I<<") to cell ("<<J<<","<<K<<")"<<std::endl;
#endif
                //find the TRSM and add cell J,K as outgoing dependency
                auto k1 = scheduling::key_t(J,I,I,Factorization::op_type::TRSM,this->iam);
                auto outtask_idx_it = std::lower_bound(task_idx.begin(), task_idx.end(), k1, key_lb_comp);
                auto taskptr = task_graph[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr); 
                auto k2 = scheduling::key_t(this->SupMembership_[std::get<4>(taskptr->_meta)-1], std::get<1>(taskptr->_meta), std::get<0>(taskptr->_meta), std::get<2>(taskptr->_meta),this->iam);
                bassert(k2 == k1 );

                auto owner = pQueryCELL(J-1,K-1)->owner;
                auto task_idx_beg = &task_idx_recvbuf[task_idx_rdispls[owner]];
                auto task_idx_end = task_idx_beg + task_idx_rsizes[owner];

                std::size_t remote_task_idx = 0;
                scheduling::key_t key(J,K,I,Factorization::op_type::UPDATE2D_COMP,owner);
                auto task_idx_it = std::lower_bound(task_idx_beg, task_idx_end, (key), key_lb_comp);
                bassert(task_idx_it != task_idx_end); 
                remote_task_idx = std::get<1>(*task_idx_it);

                taskptr->out_dependencies[owner].push_back(remote_task_idx);

#ifdef _VERBOSE_
                logfileptr->OFS()<<"        out dep added to TRSM"<<" from "<<I<<" to "<<I<<std::endl;
#endif
              }
              break;
            case Factorization::op_type::UPDATE2D_SEND_OD:
              {
                auto K = this->SupMembership_[lt_first_row-1];
#ifdef _VERBOSE_
                logfileptr->OFS()<<"UPDATE_SEND OD"<<" from "<<I<<" to "<<J<<" cell ("<<J<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;
#endif

                //find the TRSM and add cell K,J as outgoing dependency
                auto k1 = scheduling::key_t(J,I,I,Factorization::op_type::TRSM,this->iam);
                auto outtask_idx_it = std::lower_bound(task_idx.begin(), task_idx.end(), k1, key_lb_comp);
                auto taskptr = task_graph[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr); 
                auto k2 = scheduling::key_t(this->SupMembership_[std::get<4>(taskptr->_meta)-1], std::get<1>(taskptr->_meta), std::get<0>(taskptr->_meta), std::get<2>(taskptr->_meta),this->iam);
                bassert(k2 == k1 );

                auto owner = pQueryCELL(K-1,J-1)->owner;
                auto task_idx_beg = &task_idx_recvbuf[task_idx_rdispls[owner]];
                auto task_idx_end = task_idx_beg + task_idx_rsizes[owner];

                std::size_t remote_task_idx = 0;
                scheduling::key_t key(K,J,I,Factorization::op_type::UPDATE2D_COMP,owner);
                auto task_idx_it = std::lower_bound(task_idx_beg, task_idx_end, (key), key_lb_comp);
                bassert(task_idx_it != task_idx_end); 
                remote_task_idx = std::get<1>(*task_idx_it);
                taskptr->out_dependencies[owner].push_back(remote_task_idx);

                if (this->options_.decomposition == DecompositionType::LDL) {
                  key = scheduling::key_t(J,J,I,Factorization::op_type::DLT2D_COMP,owner);
                  task_idx_it = std::lower_bound(task_idx_beg, task_idx_end, (key), key_lb_comp);
                  bassert(task_idx_it != task_idx_end); 
                  remote_task_idx = std::get<1>(*task_idx_it);
                  taskptr->out_dependencies[owner].push_back(remote_task_idx);

                }

#ifdef _VERBOSE_
                logfileptr->OFS()<<"        out dep DOWN added to TRSM"<<" from "<<I<<" to "<<I<<std::endl;
#endif
              }
              break;
            case Factorization::op_type::UPDATE2D_RECV_OD:
              {
                auto K = this->SupMembership_[lt_first_row-1];
#ifdef _VERBOSE_
                logfileptr->OFS()<<"UPDATE_RECV OD"<<" from "<<I<<" to "<<J<<" cell ("<<J<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;
#endif

                //find the UPDATE2D_COMP and add cell J,I as incoming dependency
                auto k1 = scheduling::key_t(K,J,I,Factorization::op_type::UPDATE2D_COMP,this->iam);

                auto outtask_idx_it = std::lower_bound(task_idx.begin(), task_idx.end(), k1, key_lb_comp);
                auto taskptr = task_graph[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr);
                auto k2 = scheduling::key_t(this->SupMembership_[std::get<4>(taskptr->_meta)-1], std::get<1>(taskptr->_meta), std::get<0>(taskptr->_meta), std::get<2>(taskptr->_meta),this->iam);
                bassert(k2 == k1 );

                if (  pQueryCELL(J-1,I-1)->owner != pQueryCELL(K-1,J-1)->owner )
                  taskptr->in_remote_dependencies_cnt++;
                else
                  taskptr->in_local_dependencies_cnt++;


                if ( this->options_.decomposition == DecompositionType::LDL ) {
                  auto k3 = scheduling::key_t(J,J,I,Factorization::op_type::DLT2D_COMP,this->iam);
                  auto dlttask_idx_it = std::lower_bound(task_idx.begin(), task_idx.end(), k3, key_lb_comp);
                  auto dlttaskptr = task_graph[std::get<1>(*dlttask_idx_it)].get();
                  bassert(dlttaskptr!=nullptr);
                  bassert(std::get<2>(dlttaskptr->_meta)==Factorization::op_type::DLT2D_COMP);
                  if (  pQueryCELL(J-1,I-1)->owner != pQueryCELL(K-1,J-1)->owner )
                    dlttaskptr->in_remote_dependencies_cnt++;
                  else
                    dlttaskptr->in_local_dependencies_cnt++;

                  //now add a local dependency for the UPDATE2D_COMP task
                  taskptr->in_local_dependencies_cnt++;
                  auto owner = pQueryCELL(K-1,J-1)->owner;
                  std::size_t remote_task_idx = std::get<1>(*outtask_idx_it);
                  dlttaskptr->out_dependencies[owner].push_back(remote_task_idx);
                }





#ifdef _VERBOSE_
                logfileptr->OFS()<<"        in dep DOWN added to UPDATE2D_COMP"<<" from "<<I<<" to "<<J<<std::endl;
#endif
              }
              break;
            default:
              break;
          }
        }

        //fetch the addresses of diagonal cells if necessary (LDL) ?
        //Allgatherv

        if (this->options_.decomposition == DecompositionType::LDL) {
          using diag_ptr_t = std::tuple<Int,upcxx::global_ptr<char>>;
          MPI_Datatype upcxx_ptr_type;
          MPI_Type_contiguous( sizeof(diag_ptr_t), MPI_BYTE, &upcxx_ptr_type );
          MPI_Type_commit(&upcxx_ptr_type);
          std::vector<int> rsizes(this->np,0);
          std::vector<int> rdispls(this->np+1,0);

          int ssize = local_diag_pointers.size();
          MPI_Allgather(&ssize,sizeof(ssize),MPI_BYTE,rsizes.data(),sizeof(ssize),MPI_BYTE,this->workcomm_);

          rdispls[0] = 0;
          std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);
          this->diag_pointers_.resize(rdispls.back());

          //now communicate the task_idx arrays: allgatherv
          MPI_Allgatherv(local_diag_pointers.data(),local_diag_pointers.size(),upcxx_ptr_type,this->diag_pointers_.data(),rsizes.data(),rdispls.data(),upcxx_ptr_type,this->workcomm_);
          MPI_Type_free(&upcxx_ptr_type);

          std::sort(this->diag_pointers_.begin(),this->diag_pointers_.end(),[](diag_ptr_t & a, diag_ptr_t & b) { return std::get<0>(a) < std::get<0>(b) ; });

        }




        scheduler.sp_handle = this->sp_handle;
        //Now we have our local part of the task graph
        for (auto it = this->task_graph.begin(); it != this->task_graph.end(); it++) {
          auto & ptask = *it;
          auto meta = &ptask->_meta;

          auto & src_snode = std::get<0>(meta[0]);
          auto & tgt_snode = std::get<1>(meta[0]);
          auto & lt_first_row = std::get<3>(meta[0]);
          auto & facing_row = std::get<4>(meta[0]);
          auto & type = std::get<2>(meta[0]);

          auto I = src_snode;
          auto J = tgt_snode;
          auto K = this->SupMembership_[facing_row-1];

          auto remote_deps = ptask->in_remote_dependencies_cnt;
          auto local_deps = ptask->in_local_dependencies_cnt;

#ifdef _USE_PROM_AVAIL_
          ptask->in_avail_prom.require_anonymous(remote_deps);
          ptask->in_avail_counter = remote_deps;
#else
          ptask->in_avail_counter = remote_deps;
#endif
#ifdef _USE_PROM_RDY_
          ptask->in_prom.require_anonymous(local_deps + remote_deps);
          ptask->in_counter = local_deps + remote_deps;
#else
          ptask->in_counter = local_deps + remote_deps;
#endif

          auto ptr = ptask.get();
          switch(type) {
            case Factorization::op_type::FACTOR:
              {
                ptask->execute = [this,src_snode,ptr,I,J,K] () {
                  scope_timer(b,FB_FACTOR_DIAG_TASK);
                  auto ptask = ptr;

                  auto ptr_diagcell = pQueryCELL2(I-1,I-1);
                  assert(ptr_diagcell);

#ifdef SP_THREADS
                  std::thread::id tid = std::this_thread::get_id();
                  auto & tmpBuf = tmpBufs_th[tid];
#else
                  auto & tmpBuf = tmpBufs;
#endif
#ifdef _TIMING_
                  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif

                  ptr_diagcell->factorize(tmpBuf);
#ifdef _TIMING_
                  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
                  comp_fact_ticks += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
#endif
                  //Scan through all dependent tasks, set metadata of dependent tasks, including the global pointer of the dependent task. (which will be on the src processor)
                  for (auto it = ptask->out_dependencies.begin(); it!=ptask->out_dependencies.end(); it++) {
                    auto pdest = it->first;
                    auto & tgt_cells = it->second;
                    //serialize data once, and list of meta data
                    //factor is output data so it will not be deleted
#ifdef CUDA_MODE
                    bool is_gpu_block = false;
                    if (pdest != this->iam 
                        && ptr_diagcell->_nnz*sizeof(T)>symPACK::gpu_block_limit
                        && tgt_cells.size()>0) {
                            is_gpu_block = true;
                    }
#endif
                    if ( pdest != this->iam ) {
                      upcxx::rpc_ff( pdest,
#ifdef CUDA_MODE
                          [is_gpu_block]
#else 
                          [] 
#endif
                          (int sp_handle, upcxx::global_ptr<char> gptr, 
                               size_t storage_size, 
                               size_t nnz, size_t nblocks, rowind_t width, 
                               SparseTask2D::meta_t meta, upcxx::view<std::size_t> target_cells ) { 
                              //there is a map between sp_handle and task_graphs
#ifdef _TIMING_
                              gasneti_tick_t start = gasneti_ticks_now();
#endif
                              auto matptr = (symPACKMatrix2D<colptr_t,rowind_t,T> *) g_sp_handle_to_matrix[sp_handle];
                              auto I = std::get<1>(meta);
                              rowind_t fc = matptr->Xsuper_[I-1];

                              //store pointer & associated metadata somewhere
                              std::shared_ptr <SparseTask2D::data_t > data;
                              std::shared_ptr <SparseTask2D::data_t > diag_data;

                              for ( auto & tgt_cell_idx: target_cells) {
                                  auto taskptr = matptr->task_graph[tgt_cell_idx].get();

                                  if ( std::get<2>(taskptr->_meta) == Factorization::op_type::DLT2D_COMP
                                      || std::get<2>(taskptr->_meta) == Factorization::op_type::UPDATE2D_COMP
                                     ) {
                                      bassert(std::get<2>(taskptr->_meta) == Factorization::op_type::DLT2D_COMP);
                                      bassert (matptr->options_.decomposition == DecompositionType::LDL);
                                      if ( !diag_data ) {
                                        upcxx::global_ptr<char> diag_ptr = matptr->find_diag_pointer( I );
                                        bassert ( diag_ptr.where() != upcxx::rank_me() );
                                        diag_data = std::make_shared<SparseTask2D::data_t >();
                                        diag_data->in_meta = std::make_tuple(I,I,Factorization::op_type::DIAG_ENTRIES,0,0);;
                                        diag_data->size = (matptr->Xsuper_[I] - matptr->Xsuper_[I-1])*sizeof(T);
                                        diag_data->remote_gptr = diag_ptr;
                                      }

                                      taskptr->input_msg.push_back(diag_data);
                                      diag_data->target_tasks.push_back(taskptr);
                                      taskptr->avail_dep(1,matptr->scheduler);
                                  } else {
                                    
                                    if ( ! data ) {
                                      data = std::make_shared<SparseTask2D::data_t >();
                                      data->in_meta = meta;
                                      data->size = storage_size;
                                      data->remote_gptr = gptr;
#ifdef CUDA_MODE
                                      if (is_gpu_block) {
                                        data->is_gpu_block = true;
                                        data->d_size = nnz;
                                      }
#endif
                                    }
                                    taskptr->input_msg.push_back(data);
                                    data->target_tasks.push_back(taskptr);
                                    taskptr->avail_dep(1,matptr->scheduler);
                                  }
                              }


                              if ( data ) {
                                data->on_fetch_future = data->on_fetch_future.then(
                                    [fc,width,nnz,nblocks,I,matptr,data](SparseTask2D::data_t * pdata) {
#if not defined(_NO_COMPUTATION_)
                                    //create snodeBlock_t and store it in the extra_data
                                    if (matptr->options_.decomposition == DecompositionType::LDL) {
                                      pdata->extra_data = std::shared_ptr<blockCellBase_t>( (blockCellBase_t*)new snodeBlockLDL_t(I,I,pdata->landing_zone,fc,width,nnz,nblocks) );
                                    }
                                    else {
#ifdef CUDA_MODE
                                      if (data->is_gpu_block) {
                                        pdata->extra_data = std::shared_ptr<blockCellBase_t>( 
                                        (blockCellBase_t*)new snodeBlock_t(I,I,pdata->landing_zone,pdata->d_landing_zone,fc,width,nnz,nblocks) );
                                      } else {
                                        pdata->extra_data = std::shared_ptr<blockCellBase_t>( 
                                        (blockCellBase_t*)new snodeBlock_t(I,I,pdata->landing_zone,fc,width,nnz,nblocks) );
                                      }
#else
                                      pdata->extra_data = std::shared_ptr<blockCellBase_t>( (blockCellBase_t*)new snodeBlock_t(I,I,pdata->landing_zone,fc,width,nnz,nblocks) );
#endif
                                    }
#endif
                                    return upcxx::to_future(pdata);
                                    });
#ifdef _EAGER_FETCH_
                                data->allocate();
                                data->fetch();
#endif
                              }

#ifdef _EAGER_FETCH_
                              if ( diag_data ) {
                                diag_data->allocate();
                                diag_data->fetch();
                              }
#endif

#ifdef _TIMING_
                              matptr->rpc_fact_ticks += gasneti_ticks_to_ns(gasneti_ticks_now() - start);
#endif
                          }, this->sp_handle, ptr_diagcell->_gstorage, 
                             ptr_diagcell->_storage_size, ptr_diagcell->nnz(), 
                             ptr_diagcell->nblocks(), std::get<0>(ptr_diagcell->_dims) ,ptask->_meta, 
                             upcxx::make_view(tgt_cells.begin(),tgt_cells.end()));                           
                    }
                    else {
                      for ( auto & tgt_cell_idx: tgt_cells ) {
                        auto taskptr = task_graph[tgt_cell_idx].get();
                        bassert(taskptr!=nullptr); 
                        bassert(std::get<2>(taskptr->_meta)==Factorization::op_type::TRSM
                            ||  std::get<2>(taskptr->_meta) == Factorization::op_type::DLT2D_COMP
                            ||  std::get<2>(taskptr->_meta) == Factorization::op_type::UPDATE2D_COMP);
                        //mark the dependency as satisfied
                        taskptr->satisfy_dep(1,this->scheduler);
                      }
                    }
                  }

#ifdef _TIMING_
                  deps_fact_ticks += gasneti_ticks_to_ns(gasneti_ticks_now() - start);
#endif
                  //TODO what do do with the task OR return a future
                  ptask->executed = true;
                };
              }
              break;
            case Factorization::op_type::TRSM:
              {

                ptask->execute = [this,src_snode,tgt_snode,ptr,I,J,K] () {
                  scope_timer(b,FB_TRSM_TASK);
                  auto ptask = ptr;
                  auto ptr_od_cell = pQueryCELL2(K-1,I-1);
                  assert(ptr_od_cell);

                  bassert( ptr_od_cell->owner == this->iam);

                  {
#ifdef SP_THREADS
                    std::thread::id tid = std::this_thread::get_id();
                    cell_lock<snodeBlock_sptr_t> lock(ptr_od_cell);
                    auto & tmpBuf = tmpBufs_th[tid];
#else
                    auto & tmpBuf = tmpBufs;
#endif


                    //input data is one or 0 (local diagonal block)
                    bassert(ptask->input_msg.size()<=1);
                    auto ptr_diagCell = pQueryCELL(I-1,I-1).get(); 
#if not defined(_NO_COMPUTATION_)
                    if ( ptr_diagCell->owner != this->iam ) {
                      ptr_diagCell = (snodeBlock_t*)(ptask->input_msg.begin()->get()->extra_data.get());
                    }
                    bassert(ptr_diagCell!=nullptr);
#endif

#ifdef _TIMING_
                    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif
                    //TODO DEBUG
                    ptr_od_cell->trsm(ptr_diagCell,tmpBuf);
#ifdef _TIMING_
                    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
                    comp_trsm_ticks += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
#endif
                  }
#ifdef _TIMING_
                  start = gasneti_ticks_now();
#endif
                  //Iterate thru dependent tasks
                  for (auto it = ptask->out_dependencies.begin(); it!=ptask->out_dependencies.end(); it++) {
                    auto pdest = it->first;
                    auto & tgt_cells = it->second;
                    //serialize data once, and list of meta data
                    //factor is output data so it will not be deleted
                    if ( pdest != this->iam ) {
                      upcxx::rpc_ff( pdest,
                          [K,I] (int sp_handle, upcxx::global_ptr<char> gptr, size_t storage_size, 
                          size_t nnz, size_t nblocks, rowind_t width, SparseTask2D::meta_t meta, upcxx::view<std::size_t> target_cells ) { 
#ifdef _TIMING_
                              gasneti_tick_t start = gasneti_ticks_now();
#endif
                              //store pointer & associated metadata somewhere
                              auto data = std::make_shared<SparseTask2D::data_t >();
                              data->in_meta = meta;
                              data->size = storage_size;
                              data->remote_gptr = gptr;

                              //there is a map between sp_handle and task_graphs
                              auto matptr = (symPACKMatrix2D<colptr_t,rowind_t,T> *) g_sp_handle_to_matrix[sp_handle];

                              for ( auto & tgt_cell_idx: target_cells) {
                              auto taskptr = matptr->task_graph[tgt_cell_idx].get();

                              bassert( std::get<2>(taskptr->_meta)==Factorization::op_type::DLT2D_COMP
                                  || std::get<2>(taskptr->_meta)==Factorization::op_type::UPDATE2D_COMP);

                              taskptr->input_msg.push_back(data);
                              data->target_tasks.push_back(taskptr);
                              taskptr->avail_dep(1,matptr->scheduler);
                              }

                              auto owner = gptr.where();
                              rowind_t fc = matptr->Xsuper_[I-1];
                              data->on_fetch_future = data->on_fetch_future.then(
                                  [fc,matptr,width,nnz,nblocks,K,I,owner](SparseTask2D::data_t * pdata) {
#if not defined(_NO_COMPUTATION_)
                                  //create snodeBlock_t and store it in the extra_data
                                  if (matptr->options_.decomposition == DecompositionType::LDL) {
                                    pdata->extra_data = std::shared_ptr<blockCellBase_t>( (blockCellBase_t*)new snodeBlockLDL_t(K,I,pdata->landing_zone,fc,width,nnz,nblocks) );
                                  }
                                  else {
                                    pdata->extra_data = std::shared_ptr<blockCellBase_t>( (blockCellBase_t*)new snodeBlock_t(K,I,pdata->landing_zone,fc,width,nnz,nblocks) );
                                  }
                                  
                                  pdata->extra_data->owner = owner;
#endif
                                  return upcxx::to_future(pdata); 
                                  }
                                  );
#ifdef _EAGER_FETCH_
                              data->allocate();
                              data->fetch();
#endif

#ifdef _TIMING_
                              matptr->rpc_trsm_ticks += gasneti_ticks_to_ns(gasneti_ticks_now() - start);
#endif
                          }, this->sp_handle, ptr_od_cell->_gstorage,ptr_od_cell->_storage_size, ptr_od_cell->nnz(),ptr_od_cell->nblocks(), std::get<0>(ptr_od_cell->_dims),ptask->_meta, upcxx::make_view(tgt_cells.begin(),tgt_cells.end())); 
                    }
                    else {
                      std::shared_ptr<SparseTask2D::data_t> diag_data;
                      if (this->options_.decomposition == DecompositionType::LDL) {
                        auto ptr_ldl = (snodeBlockLDL_t*)(ptr_od_cell.get());
                        ptr_ldl->local_pivot = 0;
                      }
                      for ( auto & tgt_cell_idx: tgt_cells ) {
                        auto taskptr = task_graph[tgt_cell_idx].get();
                        bassert(taskptr!=nullptr); 
                        bassert( std::get<2>(taskptr->_meta)==Factorization::op_type::DLT2D_COMP
                            || std::get<2>(taskptr->_meta)==Factorization::op_type::UPDATE2D_COMP);


                        if (std::get<2>(taskptr->_meta)==Factorization::op_type::UPDATE2D_COMP &&
                            this->options_.decomposition == DecompositionType::LDL) {
                          auto ptr_ldl = (snodeBlockLDL_t*)(ptr_od_cell.get());
                          //get coordinates of target cell and check if ptr_od_cell is a pivot or a facing cell
                          auto tJ = std::get<1>(taskptr->_meta);
                          if ( ptr_ldl->i == tJ ) {
                            ptr_ldl->local_pivot++;
                          }
                        }

                        //mark the dependency as satisfied
                        taskptr->satisfy_dep(1,this->scheduler);
                      }
                    }
                  }

#ifdef _TIMING_
                  deps_trsm_ticks += gasneti_ticks_to_ns(gasneti_ticks_now() - start);
#endif

                  ptask->executed = true;

                };
              }
              break;
            case Factorization::op_type::DLT2D_COMP:
              {
                ptask->execute = [this,src_snode,ptr,I,J,K] () {
                  scope_timer(b,FB_DLT2D_TASK);
                  auto ptask = ptr;

                  bassert(this->options_.decomposition == DecompositionType::LDL);

                  T * diag_ptr = nullptr;
                  auto ptr_odCell = pQueryCELL(J-1,I-1).get(); 
                  bool odSet = false;
                  bool diagSet = false;
                  for ( auto && msg_ptr: ptask->input_msg ) {
                    if ( msg_ptr->extra_data ) {
                      if ( ptr_odCell->owner != this->iam && !odSet) {
                        odSet = true;
                        ptr_odCell = (snodeBlock_t*)(msg_ptr->extra_data.get());
                      }
#ifndef NDEBUG
                      else if (odSet && msg_ptr->extra_data.get() != ptr_odCell){
                        gdb_lock();
                      }
#endif
                    }
                    else {
                      diagSet = true;
                      diag_ptr = (T*)msg_ptr->landing_zone;
                    }
                  }



                  if ( ! diagSet) {
                    auto ptr_diagCell = pQueryCELL(I-1,I-1).get();
                    bassert(ptr_diagCell->owner == this->iam);
                    diag_ptr = (T*)dynamic_cast<snodeBlockLDL_t*>(ptr_diagCell)->GetDiag();
                  }

                  auto ptr_odCellLDL =  dynamic_cast<snodeBlockLDL_t*>(ptr_odCell);

                  if ( ptr_odCellLDL->owner == this->iam ) {
                    bassert( ptr_odCellLDL->local_pivot>0 );

                    if ( ptr_odCellLDL->_bufLDL == nullptr ) {
                      ptr_odCellLDL->computeDLT(diag_ptr);
                    }
                  }
                  else{
                    if ( ptr_odCellLDL->_bufLDL == nullptr ) {
                      ptr_odCellLDL->computeDLT(diag_ptr);
                    }
                  }

                  //Signal the LOCAL dependencies
                  for (auto it = ptask->out_dependencies.begin(); it!=ptask->out_dependencies.end(); it++) {
                    auto pdest = it->first;
                    auto & tgt_cells = it->second;
                    //serialize data once, and list of meta data
                    //factor is output data so it will not be deleted
                    bassert ( pdest == this->iam );
                    for ( auto & tgt_cell_idx: tgt_cells ) {
                      auto taskptr = task_graph[tgt_cell_idx].get();
                      bassert(taskptr!=nullptr); 
                      bassert(std::get<2>(taskptr->_meta)==Factorization::op_type::UPDATE2D_COMP);
                      //mark the dependency as satisfied
                      taskptr->satisfy_dep(1,this->scheduler);
                    }

                  }
                };
              }
              break;
            case Factorization::op_type::UPDATE2D_COMP:
              {

                ptask->execute = [this,src_snode,ptr,I,J,K] () {
                  scope_timer(b,FB_UPDATE2D_TASK);
                  auto ptask = ptr;

                  auto ptr_upd_cell = pQueryCELL2(K-1,J-1);
                  assert(ptr_upd_cell);
                  //TODO false for fan-both
                  bassert(ptr_upd_cell->owner==this->iam);

                  {
#ifdef SP_THREADS
                    cell_lock<snodeBlock_sptr_t> lock(ptr_upd_cell);
#endif
#ifdef SP_THREADS
                    std::thread::id tid = std::this_thread::get_id();
                    auto & tmpBuf = tmpBufs_th[tid];
#else
                    auto & tmpBuf = tmpBufs;
#endif

                    //input data should be at most two
                    bassert(ptask->input_msg.size()<=2 );

                    auto ptr_odCell = pQueryCELL(J-1,I-1).get(); 
                    auto ptr_facingCell = pQueryCELL(K-1,I-1).get(); 
                    bool odSet = false;
                    bool facingSet = false;
#if not defined(_NO_COMPUTATION_)
                    for ( auto && msg_ptr: ptask->input_msg ) {
                      bassert ( msg_ptr->extra_data != nullptr );
                      if ( ptr_odCell->owner != this->iam && !odSet) {
                        odSet = true;
                        ptr_odCell = (snodeBlock_t*)(msg_ptr->extra_data.get());
                      }
                      else if ( ptr_facingCell->owner != this->iam && !facingSet) {
                        facingSet = true;
                        ptr_facingCell = (snodeBlock_t*)(msg_ptr->extra_data.get());
                      }
                    }

                    bassert( !odSet || ptr_odCell != pQueryCELL(J-1,I-1).get() );
                    bassert( !facingSet || ptr_facingCell != pQueryCELL(K-1,I-1).get() );

                    bassert( ptr_facingCell->owner == this->iam || ptr_facingCell != pQueryCELL(K-1,I-1).get() );
                    bassert( ptr_odCell->owner == this->iam || ptr_odCell != pQueryCELL(J-1,I-1).get() );

                    bassert(ptr_odCell!=nullptr);
                    bassert(ptr_facingCell!=nullptr);
                    //check that they don't need to be swapped
                    if ( ptr_facingCell->i < ptr_odCell->i ) {
                      std::swap(ptr_facingCell, ptr_odCell);
                    }

#endif

#ifdef _TIMING_
                    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif

                    if (this->options_.decomposition == DecompositionType::LDL) {
#ifndef NDEBUG
#if not defined(_NO_COMPUTATION_)
                      auto ptr_od = dynamic_cast<snodeBlockLDL_t*>(ptr_odCell);
                      auto ptr_facing = dynamic_cast<snodeBlockLDL_t*>(ptr_facingCell);
                      bassert(ptr_od!=nullptr);
                      bassert(ptr_facing!=nullptr);
#endif
#endif
                      ptr_upd_cell->update(ptr_odCell,ptr_facingCell,tmpBuf);
                    }
                    else {
                      ptr_upd_cell->update(ptr_odCell,ptr_facingCell,tmpBuf);
                    }

#ifdef _TIMING_
                    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
                    comp_upd_ticks += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
#endif

                  }

                  for (auto it = ptask->out_dependencies.begin(); it!=ptask->out_dependencies.end(); it++) {
                    auto pdest = it->first;
                    auto & tgt_cells = it->second;
                    //serialize data once, and list of meta data
                    //factor is output data so it will not be deleted
                    if ( pdest != this->iam ) {
                      upcxx::rpc_ff( pdest, 
                          [K,J] (int sp_handle, upcxx::global_ptr<char> gptr, size_t storage_size, size_t nnz, size_t nblocks, rowind_t width, SparseTask2D::meta_t meta, upcxx::view<std::size_t> target_cells ) { 
#ifdef _TIMING_
                              gasneti_tick_t start = gasneti_ticks_now();
#endif
                              //store pointer & associated metadata somewhere
                              auto data = std::make_shared<SparseTask2D::data_t >();
                              data->in_meta = meta;
                              data->size = storage_size;
                              data->remote_gptr = gptr;

                              //there is a map between sp_handle and task_graphs
                              auto matptr = (symPACKMatrix2D<colptr_t,rowind_t,T> *) g_sp_handle_to_matrix[sp_handle];
                              for ( auto & tgt_cell_idx: target_cells) {
                                auto taskptr = matptr->task_graph[tgt_cell_idx].get();
                                taskptr->input_msg.push_back(data);
                                data->target_tasks.push_back(taskptr);
                                taskptr->avail_dep(1,matptr->scheduler);
                              }


                              auto I = std::get<1>(meta);
                              rowind_t fc = matptr->Xsuper_[I-1];
                              data->on_fetch_future = data->on_fetch_future.then(
                                  [fc,width,nnz,nblocks,K,J,matptr](SparseTask2D::data_t * pdata) {
#if not defined(_NO_COMPUTATION_)
                                  //create snodeBlock_t and store it in the extra_data
                                  if (matptr->options_.decomposition == DecompositionType::LDL) {
                                    pdata->extra_data = std::shared_ptr<blockCellBase_t>( (blockCellBase_t*)new snodeBlockLDL_t(K,J,pdata->landing_zone,fc,width,nnz,nblocks) );
                                  }
                                  else {
                                    pdata->extra_data = std::shared_ptr<blockCellBase_t>( (blockCellBase_t*)new snodeBlock_t(K,J,pdata->landing_zone,fc,width,nnz,nblocks) );
                                  }
#endif
                                  return upcxx::to_future(pdata);
                                  }
                                  );
                              //TODO check this
#ifdef _EAGER_FETCH_
                              data->allocate();
                              data->fetch();
#endif
#ifdef _TIMING_
                              matptr->rpc_upd_ticks += gasneti_ticks_to_ns(gasneti_ticks_now() - start);
#endif
                          }, this->sp_handle, ptr_upd_cell->_gstorage, ptr_upd_cell->_storage_size,ptr_upd_cell->nnz(),ptr_upd_cell->nblocks(), std::get<0>(ptr_upd_cell->_dims),ptask->_meta, upcxx::make_view(tgt_cells.begin(),tgt_cells.end())); 
                    }
                    else {
                      for ( auto & tgt_cell_idx: tgt_cells ) {
                        auto taskptr = task_graph[tgt_cell_idx].get();
                        bassert(taskptr!=nullptr); 
                        bassert(std::get<2>(taskptr->_meta)==Factorization::op_type::FACTOR || std::get<2>(taskptr->_meta)==Factorization::op_type::TRSM);
                        //mark the dependency as satisfied
                        taskptr->satisfy_dep(1,this->scheduler);
                      }
                    }
                  }

#ifdef _TIMING_
                  deps_upd_ticks += gasneti_ticks_to_ns(gasneti_ticks_now() - start);
#endif

                  ptask->executed = true;
                };

              }
              break;
            default:
              {
                ptask.reset(nullptr);
              }
              break;
          }
        }


        logfileptr->OFS()<<"Task graph created"<<std::endl;
      }

      //generate task graph for the solution phase
      {
        MPI_Datatype type;
        MPI_Type_contiguous( sizeof(SparseTask2D::meta_t), MPI_BYTE, &type );
        MPI_Type_commit(&type);
        vector<int> ssizes(this->np,0);
        vector<int> sdispls(this->np+1,0);
        vector<SparseTask2D::meta_t> sendbuf;

        auto supETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);

        vector<int> ssizes_dep(this->np,0);
        vector<int> sdispls_dep(this->np+1,0);
        vector<SparseTask2D::meta_t> sendbuf_dep;
        {  
          //we will need to communicate if only partial xlindx_, lindx_
          //idea: build tasklist per processor and then exchange
          //tuple is: src_snode,tgt_snode,op_type,lt_first_row = updated fc, facing_first_row = updated fr,
          std::map<Idx, std::list< SparseTask2D::meta_t> > Updates;
          std::map<Idx, std::list< SparseTask2D::meta_t> > Messages;
          std::vector<int> marker(this->np,0);

          Int numLocSnode = this->XsuperDist_[this->iam+1]-this->XsuperDist_[this->iam];
          Int firstSnode = this->XsuperDist_[this->iam];
          for (Int locsupno = 1; locsupno<this->locXlindx_.size(); ++locsupno) {
            Idx I = locsupno + firstSnode-1;
            Int first_col = this->Xsuper_[I-1];
            Int last_col = this->Xsuper_[I]-1;

            bassert( cells_.find(coord2supidx(I-1,I-1)) != cells_.end() );

            auto ptr_fcell = pQueryCELL(I-1,I-1);
            Int iOwner = ptr_fcell->owner;

            //Create the forward update task for the diagonal block
            Updates[iOwner].push_back(std::make_tuple(I,I,Factorization::op_type::FUC,I,first_col));


            Int parent = supETree.Parent(I-1);
            //if root, add a FUC send to a BUC task
            if ( parent == 0 ) {

              Messages[iOwner].push_back(std::make_tuple(I,I,Factorization::op_type::FUC_SEND,I,first_col));
              Messages[iOwner].push_back(std::make_tuple(I,I,Factorization::op_type::FUC_RECV,I,first_col));
            }

            Ptr lfi = this->locXlindx_[locsupno-1];
            Ptr lli = this->locXlindx_[locsupno]-1;
            std::list<Int> ancestor_rows;

            Int K = -1; 
            Idx K_prevSnode = (Idx)-1;
            for (Ptr K_sidx = lfi; K_sidx<=lli;K_sidx++) {
              Idx K_row = this->locLindx_[K_sidx-1];
              K = this->SupMembership_[K_row-1];

              //Split at boundary or after diagonal block
              if (K!=K_prevSnode) {
                if (K>I) {
                  ancestor_rows.push_back(K_row);
                }
              }
              K_prevSnode = K;
            }

            for (auto J_row : ancestor_rows) {
              Int J = this->SupMembership_[J_row-1];

              auto ptr_fodcell = pQueryCELL(J-1,I-1);
              Int iFODOwner = ptr_fodcell->owner;

              Messages[iOwner].push_back(std::make_tuple(J,I,Factorization::op_type::FUC_D_SEND,J_row,J_row));
              Messages[iFODOwner].push_back(std::make_tuple(J,I,Factorization::op_type::FUC_D_RECV,J_row,J_row));

              Updates[iFODOwner].push_back(std::make_tuple(J,I,Factorization::op_type::FUC,J,J_row));

              //send the contrib to the last local ancestor cell(J,last_local_ancestor)
              //if already last, send to cell (J,J) 
              auto ptr_tgtupdcell = pQueryCELL(J-1,J-1);
              Int iTgtOwner = ptr_tgtupdcell->owner;
              //TODO update this for non fan-out mapping
              Int iUpdOwner = ptr_tgtupdcell->owner;

              Messages[iFODOwner].push_back(std::make_tuple(J,I,Factorization::op_type::FUC_SEND,J,J_row));
              Messages[iUpdOwner].push_back(std::make_tuple(J,I,Factorization::op_type::FUC_RECV,J,J_row));
            }

            //Create the backward update task for the diagonal block
            Updates[iOwner].push_back(std::make_tuple(I,I,Factorization::op_type::BUC,I,first_col));

            for (auto J_row : ancestor_rows) {
              Int J = this->SupMembership_[J_row-1];

              auto ptr_fodcell = pQueryCELL(J-1,I-1);
              Int iFODOwner = ptr_fodcell->owner;

              Messages[iOwner].push_back(std::make_tuple(J,I,Factorization::op_type::BUC_D_RECV,J_row,J_row));
              Messages[iFODOwner].push_back(std::make_tuple(J,I,Factorization::op_type::BUC_D_SEND,J_row,J_row));
              Updates[iFODOwner].push_back(std::make_tuple(J,I,Factorization::op_type::BUC,J_row,J_row));

              //send the contrib to cell(J,J)
              auto ptr_srccell = pQueryCELL(J-1,J-1);
              Int iSrcOwner = ptr_srccell->owner;
              //TODO update this for non fan-out mapping
              Int iUpdOwner = ptr_srccell->owner;
              Messages[iFODOwner].push_back(std::make_tuple(J,I,Factorization::op_type::BUC_RECV,J_row,J_row));
              Messages[iUpdOwner].push_back(std::make_tuple(J,I,Factorization::op_type::BUC_SEND,J_row,J_row));
            }
          }

          //then do an alltoallv
          //compute send sizes
          for (auto itp = Updates.begin();itp!=Updates.end();itp++) {
            ssizes[itp->first] = itp->second.size();
          }

          //compute send displacements
          sdispls[0] = 0;
          std::partial_sum(ssizes.begin(),ssizes.end(),&sdispls[1]);

          //Build the contiguous array of pairs
          sendbuf.reserve(sdispls.back());

          for (auto itp = Updates.begin();itp!=Updates.end();itp++) {
            sendbuf.insert(sendbuf.end(),itp->second.begin(),itp->second.end());
          }

          //then do an alltoallv
          //compute send sizes
          for (auto itp = Messages.begin();itp!=Messages.end();itp++) {
            ssizes_dep[itp->first] = itp->second.size();
          }

          //compute send displacements
          sdispls_dep[0] = 0;
          std::partial_sum(ssizes_dep.begin(),ssizes_dep.end(),&sdispls_dep[1]);

          //Build the contiguous array of pairs
          sendbuf_dep.reserve(sdispls_dep.back());

          for (auto itp = Messages.begin();itp!=Messages.end();itp++) {
            sendbuf_dep.insert(sendbuf_dep.end(),itp->second.begin(),itp->second.end());
          }

        }

        //gather receive sizes
        vector<int> rsizes(this->np,0);
        MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,this->workcomm_);

        //compute receive displacements
        vector<int> rdispls(this->np+1,0);
        rdispls[0] = 0;
        std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);

        //Now do the alltoallv
        vector<SparseTask2D::meta_t> recvbuf(rdispls.back());
        MPI_Alltoallv(&sendbuf[0],&ssizes[0],&sdispls[0],type,&recvbuf[0],&rsizes[0],&rdispls[0],type,this->workcomm_);

        //clear the send buffer
        {
          vector<SparseTask2D::meta_t> tmp;
          sendbuf.swap( tmp );
        }


        //gather receive sizes
        vector<int> rsizes_dep(this->np,0);
        MPI_Alltoall(&ssizes_dep[0],sizeof(int),MPI_BYTE,&rsizes_dep[0],sizeof(int),MPI_BYTE,this->workcomm_);

        //compute receive displacements
        vector<int> rdispls_dep(this->np+1,0);
        rdispls_dep[0] = 0;
        std::partial_sum(rsizes_dep.begin(),rsizes_dep.end(),&rdispls_dep[1]);

        //Now do the alltoallv
        vector<SparseTask2D::meta_t> recvbuf_dep(rdispls_dep.back());
        MPI_Alltoallv(&sendbuf_dep[0],&ssizes_dep[0],&sdispls_dep[0],type,&recvbuf_dep[0],&rsizes_dep[0],&rdispls_dep[0],type,this->workcomm_);
        MPI_Type_free(&type);

        //clear the send buffer
        {
          vector<SparseTask2D::meta_t> tmp;
          sendbuf_dep.swap( tmp );
        }


        //do a top down traversal of my local task list and create the computational tasks
        this->task_graph_solve.clear();
        this->task_graph_solve.reserve(recvbuf.size());
        task_idx_solve.clear();
        task_idx_solve.reserve(recvbuf.size());

        for (auto it = recvbuf.begin();it!=recvbuf.end();it++) {
          auto & cur_op = (*it);
          auto & src_snode = std::get<0>(cur_op);
          auto & tgt_snode = std::get<1>(cur_op);
          auto & type = std::get<2>(cur_op);
          auto & last_local_ancestor = std::get<3>(cur_op);
          auto & facing_first_row = std::get<4>(cur_op);

          SparseTask2D * ptask = nullptr;
          switch(type) {
            case Factorization::op_type::FUC:
            case Factorization::op_type::BUC:
              {
                ptask = new SparseTask2D;
                size_t tuple_size = sizeof(SparseTask2D::meta_t);
                auto meta = &ptask->_meta;
                meta[0] = *it;
              }
              break;
            default:
              break;
          }

          auto I = src_snode;
          auto J = tgt_snode;

          if (ptask!=nullptr) {
            scheduling::key_t key(std::get<0>(cur_op), std::get<1>(cur_op), 0, std::get<2>(cur_op),this->iam);
            this->task_graph_solve.push_back( std::unique_ptr<SparseTask2D>(ptask) );
            task_idx_solve.push_back( std::make_tuple(key,task_idx_solve.size()));
          }

        }


        //sort task_idx_solve by keys
        std::sort(task_idx_solve.begin(),task_idx_solve.end(),[](std::tuple<scheduling::key_t, std::size_t > &a, std::tuple<scheduling::key_t, std::size_t >&b) { return std::get<0>(a) < std::get<0>(b);});


        MPI_Datatype task_idx_type;
        MPI_Type_contiguous( sizeof(std::tuple<scheduling::key_t, std::size_t >), MPI_BYTE, &task_idx_type );
        MPI_Type_commit(&task_idx_type);
        std::vector< std::tuple<scheduling::key_t, std::size_t > > task_idx_recvbuf;
        std::vector<int> task_idx_rsizes(this->np,0);
        std::vector<int> task_idx_rdispls(this->np+1,0);

        int task_idx_size = task_idx_solve.size();
        MPI_Allgather(&task_idx_size,sizeof(task_idx_size),MPI_BYTE,task_idx_rsizes.data(),sizeof(task_idx_size),MPI_BYTE,this->workcomm_);

        task_idx_rdispls[0] = 0;
        std::partial_sum(task_idx_rsizes.begin(),task_idx_rsizes.end(),&task_idx_rdispls[1]);
        task_idx_recvbuf.resize(task_idx_rdispls.back());

        //now communicate the task_idx arrays: allgatherv
        MPI_Allgatherv(task_idx_solve.data(),task_idx_size,task_idx_type,task_idx_recvbuf.data(),task_idx_rsizes.data(),task_idx_rdispls.data(),task_idx_type,this->workcomm_);
        MPI_Type_free(&task_idx_type);

        auto key_lb_comp = [](const std::tuple<scheduling::key_t, std::size_t > & a, const scheduling::key_t & key) { return std::get<0>(a) < key;};

        std::vector<int> update_right_cnt(this->nsuper+1,0);
        std::vector<int> update_up_cnt(this->nsuper+1,0);
        //now we can process the dependency tasks
        for (auto it = recvbuf_dep.begin();it!=recvbuf_dep.end();it++) {
          auto & cur_op = (*it);
          auto & src_snode = std::get<0>(cur_op);
          auto & tgt_snode = std::get<1>(cur_op);
          auto & type = std::get<2>(cur_op);
          auto & last_local_ancestor = std::get<3>(cur_op);
          auto & facing_first_row = std::get<4>(cur_op);

          auto J = src_snode;
          auto I = tgt_snode;

          switch(type) {
            case Factorization::op_type::FUC_DIAG_SEND:
              {
                Idx K = this->SupMembership_[facing_first_row-1];

                auto k1 = scheduling::key_t(I,I,0,Factorization::op_type::FUC,this->iam);
                auto outtask_idx_it = std::lower_bound(task_idx_solve.begin(), task_idx_solve.end(), k1, key_lb_comp);
                auto taskptr = task_graph_solve[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr);
                auto k2 = scheduling::key_t(std::get<0>(taskptr->_meta), std::get<1>(taskptr->_meta), 0, std::get<2>(taskptr->_meta),this->iam);
                bassert(k2 == k1 );

                auto owner = pQueryCELL(J-1,K-1)->owner;
                auto task_idx_beg = &task_idx_recvbuf[task_idx_rdispls[owner]];
                auto task_idx_end = task_idx_beg + task_idx_rsizes[owner];

                std::size_t remote_task_idx = 0;
                scheduling::key_t key(J,K,0,Factorization::op_type::FUC,owner);
                auto task_idx_it = std::lower_bound(task_idx_beg, task_idx_end, (key), key_lb_comp);
                bassert(task_idx_it != task_idx_end); 
                remote_task_idx = std::get<1>(*task_idx_it);
#ifdef _VERBOSE_
                logfileptr->OFS()<<"FUC_DIAG_SEND"<<" cell ("<<I<<","<<I<<") to cell("<<J<<","<<K<<")"<<std::endl;
#endif

                taskptr->out_dependencies[owner].push_back(remote_task_idx);
              }
              break;
            case Factorization::op_type::FUC_DIAG_RECV:
              {
                Idx K = this->SupMembership_[facing_first_row-1];
                std::swap(K,J);
                //find the FUC and add cell I,I as incoming dependency
                auto k1 = scheduling::key_t(J,K,I,Factorization::op_type::FUC,this->iam);
                auto outtask_idx_it = std::lower_bound(task_idx_solve.begin(), task_idx_solve.end(), k1, key_lb_comp);
                auto taskptr = task_graph_solve[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr);

                if (  pQueryCELL(I-1,I-1)->owner != pQueryCELL(J-1,K-1)->owner )
                  taskptr->in_remote_dependencies_cnt++;
                else
                  taskptr->in_local_dependencies_cnt++;

#ifdef _VERBOSE_
                logfileptr->OFS()<<"FUC_DIAG_RECV"<<" cell ("<<I<<","<<I<<") to cell("<<J<<","<<K<<")"<<std::endl;
#endif
              }
              break;

            case Factorization::op_type::FUC_D_SEND:
              {

                auto k1 = scheduling::key_t(I,I,0,Factorization::op_type::FUC,this->iam);
                auto outtask_idx_it = std::lower_bound(task_idx_solve.begin(), task_idx_solve.end(), k1, key_lb_comp);
                auto taskptr = task_graph_solve[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr);
#ifdef _VERBOSE_
                logfileptr->OFS()<<"FUC_D_SEND"<<" cell ("<<I<<","<<I<<") to cell("<<J<<","<<I<<")"<<std::endl;
#endif
                auto k2 = scheduling::key_t(std::get<0>(taskptr->_meta), std::get<1>(taskptr->_meta), 0, std::get<2>(taskptr->_meta),this->iam);
                bassert(k2 == k1 );
                //get index of task FUC(J,I)
                auto owner = pQueryCELL(J-1,I-1)->owner;
                scheduling::key_t key(J,I,0,Factorization::op_type::FUC,owner);
                auto task_idx_beg = &task_idx_recvbuf[task_idx_rdispls[owner]];
                auto task_idx_end = task_idx_beg + task_idx_rsizes[owner];
                auto task_idx_it = std::lower_bound(task_idx_beg, task_idx_end, key, key_lb_comp);
                bassert(task_idx_it != task_idx_end); 
                std::size_t remote_task_idx = std::get<1>(*task_idx_it);

                taskptr->out_dependencies[owner].push_back(remote_task_idx);
              }
              break;
            case Factorization::op_type::FUC_D_RECV:
              {

                //find the FUC and add cell I,I as incoming dependency
                auto k1 = scheduling::key_t(J,I,0,Factorization::op_type::FUC,this->iam);

                auto outtask_idx_it = std::lower_bound(task_idx_solve.begin(), task_idx_solve.end(), k1, key_lb_comp);
                auto taskptr = task_graph_solve[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr); 
                auto k2 = scheduling::key_t(std::get<0>(taskptr->_meta), std::get<1>(taskptr->_meta), 0, std::get<2>(taskptr->_meta),this->iam);
                bassert(k2 == k1 );
                if (  pQueryCELL(J-1,I-1)->owner != pQueryCELL(I-1,I-1)->owner )
                  taskptr->in_remote_dependencies_cnt++;
                else
                  taskptr->in_local_dependencies_cnt++;
#ifdef _VERBOSE_
                logfileptr->OFS()<<"FUC_D_RECV"<<" cell ("<<I<<","<<I<<") to cell("<<J<<","<<I<<")"<<std::endl;
#endif
              }
              break;

            case Factorization::op_type::FUC_RECV:
              {
                //TODO this might be useless
                //find the FUC and add cell J,I as incoming dependency

                auto k1 = scheduling::key_t(J,J,0,Factorization::op_type::FUC,pQueryCELL(J-1,J-1)->owner);
                //TODO keep going from here
                Int parent = supETree.Parent(I-1);
                if (parent == 0 ) {
                  bassert(J==I);
                  k1 = scheduling::key_t(J,J,0,Factorization::op_type::BUC,pQueryCELL(J-1,J-1)->owner);
                }



                auto outtask_idx_it = std::lower_bound(task_idx_solve.begin(), task_idx_solve.end(), k1, key_lb_comp);
                auto taskptr = task_graph_solve[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr);
                if (  pQueryCELL(J-1,I-1)->owner != pQueryCELL(J-1,J-1)->owner )
                  taskptr->in_remote_dependencies_cnt++;
                else
                  taskptr->in_local_dependencies_cnt++;

#ifdef _VERBOSE_
                logfileptr->OFS()<<"FUC_RECV"<<" cell ("<<J<<","<<I<<") to cell("<<J<<","<<J<<") "<<(parent==0?"BUC":"FUC")<<std::endl;
#endif
              }
              break;
            case Factorization::op_type::FUC_SEND:
              {
                //find the FUC and add cell J,J as outgoing dependency
                auto k1 = scheduling::key_t(J,I,0,Factorization::op_type::FUC,pQueryCELL(J-1,I-1)->owner);
                auto outtask_idx_it = std::lower_bound(task_idx_solve.begin(), task_idx_solve.end(), k1, key_lb_comp);
                auto taskptr = task_graph_solve[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr); 
                auto k2 = scheduling::key_t(std::get<0>(taskptr->_meta), std::get<1>(taskptr->_meta), 0, std::get<2>(taskptr->_meta),this->iam);
                bassert(k2 == k1 );

                auto owner = pQueryCELL(J-1,J-1)->owner;
                auto task_idx_beg = &task_idx_recvbuf[task_idx_rdispls[owner]];
                auto task_idx_end = task_idx_beg + task_idx_rsizes[owner];

                std::size_t remote_task_idx = 0;
                scheduling::key_t key(J,J,0,Factorization::op_type::FUC,owner);
                Int parent = supETree.Parent(I-1);
                if (parent == 0 ) {
                  key = scheduling::key_t(J,J,0,Factorization::op_type::BUC,owner);
                }
                auto task_idx_it = std::lower_bound(task_idx_beg, task_idx_end, (key), key_lb_comp);
                bassert(task_idx_it != task_idx_end); 
                remote_task_idx = std::get<1>(*task_idx_it);

                taskptr->out_dependencies[owner].push_back(remote_task_idx);

                //check if it's not in it first
                //this count should also be sent along the rpc
                update_right_cnt[J]++;

#ifdef _VERBOSE_
                logfileptr->OFS()<<"FUC_SEND"<<" cell ("<<J<<","<<I<<") to cell("<<J<<","<<J<<") "<<(parent==0?"BUC":"FUC")<<std::endl;
#endif
              }
              break;


            case Factorization::op_type::BUC_D_RECV:
              {

                //find the BUC and add cell J,I as incoming dependency
                auto k1 = scheduling::key_t(I,I,0,Factorization::op_type::BUC,pQueryCELL(I-1,I-1)->owner);
                auto outtask_idx_it = std::lower_bound(task_idx_solve.begin(), task_idx_solve.end(), k1, key_lb_comp);
                auto taskptr = task_graph_solve[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr);
                auto k2 = scheduling::key_t(std::get<0>(taskptr->_meta), std::get<1>(taskptr->_meta), 0, std::get<2>(taskptr->_meta),this->iam);
                bassert(k2 == k1 );

                if (  pQueryCELL(J-1,I-1)->owner != pQueryCELL(I-1,I-1)->owner )
                  taskptr->in_remote_dependencies_cnt++;
                else
                  taskptr->in_local_dependencies_cnt++;
#ifdef _VERBOSE_
                logfileptr->OFS()<<"BUC_D_RECV"<<" cell ("<<J<<","<<I<<") to cell("<<I<<","<<I<<")"<<std::endl;
#endif

              }
              break;
            case Factorization::op_type::BUC_D_SEND:
              {

                //find the BUC and add cell I,I as outgoing dependency
                auto k1 = scheduling::key_t(J,I,0,Factorization::op_type::BUC,pQueryCELL(J-1,I-1)->owner);
                auto outtask_idx_it = std::lower_bound(task_idx_solve.begin(), task_idx_solve.end(), k1, key_lb_comp);
                auto taskptr = task_graph_solve[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr); 
                auto k2 = scheduling::key_t(std::get<0>(taskptr->_meta), std::get<1>(taskptr->_meta), 0, std::get<2>(taskptr->_meta),this->iam);
                bassert(k2 == k1 );
                auto owner = pQueryCELL(I-1,I-1)->owner;
                //get index of task BUC(I,I)
                scheduling::key_t key(I,I,0,Factorization::op_type::BUC,owner);
                auto task_idx_beg = &task_idx_recvbuf[task_idx_rdispls[owner]];
                auto task_idx_end = task_idx_beg + task_idx_rsizes[owner];
                auto task_idx_it = std::lower_bound(task_idx_beg, task_idx_end, key, key_lb_comp);
                bassert(task_idx_it != task_idx_end); 
                std::size_t remote_task_idx = std::get<1>(*task_idx_it);

                taskptr->out_dependencies[owner].push_back(remote_task_idx);
                update_up_cnt[I]++;
#ifdef _VERBOSE_
                logfileptr->OFS()<<"BUC_D_SEND"<<" cell ("<<J<<","<<I<<") to cell("<<I<<","<<I<<")"<<std::endl;
#endif
              }
              break;


            case Factorization::op_type::BUC_SEND:
              {
                //TODO this might be useless
#ifdef _VERBOSE_
                logfileptr->OFS()<<"BUC_SEND"<<" cell ("<<J<<","<<J<<") to cell ("<<J<<","<<I<<")"<<std::endl;
#endif
                //find the BUC and add cell J,I as outgoing dependency
                auto k1 = scheduling::key_t(J,J,0,Factorization::op_type::BUC,pQueryCELL(J-1,J-1)->owner);
                auto outtask_idx_it = std::lower_bound(task_idx_solve.begin(), task_idx_solve.end(), k1, key_lb_comp);
                auto taskptr = task_graph_solve[std::get<1>(*outtask_idx_it)].get();

                auto k2 = scheduling::key_t(std::get<0>(taskptr->_meta), std::get<1>(taskptr->_meta), 0, std::get<2>(taskptr->_meta),this->iam);
                bassert( k1==k2);

                bassert(taskptr!=nullptr);

                auto owner = pQueryCELL(J-1,I-1)->owner;
                auto task_idx_beg = &task_idx_recvbuf[task_idx_rdispls[owner]];
                auto task_idx_end = task_idx_beg + task_idx_rsizes[owner];

                std::size_t remote_task_idx = 0;
                scheduling::key_t key(J,I,0,Factorization::op_type::BUC,owner);
                auto task_idx_it = std::lower_bound(task_idx_beg, task_idx_end, (key), key_lb_comp);
                bassert(task_idx_it != task_idx_end); 
                remote_task_idx = std::get<1>(*task_idx_it);

                taskptr->out_dependencies[owner].push_back(remote_task_idx);
              }
              break;
              //WE SHOULD BE RECEIVING ONLY ONCE PER REMOTE PROCESS
            case Factorization::op_type::BUC_RECV:
              {
#ifdef _VERBOSE_
                logfileptr->OFS()<<"BUC_RECV"<<" cell ("<<J<<","<<J<<") to cell ("<<J<<","<<I<<")"<<std::endl;
#endif
                //find the BUC and add cell J,J as incoming dependency
                auto k1 = scheduling::key_t(J,I,0,Factorization::op_type::BUC,pQueryCELL(J-1,I-1)->owner);
                auto outtask_idx_it = std::lower_bound(task_idx_solve.begin(), task_idx_solve.end(), k1, key_lb_comp);
                auto taskptr = task_graph_solve[std::get<1>(*outtask_idx_it)].get();
                bassert(taskptr!=nullptr); 
                auto k2 = scheduling::key_t(std::get<0>(taskptr->_meta), std::get<1>(taskptr->_meta), 0, std::get<2>(taskptr->_meta),this->iam);
                bassert(k2 == k1 );
                bassert( std::get<2>(taskptr->_meta)==Factorization::op_type::BUC);

                if (  pQueryCELL(J-1,I-1)->owner != pQueryCELL(J-1,J-1)->owner )
                  taskptr->in_remote_dependencies_cnt++;
                else
                  taskptr->in_local_dependencies_cnt++;
              }
              break;

            default:
              break;
          }
        }


#ifdef _VERBOSE_
        logfileptr->OFS()<<update_right_cnt<<std::endl;
        logfileptr->OFS()<<update_up_cnt<<std::endl;
#endif

        //Now we have our local part of the task graph
        for (auto it = this->task_graph_solve.begin(); it != this->task_graph_solve.end(); it++) {
          auto & ptask = *it;
          auto meta = &ptask->_meta;

          auto & src_snode = std::get<0>(meta[0]);
          auto & tgt_snode = std::get<1>(meta[0]);
          auto & last_local_ancestor = std::get<3>(meta[0]);
          auto & facing_row = std::get<4>(meta[0]);
          auto & type = std::get<2>(meta[0]);

          auto J = src_snode;
          auto I = tgt_snode;

          auto remote_deps = ptask->in_remote_dependencies_cnt;
          auto local_deps = ptask->in_local_dependencies_cnt;

#ifdef _USE_PROM_AVAIL_
          ptask->in_avail_prom.require_anonymous(remote_deps);
          ptask->in_avail_counter = remote_deps;
#else
          ptask->in_avail_counter = remote_deps;
#endif
#ifdef _USE_PROM_RDY_
          ptask->in_prom.require_anonymous(local_deps + remote_deps);
          ptask->in_counter = local_deps + remote_deps;
#else
          ptask->in_counter = local_deps + remote_deps;
#endif

          auto ptr = ptask.get();
#ifdef _VERBOSE_
          switch(type) {
            case Factorization::op_type::FUC:
              {
                logfileptr->OFS()<<"FUC"<<" from "<<I<<" to ("<<I<<") cell ("<<J<<","<<I<<") "<<remote_deps<<std::endl;
              }
              break;
            case Factorization::op_type::BUC:
              {
                logfileptr->OFS()<<"BUC"<<" from "<<I<<" to ("<<I<<") cell ("<<J<<","<<I<<") "<<remote_deps<<std::endl;
              }
              break;
            default:
              break;
          }
#endif

          switch(type) {
            case Factorization::op_type::FUC:
              {
                int dep_cnt = update_right_cnt[J];
                int upd_diag_cnt = 0;
                ptask->execute = [this,ptr,I,J,dep_cnt,upd_diag_cnt] () {
                  scope_timer(b,SOLVE_FUC_TASK);
                  auto ptask = ptr;
                  auto ptr_cell = pQueryCELL2(J-1,I-1);
                  auto & update_right_cnt = this->solve_data.update_right_cnt;
                  auto & contribs = this->solve_data.contribs;
                  auto rhs = this->solve_data.rhs;
     		  auto nrhs = this->solve_data.nrhs;
                  snodeBlock_sptr_t ptr_contrib = nullptr;

                  if ( I == J ) {
                    { 
#ifdef SP_THREADS
                      cell_lock<std::atomic<bool> > lock_cell(this->solve_data.contribs_lock[J]);
#endif
                      {
                        //Allocate Y(K) with the same structure as L(K,M) where M is the highest anscestor in the Etree this rank owns
                        auto & contrib_slot = contribs[J];
                        auto & rptr_contrib = std::get<1>(contrib_slot);
                        if ( ! rptr_contrib ) {
                          bassert(update_right_cnt[J] == 0);
                          rptr_contrib = std::static_pointer_cast<snodeBlock_t>(this->options_.decomposition == DecompositionType::LDL?std::make_shared<snodeBlockLDL_t>():std::make_shared<snodeBlock_t>());
                          auto ptr_test_cell = pQueryCELL(J-1,J-1);
                          bassert( ptr_test_cell->owner == this->iam );
			  //copy ptr_cell test into new rptr_contrib block
                          rptr_contrib->copy_row_structure(nrhs,(snodeBlock_t*)ptr_test_cell.get());
                          for (rowind_t row = 0; row< rptr_contrib->total_rows(); ++row) {
                            rowind_t srcRow = this->Order_.perm[rptr_contrib->first_col-1+row] -1;
                            for (rowind_t col = 0; col<nrhs;++col) {
                              rptr_contrib->_nzval[row*nrhs+col] = rhs[srcRow + col*this->iSize_];
                            }
      			  }
                        }
                        else {
                          bassert( rptr_contrib->nnz() >0 && rptr_contrib->width()>0);
                          //Add data from RHS
                          for (rowind_t row = 0; row< rptr_contrib->total_rows(); ++row) {
                            rowind_t srcRow = this->Order_.perm[rptr_contrib->first_col-1+row] -1;
                            for (rowind_t col = 0; col<nrhs;++col) {
                              rptr_contrib->_nzval[row*nrhs+col] += rhs[srcRow + col*this->iSize_];
                            }
                          }	
                        }
                        ptr_contrib = rptr_contrib;
                      }


                      //Then accumulate everything coming from children
                      for ( auto && msg_ptr: ptask->input_msg ) {
                        if ( msg_ptr->extra_data ) {
                          auto ptr_rem_contrib = (snodeBlock_t*)(msg_ptr->extra_data.get());
                          ptr_contrib->forward_update(ptr_rem_contrib);
                        }
                      }

                      //compute the product Y(J) -= L(J,J) * Y(J)
                      ptr_cell->forward_update_contrib(ptr_contrib.get());
                    } //release the lock here

                    //send the contrib down
                    for (auto it = ptask->out_dependencies.begin(); it!=ptask->out_dependencies.end(); it++) {
                      auto pdest = it->first;
                      auto & tgt_cells = it->second;
                      //serialize data once, and list of meta data
                      //diag contrib is output data so it will not be deleted
                      if ( pdest != this->iam ) {
                        upcxx::rpc_ff( pdest,  
                            [] (int sp_handle, upcxx::global_ptr<char> gptr,
				size_t storage_size, size_t nnz, size_t nblocks, rowind_t width,  SparseTask2D::meta_t meta, 
				upcxx::view<std::size_t> target_cells ) { 
                                //there is a map between sp_handle and task_graphs
                                auto matptr = (symPACKMatrix2D<colptr_t,rowind_t,T> *) g_sp_handle_to_matrix[sp_handle];
                                auto I = std::get<1>(meta);
                                rowind_t fc = matptr->Xsuper_[I-1];

                                //store pointer & associated metadata somewhere
                                std::shared_ptr <SparseTask2D::data_t > data;

                                for ( auto & tgt_cell_idx: target_cells) {
                                auto taskptr = matptr->task_graph_solve[tgt_cell_idx].get();
                                {
                                if ( ! data ) {
                                data = std::make_shared<SparseTask2D::data_t >();
                                data->in_meta = meta;
                                data->size = storage_size;
                                data->remote_gptr = gptr;
                                }

                                taskptr->input_msg.push_back(data);
                                data->target_tasks.push_back(taskptr);
                                taskptr->avail_dep(1,matptr->scheduler);
                                }
                                }


                                if ( data ) {
                                  data->on_fetch_future = data->on_fetch_future.then(
                                      [fc,width,nnz,nblocks,I,matptr](SparseTask2D::data_t * pdata) {
                                      //create snodeBlock_t and store it in the extra_data
                                      if (matptr->options_.decomposition == DecompositionType::LDL) {
                                        pdata->extra_data = std::shared_ptr<blockCellBase_t>( (blockCellBase_t*)new snodeBlockLDL_t(I,I,pdata->landing_zone,fc,width,nnz,nblocks,0) );
                                      }
                                      else {
                                        pdata->extra_data = std::shared_ptr<blockCellBase_t>( (blockCellBase_t*)new snodeBlock_t(I,I,pdata->landing_zone,fc,width,nnz,nblocks) );
                                      }
                                      return upcxx::to_future(pdata);
                                      });
                                  //TODO check this
#ifdef _EAGER_FETCH_
                                  data->allocate();
                                  data->fetch();
#endif
                                }


                            }, this->sp_handle, ptr_contrib->_gstorage, 
			       ptr_contrib->_storage_size, ptr_contrib->nnz(), ptr_contrib->nblocks(), 
			       ptr_contrib->width(), ptask->_meta, upcxx::make_view(tgt_cells.begin(),tgt_cells.end())); 
                      }
                      else {
                        for ( auto & tgt_cell_idx: tgt_cells ) {
                          auto taskptr = task_graph_solve[tgt_cell_idx].get();
                          bassert(taskptr!=nullptr); 
                          bassert(std::get<2>(taskptr->_meta)==Factorization::op_type::FUC
                              || std::get<2>(taskptr->_meta)==Factorization::op_type::BUC);
                          //mark the dependency as satisfied
                          taskptr->satisfy_dep(1,this->scheduler);
                        }
                      }
                    }
                  }
                  else {
                    bool deleteContrib = false;
                    {
#ifdef SP_THREADS
                      cell_lock<std::atomic<bool> > lock_cell(this->solve_data.contribs_lock[J]);
#endif
                      {
                        auto & contrib_slot = contribs[J];
                        auto & counter = std::get<0>(contrib_slot);

                        auto & rptr_contrib = std::get<1>(contrib_slot);
                        if ( ! rptr_contrib ) {
                          bassert(update_right_cnt[J] == 0);

                          rptr_contrib = std::static_pointer_cast<snodeBlock_t>(this->options_.decomposition == DecompositionType::LDL?std::make_shared<snodeBlockLDL_t>():std::make_shared<snodeBlock_t>());

                          rowind_t nrows = this->Xsuper_[J] - this->Xsuper_[J-1];
                          
			  rptr_contrib = std::static_pointer_cast<snodeBlock_t>(this->options_.decomposition == DecompositionType::LDL?std::make_shared<snodeBlockLDL_t>(J,J,this->Xsuper_[J-1],nrhs,nrows*nrhs,1,nrows,true):std::make_shared<snodeBlock_t>(J,J,this->Xsuper_[J-1],nrhs,nrows*nrhs,1));
                          
			  //TODO THIS IS TERRIBLE!!!
			  rptr_contrib->add_block(this->Xsuper_[J-1],nrows);
			  //TODO: add cublas scal here
                          std::fill(rptr_contrib->_nzval,rptr_contrib->_nzval+rptr_contrib->_nnz,T(0));
     			}


                        if (dep_cnt>0) {
                          auto ptr_test_cell = pQueryCELL(J-1,J-1);
                          if ( ptr_test_cell->owner != this->iam ) {
                            deleteContrib = true;

                            if (update_right_cnt[J] == 0 ) {
                              //TODO THIS IS NOT THREADSAFE, CAN THIS BE DONE BEFORE?
                              ((*this->solve_data.remoteDeallocCounter))+=dep_cnt;
                              counter+=dep_cnt;
                            }
                            else {
                            }
                          }
                        }


                        ptr_contrib = rptr_contrib;
                      }



                      bassert(ptask->input_msg.size() <= 1 );
                      //unpack the part of the solution (y) facing the diagonal block of the supernode
                      auto ptr_diagContrib = pQueryCELL(I-1,I-1); 
                      if ( ptr_diagContrib->owner == this->iam ) {
                        auto tmp2 = contribs.at(I);
                        auto tmp = contribs[I];
                        ptr_diagContrib = std::get<1>(contribs[I]);
                        bassert(ptr_diagContrib != nullptr);
                      }
                      else {
                        if ( ptask->input_msg.size() > 0 ) {
                          auto msg_ptr = *ptask->input_msg.begin();
                          bassert ( msg_ptr->extra_data != nullptr );
                          ptr_diagContrib = (msg_ptr->extra_data);
                        }
                      }

                      bassert(I<J);
                      //compute the product Y(J) -= L(J,I) * Y(I)
                      ptr_cell->forward_update_contrib(ptr_contrib.get(),ptr_diagContrib.get());
                    } //release the lock here

                    //send the contrib
                    for (auto it = ptask->out_dependencies.begin(); it!=ptask->out_dependencies.end(); it++) {
                      auto pdest = it->first;
                      auto & tgt_cells = it->second;
                      //serialize data once, and list of meta data
                      //factor is output data so it will not be deleted
                      if ( pdest != this->iam ) {
                        update_right_cnt[J]++;
                        if ( dep_cnt == update_right_cnt[J]) {
                          upcxx::rpc_ff( pdest, 
                              [dep_cnt,deleteContrib ] (int sp_handle, upcxx::global_ptr<char> gptr, 
				      			size_t storage_size, size_t nnz, size_t nblocks, rowind_t width, SparseTask2D::meta_t meta, upcxx::view<std::size_t> target_cells ) { 
                                  //there is a map between sp_handle and task_graphs
                                  auto matptr = (symPACKMatrix2D<colptr_t,rowind_t,T> *) g_sp_handle_to_matrix[sp_handle];
                                  auto I = std::get<1>(meta);
                                  auto J = std::get<0>(meta);


                                  rowind_t fc = matptr->Xsuper_[I-1];

                                  //store pointer & associated metadata somewhere
                                  std::shared_ptr <SparseTask2D::data_t > data;

                                  for ( auto & tgt_cell_idx: target_cells) {
                                  auto taskptr = matptr->task_graph_solve[tgt_cell_idx].get();
                                  if ( ! data ) {
                                  data = std::make_shared<SparseTask2D::data_t >();
                                  data->in_meta = meta;
                                  data->size = storage_size;
                                  data->remote_gptr = gptr;
				  }


                                  taskptr->input_msg.push_back(data);
                                  data->target_tasks.push_back(taskptr);
                                  taskptr->avail_dep(dep_cnt,matptr->scheduler);

                                  if (dep_cnt>1) {
                                    taskptr->satisfy_dep(dep_cnt-1,matptr->scheduler);
                                  }
                                  }


                                  if ( data ) {
                                    int owner = gptr.where();
                                    data->on_fetch_future = data->on_fetch_future.then(
                                        [fc,width,nnz,nblocks,I,matptr,sp_handle,dep_cnt,owner,J,deleteContrib](SparseTask2D::data_t * pdata) {
                                        //create snodeBlock_t and store it in the extra_data
                                        if (matptr->options_.decomposition == DecompositionType::LDL) {
                                        pdata->extra_data = std::shared_ptr<blockCellBase_t>( (blockCellBase_t*)new snodeBlockLDL_t(I,I,pdata->landing_zone,fc,width,nnz,nblocks,0) );
                                        }
                                        else {
					pdata->extra_data = std::shared_ptr<blockCellBase_t>( (blockCellBase_t*)new snodeBlock_t(I,I,pdata->landing_zone,fc,width,nnz,nblocks) );
					}

                                        //send a rpc_ff on the owner of the data to signal we have fetched it
                                        if ( deleteContrib ) {
                                        matptr->solve_data.deallocRemote(owner,sp_handle,J,dep_cnt);
                                        }

                                        return upcxx::to_future(pdata);
                                        });
                                    //TODO check this
#ifdef _EAGER_FETCH_
                                    data->allocate();
                                    data->fetch();
#endif
                                  }

                              }, this->sp_handle, ptr_contrib->_gstorage,
			 	 ptr_contrib->_storage_size, ptr_contrib->nnz(), ptr_contrib->nblocks(), ptr_contrib->width() ,ptask->_meta, upcxx::make_view(tgt_cells.begin(),tgt_cells.end())); 
                        }
                      }
                      else {
                        for ( auto & tgt_cell_idx: tgt_cells ) {
                          auto taskptr = task_graph_solve[tgt_cell_idx].get();
                          bassert(taskptr!=nullptr); 
                          bassert(std::get<2>(taskptr->_meta)==Factorization::op_type::FUC);
                          //mark the dependency as satisfied
                          taskptr->satisfy_dep(1,this->scheduler);
                        }
                      }
                    }
                  }
                  ptask->executed = true;

                };

              }
              break;
            case Factorization::op_type::BUC:
              {
                int dep_cnt = update_up_cnt[I];
                ptask->execute = [this,ptr,I,J,dep_cnt] () {
                  scope_timer(b,SOLVE_BUC_TASK);
                  auto ptask = ptr;
                  auto ptr_cell = pQueryCELL2(J-1,I-1);

                  auto & update_up_cnt = this->solve_data.update_up_cnt;
                  auto & contribs = this->solve_data.contribs;
                  auto rhs = this->solve_data.rhs;
                  auto nrhs = this->solve_data.nrhs;
                  auto ptr_tgtcell = pQueryCELL(I-1,I-1);

                  snodeBlock_sptr_t ptr_contrib = nullptr;
                  if ( I == J ) {
                    {
#ifdef SP_THREADS
                      cell_lock<std::atomic<bool> > lock_cell(this->solve_data.contribs_lock[I]);
#endif
                      {
                        bassert ( ptr_tgtcell->owner == this->iam ) ;
                        auto & contrib_slot = contribs[I];
                        auto & rptr_contrib = std::get<1>(contrib_slot);
                        bassert( rptr_contrib != nullptr );
                        bassert( rptr_contrib->nnz() >0 && rptr_contrib->width()>0);
                        ptr_contrib = rptr_contrib;
                      }  
                      if (this->options_.decomposition == DecompositionType::LDL) {
                        if ( ptr_tgtcell->owner == this->iam ) {
                          auto & tgt_ldlcell = *std::dynamic_pointer_cast<snodeBlockLDL_t>(ptr_contrib);
                          if(!tgt_ldlcell.scaled){
                            ((snodeBlockLDL_t*)ptr_tgtcell.get())->scale_contrib(&tgt_ldlcell);
                          }
                        }
                      }
                      //Then accumulate everything coming from children
                      for ( auto && msg_ptr: ptask->input_msg ) {
                        if ( msg_ptr->extra_data ) {
                          auto ptr_rem_contrib = (snodeBlock_t*)(msg_ptr->extra_data.get());
                          ptr_contrib->back_update(ptr_rem_contrib);
                        }
                      }

                      //this is a diagonal block update
                      //this could be implemented as a straight trsm
                      //compute the product Y(I) -= L(I,I)^T * Y(I)
                      ptr_cell->back_update_contrib(ptr_contrib.get());


#ifndef NDEBUG
                      if (this->options_.decomposition == DecompositionType::LDL) {
                        auto tgt_ldlcell = (snodeBlockLDL_t*)(ptr_contrib.get());
                        bassert(tgt_ldlcell->scaled);
                      }
#endif
                    }
                    //send the contrib left
                    for (auto it = ptask->out_dependencies.begin(); it!=ptask->out_dependencies.end(); it++) {
                      auto pdest = it->first;
                      auto & tgt_cells = it->second;
                      //serialize data once, and list of meta data
                      //diag contrib is output data so it will not be deleted
                      if ( pdest != this->iam ) {
                        upcxx::rpc_ff( pdest,  
                            [ ] (int sp_handle, upcxx::global_ptr<char> gptr, 
				 size_t storage_size, size_t nnz, size_t nblocks, rowind_t width, SparseTask2D::meta_t meta, upcxx::view<std::size_t> target_cells ) { 
                                //there is a map between sp_handle and task_graphs
                                auto matptr = (symPACKMatrix2D<colptr_t,rowind_t,T> *) g_sp_handle_to_matrix[sp_handle];
                                auto I = std::get<1>(meta);
                                rowind_t fc = matptr->Xsuper_[I-1];

                                //store pointer & associated metadata somewhere
                                std::shared_ptr <SparseTask2D::data_t > data;

                                for ( auto & tgt_cell_idx: target_cells) {
                                auto taskptr = matptr->task_graph_solve[tgt_cell_idx].get();


                                {
                                if ( ! data ) {
                                data = std::make_shared<SparseTask2D::data_t >();
                                data->in_meta = meta;
                                data->size = storage_size;
                                data->remote_gptr = gptr;
                                }

                                taskptr->input_msg.push_back(data);
                                data->target_tasks.push_back(taskptr);
                                taskptr->avail_dep(1,matptr->scheduler);
                                }
                                }


                                if ( data ) {
                                  data->on_fetch_future = data->on_fetch_future.then(
                                      [fc,width,nnz,nblocks,I,matptr](SparseTask2D::data_t * pdata) {
                                      //create snodeBlock_t and store it in the extra_data
                                      if (matptr->options_.decomposition == DecompositionType::LDL) {
                                      pdata->extra_data = std::shared_ptr<blockCellBase_t>( (blockCellBase_t*)new snodeBlockLDL_t(I,I,pdata->landing_zone,fc,width,nnz,nblocks,0) );
                                      }
                                      else {
                                      pdata->extra_data = std::shared_ptr<blockCellBase_t>( (blockCellBase_t*)new snodeBlock_t(I,I,pdata->landing_zone,fc,width,nnz,nblocks) );
                                      }
                                      return upcxx::to_future(pdata); 
                                      });
                                  //TODO check this
#ifdef _EAGER_FETCH_
                                  data->allocate();
                                  data->fetch();
#endif
                                }


                            },  this->sp_handle, ptr_contrib->_gstorage,
			        ptr_contrib->_storage_size, ptr_contrib->nnz(), ptr_contrib->nblocks(), ptr_contrib->width(), ptask->_meta, upcxx::make_view(tgt_cells.begin(),tgt_cells.end())); 
                      }
                      else {
                        for ( auto & tgt_cell_idx: tgt_cells ) {
                          auto taskptr = task_graph_solve[tgt_cell_idx].get();
                          bassert(taskptr!=nullptr); 
                          bassert(std::get<2>(taskptr->_meta)==Factorization::op_type::BUC);
                          //mark the dependency as satisfied
                          taskptr->satisfy_dep(1,this->scheduler);
                        }
                      }
                    }
                  }
                  else {
                    bool deleteContrib = false;
                    {
#ifdef SP_THREADS
                      cell_lock<std::atomic<bool> > lock_cell(this->solve_data.contribs_lock[I]);
#endif
                      {
                        {
                          auto & contrib_slot = contribs[I];
                          auto & counter = std::get<0>(contrib_slot);
                          auto & rptr_contrib = std::get<1>(contrib_slot);


                          if ( ptr_tgtcell->owner != this->iam ) {
                            if (update_up_cnt[I] == 0 ) {
                              if ( !rptr_contrib ) {
                                rowind_t nrows = this->Xsuper_[I] - this->Xsuper_[I-1];
                                rptr_contrib = std::static_pointer_cast<snodeBlock_t>(this->options_.decomposition == DecompositionType::LDL?std::make_shared<snodeBlockLDL_t>(I,I,this->Xsuper_[I-1],nrhs,nrows*nrhs,1,nrows,true):std::make_shared<snodeBlock_t>(I,I,this->Xsuper_[I-1],nrhs,nrows*nrhs,1));
                                //TODO THIS IS TERRIBLE!!!
                                rptr_contrib->add_block(this->Xsuper_[I-1],nrows);

                              }
                              std::fill(rptr_contrib->_nzval,rptr_contrib->_nzval+rptr_contrib->_nnz,T(0));

                            }
                            else {
                              bassert(rptr_contrib!=nullptr);
                            }

                            if (dep_cnt>0) {
                              deleteContrib = true;

                              if (deleteContrib && update_up_cnt[I] == 0 ) {
                                //TODO THIS IS NOT THREADSAFE, CAN THIS BE DONE BEFORE?
                                ((*this->solve_data.remoteDeallocCounter))+=dep_cnt;
                                counter+=dep_cnt;
                              }
                            }


                          }
                          else {
                            bassert(rptr_contrib!=nullptr);
                          }

                          bassert( rptr_contrib->nnz() >0 && rptr_contrib->width()>0);
                          ptr_contrib = rptr_contrib;
                        }

                        if (this->options_.decomposition == DecompositionType::LDL) {
                          if ( ptr_tgtcell->owner == this->iam ) {
                            auto & tgt_ldlcell = *std::dynamic_pointer_cast<snodeBlockLDL_t>(ptr_contrib);
                            if(!tgt_ldlcell.scaled){
                              ((snodeBlockLDL_t*)ptr_tgtcell.get())->scale_contrib(&tgt_ldlcell);
                            }
                          }
                        }


                        bassert(ptask->input_msg.size() <= 1 );
                        //unpack the part of the solution (y) facing the diagonal block of the supernode
                        auto ptr_diagContrib = pQueryCELL(J-1,J-1); 
                        if ( ptr_diagContrib->owner == this->iam ) {
                          ptr_diagContrib = std::get<1>(contribs[J]);
                          bassert(ptr_diagContrib != nullptr);

                          if (this->options_.decomposition == DecompositionType::LDL) {
                            auto tgt_ldlcell = (snodeBlockLDL_t*)(ptr_diagContrib.get());
                            bassert(tgt_ldlcell->scaled);
                          }


                        }
                        else {
                          if ( ptask->input_msg.size() > 0 ) {
                            auto msg_ptr = *ptask->input_msg.begin();
                            bassert ( msg_ptr->extra_data != nullptr );
                            ptr_diagContrib = (msg_ptr->extra_data);
                          }
                        }



                        //compute the product Y(I) -= L(I,J)^T * Y(J)
                        ptr_cell->back_update_contrib(ptr_contrib.get(),ptr_diagContrib.get());
                      }
                    }
                    //send the contrib up
                    for (auto it = ptask->out_dependencies.begin(); it!=ptask->out_dependencies.end(); it++) {
                      auto pdest = it->first;
                      auto & tgt_cells = it->second;
                      //serialize data once, and list of meta data
                      //factor is output data so it will not be deleted
                      if ( pdest != this->iam ) {
                        update_up_cnt[I]++;
                        if ( dep_cnt == update_up_cnt[I]) {
                          upcxx::rpc_ff( pdest, 
                              [ dep_cnt,deleteContrib ] (int sp_handle, upcxx::global_ptr<char> gptr, 
				      			size_t storage_size, size_t nnz, size_t nblocks, rowind_t width, SparseTask2D::meta_t meta, upcxx::view<std::size_t> target_cells ) { 
                                  //there is a map between sp_handle and task_graphs
                                  auto matptr = (symPACKMatrix2D<colptr_t,rowind_t,T> *) g_sp_handle_to_matrix[sp_handle];
                                  auto I = std::get<1>(meta);
                                  rowind_t fc = matptr->Xsuper_[I-1];

                                  //store pointer & associated metadata somewhere
                                  std::shared_ptr <SparseTask2D::data_t > data;

                                  for ( auto & tgt_cell_idx: target_cells) {
                                  auto taskptr = matptr->task_graph_solve[tgt_cell_idx].get();
                                  if ( ! data ) {
                                  data = std::make_shared<SparseTask2D::data_t >();
                                  data->in_meta = meta;
                                  data->size = storage_size;
                                  data->remote_gptr = gptr;
                                  }


                                  taskptr->input_msg.push_back(data);
                                  data->target_tasks.push_back(taskptr);
                                  taskptr->avail_dep(dep_cnt,matptr->scheduler);

                                  if (dep_cnt>1) {
                                    taskptr->satisfy_dep(dep_cnt-1,matptr->scheduler);
                                  }
                                  }


                                  if ( data ) {
                                    int owner = gptr.where();
                                    data->on_fetch_future = data->on_fetch_future.then(
                                        [fc,width,nnz,nblocks,I,matptr,owner,sp_handle,dep_cnt,deleteContrib](SparseTask2D::data_t * pdata) {
                                        //create snodeBlock_t and store it in the extra_data
                                        if (matptr->options_.decomposition == DecompositionType::LDL) {
                                          pdata->extra_data = std::shared_ptr<blockCellBase_t>( (blockCellBase_t*)new snodeBlockLDL_t(I,I,pdata->landing_zone,fc,width,nnz,nblocks,0) );
                                        }
                                        else {
                                          pdata->extra_data = std::shared_ptr<blockCellBase_t>( (blockCellBase_t*)new snodeBlock_t(I,I,pdata->landing_zone,fc,width,nnz,nblocks) );
                                        }


                                        //send a rpc_ff on the owner of the data to signal we have fetched it
                                        if ( deleteContrib ) {
                                        matptr->solve_data.deallocRemote(owner,sp_handle,I,dep_cnt);
                                        }

                                        return upcxx::to_future(pdata);
                                        });
#ifdef _EAGER_FETCH_
                                    data->allocate();
                                    data->fetch();
#endif
                                  }

                              }, this->sp_handle, ptr_contrib->_gstorage,
				ptr_contrib->_storage_size, ptr_contrib->nnz(), ptr_contrib->nblocks(), ptr_contrib->width(), ptask->_meta, upcxx::make_view(tgt_cells.begin(),tgt_cells.end())); 
                        }
                      }
                      else {
                        for ( auto & tgt_cell_idx: tgt_cells ) {
                          auto taskptr = task_graph_solve[tgt_cell_idx].get();
                          bassert(taskptr!=nullptr); 
                          bassert(std::get<2>(taskptr->_meta)==Factorization::op_type::BUC);
                          //mark the dependency as satisfied
                          taskptr->satisfy_dep(1,this->scheduler);
                        }
                      }
                    }
                  }
                  ptask->executed = true;
                };
              }
              break;
            default:
              break;
          }
        }


        logfileptr->OFS()<<"Solve task graph created"<<std::endl;
        
      }


    } 

  template <typename colptr_t, typename rowind_t, typename T, typename int_t>
    void symPACKMatrix2D<colptr_t,rowind_t,T,int_t>::DistributeMatrix(DistSparseMatrix<T> & pMat ) {

      std::map<Int,size_t > send_map;

      typedef std::conditional< sizeof(Idx) < sizeof(Ptr), Idx, Ptr>::type minTypeIdxPtr;
      using minType = typename std::conditional< sizeof(minTypeIdxPtr) < sizeof(T), minTypeIdxPtr, T >::type;

      size_t minSize = sizeof(minType);
      size_t IdxToMin = sizeof(Idx) / minSize;
      size_t PtrToMin = sizeof(Ptr) / minSize;
      size_t TToMin = sizeof(T) / minSize;

      Int baseval = pMat.Localg_.GetBaseval();
      Idx FirstLocalCol = pMat.Localg_.vertexDist[this->iam] + (1 - baseval); //1-based
      Idx LastLocalCol = pMat.Localg_.vertexDist[this->iam+1] + (1 - baseval); //1-based

      {
        scope_timer(a,symPACKMatrix2D::DistributeMatrix::Counting);
        for (Int I=1;I<this->Xsuper_.size();I++) {
          Idx fc = this->Xsuper_[I-1];
          Idx lc = this->Xsuper_[I]-1;
          Int iWidth = lc-fc+1;

          //post all the recv and sends
          for (Idx col = fc;col<=lc;col++) {
            //corresponding column in the unsorted matrix A
            Idx orig_col = this->Order_.perm[col-1];
            if (orig_col>= FirstLocalCol && orig_col < LastLocalCol) {
              Ptr nrows = 0;
              Idx local_col = orig_col - FirstLocalCol+1;//1 based
              for (Ptr rowidx = pMat.Localg_.colptr[local_col-1] + (1-baseval); rowidx<pMat.Localg_.colptr[local_col]+(1-baseval); ++rowidx) {
                Idx orig_row = pMat.Localg_.rowind[rowidx-1]+(1-baseval);//1-based
                Idx row = this->Order_.invp[orig_row-1];

                Int J = this->SupMembership_[row-1];

                if (row>=col) {
                  //this has to go in cell(J,I)
                  auto ptr_tgt_cell = pQueryCELL(J-1,I-1);
                  Int iDest = ptr_tgt_cell->owner;
                  send_map[iDest] += IdxToMin+PtrToMin+ (IdxToMin + TToMin);
                }
                else {
                  //this has to go in cell(I,J)
                  auto ptr_tgt_cell = pQueryCELL(I-1,J-1);
                  Int iDestJ = ptr_tgt_cell->owner;
                  //add the pair (col,row) to processor owning column row
                  send_map[iDestJ] += IdxToMin+PtrToMin+ (IdxToMin + TToMin);
                }
              }
            }
          }
        }
      }



      {
        //allocate one buffer for every remote processor
        //compute the send structures
        size_t total_send_size = 0;
        std::vector<size_t> stotcounts(this->all_np,0);
        for (auto it = send_map.begin(); it!=send_map.end();it++) {
          Int iCurDest = it->first;
          size_t & send_bytes = it->second;
          stotcounts[iCurDest] = send_bytes;
        }
        //compute send displacements
        std::vector<size_t> spositions(this->all_np+1,0);
        spositions[0] = 0;
        std::partial_sum(stotcounts.begin(),stotcounts.end(),&spositions[1]);
        total_send_size = spositions.back();

        std::vector<minType, Mallocator<minType> > sendBuffer(total_send_size);

        {
          scope_timer(a,symPACKMatrix2D::DistributeMatrix::Serializing);

          for (Int I=1;I<this->Xsuper_.size();I++) {
            Idx fc = this->Xsuper_[I-1];
            Idx lc = this->Xsuper_[I]-1;
            Int iWidth = lc-fc+1;

            //Serialize
            for (Idx col = fc;col<=lc;col++) {
              Idx orig_col = this->Order_.perm[col-1];

              if (orig_col>= FirstLocalCol && orig_col < LastLocalCol) {
                Ptr nrows = 0;
                Idx local_col = orig_col - FirstLocalCol+1;//1 based

                for (Ptr rowidx = pMat.Localg_.colptr[local_col-1]+(1-baseval); rowidx<pMat.Localg_.colptr[local_col]+(1-baseval); ++rowidx) {
                  Idx orig_row = pMat.Localg_.rowind[rowidx-1]+(1-baseval);
                  Idx row = this->Order_.invp[orig_row-1];

                  Int J = this->SupMembership_[row-1];

                  if (row>=col) {
                    //this has to go in cell(J,I)
                    auto ptr_tgt_cell = pQueryCELL(J-1,I-1);
                    bassert(ptr_tgt_cell->i==J && ptr_tgt_cell->j==I);
                    Int iDest = ptr_tgt_cell->owner;

                    T val = pMat.nzvalLocal[rowidx-1];

                    *((Idx*)&sendBuffer[spositions[iDest]]) = col;
                    spositions[iDest]+=IdxToMin;
                    *((Ptr*)&sendBuffer[spositions[iDest]]) = 1;
                    spositions[iDest]+=PtrToMin;
                    *((Idx*)&sendBuffer[spositions[iDest]]) = row;
                    spositions[iDest]+=IdxToMin;
                    *((T*)&sendBuffer[spositions[iDest]]) = val;
                    spositions[iDest]+=TToMin;
                  }
                  else {
                    //this has to go in cell(I,J)
                    auto ptr_tgt_cell = pQueryCELL(I-1,J-1);
                    bassert(ptr_tgt_cell->i==I && ptr_tgt_cell->j==J);
                    Int iDestJ = ptr_tgt_cell->owner;
                    //add the pair (col,row) to processor owning column row

                    T val = pMat.nzvalLocal[rowidx-1];

                    *((Idx*)&sendBuffer[spositions[iDestJ]]) = row;
                    spositions[iDestJ]+=IdxToMin;
                    *((Ptr*)&sendBuffer[spositions[iDestJ]]) = 1;
                    spositions[iDestJ]+=PtrToMin;
                    *((Idx*)&sendBuffer[spositions[iDestJ]]) = col;
                    spositions[iDestJ]+=IdxToMin;
                    *((T*)&sendBuffer[spositions[iDestJ]]) = val;
                    spositions[iDestJ]+=TToMin;
                  }
                }
              }
            }
          }
        }

        spositions[0] = 0;
        std::partial_sum(stotcounts.begin(),stotcounts.end(),&spositions[1]);


        size_t total_recv_size = 0;
        std::vector<minType, Mallocator<minType> > recvBuffer;
        std::function<void(std::vector<minType, Mallocator<minType> >&,size_t)> resize_lambda =
          [](std::vector<minType, Mallocator<minType> >& container, size_t sz) {
            container.resize(sz);
          };

        MPI_Datatype type;
        MPI_Type_contiguous( sizeof(minType), MPI_BYTE, &type );
        MPI_Type_commit(&type);


        mpi::Alltoallv(sendBuffer, &stotcounts[0], &spositions[0], type ,
            recvBuffer,this->fullcomm_, resize_lambda);

        total_recv_size = recvBuffer.size();

        MPI_Type_free(&type);
        //Need to parse the structure sent from the processor owning the first column of the supernode


        {
          scope_timer(a,symPACKMatrix2D::DistributeMatrix::Deserializing);
          size_t head = 0;

          while (head < total_recv_size)
          { 
            //Deserialize
            Idx col = *((Idx*)&recvBuffer[head]);
            head+=IdxToMin;
            //nrows of column col sent by processor p
            Ptr nrows = *((Ptr*)&recvBuffer[head]);
            head+=PtrToMin;
            Idx * rowind = (Idx*)(&recvBuffer[head]);
            head+=nrows*IdxToMin;
            T * nzvalA = (T*)(&recvBuffer[head]);
            head+=nrows*TToMin;

            Int I = this->SupMembership_[col-1];

            Int fc = this->Xsuper_[I-1];
            Int lc = this->Xsuper_[I]-1;
            Int iWidth = lc-fc+1;

            //Here, do a linear search instead for the blkidx
            Ptr colbeg = 1;
            Ptr colend = nrows;
            if (colbeg<=colend) {
              //sort rowind and nzvals
              std::vector<size_t> lperm = sort_permutation(&rowind[colbeg-1],&rowind[colend-1]+1,std::less<Idx>());
              apply_permutation(&rowind[colbeg-1],&rowind[colend-1]+1,lperm);
              apply_permutation(&nzvalA[colbeg-1],&nzvalA[colend-1]+1,lperm);
              Idx firstRow = rowind[colbeg-1];

              for (Ptr rowidx = colbeg; rowidx<=colend; ++rowidx) {
                Idx row = rowind[rowidx-1];

                Int J = this->SupMembership_[row-1];

                auto ptr_tgt_cell = pQueryCELL(J-1,I-1);
                assert(ptr_tgt_cell);

                auto & tgt_cell =CELL(J-1,I-1);
                bassert(this->iam == tgt_cell.owner);
                bassert(tgt_cell.i==J && tgt_cell.j==I);

                bool found = false;
                for (auto & block: tgt_cell.blocks()) {
                  //Match found if row between first and last row of block
                  if ( block.first_row <= row && row < block.first_row + tgt_cell.block_nrows(block) ) {
                    //Offset is block offset + row index in block & number of cols in tgt_cell
                    auto offset = block.offset + (row - block.first_row)*tgt_cell.width() + (col-fc);
                    tgt_cell._nzval[offset] = nzvalA[rowidx-1];
                    found = true;
                    break;

                  }  
                }

                bassert(found);

              }
              
            }
          }

        }
        

      }







    } 

  template <typename colptr_t, typename rowind_t, typename T, typename int_t>
    void symPACKMatrix2D<colptr_t,rowind_t,T,int_t>::Factorize( ) {
      using block_t = typename symPACK::symPACKMatrix2D<colptr_t, rowind_t, T, int_t>::snodeBlock_t::block_t;

#ifdef SP_THREADS
      this->scheduler.threadInitHandle_ = nullptr;
      this->scheduler.extraTaskHandle_  = nullptr;

      this->scheduler.threadInitHandle_ = [&,this]() {
        std::thread::id tid = std::this_thread::get_id();
        std::lock_guard<std::recursive_mutex> lock( this->scheduler.scheduler_mutex_);
        auto & tmpBuf = tmpBufs_th[tid];
      };

      if (Multithreading::NumThread>2) {
        this->scheduler.extraTaskHandle_ = [&,this](SparseTask2D * pTask)->bool {
          auto type = std::get<2>(pTask->_meta);

          if (type == Factorization::op_type::FACTOR) {
            return false;
          }
          else if (type == Factorization::op_type::DLT2D_COMP) {
            return false;
          }
          else {
            auto src = std::get<0>(pTask->_meta);
            auto J = std::get<1>(pTask->_meta);
            auto facing_row = std::get<4>(pTask->_meta);
            auto K = this->SupMembership_[facing_row-1];

            auto ptr_tgtcell = pQueryCELL2(K-1,J-1);
            bool exp = false;
            pTask->_lock_ptr = &ptr_tgtcell->in_use;
            if (std::atomic_compare_exchange_weak( &ptr_tgtcell->in_use, &exp, true )) {
              return false;
            }
            else {
              return true;
            }
          }
        };
      }
#endif

#ifdef _TIMING_
      std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
#endif
      this->scheduler.execute(this->task_graph,this->mem_budget);
#ifdef _TIMING_
      std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
      double execute_graph_ticks = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();
#endif

#ifdef _TIMING_
      std::stringstream sstr;
      sstr<<upcxx::rank_me()<<" "<<(double)CELL_ticks*1.0e-9<<" "<<(double)rpc_fact_ticks*1.0e-9<<" "<<(double)rpc_trsm_ticks*1.0e-9<<" "<<(double)rpc_upd_ticks*1.0e-9<<std::endl;
      sstr<<upcxx::rank_me()<<" "<<(double)deps_fact_ticks*1.0e-9<<" "<<(double)deps_trsm_ticks*1.0e-9<<" "<<(double)deps_upd_ticks*1.0e-9<<std::endl;
      sstr<<upcxx::rank_me()<<" "<<(double)comp_fact_ticks*1.0e-6<<" "<<(double)comp_trsm_ticks*1.0e-6<<" "<<(double)comp_upd_ticks*1.0e-6<<std::endl;
      sstr<<upcxx::rank_me()<<" "<<(double)execute_graph_ticks*1.0e-6<<std::endl;
      logfileptr->OFS()<<sstr.str();
#endif
    } 


  template <typename colptr_t, typename rowind_t, typename T, typename int_t>
    void symPACKMatrix2D<colptr_t,rowind_t,T,int_t>::Solve( T * RHS, int nrhs, int rhs_size,  T * Xptr ) {
      //set solve_data
      //TODO RHS should be distributed according to Xsuper
      std::vector<T> distRHS;
      this->solve_data.rhs = RHS;
      this->solve_data.nrhs = nrhs;
      this->solve_data.contribs.clear();
      this->solve_data.contribs.resize(nsuper+1);
      delete [] this->solve_data.contribs_lock;
      this->solve_data.contribs_lock = new std::atomic<bool>[this->nsuper+1];
      for ( int i = 0; i < this->nsuper+1; i++ ) { this->solve_data.contribs_lock[i] = false; } 

      this->solve_data.update_right_cnt.assign(this->nsuper+1,0);
      this->solve_data.update_up_cnt.assign(this->nsuper+1,0);


#ifdef SP_THREADS
      this->scheduler.threadInitHandle_ = nullptr;
      this->scheduler.extraTaskHandle_  = nullptr;

      this->scheduler.quiesceHandle_ = [this](){
        while ( (*this->solve_data.remoteDeallocCounter) > 0 ) { upcxx::progress(); }
      };

      this->scheduler.threadInitHandle_ = [&,this]() {
        std::thread::id tid = std::this_thread::get_id();
        std::lock_guard<std::recursive_mutex> lock( this->scheduler.scheduler_mutex_);
        auto & tmpBuf = tmpBufs_th[tid];
      };


      if (Multithreading::NumThread>2) {
        this->scheduler.extraTaskHandle_ = [&,this](SparseTask2D * pTask)->bool {
          auto type = std::get<2>(pTask->_meta);

          if (type == Factorization::op_type::BUC) {
            auto J = std::get<0>(pTask->_meta);
            auto I = std::get<1>(pTask->_meta);
            bool exp = false;
            auto & lock = this->solve_data.contribs_lock[I];
            pTask->_lock_ptr = &lock;
            if (std::atomic_compare_exchange_weak( &lock, &exp, true )) {
              return false;
            }
            else {
              return true;
            }
          }
          else {
            auto J = std::get<0>(pTask->_meta);
            auto I = std::get<1>(pTask->_meta);
            bool exp = false;
            auto & lock = this->solve_data.contribs_lock[J];
            pTask->_lock_ptr = &lock;
            if (std::atomic_compare_exchange_weak( &lock, &exp, true )) {
              return false;
            }
            else {
              return true;
            }
          }
        };
      }
#endif

      upcxx::barrier();
      this->scheduler.execute(this->task_graph_solve,this->mem_budget);
  } 

  template <typename colptr_t, typename rowind_t, typename T, typename int_t>
    void symPACKMatrix2D<colptr_t,rowind_t,T,int_t>::GetSolution(T * B, int nrhs) {
      Int n = this->iSize_;
      {
        std::fill(B,B+n*nrhs,T(0.0));

        //Gather B from everybody and put it in the original matrix order
        for ( int I = 1; I <= this->nsuper; I++ ) {
          auto ptr_tgt_cell = pQueryCELL(I-1,I-1);
          Int iOwner = ptr_tgt_cell->owner;
          T * data;
          Int snode_size = this->Xsuper_[I] - this->Xsuper_[I-1];
          Int nzcnt = snode_size * nrhs;

          if ( iOwner == this->iam ) {
            auto & contrib_slot = this->solve_data.contribs[I];
            auto & ptr_tgt_cell = std::get<1>(contrib_slot);
            bassert(ptr_tgt_cell!=nullptr);
            auto & tgt_cell = *std::dynamic_pointer_cast<snodeBlock_t>(ptr_tgt_cell);
            if (this->options_.decomposition == DecompositionType::LDL) {
              auto & tgt_ldlcell = *std::dynamic_pointer_cast<snodeBlockLDL_t>(ptr_tgt_cell);
              bassert(tgt_ldlcell.scaled);
            }


            for (auto & block: tgt_cell.blocks()) {
              T * val = &tgt_cell._nzval[block.offset];
              auto nRows = tgt_cell.block_nrows(block);
              auto row = block.first_row;
              for (auto i = 0; i< nRows; ++i) {
                for (auto j = 0; j< tgt_cell.width(); ++j) {

                  Int destRow = tgt_cell.first_col + i;
                  destRow = this->Order_.perm[destRow - 1];

                  B[ (destRow -1) +j*this->iSize_ ] = val[i*tgt_cell.width()+j];
                }
              }
            }            
          }
        }

        mpi::Allreduce((T*)MPI_IN_PLACE,&B[0],n*nrhs,MPI_SUM,this->fullcomm_);
        MPI_Barrier(this->fullcomm_);
      }
    }




  //TODO redo this
  template <typename colptr_t, typename rowind_t, typename T, typename int_t>
    inline void symPACKMatrix2D<colptr_t,rowind_t,T,int_t>::DumpMatlab() {
      logfileptr->OFS()<<"+sparse([";
      for (Int I=1;I<this->Xsuper_.size();I++) {
        Int src_first_col = this->Xsuper_[I-1];
        Int src_last_col = this->Xsuper_[I]-1;

        for (Int J=I;J<this->Xsuper_.size();J++) {
          auto idx = coord2supidx(J-1,I-1);
          if (this->cells_.find(idx)!= this->cells_.end()) {


            auto ptr_tgt_cell = pQueryCELL(J-1,I-1);
            assert(ptr_tgt_cell!=nullptr);


            Int iOwner = ptr_tgt_cell->owner;
            if ( iOwner == this->iam ) {
              auto & tgt_cell = *std::dynamic_pointer_cast<snodeBlock_t>(ptr_tgt_cell);
              for (auto & block: tgt_cell.blocks()) {
                T * val = &tgt_cell._nzval[block.offset];
                Int nRows = tgt_cell.block_nrows(block);

                Int row = block.first_row;
                for (Int i = 0; i< nRows; ++i) {
                  for (Int j = 0; j< tgt_cell.width(); ++j) {
                    logfileptr->OFS()<<row+i<<" ";
                  }
                }
              }
            }
          }
        }
      }
      logfileptr->OFS()<<"],[";
      for (Int I=1;I<this->Xsuper_.size();I++) {
        Int src_first_col = this->Xsuper_[I-1];
        Int src_last_col = this->Xsuper_[I]-1;

        for (Int J=I;J<this->Xsuper_.size();J++) {
          auto idx = coord2supidx(J-1,I-1);
          if (this->cells_.find(idx)!= this->cells_.end()) {
            auto ptr_tgt_cell = pQueryCELL(J-1,I-1);
            Int iOwner = ptr_tgt_cell->owner;
            if ( iOwner == this->iam ) {
              auto & tgt_cell = *std::dynamic_pointer_cast<snodeBlock_t>(ptr_tgt_cell);
              for (auto & block: tgt_cell.blocks()) {
                T * val = &tgt_cell._nzval[block.offset];
                Int nRows = tgt_cell.block_nrows(block);

                Int row = block.first_row;
                for (Int i = 0; i< nRows; ++i) {
                  for (Int j = 0; j< tgt_cell.width(); ++j) {
                    logfileptr->OFS()<<src_first_col+j<<" ";
                  }
                }
              }
            }
          }
        }
      }
      logfileptr->OFS()<<"],[";
      logfileptr->OFS().precision(std::numeric_limits< T >::max_digits10);
      for (Int I=1;I<this->Xsuper_.size();I++) {
        Int src_first_col = this->Xsuper_[I-1];
        Int src_last_col = this->Xsuper_[I]-1;

        for (Int J=I;J<this->Xsuper_.size();J++) {
          auto idx = coord2supidx(J-1,I-1);
          auto it = this->cells_.end();
          if ( (it = this->cells_.find(idx)) != this->cells_.end()) {
            auto ptr_tgt_cell = pQueryCELL(J-1,I-1);
            Int iOwner = ptr_tgt_cell->owner;
            if ( iOwner == this->iam ) {
              auto & tgt_cell = *std::dynamic_pointer_cast<snodeBlock_t>(ptr_tgt_cell);
              for (auto & block: tgt_cell.blocks()) {
                T * val = &tgt_cell._nzval[block.offset];
                auto nRows = tgt_cell.block_nrows(block);

                auto row = block.first_row;
                for (auto i = 0; i< nRows; ++i) {
                  for (auto j = 0; j< tgt_cell.width(); ++j) {
                    logfileptr->OFS()<<std::scientific<<"("<<row+i<<","<<src_first_col+j<<") "<<ToMatlabScalar(val[i*tgt_cell.width()+j])<<" "<<std::endl;
                  }
                }
              }
            }
          }
        }
      }
      logfileptr->OFS()<<"],"<<this->iSize_<<","<<this->iSize_<<")"<<std::endl;
    }













  namespace scheduling {
#ifdef SP_THREADS
    template <typename ttask_t , typename ttaskgraph_t >
      inline bool Scheduler2D<ttask_t,ttaskgraph_t>::assignWork(ttask_t * ptask) 
      {
        if ( ! free_workers_.empty() ) {
          auto id = free_workers_.front();
          auto & w_p = work_pointers_[id];
          bool success = w_p.set(ptask);
          if ( success ) {
            bassert(worker_loads_[id]<2);
            ++worker_loads_[id]; 
            if ( worker_loads_[id] == 2 ) {
              free_workers_.pop_front();
              bassert(std::find(this->free_workers_.begin(),this->free_workers_.end(),id)==this->free_workers_.end());
            }
          }
          return success;
        }
        else {
          return false;
        }
      }
#endif

    template <typename ttask_t , typename ttaskgraph_t >
      inline void Scheduler2D<ttask_t,ttaskgraph_t>::execute(ttaskgraph_t & task_graph, double & mem_budget )
      {
        double progress_ticks = 0;
        double progress_ticks2 = 0;
        double message_ticks = 0;
        double execute_ticks = 0;
        double while_ticks = 0;
        double ifempty_ticks = 0;
        try{
          int64_t local_task_cnt;
          {
            local_task_cnt = task_graph.size();
            for (auto it = task_graph.begin(); it != task_graph.end(); it++) {
              auto & ptask = *it;
              auto remote_deps = ptask->in_remote_dependencies_cnt;
              auto ptr = ptask.get();

#ifdef _USE_PROM_AVAIL_
              auto fut_comm = ptask->in_avail_prom.finalize();
              if (remote_deps >0 ) {
                fut_comm.then([this,ptr]() {
                    push_avail(this,ptr);
                    });
              }
#else
              if (remote_deps >0 ) {
                if (ptr->in_avail_counter == 0 ) {
                  push_avail(this,ptr);
                }
              }
#endif

#ifdef _USE_PROM_RDY_
              auto fut = ptask->in_prom.finalize();
              fut.then([this, ptr]() {
#ifdef SP_THREADS
                  if (this->extraTaskHandle_!=nullptr) {
                  bool delay = this->extraTaskHandle_(ptr);
                  if (delay) {
                  this->delayedTasks_[ptr->_lock_ptr].push_back(ptr);
                  return;
                  }
                  }
#endif
                  push_ready(this,ptr);
                  });
#else
              if (ptr->in_counter == 0 ) {
#ifdef SP_THREADS
                if (this->extraTaskHandle_!=nullptr) {
                  bool delay = this->extraTaskHandle_(ptr);
                  if (delay) {
                    this->delayedTasks_[ptr->_lock_ptr].push_back(ptr);
                  }
                  else 
                  {
                    push_ready(this,ptr);
                  }
                } else 
#endif
                {
                  push_ready(this,ptr);
                }
              }
#endif
            }
          }

#ifdef SP_THREADS
          if (Multithreading::NumThread>1) {
            std::atomic<bool> done(false);
            int nthreads = std::max(1,Multithreading::NumThread-1);
            this->work_pointers_ = std::vector<worker_ptrs<ttask_t>>(nthreads);

            std::vector<std::thread> workers;
            workers.reserve(nthreads);

            this->free_workers_.clear();
            this->worker_loads_.clear();
            this->worker_loads_.resize(nthreads,0);
            for (int count {0}; count < nthreads; count += 1) {
              bassert(std::find(this->free_workers_.begin(),this->free_workers_.end(),count)==this->free_workers_.end());
              this->free_workers_.push_back(count);
            }

            for (int count {0}; count < nthreads; count += 1) {
              workers.emplace_back(
                  [this,&done](int tid) {
                  if (this->threadInitHandle_!=nullptr) {
                  this->threadInitHandle_();
                  }
                  auto & w_p = this->work_pointers_[tid];

                  bool local_done = false;
                  int drain = 0;
                  while (true) {
                  ttask_t * t = w_p.get();
                  while (t!=nullptr) { 
                  t->execute();
                  auto lock_ptr = t->_lock_ptr;
                  t->reset();

                  int J=std::get<0>(t->_meta);
                  int I=std::get<1>(t->_meta);
                  upcxx::master_persona().lpc_ff( 
                      [this,tid,I,J,lock_ptr] () {
                      // if _lock_ptr is not null, it means that t had to acquire a lock,
                      // which means it may have to unlock other tasks
                      if ( lock_ptr ) {
                      bassert(this->extraTaskHandle_!=nullptr);
                      auto & task_slot = this->delayedTasks_[lock_ptr];
                      if ( ! task_slot.empty() ) {
                      //get the oldest delayed task waiting on that lock
                      auto nt = task_slot.front();
                      bassert(nt);
                      bassert(nt->_lock_ptr== lock_ptr);
                      //check if task needs to be delayed: this will also take the lock if no delay is required
                      bool delay = extraTaskHandle_(nt);
                      if (!delay) {
                      //remove it from the list of delayed tasks
                      task_slot.pop_front();
                      //push it in the list of ready tasks
                      push_ready(this, nt );
                      }
                      }
                      //if the list is empty, then there is no point in having it in delayed tasks anymore
                      if ( task_slot.empty() ) {
                        this->delayedTasks_.erase(lock_ptr);
                      }
                      }

                      //decrement the load of that worker
                      bassert(this->worker_loads_[tid]>=1);
                      if (--this->worker_loads_[tid]==1) {
                        bassert(std::find(this->free_workers_.begin(),this->free_workers_.end(),tid)==this->free_workers_.end());
                        this->free_workers_.push_back(tid);
                      }
                      else {
                        bassert(std::find(this->free_workers_.begin(),this->free_workers_.end(),tid)!=this->free_workers_.end());
                      }
                      });
                  t = w_p.get();
                  }

                  if (done.load(std::memory_order_acquire) ) {
                    break;
                  }
                  else {
                    sched_yield();
                  }
                  }
                  }
              , count);
            }

            while (local_task_cnt>0) {
              bool empty = false;
              while (!empty && !this->free_workers_.empty()) {
                ttask_t * ptask = nullptr;
                {
                  empty = ready_tasks.empty();
                  if (!empty) {
                    ptask = top_ready();
                  }
                }

                if (empty)
                  break;

                if (!empty) {
                  bool delay = false;
                  if (!delay) {
                    if (assignWork(ptask)) {
                      local_task_cnt--;
                      pop_ready();
                    }
                    else {
                      break;
                    }
                  }
                }
              }

              upcxx::progress();
              sched_yield();

              //handle communications
              empty = false;
              while (!empty) {
                ttask_t * ptask = nullptr;
                {
                  empty = avail_tasks.empty();

                  if (!empty) {
                    ptask = top_avail();
                    pop_avail();
                  }
                }
                if (!empty) {
                  for (auto & msg : ptask->input_msg) {
                    msg->allocate();
                    msg->fetch().then([this,ptask](incoming_data_t<ttask_t,meta_t> * pmsg) {
                        //fulfill promise by one, when this reaches 0, ptask is moved to scheduler.ready_tasks
                        ptask->satisfy_dep(1,*this);
                        });
                  }
                  upcxx::progress(upcxx::progress_level::internal);
                }
              }
            }
            upcxx::progress();
            upcxx::discharge();
            done = true;
            for (auto && thread: workers) thread.join();
            //ensure the last lpc_ff are executed;
            if ( this->quiesceHandle_ != nullptr ) this->quiesceHandle_();
            upcxx::progress();
            upcxx::barrier();
          }
          else
#endif
          {
            while (local_task_cnt>0) {
              if (!ready_tasks.empty()) {
                upcxx::progress(upcxx::progress_level::internal);
                auto ptask = top_ready();
                pop_ready();
                ptask->execute(); 
                local_task_cnt--;
              }
              upcxx::progress();
              //handle communications
              if (!avail_tasks.empty()) {
                //get all comms from a task
                auto ptask = top_avail();
                pop_avail();
                for (auto & msg : ptask->input_msg) {
                  msg->allocate();
                  msg->fetch().then([this,ptask](incoming_data_t<ttask_t,meta_t> * pmsg) {
                      //fulfill promise by one, when this reaches 0, ptask is moved to scheduler.ready_tasks
                      ptask->satisfy_dep(1,*this);
                      });
                } 
	      	upcxx::progress(upcxx::progress_level::internal);
              }
            }
	    upcxx::progress();
            upcxx::discharge();
#ifdef SP_THREADS
            if ( this->quiesceHandle_ != nullptr ) this->quiesceHandle_();
#endif
            upcxx::barrier();
          }
        }
        catch(const std::runtime_error& e) {
          std::cerr << "Runtime error: " << e.what() << '\n';
        }
      }
  }

}

#include <sympack/impl/symPACKMatrix2D_impl.hpp>

#endif //_SYMPACK_MATRIX2D_DECL_HPP_

