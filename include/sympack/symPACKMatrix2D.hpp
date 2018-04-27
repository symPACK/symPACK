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
#ifndef _SYMPACK_MATRIX2D_DECL_HPP_
#define _SYMPACK_MATRIX2D_DECL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/symPACKMatrix.hpp"




namespace symPACK{


#ifdef NEW_UPCXX

#define _DEBUG_DEPENDENCIES_

#define pCELL(a,b) (std::static_pointer_cast<snodeBlock_t>(this->cells_[coord2supidx((a),(b))]))
#define CELL(a,b) (*pCELL((a),(b)))

  class blockCellBase_t {
    protected:
      blockCellBase_t(){}
  };

namespace scheduling {

  template <typename ttask_t, typename tmeta_t>
  class incoming_data_t;
  
  using meta_t = std::tuple<Idx,Idx,Factorization::op_type,Idx,Idx>;
#ifdef _DEBUG_DEPENDENCIES_
  using depend_t = std::tuple<Int,Int,bool>; //(I,J,satisfied)
#else
  using depend_t = std::tuple<Int,Int>; //(I,J)
#endif

  template
    < typename tmeta_t = meta_t,
      typename tdepend_t = depend_t >
  class task_t{
    public:
      //cell description
      using depend_task_t = tdepend_t;
      using meta_t = tmeta_t;
      using data_t = incoming_data_t<task_t,tmeta_t>;

      meta_t _meta;

      std::function< void() > execute;

      //byte storage for tasks 
      std::deque< depend_task_t > out_dependencies;
      std::deque< depend_task_t > in_remote_dependencies;
      std::deque< depend_task_t > in_local_dependencies;

      std::deque< upcxx::global_ptr<char> > in_ptr;

      //promise to sync all outgoing RPCs
      upcxx::promise<> * out_prom;
      //promise to wait for all incoming_data_t to be created
      upcxx::promise<> * in_avail_prom;
      //promise to sync on all incoming_data_t fetch
      upcxx::promise<> * in_prom;

      int dep_count;
      bool executed;

      task_t( ):dep_count(0),executed(false){
        out_prom = new upcxx::promise<>();
        in_prom = new upcxx::promise<>();
        in_avail_prom = new upcxx::promise<>();
      }

      ~task_t(){
        bassert(out_prom->get_future().ready());
        bassert(in_prom->get_future().ready());
        bassert(in_avail_prom->get_future().ready());
        delete out_prom;
        delete in_prom;
        delete in_avail_prom;
      }

      void reset(){
        bassert(out_prom->get_future().ready());
        bassert(in_prom->get_future().ready());
        bassert(in_avail_prom->get_future().ready());
        delete out_prom;
        delete in_prom;
        delete in_avail_prom;
        out_prom = new upcxx::promise<>();
        in_prom = new upcxx::promise<>();
        in_avail_prom = new upcxx::promise<>();
      }

      //need counter here as same REMOTE cell can be input to many tasks.
      //std::deque< std::shared_ptr< blockCellBase_t > > input_data;
      std::deque< std::shared_ptr< incoming_data_t<task_t,meta_t> > > input_msg;
  };


  template <typename ttask_t = task_t<meta_t, depend_t> , typename tmeta_t = task_t<meta_t, depend_t>::meta_t>
  class incoming_data_t {
    public:
      upcxx::promise<> on_delete;
      upcxx::promise<incoming_data_t *> on_fetch;

      using task_t = ttask_t;
      using meta_t = tmeta_t;
      //this should be made generic, not specific to sparse matrices
      std::deque< task_t *> target_tasks;
      meta_t in_meta;
      //a pointer to be used by the user if he wants to attach some data
      void * extra_data;

      incoming_data_t():transfered(false),size(0),extra_data(nullptr)/*,cell(nullptr)*/,landing_zone(nullptr){};
      bool transfered;
      upcxx::global_ptr<char> remote_gptr;
      size_t size;
      char * landing_zone;
      //this has to be made generic...
      //blockCellBase_t * cell;

      void allocate(){
          if(landing_zone == nullptr)
            landing_zone = new char[size];
      }

      upcxx::future<incoming_data_t *> fetch(){
        if(!transfered){
          transfered=true;
          upcxx::rget(remote_gptr,landing_zone,size).then([this](){
              //gdb_lock();
              on_fetch.fulfill_result(this);
              return;});
        }
        return on_fetch.get_future();
      }

      ~incoming_data_t(){
        //wait on callbacks
        on_delete.finalize().wait();

        //this has to be made generic
        //if (cell) {
        //  delete cell;
        //}
        if (landing_zone) {
          delete [] landing_zone;
        }
      }
    
  };

  //template <typename ttask_t> class Scheduler2D;
  //using taskGraph2D_t = std::unordered_map< Int,  std::deque< SparseTask2D * > > ;
  template <typename bin_label_t = Int, typename ttask_t = task_t<meta_t, depend_t> >
  class task_graph_t: public std::unordered_map< bin_label_t,  std::deque< ttask_t * > > {
  };

  template <typename ttask_t = task_t<meta_t, depend_t> , typename ttaskgraph_t = task_graph_t<Int,task_t<meta_t, depend_t> > >
  class Scheduler2D{
    public:
    int sp_handle;
    std::deque<ttask_t*> ready_tasks;
    std::deque<ttask_t*> avail_tasks;
#ifdef _DEBUG_DEPENDENCIES_
    void execute(ttaskgraph_t & graph, Int * SupMembership );
#else
    void execute(ttaskgraph_t & graph );
#endif
  };

}


  //extern std::map<int, std::deque< scheduling::incoming_data_t >  > g_sp_handle_incoming;
  //extern std::map<int, taskGraph2D_t *  > g_sp_handle_to_graph;
  extern std::map<int, symPACKMatrixBase *  > g_sp_handle_to_matrix;


  template <typename colptr_t, typename rowind_t, typename T, typename int_t = int> 
    class blockCell_t: public blockCellBase_t {
      private:
        using intrank_t = upcxx::intrank_t;

      public:
        class block_t {
          public:
            //this is just an offset from the begining of the cell;
            size_t _offset;
            rowind_t first_row;
            rowind_t nrows;
        };

        int_t i;
        int_t j;
        intrank_t owner;
        rowind_t first_col;
        std::tuple<rowind_t> _dims;
        upcxx::global_ptr< char > _gstorage;
        char * _storage;
        block_t * _blocks;
        size_t _cblocks; //capacity
        size_t _nblocks;
        T* _nzval;
        size_t _cnz; //capacity
        size_t _nnz;
        size_t _storage_size;
        rowind_t _total_rows;


        blockCell_t():i(0),j(0),owner(0),
        first_col(0),_dims(std::make_tuple(0)),_total_rows(0),
        _storage(nullptr),_nzval(nullptr),_blocks(nullptr),
        _gstorage(nullptr),_nnz(0),_storage_size(0),_nblocks(0){}

        ~blockCell_t() {
          if ( _gstorage.is_null() ) {
            bassert(_storage == _gstorage.local());
            upcxx::deallocate( _gstorage );
          }
          else {
            delete [] _storage;
          }
        }

        rowind_t total_rows(){
          if(this->_nblocks>0 && _total_rows==0){
            for( int i = 0; i< _nblocks; i++){ _total_rows += _blocks[i].nrows;}
          }
          return _total_rows;
        }

        blockCell_t ( rowind_t width, size_t nzval_cnt, size_t block_cnt, bool shared_segment = true  ): blockCell_t() {
          _dims = std::make_tuple(width);
          allocate(nzval_cnt,block_cnt, shared_segment);
          initialize(nzval_cnt,block_cnt);
        }

        blockCell_t ( char * ext_storage, rowind_t width, size_t nzval_cnt, size_t block_cnt ): blockCell_t() {
          _gstorage = nullptr;
          _storage = ext_storage;
          _dims = std::make_tuple(width);
          initialize(nzval_cnt,block_cnt);
          _nnz = nzval_cnt;
          _nblocks = block_cnt;
        }

        blockCell_t ( upcxx::global_ptr<char> ext_gstorage, rowind_t width, size_t nzval_cnt, size_t block_cnt ): blockCell_t() {
          _gstorage = ext_gstorage;
          _storage = _gstorage.local();
          _dims = std::make_tuple(width);
          initialize(nzval_cnt,block_cnt);
          _nnz = nzval_cnt;
          _nblocks = block_cnt;
        }

        // Copy constructor.  
        blockCell_t ( const blockCell_t & other ): blockCell_t() {
          i           = other.i;
          j           = other.j;
          owner       = other.owner;
          first_col   = other.first_col;
          _dims       = other._dims;
          _total_rows = other._total_rows;

          allocate( other.nz_capacity(), other.block_capacity , !other._gstorage.is_null() );
          //now copy the data
          std::copy( other._storage, other._storage + other._storage_size, _storage );
          _nblocks = other._nblocks;
          _nnz = other._nnz;
        }

        // Move constructor.  
        blockCell_t ( const blockCell_t && other ): blockCell_t() {
          i           = other.i;
          j           = other.j;
          owner       = other.owner;
          first_col   = other.first_col;
          _dims       = other._dims;
          _total_rows = other._total_rows;

          _gstorage = other._gstorage;
          _storage = other._storage;
          initialize( other._cnz , other._cblocks );

          //invalidate other
          other._gstorage = nullptr;
          other._storage = nullptr;
          other._storage_size = 0;
          other._cblocks = 0;
          other._cnz = 0;
          other._nblocks = 0;
          other._nnz = 0;
          other._blocks = nullptr;
          other._nzval = nullptr;
        }




        // Copy assignment operator.  
        blockCell_t& operator=(const blockCell_t& other)  {  
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

            allocate( other._cnz, other._cblocks , !other._gstorage.is_null() );
            //now copy the data
            std::copy( other._storage, other._storage + other._storage_size, _storage );
            _nblocks = other._nblocks;
            _nnz = other._nnz;

          } 
          return *this;  
        }  

        // Move assignment operator.  
        blockCell_t& operator=(const blockCell_t&& other)  {  
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

            _gstorage = other._gstorage;
            _storage = other._storage;
            initialize( other._cnz , other._cblocks );

            //invalidate other
            other._gstorage = nullptr;
            other._storage = nullptr;
            other._storage_size = 0;
            other._cblocks = 0;
            other._cnz = 0;
            other._nblocks = 0;
            other._nnz = 0;
            other._blocks = nullptr;
            other._nzval = nullptr;        
          } 
          return *this;  
        }  




        inline int_t block_capacity() const { return _cblocks; }
        inline int_t nblocks() const { return _nblocks; }
        inline int_t nz_capacity() const { return _cnz; }
        inline int_t nnz() const { return _nnz; }

        void add_block( rowind_t first_row, rowind_t nrows ){
          bassert( nblocks() + 1 <= block_capacity() );
          bassert( nnz() + nrows*std::get<0>(_dims) <= nz_capacity() );
          auto & new_block = _blocks[_nblocks++];
          new_block.first_row = first_row;
          new_block.nrows = nrows;
          new_block.offset = _nnz;
          _nnz += nrows*std::get<0>(_dims);
        }

        void initialize ( size_t nzval_cnt, size_t block_cnt ) {
          _storage_size = nzval_cnt*sizeof(T) + block_cnt*sizeof(block_t);
          _cnz = nzval_cnt;
          _cblocks = block_cnt;
          _nnz = 0;
          _nblocks = 0;
          _blocks = reinterpret_cast<block_t*>( _storage );
          _nzval = reinterpret_cast<T*>( _blocks + _cblocks );
        }

        void allocate ( size_t nzval_cnt, size_t block_cnt, bool shared_segment ) {
          if ( shared_segment ) {
            _gstorage = upcxx::allocate<char>( nzval_cnt*sizeof(T) + block_cnt*sizeof(block_t) );
            _storage = _gstorage.local();
          }
          else {
            _gstorage = nullptr;
            _storage = new char[nzval_cnt*sizeof(T) + block_cnt*sizeof(block_t)];
          }
          initialize( nzval_cnt, block_cnt );
        }

        int factorize( TempUpdateBuffers<T> & tmpBuffers){
          scope_timer(a,blockCell_t::factorize);
#if defined(_NO_COMPUTATION_)
          return 0;
#endif
          auto snode_size = std::get<0>(_dims);
          auto diag_nzval = _nzval;
          lapack::Potrf( 'U', snode_size, diag_nzval, snode_size);
          return 0;
        }

        int trsm( const blockCell_t & diag, TempUpdateBuffers<T> & tmpBuffers){
          scope_timer(a,blockCell_t::trsm);
#if defined(_NO_COMPUTATION_)
          return 0;
#endif

          bassert(diag.nblocks()>0);
          auto diag_nzval = diag._nzval;

          auto snode_size = std::get<0>(_dims);
          auto nzblk_nzval = _nzval;
          blas::Trsm('L','U','T','N',snode_size, total_rows(), T(1.0),  diag_nzval, snode_size, nzblk_nzval, snode_size);
          return 0;
        }

        int update( /*const*/ blockCell_t & pivot, /*const*/ blockCell_t & facing, TempUpdateBuffers<T> & tmpBuffers){
          scope_timer(a,blockCell_t::update);
#if defined(_NO_COMPUTATION_)
          return 0;
#endif
          //do the owner compute update first



          {
            bassert(nblocks()>0);
            auto & first_block = _blocks[0];

            bassert(pivot.nblocks()>0);
            auto pivot_fr = pivot._blocks[0].first_row;

            bassert(facing.nblocks()>0);
            auto facing_fr = facing._blocks[0].first_row;

            auto src_snode_size = std::get<0>(pivot._dims);
            auto tgt_snode_size = std::get<0>(_dims);

            //find the first row updated by src_snode
            auto tgt_fc = pivot_fr;
            auto tgt_lc = pivot_fr+pivot._blocks[0].nrows-1;;

            int_t first_pivot_idx = 0;
            int_t last_pivot_idx = pivot.nblocks()-1;

            //determine the first column that will be updated in the target supernode
            rowind_t tgt_local_fc =  tgt_fc - this->first_col;
            rowind_t tgt_local_lc =  tgt_lc - this->first_col;

            rowind_t pivot_nrows = pivot.total_rows();
            rowind_t tgt_nrows = this->total_rows();
            rowind_t src_nrows = facing.total_rows();
            rowind_t src_lr = facing_fr + src_nrows-1;

            //condensed update width is the number of rows in the pivot block
            int_t tgt_width = pivot_nrows;

            T * pivot_nzval = pivot._nzval;
            T * facing_nzval = facing._nzval;
            T * tgt = this->_nzval;

            //Pointer to the output buffer of the GEMM
            T * buf = NULL;
            T beta = T(0);
#ifdef SP_THREADS
            tmpBuffers.tmpBuf.resize(tgt_width*src_nrows + src_snode_size*tgt_width);
#endif

            //If the target supernode has the same structure,
            //The GEMM is directly done in place
            if(src_nrows == tgt_nrows){
              buf = &tgt[tgt_local_fc];
              beta = T(1);
            }
            else{
              //Compute the update in a temporary buffer
              buf = &tmpBuffers.tmpBuf[0];
            }

            //everything is in row-major
            SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_GEMM);
            blas::Gemm('T','N',tgt_width, src_nrows,src_snode_size,
                T(-1.0),pivot_nzval,src_snode_size,
                facing_nzval,src_snode_size,beta,buf,tgt_width);
            SYMPACK_TIMER_SPECIAL_STOP(UPDATE_SNODE_GEMM);

            //If the GEMM wasn't done in place we need to aggregate the update
            //This is the assembly phase
            if(src_nrows != tgt_nrows){
              //now add the update to the target supernode
              //              if(tgt_snode_size==1){
              //                int_t rowidx = 0;
              //                int_t src_blkcnt = src_snode->NZBlockCnt();
              //                for(int_t blkidx = first_pivot_idx; blkidx < src_blkcnt; ++blkidx){
              //                  NZBlockDesc & cur_block_desc = src_snode->GetNZBlockDesc(blkidx);
              //                  int_t cur_src_nrows = src_snode->NRows(blkidx);
              //                  int_t cur_src_lr = cur_block_desc.GIndex + cur_src_nrows -1;
              //                  int_t cur_src_fr = std::max(tgt_fc, cur_block_desc.GIndex);
              //
              //                  int_t row = cur_src_fr;
              //                  while(row<=cur_src_lr){
              //                    int_t tgt_blk_idx = FindBlockrowind_t(row);
              //                    assert(tgt_blk_idx>=0);
              //                    NZBlockDesc & cur_tgt_desc = GetNZBlockDesc(tgt_blk_idx);
              //                    int_t lr = std::min(cur_src_lr,cur_tgt_desc.GIndex + NRows(tgt_blk_idx)-1);
              //                    int_t tgtOffset = cur_tgt_desc.Offset 
              //                      + (row - cur_tgt_desc.GIndex)*tgt_snode_size;
              //                    SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_ADD_SINGLE);
              //                    for(int_t cr = row ;cr<=lr;++cr){
              //                      tgt[tgtOffset + (cr - row)*tgt_snode_size] += buf[rowidx]; 
              //                      rowidx++;
              //                    }
              //                    SYMPACK_TIMER_SPECIAL_STOP(UPDATE_SNODE_ADD_SINGLE);
              //                    row += (lr-row+1);
              //                  }
              //                }
              //              }
              //              else
              {
                SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_INDEX_MAP);
                tmpBuffers.src_colindx.resize(tgt_width);
                tmpBuffers.src_to_tgt_offset.resize(src_nrows);
                colptr_t colidx = 0;
                colptr_t rowidx = 0;
                size_t offset = 0;

                block_t * tgt_ptr = _blocks;

                for ( int block_idx = 0; block_idx < facing.nblocks() ; block_idx++ ) {
                  auto & cur_block = facing._blocks[block_idx];
                  rowind_t cur_src_nrows = cur_block.nrows;
                  rowind_t cur_src_lr = cur_block.first_row + cur_src_nrows -1;
                  rowind_t cur_src_fr = cur_block.first_row;

                  //The other one MUST reside into a single block in the target
                  rowind_t row = cur_src_fr;
                  while(row<=cur_src_lr){
                    do {
                      if (tgt_ptr->first_row <= row 
                          && row< tgt_ptr->first_row + tgt_ptr->nrows ) {
                        break;
                      }
                    } while( ++tgt_ptr<_blocks+nblocks() ); 

                    bassert( tgt_ptr < _blocks+nblocks() );

                    int_t lr = std::min(cur_src_lr,tgt_ptr->first_row + tgt_ptr->nrows-1);
                    int_t tgtOffset = tgt_ptr->_offset + (row - tgt_ptr->first_row)*tgt_snode_size;

                    for(int_t cr = row ;cr<=lr;++cr){
                      if(cr<=tgt_lc){
                        tmpBuffers.src_colindx[colidx++] = cr;
                      }
                      offset+=tgt_width;
                      tmpBuffers.src_to_tgt_offset[rowidx] = tgtOffset + (cr - row)*tgt_snode_size;
                      rowidx++;
                    }
                    row += (lr-row+1);
                  }
                }
                //Multiple cases to consider
                SYMPACK_TIMER_SPECIAL_STOP(UPDATE_SNODE_INDEX_MAP);

                SYMPACK_TIMER_SPECIAL_START(UPDATE_SNODE_ADD);
                if(first_pivot_idx==last_pivot_idx){
                  // Updating contiguous columns
                  rowind_t tgt_offset = (tgt_fc - this->first_col);
                  for (rowind_t rowidx = 0; rowidx < src_nrows; ++rowidx) {
                    T * A = &buf[rowidx*tgt_width];
                    T * B = &tgt[tmpBuffers.src_to_tgt_offset[rowidx] + tgt_offset];
#pragma unroll
                    for(rowind_t i = 0; i < tgt_width; ++i){ B[i] += A[i]; }
                  }
                }
                else{
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
    };



  template <typename colptr_t, typename rowind_t, typename T> 
    class symPACKMatrix2D: public symPACKMatrixMeta<T>{
      using snodeBlock_t = blockCell_t<colptr_t,rowind_t, T>;
      using SparseTask2D = scheduling::task_t<scheduling::meta_t, scheduling::depend_t>;

      public:
        int nsuper;
      Int coord2supidx(Int i, Int j){ 
        Int val =  (i-j) + j*(nsuper+1) - (j+1)*(j)/2;  
        bassert(val>=0);
        return val;
      };
      //TODO implement this
      //auto idx2coord = [nsuper](Int idx){ return 0; };

        SparseTask2D * task_lookup(Int cellI_, Int cellJ_, Int I_, Int J_, Int K_, Factorization::op_type type_){
          static std::unordered_map< Factorization::op_type , SparseTask2D * > last_ptr;


          //first check last_ptr
          auto last_it = last_ptr.find(type_);
          if ( last_it != last_ptr.end() ) {
            //auto meta = reinterpret_cast<SparseTask2D::meta_t*>( last_it->second->meta.data() );
            auto meta = &last_it->second->_meta;
            auto & src_snode = std::get<0>(meta[0]);
            auto & tgt_snode = std::get<1>(meta[0]);
            auto & facing_row = std::get<4>(meta[0]);
            auto K = this->SupMembership_[facing_row-1];
            auto & type = std::get<2>(meta[0]);

            bassert(type==type_);

            if( (src_snode == I_) && (tgt_snode == J_) && ( K == K_) ){
              return last_it->second;
            }
          }

          auto idx = coord2supidx(cellI_-1,cellJ_-1);
          auto it = std::find_if(this->task_graph[idx].begin(),this->task_graph[idx].end(),
              [I_,J_,K_,type_,this]( SparseTask2D * ptask)->bool{
              auto meta = &ptask->_meta;
              //auto meta = reinterpret_cast<SparseTask2D::meta_t*>(ptask->meta.data());
              auto & src_snode = std::get<0>(meta[0]);
              auto & tgt_snode = std::get<1>(meta[0]);
              auto & facing_row = std::get<4>(meta[0]);
              auto K = this->SupMembership_[facing_row-1];
              auto & type = std::get<2>(meta[0]);
              return (src_snode == I_) && (tgt_snode == J_) && (K == K_) && (type == type_) ;
              });

          if (it!=this->task_graph[idx].end()){
            //backup
            last_ptr[type_] = *it;
            return *it;
          }
          else{
            return (SparseTask2D*)nullptr;
          }
        };






      protected:
        //std::unordered_map<Int, snodeBlock_t > cells_;
        std::unordered_map<Int, std::shared_ptr<blockCellBase_t> > cells_;


        //TODO ideally, this should be DistSparseMatrixGraph<colptr_t,rowind_t>
#ifdef SP_THREADS
        std::map<std::thread::id,TempUpdateBuffers<T> > tmpBufs_th;
#else
        TempUpdateBuffers<T> tmpBufs;
#endif


#ifndef NO_MPI
#endif
      public:
        using TaskGraph2D = scheduling::task_graph_t<Int,SparseTask2D>; 
        TaskGraph2D task_graph;
        scheduling::Scheduler2D<SparseTask2D,TaskGraph2D> scheduler;


        symPACKMatrix2D();
        ~symPACKMatrix2D();



        void Init(symPACKOptions & options );
        void Init(DistSparseMatrix<T> & pMat, symPACKOptions & options );

        void SymbolicFactorization(DistSparseMatrix<T> & pMat);
        void DistributeMatrix(DistSparseMatrix<T> & pMat);


        void Factorize();



    };

  template <typename colptr_t, typename rowind_t, typename T>
    symPACKMatrix2D<colptr_t,rowind_t,T>::symPACKMatrix2D():
      symPACKMatrixMeta<T>()
  {
#ifndef NO_MPI
#endif
    g_sp_handle_to_matrix[this->sp_handle] = this;
  }

  template <typename colptr_t, typename rowind_t, typename T>
    symPACKMatrix2D<colptr_t,rowind_t,T>::~symPACKMatrix2D(){
    }

  template <typename colptr_t, typename rowind_t, typename T>
    void symPACKMatrix2D<colptr_t,rowind_t,T>::Init(symPACKOptions & options ){
      scope_timer(a,symPACKMatrix2D::Init);
      this->options_ = options;
      logfileptr->verbose = this->options_.verbose>0;

#ifndef NO_MPI
      this->all_np = 0;
      MPI_Comm_size(this->options_.MPIcomm,&this->all_np);
      MPI_Comm_rank(this->options_.MPIcomm,&this->iam);
#endif

#ifdef NEW_UPCXX
      this->iam = upcxx::rank_me();
      this->all_np = upcxx::rank_n();
#endif

      this->np = this->options_.used_procs(this->all_np);


#ifndef NO_MPI
      MPI_Comm_split(this->options_.MPIcomm,this->iam<this->np,this->iam,&this->workcomm_);

      //do another split to contain P0 and all the non working processors
      if(this->all_np!=this->np){
        this->non_workcomm_ = MPI_COMM_NULL;
        MPI_Comm_split(this->options_.MPIcomm,this->iam==0||this->iam>=this->np,this->iam==0?0:this->iam-this->np+1,&this->non_workcomm_);
      }
#else
      //    upcxx::team_all.split(this->iam<this->np,new_rank, this->team_);
#endif

    } 

  template <typename colptr_t, typename rowind_t, typename T>
    void symPACKMatrix2D<colptr_t,rowind_t,T>::Init(DistSparseMatrix<T> & pMat,symPACKOptions & options ){
    } 

  template <typename colptr_t, typename rowind_t, typename T>
    void symPACKMatrix2D<colptr_t,rowind_t,T>::SymbolicFactorization(DistSparseMatrix<T> & pMat ){
      scope_timer(a,symPACKMatrix2D::SymbolicFactorization);

      //This has to be declared here to be able to debug ...
      std::vector<int, Mallocator<int> > xadj;
      std::vector<int, Mallocator<int> > adj;
      Idx row;
      int fc,lc,colbeg,colend,col;
      Ptr supbeg,supend,rowidx;
      Int I;

#ifndef NO_MPI
      if(this->fullcomm_!=MPI_COMM_NULL){
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

      SparseMatrixGraph * sgraph = NULL;
      {
        {
          double timeSta = get_time();
          scope_timer(c,symPACKMatrix2D::ordering);

          if(this->options_.orderingStr=="MMD"){
            this->options_.ordering = symPACK::MMD;
          }
          else if(this->options_.orderingStr=="RCM"){
            this->options_.ordering = symPACK::RCM;
          }
          else if(this->options_.orderingStr=="AMD"){
            this->options_.ordering = symPACK::AMD;
          }
#ifdef USE_METIS
          else if(this->options_.orderingStr=="METIS"){
            this->options_.ordering = symPACK::METIS;
          }
#endif
#ifdef USE_SCOTCH
          else if(this->options_.orderingStr=="SCOTCH"){
            this->options_.ordering = symPACK::SCOTCH;
          }
#endif
#ifdef USE_PARMETIS
          else if(this->options_.orderingStr=="PARMETIS"){
            this->options_.ordering = symPACK::PARMETIS;
          }
#endif
#ifdef USE_PTSCOTCH
          else if(this->options_.orderingStr=="PTSCOTCH"){
            this->options_.ordering = symPACK::PTSCOTCH;
          }
#endif
          else if(this->options_.orderingStr=="NATURAL"){
            this->options_.ordering = symPACK::NATURAL;
          }
          else if(this->options_.orderingStr=="USER"){
            if(this->options_.perm ==NULL){
              throw std::logic_error( "When using USER, symPACKOptions.perm must be provided.\n" );
            }
            else{
              this->options_.ordering = symPACK::USER;
            }
          }
          else{
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
          switch(this->options_.ordering){
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
                for(int i=0;i<this->Order_.perm.size();i++){this->Order_.perm[i] = this->options_.perm[i] - *baseval + 1; /*1 based*/}
                this->Order_.invp.resize(this->iSize_);
                for(int i=0;i<this->Order_.perm.size();i++){this->Order_.invp[this->Order_.perm[i]-1] = i;} 
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
              for(int i=0;i<this->Order_.perm.size();i++){this->Order_.perm[i]=i+1;} 
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
          if(this->iam==0 && this->options_.verbose){
            std::cout<<"Ordering time: "<<timeStop - timeSta<<std::endl;
          }
          logfileptr->OFS()<<"Ordering done"<<std::endl;
        }
      }

      if(this->options_.dumpPerm>0){
        logfileptr->OFS()<<"perm = [";
        for (auto i : this->Order_.perm){ 
          logfileptr->OFS()<<i<<" ";
        }
        logfileptr->OFS()<<"]"<<std::endl;
      }

      std::vector<Int> cc;
      std::vector<Int> rc;


      if(sgraph==NULL){ 
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
        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Column count (distributed) construction time: "<<timeStop_cc - timeSta_cc<<std::endl;
        }

        if(this->options_.ordering != NATURAL){
          this->ETree_.SortChildren(cc,this->Order_);
          if(this->options_.dumpPerm>0){
            logfileptr->OFS()<<"perm = "<<this->Order_.perm<<std::endl;
          }
        }

        double timeStop = get_time();
        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Elimination tree construction time: "<<timeStop - timeSta<<std::endl;
        }
      }
      else
      {
        //gather the graph if necessary to build the elimination tree
        if(sgraph==NULL){ 
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
        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Column count (gather + serial + bcast) construction time: "<<timeStop_cc - timeSta_cc<<std::endl;
        }

        if(this->options_.ordering != NATURAL){
          this->ETree_.SortChildren(cc,this->Order_);
          if(this->options_.dumpPerm>0){
            logfileptr->OFS()<<"perm = "<<this->Order_.perm<<std::endl;
          }
        }

        double timeStop = get_time();
        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Elimination tree construction time: "<<timeStop - timeSta<<std::endl;
        }
      }

      //compute some statistics
      if(this->iam==0 && this->options_.verbose){
        double flops = 0.0;
        int64_t NNZ = 0;
        for(Int i = 0; i<cc.size();++i){
          flops+= (double)pow((double)cc[i],2.0);
          NNZ+=cc[i];
        }
        std::cout<<"Flops: "<<flops<<std::endl;
        std::cout<<"NNZ in L factor: "<<NNZ<<std::endl;
      }

      //get rid of the sequential graph
      if(sgraph!=NULL && this->options_.order_refinement_str.substr(0,4) != "TSPB"){delete sgraph;}

      { 
        double timeSta = get_time();
        this->findSupernodes(this->ETree_,this->Order_,cc,this->SupMembership_,this->Xsuper_,this->options_.relax.maxSize);
        logfileptr->OFS()<<"Supernodes found"<<std::endl;

        if(this->options_.relax.nrelax0>0)
        {
          this->relaxSupernodes(this->ETree_, cc,this->SupMembership_, this->Xsuper_, this->options_.relax );
          logfileptr->OFS()<<"Relaxation done"<<std::endl;
        }

        //modify this->np since it cannot be greater than the number of supernodes
        this->np = std::min(this->np,this->options_.used_procs(this->Xsuper_.size()-1)); 

        if(this->workcomm_!=MPI_COMM_NULL){
          MPI_Comm_free(&this->workcomm_);
        }
        MPI_Comm_split(this->options_.MPIcomm,this->iam<this->np,this->iam,&this->workcomm_);
        this->group_.reset( new RankGroup( this->workcomm_ ) ); 

        //do another split to contain P0 and all the non working processors
        if(this->all_np!=this->np){
          if(this->non_workcomm_!=MPI_COMM_NULL){
            MPI_Comm_free(&this->non_workcomm_);
          }
          this->non_workcomm_ = MPI_COMM_NULL;
          MPI_Comm_split(this->options_.MPIcomm,this->iam==0||this->iam>=this->np,this->iam==0?0:this->iam-this->np+1,&this->non_workcomm_);
        }

        //Compute this->XsuperDist_
        std::vector<Idx> newVertexDist;
        {
          Idx supPerProc = std::max((size_t)1,(this->Xsuper_.size()-1) / this->np);
          this->XsuperDist_.resize(this->all_np+1,0);
          this->XsuperDist_[this->np] = this->Xsuper_.size();
          for(int p = 0; p<this->np; p++){
            this->XsuperDist_[p]= std::min(this->XsuperDist_[this->np],(Int)(p*supPerProc+1));
          }
          for(int p = this->np+1; p<this->all_np; p++){
            this->XsuperDist_[p]= this->XsuperDist_[p-1];
          }
          this->XsuperDist_[this->all_np] = this->Xsuper_.size();


          newVertexDist.resize(this->all_np+1,0);
          newVertexDist[this->all_np] = this->iSize_+1;
          for(int p = 0; p < this->all_np; p++){
            Int S = this->XsuperDist_[p];
            newVertexDist[p] = this->Xsuper_[S-1];
          }
        }

        //      //now create the mapping
        //      if(this->Mapping_ != NULL){
        //        delete this->Mapping_;
        //      }
        //
        //      Int pmapping = this->np;
        //      Int pr = (Int)sqrt((double)pmapping);
        //      if(this->options_.mappingTypeStr ==  "MODWRAP2DTREE"){
        //        this->Mapping_ = new ModWrap2DTreeMapping(pmapping, pr, pr, 1);
        //      }
        //      else if(this->options_.mappingTypeStr ==  "WRAP2D"){
        //        this->Mapping_ = new Wrap2D(pmapping, pr, pr, 1);
        //      }
        //      else if(this->options_.mappingTypeStr ==  "WRAP2DFORCED"){
        //        this->Mapping_ = new Wrap2DForced(pmapping, pr, pr, 1);
        //      }
        //      else if(this->options_.mappingTypeStr ==  "MODWRAP2DNS"){
        //        this->Mapping_ = new Modwrap2DNS(pmapping, pr, pr, 1);
        //      }
        //      else if(this->options_.mappingTypeStr ==  "ROW2D"){
        //        this->Mapping_ = new Row2D(pmapping, pmapping, pmapping, 1);
        //      }
        //      else if(this->options_.mappingTypeStr ==  "COL2D"){
        //        this->Mapping_ = new Col2D(pmapping, pmapping, pmapping, 1);
        //      }
        //      else{
        //        this->Mapping_ = new Modwrap2D(pmapping, pr, pr, 1);
        //      }


        double timeStaSymb = get_time();
        this->symbolicFactorizationRelaxedDist(cc);

        double timeStopSymb = get_time();
        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Symbolic factorization time: "<<timeStopSymb - timeStaSymb<<std::endl;
        }
        logfileptr->OFS()<<"Symbfact done"<<std::endl;


        if(this->options_.order_refinement_str != "NO") {
          double timeSta = get_time();
          if(this->options_.order_refinement_str == "SET"){ 
            this->refineSupernodes(3,1,&pMat);
          }
          else if(this->options_.order_refinement_str == "SET10"){ 
            this->refineSupernodes(1,0,&pMat);
          }
          else if(this->options_.order_refinement_str == "SET11"){ 
            this->refineSupernodes(1,1,&pMat);
          }
          else if(this->options_.order_refinement_str == "SET20"){ 
            this->refineSupernodes(2,0,&pMat);
          }
          else if(this->options_.order_refinement_str == "SET21"){ 
            this->refineSupernodes(2,1,&pMat);
          }
          else if(this->options_.order_refinement_str == "SET30"){ 
            this->refineSupernodes(3,0,&pMat);
          }
          else if(this->options_.order_refinement_str == "SET31"){ 
            this->refineSupernodes(3,1,&pMat);
          }
          else if(this->options_.order_refinement_str == "TSP"){
            auto SupETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);

            std::vector<Ptr> xlindx;
            std::vector<Idx> lindx;
            this->gatherLStructure(xlindx, lindx);

            if(this->iam==0){
              TSP::SymbolMatrix * symbmtx = TSP::GetPastixSymbolMatrix(this->Xsuper_,this->SupMembership_, xlindx, lindx);
              TSP::Order * psorder = TSP::GetPastixOrder(symbmtx,this->Xsuper_, SupETree, &this->Order_.perm[0], &this->Order_.invp[0]);

              double timeSta = get_time();
              TSP::symbolReordering( symbmtx, psorder, 0, std::numeric_limits<int>::max(), 0 );
              double timeStop = get_time();
              if(this->iam==0 && this->options_.verbose){
                std::cout<<"TSP reordering done in "<<timeStop-timeSta<<std::endl;
              }

              //overwrite order
              for(int i = 0; i < this->Order_.perm.size(); ++i){
                this->Order_.perm[i] = psorder->peritab[i]+1;
                this->Order_.invp[i] = psorder->permtab[i]+1;
              }
            }

            MPI_Bcast(this->Order_.perm.data(),this->Order_.perm.size()*sizeof(Int),MPI_BYTE,0,this->fullcomm_);
            MPI_Bcast(this->Order_.invp.data(),this->Order_.invp.size()*sizeof(Int),MPI_BYTE,0,this->fullcomm_);
          } 
          else if(this->options_.order_refinement_str.substr(0,4) == "TSPB"){

            if(sgraph==NULL){ 
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

              if(this->iam==0){
                ixlindx.resize(xlindx.size());
                for(int i = 0;i<xlindx.size();i++){
                  ixlindx[i] = xlindx[i];
                }
                ilindx.resize(lindx.size());
                for(int i = 0;i<lindx.size();i++){
                  ilindx[i] = lindx[i];
                }
              }
            }

            if(this->iam==0){
              bassert(sgraph!=nullptr);
              int neqns = this->iSize_;

              int nofsub =ilindx.size();
              nsuper = xsuper.size()-1;


              new_invp.assign(neqns,0);


              std::vector<int>  new_perm(neqns,0);

              for(size_t i =0; i<new_perm.size(); i++){ new_invp[i] = this->Order_.invp[i];}
              for(size_t i =0; i<new_perm.size(); i++){ new_perm[i] = this->Order_.perm[i];}

              int supsiz = 0;
              for(I = 1; I <= nsuper; I++){
                Int fc = this->Xsuper_[I-1];
                Int lc = this->Xsuper_[I]-1;
                supsiz = std::max(supsiz,lc-fc+1);
              }

              sgraph->SetKeepDiag(0);
              xadj.resize(sgraph->colptr.size());
              adj.resize(sgraph->rowind.size());
              int nadj = adj.size();
              for(size_t i = 0; i< sgraph->colptr.size(); i++){ xadj[i] = int(sgraph->colptr[i]); }
              for(size_t i = 0; i< sgraph->rowind.size(); i++){ adj[i] = int(sgraph->rowind[i]); }

              if(sgraph!=NULL){delete sgraph;}

              std::vector<int, Mallocator<int> > etpar(neqns);
              for(int i = 0; i<neqns; i++){ etpar[i] = tree.PostParent(i); }


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
              if(this->options_.order_refinement_str == "TSPB"){
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
              else{
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
              if(this->iam==0 && this->options_.verbose){
                std::cout<<"TSPB reordering done in "<<timeStop-timeSta<<std::endl;
              }


              //this->Order_.Compose(new_invp);
              for(size_t i =0; i<new_perm.size(); i++){ this->Order_.invp[i] = new_invp[i];}
              for(size_t i =0; i<new_perm.size(); i++){ this->Order_.perm[i] = new_perm[i];}
            }

            // broadcast invp
            Int N = aOrder.invp.size();
            MPI_Bcast(&aOrder.invp[0],N*sizeof(Int),MPI_BYTE,0,this->fullcomm_);
            MPI_Bcast(&aOrder.perm[0],N*sizeof(Int),MPI_BYTE,0,this->fullcomm_);

          }

          double timeStop = get_time();

          if(this->iam==0 && this->options_.verbose){
            std::cout<<"Supernode reordering done in "<<timeStop-timeSta<<std::endl;
          }

          {
            double timeSta = get_time();
            this->symbolicFactorizationRelaxedDist(cc);
            double timeStop = get_time();
            if(this->iam==0 && this->options_.verbose){
              std::cout<<"Symbolic factorization time: "<<timeStop - timeSta<<std::endl;
            }
            logfileptr->OFS()<<"Symbfact done"<<std::endl;
          }
        }

        double timeStop = get_time();
        if(this->iam==0 && this->options_.verbose){
          std::cout<<"Total symbolic factorization time: "<<timeStop - timeSta<<std::endl;
        }

        //Print statistics
        if(this->options_.print_stats){
          OrderStats stats;
          stats.get(this->Xsuper_, this->XsuperDist_, this->locXlindx_, this->locLindx_, this->fullcomm_);
          if (this->iam==0){
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

      //compute a mapping
      // Do a 2D cell-cyclic distribution for now.
      nsuper = this->Xsuper_.size()-1;
      auto iam = this->iam;
      auto np = this->np;
      Int npcol = std::floor(std::sqrt(np));
      Int nprow = std::floor(np / npcol);


#if 1    
      {
        using cell_tuple_t = std::tuple<Idx,Idx>;
        std::list< cell_tuple_t> cellsToSend;

        Int numLocSnode = this->XsuperDist_[this->iam+1]-this->XsuperDist_[this->iam];
        Int firstSnode = this->XsuperDist_[this->iam];
        for(Int locsupno = 1; locsupno<this->locXlindx_.size(); ++locsupno){
          Idx I = locsupno + firstSnode-1;
          Int first_col = this->Xsuper_[I-1];
          Int last_col = this->Xsuper_[I]-1;

          Ptr lfi = this->locXlindx_[locsupno-1];
          Ptr lli = this->locXlindx_[locsupno]-1;

          Int K = -1; 
          Idx K_prevSnode = -1;
          for(Ptr K_sidx = lfi; K_sidx<=lli;K_sidx++){
            Idx K_row = this->locLindx_[K_sidx-1];
            K = this->SupMembership_[K_row-1];

            //Split at boundary or after diagonal block
            if(K!=K_prevSnode){
              if(K>=I){
                cellsToSend.push_back(std::make_tuple(K,I));
              }
            }
            K_prevSnode = K;
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

        for(auto it = recvbuf.begin();it!=recvbuf.end();it++){
          auto & cur_cell = (*it);
          Int i = std::get<0>(cur_cell);
          Int j = std::get<1>(cur_cell);
          auto pcol = j % npcol;
          //auto idx = coord2supidx(i,j);
          auto idx = coord2supidx(i-1,j-1);
          auto prow = i % nprow;
          auto p = pcol*nprow+prow;
          //cells_[idx].owner = p;
          //cells_[idx].i = i;
          //cells_[idx].j = j;
          auto sptr = std::make_shared<snodeBlock_t>();
          cells_[idx] = std::static_pointer_cast<blockCellBase_t>(sptr);
          //logfileptr->OFS()<<idx<<" "<<1<<" "<<1<<" = "<<cells_[1]<<std::endl;
          //logfileptr->OFS()<<i<<" "<<j<<" = "<<sptr<<" "<<std::static_pointer_cast<blockCellBase_t>(sptr)<<" "<<cells_[idx]<<std::endl;
          sptr->owner = p;
          sptr->i = i;
          sptr->j = j;
        }
        //logfileptr->OFS()<<1<<" "<<1<<" = "<<cells_[1]<<std::endl;
      }

#else

      cells_.resize(nsuper*(nsuper-1)/2.0+nsuper);

      //should create cell only if non empty
      for(Int j = 0; j < nsuper; j++){
        auto pcol = j % npcol;
        for(Int i = j; i < nsuper; i++){
          auto idx = coord2supidx(i,j);
          auto prow = i % nprow;
          auto p = pcol*nprow+prow;
          //cells_[idx].owner = p;
          //cells_[idx].i = i;
          //cells_[idx].j = j;
          auto sptr = std::make_shared<snodeBlock_t>();
          cells_[idx] = std::static_pointer_cast<blockCellBase_t>(sptr);
          sptr->owner = p;
          sptr->i = i;
          sptr->j = j;

        }
      }
#endif

      if(iam==0){
        std::cout<<"#supernodes = "<<nsuper<<" #cells = "<<cells_.size()<< " which is "<<cells_.size()*sizeof(snodeBlock_t)<<" bytes"<<std::endl;
      }

      //generate task graph
      {
        auto supETree = this->ETree_.ToSupernodalETree(this->Xsuper_,this->SupMembership_,this->Order_);

        //we will need to communicate if only partial xlindx_, lindx_
        //idea: build tasklist per processor and then exchange
        //tuple is: src_snode,tgt_snode,op_type,lt_first_row = updated fc, facing_first_row = updated fr,
        std::map<Idx, std::list< SparseTask2D::meta_t> > Updates;
        std::vector<int> marker(this->np,0);

        std::vector< std::map< Idx, std::pair<Idx,std::set<Idx> > > > updCnt(this->Xsuper_.size() -1 );

        Int numLocSnode = this->XsuperDist_[this->iam+1]-this->XsuperDist_[this->iam];
        Int firstSnode = this->XsuperDist_[this->iam];
        for(Int locsupno = 1; locsupno<this->locXlindx_.size(); ++locsupno){
          Idx I = locsupno + firstSnode-1;
          Int first_col = this->Xsuper_[I-1];
          Int last_col = this->Xsuper_[I]-1;

          //logfileptr->OFS()<<I<<" "<<I<<" = "<<cells_[coord2supidx(I-1,I-1)]<<std::endl;
          bassert( cells_[coord2supidx(I-1,I-1)] != nullptr );

          auto & fcell = CELL(I-1,I-1);//*std::static_pointer_cast<snodeBlock_t>(cells_[coord2supidx(I-1,I-1)]);
          //auto & fcell = cells_[coord2supidx(I-1,I-1)];
          Int iOwner = fcell.owner;

          //Create the factor task on the owner
          Updates[iOwner].push_back(std::make_tuple(I,I,Factorization::op_type::FACTOR,first_col,first_col));


          Ptr lfi = this->locXlindx_[locsupno-1];
          Ptr lli = this->locXlindx_[locsupno]-1;
          std::list<Int> ancestor_rows;

          Int K = -1; 
          Idx K_prevSnode = -1;
          for(Ptr K_sidx = lfi; K_sidx<=lli;K_sidx++){
            Idx K_row = this->locLindx_[K_sidx-1];
            K = this->SupMembership_[K_row-1];

            //Split at boundary or after diagonal block
            if(K!=K_prevSnode){
              if(K>I){
                ancestor_rows.push_back(K_row);
              }
            }
            K_prevSnode = K;
          }

          for(auto J_row : ancestor_rows){
            Int J = this->SupMembership_[J_row-1];

            auto & fodcell = CELL(J-1,I-1);//*std::static_pointer_cast<snodeBlock_t>(cells_[coord2supidx(J-1,I-1)]);
            //auto & fodcell = cells_[coord2supidx(J-1,I-1)];
            Int iFODOwner = fodcell.owner;

            //          if(iOwner!=iFODOwner){
            Updates[iOwner].push_back(std::make_tuple(I,I,Factorization::op_type::TRSM_SEND,J_row,J_row));
            Updates[iFODOwner].push_back(std::make_tuple(I,I,Factorization::op_type::TRSM_RECV,J_row,J_row));
            //          }

            Updates[iFODOwner].push_back(std::make_tuple(I,I,Factorization::op_type::TRSM,J_row,J_row));

            //add the UPDATE_SEND from FODOwner tasks
            for(auto K_row : ancestor_rows){
              K = this->SupMembership_[K_row-1];

              if(K>=J){
                auto & tgtupdcell = CELL(K-1,J-1);//*std::static_pointer_cast<snodeBlock_t>(cells_[coord2supidx(K-1,J-1)]);
                //auto & tgtupdcell = cells_[coord2supidx(K-1,J-1)];
                Int iTgtOwner = tgtupdcell.owner;
                //TODO update this for non fan-out mapping
                Int iUpdOwner = tgtupdcell.owner;

                //              if(iFODOwner!=iUpdOwner){
                //cell(J,I) to cell(K,J)
                Updates[iFODOwner].push_back(std::make_tuple(I,J,Factorization::op_type::UPDATE2D_SEND_OD,K_row,J_row));
                Updates[iUpdOwner].push_back(std::make_tuple(I,J,Factorization::op_type::UPDATE2D_RECV_OD,K_row,J_row));
                //              }
              }

              if(K<=J){
                auto & tgtupdcell = CELL(J-1,K-1);//*std::static_pointer_cast<snodeBlock_t>(cells_[coord2supidx(J-1,K-1)]);
                //auto & tgtupdcell = cells_[coord2supidx(J-1,K-1)];
                Int iTgtOwner = tgtupdcell.owner;
                //TODO update this for non fan-out mapping
                Int iUpdOwner = tgtupdcell.owner;

                auto & facingcell = CELL(J-1,I-1);//*std::static_pointer_cast<snodeBlock_t>(cells_[coord2supidx(J-1,I-1)]);
                //auto & facingcell = cells_[coord2supidx(J-1,I-1)];
                Int iFacingOwner = facingcell.owner;

                //sender point of view
                //cell(J,I) to cell(J,K)
                //              if(iFacingOwner!=iUpdOwner){
                Updates[iFacingOwner].push_back(std::make_tuple(I,K,Factorization::op_type::UPDATE2D_SEND,K_row,J_row));
                Updates[iUpdOwner].push_back(std::make_tuple(I,K,Factorization::op_type::UPDATE2D_RECV,K_row,J_row));
                //              }
              }
            }

            for(auto K_row : ancestor_rows){
              K = this->SupMembership_[K_row-1];
              if(K>=J){
                auto & tgtupdcell = CELL(K-1,J-1);//*std::static_pointer_cast<snodeBlock_t>(cells_[coord2supidx(K-1,J-1)]);
                //auto & tgtupdcell = cells_[coord2supidx(K-1,J-1)];
                Int iTgtOwner = tgtupdcell.owner;
                //TODO update this for non fan-out mapping
                Int iUpdOwner = tgtupdcell.owner;

                //update on cell(K,J)
                Updates[iUpdOwner].push_back(std::make_tuple(I,J,Factorization::op_type::UPDATE2D_COMP,J_row,K_row));

                //              if(iUpdOwner!=iTgtOwner){
                //cell(K,J) to cell (K,J)
                Updates[iUpdOwner].push_back(std::make_tuple(I,J,Factorization::op_type::AGGREGATE2D_SEND,J_row,K_row));
                Updates[iTgtOwner].push_back(std::make_tuple(I,J,Factorization::op_type::AGGREGATE2D_RECV,J_row,K_row));
                //              }
              }
            }
          }
        }

        MPI_Datatype type;
        MPI_Type_contiguous( sizeof(SparseTask2D::meta_t), MPI_BYTE, &type );
        MPI_Type_commit(&type);

        //then do an alltoallv
        //compute send sizes
        vector<int> ssizes(this->np,0);
        for(auto itp = Updates.begin();itp!=Updates.end();itp++){
          ssizes[itp->first] = itp->second.size();//*sizeof(std::pair<Idx,Idx>);
        }

        //compute send displacements
        vector<int> sdispls(this->np+1,0);
        sdispls[0] = 0;
        std::partial_sum(ssizes.begin(),ssizes.end(),&sdispls[1]);

        //Build the contiguous array of pairs
        vector<SparseTask2D::meta_t> sendbuf;
        sendbuf.reserve(sdispls.back());

        for(auto itp = Updates.begin();itp!=Updates.end();itp++){
          sendbuf.insert(sendbuf.end(),itp->second.begin(),itp->second.end());
        }

        //gather receive sizes
        vector<int> rsizes(this->np,0);
        MPI_Alltoall(&ssizes[0],sizeof(int),MPI_BYTE,&rsizes[0],sizeof(int),MPI_BYTE,this->workcomm_);

        //compute receive displacements
        vector<int> rdispls(this->np+1,0);
        rdispls[0] = 0;
        std::partial_sum(rsizes.begin(),rsizes.end(),&rdispls[1]);

        //Now do the alltoallv
        vector<SparseTask2D::meta_t> recvbuf(rdispls.back());///sizeof(std::pair<Idx,Idx>));
        MPI_Alltoallv(&sendbuf[0],&ssizes[0],&sdispls[0],type,&recvbuf[0],&rsizes[0],&rdispls[0],type,this->workcomm_);

        MPI_Type_free(&type);


        //do a top down traversal of my local task list
        for(auto it = recvbuf.begin();it!=recvbuf.end();it++){
          auto & cur_op = (*it);
          auto & src_snode = std::get<0>(cur_op);
          auto & tgt_snode = std::get<1>(cur_op);
          auto & type = std::get<2>(cur_op);
          auto & lt_first_row = std::get<3>(cur_op);
          auto & facing_first_row = std::get<4>(cur_op);

          SparseTask2D * ptask = nullptr;
          switch(type){
            case Factorization::op_type::TRSM:
            case Factorization::op_type::FACTOR:
            case Factorization::op_type::UPDATE2D_COMP:
              {
                ptask = new SparseTask2D;
                size_t tuple_size = sizeof(SparseTask2D::meta_t);
                //ptask->meta.resize(tuple_size);
                //auto meta = reinterpret_cast<SparseTask2D::meta_t*>(ptask->meta.data());
                auto meta = &ptask->_meta;
                meta[0] = *it;
                //backup for convenience
                //ptask->_meta = meta;
              }
              break;
            default:
              break;
          }


          auto I = src_snode;
          auto J = tgt_snode;
          switch(type){
            case Factorization::op_type::TRSM:
              {
                auto J = this->SupMembership_[facing_first_row-1];
                logfileptr->OFS()<<"TRSM"<<" from "<<I<<" to ("<<I<<") cell ("<<J<<","<<I<<")"<<std::endl;
                this->task_graph[coord2supidx(J-1,I-1)].push_back(ptask);
              }
              break;
            case Factorization::op_type::FACTOR:
              {
                logfileptr->OFS()<<"FACTOR"<<" cell ("<<J<<","<<I<<")"<<std::endl;
                this->task_graph[coord2supidx(I-1,I-1)].push_back(ptask);
              }
              break;
            case Factorization::op_type::UPDATE2D_COMP:
              {
                auto K = this->SupMembership_[facing_first_row-1];
                logfileptr->OFS()<<"UPDATE"<<" from "<<I<<" to "<<J<<" facing cell ("<<K<<","<<I<<") and cell ("<<J<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;
                this->task_graph[coord2supidx(K-1,J-1)].push_back(ptask);
              }
              break;
            default:
              break;
          }
        }


        //process dependency tasks
        for(auto it = recvbuf.begin();it!=recvbuf.end();it++){
          auto & cur_op = (*it);
          auto & src_snode = std::get<0>(cur_op);
          auto & tgt_snode = std::get<1>(cur_op);
          auto & type = std::get<2>(cur_op);
          auto & lt_first_row = std::get<3>(cur_op);
          auto & facing_first_row = std::get<4>(cur_op);


          auto I = src_snode;
          auto J = tgt_snode;

          switch(type){
            case Factorization::op_type::TRSM_SEND:
              {
                auto J = this->SupMembership_[facing_first_row-1];

                auto taskptr = task_lookup(I,I,I,I,I,Factorization::op_type::FACTOR);
                bassert(taskptr!=nullptr); 
                logfileptr->OFS()<<"TRSM_SEND"<<" from "<<I<<" to "<<I<<" cell ("<<J<<","<<I<<")"<<std::endl;
#ifdef _DEBUG_DEPENDENCIES_
                taskptr->out_dependencies.push_back( std::make_tuple(J,I,false) );
#else
                taskptr->out_dependencies.push_back( std::make_tuple(J,I) );
#endif
                logfileptr->OFS()<<"        out dep added to FACTOR"<<" from "<<I<<" to "<<I<<std::endl;
              }
              break;
            case Factorization::op_type::TRSM_RECV:
              {
                auto J = this->SupMembership_[facing_first_row-1];

                //find the TRSM and add cell I,I as incoming dependency
                auto taskptr = task_lookup(J,I,I,I,J,Factorization::op_type::TRSM);
                bassert(taskptr!=nullptr); 
                logfileptr->OFS()<<"TRSM_RECV"<<" from "<<I<<" to "<<I<<" cell ("<<J<<","<<I<<")"<<std::endl;
#ifdef _DEBUG_DEPENDENCIES_
                if (  CELL(J-1,I-1).owner != CELL(I-1,I-1).owner )
                  taskptr->in_remote_dependencies.push_back( std::make_tuple(I,I,false) );
                else
                  taskptr->in_local_dependencies.push_back( std::make_tuple(I,I,false) );
#else
                if (  CELL(J-1,I-1).owner != CELL(I-1,I-1).owner )
                  taskptr->in_remote_dependencies.push_back( std::make_tuple(I,I) );
                else
                  taskptr->in_local_dependencies.push_back( std::make_tuple(I,I) );
#endif
                logfileptr->OFS()<<"        in dep added to TRSM"<<" from "<<I<<" to "<<I<<std::endl;
              }
              break;
            case Factorization::op_type::AGGREGATE2D_RECV:
              {
                auto K = this->SupMembership_[facing_first_row-1];

                logfileptr->OFS()<<"AGGREGATE2D_RECV"<<" from "<<I<<" to "<<J<<" cell ("<<K<<","<<J<<") to cell("<<K<<","<<J<<")"<<std::endl;
                //find the FACTOR or TRSM and add cell K,J as incoming dependency
                if(K==J){
                  auto taskptr = task_lookup(K,J,J,J,J,Factorization::op_type::FACTOR);
                  bassert(taskptr!=nullptr); 

#ifdef _DEBUG_DEPENDENCIES_
                  taskptr->in_local_dependencies.push_back( std::make_tuple(K,J,false) );
#else
                  taskptr->in_local_dependencies.push_back( std::make_tuple(K,J) );
#endif
                  logfileptr->OFS()<<"        in dep added to FACTOR"<<" from "<<J<<" to "<<J<<std::endl;
                }
                else{
                  auto taskptr = task_lookup(K,J,J,J,K,Factorization::op_type::TRSM);
                  bassert(taskptr!=nullptr); 
#ifdef _DEBUG_DEPENDENCIES_
                  taskptr->in_local_dependencies.push_back( std::make_tuple(K,J,false) );
#else
                  taskptr->in_local_dependencies.push_back( std::make_tuple(K,J) );
#endif
                  logfileptr->OFS()<<"        in dep added to TRSM"<<" from "<<J<<" to "<<J<<std::endl;
                }
              }
              break;
            case Factorization::op_type::AGGREGATE2D_SEND:
              {
                //TODO this might be useless
                auto K = this->SupMembership_[facing_first_row-1];

                logfileptr->OFS()<<"AGGREGATE2D_SEND"<<" from "<<I<<" to "<<J<<" cell ("<<K<<","<<J<<") to cell("<<K<<","<<J<<")"<<std::endl;
                //find the UPDATE2D_COMP and add cell K,J as outgoing dependency
                auto taskptr = task_lookup(K,J,I,J,K,Factorization::op_type::UPDATE2D_COMP);
                bassert(taskptr!=nullptr); 

#ifdef _DEBUG_DEPENDENCIES_
                taskptr->out_dependencies.push_back( std::make_tuple(K,J,false) );
#else
                taskptr->out_dependencies.push_back( std::make_tuple(K,J) );
#endif
                logfileptr->OFS()<<"        out dep added to UPDATE_COMP"<<" from "<<I<<" to "<<J<<std::endl;
              }
              break;

            case Factorization::op_type::UPDATE2D_RECV:
              {
                //TODO this might be useless
                Idx K = this->SupMembership_[facing_first_row-1];
                std::swap(K,J);
                logfileptr->OFS()<<"UPDATE_RECV"<<" from "<<I<<" to "<<K<<" cell ("<<J<<","<<I<<") to cell ("<<J<<","<<K<<")"<<std::endl;
                //find the UPDATE2D_COMP and add cell J,I as incoming dependency
                auto taskptr = task_lookup(J,K,I,K,J,Factorization::op_type::UPDATE2D_COMP);
                bassert(taskptr!=nullptr);
#ifdef _DEBUG_DEPENDENCIES_
                if (  CELL(J-1,I-1).owner != CELL(J-1,K-1).owner )
                  taskptr->in_remote_dependencies.push_back( std::make_tuple(J,I,false) );
                else
                  taskptr->in_local_dependencies.push_back( std::make_tuple(J,I,false) );
#else
                if (  CELL(J-1,I-1).owner != CELL(J-1,K-1).owner )
                  taskptr->in_remote_dependencies.push_back( std::make_tuple(J,I) );
                else
                  taskptr->in_local_dependencies.push_back( std::make_tuple(J,I) );
#endif

                logfileptr->OFS()<<"        in dep added to UPDATE_COMP"<<" from "<<I<<" to "<<K<<std::endl;
              }
              break;
            case Factorization::op_type::UPDATE2D_SEND:
              {
                Idx K = this->SupMembership_[facing_first_row-1];
                std::swap(K,J);
                logfileptr->OFS()<<"UPDATE_SEND"<<" from "<<I<<" to "<<K<<" cell ("<<J<<","<<I<<") to cell ("<<J<<","<<K<<")"<<std::endl;
                //find the TRSM and add cell J,K as outgoing dependency
                auto taskptr = task_lookup(J,I,I,I,J,Factorization::op_type::TRSM);
                bassert(taskptr!=nullptr); 
#ifdef _DEBUG_DEPENDENCIES_
                taskptr->out_dependencies.push_back( std::make_tuple(J,K,false) );
#else
                taskptr->out_dependencies.push_back( std::make_tuple(J,K) );
#endif
                logfileptr->OFS()<<"        out dep added to TRSM"<<" from "<<I<<" to "<<I<<std::endl;
              }
              break;
            case Factorization::op_type::UPDATE2D_SEND_OD:
              {
                auto K = this->SupMembership_[lt_first_row-1];
                logfileptr->OFS()<<"UPDATE_SEND OD"<<" from "<<I<<" to "<<J<<" cell ("<<J<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;

                //find the TRSM and add cell K,J as outgoing dependency
                auto taskptr = task_lookup(J,I,I,I,J,Factorization::op_type::TRSM);
                bassert(taskptr!=nullptr); 
#ifdef _DEBUG_DEPENDENCIES_
                taskptr->out_dependencies.push_back( std::make_tuple(K,J,false) );
#else
                taskptr->out_dependencies.push_back( std::make_tuple(K,J) );
#endif
                logfileptr->OFS()<<"        out dep DOWN added to TRSM"<<" from "<<I<<" to "<<I<<std::endl;
              }
              break;
            case Factorization::op_type::UPDATE2D_RECV_OD:
              {
                auto K = this->SupMembership_[lt_first_row-1];
                logfileptr->OFS()<<"UPDATE_RECV OD"<<" from "<<I<<" to "<<J<<" cell ("<<J<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;

                //find the UPDATE2D_COMP and add cell J,I as incoming dependency
                auto taskptr = task_lookup(K,J,I,J,K,Factorization::op_type::UPDATE2D_COMP);
                bassert(taskptr!=nullptr);
#ifdef _DEBUG_DEPENDENCIES_
                if (  CELL(J-1,I-1).owner != CELL(K-1,J-1).owner )
                  taskptr->in_remote_dependencies.push_back( std::make_tuple(J,I,false) );
                else
                  taskptr->in_local_dependencies.push_back( std::make_tuple(J,I,false) );
#else
                if (  CELL(J-1,I-1).owner != CELL(K-1,J-1).owner )
                  taskptr->in_remote_dependencies.push_back( std::make_tuple(J,I) );
                else
                  taskptr->in_local_dependencies.push_back( std::make_tuple(J,I) );
#endif

                logfileptr->OFS()<<"        in dep DOWN added to UPDATE2D_COMP"<<" from "<<I<<" to "<<J<<std::endl;
              }
              break;
            default:
              break;
          }

        }


        scheduler.sp_handle = this->sp_handle;

        //Now we have our local part of the task graph
        for(auto it = this->task_graph.begin(); it != this->task_graph.end(); it++){
          auto idx = it->first;
          auto tasklist = it->second;

          for( auto && ptask: tasklist){
            //auto meta = reinterpret_cast<SparseTask2D::meta_t*>(ptask->meta.data());
            auto meta = &ptask->_meta;
            auto & src_snode = std::get<0>(meta[0]);
            auto & tgt_snode = std::get<1>(meta[0]);
            auto & lt_first_row = std::get<3>(meta[0]);
            auto & facing_row = std::get<4>(meta[0]);
            auto & type = std::get<2>(meta[0]);

            auto I = src_snode;
            auto J = tgt_snode;
            auto K = this->SupMembership_[facing_row-1];

            switch(type){
              case Factorization::op_type::TRSM:
                {
                  logfileptr->OFS()<<"TRSM"<<" from "<<I<<" to ("<<I<<") cell ("<<K<<","<<I<<")"<<std::endl;
                }
                break;
              case Factorization::op_type::FACTOR:
                {
                  logfileptr->OFS()<<"FACTOR"<<" cell ("<<I<<","<<I<<")"<<std::endl;
                }
                break;
              case Factorization::op_type::UPDATE2D_COMP:
                {
                  logfileptr->OFS()<<"UPDATE"<<" from "<<I<<" to "<<J<<" facing cell ("<<K<<","<<I<<") and cell ("<<J<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;
                }
                break;
              default:
                delete ptask;
                break;
            }

            logfileptr->OFS()<<"      local input dependencies: ";
            for(auto &&tpl : ptask->in_local_dependencies){
              logfileptr->OFS()<<"("<<std::get<0>(tpl)<<","<<std::get<1>(tpl)<<") ";
            }
            logfileptr->OFS()<<std::endl;
            logfileptr->OFS()<<"      remote input dependencies: ";
            for(auto &&tpl : ptask->in_remote_dependencies){
              logfileptr->OFS()<<"("<<std::get<0>(tpl)<<","<<std::get<1>(tpl)<<") ";
            }
            logfileptr->OFS()<<std::endl;
            logfileptr->OFS()<<"      output dependencies: ";
            for(auto &&tpl : ptask->out_dependencies){
              logfileptr->OFS()<<"("<<std::get<0>(tpl)<<","<<std::get<1>(tpl)<<") ";
            }
            logfileptr->OFS()<<std::endl;

//            ptask->dep_count = ptask->in_dependencies.size();      

            auto remote_deps = ptask->in_remote_dependencies.size();
            auto local_deps = ptask->in_local_dependencies.size();

            ptask->in_prom->require_anonymous(local_deps + remote_deps);
            ptask->in_avail_prom->require_anonymous(remote_deps);

            if (remote_deps >0 ) {
              auto fut_comm = ptask->in_avail_prom->get_future();
              fut_comm.then([this,ptask](){
                  //TODO enqueue ptask somewhere
                  //gdb_lock();
                  this->scheduler.avail_tasks.push_back(ptask);
                  });
            }
            //fulfill the promise by one to "unlock" that promis
            ptask->in_avail_prom->fulfill_anonymous(1);

            //do not fulfill it to lock the tasks
            auto fut = ptask->in_prom->get_future();
            fut.then([this,ptask](){
                //TODO enqueue ptask somewhere
                //gdb_lock();
                this->scheduler.ready_tasks.push_back(ptask);
                     });
          }
        }

        for(auto it = this->task_graph.begin(); it != this->task_graph.end(); it++){
          auto idx = it->first;
          auto tasklist = it->second;

          for( auto && ptask: tasklist){
            //auto meta = reinterpret_cast<SparseTask2D::meta_t*>(ptask->meta.data());
            auto meta = &ptask->_meta;
            auto & src_snode = std::get<0>(meta[0]);
            auto & tgt_snode = std::get<1>(meta[0]);
            auto & lt_first_row = std::get<3>(meta[0]);
            auto & facing_row = std::get<4>(meta[0]);
            auto & type = std::get<2>(meta[0]);

            auto I = src_snode;
            auto J = tgt_snode;
            auto K = this->SupMembership_[facing_row-1];

            switch(type){
              case Factorization::op_type::FACTOR:
                {

                  ptask->execute = [this,src_snode,ptask,I,J,K] () {
                    scope_timer(b,FB_FACTOR_DIAG_TASK);
                    logfileptr->OFS()<<"Exec FACTOR"<<" cell ("<<I<<","<<I<<")"<<std::endl;

                    auto & diagcell = CELL(I-1,I-1);
                    bassert( diagcell.owner == this->iam);

#ifdef SP_THREADS
                    std::thread::id tid = std::this_thread::get_id();
                    auto & tmpBuf = tmpBufs_th[tid];
#else
                    auto & tmpBuf = tmpBufs;
#endif
                    //diagcell.factorize(tmpBuf);

                    std::unordered_map<Int,std::list<SparseTask2D::depend_task_t> > data_to_send;
                    for (auto &tpl : ptask->out_dependencies) {
                      auto tgt_i = std::get<0>(tpl);
                      auto tgt_j = std::get<1>(tpl);
                      auto & tgt_cell = CELL(tgt_i-1,tgt_j-1);
                      //add this cell to the list of cells depending on diagcell for the TRSM task(s)
                      data_to_send[tgt_cell.owner].push_back(tpl);
                    }

                    //send factor and update local tasks
                    //ptask->out_prom->require_anonymous(data_to_send.size());
                    for (auto it = data_to_send.begin(); it!=data_to_send.end(); it++) {
                      auto pdest = it->first;
                      auto & tgt_cells = it->second;
                      //serialize data once, and list of meta data
                      //factor is output data so it will not be deleted
                      if ( pdest != this->iam ) {
                        auto cxs = upcxx::source_cx::as_buffered() | upcxx::source_cx::as_promise(*ptask->out_prom);
#if 1
                        upcxx::rpc_ff( pdest, /*cxs,*/ 
                            [ ] (int sp_handle, upcxx::global_ptr<char> gptr, size_t storage_size, size_t nnz, size_t nblocks, rowind_t width, SparseTask2D::meta_t meta, upcxx::view<SparseTask2D::depend_task_t> target_cells ) { 
                            //gdb_lock();
                            //store pointer & associated metadata somewhere
                            auto data = std::make_shared<SparseTask2D::data_t >();
                            data->in_meta = meta;
                            data->size = storage_size;
                            data->remote_gptr = gptr;

                            //there is a map between sp_handle and task_graphs
                            auto matptr = (symPACKMatrix2D<colptr_t,rowind_t,T> *) g_sp_handle_to_matrix[sp_handle];
                            for ( auto & tgt_cell: target_cells) {
                              auto tgt_i = std::get<0>(tgt_cell);
                              auto tgt_j = std::get<1>(tgt_cell);
                              auto taskptr = matptr->task_lookup(tgt_i,tgt_j,tgt_j,tgt_j,tgt_i,Factorization::op_type::TRSM);
                              taskptr->input_msg.push_back(data);
                              data->target_tasks.push_back(taskptr);
                              taskptr->in_avail_prom->fulfill_anonymous(1);

#ifdef _DEBUG_DEPENDENCIES_
                              auto it = std::find(taskptr->in_remote_dependencies.begin(),taskptr->in_remote_dependencies.end(),std::make_tuple(tgt_i,tgt_j,false));
                              if (it!=taskptr->in_remote_dependencies.end()) {
                                std::get<2>(*it) = true;
                              }
#endif


                            }


                            data->on_fetch.get_future().then(
                                [width,nnz,nblocks](SparseTask2D::data_t * pdata){
                                  //create snodeBlock_t and store it in the extra_data
                                  pdata->extra_data = (void*)new snodeBlock_t(pdata->landing_zone,width,nnz,nblocks);
                                }
                                );

                            data->on_delete.get_future().then(
                                [data](){
                                  delete (snodeBlock_t*)data->extra_data;
                                }
                                );




                            }, this->sp_handle, diagcell._gstorage, diagcell._storage_size, diagcell.nnz(), diagcell.nblocks(), std::get<0>(diagcell._dims) ,ptask->_meta, upcxx::make_view(tgt_cells.begin(),tgt_cells.end())); 
#endif
                      }
                      else {
                        for ( auto & tgt_cell: tgt_cells ) {
                          auto tgt_i = std::get<0>(tgt_cell);
                          auto tgt_j = std::get<1>(tgt_cell);
                          bassert(I==tgt_j);
                          auto taskptr = task_lookup(tgt_i,tgt_j,tgt_j,tgt_j,tgt_i,Factorization::op_type::TRSM);
                          bassert(taskptr!=nullptr); 
                          //mark the dependency as satisfied
#ifdef _DEBUG_DEPENDENCIES_
                              auto it = std::find(taskptr->in_local_dependencies.begin(),taskptr->in_local_dependencies.end(),std::make_tuple(tgt_i,tgt_j,false));
                              if (it!=taskptr->in_local_dependencies.end()) {
                                std::get<2>(*it) = true;
                              }
#endif


                          taskptr->in_prom->fulfill_anonymous(1);
                          //taskptr->in_prom->fulfill_anonymous(1);
//                          taskptr->dep_count--;

                        }
                        //ptask->out_prom->fulfill_anonymous(1);
                      }
                    }

                    //the task pointed by ptask can be deleted when all outgoing communications have been performed.
                    auto fut = ptask->out_prom->finalize();
                    fut.wait();
                    //fut.then( [ptask]() { delete ptask; } );
                    //
                    //
                    //TODO what do do with the task OR return a future
                    ptask->executed = true;
                  };
                }
                break;
              case Factorization::op_type::TRSM:
                {

                  ptask->execute = [this,src_snode,tgt_snode,ptask,I,J,K] () {
                    scope_timer(b,FB_TRSM_TASK);

                    logfileptr->OFS()<<"Exec TRSM"<<" from "<<I<<" to ("<<I<<") cell ("<<K<<","<<I<<")"<<std::endl;

                    auto & od_cell = CELL(K-1,I-1);//*std::static_pointer_cast<snodeBlock_t>(cells_[coord2supidx(tgt_snode-1,src_snode-1)]);
                    //auto & od_cell = cells_[coord2supidx(tgt_snode-1,src_snode-1)];
                    bassert( od_cell.owner == this->iam);

#ifdef SP_THREADS
                    std::thread::id tid = std::this_thread::get_id();
                    auto & tmpBuf = tmpBufs_th[tid];
#else
                    auto & tmpBuf = tmpBufs;
#endif

                    //input data is one or 0 (local diagonal block)
                    bassert(ptask->input_msg.size()<=1);

                    auto ptr_diagCell = pCELL(I-1,I-1).get(); 
                    if ( ptr_diagCell->owner != this->iam ){
                      ptr_diagCell = (snodeBlock_t*)(ptask->input_msg.begin()->get()->extra_data);
                    }
                    bassert(ptr_diagCell!=NULL);


                    //od_cell.trsm(*ptr_diagCell,tmpBuf);

                    //get rid of the shared_ptr
                    ptask->input_msg.clear();

                    //send this to the all facing blocks (tgt_snode-1,*) and off diagonal blocks (*,tgt_snode-1)
                    std::unordered_map<Int,std::list<SparseTask2D::depend_task_t> > data_to_send;
                    for (auto &tpl : ptask->out_dependencies) {
                      auto tgt_i = std::get<0>(tpl);
                      auto tgt_j = std::get<1>(tpl);
                      auto & tgt_cell = CELL(tgt_i-1,tgt_j-1);//*std::static_pointer_cast<snodeBlock_t>(cells_[coord2supidx(tgt_i-1,tgt_j-1)]);
                      data_to_send[tgt_cell.owner].push_back(tpl);
                    }

                    //send factor and update local tasks
//                    ptask->out_prom->require_anonymous(data_to_send.size());
                    for (auto it = data_to_send.begin(); it!=data_to_send.end(); it++) {
                      auto pdest = it->first;
                      auto & tgt_cells = it->second;
                      //serialize data once, and list of meta data
                      //factor is output data so it will not be deleted
                      if ( pdest != this->iam ) {
                        auto cxs = upcxx::source_cx::as_buffered() | upcxx::source_cx::as_promise(*ptask->out_prom);
#if 1
                        upcxx::rpc_ff( pdest, //cxs, 
                            [I] (int sp_handle, upcxx::global_ptr<char> gptr, size_t storage_size, size_t nnz, size_t nblocks, rowind_t width, SparseTask2D::meta_t meta, upcxx::view<SparseTask2D::depend_task_t> target_cells ) { 
//                            gdb_lock();
                            //store pointer & associated metadata somewhere
                            auto data = std::make_shared<SparseTask2D::data_t >();
                            data->in_meta = meta;
                            data->size = storage_size;
                            data->remote_gptr = gptr;
                            //there is a map between sp_handle and task_graphs
                            auto matptr = (symPACKMatrix2D<colptr_t,rowind_t,T> *) g_sp_handle_to_matrix[sp_handle];
                            for ( auto & tgt_cell: target_cells) {
                              auto tgt_i = std::get<0>(tgt_cell);
                              auto tgt_j = std::get<1>(tgt_cell);
                              auto taskptr = matptr->task_lookup(tgt_i,tgt_j,I,tgt_j,tgt_i,Factorization::op_type::UPDATE2D_COMP);
                              taskptr->input_msg.push_back(data);
                              data->target_tasks.push_back(taskptr);
                              taskptr->in_avail_prom->fulfill_anonymous(1);

#ifdef _DEBUG_DEPENDENCIES_
                              auto it = std::find(taskptr->in_remote_dependencies.begin(),taskptr->in_remote_dependencies.end(),std::make_tuple(tgt_i,tgt_j,false));
                              if (it!=taskptr->in_remote_dependencies.end()) {
                                std::get<2>(*it) = true;
                              }
#endif


                            }

                            data->on_fetch.get_future().then(
                                [width,nnz,nblocks](SparseTask2D::data_t * pdata){
                                  //create snodeBlock_t and store it in the extra_data
                                  pdata->extra_data = (void*)new snodeBlock_t(pdata->landing_zone,width,nnz,nblocks);
                                }
                                );

                            data->on_delete.get_future().then(
                                [data](){
                                  delete (snodeBlock_t*)data->extra_data;
                                }
                                );




                            }, this->sp_handle, od_cell._gstorage,od_cell._storage_size, od_cell.nnz(),od_cell.nblocks(), std::get<0>(od_cell._dims),ptask->_meta, upcxx::make_view(tgt_cells.begin(),tgt_cells.end())); 
#endif
                      }
                      else {
                        for ( auto & tgt_cell: tgt_cells ) {
                          auto tgt_i = std::get<0>(tgt_cell);
                          auto tgt_j = std::get<1>(tgt_cell);

                          bassert(K==tgt_i || K==tgt_j); // K is either facing or pivot row

                          auto taskptr = task_lookup(tgt_i,tgt_j,I,tgt_j,tgt_i,Factorization::op_type::UPDATE2D_COMP);
                          bassert(taskptr!=nullptr); 
                          //mark the dependency as satisfied
                          logfileptr->OFS()<<"UPDATE"<<" from "<<I<<" to "<<tgt_j<<" facing cell ("<<K<<","<<I<<") and cell ("<<tgt_j<<","<<I<<") to cell ("<<K<<","<<tgt_j<<")"<<std::endl;
                          //logfileptr->OFS()<<"      input dependencies: ";
                          //for(auto &&tpl : taskptr->in_dependencies){
                          //  logfileptr->OFS()<<"("<<std::get<0>(tpl)<<","<<std::get<1>(tpl)<<") ";
                          //}
                          //
                      
#ifdef _DEBUG_DEPENDENCIES_
                              auto it = std::find(taskptr->in_local_dependencies.begin(),taskptr->in_local_dependencies.end(),std::make_tuple(tgt_i,tgt_j,false));
                              if (it!=taskptr->in_local_dependencies.end()) {
                                std::get<2>(*it) = true;
                              }
#endif


                          taskptr->in_prom->fulfill_anonymous(1);
                          //taskptr->in_prom->fulfill_anonymous(1);
//                          taskptr->dep_count--;

                        }
//                        ptask->out_prom->fulfill_anonymous(1);
                      }
                    }

                    //the task pointed by ptask can be deleted when all outgoing communications have been performed.
                    auto fut = ptask->out_prom->finalize();
                    fut.wait();

                    ptask->executed = true;

                  };
                }
                break;
              case Factorization::op_type::UPDATE2D_COMP:
                {

                  ptask->execute = [this,src_snode,ptask,I,J,K] () {
                    scope_timer(b,FB_UPDATE2D_TASK);

                    logfileptr->OFS()<<"Exec UPDATE"<<" from "<<I<<" to "<<J<<" facing cell ("<<K<<","<<I<<") and cell ("<<J<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;
                    auto & upd_cell = CELL(K-1,J-1);

#ifdef SP_THREADS
                    std::thread::id tid = std::this_thread::get_id();
                    auto & tmpBuf = tmpBufs_th[tid];
#else
                    auto & tmpBuf = tmpBufs;
#endif

                    //TODO do something
                    //input data should be at most two

                    bassert(ptask->input_msg.size()<=2);

                    auto ptr_odCell = pCELL(J-1,I-1).get(); 
                    if ( ptr_odCell->owner != this->iam ){
                      ptr_odCell = (snodeBlock_t*)(ptask->input_msg.begin()->get()->extra_data);
                    }
                    bassert(ptr_odCell!=NULL);

                    auto ptr_facingCell = pCELL(K-1,I-1).get(); 
                    if ( ptr_facingCell->owner != this->iam ){
                      ptr_facingCell = (snodeBlock_t*)(ptask->input_msg.rbegin()->get()->extra_data);
                    }
                    bassert(ptr_facingCell!=NULL);

                    //upd_cell.update(*ptr_odCell,*ptr_facingCell,tmpBuf);

                    //get rid of the shared_ptr
                    ptask->input_msg.clear();

                    //send this to all facing blocks (tgt_snode-1,*) and off diagonal blocks (*,tgt_snode-1)
                    std::unordered_map<Int,std::list<SparseTask2D::depend_task_t> > data_to_send;
                    for (auto &tpl : ptask->out_dependencies) {
                      auto tgt_i = std::get<0>(tpl);
                      auto tgt_j = std::get<1>(tpl);
                      auto & tgt_cell = CELL(tgt_i-1,tgt_j-1);
                      data_to_send[tgt_cell.owner].push_back(tpl);
                    }

                    //send factor and update local tasks
//                    ptask->out_prom->require_anonymous(data_to_send.size());
                    for (auto it = data_to_send.begin(); it!=data_to_send.end(); it++) {
                      auto pdest = it->first;
                      auto & tgt_cells = it->second;
                      //serialize data once, and list of meta data
                      //factor is output data so it will not be deleted
                      if ( pdest != this->iam ) {
                        auto cxs = upcxx::source_cx::as_buffered() | upcxx::source_cx::as_promise(*ptask->out_prom);
#if 1
                        upcxx::rpc_ff( pdest, //cxs, 
                            [ ] (int sp_handle, upcxx::global_ptr<char> gptr, size_t storage_size, size_t nnz, size_t nblocks, rowind_t width, SparseTask2D::meta_t meta, upcxx::view<SparseTask2D::depend_task_t> target_cells ) { 
//                            gdb_lock();
                            //store pointer & associated metadata somewhere
                            auto data = std::make_shared<SparseTask2D::data_t >();
                            //data->target_cells.insert(data->target_cells.begin(),
                            //    target_cells.begin(),target_cells.end());
                            data->in_meta = meta;
                            data->size = storage_size;
                            data->remote_gptr = gptr;

                            //there is a map between sp_handle and task_graphs
                            auto matptr = (symPACKMatrix2D<colptr_t,rowind_t,T> *) g_sp_handle_to_matrix[sp_handle];
                            for ( auto & tgt_cell: target_cells) {
                              auto tgt_i = std::get<0>(tgt_cell);
                              auto tgt_j = std::get<1>(tgt_cell);
                              auto K = tgt_i;
                              auto J = tgt_j;
                              auto taskptr = matptr->task_lookup(tgt_i,tgt_j,J,J,K,
                                (tgt_i==tgt_j?Factorization::op_type::FACTOR:Factorization::op_type::TRSM));
                              taskptr->input_msg.push_back(data);
                              data->target_tasks.push_back(taskptr);
                              taskptr->in_avail_prom->fulfill_anonymous(1);

#ifdef _DEBUG_DEPENDENCIES_
                              auto it = std::find(taskptr->in_remote_dependencies.begin(),taskptr->in_remote_dependencies.end(),std::make_tuple(tgt_i,tgt_j,false));
                              if (it!=taskptr->in_remote_dependencies.end()) {
                                std::get<2>(*it) = true;
                              }
#endif


                            }


                            data->on_fetch.get_future().then(
                                [width,nnz,nblocks](SparseTask2D::data_t * pdata){
                                  //create snodeBlock_t and store it in the extra_data
                                  pdata->extra_data = (void*)new snodeBlock_t(pdata->landing_zone,width,nnz,nblocks);
                                }
                                );

                            data->on_delete.get_future().then(
                                [data](){
                                  delete (snodeBlock_t*)data->extra_data;
                                }
                                );




                            }, this->sp_handle, upd_cell._gstorage, upd_cell._storage_size,upd_cell.nnz(),upd_cell.nblocks(), std::get<0>(upd_cell._dims),ptask->_meta, upcxx::make_view(tgt_cells.begin(),tgt_cells.end())); 
#endif
                      }
                      else {
                        for ( auto & tgt_cell: tgt_cells ) {
                          auto tgt_i = std::get<0>(tgt_cell);
                          auto tgt_j = std::get<1>(tgt_cell);

                          bassert(K==tgt_i && J==tgt_j);
                          auto taskptr = task_lookup(tgt_i,tgt_j,J,J,K,
                              (tgt_i==tgt_j?Factorization::op_type::FACTOR:Factorization::op_type::TRSM));
                          bassert(taskptr!=nullptr); 
                          //mark the dependency as satisfied
                          //taskptr->dep_count--;


#ifdef _DEBUG_DEPENDENCIES_
                              auto it = std::find(taskptr->in_local_dependencies.begin(),taskptr->in_local_dependencies.end(),std::make_tuple(tgt_i,tgt_j,false));
                              if (it!=taskptr->in_local_dependencies.end()) {
                                std::get<2>(*it) = true;
                              }
#endif




                          taskptr->in_prom->fulfill_anonymous(1);
                          //taskptr->in_prom->fulfill_anonymous(1);

                        }
//                        ptask->out_prom->fulfill_anonymous(1);
                      }
                    }

                    //the task pointed by ptask can be deleted when all outgoing communications have been performed.
                    auto fut = ptask->out_prom->finalize();
                    fut.wait();

                    ptask->executed = true;
                  };




                }
                break;
              default:
                delete ptask;
                break;
            }
          }
        }


      }



      //redistribute the matrix

    } 

  template <typename colptr_t, typename rowind_t, typename T>
    void symPACKMatrix2D<colptr_t,rowind_t,T>::DistributeMatrix(DistSparseMatrix<T> & pMat ){
    } 

  template <typename colptr_t, typename rowind_t, typename T>
    void symPACKMatrix2D<colptr_t,rowind_t,T>::Factorize( ){
      //bool finished = false;
      //while (!finished) {
      //  finished = true;
      //  for (auto it = this->task_graph.begin(); it != this->task_graph.end(); it++) {
      //    auto idx = it->first;
      //    auto tasklist = it->second;

      //    for ( auto && ptask: tasklist) {
      //      if (ptask->dep_count==0) {
      //        if (!ptask->executed) {
      //          ptask->execute();
      //        }
      //      }
      //      else {
      //        finished = false;
      //      }
      //    }
      //  }
      //}
#ifdef _DEBUG_DEPENDENCIES_
      this->scheduler.execute(this->task_graph,this->SupMembership_.data());
#else
      this->scheduler.execute(this->task_graph);
#endif
    } 


namespace scheduling {

  template <typename ttask_t , typename ttaskgraph_t >
#ifdef _DEBUG_DEPENDENCIES_
  inline void Scheduler2D<ttask_t,ttaskgraph_t>::execute(ttaskgraph_t & task_graph, Int * SupMembership )
#else
  inline void Scheduler2D<ttask_t,ttaskgraph_t>::execute(ttaskgraph_t & task_graph )
#endif
  {
//    if(upcxx::rank_me()==0){gdb_lock();}
    //fulfill_anonymous each in_prom to launch the computations
    int64_t local_task_cnt = 0;
    for (auto it = task_graph.begin(); it != task_graph.end(); it++) {
      auto idx = it->first;
      auto tasklist = it->second;
      for ( auto && ptask: tasklist) {
        local_task_cnt++;
        //this finalize the in_prom promise, trigering tasks without input dependencies 
        ptask->in_prom->fulfill_anonymous(1);
      }
    }

    while(local_task_cnt>0){
      if (!ready_tasks.empty()) {
        auto ptask = ready_tasks.front();
        ready_tasks.pop_front();
//        bassert(ptask->dep_count==0 && ptask->executed==false);
        ptask->execute(); 
        local_task_cnt--;
      }

      //handle communications
      if (!avail_tasks.empty()) {
        //get all comms from a task
        auto ptask = avail_tasks.front();
        avail_tasks.pop_front();
        for(auto & msg : ptask->input_msg){
          msg->allocate();


          msg->fetch().then([ptask](incoming_data_t<ttask_t,meta_t> * pmsg){
              //fulfill promise by one, when this reaches 0, ptask is moved to scheduler.ready_tasks
              ptask->in_prom->fulfill_anonymous(1);

              //for(auto ptask: pmsg->target_tasks){
                //fulfill promise by one, when this reaches 0, ptask is moved to scheduler.ready_tasks
              //  ptask->in_prom->fulfill_anonymous(1);
              //}

              });
        } 
        //fulfill by 1 to finalize the promise 
        //ptask->in_prom->fulfill_anonymous(1);
      }
      
      
      upcxx::progress();

#ifdef _DEBUG_DEPENDENCIES_
        for(auto it = task_graph.begin(); it != task_graph.end(); it++){
          auto idx = it->first;
          auto tasklist = it->second;

          for( auto && ptask: tasklist){
            if (!ptask->executed){

              auto meta = &ptask->_meta;
              auto & src_snode = std::get<0>(meta[0]);
              auto & tgt_snode = std::get<1>(meta[0]);
              auto & lt_first_row = std::get<3>(meta[0]);
              auto & facing_row = std::get<4>(meta[0]);
              auto & type = std::get<2>(meta[0]);

              auto I = src_snode;
              auto J = tgt_snode;
              auto K = SupMembership[facing_row-1];

              switch(type){
                case Factorization::op_type::TRSM:
                  {
                    logfileptr->OFS()<<"TRSM"<<" from "<<I<<" to ("<<I<<") cell ("<<K<<","<<I<<")"<<std::endl;
                  }
                  break;
                case Factorization::op_type::FACTOR:
                  {
                    logfileptr->OFS()<<"FACTOR"<<" cell ("<<I<<","<<I<<")"<<std::endl;
                  }
                  break;
                case Factorization::op_type::UPDATE2D_COMP:
                  {
                    logfileptr->OFS()<<"UPDATE"<<" from "<<I<<" to "<<J<<" facing cell ("<<K<<","<<I<<") and cell ("<<J<<","<<I<<") to cell ("<<K<<","<<J<<")"<<std::endl;
                  }
                  break;
                default:
                  delete ptask;
                  break;
              }




              logfileptr->OFS()<<"      local input dependencies: ";
              for(auto &&tpl : ptask->in_local_dependencies){
                logfileptr->OFS()<<"("<<std::get<0>(tpl)<<","<<std::get<1>(tpl)<<","<<(std::get<2>(tpl)?"yes":"no")<<") ";
              }
              logfileptr->OFS()<<std::endl;
              logfileptr->OFS()<<"      remote input dependencies: ";
              for(auto &&tpl : ptask->in_remote_dependencies){
                logfileptr->OFS()<<"("<<std::get<0>(tpl)<<","<<std::get<1>(tpl)<<","<<(std::get<2>(tpl)?"yes":"no")<<") ";
              }
            }
          }

        }
#endif

    }

//    //fulfill_anonymous each in_prom to launch the computations
//    int64_t local_task_cnt = 0;
//    for (auto it = task_graph.begin(); it != task_graph.end(); it++) {
//      auto idx = it->first;
//      auto tasklist = it->second;
//      for ( auto && ptask: tasklist) {
//        local_task_cnt++;
//        ptask->in_prom->fulfill_anonymous(1);
//      }
//    }
//   
//    std::deque<incoming_data_t> & comm_queue = g_sp_handle_incoming[this->sp_handle];
//
//    while(local_task_cnt>0){
//      if (!ready_tasks.empty()) {
//        auto ptask = ready_tasks.front();
//        ready_tasks.pop_front();
//        bassert(ptask->dep_count==0 && ptask->executed==false);
//        ptask->execute(); 
//      }
//
//      //handle communications
//      if (!avail_tasks.empty()) {
//        //get all comms from a task
//        auto ptask = avail_tasks.front();
//        avail_tasks.pop_front();
//        
//      }
//    }
  }


//  inline void taskGraph2D_t::execute( Scheduler2D & scheduler ){
//    //fulfill_anonymous each in_prom to launch the computations
//    int64_t local_task_cnt = 0;
//    for (auto it = begin(); it != end(); it++) {
//      auto idx = it->first;
//      auto tasklist = it->second;
//      for ( auto && ptask: tasklist) {
//        local_task_cnt++;
//        ptask->in_prom->fulfill_anonymous(1);
//      }
//    }
//
//    while(local_task_cnt>0){
//      if (!scheduler.ready_tasks.empty()) {
//        auto ptask = scheduler.ready_tasks.front();
//        scheduler.ready_tasks.pop_front();
//        bassert(ptask->dep_count==0 && ptask->executed==false);
//        ptask->execute(); 
//      }
//
//      //handle communications
//      if (!scheduler.avail_tasks.empty()) {
//        //get all comms from a task
//        auto ptask = scheduler.avail_tasks.front();
//        scheduler.avail_tasks.pop_front();
//        for(auto & msg : ptask->input_msg){
//          msg->allocate();
//          msg->fetch().then([](incoming_data_t * pmsg){
//              if (!pmsg->cell) {
//                //call constructor of snodeBlock_t
//              }
//              for(auto ptask: pmsg->target_tasks){
//                //fulfill promise by one, when this reaches 0, ptask is moved to scheduler.read_tasks
//                ptask->in_prom->fulfill_anonymous(1);
//              }
//
//              });
//        } 
//      }
//    }
//  }
}


#endif

}

#endif //_SYMPACK_MATRIX2D_DECL_HPP_

