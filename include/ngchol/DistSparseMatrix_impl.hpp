/// @file sparse_matrix_impl.hpp
/// @brief Sparse matrix and Distributed sparse matrix in compressed
/// column format.
/// @author Lin Lin
/// @date 2012-11-10
#ifndef _DIST_SPARSE_MATRIX_IMPL_HPP_
#define _DIST_SPARSE_MATRIX_IMPL_HPP_

#include "ngchol/Environment.hpp"
#include "ngchol/DistSparseMatrix.hpp"
#include "ngchol/utility.hpp"
#include "ngchol/SparseMatrixStructure.hpp"

#include <vector>

#ifdef UPCXX
#include <upcxx.h>
#endif
namespace LIBCHOLESKY{






  template <class F> SparseMatrixStructure  DistSparseMatrix<F>::GetLocalStructure() const {
    return Local_;
  }

  template <class F> SparseMatrixStructure DistSparseMatrix<F>::GetGlobalStructure(){
    if(!globalAllocated){
      Local_.ToGlobal(Global_);
      globalAllocated = true;
    }
    return Global_;
  }



  template <class F> void DistSparseMatrix<F>::CopyData(const int n, const int nnz, const int * colptr, const int * rowidx, const F * nzval ){
    int np;
    int iam;

      MPI_Comm_size(comm,&np);
      MPI_Comm_rank(comm, &iam);

    //fill global structure info as we have it directly
    this->size = n; 
    this->nnz = nnz; 
    this->Global_.size = this->size;
    this->Global_.nnz = this->nnz;
    this->Global_.colptr.Resize(n+1);
    std::copy(colptr,colptr+n+1,this->Global_.colptr.Data());
    this->Global_.rowind.Resize(nnz+1);
    std::copy(rowidx,rowidx+nnz+1,this->Global_.rowind.Data());

    //move to 1 based indices
    for(int i=0;i<this->Global_.colptr.m();i++){ ++this->Global_.colptr[i]; }
    for(int i=0;i<this->Global_.rowind.m();i++){ ++this->Global_.rowind[i]; }

    this->globalAllocated = true;

    //Compute local structure info
    // Compute the number of columns on each processor
    IntNumVec numColLocalVec(np);
    Int numColLocal, numColFirst;
    numColFirst = this->size / np;
    SetValue( numColLocalVec, numColFirst );
    numColLocalVec[np-1] = this->size - numColFirst * (np-1) ;  // Modify the last entry	
    numColLocal = numColLocalVec[iam];

    this->Local_.colptr.Resize( numColLocal + 1 );

    for( Int i = 0; i < numColLocal+1; i++ ){
      this->Local_.colptr[i] = this->Global_.colptr[iam * numColFirst+i] - this->Global_.colptr[iam * numColFirst] + 1;
    }

    this->Local_.size = this->size;
    // Calculate nnz_loc on each processor
    this->Local_.nnz = this->Local_.colptr[numColLocal] - this->Local_.colptr[0];

    // Resize rowind and nzval appropriately 
    this->Local_.rowind.Resize( this->Local_.nnz );
    this->nzvalLocal.Resize ( this->Local_.nnz );

    //Read my row indices
    Int prevRead = 0;
    Int numRead = 0;
    for( Int ip = 0; ip <iam; ip++ ){	
      prevRead += this->Global_.colptr[ip*numColFirst + numColLocalVec[ip]]
        - this->Global_.colptr[ip*numColFirst];
    }

    numRead = this->Global_.colptr[iam*numColFirst + numColLocalVec[iam]] - this->Global_.colptr[iam*numColFirst];
    std::copy(&this->Global_.rowind[prevRead],&this->Global_.rowind[prevRead+numRead],this->Local_.rowind.Data());

      std::copy(&((const F*)nzval)[prevRead],&((const F*)nzval)[prevRead+numRead],this->nzvalLocal.Data());
    //copy appropriate nnz values
//    if(cscptr->value_type == REAL){
//      std::copy(&((const double*)cscptr->values)[prevRead],&((const double*)cscptr->values)[prevRead+numRead],this->nzvalLocal.Data());
//    }
//    else if(cscptr->value_type == COMPLEX){
//      std::copy(&((const double*)cscptr->values)[prevRead],&((const double*)cscptr->values)[prevRead+numRead],this->nzvalLocal.Data());
//    }
  }


  template <class F> DistSparseMatrix<F>::DistSparseMatrix(const int n, const int nnz, const int * colptr, const int * rowidx, const F * nzval ,MPI_Comm oComm ):comm(oComm){
    globalAllocated=false;
    this->CopyData(n,nnz,colptr,rowidx,nzval);
  }



}


#endif // _DIST_SPARSE_MATRIX_IMPL_HPP_
