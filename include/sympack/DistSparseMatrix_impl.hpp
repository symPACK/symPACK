/// @file sparse_matrix_impl.hpp
/// @brief Sparse matrix and Distributed sparse matrix in compressed
/// column format.
/// @author Lin Lin
/// @date 2012-11-10
#ifndef _DIST_SPARSE_MATRIX_IMPL_HPP_
#define _DIST_SPARSE_MATRIX_IMPL_HPP_

#include "sympack/Environment.hpp"
#include "sympack/DistSparseMatrix.hpp"
#include "sympack/utility.hpp"
#include "sympack/SparseMatrixStructure.hpp"

#include <vector>

#ifdef UPCXX
#include <upcxx.h>
#endif
namespace SYMPACK{






  template <class F> SparseMatrixStructure  DistSparseMatrix<F>::GetLocalStructure() const {
    return Local_;
  }

  template <class F> SparseMatrixStructure DistSparseMatrix<F>::GetGlobalStructure(){
    if(!globalAllocated){
      Local_.ToGlobal(Global_,comm);
      globalAllocated = true;
    }
    return Global_;
  }


  template <class F> template <typename T> void DistSparseMatrix<F>::ConvertData(const int n, const int nnz, const int * colptr, const int * rowidx, const T * nzval ,bool onebased){
    int np;
    int iam;

    MPI_Comm_size(comm,&np);
    MPI_Comm_rank(comm, &iam);

    this->size = n; 
    this->nnz = nnz; 

    //fill global structure info as we have it directly
    this->Global_.size = this->size;
    this->Global_.nnz = this->nnz;
    this->Global_.colptr.resize(n+1);
    std::copy(colptr,colptr+n+1,&this->Global_.colptr[0]);
    this->Global_.rowind.resize(nnz+1);
    std::copy(rowidx,rowidx+nnz+1,&this->Global_.rowind[0]);

    if(!onebased){
      //move to 1 based indices
      for(int i=0;i<this->Global_.colptr.size();i++){ ++this->Global_.colptr[i]; }
      for(int i=0;i<this->Global_.rowind.size();i++){ ++this->Global_.rowind[i]; }
    }

    this->globalAllocated = true;

    //Compute local structure info
    // Compute the number of columns on each processor
    Idx numColLocal, numColFirst;
    numColFirst = this->size / np;
    SYMPACK::vector<Idx> numColLocalVec(np,numColFirst );
    numColLocalVec[np-1] = this->size - numColFirst * (np-1) ;  // Modify the last entry	
    numColLocal = numColLocalVec[iam];

    this->Local_.colptr.resize( numColLocal + 1 );

    for( Idx i = 0; i < numColLocal+1; i++ ){
      this->Local_.colptr[i] = this->Global_.colptr[iam * numColFirst+i] - this->Global_.colptr[iam * numColFirst] + 1;
    }

    this->Local_.size = this->size;
    // Calculate nnz_loc on each processor
    this->Local_.nnz = this->Local_.colptr[numColLocal] - this->Local_.colptr[0];
    this->Local_.colptr.back() = this->Local_.nnz +1;

    // resize rowind and nzval appropriately 
    this->Local_.rowind.resize( this->Local_.nnz );
    this->nzvalLocal.resize ( this->Local_.nnz );

    //Read my row indices
    Int prevRead = 0;
    Int numRead = 0;
    for( Int ip = 0; ip <iam; ip++ ){	
      prevRead += this->Global_.colptr[ip*numColFirst + numColLocalVec[ip]]
        - this->Global_.colptr[ip*numColFirst];
    }

    numRead = this->Global_.colptr[iam*numColFirst + numColLocalVec[iam]] - this->Global_.colptr[iam*numColFirst];
    std::copy(&this->Global_.rowind[prevRead],&this->Global_.rowind[prevRead+numRead],&this->Local_.rowind[0]);

    for(int i = 0; i<numRead;++i){
      this->nzvalLocal[i] = (F)((const T*)nzval)[prevRead+i];
    }
    //copy appropriate nnz values
    //    if(cscptr->value_type == REAL){
    //      std::copy(&((const double*)cscptr->values)[prevRead],&((const double*)cscptr->values)[prevRead+numRead],this->nzvalLocal.Data());
    //    }
    //    else if(cscptr->value_type == COMPLEX){
    //      std::copy(&((const double*)cscptr->values)[prevRead],&((const double*)cscptr->values)[prevRead+numRead],this->nzvalLocal.Data());
    //    }

#if 1
{
vector<Ptr> dummy;
this->Global_.colptr.swap(dummy);
}
{
vector<Idx> dummy;
this->Global_.rowind.swap(dummy);
}
this->Global_.nnz=-1;
this->Global_.size=1;
    this->globalAllocated = false;

//Local_.ToGlobal(Global_,comm);

#endif



  }



  template <class F> void DistSparseMatrix<F>::CopyData(const int n, const int nnz, const int * colptr, const int * rowidx, const F * nzval ,bool onebased){
    int np;
    int iam;

    MPI_Comm_size(comm,&np);
    MPI_Comm_rank(comm, &iam);

    //fill global structure info as we have it directly
    this->size = n; 
    this->nnz = nnz; 
    this->Global_.size = this->size;
    this->Global_.nnz = this->nnz;
    this->Global_.colptr.resize(n+1);
    std::copy(colptr,colptr+n+1,&this->Global_.colptr[0]);
    this->Global_.rowind.resize(nnz+1);
    std::copy(rowidx,rowidx+nnz+1,&this->Global_.rowind[0]);

    if(!onebased){
      //move to 1 based indices
      for(int i=0;i<this->Global_.colptr.size();i++){ ++this->Global_.colptr[i]; }
      for(int i=0;i<this->Global_.rowind.size();i++){ ++this->Global_.rowind[i]; }
    }

    this->globalAllocated = true;

    //Compute local structure info
    // Compute the number of columns on each processor
    SYMPACK::vector<Int> numColLocalVec(np);
    Int numColLocal, numColFirst;
    numColFirst = this->size / np;
    std::fill(numColLocalVec.begin(),numColLocalVec.end(),numColFirst);

    numColLocalVec[np-1] = this->size - numColFirst * (np-1) ;  // Modify the last entry	
    numColLocal = numColLocalVec[iam];

    this->Local_.colptr.resize( numColLocal + 1 );

    for( Int i = 0; i < numColLocal+1; i++ ){
      this->Local_.colptr[i] = this->Global_.colptr[iam * numColFirst+i] - this->Global_.colptr[iam * numColFirst] + 1;
    }

    this->Local_.size = this->size;
    // Calculate nnz_loc on each processor
    this->Local_.nnz = this->Local_.colptr[numColLocal] - this->Local_.colptr[0];
    this->Local_.colptr.back() = this->Local_.nnz +1;

    // resize rowind and nzval appropriately 
    this->Local_.rowind.resize( this->Local_.nnz );
    this->nzvalLocal.resize ( this->Local_.nnz );

    //Read my row indices
    Int prevRead = 0;
    Int numRead = 0;
    for( Int ip = 0; ip <iam; ip++ ){	
      prevRead += this->Global_.colptr[ip*numColFirst + numColLocalVec[ip]]
        - this->Global_.colptr[ip*numColFirst];
    }

    numRead = this->Global_.colptr[iam*numColFirst + numColLocalVec[iam]] - this->Global_.colptr[iam*numColFirst];
    std::copy(&this->Global_.rowind[prevRead],&this->Global_.rowind[prevRead+numRead],&this->Local_.rowind[0]);

    std::copy(&((const F*)nzval)[prevRead],&((const F*)nzval)[prevRead+numRead],&this->nzvalLocal[0]);
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


  template <class F> void DistSparseMatrix<F>::Dump() const{
    Int baseval = 1;
    Int numColFirst = size / np;
    Int firstLocCol = iam * numColFirst+1;

    for(Idx locCol = 0 ; locCol< Local_.colptr.size()-1; locCol++){
      Idx col = locCol + firstLocCol;
      Ptr colbeg = Local_.colptr[locCol]-baseval; //now 0 based
      Ptr colend = Local_.colptr[locCol+1]-baseval; // now 0 based 

      logfileptr->OFS()<<col<<": ";
      for(Ptr pos = colbeg; pos<colend; pos++){
        Idx row = Local_.rowind[pos];
        logfileptr->OFS()<<row<<" ";
      }
      logfileptr->OFS()<<endl;
    }
  }

//chunk of code that might be useful to expand a matrix to unsymmetric storage
#if 0
    TIMER_START(EXPANDING_MATRIX);
    DistSparseMatrix<T> ExpA(CommEnv_->MPI_GetComm());
    {
      double timeSta = get_time();
      ExpA.size = pMat.size;
      ExpA.nnz = 2*pMat.nnz-pMat.size;
      Idx numColFirst = std::max(1,ExpA.size / np);
      SYMPACK::vector<Idx> numColLocalVec(np,numColFirst);
      numColLocalVec[np-1] = ExpA.size - numColFirst * (np-1);  // Modify the last entry	
      Idx numColLocal = numColLocalVec[iam];
      //Expand A to asymmetric storage
      Idx localFirstCol = iam*numColFirst+1;
      Idx localLastCol = localFirstCol+numColLocal-1;

      //Initialize the Local structure
      ExpA.Local_.size = ExpA.size;
      ExpA.Local_.colptr.resize(numColLocal+1);

      for( Idx i = 0; i < numColLocal + 1; i++ ){
        ExpA.Local_.colptr[i] = Global_->expColptr[iam * numColFirst+i] - Global_->expColptr[iam * numColFirst] + 1;
      }

      ExpA.Local_.nnz = ExpA.Local_.colptr[numColLocal] - ExpA.Local_.colptr[0];

      ExpA.Local_.rowind.resize(ExpA.Local_.nnz);

      Ptr globalColBeg = Global_->expColptr[iam*numColFirst];
      std::copy(&Global_->expRowind[globalColBeg-1],&Global_->expRowind[globalColBeg-1]+ExpA.Local_.nnz,&ExpA.Local_.rowind[0]);
      ExpA.nzvalLocal.resize(ExpA.Local_.nnz);

#ifdef _DEBUG_
      logfileptr->OFS()<<"Global_.colptr: "<<Global_->colptr<<endl;
      logfileptr->OFS()<<"Global_.rowind: "<<Global_->rowind<<endl;
      logfileptr->OFS()<<"Global_.expColptr: "<<Global_->expColptr<<endl;
      logfileptr->OFS()<<"Global_.expRowind: "<<Global_->expRowind<<endl;
      logfileptr->OFS()<<"ExpA.colptr: "<<ExpA.Local_.colptr<<endl;
      logfileptr->OFS()<<"ExpA.rowind: "<<ExpA.Local_.rowind<<endl;
      logfileptr->OFS()<<"pMat.colptr: "<<pMat.Local_.colptr<<endl;
      logfileptr->OFS()<<"pMat.rowind: "<<pMat.Local_.rowind<<endl;
#endif

      SYMPACK::vector<Ptr> localDestColHead = ExpA.Local_.colptr;

      SYMPACK::vector<T> recvNzval;
      SYMPACK::vector<Ptr> recvColptr;
      SYMPACK::vector<Idx> recvRowind;
      for(Int proc = 0; proc<np; ++proc){
        //communicate colptr
        //Broadcast size
        Int size = iam==proc?pMat.Local_.colptr.size():0;
        MPI_Bcast((void*)&size,sizeof(size),MPI_BYTE,proc,ExpA.comm);
        //Broadcast colptr
        if(iam==proc){
          MPI_Bcast((void*)&pMat.Local_.colptr[0],size*sizeof(Ptr),MPI_BYTE,proc,ExpA.comm);
        }
        else{
          recvColptr.resize(size);
          MPI_Bcast((void*)&recvColptr[0],size*sizeof(Ptr),MPI_BYTE,proc,ExpA.comm);
        }

        //communicate rowind
        size = iam==proc?pMat.Local_.rowind.size():0;
        MPI_Bcast(&size,sizeof(size),MPI_BYTE,proc,ExpA.comm);
        //Broadcast rowind
        if(iam==proc){
          MPI_Bcast((void*)&pMat.Local_.rowind[0],size*sizeof(Idx),MPI_BYTE,proc,ExpA.comm);
        }
        else{
          recvRowind.resize(size);
          MPI_Bcast((void*)&recvRowind[0],size*sizeof(Idx),MPI_BYTE,proc,ExpA.comm);
        }


        //communicate nzvalLocal
        //Broadcast size
        size = iam==proc?pMat.nzvalLocal.size():0;
        MPI_Bcast((void*)&size,sizeof(size),MPI_BYTE,proc,ExpA.comm);
        //Broadcast nzvalLocal
        if(iam==proc){
          MPI_Bcast((void*)&pMat.nzvalLocal[0],size*sizeof(T),MPI_BYTE,proc,ExpA.comm);
        }
        else{
          recvNzval.resize(size);
          MPI_Bcast((void*)&recvNzval[0],size*sizeof(T),MPI_BYTE,proc,ExpA.comm);
        }

        int colptrSize = iam==proc?pMat.Local_.colptr.size()-1:recvColptr.size()-1;
        const Ptr * colptr = iam==proc?&pMat.Local_.colptr[0]:&recvColptr[0];
        const Idx * rowind = iam==proc?&pMat.Local_.rowind[0]:&recvRowind[0];
        const T * nzval = iam==proc?&pMat.nzvalLocal[0]:&recvNzval[0];
        //Parse the received data and copy it in the ExpA.nzvalLocal
        Int recvFirstCol = proc*numColFirst+1;
        Int recvLastCol = recvFirstCol+numColLocalVec[proc]-1;
        for(Idx colIdx = 1; colIdx<=colptrSize;++colIdx){
          Idx col = recvFirstCol + colIdx-1;
          Ptr colBeg = colptr[colIdx-1];
          Ptr colEnd = colptr[colIdx]-1;
          for(Ptr rowIdx = colBeg; rowIdx<=colEnd;++rowIdx){
            Idx row = rowind[rowIdx-1];
            if(row>=localFirstCol && row<=localLastCol){
              //copy only upper triangular part
              if(col<row){
                Idx localCol = row - localFirstCol +1;
                Ptr & localDestColIdx = localDestColHead[localCol-1];
                ExpA.nzvalLocal[localDestColIdx-1] = nzval[rowIdx-1];
                localDestColIdx++;
              }
            }
          }
        }
      }
      //copy upper triangular part
      for(Idx colIdx = 1; colIdx<pMat.Local_.colptr.size();++colIdx){
        Ptr & localDestColIdx = localDestColHead[colIdx-1];
        Ptr colBeg = pMat.Local_.colptr[colIdx-1];
        Ptr colEnd = pMat.Local_.colptr[colIdx]-1;
        for(Ptr rowIdx = colBeg; rowIdx<=colEnd;++rowIdx){
          Idx row = pMat.Local_.rowind[rowIdx-1];
          ExpA.nzvalLocal[localDestColIdx-1] = pMat.nzvalLocal[rowIdx-1];
          localDestColIdx++;
        }
      }
      double timeEnd = get_time();


//#ifdef _DEBUG_
//      logfileptr->OFS()<<"ExpA.nzvalLocal: "<<ExpA.nzvalLocal<<endl;
//      logfileptr->OFS()<<"pMat.nzvalLocal: "<<pMat.nzvalLocal<<endl;
//#endif

      //now we can do the copy by doing sends of whole columns to the dest processor
      if(iam==0){
        cout<<"Time for expanding matrix to asymmetric: "<<timeEnd-timeSta<<endl;
      }
    }

    TIMER_STOP(EXPANDING_MATRIX);
#endif





}


#endif // _DIST_SPARSE_MATRIX_IMPL_HPP_
