#include "sympack/SparseMatrixStructure.hpp"
#include "sympack/DistSparseMatrixGraph.hpp"
#include "sympack/ETree.hpp"
#include "sympack/utility.hpp"
#include <limits>       // std::numeric_limits

#include <iterator>
#include <set>
#include <list>
#include <vector>
#include <algorithm>

#define USE_REDUCE

namespace SYMPACK{


  //special reduction: first element contains the max local sum
  void _PtrSum( void *in, void *inout, int *len, MPI_Datatype *dptr ) 
  { 
    int i; 

    Ptr * pinout = (Ptr*)inout;
    Ptr * pin = (Ptr*)in;
    pinout[0] = max(pinout[0],pin[0]);
    //logfileptr->OFS()<<"REDUCE max is "<<pinout[0]<<" vs (max) "<<pin[0]<<endl;
#pragma unroll
    for (i=1; i< *len; ++i) { 
      //logfileptr->OFS()<<"REDUCE elem is "<<pin[i]<<endl;
      pinout[i] += pin[i];
    } 
  } 


  SparseMatrixStructure::SparseMatrixStructure(const DistSparseMatrixGraph & G){
    bIsExpanded = G.IsExpanded();
    baseval = G.baseval;
    keepDiag = 1;
    sorted = 0;
    //sorted = G.sorted;
    size = G.size;

    int np;
    int iam;

    int ismpi=0;
    MPI_Initialized( &ismpi);

    //FIXME needs to be passed as an argument ?
    //MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm comm = G.comm;
    int isnull= (comm == MPI_COMM_NULL);
    //    logfileptr->OFS()<<ismpi<<std::endl;
    //    logfileptr->OFS()<<comm<<std::endl;
    //    logfileptr->OFS()<<MPI_COMM_NULL<<std::endl;
    //    logfileptr->OFS()<<isnull<<std::endl;

    if(ismpi && isnull==0){
      MPI_Comm_size(comm,&np);
      MPI_Comm_rank(comm, &iam);
    }
    else{
      //throw an exception
      throw std::logic_error("MPI needs to be available.");
    }

    Ptr * pcolptr = NULL;
    Idx * prowind = NULL;

    if(bIsExpanded){
      expColptr.resize(size+1);
      pcolptr = &expColptr[0];
    }
    else{
      colptr.resize(size+1);
      pcolptr = &colptr[0];
    }


    MPI_Datatype MPI_SYMPACK_PTR; 
    MPI_Datatype MPI_SYMPACK_IDX; 
    MPI_Type_contiguous( sizeof(Ptr), MPI_BYTE, &MPI_SYMPACK_PTR );
    MPI_Type_contiguous( sizeof(Idx), MPI_BYTE, &MPI_SYMPACK_IDX );
    MPI_Type_commit( &MPI_SYMPACK_PTR ); 
    MPI_Type_commit( &MPI_SYMPACK_IDX ); 

    /* Allgatherv for row indices. */ 
    SYMPACK::vector<int> prevnz(np);
    SYMPACK::vector<int> rcounts(np);
    Ptr lnnz = G.LocalEdgeCount();
    MPI_Allgather(&lnnz, 1, MPI_SYMPACK_PTR, &rcounts[0], 1, MPI_SYMPACK_PTR, comm);

    prevnz[0] = 0; 
    for (Int i = 0; i < np-1; ++i) { prevnz[i+1] = prevnz[i] + rcounts[i]; } 

    if(G.keepDiag){
      this->nnz = 0;
    }
    else{
      this->nnz = this->size;
    }
    for (Int i = 0; i < np; ++i) { this->nnz += rcounts[i]; } 

    if(bIsExpanded){
      expRowind.resize(this->nnz);
      prowind = &expRowind[0];
    }
    else{
      rowind.resize(this->nnz);
      prowind = &rowind[0];
    }

    SYMPACK::vector<int> rdispls = prevnz;
    MPI_Allgatherv(&G.rowind[0], lnnz, MPI_SYMPACK_IDX, &prowind[0],&rcounts[0], &rdispls[0], MPI_BYTE, comm); 
    
    /* Allgatherv for colptr */
    // Compute the number of columns on each processor
    Int numColFirst = std::max(1,size / np);
    std::fill(rcounts.begin(),rcounts.end(),numColFirst);
    rcounts[np-1] = (size - numColFirst * (np-1));  // Modify the last entry     

    rdispls[0] = 0;
    for (Int i = 0; i < np-1; ++i) { rdispls[i+1] = rdispls[i] + rcounts[i]; } 

    Idx locN = G.LocalVertexCount();
    MPI_Allgatherv(&G.colptr[0], locN, MPI_SYMPACK_PTR, &pcolptr[0], &rcounts[0], &rdispls[0], MPI_SYMPACK_PTR, comm);

    /* Recompute column pointers. */
    for (Int p = 1; p < np; p++) {
      Int idx = rdispls[p];
      for (Int j = 0; j < rcounts[p]; ++j) pcolptr[idx++] += prevnz[p];
    }
    pcolptr[this->size]= this->nnz+baseval;

    //add diagonal entries if necessary
    if(!G.keepDiag){
      for(Idx col = this->size-1; col>= 0; col--){
        Ptr colbeg = pcolptr[col]-baseval;//0 based
        Ptr colend = pcolptr[col+1]-baseval;//0 based
        //shift this by col+1;
        assert(colend+col+1<=this->nnz);
        std::copy_backward(&prowind[0] + colbeg,&prowind[0]+colend,&prowind[0] + colend + col + 1);
        //add diagonal entry
        prowind[colbeg+col] = col + baseval;
      }
      //recompute column pointers
      for(Idx col = 0; col < this->size;col++){ pcolptr[col] += col; }
    }

    bIsGlobal = true;

    MPI_Type_free(&MPI_SYMPACK_PTR);
    MPI_Type_free(&MPI_SYMPACK_IDX);
  }




  SparseMatrixStructure::SparseMatrixStructure(){
    bIsGlobal = false;
    bIsExpanded = false;
    baseval = 1;
    keepDiag = 1;
    sorted = 1;
    nnz = 0;
    size = 0;
  }

  void SparseMatrixStructure::ClearExpandedSymmetric(){
    {
      SYMPACK::vector<Ptr> dummy;
      expColptr.swap(dummy);
    }
    {
      SYMPACK::vector<Idx> dummy;
      expRowind.swap(dummy);
    }
    bIsExpanded=false;
  }





  void SparseMatrixStructure::ExpandSymmetric(MPI_Comm * workcomm){
    TIMER_START(EXPAND);
    if(!bIsExpanded){
      if(bIsGlobal){

      //code from sparsematrixconverter
      /* set-up */

      SYMPACK::vector<Ptr> cur_col_nnz(size);
      SYMPACK::vector<Ptr> new_col_nnz(size);
//gdb_lock(0);
      /*
       * Scan A and count how many new non-zeros we'll need to create.
       *
       * Post:
       *   cur_col_nnz[i] == # of non-zeros in col i of the original symmetric 
       *                     matrix.
       *   new_col_nnz[i] == # of non-zeros to be stored in col i of the final 
       *                     expanded matrix.
       *   new_nnz == total # of non-zeros to be stored in the final expanded 
       *              matrix.
       */


//logfileptr->OFS()<<colptr<<endl;
//logfileptr->OFS()<<rowind<<endl;


//gdb_lock(0);
      Int new_nnz = 0;
      for (Int i = 0; i < size; i++) 
      {    
        cur_col_nnz[i] = colptr[i+1] - colptr[i];
        new_col_nnz[i] = cur_col_nnz[i];
        new_nnz += new_col_nnz[i];
      }    

      for (Int i = 0; i < size; i++) 
      {    
        Int k;
        for (k = colptr[i]; k < colptr[i+1]; k++) 
        {
          Int j = rowind[k-1]-1;
          if (j != i)
          {
            new_col_nnz[j]++;
            new_nnz++;
          }
        }
      }

      /*
       *  Initialize row pointers in expanded matrix.
       *
       *  Post:
       *    new_colptr initialized to the correct, final values.
       *    new_col_nnz[i] reset to be equal to cur_col_nnz[i].
       */

      expColptr.resize(size+1);
      expColptr[0] = 1;
      for (Int i = 1; i <= size; i++)
      {
        expColptr[i] = expColptr[i-1] + new_col_nnz[i-1];
        new_col_nnz[i-1] = cur_col_nnz[i-1];
      }
      
      for (Int i = 0; i < size; i++){
        new_col_nnz[i] = 0;
      }
      //expColptr[size] = new_nnz;

      expRowind.resize(new_nnz,-1);

      /*
       *  Complete expansion of A to full storage.
       *
       *  Post:
       *    (new_colptr, new_rowind, new_values) is the full-storage equivalent of A.
       *    new_col_nnz[i] == # of non-zeros in col i of A.
       */



      for (Int i = 0; i < size; i++)
      {
#if 0
        Int cur_nnz = cur_col_nnz[i];
        Int k_cur   = colptr[i] -1;
        Int k_new   = expColptr[i] -1;


        /* copy current non-zeros from old matrix to new matrix */
        std::copy(rowind.Data() + k_cur, rowind.Data() + k_cur + cur_nnz , expRowind.Data() + k_new);

        /* fill in the symmetric "missing" values */
        while (k_cur < colptr[i+1]-1)
        {
          /* non-zero of original matrix */
          Int j = rowind[k_cur]-1;

          if (j != i)  /* if not a non-diagonal element ... */
          {
            /* position of this transposed element in new matrix */
            k_new = expColptr[j]-1 + new_col_nnz[j];

            /* store */
            expRowind[k_new] = i+1;
            /*  update so next element stored at row j will appear
             *  at the right place.  */
            new_col_nnz[j]++;
          }

          k_cur++;
        }
#else
        Int cur_nnz = cur_col_nnz[i];
        Int k_cur   = colptr[i] -1;
        Int k_new;



        /* fill in the symmetric "missing" values */
        while (k_cur < colptr[i+1]-1)
        {
          /* non-zero of original matrix */
          Int j = rowind[k_cur]-1;

          if (j != i)  /* if not a non-diagonal element ... */
          {
            /* position of this transposed element in new matrix */
            k_new = expColptr[j]-1 + new_col_nnz[j];

            /* store */
            expRowind[k_new] = i+1;
            /*  update so next element stored at row j will appear
             *  at the right place.  */
            new_col_nnz[j]++;
          }

          k_cur++;
        }


        /* position of this transposed element in new matrix */
        k_new = expColptr[i] -1 + new_col_nnz[i];

        k_cur   = colptr[i] -1;
        /* copy current non-zeros from old matrix to new matrix */
        std::copy(&rowind[0] + k_cur, &rowind[0] + k_cur + cur_nnz , &expRowind[0] + k_new);
#endif
      }
      //expRowind[new_nnz-1] = size;

//logfileptr->OFS()<<expColptr<<endl;
//logfileptr->OFS()<<expRowind<<endl;




      }
      else{
        throw std::logic_error("Please use DistSparseMatrixGraph instead.");
        //logfileptr->OFS()<<"N is"<<this->size<<endl;
        //logfileptr->OFS()<<"baseval is"<<this->baseval<<endl;
        //logfileptr->OFS()<<"keepDiag is"<<this->keepDiag<<endl;
        //logfileptr->OFS()<<"sorted is"<<this->sorted<<endl;
        //int baseval =1;

        //this would be declared by sympack
        MPI_Op MPI_SYMPACK_PTR_SUM; 
        MPI_Datatype MPI_SYMPACK_PTR; 

        /* explain to MPI how type Complex is defined 
         */ 
        MPI_Type_contiguous( sizeof(Ptr), MPI_BYTE, &MPI_SYMPACK_PTR ); 
        MPI_Type_commit( &MPI_SYMPACK_PTR ); 
        /* create the complex-product user-op 
         */ 
        MPI_Op_create( _PtrSum, true, &MPI_SYMPACK_PTR_SUM ); 

        int mpirank,mpisize;
        int gmpirank,gmpisize;
        MPI_Comm expandComm = *workcomm;
        Idx N = size; 

        Ptr * pcolptr = &this->colptr[0];
        Idx * prowind = &this->rowind[0];
        Idx locColCnt = this->colptr.size()-1;
        Ptr locNNZ = this->rowind.size();

        MPI_Comm_size(expandComm,&gmpisize);
        MPI_Comm_rank(expandComm,&gmpirank);

        MPI_Comm splitcomm;
        MPI_Comm_split(expandComm,locColCnt>0,gmpirank,&splitcomm);


        if(locColCnt>0){

          MPI_Comm_size(splitcomm,&mpisize);
          MPI_Comm_rank(splitcomm,&mpirank);

          Idx colPerProc = N / mpisize;
          //first generate the expanded distributed structure
          {
            Int firstLocCol = mpirank>0?(mpirank)*colPerProc:0; //0 based
            Int maxLocN = max(locColCnt,N-(mpisize-1)*colPerProc); // can be 0

            //SYMPACK::vector<Ptr> extra_nnz_percol(maxLocN,0);
            SYMPACK::vector<Ptr> remote_colptr(maxLocN+1);
            SYMPACK::vector<Idx> remote_rowind;
            SYMPACK::vector<Ptr> remote_rowindPos(maxLocN+1);
            //SYMPACK::vector<Idx> extra_rowind(maxLocN,0);
            SYMPACK::vector<Ptr> curPos(locColCnt);
            SYMPACK::vector<Ptr> prevPos(locColCnt);

            std::copy(pcolptr,pcolptr+locColCnt,curPos.begin());



            for(Int prow = 0; prow<mpisize; prow++){
              Int firstRemoteCol = prow>0?(prow)*colPerProc:0; // 0 based
              Int pastLastRemoteCol = prow==mpisize-1?N:(prow+1)*colPerProc; // 0 based
              Ptr remColCnt = pastLastRemoteCol - firstRemoteCol;
              Ptr maxExtraNNZ = 0; 
              Ptr extraNNZ = 0;
              //receive from all the previous processor
              //receive extra nnzcnt...

              std::fill(remote_colptr.begin(),remote_colptr.end(),0);
              if(mpirank<=prow){
                //backup the position in each column
                std::copy(curPos.begin(),curPos.end(),prevPos.begin());
                //use remote_colptr to store the number of nnz per col first
                //loop through my local columns and figure out the extra nonzero per col on remote processors
                for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
                  Idx col = firstLocCol + locCol;  // 0 based
                  Ptr colbeg = curPos[locCol]-baseval; //now 0 based
                  Ptr colend = pcolptr[locCol+1]-baseval; // now 0 based
                  for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
                    Idx row = prowind[rptr]-baseval; //0 based
                    assert(row>=firstRemoteCol);
                    if(row>col){
                      //this goes onto prow
                      if(row<pastLastRemoteCol){
                        //this is shifted by one to compute the prefix sum without copying afterwards
                        remote_colptr[row-firstRemoteCol+1]++;
                        extraNNZ++;
                      }
                      else{
                        break;
                      }
                    }
                    curPos[locCol]++; // baseval based
                  }
                }

#ifdef DEBUG
                {
                  Ptr checkNNZ = 0;
                  for(Idx col = 1;col<=remColCnt;col++){
                    checkNNZ += remote_colptr[col];
                  }
                  assert(extraNNZ==checkNNZ);
                }
#endif

                //put the local sum into the first element of the array
                remote_colptr[0] = 0;
                for(Idx p = 1; p<remColCnt+1;p++){ remote_colptr[0] += remote_colptr[p];}

                if(mpirank==prow){
                  //we can now receive the number of NZ per col into the expanded pcolptr
                  expColptr.resize(remColCnt+1,0);
                  //this array will contain the max element in our custom reduce
                  expColptr[0] = 0;
                  TIMER_START(REDUCE);
#ifdef USE_REDUCE
                  MPI_Reduce(&remote_colptr[0],&expColptr[0],remColCnt+1,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_SUM,prow,splitcomm);
#else
                  {
                    SYMPACK::vector<Ptr> recv(remColCnt+1);
                    //first copy remote_colptr in expColptr
                    int len = remColCnt+1;
                    _PtrSum(&remote_colptr[0],&expColptr[0],&len,NULL);
                    for(Int pcol = 0; pcol<prow; pcol++){
                      MPI_Status status;
                      MPI_Recv(&recv[0],(remColCnt+1)*sizeof(Ptr),MPI_BYTE,MPI_ANY_SOURCE,mpisize+prow,splitcomm,&status);
                      _PtrSum(&recv[0],&expColptr[0],&len,NULL);
                    }
                  }
#endif
                  TIMER_STOP(REDUCE);

                  maxExtraNNZ = expColptr[0];
                  remote_rowind.resize(maxExtraNNZ);

                  for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
                    Ptr colbeg = pcolptr[locCol]-baseval; //now 0 based
                    Ptr colend = pcolptr[locCol+1]-baseval; // now 0 based
                    expColptr[locCol+1] += colend - colbeg - 1 + keepDiag;
                    //At this point expColptr[locCol+1] contains the total number of NNZ of locCol 
                    //we can now compute the expanded pcolptr
                  }

                  expColptr[0] = baseval;
                  for(Idx col = 1;col<remColCnt+1;col++){ expColptr[col]+=expColptr[col-1]; }
                }
                else{
                  TIMER_START(REDUCE);
#ifdef USE_REDUCE
                  MPI_Reduce(&remote_colptr[0],NULL,remColCnt+1,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_SUM,prow,splitcomm);
#else
                  {
                    MPI_Send(&remote_colptr[0],(remColCnt+1)*sizeof(Ptr),MPI_BYTE,prow,mpisize+prow,splitcomm);
                  }
#endif
                  TIMER_STOP(REDUCE);
                  remote_rowind.resize(extraNNZ);
                }



                /**************     Compute remote_colptr from the local nnz per col ***********/
                //compute a prefix sum of the nnz per column to get the new pcolptr
                remote_colptr[0] = baseval;
                for(Idx col = 1;col<=remColCnt;col++){ remote_colptr[col]+=remote_colptr[col-1]; }

                /**************     Fill remote_rowind now ****************/
                //make a copy of pcolptr in order to track the current position in each column
                std::copy(&remote_colptr[0],&remote_colptr[0]+remColCnt+1,remote_rowindPos.begin());
                //loop through my local columns and figure out the extra nonzero per col on remote processors
                for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
                  Idx col = firstLocCol + locCol;  // 0 based
                  Ptr colbeg = prevPos[locCol]-baseval; //now 0 based
                  Ptr colend = curPos[locCol]-baseval; // now 0 based

                  for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
                    Idx row = prowind[rptr]-baseval; //0 based
                    assert(row>=firstRemoteCol && row<pastLastRemoteCol);
                    if(row>col){
                      //this goes onto prow
                      Idx locRow = row - firstRemoteCol;      
                      remote_rowind[ remote_rowindPos[locRow]++ - baseval ] = col + baseval;
                    }
                  }
                }

#ifdef DEBUG
                logfileptr->OFS()<<"remote_colptr ";
                for(Idx col = 0;col<=remColCnt;col++){
                  logfileptr->OFS()<<remote_colptr[col]<<" ";
                }
                logfileptr->OFS()<<endl;

                logfileptr->OFS()<<"remote_rowind ";
                for(Ptr col = 0;col<extraNNZ;col++){
                  logfileptr->OFS()<<remote_rowind[col]<<" ";
                }
                logfileptr->OFS()<<endl;
#endif

                if(prow==mpirank){
#ifdef DEBUG
                  logfileptr->OFS()<<"expColptr "<<expColptr<<endl;
                  //logfileptr->OFS()<<"true expColptr "<<Global.expColptr<<endl;
#endif
                  expRowind.resize(expColptr.back()-baseval);//totalExtraNNZ+locNNZ-(1-keepDiag)*locColCnt);
                  std::copy(expColptr.begin(),expColptr.end(),remote_rowindPos.begin()); 
                  for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
                    Idx col = firstLocCol + locCol;  // 0 based
                    Ptr colbeg = pcolptr[locCol]-baseval; //now 0 based
                    Ptr colend = pcolptr[locCol+1]-baseval; // now 0 based
                    Ptr & pos = remote_rowindPos[locCol];

                    //copy the local lower triangular part into the expanded structure
                    for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
                      Idx row = prowind[rptr]-baseval; //0 based
                      if(col!=row || keepDiag){
                        expRowind[pos++ - baseval] = row + baseval;  
                      }
                    }

                    //copy the local extra NNZ into the expanded structure
                    colbeg = remote_colptr[locCol]-baseval; //now 0 based
                    colend = remote_colptr[locCol+1]-baseval; // now 0 based 
                    for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
                      Idx row = remote_rowind[rptr]-baseval; //0 based
                      expRowind[pos++ - baseval] = row + baseval;  
                    }
                  }

                  TIMER_START(RECV);
                  for(Int pcol = 0; pcol<prow; pcol++){
                    //Use an MPI_Gatherv instead ? >> memory usage : p * n/p
                    //Do mpi_recv from any, anytag for pcolptr and then do the matching rowind ?

                    //logfileptr->OFS()<<"P"<<mpirank<<" receives pcolptr from P"<<pcol<<endl;
                    //receive colptrs...
                    MPI_Status status;
                    MPI_Recv(&remote_colptr[0],(remColCnt+1)*sizeof(Ptr),MPI_BYTE,MPI_ANY_SOURCE,prow,splitcomm,&status);

                    //logfileptr->OFS()<<"P"<<mpirank<<" receives rowind from P"<<pcol<<endl;
                    //receive rowinds...
                    MPI_Recv(&remote_rowind[0],maxExtraNNZ*sizeof(Idx),MPI_BYTE,status.MPI_SOURCE,prow,splitcomm,MPI_STATUS_IGNORE);

                    TIMER_START(PROCESSING_RECV_DATA);
                    //logfileptr->OFS()<<"P"<<mpirank<<" done receiving from P"<<pcol<<endl;
                    for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
                      Idx col = firstLocCol + locCol;  // 0 based
                      //copy the extra NNZ into the expanded structure
                      Ptr colbeg = remote_colptr[locCol]-baseval; //now 0 based
                      Ptr colend = remote_colptr[locCol+1]-baseval; // now 0 based 
                      Ptr & pos = remote_rowindPos[locCol];
                      for(Ptr rptr = colbeg ; rptr< colend ; rptr++ ){
                        Idx row = remote_rowind[rptr]-baseval; //0 based
                        //perColRowind[locCol].push_back(row+baseval);
                        expRowind[pos++ - baseval] = row + baseval;  
                      }
                    }
                    TIMER_STOP(PROCESSING_RECV_DATA);
                  }
                  TIMER_STOP(RECV);

                  if(sorted){ 
                    for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
                      Ptr colbeg = expColptr[locCol]-baseval; //now 0 based
                      Ptr colend = expColptr[locCol+1]-baseval; // now 0 based 
                      sort(&expRowind[0]+colbeg,&expRowind[0]+colend,std::less<Ptr>());
                    }
                  }

#ifdef DEBUG
                  logfileptr->OFS()<<"expRowind "<<expRowind<<endl;
                  //logfileptr->OFS()<<"true expRowind "<<Global.expRowind<<endl;
#endif

                }
                else{
                  TIMER_START(SEND);
                  //MPI_Send
                  //            logfileptr->OFS()<<"P"<<mpirank<<" sends pcolptr to P"<<prow<<endl;
                  MPI_Send(&remote_colptr[0],(remColCnt+1)*sizeof(Ptr),MPI_BYTE,prow,prow,splitcomm);
                  //            logfileptr->OFS()<<"P"<<mpirank<<" sends rowind to P"<<prow<<endl;
                  MPI_Send(&remote_rowind[0],extraNNZ*sizeof(Idx),MPI_BYTE,prow,prow,splitcomm);
                  //            logfileptr->OFS()<<"P"<<mpirank<<" done sending to P"<<prow<<endl;
                  TIMER_STOP(SEND);
                }

              }
              else{
                TIMER_START(REDUCE);
#ifdef USE_REDUCE
                MPI_Reduce(&remote_colptr[0],NULL,remColCnt+1,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_SUM,prow,splitcomm);
#endif
                TIMER_STOP(REDUCE);
              }

            }
          }
        }
        else{
          expColptr.resize(0);
          expRowind.resize(0);
        }

        MPI_Op_free(&MPI_SYMPACK_PTR_SUM);
        MPI_Type_free(&MPI_SYMPACK_PTR);

        nnz = expRowind.size();



      }
      bIsExpanded =true;
    }
    TIMER_STOP(EXPAND);

  }


  void SparseMatrixStructure::ToGlobal(SparseMatrixStructure & pGlobal,MPI_Comm & comm){

    TIMER_START(ToGlobalStructure);
    // test if structure hasn't been allocated yet
    if(bIsGlobal){
      pGlobal = *this;
    }
    else{


      int np;
      int iam;

      int ismpi=0;
      MPI_Initialized( &ismpi);

      //FIXME needs to be passed as an argument ?
      //MPI_Comm comm = MPI_COMM_WORLD;

      int isnull= (comm == MPI_COMM_NULL);
      //    logfileptr->OFS()<<ismpi<<std::endl;
      //    logfileptr->OFS()<<comm<<std::endl;
      //    logfileptr->OFS()<<MPI_COMM_NULL<<std::endl;
      //    logfileptr->OFS()<<isnull<<std::endl;

      if(ismpi && isnull==0){
        MPI_Comm_size(comm,&np);
        MPI_Comm_rank(comm, &iam);
      }
      else{
        //throw an exception
        throw std::logic_error("MPI needs to be available.");
      }


      pGlobal.bIsGlobal = true;
      pGlobal.bIsExpanded = false;
      pGlobal.size = size;
      pGlobal.colptr.resize(size+1);

      {
        /* Allgatherv for row indices. */ 
        SYMPACK::vector<int> prevnz(np);
        SYMPACK::vector<int> rcounts(np);
        MPI_Allgather(&nnz, sizeof(nnz), MPI_BYTE, &rcounts[0], sizeof(nnz), MPI_BYTE, comm);

        prevnz[0] = 0;
        for (Int i = 0; i < np-1; ++i) { prevnz[i+1] = prevnz[i] + rcounts[i]; } 

        pGlobal.nnz = 0;
        for (Int i = 0; i < np; ++i) { pGlobal.nnz += rcounts[i]; } 
        pGlobal.rowind.resize(pGlobal.nnz);

        SYMPACK::vector<int> rdispls = prevnz;
        for (Int i = 0; i < np; ++i) { rcounts[i] *= sizeof(Idx); } 
        for (Int i = 0; i < np; ++i) { rdispls[i] *= sizeof(Idx); } 

        //    logfileptr->OFS()<<"Global nnz is "<<pGlobal.nnz<<std::endl;

        MPI_Allgatherv(&rowind[0], nnz*sizeof(Idx), MPI_BYTE, &pGlobal.rowind[0],&rcounts[0], &rdispls[0], MPI_BYTE, comm); 

        //    logfileptr->OFS()<<"Global rowind is "<<pGlobal.rowind<<std::endl;

        /* Allgatherv for colptr */
        // Compute the number of columns on each processor
        Int numColFirst = std::max(1,size / np);
        std::fill(rcounts.begin(),rcounts.end(),numColFirst*sizeof(Ptr));
        rcounts[np-1] = (size - numColFirst * (np-1))*sizeof(Ptr);  // Modify the last entry     

        rdispls[0] = 0;
        for (Int i = 0; i < np-1; ++i) { rdispls[i+1] = rdispls[i] + rcounts[i]; } 

        MPI_Allgatherv(&colptr[0], (colptr.size()-1)*sizeof(Ptr), MPI_BYTE, &pGlobal.colptr[0],
            &rcounts[0], &rdispls[0], MPI_BYTE, comm);

        for (Int i = 0; i < np; ++i) { rcounts[i] /= sizeof(Ptr); } 
        for (Int i = 0; i < np; ++i) { rdispls[i] /= sizeof(Ptr); } 
        /* Recompute column pointers. */
        for (Int p = 1; p < np; p++) {
          Int idx = rdispls[p];
          for (Int j = 0; j < rcounts[p]; ++j) pGlobal.colptr[idx++] += prevnz[p];
        }
        pGlobal.colptr[pGlobal.size]= pGlobal.nnz+1;
      }

      //    logfileptr->OFS()<<"Global colptr is "<<pGlobal.colptr<<std::endl;
      if(bIsExpanded){
        pGlobal.bIsExpanded = true;
        pGlobal.expColptr.resize(size+1);


        /* Allgatherv for row indices. */ 
        SYMPACK::vector<int> prevnz(np);
        SYMPACK::vector<int> rcounts(np);
        Ptr gnnz = expRowind.size();
        MPI_Allgather(&gnnz, sizeof(gnnz), MPI_BYTE, &rcounts[0], sizeof(gnnz), MPI_BYTE, comm);

        prevnz[0] = 0;
        for (Int i = 0; i < np-1; ++i) { prevnz[i+1] = prevnz[i] + rcounts[i]; } 

        Ptr totnnz=0;
        for (Int i = 0; i < np; ++i) { totnnz += rcounts[i]; } 
        pGlobal.expRowind.resize(totnnz);

        SYMPACK::vector<int> rdispls = prevnz;
        for (Int i = 0; i < np; ++i) { rcounts[i] *= sizeof(Idx); } 
        for (Int i = 0; i < np; ++i) { rdispls[i] *= sizeof(Idx); } 

        MPI_Allgatherv(&expRowind[0], gnnz*sizeof(Idx), MPI_BYTE, &pGlobal.expRowind[0],&rcounts[0], &rdispls[0], MPI_BYTE, comm); 

        /* Allgatherv for colptr */
        // Compute the number of columns on each processor
        Int numColFirst = std::max(1,size / np);
        std::fill(rcounts.begin(),rcounts.end(),numColFirst*sizeof(Ptr));
        rcounts[np-1] = (size - numColFirst * (np-1))*sizeof(Ptr);  // Modify the last entry     

        rdispls[0] = 0;
        for (Int i = 0; i < np-1; ++i) { rdispls[i+1] = rdispls[i] + rcounts[i]; } 

        MPI_Allgatherv(&expColptr[0], (expColptr.size()-1)*sizeof(Ptr), MPI_BYTE, &pGlobal.expColptr[0],
            &rcounts[0], &rdispls[0], MPI_BYTE, comm);

        for (Int i = 0; i < np; ++i) { rcounts[i] /= sizeof(Ptr); } 
        for (Int i = 0; i < np; ++i) { rdispls[i] /= sizeof(Ptr); } 
        /* Recompute column pointers. */
        for (Int p = 1; p < np; p++) {
          Int idx = rdispls[p];
          for (Int j = 0; j < rcounts[p]; ++j) pGlobal.expColptr[idx++] += prevnz[p];
        }
        pGlobal.expColptr[pGlobal.size]= totnnz+1;


      }

    }

    TIMER_STOP(ToGlobalStructure);
  }

  void SparseMatrixStructure::GetLColRowCount(ETree & tree, Ordering & aOrder, SYMPACK::vector<Int> & cc, SYMPACK::vector<Int> & rc){
    //The tree need to be postordered
    if(!tree.IsPostOrdered()){
      tree.PostOrderTree(aOrder);
    }

    ExpandSymmetric();

    TIMER_START(GetColRowCount_Classic);
    cc.resize(size);
    rc.resize(size);

    SYMPACK::vector<Idx> level(size+1);
    SYMPACK::vector<Int> weight(size+1);
    SYMPACK::vector<Idx> fdesc(size+1);
    SYMPACK::vector<Idx> nchild(size+1);
    SYMPACK::vector<Idx> set(size);
    SYMPACK::vector<Idx> prvlf(size);
    SYMPACK::vector<Idx> prvnbr(size);


    Idx xsup = 1;
    level[0] = 0;
    for(Idx k = size; k>=1; --k){
      rc[k-1] = 1;
      cc[k-1] = 0;
      set[k-1] = k;
      prvlf[k-1] = 0;
      prvnbr[k-1] = 0;
      level[k] = level[tree.PostParent(k-1)] + 1;
      weight[k] = 1;
      fdesc[k] = k;
      nchild[k] = 0;
    }

    nchild[0] = 0;
    fdesc[0] = 0;
    for(Idx k =1; k<size; ++k){
      Idx parent = tree.PostParent(k-1);
      weight[parent] = 0;
      ++nchild[parent];
      Idx ifdesc = fdesc[k];
      if  ( ifdesc < fdesc[parent] ) {
        fdesc[parent] = ifdesc;
      }
    }

    for(Idx lownbr = 1; lownbr<=size; ++lownbr){
      Int lflag = 0;
      Idx ifdesc = fdesc[lownbr];
      Idx oldnbr = aOrder.perm[lownbr-1];
      Ptr jstrt = expColptr[oldnbr-1];
      Ptr jstop = expColptr[oldnbr] - 1;


      //           -----------------------------------------------
      //           for each ``high neighbor'', hinbr of lownbr ...
      //           -----------------------------------------------
      for(Ptr j = jstrt; j<=jstop;++j){
        Idx hinbr = expRowind[j-1];
        hinbr = aOrder.invp[hinbr-1];
        if  ( hinbr > lownbr )  {
          if  ( ifdesc > prvnbr[hinbr-1] ) {
            //                       -------------------------
            //                       increment weight[lownbr].
            //                       -------------------------
            ++weight[lownbr];
            Idx pleaf = prvlf[hinbr-1];
            //                       -----------------------------------------
            //                       if hinbr has no previous ``low neighbor'' 
            //                       then ...
            //                       -----------------------------------------
            if  ( pleaf == 0 ) {
              //                           -----------------------------------------
              //                           ... accumulate lownbr-->hinbr path length 
              //                               in rowcnt[hinbr].
              //                           -----------------------------------------
              rc[hinbr-1] += level[lownbr] - level[hinbr];
            }
            else{
              //                           -----------------------------------------
              //                           ... otherwise, lca <-- find[pleaf], which 
              //                               is the least common ancestor of pleaf 
              //                               and lownbr.
              //                               (path halving.)
              //                           -----------------------------------------
              Idx last1 = pleaf;
              Idx last2 = set[last1-1];
              Idx lca = set[last2-1];
              while(lca != last2){
                set[last1-1] = lca;
                last1 = lca;
                last2 = set[last1-1];
                lca = set[last2-1];
              }
              //                           -------------------------------------
              //                           accumulate pleaf-->lca path length in 
              //                           rowcnt[hinbr].
              //                           decrement weight(lca).
              //                           -------------------------------------
              rc[hinbr-1] += level[lownbr] - level[lca];
              --weight[lca];
            }
            //                       ----------------------------------------------
            //                       lownbr now becomes ``previous leaf'' of hinbr.
            //                       ----------------------------------------------
            prvlf[hinbr-1] = lownbr;
            lflag = 1;
          }
          //                   --------------------------------------------------
          //                   lownbr now becomes ``previous neighbor'' of hinbr.
          //                   --------------------------------------------------
          prvnbr[hinbr-1] = lownbr;
        }
      }
      //           ----------------------------------------------------
      //           decrement weight ( parent[lownbr] ).
      //           set ( p[lownbr] ) <-- set ( p[lownbr] ) + set[xsup].
      //           ----------------------------------------------------
      Idx parent = tree.PostParent(lownbr-1);
      --weight[parent];


      //merge the sets
      if  ( lflag == 1  || nchild[lownbr] >= 2 ) {
        xsup = lownbr;
      }
      set[xsup-1] = parent;
    }



#ifdef _DEBUG_
    logfileptr->OFS()<<"deltas "<<weight<<std::endl;
#endif

    for(Int k = 1; k<=size; ++k){
      Int temp = cc[k-1] + weight[k];
      cc[k-1] = temp;
      Int parent = tree.PostParent(k-1);
      if  ( parent != 0 ) {
        cc[parent-1] += temp;
      }
    }

    //      logfileptr->OFS()<<"column counts "<<cc<<std::endl;

    TIMER_STOP(GetColRowCount_Classic);
  }



  void SparseMatrixStructure::FindSupernodes(ETree& tree, Ordering & aOrder, SYMPACK::vector<Int> & cc,SYMPACK::vector<Int> & supMembership, SYMPACK::vector<Int> & xsuper, Int maxSize ){
    TIMER_START(FindSupernodes);

    if(!bIsGlobal){
      throw std::logic_error( "SparseMatrixStructure must be global in order to call FindSupernodes\n" );
    }

    supMembership.resize(size);

    Int nsuper = 1;
    Int supsize = 1;
    supMembership[0] = 1;

    for(Int i =2; i<=size;i++){
      Int prev_parent = tree.PostParent(i-2);
      if(prev_parent == i){
        if(cc[i-2] == cc[i-1]+1 ) {
          if(supsize<maxSize || maxSize==-1){
            ++supsize;
            supMembership[i-1] = nsuper;
            continue;
          }
        }
      }

      nsuper++;
      supsize = 1;
      supMembership[i-1] = nsuper;
    }

    xsuper.resize(nsuper+1);
    Int lstsup = nsuper+1;
    for(Int i = size; i>=1;--i){
      Int ksup = supMembership[i-1];
      if(ksup!=lstsup){
        xsuper[lstsup-1] = i + 1; 
      }
      lstsup = ksup;
    }
    xsuper[0]=1;
    TIMER_STOP(FindSupernodes);
  }

#ifdef REFINE_SNODE
  //EXPERIMENTAL STUFF
  typedef std::set<Int> nodeset;
  typedef std::list<nodeset*> partitions;
  typedef SYMPACK::vector<nodeset*> vecset;
  void SparseMatrixStructure::RefineSupernodes(ETree& tree, Ordering & aOrder, SYMPACK::vector<Int> & supMembership, SYMPACK::vector<Int> & xsuper, PtrVec & xlindx, IdxVec & lindx, SYMPACK::vector<Int> & perm){

    perm.resize(size);

    Int nsuper = xsuper.size()-1;
    //Using std datatypes first

    vecset adj(nsuper);
    vecset snodes(nsuper);

    partitions L;

    SYMPACK::vector<Int> origPerm(size);

    //init L with curent supernodal partition
    Int pos = 1;
    for(Int i = 1; i<=nsuper; ++i){
      adj[i-1] = new nodeset();
      Ptr fi = xlindx[i-1];
      Ptr li = xlindx[i]-1;
      for(Ptr idx = fi; idx<=li;idx++){
        adj[i-1]->insert(lindx[idx-1]);
      }

      L.push_back(new nodeset());
      snodes[i-1] = L.back();
      nodeset * curL = L.back();
      Int fc = xsuper[i-1];
      Int lc = xsuper[i]-1;
      for(Int node = fc; node<=lc;++node){
        curL->insert(node);
        origPerm[node-1] = pos++;
      }
    }

    Int K = nsuper;
    partitions::reverse_iterator Kit;
    for(Kit = L.rbegin();Kit != L.rend(); ++Kit){
      logfileptr->OFS()<<"Looking at snode "<<K<<std::endl;

      assert( snodes[K-1] == *Kit);

      //        logfileptr->OFS()<<"Adj is "<<*adj[K-1]<<std::endl;

      partitions::reverse_iterator Jit;
      //    partitions tmp;
      Jit = L.rbegin();
      Int count = L.size() - K;
      while(count>0){
        //        logfileptr->OFS()<<"L is "<<**Jit<<std::endl;
        //for(Jit = L.rbegin();Jit!=Kit;++Jit){
        nodeset * inter = new nodeset();

        std::set_intersection((*Jit)->begin(),(*Jit)->end(),
            adj[K-1]->begin(),adj[K-1]->end(),
            std::inserter(*inter, inter->begin()));
        //        logfileptr->OFS()<<"Intersect is "<<*inter<<std::endl;

        if(inter->size()>0){
          nodeset * diff = new nodeset();
          std::set_difference((*Jit)->begin(),(*Jit)->end(),
              adj[K-1]->begin(),adj[K-1]->end(),
              std::inserter(*diff, diff->begin()));

          //        logfileptr->OFS()<<"Diff is "<<*diff<<std::endl;

          if(diff->size()>0){
            //          tmp.push_back(diff);
            //          tmp.push_back(inter);
            //replace Jit by inter and diff
            (*Jit)->swap(*diff);
            L.insert(Jit.base(),inter);
            delete diff;
            //          (*Jit)->swap(*inter);
            //          L.insert(Jit.base(),diff);
            //          delete inter;
            ++Jit; 
            ++Jit; 
          }
          else{
            //          tmp.push_back(*Jit);
            delete diff;
          }
        }
        else{
          //        tmp.push_back(*Jit);
          delete inter;
        }
        //}
        ++Jit;
        --count;
      }


      partitions::iterator it;
      for(it = L.begin();it != L.end(); ++it){
        logfileptr->OFS()<<"[ ";
        for(nodeset::iterator nit = (*it)->begin();nit != (*it)->end(); ++nit){
          logfileptr->OFS()<<*nit<<" ";
        }
        logfileptr->OFS()<<"] ";
      }
      logfileptr->OFS()<<std::endl;
      --K;
    }

    partitions::iterator it;
    for(it = L.begin();it != L.end(); ++it){
      logfileptr->OFS()<<"[ ";
      for(nodeset::iterator nit = (*it)->begin();nit != (*it)->end(); ++nit){
        logfileptr->OFS()<<*nit<<" ";
      }
      logfileptr->OFS()<<"] ";
    }
    logfileptr->OFS()<<std::endl;


    //construct perm
    pos = 1;
    for(it = L.begin();it != L.end(); ++it){
      for(nodeset::iterator nit = (*it)->begin();nit != (*it)->end(); ++nit){
        perm[*nit-1] = pos++;
      }
    }

    logfileptr->OFS()<<"Orig col order "<<origPerm<<std::endl;
    logfileptr->OFS()<<"Refined col order "<<perm<<std::endl;

    //  SYMPACK::vector<Int> lindxTemp = lindx;
    //  //change lindx to reflect the new ordering
    //  for(Int i = 1; i<=nsuper; ++i){
    //    Int fi = xlindx(i-1);
    //    Int li = xlindx(i)-1;
    //    for(Int idx = fi; idx<=li;idx++){
    //       lindx[idx-1] = perm[lindxTemp[idx-1]-1];
    //    }
    //  }
    //
    //    logfileptr->OFS()<<"Previous lindx "<<lindxTemp<<std::endl;
    //    logfileptr->OFS()<<"Refined lindx "<<lindx<<std::endl;


    for(it = L.begin();it != L.end(); ++it){
      delete (*it);
    }

    for(Int i = 1; i<=nsuper; ++i){
      delete adj[i-1];
    }

    aOrder.Compose(perm);

  }
#endif


  void SparseMatrixStructure::RelaxSupernodes(ETree& tree, SYMPACK::vector<Int> & cc,SYMPACK::vector<Int> & supMembership, SYMPACK::vector<Int> & xsuper, RelaxationParameters & params ){

    Int nsuper = xsuper.size()-1;

    DisjointSet sets;
    sets.Initialize(nsuper);
    SYMPACK::vector<Int> ncols(nsuper);
    SYMPACK::vector<Int> zeros(nsuper);
    SYMPACK::vector<Int> newCC(nsuper);
    for(Int ksup=nsuper;ksup>=1;--ksup){
      Int cset = sets.makeSet(ksup);
      sets.Root(cset-1)=ksup;

      Int fstcol = xsuper[ksup-1];
      Int lstcol = xsuper[ksup]-1;
      Int width = lstcol - fstcol +1;
      Int length = cc[fstcol-1];
      ncols[ksup-1] = width;
      zeros[ksup-1] = 0;
      newCC[ksup-1] = length;
    }



    for(Int ksup=nsuper;ksup>=1;--ksup){
      Int fstcol = xsuper[ksup-1];
      Int lstcol = xsuper[ksup]-1;
      Int width = ncols[ksup-1];
      Int length = cc[fstcol-1];

      Int parent_fstcol = tree.PostParent(lstcol-1);
      if(parent_fstcol!=0){
        Int parent_snode = supMembership[parent_fstcol-1];
        Int pset = sets.find(parent_snode);
        parent_snode = sets.Root(pset-1);

        bool merge = (parent_snode == ksup+1);


        if(merge){
          Int parent_width = ncols[parent_snode-1];

          Int parent_fstcol = xsuper[parent_snode-1];
          Int parent_lstcol = xsuper[parent_snode]-1;
          Int totzeros = zeros[parent_snode-1];
          Int fused_cols = width + parent_width;

          merge = false;
          if(fused_cols <= params.nrelax0){
            merge = true;
          }
          else if(fused_cols <=params.maxSize){
            double child_lnz = cc[fstcol-1];
            double parent_lnz = cc[parent_fstcol-1];
            double xnewzeros = width * (parent_lnz + width  - child_lnz);

            if(xnewzeros == 0){
              merge = true;
            }
            else{
              //all these values are the values corresponding to the merged snode
              double xtotzeros = (double)totzeros + xnewzeros;
              double xfused_cols = (double) fused_cols;
              //new number of nz
              double xtotsize = (xfused_cols * (xfused_cols+1)/2) + xfused_cols * (parent_lnz - parent_width);
              //percentage of explicit zeros
              double z = xtotzeros / xtotsize;

              Int totsize = (fused_cols * (fused_cols+1)/2) + fused_cols * ((Int)parent_lnz - parent_width);
              totzeros += (Int)xnewzeros;

              merge = ((fused_cols <= params.nrelax1 && z < params.zrelax0) 
                  || (fused_cols <= params.nrelax2 && z < params.zrelax1)
                  || (z<params.zrelax2)) &&
                (xtotsize < std::numeric_limits<Int>::max() / sizeof(double));
            }

          }

          if(merge){
            //              std::cout<<"merge "<<ksup<<" and "<<parent_snode<<std::endl;
            ncols[ksup-1] += ncols[parent_snode-1]; 
            zeros[ksup-1] = totzeros;
            newCC[ksup-1] = width + newCC[parent_snode-1];
            sets.Union(ksup,parent_snode,ksup);
          }
        } 

      }
    }

    SYMPACK::vector<Int> relXSuper(nsuper+1);
    Int nrSuper = 0;
    for(Int ksup=1;ksup<=nsuper;++ksup){
      Int kset = sets.find(ksup);
      if(ksup == sets.Root(kset-1)){
        Int fstcol = xsuper[ksup-1];
        relXSuper[nrSuper] = fstcol;
        newCC[nrSuper] = newCC[ksup-1];
        ++nrSuper;
      }
    }
    relXSuper[nrSuper] = xsuper[nsuper];
    relXSuper.resize(nrSuper+1);

    for(Int ksup=1;ksup<=nrSuper;++ksup){
      Int fstcol = relXSuper[ksup-1];
      Int lstcol = relXSuper[ksup]-1;
      for(Int col = fstcol; col<=lstcol;++col){
        supMembership[col-1] = ksup;
        cc[col-1] = newCC[ksup-1] + col-fstcol;
      }
    }

    xsuper = relXSuper;
    ///      //adjust the column counts
    ///      for(Int col=i-2;col>=i-supsize;--col){
    ///        cc[col-1] = cc[col]+1;
    ///      }


  }


//  void SparseMatrixStructure::SymbolicFactorizationRelaxed_upcxx(ETree& tree,Ordering & aOrder, const SYMPACK::vector<Int> & cc,const SYMPACK::vector<Int> & xsuper,const SYMPACK::vector<Int> & SupMembership, upcxx::shared_array<Ptr> & xlindx, upcxx::shared_array<Idx> & lindx){
//    TIMER_START(SymbolicFactorization);
//
//
//    if(!bIsGlobal){
//      throw std::logic_error( "SparseMatrixStructure must be global in order to call SymbolicFactorization\n" );
//    }
//
//    Int nsuper = xsuper.size()-1;
//
//
//
//
//    Ptr nzbeg = 0;
//    //nzend points to the last used slot in lindx
//    Ptr nzend = 0;
//
//    //tail is the end of list indicator (in rchlnk, not mrglnk)
//    Idx tail = size +1;
//
//    Idx head = 0;
//
//    //Array of length nsuper containing the children of 
//    //each supernode as a linked list
//    SYMPACK::vector<Idx> mrglnk(nsuper,0);
//
//    //Array of length n+1 containing the current linked list 
//    //of merged indices (the "reach" set)
//    SYMPACK::vector<Idx> rchlnk(size+1);
//
//    //Array of length n used to mark indices as they are introduced
//    // into each supernode's index set
//    SYMPACK::vector<Int> marker(size,0);
//
//
//    xlindx.init(nsuper+1);
//
//    //Compute the sum of the column count and resize lindx accordingly
//    Ptr nofsub = 1;
//    for(Int i =0; i<cc.size();++i){
//      nofsub+=cc[i];
//    }
//
//    lindx.init(nofsub);
//
//
//    Ptr point = 1;
//    for(Int ksup = 1; ksup<=nsuper; ++ksup){
//      Int fstcol = xsuper[ksup-1];
//      xlindx[ksup-1] = point;
//      point += cc[fstcol-1]; 
//    } 
//    xlindx[nsuper] = point;
//
//
//
//    for(Int ksup = 1; ksup<=nsuper; ++ksup){
//      Int fstcol = xsuper[ksup-1];
//      Int lstcol = xsuper[ksup]-1;
//      Int width = lstcol - fstcol +1;
//      Int length = cc[fstcol-1];
//      Ptr knz = 0;
//      rchlnk[head] = tail;
//      Int jsup = mrglnk[ksup-1];
//
//      //If ksup has children in the supernodal e-tree
//      if(jsup>0){
//        //copy the indices of the first child jsup into 
//        //the linked list, and mark each with the value 
//        //ksup.
//        Int jwidth = xsuper[jsup]-xsuper[jsup-1];
//        Ptr jnzbeg = xlindx[jsup-1] + jwidth;
//        Ptr jnzend = xlindx[jsup] -1;
//        for(Ptr jptr = jnzend; jptr>=jnzbeg; --jptr){
//          Idx newi = lindx[jptr-1];
//          ++knz;
//          marker[newi-1] = ksup;
//          rchlnk[newi] = rchlnk[head];
//          rchlnk[head] = newi;
//        }
//
//        //for each subsequent child jsup of ksup ...
//        jsup = mrglnk[jsup-1];
//        while(jsup!=0 && knz < length){
//          //merge the indices of jsup into the list,
//          //and mark new indices with value ksup.
//
//          jwidth = xsuper[jsup]-xsuper[jsup-1];
//          jnzbeg = xlindx[jsup-1] + jwidth;
//          jnzend = xlindx[jsup] -1;
//          Int nexti = head;
//          for(Ptr jptr = jnzbeg; jptr<=jnzend; ++jptr){
//            Idx newi = lindx[jptr-1];
//            Idx i;
//            do{
//              i = nexti;
//              nexti = rchlnk[i];
//            }while(newi > nexti);
//
//            if(newi < nexti){
//#ifdef _DEBUG_
//              logfileptr->OFS()<<jsup<<" is a child of "<<ksup<<" and "<<newi<<" is inserted in the structure of "<<ksup<<std::endl;
//#endif
//              ++knz;
//              rchlnk[i] = newi;
//              rchlnk[newi] = nexti;
//              marker[newi-1] = ksup;
//              nexti = newi;
//            }
//          }
//          jsup = mrglnk[jsup-1];
//        }
//      }
//
//      //structure of a(*,fstcol) has not been examined yet.  
//      //"sort" its structure into the linked list,
//      //inserting only those indices not already in the
//      //list.
//      if(knz < length){
//        for(Int row = fstcol; row<=lstcol; ++row){
//          Idx newi = row;
//          if(newi > fstcol && marker[newi-1] != ksup){
//            //position and insert newi in list and
//            // mark it with kcol
//            Idx nexti = head;
//            Idx i;
//            do{
//              i = nexti;
//              nexti = rchlnk[i];
//            }while(newi > nexti);
//            ++knz;
//            rchlnk[i] = newi;
//            rchlnk[newi] = nexti;
//            marker[newi-1] = ksup;
//          }
//        }
//
//
//        for(Int col = fstcol; col<=lstcol; ++col){
//          Int node = aOrder.perm[col-1];
//
//          Ptr knzbeg = expColptr[node-1];
//          Ptr knzend = expColptr[node]-1;
//          for(Ptr kptr = knzbeg; kptr<=knzend;++kptr){
//            Idx newi = expRowind[kptr-1];
//            newi = aOrder.invp[newi-1];
//
//            if(newi > fstcol && marker[newi-1] != ksup){
//              //position and insert newi in list and
//              // mark it with kcol
//              Idx nexti = head;
//              Idx i;
//              do{
//                i = nexti;
//                nexti = rchlnk[i];
//              }while(newi > nexti);
//              ++knz;
//              rchlnk[i] = newi;
//              rchlnk[newi] = nexti;
//              marker[newi-1] = ksup;
//            }
//          }
//        }
//
//      } 
//
//      //if ksup has no children, insert fstcol into the linked list.
//      if(rchlnk[head] != fstcol){
//        rchlnk[fstcol] = rchlnk[head];
//        rchlnk[head] = fstcol;
//        ++knz;
//      }
//
//      assert(knz == cc[fstcol-1]);
//
//
//      //copy indices from linked list into lindx(*).
//      nzbeg = nzend+1;
//      nzend += knz;
//      assert(nzend+1 == xlindx[ksup]);
//      Idx i = head;
//      for(Ptr kptr = nzbeg; kptr<=nzend;++kptr){
//        i = rchlnk[i];
//        lindx[kptr-1] = i;
//      } 
//
//      //if ksup has a parent, insert ksup into its parent's 
//      //"merge" list.
//      if(length > width){
//        Idx pcol = lindx[xlindx[ksup-1] + width -1];
//        Int psup = SupMembership[pcol-1];
//        mrglnk[ksup-1] = mrglnk[psup-1];
//        mrglnk[psup-1] = ksup;
//      }
//    }
//
//    lindx.resize(nzend+1);
//
//    TIMER_STOP(SymbolicFactorization);
//  }


  void SparseMatrixStructure::SymbolicFactorizationRelaxedDist(ETree& tree,Ordering & aOrder, const SYMPACK::vector<Int> & cc,const SYMPACK::vector<Int> & xsuper,const SYMPACK::vector<Int> & SupMembership, PtrVec & xlindx, IdxVec & lindx,MPI_Comm & comm){
    TIMER_START(SymbolicFactorization);


    if(!bIsGlobal){
      throw std::logic_error( "SparseMatrixStructure must be global in order to call SymbolicFactorization\n" );
    }

    Int nsuper = xsuper.size()-1;

    std::list<MPI_Request> mpirequests;


    Ptr nzbeg = 0;
    //nzend points to the last used slot in lindx
    Ptr nzend = 0;

    //tail is the end of list indicator (in rchlnk, not mrglnk)
    Idx tail = size +1;
    Idx head = 0;

    //Array of length nsuper containing the children of 
    //each supernode as a linked list
    SYMPACK::vector<Idx> mrglnk(nsuper,0);

    //Array of length n+1 containing the current linked list 
    //of merged indices (the "reach" set)
    SYMPACK::vector<Idx> rchlnk(size+1);

    //Array of length n used to mark indices as they are introduced
    // into each supernode's index set
    SYMPACK::vector<Int> marker(size,0);

    Int nsuperLocal = (iam!=np-1)?nsuper/np:nsuper-iam*(nsuper/np);
    Int firstSnode = iam*(nsuper/np)+1;
    Int lastSnode = firstSnode + nsuperLocal-1;
    xlindx.resize(nsuperLocal+1);

    //Compute the sum of the column count and resize lindx accordingly
    //nofsub will be the local nnz now
    Ptr nofsub = 1;
    for(Int ksup = 1; ksup<=nsuper; ++ksup){
      if(ksup>=firstSnode && ksup<=lastSnode){
        Int fstcol = xsuper[ksup-1];
        nofsub+=cc[fstcol-1];
      }
    }
    lindx.resize(nofsub);

    Ptr point = 1;
    for(Int ksup = 1; ksup<=nsuper; ++ksup){
      if(ksup>=firstSnode && ksup<=lastSnode){
        Int locSnode = ksup - firstSnode +1;
        Int fstcol = xsuper[ksup-1];
        xlindx[locSnode-1] = point;
        point += cc[fstcol-1]; 
      }
    } 
    xlindx[nsuperLocal] = point;
    //logfileptr->OFS()<<xlindx<<endl;

    SYMPACK::vector<Ptr> recvXlindx;
    SYMPACK::vector<Idx> recvLindx;

    if(iam>0){
      //build the mrglnk array


      for(Int ksup = 0; ksup<firstSnode; ++ksup){
        Int fstcol = xsuper[ksup-1];
        Int lstcol = xsuper[ksup]-1;
        Int width = lstcol - fstcol +1;
        Int length = cc[fstcol-1];

        //if ksup has a parent, insert ksup into its parent's 
        //"merge" list.
        if(length > width){
        Idx pcol = tree.PostParent(fstcol+width-1-1);  
          Int psup = SupMembership[pcol-1];
          mrglnk[ksup-1] = mrglnk[psup-1];
          mrglnk[psup-1] = ksup;
        }
      }





      //recv from iam-1
//      MPI_Recv(&rchlnk[0],rchlnk.size()*sizeof(Idx),MPI_BYTE,iam-1,(iam-1),comm,MPI_STATUS_IGNORE);
//      MPI_Recv(&mrglnk[0],mrglnk.size()*sizeof(Idx),MPI_BYTE,iam-1,(iam-1),comm,MPI_STATUS_IGNORE);
    }

//logfileptr->OFS()<<"P"<<iam<<" is active"<<endl;
//logfileptr->OFS()<<mrglnk<<endl;
//logfileptr->OFS()<<rchlnk<<endl;

    for(Int locksup = 1; locksup<=nsuperLocal; ++locksup){
      Int ksup = locksup + firstSnode - 1;
      Int fstcol = xsuper[ksup-1];
      Int lstcol = xsuper[ksup]-1;
      Int width = lstcol - fstcol +1;
      Int length = cc[fstcol-1];
      Ptr knz = 0;
      rchlnk[head] = tail;

      //If ksup has children in the supernodal e-tree
      Int jsup = mrglnk[ksup-1];
      if(jsup>0){
        //copy the indices of the first child jsup into 
        //the linked list, and mark each with the value 
        //ksup.
        Int parentJ = ksup;


        do{

//logfileptr->OFS()<<jsup<<endl;
          Int jwidth = xsuper[jsup]-xsuper[jsup-1];

          Ptr * jxlindx = NULL;
          Idx * jlindx = NULL;
          Int locjsup = -1;
          if(jsup>=firstSnode && jsup<=lastSnode){
            locjsup = jsup - firstSnode +1;
            jxlindx = &xlindx[0];
            jlindx = &lindx[0];
          }
          else{
            MPI_Status status;
            recvLindx.resize(size);
            //receive jsup lindx
            Int psrc = min( (Int)np-1, (jsup-1) / (nsuper/np) );
            //logfileptr->OFS()<<"trying to recv "<<jsup<<" max "<<size*sizeof(Idx)<<" bytes"<<" from P"<<psrc<<endl;
            MPI_Recv(&recvLindx[0],recvLindx.size()*sizeof(Idx),MPI_BYTE,psrc,jsup+np,comm,&status);
            //get actual number of received elements
            int count = 0;
            MPI_Get_count(&status,MPI_BYTE,&count);
            count/=sizeof(Idx);

            //compute jsup xlindx
            recvXlindx.resize(2);
            recvXlindx[0] = 1;
            recvXlindx[1] = count +1; 

            locjsup = 1;
            jxlindx = &recvXlindx[0];
            jlindx = &recvLindx[0];
          }

          Ptr jnzbeg = jxlindx[locjsup-1] + jwidth;
          Ptr jnzend = jxlindx[locjsup] -1;
          if(parentJ == ksup){
            for(Ptr jptr = jnzend; jptr>=jnzbeg; --jptr){
              Idx newi = jlindx[jptr-1];
              ++knz;
              marker[newi-1] = ksup;
              rchlnk[newi] = rchlnk[head];
              rchlnk[head] = newi;
            }
          }
          else{
            Int nexti = head;
            for(Ptr jptr = jnzbeg; jptr<=jnzend; ++jptr){
              Idx newi = jlindx[jptr-1];
              Idx i;
              do{
                i = nexti;
                nexti = rchlnk[i];
              }while(newi > nexti);

              if(newi < nexti){
                ++knz;
                rchlnk[i] = newi;
                rchlnk[newi] = nexti;
                marker[newi-1] = ksup;
                nexti = newi;
              }
            }
          }

          parentJ = jsup;
          jsup = mrglnk[jsup-1];
        } while(jsup!=0 && knz < length);
      }

      //structure of a(*,fstcol) has not been examined yet.  
      //"sort" its structure into the linked list,
      //inserting only those indices not already in the
      //list.
      if(knz < length){

        //loop on local columns instead for LOCAL EXPANDED structure

        for(Int row = fstcol; row<=lstcol; ++row){
          Idx newi = row;
          if(newi > fstcol && marker[newi-1] != ksup){
            //position and insert newi in list and
            // mark it with kcol
            Idx nexti = head;
            Idx i;
            do{
              i = nexti;
              nexti = rchlnk[i];
            }while(newi > nexti);
            ++knz;
            rchlnk[i] = newi;
            rchlnk[newi] = nexti;
            marker[newi-1] = ksup;
          }
        }


        for(Int col = fstcol; col<=lstcol; ++col){
          Int node = aOrder.perm[col-1];

          Ptr knzbeg = expColptr[node-1];
          Ptr knzend = expColptr[node]-1;
          for(Ptr kptr = knzbeg; kptr<=knzend;++kptr){
            Idx newi = expRowind[kptr-1];
            newi = aOrder.invp[newi-1];

            if(newi > fstcol && marker[newi-1] != ksup){
              //position and insert newi in list and
              // mark it with kcol
              Idx nexti = head;
              Idx i;
              do{
                i = nexti;
                nexti = rchlnk[i];
              }while(newi > nexti);
              ++knz;
              rchlnk[i] = newi;
              rchlnk[newi] = nexti;
              marker[newi-1] = ksup;
            }
          }
        }

      } 

      //if ksup has no children, insert fstcol into the linked list.
      if(rchlnk[head] != fstcol){
        rchlnk[fstcol] = rchlnk[head];
        rchlnk[head] = fstcol;
        ++knz;
      }

      assert(knz == cc[fstcol-1]);


      //copy indices from linked list into lindx(*).
      nzbeg = nzend+1;
      nzend += knz;
      assert(nzend+1 == xlindx[locksup]);
      Idx i = head;
      for(Ptr kptr = nzbeg; kptr<=nzend;++kptr){
        i = rchlnk[i];
        lindx[kptr-1] = i;
      } 

      //if ksup has a parent, insert ksup into its parent's 
      //"merge" list.
      if(length > width){
        Idx pcol = tree.PostParent(fstcol+width-1-1);  
//        Idx pcol = lindx[xlindx[locksup-1] + width -1];
//logfileptr->OFS()<<pcol2<<" vs "<<pcol<<endl;
//assert(pcol==pcol2);
        Int psup = SupMembership[pcol-1];
        mrglnk[ksup-1] = mrglnk[psup-1];
        mrglnk[psup-1] = ksup;


        //send L asap
        Int pdest = min( (Int)np-1, (psup-1) / (nsuper/np) );
        //if remote
        if(pdest!=iam){
                mpirequests.push_back(MPI_REQUEST_NULL);
                MPI_Request & request = mpirequests.back();
                //logfileptr->OFS()<<"sending "<<lengthj*sizeof(Idx)<<" bytes of "<<jsup<<" to P"<<pdest<<" for "<<ksup<<endl;
                MPI_Isend(&lindx[xlindx[locksup-1]-1],length*sizeof(Idx),MPI_BYTE,pdest,ksup+np,comm,&request);
        }
      }

    }
    if(np>1 && iam<np-1){

      //send to iam+1
      //marker, rchlnk, mrglnk
      //MPI_Send(&rchlnk[0],rchlnk.size()*sizeof(Idx),MPI_BYTE,iam+1,iam,comm);
      //MPI_Send(&mrglnk[0],mrglnk.size()*sizeof(Idx),MPI_BYTE,iam+1,iam,comm);


#if 0
      Int nextFirstSnode = (iam+1)*(nsuper/np)+1;
      for(Int ksup = nextFirstSnode; ksup<=nsuper; ++ksup){
        Int fstcol = xsuper[ksup-1];
        Int lstcol = xsuper[ksup]-1;
        Int width = lstcol - fstcol +1;
        Int length = cc[fstcol-1];

        //for all children

        //If ksup has children in the supernodal e-tree
        Int jsup = mrglnk[ksup-1];
        if(jsup>0){
          do{
            if(jsup>=firstSnode && jsup<=lastSnode){
              Int locjsup = jsup - firstSnode + 1;
              Int pdest = min( (Int)np-1, (ksup-1) / (nsuper/np) );
              //if remote
              if(pdest!=iam){
                Int lengthj = xlindx[locjsup] - xlindx[locjsup-1];
                mpirequests.push_back(MPI_REQUEST_NULL);
                MPI_Request & request = mpirequests.back();
                //logfileptr->OFS()<<"sending "<<lengthj*sizeof(Idx)<<" bytes of "<<jsup<<" to P"<<pdest<<" for "<<ksup<<endl;
                MPI_Isend(&lindx[xlindx[locjsup-1]-1],lengthj*sizeof(Idx),MPI_BYTE,pdest,jsup+np,comm,&request);
             }
            }
            //for each subsequent child jsup of ksup ...
            jsup = mrglnk[jsup-1];
          }while(jsup>0);
        }
      }
#endif




      //wait all on every requests
      //for(auto it = mpirequests.begin();it!=mpirequests.end();it++){
      //  MPI_Request & request = *it;
      //  MPI_Wait(&request,MPI_STATUS_IGNORE);
      //}
    }


    lindx.resize(nzend+1);
    MPI_Barrier(comm);

    TIMER_STOP(SymbolicFactorization);
  }




  void SparseMatrixStructure::SymbolicFactorizationRelaxed(ETree& tree,Ordering & aOrder, const SYMPACK::vector<Int> & cc,const SYMPACK::vector<Int> & xsuper,const SYMPACK::vector<Int> & SupMembership, PtrVec & xlindx, IdxVec & lindx){
    TIMER_START(SymbolicFactorization);


    if(!bIsGlobal){
      throw std::logic_error( "SparseMatrixStructure must be global in order to call SymbolicFactorization\n" );
    }

    Int nsuper = xsuper.size()-1;




    Ptr nzbeg = 0;
    //nzend points to the last used slot in lindx
    Ptr nzend = 0;

    //tail is the end of list indicator (in rchlnk, not mrglnk)
    Idx tail = size +1;

    Idx head = 0;

    //Array of length nsuper containing the children of 
    //each supernode as a linked list
    SYMPACK::vector<Idx> mrglnk(nsuper,0);

    //Array of length n+1 containing the current linked list 
    //of merged indices (the "reach" set)
    SYMPACK::vector<Idx> rchlnk(size+1);

    //Array of length n used to mark indices as they are introduced
    // into each supernode's index set
    SYMPACK::vector<Int> marker(size,0);



    xlindx.resize(nsuper+1);

    //Compute the sum of the column count and resize lindx accordingly
    Ptr nofsub = 1;
    for(Int ksup = 1; ksup<=nsuper; ++ksup){
      Int fstcol = xsuper[ksup-1];
      nofsub+=cc[fstcol-1];
    }

//    for(Int i =0; i<cc.size();++i){
//      nofsub+=cc[i];
//    }

    lindx.resize(nofsub);


    Ptr point = 1;
    for(Int ksup = 1; ksup<=nsuper; ++ksup){
      Int fstcol = xsuper[ksup-1];
      xlindx[ksup-1] = point;
      point += cc[fstcol-1]; 
    } 
    xlindx[nsuper] = point;
    //logfileptr->OFS()<<xlindx<<endl;



    for(Int ksup = 1; ksup<=nsuper; ++ksup){
      Int fstcol = xsuper[ksup-1];
      Int lstcol = xsuper[ksup]-1;
      Int width = lstcol - fstcol +1;
      Int length = cc[fstcol-1];
      Ptr knz = 0;
      rchlnk[head] = tail;
      Int jsup = mrglnk[ksup-1];


      //If ksup has children in the supernodal e-tree
      if(jsup>0){
        //copy the indices of the first child jsup into 
        //the linked list, and mark each with the value 
        //ksup.
        Int jwidth = xsuper[jsup]-xsuper[jsup-1];
        Ptr jnzbeg = xlindx[jsup-1] + jwidth;
        Ptr jnzend = xlindx[jsup] -1;
        for(Ptr jptr = jnzend; jptr>=jnzbeg; --jptr){
          Idx newi = lindx[jptr-1];
          ++knz;
          marker[newi-1] = ksup;
          rchlnk[newi] = rchlnk[head];
          rchlnk[head] = newi;
        }

        //for each subsequent child jsup of ksup ...
        jsup = mrglnk[jsup-1];
        while(jsup!=0 && knz < length){
          //merge the indices of jsup into the list,
          //and mark new indices with value ksup.

          jwidth = xsuper[jsup]-xsuper[jsup-1];
          jnzbeg = xlindx[jsup-1] + jwidth;
          jnzend = xlindx[jsup] -1;
          Int nexti = head;
          for(Ptr jptr = jnzbeg; jptr<=jnzend; ++jptr){
            Idx newi = lindx[jptr-1];
            Idx i;
            do{
              i = nexti;
              nexti = rchlnk[i];
            }while(newi > nexti);

            if(newi < nexti){
#ifdef _DEBUG_
              logfileptr->OFS()<<jsup<<" is a child of "<<ksup<<" and "<<newi<<" is inserted in the structure of "<<ksup<<std::endl;
#endif
              ++knz;
              rchlnk[i] = newi;
              rchlnk[newi] = nexti;
              marker[newi-1] = ksup;
              nexti = newi;
            }
          }
          jsup = mrglnk[jsup-1];
        }
      }

      //structure of a(*,fstcol) has not been examined yet.  
      //"sort" its structure into the linked list,
      //inserting only those indices not already in the
      //list.
      if(knz < length){
        for(Int row = fstcol; row<=lstcol; ++row){
          Idx newi = row;
          if(newi > fstcol && marker[newi-1] != ksup){
            //position and insert newi in list and
            // mark it with kcol
            Idx nexti = head;
            Idx i;
            do{
              i = nexti;
              nexti = rchlnk[i];
            }while(newi > nexti);
            ++knz;
            rchlnk[i] = newi;
            rchlnk[newi] = nexti;
            marker[newi-1] = ksup;
          }
        }


        for(Int col = fstcol; col<=lstcol; ++col){
          Int node = aOrder.perm[col-1];

          Ptr knzbeg = expColptr[node-1];
          Ptr knzend = expColptr[node]-1;
          for(Ptr kptr = knzbeg; kptr<=knzend;++kptr){
            Idx newi = expRowind[kptr-1];
            newi = aOrder.invp[newi-1];

            if(newi > fstcol && marker[newi-1] != ksup){
              //position and insert newi in list and
              // mark it with kcol
              Idx nexti = head;
              Idx i;
              do{
                i = nexti;
                nexti = rchlnk[i];
              }while(newi > nexti);
              ++knz;
              rchlnk[i] = newi;
              rchlnk[newi] = nexti;
              marker[newi-1] = ksup;
            }
          }
        }

      } 

      //if ksup has no children, insert fstcol into the linked list.
      if(rchlnk[head] != fstcol){
        rchlnk[fstcol] = rchlnk[head];
        rchlnk[head] = fstcol;
        ++knz;
      }

      assert(knz == cc[fstcol-1]);


      //copy indices from linked list into lindx(*).
      nzbeg = nzend+1;
      nzend += knz;
      assert(nzend+1 == xlindx[ksup]);
      Idx i = head;
      for(Ptr kptr = nzbeg; kptr<=nzend;++kptr){
        i = rchlnk[i];
        lindx[kptr-1] = i;
      } 

      //if ksup has a parent, insert ksup into its parent's 
      //"merge" list.
      if(length > width){
        Idx pcol = lindx[xlindx[ksup-1] + width -1];
        Int psup = SupMembership[pcol-1];
        mrglnk[ksup-1] = mrglnk[psup-1];
        mrglnk[psup-1] = ksup;
      }
    }

    lindx.resize(nzend+1);

    TIMER_STOP(SymbolicFactorization);
  }

  void SparseMatrixStructure::SymbolicFactorization(ETree& tree, Ordering & aOrder, const SYMPACK::vector<Int> & cc,const SYMPACK::vector<Int> & xsuper,const SYMPACK::vector<Int> & SupMembership, PtrVec & xlindx, IdxVec & lindx){
    TIMER_START(SymbolicFactorization);


    if(!bIsGlobal){
      throw std::logic_error( "SparseMatrixStructure must be global in order to call SymbolicFactorization\n" );
    }

    Int nsuper = xsuper.size()-1;

    Ptr nzbeg = 0;
    //nzend points to the last used slot in lindx
    Ptr nzend = 0;

    //tail is the end of list indicator (in rchlnk, not mrglnk)
    Idx tail = size +1;

    Idx head = 0;

    //Array of length nsuper containing the children of 
    //each supernode as a linked list
    SYMPACK::vector<Idx> mrglnk(nsuper,0);

    //Array of length n+1 containing the current linked list 
    //of merged indices (the "reach" set)
    SYMPACK::vector<Idx> rchlnk(size+1);

    //Array of length n used to mark indices as they are introduced
    // into each supernode's index set
    SYMPACK::vector<Int> marker(size,0);



    xlindx.resize(nsuper+1);

    //Compute the sum of the column count and resize lindx accordingly
    //Ptr nofsub = 1;
    //for(Int i =0; i<cc.size();++i){
    //  nofsub+=cc[i];
    //}

    Ptr nofsub = 1;
    for(Int ksup = 1; ksup<=nsuper; ++ksup){
      Int fstcol = xsuper[ksup-1];
      nofsub+=cc[fstcol-1];
    }
    lindx.resize(nofsub);


    Ptr point = 1;
    for(Int ksup = 1; ksup<=nsuper; ++ksup){
      Int fstcol = xsuper[ksup-1];
      xlindx[ksup-1] = point;
      point += cc[fstcol-1]; 
    } 
    xlindx[nsuper] = point;

    for(Int ksup = 1; ksup<=nsuper; ++ksup){
      Int fstcol = xsuper[ksup-1];
      Int lstcol = xsuper[ksup]-1;
      Int width = lstcol - fstcol +1;
      Int length = cc[fstcol-1];
      Idx knz = 0;
      rchlnk[head] = tail;
      Int jsup = mrglnk[ksup-1];

      //If ksup has children in the supernodal e-tree
      if(jsup>0){
        //copy the indices of the first child jsup into 
        //the linked list, and mark each with the value 
        //ksup.
        Int jwidth = xsuper[jsup]-xsuper[jsup-1];
        Ptr jnzbeg = xlindx[jsup-1] + jwidth;
        Ptr jnzend = xlindx[jsup] -1;
        for(Ptr jptr = jnzend; jptr>=jnzbeg; --jptr){
          Idx newi = lindx[jptr-1];
          ++knz;
          marker[newi-1] = ksup;
          rchlnk[newi] = rchlnk[head];
          rchlnk[head] = newi;
        }

        //for each subsequent child jsup of ksup ...
        jsup = mrglnk[jsup-1];
        while(jsup!=0 && knz < length){
          //merge the indices of jsup into the list,
          //and mark new indices with value ksup.

          jwidth = xsuper[jsup]-xsuper[jsup-1];
          jnzbeg = xlindx[jsup-1] + jwidth;
          jnzend = xlindx[jsup] -1;
          Idx nexti = head;
          for(Ptr jptr = jnzbeg; jptr<=jnzend; ++jptr){
            Idx newi = lindx[jptr-1];
            Idx i;
            do{
              i = nexti;
              nexti = rchlnk[i];
            }while(newi > nexti);

            if(newi < nexti){
#ifdef _DEBUG_
              logfileptr->OFS()<<jsup<<" is a child of "<<ksup<<" and "<<newi<<" is inserted in the structure of "<<ksup<<std::endl;
#endif
              ++knz;
              rchlnk[i] = newi;
              rchlnk[newi] = nexti;
              marker[newi-1] = ksup;
              nexti = newi;
            }
          }
          jsup = mrglnk[jsup-1];
        }
      }

      //structure of a(*,fstcol) has not been examined yet.  
      //"sort" its structure into the linked list,
      //inserting only those indices not already in the
      //list.
      if(knz < length){
        Idx node = aOrder.perm[fstcol-1];
        Ptr knzbeg = expColptr[node-1];
        Ptr knzend = expColptr[node]-1;
        for(Ptr kptr = knzbeg; kptr<=knzend;++kptr){
          Idx newi = expRowind[kptr-1];
          newi = aOrder.invp[newi-1];
          if(newi > fstcol && marker[newi-1] != ksup){
            //position and insert newi in list and
            // mark it with kcol
            Idx nexti = head;
            Idx i;
            do{
              i = nexti;
              nexti = rchlnk[i];
            }while(newi > nexti);
            ++knz;
            rchlnk[i] = newi;
            rchlnk[newi] = nexti;
            marker[newi-1] = ksup;
          }
        }
      }

      //if ksup has no children, insert fstcol into the linked list.
      if(rchlnk[head] != fstcol){
        rchlnk[fstcol] = rchlnk[head];
        rchlnk[head] = fstcol;
        ++knz;
      }

      assert(knz == cc[fstcol-1]);


      //copy indices from linked list into lindx(*).
      nzbeg = nzend+1;
      nzend += knz;
      assert(nzend+1 == xlindx[ksup]);
      Int i = head;
      for(Int kptr = nzbeg; kptr<=nzend;++kptr){
        i = rchlnk[i];
        lindx[kptr-1] = i;
      } 

      //if ksup has a parent, insert ksup into its parent's 
      //"merge" list.
      if(length > width){
        Int pcol = lindx[xlindx[ksup-1] + width -1];
        Int psup = SupMembership[pcol-1];
        mrglnk[ksup-1] = mrglnk[psup-1];
        mrglnk[psup-1] = ksup;
      }
    }

    lindx.resize(nzend+1);

    TIMER_STOP(SymbolicFactorization);
  }







  //////DEPRECATED METHODS
  //////NOT WORKING
  ////  void SparseMatrixStructure::SymbolicFactorizationDEPRECATED(ETree& tree,const SYMPACK::vector<Int> & cc,const SYMPACK::vector<Int> & xsuper, SYMPACK::vector<Int> & xlindx, SYMPACK::vector<Int> & lindx){
  ////    TIMER_START(SymbolicFactorization);
  ////
  ////
  ////    if(!bIsGlobal){
  ////			throw std::logic_error( "SparseMatrixStructure must be global in order to call SymbolicFactorization\n" );
  ////    }
  ////
  ////    SYMPACK::vector<std::set<Int> > sets;
  ////    sets.resize(xsuper.size(),std::set<Int>());
  ////
  ////    SYMPACK::vector<SYMPACK::vector<Int> > LIs;
  ////    LIs.resize(xsuper.size());
  ////
  ////    Int lindxCnt = 0;
  ////    for(Int I=1;I<xsuper.size();I++){
  ////      Int fi = tree.FromPostOrder(xsuper(I-1));
  ////      Int width = xsuper(I)-xsuper(I-1);
  ////      Int length = cc(fi-1);
  ////
  ////      SYMPACK::vector<Int> & LI = LIs[I-1];
  ////
  ////
  ////      //Initialize LI with nnz struct of A_*fi
  ////      Int begin = colptr(fi-1);
  ////      Int end = colptr(fi);
  ////
  ////      Int * start = &rowind(begin-1); 
  ////      Int * stop = (end>=rowind.size())?&rowind(rowind.size()-1)+1:&rowind(end-1); 
  ////      //find the diagonal block
  ////      start=std::find(start,stop,fi);
  ////
  ////
  ////      LI.Resize(stop-start);
  ////
  ////      std::copy(start,stop,&LI(0));
  ////
  ////      //     logfileptr->OFS()<<"L"<<I<<"<- A_*,fi: ";
  ////      //     for(int i=0;i<LI.size();i++){logfileptr->OFS()<<LI(i)<< " ";}
  ////      //     logfileptr->OFS()<<std::endl;
  ////
  ////      LI = tree.ToPostOrder(LI);
  ////
  ////
  ////      //     logfileptr->OFS()<<"PO L"<<I<<"<- A_*,fi: ";
  ////      //     for(int i=0;i<LI.size();i++){logfileptr->OFS()<<LI(i)<< " ";}
  ////      //     logfileptr->OFS()<<std::endl;
  ////
  ////
  ////
  ////
  ////      std::set<Int> & SI = sets[I-1];
  ////      for(std::set<Int>::iterator it = SI.begin(); it!=SI.end(); it++){
  ////        Int K = *it;
  ////        SYMPACK::vector<Int> & LK = LIs[K-1];
  ////
  ////        //        logfileptr->OFS()<<"merging "<<I<<" with "<<K<<std::endl;
  ////        //LI = LI U LK \ K
  ////        SYMPACK::vector<Int> Ltmp(LI.size()+LK.size()-1);
  ////        std::copy(&LI(0),&LI(LI.size()-1)+1,&Ltmp(0));
  ////
  ////
  ////        if(LK.size()>1){
  ////
  ////          //Be careful to not insert duplicates !
  ////          Int firstidx =1;
  ////          for(Int i =1;i<LK.size();i++){
  ////            if(LK[i]>fi ){
  ////              firstidx = i+1;
  ////              break;
  ////            }
  ////          }
  ////          Int * end = std::set_union(&LI(0),&LI(LI.size()-1)+1,&LK(firstidx-1),&LK(LK.size()-1)+1,&Ltmp(0));
  ////          Ltmp.Resize(end - &Ltmp(0));
  ////          LI = Ltmp;
  ////        }
  ////
  ////
  ////      }
  ////
  ////
  ////      lindxCnt += LI.size();
  ////#ifdef _DEBUG_
  ////        logfileptr->OFS()<<"I:"<<I<<std::endl<<"LI:"<<LI<<std::endl;
  ////#endif
  ////
  ////      if(length>width  ){
  ////#ifdef _DEBUG_
  ////        logfileptr->OFS()<<"I:"<<I<<std::endl<<"LI:"<<LI<<std::endl;
  ////#endif
  ////        Int i = LI(width);
  ////        //        logfileptr->OFS()<<I<<" : col "<<i<<" is the next to be examined. width="<<width<<" length="<<length<<std::endl;
  ////
  ////        Int J = I+1;
  ////        for(J = I+1;J<=xsuper.size()-1;J++){
  ////          Int fc = xsuper(J-1);
  ////          Int lc = xsuper(J)-1;
  ////          //          logfileptr->OFS()<<"FC = "<<fc<<" vs "<<i<<std::endl;
  ////          //          logfileptr->OFS()<<"LC = "<<lc<<" vs "<<i<<std::endl;
  ////          if(fc <=i && lc >= i){
  ////            //            logfileptr->OFS()<<I<<" : col "<<i<<" found in snode "<<J<<std::endl;
  ////            break;
  ////          }
  ////        } 
  ////
  ////        //        logfileptr->OFS()<<I<<" : col "<<i<<" is in snode "<<J<<std::endl;
  ////        std::set<Int> & SJ = sets[J-1];
  ////        SJ.insert(I);
  ////        //        logfileptr->OFS()<<"S"<<J<<" U {"<<I<<"}"<<std::endl; 
  ////
  ////      }
  ////
  ////    }  
  ////
  ////    Int nsuper = xsuper.size()-1;
  ////    Int totNnz = 1;
  ////    for(Int I=1;I<xsuper.size();I++){
  ////      Int fc = xsuper(I-1);
  ////      Int lc = xsuper(I)-1;
  ////      for(Int i=fc;i<=lc;i++){
  ////        totNnz+=cc(i-1);
  ////      }
  ////
  ////    }
  ////
  ////    lindx.Resize(lindxCnt);
  ////    xlindx.Resize(nsuper+1);
  ////    Int head = 1;
  ////
  ////    //    logfileptr->OFS()<<"There are "<<lindxCnt<<" slots in lindx"<<std::endl;
  ////    for(Int I=1;I<=nsuper;I++){
  ////      Int fi = tree.FromPostOrder(xsuper(I-1));
  ////      SYMPACK::vector<Int> & LI = LIs[I-1];
  ////      xlindx(I-1)=head;
  ////
  ////
  ////      //        logfileptr->OFS()<<"PO L"<<I<<":";
  ////      //        for(int i=0;i<LI.size();i++){logfileptr->OFS()<<LI(i)<< " ";}
  ////      //        logfileptr->OFS()<<std::endl;
  ////
  ////      //      logfileptr->OFS()<<"Copying "<<LI.size()<<" elem into lindx("<<head-1<<")"<<std::endl;
  ////      for(Int i=0;i<LI.size();++i){
  ////        if(LI(i)!=0){
  ////          lindx(head-1) = LI(i); 
  ////          head++;
  ////        }
  ////      }
  ////      //std::copy(&LI(0),LI.Data()+LI.size(),&(lindx(head-1)));
  ////      //head+=LI.size();//cc(fi-1);
  ////    }
  ////    //    lindx.Resize(head-1);
  ////    xlindx(nsuper) = head;
  ////
  ////    TIMER_STOP(SymbolicFactorization);
  ////  }
  ////
  //////FIXME correct these methods
  //////Return the row structure in the permuted matrix
  ////void SparseMatrixStructure::GetARowStruct(const ETree & etree, const Int iPORow, SYMPACK::vector<Int> & rowStruct){
  ////  TIMER_START(SparseMatrixStructure::GetARowStruct);
  ////    if(!bIsGlobal){
  ////			throw std::logic_error( "SparseMatrixStructure must be global in order to call GetARowStruct\n" );
  ////    }
  ////
  ////  for(Int iPOCurCol = 1; iPOCurCol<iPORow;++iPOCurCol){
  ////    Int iCurCol = etree.FromPostOrder(iPOCurCol);
  ////
  ////    Int iFirstRowPtr = colptr(iCurCol-1);
  ////    Int iLastRowPtr = colptr(iCurCol)-1;
  ////
  //////    logfileptr->OFS()<<iFirstRowPtr<<" to "<<iLastRowPtr<<std::endl;
  ////    for(Int iCurRowPtr=iFirstRowPtr ;iCurRowPtr<=iLastRowPtr;++iCurRowPtr){
  ////      Int iPOCurRow = etree.ToPostOrder(rowind(iCurRowPtr-1));
  ////      if(iPOCurRow == iPORow){
  ////        rowStruct.push_back(iPOCurCol);
  //////        logfileptr->OFS()<<"A("<<iPORow<<","<<iPOCurCol<<") is nz "<<std::endl;
  ////      }
  ////
  ////      if(iPOCurRow >= iPORow){
  ////        break;
  ////      }
  ////    }
  ////  }
  ////  TIMER_STOP(SparseMatrixStructure::GetARowStruct);
  ////}
  ////
  ////void SparseMatrixStructure::GetLRowStruct(const ETree & etree, const Int iPORow, const SYMPACK::vector<Int> & ARowStruct, std::set<Int> & LRowStruct){
  ////
  ////  TIMER_START(SparseMatrixStructure::GetLRowStruct);
  ////    if(!bIsGlobal){
  ////			throw std::logic_error( "SparseMatrixStructure must be global in order to call GetLRowStruct\n" );
  ////    }
  //////  LRowStruct.clear();
  ////  for(Int i = 0; i<ARowStruct.size();++i){
  ////    Int iCurNode = ARowStruct[i];
  ////    //tracing from iCurRow to iRow;
  //////    logfileptr->OFS()<<"Starting from node "<<iCurNode<<std::endl;
  ////    if(iCurNode==iPORow){
  //////      logfileptr->OFS()<<"Adding node "<<iCurNode<<std::endl;
  ////      LRowStruct.insert(iCurNode);
  ////    }
  ////    else{
  ////      while(iCurNode != iPORow && etree.PostParent(iCurNode-1) != 0){
  //////        logfileptr->OFS()<<"Adding node "<<iCurNode<<std::endl;
  ////        LRowStruct.insert(iCurNode);
  ////        iCurNode = etree.PostParent(iCurNode-1);
  //////        logfileptr->OFS()<<"Parent is "<<iCurNode<<std::endl;
  ////      }
  ////    }
  ////  } 
  ////  TIMER_STOP(SparseMatrixStructure::GetLRowStruct);
  ////}
  ////
  ////void SparseMatrixStructure::GetSuperARowStruct(const ETree & etree, const SYMPACK::vector<Int> & Xsuper, const SYMPACK::vector<Int> & SupMembership, const Int iSupNo, SYMPACK::vector<Int> & SuperRowStruct){
  ////  TIMER_START(SpStruct_GetSuperARowStruct);
  //////  TIMER_START(SparseMatrixStructure::GetSuperARowStruct);
  ////    if(!bIsGlobal){
  ////			throw std::logic_error( "SparseMatrixStructure must be global in order to call GetARowStruct\n" );
  ////    }
  ////
  ////  Int first_col = Xsuper[iSupNo-1];
  ////  Int last_col = Xsuper[iSupNo]-1;
  ////
  ////  for(Int iPOCurCol = 1; iPOCurCol<first_col;++iPOCurCol){
  ////    Int iCurCol = etree.FromPostOrder(iPOCurCol);
  ////    Int iCurSupno = SupMembership(iPOCurCol-1);
  ////    
  ////
  ////    Int iFirstRowPtr = colptr(iCurCol-1);
  ////    Int iLastRowPtr = colptr(iCurCol)-1;
  ////
  //////    logfileptr->OFS()<<iFirstRowPtr<<" to "<<iLastRowPtr<<std::endl;
  ////    for(Int iCurRowPtr=iFirstRowPtr ;iCurRowPtr<=iLastRowPtr;++iCurRowPtr){
  ////      Int iPOCurRow = etree.ToPostOrder(rowind(iCurRowPtr-1));
  ////      if(iPOCurRow >= first_col && iPOCurRow <= last_col){
  //////        SuperRowStruct.push_back(iPOCurCol);
  ////        SuperRowStruct.push_back(iCurSupno);
  //////        logfileptr->OFS()<<"A("<<iPORow<<","<<iPOCurCol<<") is nz "<<std::endl;
  ////      }
  ////
  ////      if(iPOCurRow > last_col){
  ////        break;
  ////      }
  ////    }
  ////  }
  //////  TIMER_STOP(SparseMatrixStructure::GetSuperARowStruct);
  ////  TIMER_STOP(SpStruct_GetSuperARowStruct);
  ////}
  ////
  ////void SparseMatrixStructure::GetSuperLRowStruct(const ETree & etree, const SYMPACK::vector<Int> & Xsuper, const SYMPACK::vector<Int> & SupMembership, const Int iSupNo, std::set<Int> & SuperLRowStruct){
  ////
  ////  TIMER_START(SpStruct_GetSuperLRowStruct);
  ////    if(!bIsGlobal){
  ////			throw std::logic_error( "SparseMatrixStructure must be global in order to call GetSuperLRowStruct\n" );
  ////    }
  ////
  //////  SuperLRowStruct.clear();
  ////
  ////    //Get A row struct
  ////    SYMPACK::vector<Int> SuperARowStruct;
  ////    GetSuperARowStruct(etree, Xsuper,SupMembership, iSupNo, SuperARowStruct);
  ////
  ////
  ////
  ////#ifdef _DEBUG_
  ////      logfileptr->OFS()<<"Row structure of A of Supernode "<<iSupNo<<" is ";
  ////      for(SYMPACK::vector<Int>::iterator it = SuperARowStruct.begin(); it != SuperARowStruct.end(); ++it){
  ////        logfileptr->OFS()<<*it<<" ";
  ////      }
  ////      logfileptr->OFS()<<std::endl;
  ////#endif
  ////
  ////  Int first_col = Xsuper[iSupNo-1];
  ////  Int last_col = Xsuper[iSupNo]-1;
  ////
  ////  for(Int i = 0; i<SuperARowStruct.size();++i){
  ////    Int iCurNode = SuperARowStruct[i];
  ////    //tracing from iCurRow to iRow;
  ////#ifdef _DEBUG_
  ////    logfileptr->OFS()<<"Starting from node "<<iCurNode<<std::endl;
  ////#endif
  ////    //if(iCurNode==iPORow){
  ////    if(iCurNode >= first_col && iCurNode <= last_col){
  ////#ifdef _DEBUG_
  ////      logfileptr->OFS()<<"Adding node "<<iCurNode<<std::endl;
  ////#endif
  ////      SuperLRowStruct.insert(iCurNode);
  ////    }
  ////    else{
  ////      while( iCurNode != first_col && etree.PostParent(iCurNode-1) != 0){
  //////      while( !(iCurNode >= first_col && iCurNode <= last_col) && etree.PostParent(iCurNode-1) != 0){
  ////#ifdef _DEBUG_
  ////        logfileptr->OFS()<<"Adding node "<<iCurNode<<std::endl;
  ////#endif
  ////        SuperLRowStruct.insert(iCurNode);
  ////        iCurNode = etree.PostParent(iCurNode-1);
  ////#ifdef _DEBUG_
  ////        logfileptr->OFS()<<"Parent is "<<iCurNode<<std::endl;
  ////#endif
  ////      }
  ////    }
  ////  } 
  ////  TIMER_STOP(SpStruct_GetSuperLRowStruct);
  ////}
  ////
  ////  void SparseMatrixStructure::GetLColRowCountDEPRECATED(ETree & tree, SYMPACK::vector<Int> & cc, SYMPACK::vector<Int> & rc){
  ////
  ////
  ////    if(!bIsGlobal){
  ////			throw std::logic_error( "SparseMatrixStructure must be global in order to call GetLColRowCount\n" );
  ////    }
  ////
  ////    //The tree need to be postordered
  ////    if(!tree.IsPostOrdered()){
  ////      tree.PostOrderTree();
  ////    }
  ////
  ////    ExpandSymmetric();
  ////
  ////    TIMER_START(GetColRowCount);
  //////    TIMER_START(Initialize_Data);
  ////    //cc first contains the delta
  ////    cc.Resize(size);
  ////    //Compute size of subtrees
  ////    SYMPACK::vector<Int> treeSize(size);
  ////    SetValue(treeSize,I_ONE);
  ////
  ////
  ////
  ////    SYMPACK::vector<Int> level(size);
  ////    level(size-1)=1;
  ////    for(Int vertex = 1; vertex<=size-1; vertex++){
  ////      Int curParent = tree.PostParent(vertex-1);
  ////      if(curParent!=0){
  ////        treeSize(curParent-1)+=treeSize(vertex-1);
  ////      }
  ////      if(treeSize(vertex-1)==1){
  ////        cc(vertex-1)=1;
  ////      }
  ////      else{
  ////        cc(vertex-1)= 0;
  ////      }
  ////    }
  ////
  ////    for(Int vertex = size-1; vertex>=1; --vertex){
  ////      Int curParent = tree.PostParent(vertex-1);
  ////      if(curParent!=0){
  ////        level(vertex-1) = level(curParent-1)+1;
  ////      }
  ////      else{
  ////        level(vertex-1) = 0;
  ////      }
  ////    }
  ////
  ////
  ////    if(treeSize(size-1)==1){
  ////      cc(size-1)=1;
  ////    }
  ////    else{
  ////      cc(size-1)= 0 ;
  ////    }
  ////
  ////
  ////
  ////    SYMPACK::vector<Int> prevLeaf(size);
  ////    SetValue(prevLeaf,I_ZERO);
  ////    SYMPACK::vector<Int> prevNz(size);
  ////    SetValue(prevNz,I_ZERO);
  ////
  ////    rc.Resize(size);
  ////    SetValue(rc,I_ONE);
  ////
  ////    DisjointSet sets;
  ////    sets.Initialize(size);
  ////    for(Int vertex = 1; vertex<=size; vertex++){
  ////      Int cset = sets.makeSet (vertex);
  ////      sets.Root(cset-1)=vertex;
  ////    }
  ////
  ////
  //////    TIMER_STOP(Initialize_Data);
  ////
  //////    TIMER_START(Compute_Col_Row_Count);
  ////    for(Int col=1; col<size; col++){
  ////      Int cset;
  ////
  ////      Int colPar = tree.PostParent(col-1);
  ////      if (col<size && colPar!=0){
  ////        cc(colPar-1)--;
  ////      }
  ////
  ////      Int oCol = tree.FromPostOrder(col);
  ////      for (Int i = expColptr(oCol-1); i < expColptr(oCol); i++) {
  ////        Int row = tree.ToPostOrder(expRowind(i-1));
  ////        if (row > col){
  ////          Int k = prevNz(row-1);
  ////
  ////
  ////#ifdef _DEBUG_
  ////          logfileptr->OFS()<<"prevNz("<<row<<")="<<k<<" vs "<< col - treeSize(col-1) +1<<std::endl;
  ////#endif
  ////          if(k< col - treeSize(col-1) +1){
  ////#ifdef _DEBUG_
  ////            logfileptr->OFS()<<"Vertex "<<col<<" is a leaf of Tr["<<row<<"]"<<std::endl;
  ////#endif
  ////            cc(col-1)++;
  ////
  ////            Int p = prevLeaf(row-1);
  ////            if(p==0){
  ////              rc(row-1)+=level(col-1)-level(row-1);
  ////            }
  ////            else {
  //////              TIMER_START(Get_LCA);
  ////              Int pset = sets.find(p);
  ////              Int q = sets.Root(pset-1);
  //////              TIMER_STOP(Get_LCA);
  ////
  ////#ifdef _DEBUG_
  ////              logfileptr->OFS()<<"Vertex "<<q<<" is the LCA of "<<p<<" and "<< col<<std::endl;
  ////#endif
  ////              rc(row-1)+= level(col-1) - level(q-1);
  ////              cc(q-1)--;
  ////
  ////            }
  ////            prevLeaf(row-1)=col;
  ////          }
  ////#ifdef _DEBUG_
  ////          else{
  ////            logfileptr->OFS()<<"Vertex "<<col<<" is an internal vertex of Tr["<<row<<"]"<<std::endl;
  ////          }
  ////#endif
  ////          prevNz(row-1)=col;
  ////        }
  ////      }
  ////
  //////      TIMER_START(Merge_LCA);
  ////      //merge col and parent sets (for lca computation)
  ////      if (colPar!=0){
  ////        sets.Union(col,colPar);
  ////      }
  //////      TIMER_STOP(Merge_LCA);
  ////
  ////    }
  //////    TIMER_STOP(Compute_Col_Row_Count);
  ////
  ////
  ////#ifdef _DEBUG_
  ////        logfileptr->OFS()<<"Deltas "<<cc.size()<<std::endl;
  ////        for(Int i = 0; i<cc.size();i++){
  ////          logfileptr->OFS()<<cc(i)<<" ";
  ////        }
  ////        logfileptr->OFS()<<std::endl;
  ////#endif
  ////
  ////
  ////
  ////
  ////
  ////    //convert delta to col count
  ////    for(Int col=1; col<size; col++){
  ////      Int parent = tree.PostParent(col-1);
  ////      if(parent!=0){
  ////        cc(parent-1)+= cc(col-1);
  ////      }
  ////    }
  ////
  ////
  ////
  ////
  ////    TIMER_STOP(GetColRowCount);
  ////  }



}
