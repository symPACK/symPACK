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

namespace symPACK{


  //special reduction: first element contains the max local sum
  void _PtrSum( void *in, void *inout, int *len, MPI_Datatype *dptr ) 
  { 
    int i; 

    Ptr * pinout = (Ptr*)inout;
    Ptr * pin = (Ptr*)in;
    pinout[0] = std::max(pinout[0],pin[0]);
    //logfileptr->OFS()<<"REDUCE max is "<<pinout[0]<<" vs (max) "<<pin[0]<<std::endl;
#pragma unroll
    for (i=1; i< *len; ++i) { 
      //logfileptr->OFS()<<"REDUCE elem is "<<pin[i]<<std::endl;
      pinout[i] += pin[i];
    } 
  } 


  SparseMatrixStructure::SparseMatrixStructure(const DistSparseMatrixGraph & G){
    bIsExpanded = G.IsExpanded();
    baseval = G.GetBaseval();
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

    MPI_Comm comm = G.GetComm();
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
    std::vector<int> prevnz(np);
    std::vector<int> rcounts(np);
    Ptr lnnz = G.LocalEdgeCount();
    MPI_Allgather(&lnnz, 1, MPI_SYMPACK_PTR, &rcounts[0], 1, MPI_SYMPACK_PTR, comm);

    prevnz[0] = 0; 
    for (Int i = 0; i < np-1; ++i) { prevnz[i+1] = prevnz[i] + rcounts[i]; } 

    if(G.GetKeepDiag()){
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

    std::vector<int> rdispls = prevnz;
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
    if(!G.GetKeepDiag()){
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
      std::vector<Ptr> dummy;
      expColptr.swap(dummy);
    }
    {
      std::vector<Idx> dummy;
      expRowind.swap(dummy);
    }
    bIsExpanded=false;
  }





  void SparseMatrixStructure::ExpandSymmetric(MPI_Comm * workcomm){
    SYMPACK_TIMER_START(EXPAND);
    if(!bIsExpanded){
      if(bIsGlobal){

      //code from sparsematrixconverter
      /* set-up */

      std::vector<Ptr> cur_col_nnz(size);
      std::vector<Ptr> new_col_nnz(size);
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


//logfileptr->OFS()<<colptr<<std::endl;
//logfileptr->OFS()<<rowind<<std::endl;


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

//logfileptr->OFS()<<expColptr<<std::endl;
//logfileptr->OFS()<<expRowind<<std::endl;




      }
      else{
        throw std::logic_error("Please use DistSparseMatrixGraph instead.");
        //logfileptr->OFS()<<"N is"<<this->size<<std::endl;
        //logfileptr->OFS()<<"baseval is"<<this->baseval<<std::endl;
        //logfileptr->OFS()<<"keepDiag is"<<this->keepDiag<<std::endl;
        //logfileptr->OFS()<<"sorted is"<<this->sorted<<std::endl;
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
            Int maxLocN = std::max(locColCnt,N-(mpisize-1)*colPerProc); // can be 0

            //std::vector<Ptr> extra_nnz_percol(maxLocN,0);
            std::vector<Ptr> remote_colptr(maxLocN+1);
            std::vector<Idx> remote_rowind;
            std::vector<Ptr> remote_rowindPos(maxLocN+1);
            //std::vector<Idx> extra_rowind(maxLocN,0);
            std::vector<Ptr> curPos(locColCnt);
            std::vector<Ptr> prevPos(locColCnt);

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
                  SYMPACK_TIMER_START(REDUCE);
#ifdef USE_REDUCE
                  MPI_Reduce(&remote_colptr[0],&expColptr[0],remColCnt+1,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_SUM,prow,splitcomm);
#else
                  {
                    std::vector<Ptr> recv(remColCnt+1);
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
                  SYMPACK_TIMER_STOP(REDUCE);

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
                  SYMPACK_TIMER_START(REDUCE);
#ifdef USE_REDUCE
                  MPI_Reduce(&remote_colptr[0],NULL,remColCnt+1,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_SUM,prow,splitcomm);
#else
                  {
                    MPI_Send(&remote_colptr[0],(remColCnt+1)*sizeof(Ptr),MPI_BYTE,prow,mpisize+prow,splitcomm);
                  }
#endif
                  SYMPACK_TIMER_STOP(REDUCE);
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
                logfileptr->OFS()<<std::endl;

                logfileptr->OFS()<<"remote_rowind ";
                for(Ptr col = 0;col<extraNNZ;col++){
                  logfileptr->OFS()<<remote_rowind[col]<<" ";
                }
                logfileptr->OFS()<<std::endl;
#endif

                if(prow==mpirank){
#ifdef DEBUG
                  logfileptr->OFS()<<"expColptr "<<expColptr<<std::endl;
                  //logfileptr->OFS()<<"true expColptr "<<Global.expColptr<<std::endl;
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

                  SYMPACK_TIMER_START(RECV);
                  for(Int pcol = 0; pcol<prow; pcol++){
                    //Use an MPI_Gatherv instead ? >> memory usage : p * n/p
                    //Do mpi_recv from any, anytag for pcolptr and then do the matching rowind ?

                    //logfileptr->OFS()<<"P"<<mpirank<<" receives pcolptr from P"<<pcol<<std::endl;
                    //receive colptrs...
                    MPI_Status status;
                    MPI_Recv(&remote_colptr[0],(remColCnt+1)*sizeof(Ptr),MPI_BYTE,MPI_ANY_SOURCE,prow,splitcomm,&status);

                    //logfileptr->OFS()<<"P"<<mpirank<<" receives rowind from P"<<pcol<<std::endl;
                    //receive rowinds...
                    MPI_Recv(&remote_rowind[0],maxExtraNNZ*sizeof(Idx),MPI_BYTE,status.MPI_SOURCE,prow,splitcomm,MPI_STATUS_IGNORE);

                    SYMPACK_TIMER_START(PROCESSING_RECV_DATA);
                    //logfileptr->OFS()<<"P"<<mpirank<<" done receiving from P"<<pcol<<std::endl;
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
                    SYMPACK_TIMER_STOP(PROCESSING_RECV_DATA);
                  }
                  SYMPACK_TIMER_STOP(RECV);

                  if(sorted){ 
                    for(Idx locCol = 0 ; locCol< locColCnt; locCol++){
                      Ptr colbeg = expColptr[locCol]-baseval; //now 0 based
                      Ptr colend = expColptr[locCol+1]-baseval; // now 0 based 
                      sort(&expRowind[0]+colbeg,&expRowind[0]+colend,std::less<Ptr>());
                    }
                  }

#ifdef DEBUG
                  logfileptr->OFS()<<"expRowind "<<expRowind<<std::endl;
                  //logfileptr->OFS()<<"true expRowind "<<Global.expRowind<<std::endl;
#endif

                }
                else{
                  SYMPACK_TIMER_START(SEND);
                  //MPI_Send
                  //            logfileptr->OFS()<<"P"<<mpirank<<" sends pcolptr to P"<<prow<<std::endl;
                  MPI_Send(&remote_colptr[0],(remColCnt+1)*sizeof(Ptr),MPI_BYTE,prow,prow,splitcomm);
                  //            logfileptr->OFS()<<"P"<<mpirank<<" sends rowind to P"<<prow<<std::endl;
                  MPI_Send(&remote_rowind[0],extraNNZ*sizeof(Idx),MPI_BYTE,prow,prow,splitcomm);
                  //            logfileptr->OFS()<<"P"<<mpirank<<" done sending to P"<<prow<<std::endl;
                  SYMPACK_TIMER_STOP(SEND);
                }

              }
              else{
                SYMPACK_TIMER_START(REDUCE);
#ifdef USE_REDUCE
                MPI_Reduce(&remote_colptr[0],NULL,remColCnt+1,MPI_SYMPACK_PTR,MPI_SYMPACK_PTR_SUM,prow,splitcomm);
#endif
                SYMPACK_TIMER_STOP(REDUCE);
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
    SYMPACK_TIMER_STOP(EXPAND);

  }


  void SparseMatrixStructure::ToGlobal(SparseMatrixStructure & pGlobal,MPI_Comm & comm){

    SYMPACK_TIMER_START(ToGlobalStructure);
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
        std::vector<int> prevnz(np);
        std::vector<int> rcounts(np);
        MPI_Allgather(&nnz, sizeof(nnz), MPI_BYTE, &rcounts[0], sizeof(nnz), MPI_BYTE, comm);

        prevnz[0] = 0;
        for (Int i = 0; i < np-1; ++i) { prevnz[i+1] = prevnz[i] + rcounts[i]; } 

        pGlobal.nnz = 0;
        for (Int i = 0; i < np; ++i) { pGlobal.nnz += rcounts[i]; } 
        pGlobal.rowind.resize(pGlobal.nnz);

        std::vector<int> rdispls = prevnz;
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
        std::vector<int> prevnz(np);
        std::vector<int> rcounts(np);
        Ptr gnnz = expRowind.size();
        MPI_Allgather(&gnnz, sizeof(gnnz), MPI_BYTE, &rcounts[0], sizeof(gnnz), MPI_BYTE, comm);

        prevnz[0] = 0;
        for (Int i = 0; i < np-1; ++i) { prevnz[i+1] = prevnz[i] + rcounts[i]; } 

        Ptr totnnz=0;
        for (Int i = 0; i < np; ++i) { totnnz += rcounts[i]; } 
        pGlobal.expRowind.resize(totnnz);

        std::vector<int> rdispls = prevnz;
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

    SYMPACK_TIMER_STOP(ToGlobalStructure);
  }

}
