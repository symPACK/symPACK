/// @file FBMatrix.cpp
/// @brief Matrix structure and methods for the FAN-BOTH method.
/// @author Mathias Jacquelin
/// @date 2014-01-21

#include "FBMatrix.hpp"

#include  "utility.hpp"
#include  "blas.hpp"
#include  "lapack.hpp"

#include  "LogFile.hpp"
#include  "upcxx_additions.hpp"




namespace LIBCHOLESKY{


  FBMatrix::FBMatrix():blksize(1),n(0){
  }

  FBMatrix::~FBMatrix(){

    AchunkLower.clear();
    WLower.clear();
//    for(int i =0; i<AchunkLower.size();++i){
//      if(AchunkLower[i]!=NULL){
//        delete AchunkLower[i];
//      }
//    }
//    for(int i =0; i<WLower.size();++i){
//      if(WLower[i]!=NULL){
//        delete WLower[i];
//      }
//    }
  }

  void FBMatrix::Allocate(Int np, Int pn, Int pblksize){
    this->np = np;
    n=pn;
    blksize=pblksize;

    //determine pcol and prow
    Int firstproc = MAP(0,0);
    Int j=0;
    Int curproc = -1;
    pcol =0;
    do{
      ++j;
      if(j<pn){
        curproc=MAP(j,j);
//        logfileptr->OFS()<<"firstproc = "<<firstproc<<" vs "<<curproc<<endl;
        ++pcol;
      }
    }while(firstproc!=curproc && j<pn);

//    logfileptr->OFS()<<"PCOL IS "<<pcol<<endl;   
//    exit(0);

    Int totBlk = n/blksize;
    Int remaining = n-totBlk*blksize;
    if(remaining>0) totBlk++;

    Int localBlk = totBlk/pcol;
    Int additionalBlock = totBlk - pcol*localBlk;
    if(iam<additionalBlock){ localBlk++; }


    Int prevBlk = (iam) +pcol*(localBlk-1); 
    Int chksize = (localBlk-1)*blksize + min(n-(prevBlk)*blksize,blksize);

    Int lowerChkSize = (localBlk-1)*(n-iam*blksize)*blksize - blksize*blksize*pcol*(localBlk-1)*localBlk/2 + (n-(prevBlk)*blksize)*min(n-(prevBlk)*blksize,blksize); 

    Int numStep = ceil((double)n/(double)blksize);
    WLower.resize(numStep,NULL); //, DblNumMat(0,0) );
    AchunkLower.resize(localBlk,NULL);//, DblNumMat(0,0) );

    //WARNING this function needs to be completed by specialized classes
  } 

  void FBMatrix::Aggregate(Int j, DblNumMat &DistW){
    TIMER_START(Aggregate);
    //j is local
    Int local_j = this->global_col_to_local(j);
    DblNumMat & LocalChunk = *AchunkLower[local_j/blksize];

    Int jb = min(blksize, n-j);
    //aggregate previous updates

    //assert(DistW.Size()==LocalChunk.Size());
    blas::Axpy(LocalChunk.Size(), -1.0, &DistW(0,0),1,&LocalChunk(0,0),1);

    TIMER_STOP(Aggregate);
  }

  void FBMatrix::Factor(Int j){
    TIMER_START(Factor);
    Int jb = min(blksize, n-j);
    //j is local
    Int local_j = global_col_to_local(j);
    DblNumMat & LocalChunk = *AchunkLower[local_j/blksize];

    lapack::Potrf( 'L', jb, &LocalChunk(0,0 ), n-j);
    if(n-j-jb>0){
      blas::Trsm('R','L','T','N',n-j-jb,jb, 1.0, &LocalChunk(0,0), n-j, &LocalChunk(jb,0), n-j);
    }
    TIMER_STOP(Factor);
  }

  void FBMatrix::Update(Int j, Int i, DblNumMat & Factor){
    TIMER_START(Update);
    DblNumMat & WChunk = *WLower[i/blksize];
    //i is local, j is global
    Int jb = min(blksize, n-j);
    Int jbi = min(blksize, n-i);

    Int mf = Factor.m();
    Int nf = Factor.n();

    //call dgemm
    blas::Gemm('N','T',n-i,jbi,jb,1.0,&Factor(i-j,0),n-j,&Factor(i-j,0),n-j,1.0,&WChunk(0,0),n-i);
    TIMER_STOP(Update);
  }


  void FBMatrix::WaitFactorization(){     
  }

}






////#ifndef ADVANCE_COMM
////#define ADVANCE_COMM 1
////#endif
////
////
////#define TAG_FACTOR 0
////#define TAG_AGGREGATE 1
////
////
//////#define _DEBUG_
////
////namespace LIBCHOLESKY{
////
////
////
////#ifndef UPCXX
////  MPIGrid::MPIGrid( MPI_Comm Bcomm, int nprow, int npcol )
////  {
////    Int info;
////    MPI_Initialized( &info );
////    if( !info ){
////      throw std::logic_error( "MPI has not been initialized." );
////    }
////    MPI_Group  comm_group;
////    MPI_Comm_group( Bcomm, &comm_group );
////    MPI_Comm_create( Bcomm, comm_group, &comm );
////
////    MPI_Comm_rank( comm, &mpirank );
////    MPI_Comm_size( comm, &mpisize );
////    if( mpisize != nprow * npcol ){
////      throw std::logic_error( "mpisize != nprow * npcol." ); 
////    }
////
////    numProcRow = nprow;
////    numProcCol = npcol;
////
////    Int myrow = mpirank / npcol;
////    Int mycol = mpirank % npcol;
////
////    MPI_Comm_split( comm, myrow, mycol, &rowComm );
////    MPI_Comm_split( comm, mycol, myrow, &colComm );
////
////    MPI_Group_free( &comm_group );
////
////    return ;
////  } 		// -----  end of method MPIGrid::MPIGrid  ----- 
////
////
////  MPIGrid::~MPIGrid	(  )
////  {
////    MPI_Comm_free( &rowComm );
////    MPI_Comm_free( &colComm ); 
////    MPI_Comm_free( &comm );
////
////    return ;
////  } 		// -----  end of method MPIGrid::~MPIGrid  ----- 
////
////
////
////#endif
////
////
////
////
////
////
////#ifdef UPCXX
////  volatile int fact_finished = 0;
////  void signal_exit_am()
////  {
////    fact_finished = 1;
////  }
////
////  void signal_exit()
////  {
////    for (int i = 0; i < THREADS; i++) {
////      async(i)(signal_exit_am);
////    }
////  }
////#endif
////
////
////  FBMatrix::FBMatrix():prefetch(0),blksize(1),n(0),outstdAggreg(0){
////  }
////
////  FBMatrix::~FBMatrix(){
////  }
////
////#ifdef UPCXX
////  void FBMatrix::Initialize(upcxx::shared_array<upcxx::global_ptr<FBMatrix> > * RemoteObjPtr){
////    np=THREADS;
////    iam=MYTHREAD;
////    TIMER_START(RemotePtr_fetch);
////    RemoteObjPtrs.resize(RemoteObjPtr->size());
////    for(int i =0; i< RemoteObjPtr->size();i++){
////      RemoteObjPtrs[i]=(*RemoteObjPtr)[i];
////      logfileptr->OFS()<<RemoteObjPtrs[i]<<endl;
////    }
////    TIMER_STOP(RemotePtr_fetch);
////  }
////#else
////  void FBMatrix::Initialize(MPIGrid & grid ){
////    mpigrid = &grid; 
////    iam=grid.mpirank;
////    np=grid.mpisize;
////  }
////#endif
////
////  void FBMatrix::Allocate(Int pn, Int pblksize){
////    n=pn;
////    blksize=pblksize;
////
////    Int totBlk = n/blksize;
////    Int remaining = n-totBlk*blksize;
////    if(remaining>0) totBlk++;
////
////    Int localBlk = totBlk/pcol;
////    Int additionalBlock = totBlk - pcol*localBlk;
////    if(iam<additionalBlock){ localBlk++; }
////
////
////    Int prevBlk = (iam) +pcol*(localBlk-1); 
////    Int chksize = (localBlk-1)*blksize + min(n-(prevBlk)*blksize,blksize);
////
////    Int lowerChkSize = (localBlk-1)*(n-iam*blksize)*blksize - blksize*blksize*pcol*(localBlk-1)*localBlk/2 + (n-(prevBlk)*blksize)*min(n-(prevBlk)*blksize,blksize); 
////
////    Int numStep = ceil((double)n/(double)blksize);
////    WLower.resize(numStep, DblNumMat(0,0) );
////    AchunkLower.resize(localBlk, DblNumMat(0,0) );    
////    for(Int j = 0; j<n;j+=blksize){ 
////      Int local_j = global_col_to_local(j);
////      Int jb = min(blksize, n-j);
////      if(iam==MAP(j,j)){
////        AchunkLower[local_j/blksize].Resize(n-j,jb);
////      }
////      WLower[j/blksize].Resize(n-j,jb);
////      SetValue(WLower[j/blksize],0.0);
////    }
////  }
////
////  void FBMatrix::Distribute( DblNumMat & Aorig){
////    //    SetValue(Achunk,(double)iam);
////#ifdef UPCXX
////    for(Int j = 0; j<n;j+=blksize){ 
////      Int local_j = global_col_to_local(j);
////      Int jb = min(blksize, n-j);
////      if(iam==MAP(j,j)){
////        global_ptr<double> Aorigptr(&Aorig(j,j));
////        global_ptr<double> Adestptr(&AchunkLower[local_j/blksize](0,0));
////#ifdef ASYNC_COPY
////        upcxx::ldacopy_async(n-j,jb,Aorigptr,n,Adestptr,n-j);
////#else
////        upcxx::ldacopy(n-j,jb,Aorigptr,n,Adestptr,n-j);
////#endif
////      }
////    }
////
////#ifdef ASYNC_COPY
////    upcxx::async_copy_fence();
////#endif
////#else
////    for(Int j = 0; j<n;j+=blksize){ 
////      Int local_j = global_col_to_local(j);
////      Int jb = min(blksize, n-j);
////      if(iam==MAP(j,j)){
////        double * Aorigptr = &Aorig(j,j);
////        double * Adestptr = &AchunkLower[local_j/blksize](0,0);
////
////        lapack::Lacpy( 'N', n-j, jb, Aorigptr, n,	Adestptr, n-j	);
////      }
////    }
////#endif
////    Int numStep = ceil((double)n/(double)blksize);
////
////    if(iam==0){
////      cout<<"Factorization will require "<<numStep<<" steps"<<endl;
////    }
////
////    // TODO this has to be revised
////
////    //Allocate the lock vector
////    AggLock.Resize(numStep);
////    for(int js=0; js<numStep;js++){
////      int j = min(js*blksize,n-1);
////      Int jb = min(blksize, n-j);
////
////      IntNumVec procinvolved(np);
////      SetValue(procinvolved, I_ZERO);
////      Int numProc = 0;
////      for(int prevj = j-blksize;prevj>=0;prevj-=blksize){
////        if(procinvolved(MAP(j,prevj))==0){
////          procinvolved(MAP(j,prevj))=1;
////          numProc++;
////        }
////        if(numProc==np){
////          break;
////        }
////      }
////
////
////      AggLock[js]=numProc;
////    }
////  }
////
////  void FBMatrix::Gather( DblNumMat & Adest){
////#ifndef UPCXX
////    if(iam==0)
////#endif
////    {
////      Adest.Resize(n,n);
////      SetValue(Adest,0.0);
////    }
////    
////#ifdef UPCXX
////    for(Int j = 0; j<n;j+=blksize){ 
////      Int local_j = global_col_to_local(j);
////      Int jb = min(blksize, n-j);
////      Int target = MAP(j,j);
////      global_ptr<double> Adestptr(&Adest(j,j));
////      upcxx::async(target)(Gathercopy,RemoteObjPtrs[target],Adestptr,j);
////    }
////    upcxx::wait();
////#else
////    DblNumMat recvBuf;
////    if(iam==0){
////      recvBuf.Resize(n,blksize);
////    }
////
////    for(Int j = 0; j<n;j+=blksize){ 
////      Int local_j = global_col_to_local(j);
////      Int jb = min(blksize, n-j);
////      Int target = MAP(j,j);
////
////      if(iam==0){
////        if(target==0){
////          //do a local copy
////#ifdef _DEBUG_
////          logfileptr->OFS()<<"Copying local column"<<j<<endl;
////#endif
////          DblNumMat & LocalSrcChunk = AchunkLower[local_j/blksize];
////          lapack::Lacpy( 'N', n-j, jb, LocalSrcChunk.Data(), n-j,	&Adest(j,j), n	);
////        }
////        else{
////#ifdef _DEBUG_
////logfileptr->OFS()<<"Receiving column"<<j<<" from P"<<target<<endl;
////#endif
////          //recv in the temp buffer
////          MPI_Recv(recvBuf.Data(),(n-j)*jb*sizeof(recvBuf(0,0)),MPI_BYTE,target,TAG_FACTOR,mpigrid->comm,MPI_STATUS_IGNORE);
////          //do the final lda_cpy to the output matrix
////#ifdef _DEBUG_
////logfileptr->OFS()<<"Copying remote column"<<j<<endl;
////#endif
////          lapack::Lacpy( 'N', n-j, jb, recvBuf.Data(), n-j,	&Adest(j,j), n	);
////        }
////      }
////      else{
////        if(iam==target){
////#ifdef _DEBUG_
////logfileptr->OFS()<<"Sending column"<<j<<" to P0"<<endl;
////#endif
////          DblNumMat & LocalSrcChunk = AchunkLower[local_j/blksize];
////          MPI_Send(LocalSrcChunk.Data(),LocalSrcChunk.ByteSize(),MPI_BYTE,0,TAG_FACTOR,mpigrid->comm);
////        }
////      }
////    }
////#endif
////  }
////
////
////
////  void FBMatrix::Aggregate(Int j, DblNumMat &DistW){
////    TIMER_START(Aggregate);
////    //j is local
////    Int local_j = this->global_col_to_local(j);
////    DblNumMat & LocalChunk = AchunkLower[local_j/blksize];
////
////    Int jb = min(blksize, n-j);
////    //aggregate previous updates
////
////    //assert(DistW.Size()==LocalChunk.Size());
////    blas::Axpy(LocalChunk.Size(), -1.0, &DistW(0,0),1,&LocalChunk(0,0),1);
////
////    TIMER_STOP(Aggregate);
////  }
////
////
////
////
////  void FBMatrix::Factor(Int j){
////    TIMER_START(Factor);
////    Int jb = min(blksize, n-j);
////    //j is local
////    Int local_j = global_col_to_local(j);
////    DblNumMat & LocalChunk = AchunkLower[local_j/blksize];
////
////    lapack::Potrf( 'L', jb, &LocalChunk(0,0 ), n-j);
////    if(n-j-jb>0){
////      blas::Trsm('R','L','T','N',n-j-jb,jb, 1.0, &LocalChunk(0,0), n-j, &LocalChunk(jb,0), n-j);
////    }
////    TIMER_STOP(Factor);
////  }
////
////  void FBMatrix::Update(Int j, Int i, DblNumMat & Factor){
////    TIMER_START(Update);
////    DblNumMat & WChunk = WLower[i/blksize];
////    //i is local, j is global
////    Int jb = min(blksize, n-j);
////    Int jbi = min(blksize, n-i);
////
////    //call dgemm
////    blas::Gemm('N','T',n-i,jbi,jb,1.0,&Factor(i-j,0),n-j,&Factor(i-j,0),n-j,1.0,&WChunk(0,0),n-i);
////    TIMER_STOP(Update);
////  }
////
////
////
////
////  void FBMatrix::WaitFactorization(){     
////
////#ifdef UPCXX
////    //while(!fact_finished){  upcxx::progress(); };
////    while(!fact_finished){  upcxx::advance(1, 99); };
////    if(MYTHREAD==MAP(n-1,n-1)){
////      TIMER_STOP(SYNC_END);
////    }
////#endif
////  }
////
////#ifdef UPCXX
////  void FBMatrix::NumericalFactorization(){
////    if(iam==MAP(0,0)){
////      upcxx::async(iam)(Factor_Async,RemoteObjPtrs[iam],0);
////    }
////  }
////#else
////  void FBMatrix::NumericalFactorization(){
////      Int numStep = ceil((double)n/(double)blksize);
////
////      std::vector<MPI_Request> Factor_Send_Requests(np,MPI_REQUEST_NULL);
////      std::vector<MPI_Request> Aggregate_Send_Requests(np,MPI_REQUEST_NULL);
////
////      for(int js=0; js<numStep;js++){
////        int j = min(js*blksize,n-1);
////        Int jb = min(blksize, n-j+1);
////        Int local_j = global_col_to_local(j);
////
////        //do the aggregation and factorization
////        if(iam==MAP(j,j) ){
////
////          //aggregate previous updates
////          //Do I have local aggregates ?
////          for(int i=(js-1)*blksize; i>=max(0,(js-np))*blksize; i-=blksize){
////
////            #ifdef _DEBUG_
////               logfileptr->OFS()<<i<<endl;
////            #endif
////            if(iam==MAP(j,i)){
////              //This is a local aggregate
////              DblNumMat & RemoteAggregate = WLower[j/blksize];
////
////                #ifdef _DEBUG_
////                  logfileptr->OFS()<<"Aggregating local to column "<<j<<endl;
////                #endif
////              Aggregate(j, RemoteAggregate);
////              --AggLock[j/blksize];
////              break;
////            }
////          }
////          //Process remote aggregates
////          if(AggLock[j/blksize]>0){
////            DblNumMat RemoteAggregate(n-j,jb);
////            while(AggLock[j/blksize]>0){
////
////                #ifdef _DEBUG_
////                  logfileptr->OFS()<<"Receiving aggregate"<<endl;
////                #endif
////              MPI_Recv(RemoteAggregate.Data(),RemoteAggregate.ByteSize(),MPI_BYTE,MPI_ANY_SOURCE,TAG_AGGREGATE,mpigrid->comm,MPI_STATUS_IGNORE);
////
////                #ifdef _DEBUG_
////                  logfileptr->OFS()<<"Aggregating to column "<<j<<endl;
////                #endif
////              Aggregate(j, RemoteAggregate);
////              --AggLock[j/blksize];
////            }
////          }
////
////
////                #ifdef _DEBUG_
////                  logfileptr->OFS()<<"Factoring column "<<j<<endl;
////                #endif
////          //Factor current column 
////          Factor(j);
////
////          //send factor
////          TIMER_START(Send_Factor);
////          for(Int target = 0; target<np;target++){
////            for(Int i = j+blksize; i<n; i+=blksize){
////              if(MAP(i,j)==target && target!=iam){
////                DblNumMat & Achk = AchunkLower[local_j/blksize];
////
////                #ifdef _DEBUG_
////                  logfileptr->OFS()<<"Launching update from "<<j<<" on P"<<target<<endl;
////                #endif
//////                MPI_Send(&Achk(0,0),Achk.ByteSize(),MPI_BYTE,target,TAG_FACTOR, mpigrid->comm );
////                MPI_Wait(&Factor_Send_Requests[target],MPI_STATUS_IGNORE);
////                MPI_Isend(&Achk(0,0),Achk.ByteSize(),MPI_BYTE,target,TAG_FACTOR, mpigrid->comm,&Factor_Send_Requests[target] );
////                break;
////              }
////            }
////          }
////          TIMER_STOP(Send_Factor);
////
////        }
////
////        //MPI_Barrier(mpigrid->comm);
////
////        //do the updates
////
////        DblNumMat RemoteFactorBuf(n-j,blksize);
////        bool have_factor = false;
////
////        for(int i=j+jb;i<n ;i+=jb){
////          //compute the update
////          if(iam==MAP(i,j)){
////            DblNumMat * RemoteFactorPtr;
////            if(iam!=MAP(j,j) ){
////              RemoteFactorPtr = & RemoteFactorBuf;
////              if(!have_factor){
////#ifdef _DEBUG_
////                logfileptr->OFS()<<"Receiving factor "<<j<<endl;
////#endif
////                //Receive factor
////                MPI_Recv(RemoteFactorPtr->Data(),RemoteFactorPtr->ByteSize(),MPI_BYTE,MAP(j,j),TAG_FACTOR,mpigrid->comm,MPI_STATUS_IGNORE);
////                have_factor=true;
////              }
////            }
////            else{
////              RemoteFactorPtr = & AchunkLower[local_j/blksize];
////            }
////
////                #ifdef _DEBUG_
////                  logfileptr->OFS()<<"Updating "<<i<<" with "<<j<<endl;
////                #endif
////            //Do the update
////            Update(j, i, *RemoteFactorPtr);
////
////            TIMER_START(Launch_Aggregates);
////            //is it my last update to column i ?
////            bool last_upd = true;
////            for(Int jj = j+blksize; jj<i;jj+=blksize){
////              if(iam==MAP(i,jj)){
////                last_upd=false;
////                break;
////              }
////            }
////
////            //if this is the last update, send aggregate
////            if(last_upd){
////              Int target = MAP(i,i);
////#ifdef _DEBUG_
////logfileptr->OFS()<<"I did my last update to "<<i<<endl;
////#endif
////
////              if(target!=iam){
////#ifdef _DEBUG_
////              logfileptr->OFS()<<"Updating from "<<j<<", sending aggregate to P"<<target<<" for column "<<i<<endl;
////#endif
//////              MPI_Send(WLower[i/blksize].Data(),WLower[i/blksize].ByteSize(),MPI_BYTE,target,TAG_AGGREGATE,mpigrid->comm);
////                MPI_Wait(&Aggregate_Send_Requests[target],MPI_STATUS_IGNORE);
////                MPI_Isend(WLower[i/blksize].Data(),WLower[i/blksize].ByteSize(),MPI_BYTE,target,TAG_AGGREGATE,mpigrid->comm, &Aggregate_Send_Requests[target]);
////              }
////            }
////            TIMER_STOP(Launch_Aggregates);
////
////
////
////          }
////        }
////
////        //MPI_Barrier(mpigrid->comm);
////
////      } 
////      
////      MPI_Barrier(mpigrid->comm);
////
////  }
////
////#endif
////
////#ifdef UPCXX
////  void Gathercopy(upcxx::global_ptr<FBMatrix> Objptr, global_ptr<double> Adest, Int j){
////    assert(Objptr.tid()==MYTHREAD);
////    FBMatrix & A = *Objptr;
////    Int local_j = A.global_col_to_local(j);
////    Int jb = min(A.blksize, A.n-j);
////    DblNumMat & LocalChunk = A.AchunkLower[local_j/A.blksize];
////    global_ptr<double> Asrc = &LocalChunk(0,0);
////#ifdef ASYNC_COPY
////    upcxx::ldacopy_async(A.n-j,jb,Asrc,A.n-j,Adest,A.n);
////    upcxx::async_copy_fence();
////#else
////    upcxx::ldacopy(A.n-j,jb,Asrc,A.n-j,Adest,A.n);
////#endif
////
////  }
////
////
////  void Factor_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j){
////    TIMER_START(Factor_Async);
////#ifdef _ASSERT_
////    assert(Aptr.tid()==MYTHREAD);
////#endif
////    FBMatrix & A = *Aptr;
////
////#ifdef _DEBUG_
////    logfileptr->OFS()<<"********************************************"<<endl;
////    //Factor the column
////    logfileptr->OFS()<<"Factoring column "<<j<<endl;
////#endif
////
////    A.Factor(j);
////
////#ifdef _DEBUG_
////    logfileptr->OFS()<<"Factoring column done "<<j<<endl;
////
////    Int numStep = ceil((double)A.n/(double)A.blksize);
////    Int curStep = j/A.blksize;
////
////    logfileptr->OFS()<<"STEP "<<curStep<<"/"<<numStep<<endl;
////#endif
////    //Launch the updates
////    Int local_j = A.global_col_to_local(j);
////    //TODO CHECK THIS: WE NEED TO LAUNCH THE UPDATES ON ALL PROCESSORS
////    TIMER_START(Launch_Updates);
////    for(Int target = 0; target<A.np;target++){
////      for(Int i = j+A.blksize; i<A.n; i+=A.blksize){
////        if(A.MAP(i,j)==target){
////          global_ptr<double> AchkPtr(&A.AchunkLower[local_j/A.blksize](0,0));
////
////#ifdef _DEBUG_
////          logfileptr->OFS()<<"Launching update from "<<j<<" on P"<<target<<endl;
////#endif
////          upcxx::async(target)(Update_Async,A.RemoteObjPtrs[target],j,AchkPtr);
//////          upcxx::advance(0,ADVANCE_COMM);
////          break;
////        }
////      }
////    }
////    TIMER_STOP(Launch_Updates);
////
////
////
////    if(j+A.blksize>=A.n){
////#ifdef _DEBUG_
////      logfileptr->OFS()<<"Unlocking quit"<<endl;
////#endif
////
////      TIMER_START(SYNC_END);
////      signal_exit();
////
////    }
////
////#ifdef _DEBUG_
////    logfileptr->OFS()<<"Quitting Factor_Async"<<endl;
////    logfileptr->OFS()<<"********************************************"<<endl;
////#endif
////    TIMER_STOP(Factor_Async);
////  }
////
////
////
////
////
////
////
////
////
////#ifndef ASYNC_COPIES
////  void Update_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, upcxx::global_ptr<double> remoteFactorPtr){
////
////    TIMER_START(Update_Async);
////#ifdef _ASSERT_
////    assert(Aptr.tid()==MYTHREAD);
////#endif
////    FBMatrix & A = *Aptr;
////
////
////#ifdef _DEBUG_    
////    logfileptr->OFS()<<"Fetching factor "<<j<<" from P"<<remoteFactorPtr.tid()<<endl;
////#endif
////    DblNumMat RemoteFactor(A.n-j,A.blksize);
////    TIMER_START(Fetching_factor);
////    upcxx::copy(remoteFactorPtr,RemoteFactor.GData(),RemoteFactor.Size());
////    TIMER_STOP(Fetching_factor);
////
////#ifdef _DEBUG_    
////        logfileptr->OFS()<<"Fetched factor "<<j<<endl;
////#endif
////
////    TIMER_START(Update_async_loop);
////    //do all my updates with i>j
////    for(Int i = j+A.blksize; i<A.n;i+=A.blksize){
////      if(A.iam==A.MAP(i,j)){
////#ifdef _DEBUG_
////        //        logfileptr->OFS()<<"Updating column "<<i<<" with columns from P"<<remoteFactorPtr.tid()<<endl;
////#endif
////        A.Update(j, i, RemoteFactor);
////
////
////
////        TIMER_START(Launch_Aggregates);
////        //is it my last update to column i ?
////        bool last_upd = true;
////        for(Int jj = j+A.blksize; jj<i;jj+=A.blksize){
////          if(A.iam==A.MAP(i,jj)){
////            //logfileptr->OFS()<<j<<" I also have to update column "<<jj<<endl;
////            last_upd=false;
////            break;
////          }
////        }
////        if(last_upd){
////          Int target = A.MAP(i,i);
////#ifdef _DEBUG_
////          logfileptr->OFS()<<"Updating from "<<j<<", launching aggregate on P"<<target<<" for column "<<i<<endl;
////#endif
////          global_ptr<double> WbufPtr(&A.WLower[i/A.blksize](0,0));
////          upcxx::async(target)(Aggregate_Async,A.RemoteObjPtrs[target],i,WbufPtr);
////        }
////        TIMER_STOP(Launch_Aggregates);
////
////
////      }
////    }
////    TIMER_STOP(Update_async_loop);
////
////
////#ifdef _DEBUG_
////    logfileptr->OFS()<<"Quitting Update_Async"<<endl<<endl;
////#endif
////    TIMER_STOP(Update_Async);
////  }
////
////
////
////
////  void Aggregate_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, global_ptr<double> remoteAggregatePtr){
////    TIMER_START(Aggregate_Async);
////#ifdef _ASSERT_
////    assert(Aptr.tid()==MYTHREAD);
////#endif
////    FBMatrix & A = *Aptr;
////    //fetch data
////#ifdef _DEBUG_    
////    logfileptr->OFS()<<"Aggregating Fetching data from P"<<remoteAggregatePtr.tid()<<endl;
////#endif
////    Int jb = min(A.blksize, A.n-j);
////
////    DblNumMat RemoteAggregate(A.n-j,jb); 
////    TIMER_START(Fetching_Aggregate);
////    upcxx::copy(remoteAggregatePtr,RemoteAggregate.GData(),RemoteAggregate.Size());
////    TIMER_STOP(Fetching_Aggregate);
////
////
////#ifdef _DEBUG_    
////    logfileptr->OFS()<<"Aggregating update to "<<j<<endl;
////#endif
////    A.Aggregate(j, RemoteAggregate);
////
////#ifdef _DEBUG_    
////    logfileptr->OFS()<<"Done Aggregating update to "<<j<<endl;
////#endif
////
////    A.AggLock[j/A.blksize]--;
////#ifdef _DEBUG_    
////    logfileptr->OFS()<<"Still waiting for "<<A.AggLock[j/A.blksize]<<" aggregates"<<endl;
////#endif
////    if(A.AggLock[j/A.blksize]==0){
////      //launch the factorization of column j
////#ifdef _DEBUG_    
////      logfileptr->OFS()<<"Col "<<j<<" updated. Can be factored "<<endl;
////#endif
////
////      upcxx::async(A.iam)(Factor_Async,Aptr,j);
////    } 
////#ifdef _DEBUG_
////    logfileptr->OFS()<<"Quitting Aggregate_Async"<<endl<<endl;
////#endif
////    TIMER_STOP(Aggregate_Async);
////  }
////
////
////
////
////
////
////
////#else
////
////
////
////  void Update_Compute_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, DblNumMat * remoteFactorPtr, upcxx::event * async_copy_event){
////
////    FBMatrix & A = *Aptr;
////
////    //sync with the async_copy
////    async_copy_event->wait();
////    delete async_copy_event;
////    TIMER_STOP(Fetching_factor);
////
////    A.outstdUpdate--;
////
////
////
////#ifdef _DEBUG_    
////        logfileptr->OFS()<<"Fetched factor "<<j<<endl;
////#endif
////
////    TIMER_START(Update_async_loop);
////    //do all my updates with i>j
////    for(Int i = j+A.blksize; i<A.n;i+=A.blksize){
////      if(A.iam==A.MAP(i,j)){
////#ifdef _DEBUG_
////        //        logfileptr->OFS()<<"Updating column "<<i<<" with columns from P"<<remoteFactorPtr.tid()<<endl;
////#endif
////        A.Update(j, i, *remoteFactorPtr);
////
////
////
////        TIMER_START(Launch_Aggregates);
////        //is it my last update to column i ?
////        bool last_upd = true;
////        for(Int jj = j+A.blksize; jj<i;jj+=A.blksize){
////          if(A.iam==A.MAP(i,jj)){
////            //logfileptr->OFS()<<j<<" I also have to update column "<<jj<<endl;
////            last_upd=false;
////            break;
////          }
////        }
////        if(last_upd){
////          Int target = A.MAP(i,i);
////#ifdef _DEBUG_
////          logfileptr->OFS()<<"Updating from "<<j<<", launching aggregate on P"<<target<<" for column "<<i<<endl;
////#endif
////          global_ptr<double> WbufPtr(&A.WLower[i/A.blksize](0,0));
////          upcxx::async(target)(Aggregate_Async,A.RemoteObjPtrs[target],i,WbufPtr);
//////          upcxx::advance(0,ADVANCE_COMM);
////          
////        }
////        TIMER_STOP(Launch_Aggregates);
////
////
////      }
////    }
////    TIMER_STOP(Update_async_loop);
////
////
////    //The factor can now be freed
////    delete remoteFactorPtr;
////    
////
////#ifdef _DEBUG_
////    logfileptr->OFS()<<"Quitting Update_Async"<<endl<<endl;
////#endif
//////    TIMER_STOP(Update_Async);
////  }
////
////
////
////
////
////  void Update_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, upcxx::global_ptr<double> remoteFactorPtr){
////
//////    TIMER_START(Update_Async);
////#ifdef _ASSERT_
////    assert(Aptr.tid()==MYTHREAD);
////#endif
////    FBMatrix & A = *Aptr;
////
////
////#ifdef _DEBUG_    
////    logfileptr->OFS()<<"Fetching factor "<<j<<" from P"<<remoteFactorPtr.tid()<<endl;
////#endif
////    DblNumMat * RemoteFactor = new DblNumMat(A.n-j,A.blksize);
////    upcxx::event * async_copy_event = new upcxx::event;
////    TIMER_START(Fetching_factor);
////    upcxx::async_copy(remoteFactorPtr,RemoteFactor->GData(),RemoteFactor->Size(),async_copy_event);
////
////    A.outstdUpdate++;
////    if(A.outstdUpdate<A.prefetch){
////      //add the function to the async queue
////      upcxx::async(A.iam)(Update_Compute_Async,Aptr,j,RemoteFactor, async_copy_event);
////    }
////    else{
////      //call the function inline
////      Update_Compute_Async(Aptr, j, RemoteFactor, async_copy_event);
////    }
////
////
////}
////
////
////
////  void Aggregate_Compute_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, DblNumMat * remoteAggregatePtr, upcxx::event * async_copy_event){
////
////#ifdef _ASSERT_
////    assert(Aptr.tid()==MYTHREAD);
////#endif
////    FBMatrix & A = *Aptr;
////
////    //sync with the async_copy
////    async_copy_event->wait();
////    delete async_copy_event;
////    TIMER_STOP(Fetching_Aggregate);
////
////    A.outstdAggreg--;
////#ifdef _DEBUG_    
////    logfileptr->OFS()<<"Aggregating update to "<<j<<endl;
////#endif
////    A.Aggregate(j, *remoteAggregatePtr);
////    //The aggregate can now be freed
////    delete remoteAggregatePtr;
////
////#ifdef _DEBUG_    
////    logfileptr->OFS()<<"Done Aggregating update to "<<j<<endl;
////#endif
////
////    A.AggLock[j/A.blksize]--;
////#ifdef _DEBUG_    
////    logfileptr->OFS()<<"Still waiting for "<<A.AggLock[j/A.blksize]<<" aggregates"<<endl;
////#endif
////    if(A.AggLock[j/A.blksize]==0){
////      //launch the factorization of column j
////#ifdef _DEBUG_    
////      logfileptr->OFS()<<"Col "<<j<<" updated. Can be factored "<<endl;
////#endif
////
////      upcxx::async(A.iam)(Factor_Async,Aptr,j);
////    } 
////#ifdef _DEBUG_
////    logfileptr->OFS()<<"Quitting Aggregate_Async"<<endl<<endl;
////#endif
//////    TIMER_STOP(Aggregate_Async);
////  }
////
////  void Aggregate_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, global_ptr<double> remoteAggregatePtr){
//////    TIMER_START(Aggregate_Async);
////#ifdef _ASSERT_
////    assert(Aptr.tid()==MYTHREAD);
////#endif
////    FBMatrix & A = *Aptr;
////    //fetch data
////#ifdef _DEBUG_    
////    logfileptr->OFS()<<"Aggregating Fetching data from P"<<remoteAggregatePtr.tid()<<endl;
////#endif
////    Int jb = min(A.blksize, A.n-j);
////
////    DblNumMat * RemoteAggregate =  new DblNumMat(A.n-j,jb);
////
////    upcxx::event * async_copy_event = new upcxx::event;
////    TIMER_START(Fetching_Aggregate);
////    upcxx::async_copy(remoteAggregatePtr,RemoteAggregate->GData(),RemoteAggregate->Size(),async_copy_event);
////    
////    A.outstdAggreg++;
////    if(A.outstdAggreg<A.prefetch){
////      //add the function to the async queue
////      upcxx::async(A.iam)(Aggregate_Compute_Async,Aptr,j,RemoteAggregate, async_copy_event);
////    }
////    else{
////      //call the function inline
////      Aggregate_Compute_Async(Aptr, j, RemoteAggregate, async_copy_event);
////    }
////}
////
////
////
////#endif
////
////
////#endif
////}


