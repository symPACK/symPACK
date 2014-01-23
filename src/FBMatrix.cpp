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

//#define _DEBUG_

namespace LIBCHOLESKY{

  volatile int fact_finished = 0;
  void signal_exit_am()
  {
    fact_finished = 1;
  }

  void signal_exit()
  {
    for (int i = 0; i < THREADS; i++) {
      async(i)(signal_exit_am);
    }
  }



  FBMatrix::FBMatrix():blksize(1),n(0){
    np=THREADS;
    iam=MYTHREAD;

  }

  FBMatrix::~FBMatrix(){
  }


  void FBMatrix::Initialize(upcxx::shared_array<upcxx::global_ptr<FBMatrix> > * RemoteObjPtr){
    TIMER_START(RemotePtr_fetch);
    RemoteObjPtrs.resize(RemoteObjPtr->size());
    for(int i =0; i< RemoteObjPtr->size();i++){
      RemoteObjPtrs[i]=(*RemoteObjPtr)[i];
      logfileptr->OFS()<<RemoteObjPtrs[i]<<endl;
    }
    TIMER_STOP(RemotePtr_fetch);
  }


  void FBMatrix::Allocate(Int pn, Int pblksize){
    n=pn;
    blksize=pblksize;

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
    WLower.resize(numStep, DblNumMat(0,0) );
    AchunkLower.resize(localBlk, DblNumMat(0,0) );    
    for(Int j = 0; j<n;j+=blksize){ 
      Int local_j = global_col_to_local(j);
      Int jb = min(blksize, n-j);
      if(iam==MAP(j,j)){
        AchunkLower[local_j/blksize].Resize(n-j,jb);
      }
      WLower[j/blksize].Resize(n-j,jb);
      SetValue(WLower[j/blksize],0.0);
    }
  }

  void FBMatrix::Distribute( DblNumMat & Aorig){
    //    SetValue(Achunk,(double)iam);
    for(Int j = 0; j<n;j+=blksize){ 
      Int local_j = global_col_to_local(j);
      Int jb = min(blksize, n-j);
      if(iam==MAP(j,j)){
        global_ptr<double> Aorigptr(&Aorig(j,j));
        global_ptr<double> Adestptr(&AchunkLower[local_j/blksize](0,0));
#ifdef ASYNC_COPY
        upcxx::ldacopy_async(n-j,jb,Aorigptr,n,Adestptr,n-j);
#else
        upcxx::ldacopy(n-j,jb,Aorigptr,n,Adestptr,n-j);
#endif
      }
    }

#ifdef ASYNC_COPY
    upcxx::async_copy_fence();
#endif

    Int numStep = ceil((double)n/(double)blksize);

    if(iam==0){
      cout<<"Factorization will require "<<numStep<<" steps"<<endl;
    }

    // TODO this has to be revised

    //Allocate the lock vector
    AggLock.Resize(numStep);
    for(int js=0; js<numStep;js++){
      int j = min(js*blksize,n-1);
      Int jb = min(blksize, n-j);

      IntNumVec procinvolved(np);
      SetValue(procinvolved, I_ZERO);
      Int numProc = 0;
      for(int prevj = j-blksize;prevj>=0;prevj-=blksize){
        if(procinvolved(MAP(j,prevj))==0){
          procinvolved(MAP(j,prevj))=1;
          numProc++;
        }
        if(numProc==np){
          break;
        }
      }


      AggLock[js]=numProc;
    }
  }

  void FBMatrix::Gather( DblNumMat & Adest){
    Adest.Resize(n,n);
    SetValue(Adest,0.0);

    for(Int j = 0; j<n;j+=blksize){ 
      Int local_j = global_col_to_local(j);
      Int jb = min(blksize, n-j);
      Int target = MAP(j,j);
      global_ptr<double> Adestptr(&Adest(j,j));
      upcxx::async(target)(Gathercopy,RemoteObjPtrs[target],Adestptr,j);
    }
    upcxx::wait();
  }



  void FBMatrix::Aggregate(Int j, DblNumMat &DistW){
    TIMER_START(Aggregate);
    //j is local
    Int local_j = this->global_col_to_local(j);
    DblNumMat & LocalChunk = AchunkLower[local_j/blksize];

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
    DblNumMat & LocalChunk = AchunkLower[local_j/blksize];

    lapack::Potrf( 'L', jb, &LocalChunk(0,0 ), n-j);
    if(n-j-jb>0){
      blas::Trsm('R','L','T','N',n-j-jb,jb, 1.0, &LocalChunk(0,0), n-j, &LocalChunk(jb,0), n-j);
    }
    TIMER_STOP(Factor);
  }

  void FBMatrix::Update(Int j, Int i, DblNumMat & Factor){
    TIMER_START(Update);
    DblNumMat & WChunk = WLower[i/blksize];
    //i is local, j is global
    Int jb = min(blksize, n-j);
    Int jbi = min(blksize, n-i);

    //call dgemm
    blas::Gemm('N','T',n-i,jbi,jb,1.0,&Factor(i-j,0),n-j,&Factor(i-j,0),n-j,1.0,&WChunk(0,0),n-i);
    TIMER_STOP(Update);
  }




  void FBMatrix::WaitFactorization(){     
    while(!fact_finished){  upcxx::progress(); };

    if(MYTHREAD==MAP(n-1,n-1)){
      TIMER_STOP(SYNC_END);
    }
  }


  void FBMatrix::NumericalFactorization(){
    if(iam==MAP(0,0)){
      upcxx::async(iam)(Factor_Async,RemoteObjPtrs[iam],0);
    }
  }



  void Gathercopy(upcxx::global_ptr<FBMatrix> Objptr, global_ptr<double> Adest, Int j){
    assert(Objptr.tid()==MYTHREAD);
    FBMatrix & A = *Objptr;
    Int local_j = A.global_col_to_local(j);
    Int jb = min(A.blksize, A.n-j);
    DblNumMat & LocalChunk = A.AchunkLower[local_j/A.blksize];
    global_ptr<double> Asrc = &LocalChunk(0,0);
#ifdef ASYNC_COPY
    upcxx::ldacopy_async(A.n-j,jb,Asrc,A.n-j,Adest,A.n);
    upcxx::async_copy_fence();
#else
    upcxx::ldacopy(A.n-j,jb,Asrc,A.n-j,Adest,A.n);
#endif

  }











  void Factor_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j){
    TIMER_START(Factor_Async);
#ifdef _ASSERT_
    assert(Aptr.tid()==MYTHREAD);
#endif
    FBMatrix & A = *Aptr;

#ifdef _DEBUG_
    logfileptr->OFS()<<"********************************************"<<endl;
    //Factor the column
    logfileptr->OFS()<<"Factoring column "<<j<<endl;
#endif

    A.Factor(j);

#ifdef _DEBUG_
    logfileptr->OFS()<<"Factoring column done "<<j<<endl;

    Int numStep = ceil((double)A.n/(double)A.blksize);
    Int curStep = j/A.blksize;

    logfileptr->OFS()<<"STEP "<<curStep<<"/"<<numStep<<endl;
#endif
    //Launch the updates
    Int local_j = A.global_col_to_local(j);
    //TODO CHECK THIS: WE NEED TO LAUNCH THE UPDATES ON ALL PROCESSORS
    TIMER_START(Launch_Updates);
    for(Int target = 0; target<A.np;target++){
      for(Int i = j+A.blksize; i<A.n; i+=A.blksize){
        if(A.MAP(i,j)==target){
          global_ptr<double> AchkPtr(&A.AchunkLower[local_j/A.blksize](0,0));

#ifdef _DEBUG_
          logfileptr->OFS()<<"Launching update from "<<j<<" on P"<<target<<endl;
#endif
          upcxx::async(target)(Update_Async,A.RemoteObjPtrs[target],j,AchkPtr);
          break;
        }
      }
    }
    TIMER_STOP(Launch_Updates);



    if(j+A.blksize>=A.n){
#ifdef _DEBUG_
      logfileptr->OFS()<<"Unlocking quit"<<endl;
#endif

      TIMER_START(SYNC_END);
      signal_exit();

    }

#ifdef _DEBUG_
    logfileptr->OFS()<<"Quitting Factor_Async"<<endl;
    logfileptr->OFS()<<"********************************************"<<endl;
#endif
    TIMER_STOP(Factor_Async);
  }


  void Aggregate_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, global_ptr<double> remoteAggregatePtr){
    TIMER_START(Aggregate_Async);
#ifdef _ASSERT_
    assert(Aptr.tid()==MYTHREAD);
#endif
    FBMatrix & A = *Aptr;
    //fetch data
#ifdef _DEBUG_    
    logfileptr->OFS()<<"Aggregating Fetching data from P"<<remoteAggregatePtr.tid()<<endl;
#endif
    Int jb = min(A.blksize, A.n-j);

    DblNumMat RemoteAggregate(A.n-j,jb); 
    TIMER_START(Fetching_Aggregate);
    upcxx::copy(remoteAggregatePtr,RemoteAggregate.GData(),RemoteAggregate.Size());
    TIMER_STOP(Fetching_Aggregate);


#ifdef _DEBUG_    
    logfileptr->OFS()<<"Aggregating update to "<<j<<endl;
#endif
    A.Aggregate(j, RemoteAggregate);

#ifdef _DEBUG_    
    logfileptr->OFS()<<"Done Aggregating update to "<<j<<endl;
#endif

    A.AggLock[j/A.blksize]--;
#ifdef _DEBUG_    
    logfileptr->OFS()<<"Still waiting for "<<A.AggLock[j/A.blksize]<<" aggregates"<<endl;
#endif
    if(A.AggLock[j/A.blksize]==0){
      //launch the factorization of column j
#ifdef _DEBUG_    
      logfileptr->OFS()<<"Col "<<j<<" updated. Can be factored "<<endl;
#endif

      upcxx::async(A.iam)(Factor_Async,Aptr,j);
    } 
#ifdef _DEBUG_
    logfileptr->OFS()<<"Quitting Aggregate_Async"<<endl<<endl;
#endif
    TIMER_STOP(Aggregate_Async);
  }


  void Update_Async(upcxx::global_ptr<FBMatrix> Aptr, Int j, upcxx::global_ptr<double> remoteFactorPtr){

    TIMER_START(Update_Async);
#ifdef _ASSERT_
    assert(Aptr.tid()==MYTHREAD);
#endif
    FBMatrix & A = *Aptr;


#ifdef _DEBUG_    
    logfileptr->OFS()<<"Fetching factor "<<j<<" from P"<<remoteFactorPtr.tid()<<endl;
#endif
    DblNumMat RemoteFactor(A.n-j,A.blksize);
    TIMER_START(Fetching_factor);
    upcxx::copy(remoteFactorPtr,RemoteFactor.GData(),RemoteFactor.Size());
    TIMER_STOP(Fetching_factor);

#ifdef _DEBUG_    
    //    logfileptr->OFS()<<"Fetched factor from P"<<remoteFactorPtr.tid()<<endl;
#endif

    TIMER_START(Update_async_loop);
    //do all my updates with i>j
    for(Int i = j+A.blksize; i<A.n;i+=A.blksize){
      if(A.iam==A.MAP(i,j)){
#ifdef _DEBUG_
        //        logfileptr->OFS()<<"Updating column "<<i<<" with columns from P"<<remoteFactorPtr.tid()<<endl;
#endif
        A.Update(j, i, RemoteFactor);



        TIMER_START(Launch_Aggregates);
        //is it my last update to column i ?
        bool last_upd = true;
        for(Int jj = j+A.blksize; jj<i;jj+=A.blksize){
          if(A.iam==A.MAP(i,jj)){
            //logfileptr->OFS()<<j<<" I also have to update column "<<jj<<endl;
            last_upd=false;
            break;
          }
        }
        if(last_upd){
          Int target = A.MAP(i,i);
#ifdef _DEBUG_
          logfileptr->OFS()<<"Updating from "<<j<<", launching aggregate on P"<<target<<" for column "<<i<<endl;
#endif
          global_ptr<double> WbufPtr(&A.WLower[i/A.blksize](0,0));
          upcxx::async(target)(Aggregate_Async,A.RemoteObjPtrs[target],i,WbufPtr);
        }
        TIMER_STOP(Launch_Aggregates);


      }
    }
    TIMER_STOP(Update_async_loop);

#ifdef _DEBUG_
    logfileptr->OFS()<<"Quitting Update_Async"<<endl<<endl;
#endif
    TIMER_STOP(Update_Async);
  }









}

