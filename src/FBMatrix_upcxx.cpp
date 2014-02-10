/// @file FBMatrix.cpp
/// @brief Matrix structure and methods for the FAN-BOTH method.
/// @author Mathias Jacquelin
/// @date 2014-01-21

#include "FBMatrix_upcxx.hpp"

#include <upcxx.h>
#include  "NumMat_upcxx.hpp"

#include  "utility.hpp"
#include  "blas.hpp"
#include  "lapack.hpp"


#include  "LogFile.hpp"
#include  "upcxx_additions.hpp"

#ifndef ADVANCE_COMM
#define ADVANCE_COMM 1
#endif




//#define _DEBUG_



namespace LIBCHOLESKY{




  void Gathercopy(upcxx::global_ptr<FBMatrix_upcxx> Objptr, global_ptr<double> Adest, Int j);

  void Update_Async(upcxx::global_ptr<FBMatrix_upcxx> Aptr, Int j, upcxx::global_ptr<double> remoteFactorPtr);
  void Aggregate_Async(upcxx::global_ptr<FBMatrix_upcxx> Aptr, Int j, global_ptr<double> remoteAggregatePtr);

  void Factor_Async(upcxx::global_ptr<FBMatrix_upcxx> Aptr, Int j);





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

  
  FBMatrix_upcxx::FBMatrix_upcxx():prefetch(0),outstdAggreg(0),FBMatrix(){
  }


//  FBMatrix_upcxx::~FBMatrix_upcxx(){
//    for(int i =0; i<AchunkLower.size();++i){
//      if(AchunkLower[i]!=NULL){
//        delete (DblNumMat_upcxx *)(this->AchunkLower[i]);
//      }
//    }
//    for(int i =0; i<WLower.size();++i){
//      if(WLower[i]!=NULL){
//        delete (DblNumMat_upcxx *)(this->WLower[i]);
//      }
//    }
//    AchunkLower.clear();
//    WLower.clear();
//  }



  void FBMatrix_upcxx::Initialize(upcxx::shared_array<upcxx::global_ptr<FBMatrix_upcxx> > * RemoteObjPtr){
    np=THREADS;
    iam=MYTHREAD;
    TIMER_START(RemotePtr_fetch);
    RemoteObjPtrs.resize(RemoteObjPtr->size());
    for(int i =0; i< RemoteObjPtr->size();i++){
      RemoteObjPtrs[i]=(*RemoteObjPtr)[i];
      logfileptr->OFS()<<RemoteObjPtrs[i]<<endl;
    }
    TIMER_STOP(RemotePtr_fetch);
  }


  void FBMatrix_upcxx::Distribute( DblNumMat & Aorig){
    //    SetValue(Achunk,(double)iam);
    for(Int j = 0; j<n;j+=blksize){ 
      Int local_j = global_col_to_local(j);
      Int jb = min(blksize, n-j);
      if(iam==MAP(j,j)){
        global_ptr<double> Aorigptr(&Aorig(j,j));
        global_ptr<double> Adestptr(&(*AchunkLower[local_j/blksize])(0,0));
#ifdef ASYNC_COPY
        upcxx::ldacopy_async(n-j,jb,Aorigptr,n,Adestptr,n-j);
#else
//        logfileptr->OFS()<<"Copying from "<<Aorigptr<<" to "<<Adestptr<<endl;
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






 void FBMatrix_upcxx::Allocate(Int np,Int pn, Int pblksize){
    FBMatrix::Allocate(np,pn,pblksize);
    for(Int j = 0; j<n;j+=blksize){ 
      Int local_j = global_col_to_local(j);
      Int jb = min(blksize, n-j);
      if(iam==MAP(j,j)){

        logfileptr->OFS()<<"Resizing factor buffer for column :"<<j<<endl;

        AchunkLower[local_j/blksize] = (DblNumMat*)new DblNumMat_upcxx(n-j,jb);
      }

      logfileptr->OFS()<<"Resizing update buffer for column :"<<j<<endl;
      WLower[j/blksize] = (DblNumMat*)new DblNumMat_upcxx(n-j,jb);
      SetValue(*WLower[j/blksize],0.0);
    }
    logfileptr->OFS()<<"Done"<<endl;
  }







  void FBMatrix_upcxx::Gather( DblNumMat_upcxx & Adest){
      Adest.Resize(n,n);
      SetValue(Adest,0.0);
    
//      logfileptr->OFS()<<Adest.m()<<" "<<Adest.n()<<std::endl;

    for(Int j = 0; j<n;j+=blksize){ 
      Int local_j = global_col_to_local(j);
      Int jb = min(blksize, n-j);
      Int target = MAP(j,j);
      global_ptr<double> Adestptr(&Adest(j,j));
      upcxx::async(target)(Gathercopy,RemoteObjPtrs[target],Adestptr,j);
    }
    upcxx::wait();
  }




  void FBMatrix_upcxx::WaitFactorization(){     

    //while(!fact_finished){  upcxx::progress(); };
    while(!fact_finished){  upcxx::advance(1, 99); };
    if(MYTHREAD==MAP(n-1,n-1)){
      TIMER_STOP(SYNC_END);
    }
  }

  void FBMatrix_upcxx::NumericalFactorization(){
    if(iam==MAP(0,0)){
      upcxx::async(iam)(Factor_Async,RemoteObjPtrs[iam],0);
    }
  }

  void Gathercopy(upcxx::global_ptr<FBMatrix_upcxx> Objptr, global_ptr<double> Adest, Int j){
    assert(Objptr.tid()==MYTHREAD);
    FBMatrix_upcxx & A = *Objptr;
    Int local_j = A.global_col_to_local(j);
    Int jb = min(A.blksize, A.n-j);
    DblNumMat & LocalChunk = *A.AchunkLower[local_j/A.blksize];
    global_ptr<double> Asrc = &LocalChunk(0,0);
#ifdef ASYNC_COPY
    upcxx::ldacopy_async(A.n-j,jb,Asrc,A.n-j,Adest,A.n);
    upcxx::async_copy_fence();
#else
    upcxx::ldacopy(A.n-j,jb,Asrc,A.n-j,Adest,A.n);
#endif

  }


  void Factor_Async(upcxx::global_ptr<FBMatrix_upcxx> Aptr, Int j){
    TIMER_START(Factor_Async);
#ifdef _ASSERT_
    assert(Aptr.tid()==MYTHREAD);
#endif
    FBMatrix_upcxx & A = *Aptr;

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
          //global_ptr<double> AchkPtr(&(*A.AchunkLower[local_j/A.blksize])(0,0));
          global_ptr<double> AchkPtr =  ((DblNumMat_upcxx *)A.AchunkLower[local_j/A.blksize])->GData();

#ifdef _DEBUG_
          logfileptr->OFS()<<"Launching update from "<<j<<" on P"<<target<<endl;
#endif
          upcxx::async(target)(Update_Async,A.RemoteObjPtrs[target],j,AchkPtr);
//          upcxx::advance(0,ADVANCE_COMM);
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









#ifndef ASYNC_COPIES
  void Update_Async(upcxx::global_ptr<FBMatrix_upcxx> Aptr, Int j, upcxx::global_ptr<double> remoteFactorPtr){

    TIMER_START(Update_Async);
#ifdef _ASSERT_
    assert(Aptr.tid()==MYTHREAD);
#endif
    FBMatrix_upcxx & A = *Aptr;


#ifdef _DEBUG_    
    logfileptr->OFS()<<"Fetching factor "<<j<<" from P"<<remoteFactorPtr.tid()<<endl;
#endif
    DblNumMat_upcxx RemoteFactor(A.n-j,A.blksize);
    TIMER_START(Fetching_factor);
    upcxx::copy(remoteFactorPtr,RemoteFactor.GData(),RemoteFactor.Size());
    TIMER_STOP(Fetching_factor);

#ifdef _DEBUG_    
        logfileptr->OFS()<<"Fetched factor "<<j<<endl;
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
          global_ptr<double> WbufPtr = ((DblNumMat_upcxx *)A.WLower[i/A.blksize])->GData();
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




  void Aggregate_Async(upcxx::global_ptr<FBMatrix_upcxx> Aptr, Int j, global_ptr<double> remoteAggregatePtr){
    TIMER_START(Aggregate_Async);
#ifdef _ASSERT_
    assert(Aptr.tid()==MYTHREAD);
#endif
    FBMatrix_upcxx & A = *Aptr;
    //fetch data
#ifdef _DEBUG_    
    logfileptr->OFS()<<"Aggregating Fetching data from P"<<remoteAggregatePtr.tid()<<endl;
#endif
    Int jb = min(A.blksize, A.n-j);

    DblNumMat_upcxx RemoteAggregate(A.n-j,jb); 
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







#else



  void Update_Compute_Async(upcxx::global_ptr<FBMatrix_upcxx> Aptr, Int j, DblNumMat * remoteFactorPtr, upcxx::event * async_copy_event){

    TIMER_START(Update_Async_compute);
    FBMatrix_upcxx & A = *Aptr;

    if(async_copy_event!=NULL){
      //sync with the async_copy
      async_copy_event->wait();
      delete async_copy_event;
    }
    TIMER_STOP(Fetching_factor);

    A.outstdUpdate--;



#ifdef _DEBUG_    
        logfileptr->OFS()<<"Fetched factor "<<j<<endl;
#endif

    TIMER_START(Update_async_loop);
    //do all my updates with i>j
    for(Int i = j+A.blksize; i<A.n;i+=A.blksize){
      if(A.iam==A.MAP(i,j)){
#ifdef _DEBUG_
        //        logfileptr->OFS()<<"Updating column "<<i<<" with columns from P"<<remoteFactorPtr.tid()<<endl;
#endif
        A.Update(j, i, *remoteFactorPtr);



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
          global_ptr<double> WbufPtr(&(*A.WLower[i/A.blksize])(0,0));
          upcxx::async(target)(Aggregate_Async,A.RemoteObjPtrs[target],i,WbufPtr);
//          upcxx::advance(0,ADVANCE_COMM);
          
        }
        TIMER_STOP(Launch_Aggregates);


      }
    }
    TIMER_STOP(Update_async_loop);


    //The factor can now be freed
    delete remoteFactorPtr;
    

//#ifdef _DEBUG_
    logfileptr->OFS()<<"Quitting Update_Async_compute"<<endl<<endl;
//#endif
    TIMER_STOP(Update_Async_compute);
  }





  void Update_Async(upcxx::global_ptr<FBMatrix_upcxx> Aptr, Int j, upcxx::global_ptr<double> remoteFactorPtr){

    TIMER_START(Update_Async_fetch);
#ifdef _ASSERT_
    assert(Aptr.tid()==MYTHREAD);
#endif
    FBMatrix_upcxx & A = *Aptr;


#ifdef _DEBUG_    
    logfileptr->OFS()<<"Fetching factor "<<j<<" from P"<<remoteFactorPtr.tid()<<endl;
#endif
    DblNumMat * RemoteFactor = new DblNumMat(A.n-j,A.blksize);

    if(A.outstdUpdate+1<A.prefetch){
      upcxx::event * async_copy_event = new upcxx::event;
      TIMER_START(Fetching_factor);
      upcxx::async_copy(remoteFactorPtr,RemoteFactor->GData(),RemoteFactor->Size(),async_copy_event);

      A.outstdUpdate++;
      //add the function to the async queue
      upcxx::async(A.iam)(Update_Compute_Async,Aptr,j,RemoteFactor, async_copy_event);
    }
    else{
      TIMER_START(Fetching_factor);
      upcxx::copy(remoteFactorPtr,RemoteFactor->GData(),RemoteFactor->Size());

      //call the function inline
      Update_Compute_Async(Aptr, j, RemoteFactor, NULL);
    }


    logfileptr->OFS()<<"Quitting Update_Async_fetch"<<endl<<endl;
    TIMER_STOP(Update_Async_fetch);
}



  void Aggregate_Compute_Async(upcxx::global_ptr<FBMatrix_upcxx> Aptr, Int j, DblNumMat * remoteAggregatePtr, upcxx::event * async_copy_event){

#ifdef _ASSERT_
    assert(Aptr.tid()==MYTHREAD);
#endif
    FBMatrix_upcxx & A = *Aptr;

    //sync with the async_copy
    async_copy_event->wait();
    delete async_copy_event;
    TIMER_STOP(Fetching_Aggregate);

    A.outstdAggreg--;
#ifdef _DEBUG_    
    logfileptr->OFS()<<"Aggregating update to "<<j<<endl;
#endif
    A.Aggregate(j, *remoteAggregatePtr);
    //The aggregate can now be freed
    delete remoteAggregatePtr;

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
//    TIMER_STOP(Aggregate_Async);
  }

  void Aggregate_Async(upcxx::global_ptr<FBMatrix_upcxx> Aptr, Int j, global_ptr<double> remoteAggregatePtr){
//    TIMER_START(Aggregate_Async);
#ifdef _ASSERT_
    assert(Aptr.tid()==MYTHREAD);
#endif
    FBMatrix_upcxx & A = *Aptr;
    //fetch data
#ifdef _DEBUG_    
    logfileptr->OFS()<<"Aggregating Fetching data from P"<<remoteAggregatePtr.tid()<<endl;
#endif
    Int jb = min(A.blksize, A.n-j);

    DblNumMat * RemoteAggregate =  new DblNumMat(A.n-j,jb);

    upcxx::event * async_copy_event = new upcxx::event;
    TIMER_START(Fetching_Aggregate);
    upcxx::async_copy(remoteAggregatePtr,RemoteAggregate->GData(),RemoteAggregate->Size(),async_copy_event);
    
    A.outstdAggreg++;
    if(A.outstdAggreg<A.prefetch){
      //add the function to the async queue
      upcxx::async(A.iam)(Aggregate_Compute_Async,Aptr,j,RemoteAggregate, async_copy_event);
    }
    else{
      //call the function inline
      Aggregate_Compute_Async(Aptr, j, RemoteAggregate, async_copy_event);
    }
}



#endif







}

