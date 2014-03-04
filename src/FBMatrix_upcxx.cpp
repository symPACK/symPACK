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

#include "timer.hpp"


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

    upcxx::barrier();
  }

  void signal_exit()
  {
    for (int i = 0; i < THREADS; i++) {
      async(i)(signal_exit_am);
    }
  }

  
  FBMatrix_upcxx::FBMatrix_upcxx():prefetch(0),outstdAggreg(0),outstdUpdate(0),FBMatrix(){
  }



  void FBMatrix_upcxx::Initialize(upcxx::shared_array<upcxx::global_ptr<FBMatrix_upcxx> > * RemoteObjPtr){
    aggregate_comm_time=0.0;
    factor_comm_time=0.0;
    np=THREADS;
    iam=MYTHREAD;
    TIMER_START(RemotePtr_fetch);
    RemoteObjPtrs.resize(RemoteObjPtr->size());
    for(int i =0; i< RemoteObjPtr->size();i++){
      RemoteObjPtrs[i]=(*RemoteObjPtr)[i];
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

//        logfileptr->OFS()<<"Resizing factor buffer for column :"<<j<<endl;

        AchunkLower[local_j/blksize] = (DblNumMat*)new DblNumMat_upcxx(n-j,jb);
      }

//      logfileptr->OFS()<<"Resizing update buffer for column :"<<j<<endl;
      WLower[j/blksize] = (DblNumMat*)new DblNumMat_upcxx(n-j,jb);
      SetValue(*WLower[j/blksize],0.0);
    }
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

//    while(!fact_finished){  upcxx::progress(); };
    int loc_finished = fact_finished;
    double time_sta_overall, time_end_overall, time_overall =0.0;
    double time_sta_adv, time_end_adv, time_adv =0.0;
    double time_sta_signal, time_end_signal, time_signal =0.0;

      time_sta_overall = get_time();
//    Int cnt = 0;
    while( !fact_finished ){
//      cnt++;
//      int i= 42;
//      time_sta_adv = get_time();
//      TIMER_START(MAIN_ADVANCE);
      upcxx::advance(); 
//      TIMER_STOP(MAIN_ADVANCE);
//      time_end_adv = get_time();
//      time_adv += time_end_adv - time_sta_adv;
//      time_sta_signal = get_time();
//      TIMER_START(WHILE_TEST);
//      loc_finished =fact_finished;
//      TIMER_STOP(WHILE_TEST)
//      time_end_signal = get_time();
//      time_signal += time_end_signal - time_sta_signal;
    }

    time_end_overall = get_time();
    time_overall += time_end_overall - time_sta_overall;

//    Int cnt2 = cnt;
//    while( cnt-->0 ){
//      int i= 42;
//      time_sta_signal = get_time();
//      TIMER_START(WHILE_TEST);
//      loc_finished =i;
//      TIMER_STOP(WHILE_TEST)
//      time_end_signal = get_time();
//      time_signal += time_end_signal - time_sta_signal;
//    }

    
//    logfileptr->OFS()<<"Loops: "<<cnt2<<endl; 
    logfileptr->OFS()<<"FANBOTH_WAITFACT: "<<time_overall<<endl; 
    logfileptr->OFS()<<"MAIN_ADVANCE: "<<time_adv<<endl; 
    logfileptr->OFS()<<"WHILE_TEST: "<<time_signal<<endl; 

//    while(!fact_finished){  upcxx::advance(); };
//    if(MYTHREAD==MAP(n-1,n-1)){
//      TIMER_STOP(SYNC_END);
//    }

  }

  void FBMatrix_upcxx::NumericalFactorization(){
    if(iam==MAP(0,0)){
      upcxx::async(iam)(Factor_Async,RemoteObjPtrs[iam],0);
    }
  }


typedef upcxx::global_ptr<double> dgptr;

bool isNonNull (dgptr &  i) {
  return ( static_cast<double *>(i)!=dgptr(NULL));
}



void FBMatrix_upcxx::NumericalFactorizationLoop(){
  Int numStep = ceil((double)n/(double)blksize);

  DblNumMat_upcxx RemoteAggregate(n,blksize);
  DblNumMat_upcxx RemoteFactorBuf(n,blksize);


  //    upcxx::shared_array< dgptr, np > ArrAggregateRdy2;
  //    ArrAggregateRdy2.init(np*np);
  //    for(int i=0;i<np;i++){
  //       ArrAggregateRdy2[iam*np+i]=dgptr(NULL);
  //    }

  upcxx::shared_array< dgptr, 1 > ArrFactorRdy2;
  ArrFactorRdy2.init(np);
  ArrFactorRdy2[iam]=dgptr(NULL);

  NumMat_upcxx< dgptr > RemoteAggregateRdy(np,1);
  SetValue(RemoteAggregateRdy, dgptr(NULL));
  upcxx::shared_array< upcxx::global_ptr< dgptr > > ArrAggregateRdy;
  ArrAggregateRdy.init(np);
  ArrAggregateRdy[iam] = RemoteAggregateRdy.GData();

  upcxx::barrier();
  upcxx::wait();



  for(int js=0; js<numStep;js++){
    int j = min(js*blksize,n-1);
    Int jb = min(blksize, n-j);
    Int local_j = global_col_to_local(j);

    //do the aggregation and factorization
    if(iam==MAP(j,j) ){

      //aggregate previous updates
      if(AggLock[j/blksize]>0){
        dgptr * last = RemoteAggregateRdy.Data()+RemoteAggregateRdy.Size();
        dgptr * pos = last;
        while(AggLock[j/blksize]>0){
          //poll for available data
          do{ /*upcxx::wait();*/  pos =std::find_if(RemoteAggregateRdy.Data(),last,isNonNull); }while(pos==last);
#ifdef _DEBUG_
          logfileptr->OFS()<<"Receiving aggregate"<<endl;
#endif
          TIMER_START(Aggregate_Recv);
          upcxx::copy<double>(*pos,RemoteAggregate.GData(),(n-j)*jb);
          TIMER_STOP(Aggregate_Recv);
          *pos = dgptr(NULL);


          /************** Code for shared_array with block size np ****************/
          ///              int pos=-1;
          ///              bool found = false;
          ///              do{
          ///                  //upcxx::wait();
          ///                  for(int i=0;i<np;i++){
          ///                    dgptr cur = (dgptr)ArrAggregateRdy2[iam*np+i];
          ///                    if((double*)cur!=NULL){
          ///                      pos = iam*np+i;
          ///                      found = true;
          ///                      break;
          ///                    }
          ///                  }
          ///              }while(! found);
          ///
          ///#ifdef _DEBUG_
          ///          logfileptr->OFS()<<"Receiving aggregate"<<endl;
          ///#endif
          ///
          ///              TIMER_START(Aggregate_Recv);
          ///              upcxx::copy<double>(ArrAggregateRdy2[pos],RemoteAggregate.GData(),(n-j)*jb);
          ///              TIMER_STOP(Aggregate_Recv);
          ///              ArrAggregateRdy2[pos] = dgptr(NULL);
          /********** end of code for shared_array with block size np *************/

#ifdef _DEBUG_
          logfileptr->OFS()<<"Aggregating to column "<<j<<endl;
#endif
          TIMER_START(Remote_Aggregate);
          Aggregate(j, RemoteAggregate);
          TIMER_STOP(Remote_Aggregate);
          --AggLock[j/blksize];
        }
      }


#ifdef _DEBUG_
      logfileptr->OFS()<<"Factoring column "<<j<<endl;
#endif
      //Factor current column 
      Factor(j);

      //send factor
      TIMER_START(Send_Factor);

      for(Int target = 0; target<np;target++){
        for(Int i = j+blksize; i<n; i+=blksize){
          if(MAP(i,j)==target){
            DblNumMat_upcxx & Achk = *(DblNumMat_upcxx*)AchunkLower[local_j/blksize];

#ifdef _DEBUG_
            logfileptr->OFS()<<"Launching update from "<<j<<" on P"<<target<<endl;
#endif

            //push the pointer to the factor on P(i,j)
            dgptr dest = ArrFactorRdy2[target];
            dgptr src(Achk.GData());
            ArrFactorRdy2[target] = src;
            break;
          }
        }
      }
      TIMER_STOP(Send_Factor);
    }

    //do the updates

    bool have_factor = false;
    //foreach processor, do I have an update from column j to any column ?
    for(int i=j+jb;i<n ;i+=jb){
      //compute the update
      if(iam==MAP(i,j)){
        DblNumMat_upcxx * RemoteFactorPtr;



        RemoteFactorPtr = & RemoteFactorBuf;
        if(!have_factor){
          //poll for available data
          while((double*)ArrFactorRdy2[iam]==NULL);

#ifdef _DEBUG_
          logfileptr->OFS()<<"Receiving factor "<<j<<endl;
#endif

          //Receive factor
          TIMER_START(Factor_Recv);
          upcxx::copy<double>(ArrFactorRdy2[iam],RemoteFactorPtr->GData(),(n-j)*jb);
          TIMER_STOP(Factor_Recv);
          ArrFactorRdy2[iam] = dgptr(NULL);
          have_factor=true;
        }

#ifdef _DEBUG_
        logfileptr->OFS()<<"Updating "<<i<<" with "<<j<<endl;
#endif
        //Do the update
        Update(j, i, *RemoteFactorPtr);

        //logfileptr->OFS()<<"Done Updating "<<i<<" with "<<j<<endl;
        TIMER_START(Launch_Aggregates);
        //if this is the last update to column i, send aggregate
        if(lastUpdate(j,i)){

          Int target = MAP(i,i);
#ifdef _DEBUG_
          logfileptr->OFS()<<"I did my last update to "<<i<<endl;
#endif

          if(1){
#ifdef _DEBUG_
            logfileptr->OFS()<<"Updating from "<<j<<", sending aggregate to P"<<target<<" for column "<<i<<endl;
#endif

            //push the pointer to the aggregate vector on P(i,i)
            global_ptr< dgptr > dest = (global_ptr< dgptr >)ArrAggregateRdy[target]+ iam;
            global_ptr< dgptr > src( &(((DblNumMat_upcxx *)WLower[i/blksize])->GData()) );
            upcxx::copy< dgptr >(src,dest,1);

            /************** Code for shared_array with block size np ****************/
            ///                dgptr src = (((DblNumMat_upcxx *)WLower[i/blksize])->GData());
            ///                ArrAggregateRdy2[target*np+iam] = src;
            ///                upcxx::wait();
            /********** end of code for shared_array with block size np *************/

          }
        }
        TIMER_STOP(Launch_Aggregates);
      }
    }
  } 

  fact_finished = 1;
  upcxx::barrier();

}

  //Check if j is the last column updating i
  bool FBMatrix_upcxx::lastUpdate(Int j, Int i){
    //is it my last update to column i ?
    bool last_upd = true;
    for(Int jj = j+blksize; jj<i;jj+=blksize){
      //logfileptr->OFS()<<"Next update to "<<i<<" looking at "<<jj<<endl;
      if(iam==MAP(i,jj)){
        //logfileptr->OFS()<<jj<<" updates "<<i<<endl;
        last_upd=false;
        break;
      }
    }
    return last_upd;
  }






//  void FBMatrix_upcxx::NumericalFactorizationEvent(){
//      Int numStep = ceil((double)n/(double)blksize);
//
//
//      upcxx::event * event_factos[numStep];
//      upcxx::event * event_aggregates[numStep];
//      upcxx::event * event_updates[numStep];
//
//      for(int i =0; i<numStep;i++){
//        event_factos[i] = new upcxx::event();
//        event_aggregates[i] = new upcxx::event();
//        event_updates[i] = new upcxx::event();
//      }
//
//      for(int js=0; js<numStep;js++){
//        int j = min(js*blksize,n-1);
//        Int jb = min(blksize, n-j);
//        Int local_j = global_col_to_local(j);
//
//          if(js==0){ 
//            upcxx::Async<MAP(j,j),&event_factos[js]>(Facto_Async2,,,)
//          }
//          else{
//            upcxx::Async_after<MAP(j,j),event_updates[js],&event_fggregates[js]>(Aggregate_Async2,,,)
//            upcxx::Async_after<MAP(j,j),event_aggregates[js],&event_factos[js]>(Facto_Async2,,,)
//          }
//
//
//          for(Int target = 0; target<np;target++){
//            for(Int i = j+blksize; i<n; i+=blksize){
//               Int is = i/blksize;
//                //is it possible to have JOIN ?
//               upcxx::Async_after<MAP(i,j),event_factos[js],&event_update[is]>(Update_Async2,,,)
//            }
//          }
//      }
//  }









  void Gathercopy(upcxx::global_ptr<FBMatrix_upcxx> Objptr, global_ptr<double> Adest, Int j){
    assert(Objptr.tid()==MYTHREAD);
    FBMatrix_upcxx & A = *Objptr;
    Int local_j = A.global_col_to_local(j);
    Int jb = min(A.blksize, A.n-j);
    DblNumMat_upcxx & LocalChunk = *(DblNumMat_upcxx*)A.AchunkLower[local_j/A.blksize];
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
    Int async_created = 0;
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
          if(target!=A.iam){
            async_created++;
          }
//int num_out = upcxx::advance_out_task_queue(out_task_queue, ADVANCE_COMM);
          break;
        }
      }
    }
    TIMER_STOP(Launch_Updates);

    TIMER_START(ADVANCE_UPDATES);
    upcxx::advance(0,async_created);
    TIMER_STOP(ADVANCE_UPDATES);

    if(j+A.blksize>=A.n){
#ifdef _DEBUG_
      logfileptr->OFS()<<"Unlocking quit"<<endl;
#endif

//      TIMER_START(SYNC_END);
      signal_exit();

    }

#ifdef _DEBUG_
    logfileptr->OFS()<<"Quitting Factor_Async"<<endl;
    logfileptr->OFS()<<"********************************************"<<endl;
#endif
    TIMER_STOP(Factor_Async);
  }









#ifndef ASYNC_UPDATES
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
        if(A.lastUpdate(j,i)){
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









#else



  void Update_Compute_Async(upcxx::global_ptr<FBMatrix_upcxx> Aptr, Int j, DblNumMat_upcxx * remoteFactorPtr, upcxx::event * async_copy_event, double tstart){

#ifdef _DEBUG_
    logfileptr->OFS()<<"Entering Update_Compute_Async"<<endl<<endl;
#endif
    TIMER_START(Update_Async_compute);
    FBMatrix_upcxx & A = *Aptr;

//#ifdef NO_ASYNC_COPY
      //upcxx::async_copy_fence();
//#else
    if(async_copy_event!=NULL){
#ifdef ASYNC_COPY_FENCE
      upcxx::async_copy_fence();
#else
      //sync with the async_copy
#ifndef CALLBACK
      async_copy_event->wait();
#endif

      assert(async_copy_event->isdone());

#endif 
      A.outstdUpdate--;
#ifdef _DEBUG_PREFETCH_    
      logfileptr->OFS()<<"DONE Outstanding updates: "<<A.outstdUpdate<<"/"<<A.prefetch<<endl;
#endif
    }
    TIMER_STOP(Fetching_Factor);
//#endif

    double tstop = get_time();
    A.factor_comm_time += tstop - tstart;


#ifdef _DEBUG_    
        logfileptr->OFS()<<"Fetched factor "<<j<<endl;
#endif

    TIMER_START(Update_async_loop);
    //do all my updates with i>j
    Int async_created = 0;
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
          
          if(target!=A.iam){
            async_created++;
          }
//          int num_out = upcxx::advance_out_task_queue(out_task_queue, ADVANCE_COMM);
          
        }
        TIMER_STOP(Launch_Aggregates);


      }
    }
    TIMER_STOP(Update_async_loop);

    TIMER_START(ADVANCE_AGGREGATES);
    upcxx::advance(0,async_created);
    TIMER_STOP(ADVANCE_AGGREGATES);

    //The factor can now be freed
    delete remoteFactorPtr;
    

//#ifdef _DEBUG_
//    logfileptr->OFS()<<"Quitting Update_Async_compute"<<endl<<endl;
//#endif

#ifdef _DEBUG_
    logfileptr->OFS()<<"Quitting Update_Compute_Async"<<endl<<endl;
#endif
    TIMER_STOP(Update_Async_compute);
  }





  void Update_Async(upcxx::global_ptr<FBMatrix_upcxx> Aptr, Int j, upcxx::global_ptr<double> remoteFactorPtr){

#ifdef _DEBUG_
    logfileptr->OFS()<<"Entering Update_Async (fetch)"<<endl<<endl;
#endif
    TIMER_START(Update_Async_fetch);
#ifdef _ASSERT_
    assert(Aptr.tid()==MYTHREAD);
#endif
    FBMatrix_upcxx & A = *Aptr;


#ifdef _DEBUG_    
    logfileptr->OFS()<<"Fetching factor "<<j<<" from P"<<remoteFactorPtr.tid()<<endl;
#endif
    DblNumMat_upcxx * RemoteFactor = new DblNumMat_upcxx(A.n-j,A.blksize);
//    logfileptr->OFS()<<"Done"<<endl;

    if(A.outstdUpdate+1<=A.prefetch){

    double tstart = get_time();
      TIMER_START(Fetching_Factor);
#ifdef NO_ASYNC_COPY
      upcxx::event * async_copy_event = NULL;
      upcxx::copy(remoteFactorPtr,RemoteFactor->GData(),RemoteFactor->Size());
#else
      upcxx::event * async_copy_event = new upcxx::event;
#ifdef ASYNC_COPY_FENCE
      upcxx::async_copy(remoteFactorPtr,RemoteFactor->GData(),RemoteFactor->Size());
#else

#ifdef CALLBACK
#ifdef _DEBUG_PREFETCH_    
      logfileptr->OFS()<<"USING ASYNC_AFTER: "<<A.outstdUpdate+1<<"/"<<A.prefetch<<" async copy"<<endl;
#endif
      upcxx::async_after(A.iam,async_copy_event)(Update_Compute_Async,Aptr,j,RemoteFactor, async_copy_event,tstart);
      upcxx::async_copy(remoteFactorPtr,RemoteFactor->GData(),RemoteFactor->Size(),async_copy_event);

#else
#ifdef _DEBUG_PREFETCH_    
      logfileptr->OFS()<<"Outstanding updates: "<<A.outstdUpdate+1<<"/"<<A.prefetch<<" async copy"<<endl;
#endif
      upcxx::async_copy(remoteFactorPtr,RemoteFactor->GData(),RemoteFactor->Size(),async_copy_event);
#endif



#endif
#endif

      A.outstdUpdate++;
      //add the function to the async queue
#ifndef CALLBACK
      upcxx::async(A.iam)(Update_Compute_Async,Aptr,j,RemoteFactor, async_copy_event,tstart);
#endif
    }
    else{
#ifdef _DEBUG_PREFETCH_    
      logfileptr->OFS()<<"Outstanding updates: "<<A.outstdUpdate<<"/"<<A.prefetch<<" too many > blocking copy"<<endl;
#endif
      double tstart = get_time();
      TIMER_START(Fetching_Factor);
      upcxx::copy(remoteFactorPtr,RemoteFactor->GData(),RemoteFactor->Size());

      //call the function inline
      Update_Compute_Async(Aptr, j, RemoteFactor, NULL,tstart);
    }


//    logfileptr->OFS()<<"Quitting Update_Async_fetch"<<endl<<endl;
    TIMER_STOP(Update_Async_fetch);
#ifdef _DEBUG_
    logfileptr->OFS()<<"Quitting Update_Async (fetch)"<<endl<<endl;
#endif
}

#endif


#ifndef ASYNC_AGGREGATES


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

  void Aggregate_Compute_Async(upcxx::global_ptr<FBMatrix_upcxx> Aptr, Int j, DblNumMat_upcxx * remoteAggregatePtr, upcxx::event * async_copy_event, double tstart){

#ifdef _DEBUG_
    logfileptr->OFS()<<"Entering Aggregate_Compute_Async"<<endl<<endl;
#endif
    TIMER_START(Aggregate_Async_compute);

#ifdef _ASSERT_
    assert(Aptr.tid()==MYTHREAD);
#endif
    FBMatrix_upcxx & A = *Aptr;

//#ifdef NO_ASYNC_COPY
      //upcxx::async_copy_fence();
//#else
    if(async_copy_event!=NULL){
#ifdef ASYNC_COPY_FENCE
      upcxx::async_copy_fence();
#else
      //sync with the async_copy
#ifndef CALLBACK
      async_copy_event->wait();
#endif
      assert(async_copy_event->isdone());

#endif      
      delete async_copy_event;
      A.outstdAggreg--;
#ifdef _DEBUG_PREFETCH_    
      logfileptr->OFS()<<"DONE Outstanding aggregates: "<<A.outstdAggreg<<"/"<<A.prefetch<<endl;
#endif
    }
//#endif
    TIMER_STOP(Fetching_Aggregate);
    double tstop = get_time();
    A.aggregate_comm_time += tstop - tstart;


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
    logfileptr->OFS()<<"Quitting Aggregate_Compute_Async"<<endl<<endl;
#endif
//    TIMER_STOP(Aggregate_Async);
    TIMER_STOP(Aggregate_Async_compute);
  }

  void Aggregate_Async(upcxx::global_ptr<FBMatrix_upcxx> Aptr, Int j, global_ptr<double> remoteAggregatePtr){

#ifdef _DEBUG_
    logfileptr->OFS()<<"Entering Aggregate_Async (Fetch)"<<endl<<endl;
#endif
    TIMER_START(Aggregate_Async_fetch);
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

    DblNumMat_upcxx * RemoteAggregate =  new DblNumMat_upcxx(A.n-j,jb);

    

    if(A.outstdAggreg+1<=A.prefetch){
    double tstart = get_time();
    TIMER_START(Fetching_Aggregate);
#ifdef NO_ASYNC_COPY
    upcxx::event * async_copy_event = NULL;
    upcxx::copy(remoteAggregatePtr,RemoteAggregate->GData(),RemoteAggregate->Size());
#else
    upcxx::event * async_copy_event = new upcxx::event;
#ifdef ASYNC_COPY_FENCE
    upcxx::async_copy(remoteAggregatePtr,RemoteAggregate->GData(),RemoteAggregate->Size());
#else

#ifdef CALLBACK
#ifdef _DEBUG_PREFETCH_    
      logfileptr->OFS()<<"USING ASYNC_AFTER: "<<A.outstdAggreg+1<<"/"<<A.prefetch<<" async copy"<<endl;
#endif
  upcxx::async_after(A.iam,async_copy_event)(Aggregate_Compute_Async,Aptr,j,RemoteAggregate, async_copy_event,tstart);
  upcxx::async_copy(remoteAggregatePtr,RemoteAggregate->GData(),RemoteAggregate->Size(),async_copy_event);
#else
#ifdef _DEBUG_PREFETCH_    
      logfileptr->OFS()<<"Outstanding aggregates: "<<A.outstdAggreg+1<<"/"<<A.prefetch<<" async copy"<<endl;
#endif
    upcxx::async_copy(remoteAggregatePtr,RemoteAggregate->GData(),RemoteAggregate->Size(),async_copy_event);
#endif

#endif
#endif

      A.outstdAggreg++;

#ifndef CALLBACK
      //add the function to the async queue
      upcxx::async(A.iam)(Aggregate_Compute_Async,Aptr,j,RemoteAggregate, async_copy_event,tstart);
#endif
    }
    else{
#ifdef _DEBUG_PREFETCH_    
      logfileptr->OFS()<<"Outstanding aggregates: "<<A.outstdAggreg<<"/"<<A.prefetch<<" too many > blocking copy"<<endl;
#endif
    double tstart = get_time();
      TIMER_START(Fetching_factor);
    upcxx::copy(remoteAggregatePtr,RemoteAggregate->GData(),RemoteAggregate->Size());

      //call the function inline
      Aggregate_Compute_Async(Aptr, j, RemoteAggregate, NULL,tstart);
    }

#ifdef _DEBUG_
    logfileptr->OFS()<<"Quitting Aggregate_Async (Fetch)"<<endl<<endl;
#endif
    TIMER_STOP(Aggregate_Async_fetch);
}



#endif







}

