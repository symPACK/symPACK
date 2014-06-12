/// @file FBMatrix.cpp
/// @brief Matrix structure and methods for the FAN-BOTH method.
/// @author Mathias Jacquelin
/// @date 2014-01-21

#include "FBMatrix_mpi.hpp"

#include  "utility.hpp"
#include  "blas.hpp"
#include  "lapack.hpp"

#include  "LogFile.hpp"


#define TAG_FACTOR 0
#define TAG_AGGREGATE 1
#define TAG_AGGREGATE_SEND i
#define TAG_AGGREGATE_RECV j


namespace LIBCHOLESKY{



  MPIGrid::MPIGrid( MPI_Comm Bcomm, int nprow, int npcol )
  {
    Int info;
    MPI_Initialized( &info );
    if( !info ){
      throw std::logic_error( "MPI has not been initialized." );
    }
    MPI_Group  comm_group;
    MPI_Comm_group( Bcomm, &comm_group );
    MPI_Comm_create( Bcomm, comm_group, &comm );

    MPI_Comm_rank( comm, &mpirank );
    MPI_Comm_size( comm, &mpisize );
    if( mpisize != nprow * npcol ){
      throw std::logic_error( "mpisize != nprow * npcol." ); 
    }

    numProcRow = nprow;
    numProcCol = npcol;

    Int myrow = mpirank / npcol;
    Int mycol = mpirank % npcol;

    MPI_Comm_split( comm, myrow, mycol, &rowComm );
    MPI_Comm_split( comm, mycol, myrow, &colComm );

    MPI_Group_free( &comm_group );

    return ;
  } 		// -----  end of method MPIGrid::MPIGrid  ----- 


  MPIGrid::~MPIGrid	(  )
  {
    MPI_Comm_free( &rowComm );
    MPI_Comm_free( &colComm ); 
    MPI_Comm_free( &comm );

    return ;
  } 		// -----  end of method MPIGrid::~MPIGrid  ----- 

  void FBMatrix_mpi::Initialize(MPIGrid & grid ){
    mpigrid = &grid; 
    iam=grid.mpirank;
    np=grid.mpisize;
#ifdef DRAW_GRAPH
    if(iam==0){
      char suffix[50];
      sprintf(suffix,".dot");
      graphfileptr = new LogFile("graph",suffix);
      graphfileptr->OFS()<<"digraph fanboth{"<<endl;
    }
#endif
  }

  void FBMatrix_mpi::Distribute( DblNumMat & Aorig){
    //    SetValue(Achunk,(double)iam);
    for(Int j = 0; j<n;j+=blksize){ 
      Int local_j = global_col_to_local(j);
      Int jb = min(blksize, n-j);
      if(iam==MAP(j,j)){
        double * Aorigptr = &Aorig(j,j);
        double * Adestptr = &(*AchunkLower[local_j/blksize])(0,0);

        lapack::Lacpy( 'N', n-j, jb, Aorigptr, n,	Adestptr, n-j	);
      }
    }
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



 void FBMatrix_mpi::Allocate(Int & np,Int pn, Int pblksize){
    FBMatrix::Allocate(np,pn,pblksize);
    for(Int j = 0; j<n;j+=blksize){ 
      Int local_j = global_col_to_local(j);
      Int jb = min(blksize, n-j);
      if(iam==MAP(j,j)){
#ifdef _DEBUG_
        logfileptr->OFS()<<"Resizing factor buffer for column :"<<j<<endl;
#endif
        AchunkLower[local_j/blksize] = (DblNumMat*)new DblNumMat(n-j,jb);
      }

#ifdef _DEBUG_
      logfileptr->OFS()<<"Resizing update buffer for column :"<<j<<endl;
#endif
      WLower[j/blksize] = (DblNumMat*)new DblNumMat(n-j,jb);
      SetValue(*WLower[j/blksize],0.0);
    }
  }





  void FBMatrix_mpi::Gather( DblNumMat & Adest){
    if(iam==0)
    {
      Adest.Resize(n,n);
      SetValue(Adest,0.0);
    }
    
    DblNumMat recvBuf;
    if(iam==0){
      recvBuf.Resize(n,blksize);
    }

    for(Int j = 0; j<n;j+=blksize){ 
      Int local_j = global_col_to_local(j);
      Int jb = min(blksize, n-j);
      Int target = MAP(j,j);

      if(iam==0){
        if(target==0){
          //do a local copy
#ifdef _DEBUG_
          logfileptr->OFS()<<"Copying local column"<<j<<endl;
#endif
          DblNumMat & LocalSrcChunk = *AchunkLower[local_j/blksize];
          lapack::Lacpy( 'N', n-j, jb, LocalSrcChunk.Data(), n-j,	&Adest(j,j), n	);
        }
        else{
#ifdef _DEBUG_
logfileptr->OFS()<<"Receiving column"<<j<<" from P"<<target<<endl;
#endif
          //recv in the temp buffer
          MPI_Recv(recvBuf.Data(),(n-j)*jb*sizeof(recvBuf(0,0)),MPI_BYTE,target,TAG_FACTOR,mpigrid->comm,MPI_STATUS_IGNORE);
          //do the final lda_cpy to the output matrix
#ifdef _DEBUG_
logfileptr->OFS()<<"Copying remote column"<<j<<endl;
#endif
          lapack::Lacpy( 'N', n-j, jb, recvBuf.Data(), n-j,	&Adest(j,j), n	);
        }
      }
      else{
        if(iam==target){
#ifdef _DEBUG_
logfileptr->OFS()<<"Sending column"<<j<<" to P0"<<endl;
#endif
          DblNumMat & LocalSrcChunk = *AchunkLower[local_j/blksize];
          MPI_Send(LocalSrcChunk.Data(),LocalSrcChunk.ByteSize(),MPI_BYTE,0,TAG_FACTOR,mpigrid->comm);
        }
      }
    }
  }


#ifdef DRAW_GRAPH

  FBMatrix_mpi::~FBMatrix_mpi(){
    if(iam==0){
      graphfileptr->OFS()<<"}"<<endl;
      delete graphfileptr;
    }
  }

  void FBMatrix_mpi::Draw_Graph(){

    if(iam==0){
      Int numStep = ceil((double)n/(double)blksize);

      for(int js=0; js<numStep;js++){
        int j = min(js*blksize,n-1);
        Int jb = min(blksize, n-j);
        Int local_j = global_col_to_local(j);


        for(int curproc = 0; curproc<mpigrid->mpisize;curproc++){
          int iam = (curproc + MAP(j,j)) % mpigrid->mpisize ;
          //do the aggregation and factorization
          if(iam==MAP(j,j) ){

            //aggregate previous updates
            if(j>0){
              graphfileptr->OFS()<<"a_"<<j<<" [label=\"A("<<j<<","<<iam<<")\", style=filled, fillcolor=\""<<(double)iam/(double)np<<" 0.5 0.8\" ];"<<endl;
            }



            for(Int i = j-blksize; i>=0; i-=blksize){
              Int target = MAP(j,i);
              graphfileptr->OFS()<<"u_"<<j<<"_"<<i<<"_"<<target<<" -> a_"<<j<<";"<<endl;
            }






            //Factor current column 
            graphfileptr->OFS()<<"f_"<<j<<" [label=\"F("<<j<<","<<iam<<")\", style=filled, fillcolor=\""<<(double)iam/(double)np<<" 0.5 0.8\" ];"<<endl;
            if(j>0){
              graphfileptr->OFS()<<"a_"<<j<<" -> f_"<<j<<";"<<endl;
            }
            //send factor
            TIMER_START(Send_Factor);
            for(Int target = 0; target<np;target++){
              for(Int i = j+blksize; i<n; i+=blksize){
                if(MAP(i,j)==target ){
                  graphfileptr->OFS()<<"u_"<<i<<"_"<<j<<"_"<<target<<" [label=\"U("<<i<<","<<j<<","<<target<<")\", style=filled, fillcolor=\""<<(double)target/(double)np<<" 0.5 0.8\" ];"<<endl;
                  graphfileptr->OFS()<<"f_"<<j<<" -> u_"<<i<<"_"<<j<<"_"<<target<<";"<<endl;
                  //break;
                }
              }
            }
            TIMER_STOP(Send_Factor);

          }
        } 
      }
    }
  }
#endif



  void FBMatrix_mpi::NumericalFactorization(){
      Int numStep = ceil((double)n/(double)blksize);

      std::vector<MPI_Request> Factor_Send_Requests(np,MPI_REQUEST_NULL);
      std::vector<MPI_Request> Aggregate_Send_Requests(np,MPI_REQUEST_NULL);

      DblNumMat RemoteAggregate(n,blksize);
      DblNumMat RemoteFactorBuf(n,blksize);
      for(int js=0; js<numStep;js++){
        int j = min(js*blksize,n-1);
        Int jb = min(blksize, n-j);
        Int local_j = global_col_to_local(j);

        //do the aggregation and factorization
        if(iam==MAP(j,j) ){

          //aggregate previous updates
          //Do I have local aggregates ?
          for(int i=(js-1)*blksize; i>=max(0,(js-np))*blksize; i-=blksize){

            #ifdef _DEBUG_
                  logfileptr->OFS()<<"Searching for local aggregates from "<<i<<" to "<<j<<endl;
                  logfileptr->OFS()<<i<<" is on P"<<MAP(j,i)<<endl;
            #endif
      
            if(iam==MAP(j,i)){
              //This is a local aggregate
              DblNumMat & LocalAggregate = *WLower[j/blksize];

                #ifdef _DEBUG_
                  logfileptr->OFS()<<"Aggregating local to column "<<j<<endl;
                #endif
              TIMER_START(Local_Aggregate);
              Aggregate(j, LocalAggregate);
              TIMER_STOP(Local_Aggregate);
              --AggLock[j/blksize];
              break;
            }
          }
          //Process remote aggregates
          if(AggLock[j/blksize]>0){
            while(AggLock[j/blksize]>0){

                #ifdef _DEBUG_
                  logfileptr->OFS()<<"Receiving aggregate"<<endl;
               #endif
              TIMER_START(Aggregate_Recv);
              MPI_Recv(RemoteAggregate.Data(),(n-j)*jb*sizeof(double),MPI_BYTE,MPI_ANY_SOURCE,TAG_AGGREGATE_RECV,mpigrid->comm,MPI_STATUS_IGNORE);
              TIMER_STOP(Aggregate_Recv);

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
              if(MAP(i,j)==target && target!=iam){
                DblNumMat & Achk = *AchunkLower[local_j/blksize];

                #ifdef _DEBUG_
                  logfileptr->OFS()<<"Launching update from "<<j<<" on P"<<target<<endl;
                #endif
//                MPI_Send(&Achk(0,0),Achk.ByteSize(),MPI_BYTE,target,TAG_FACTOR, mpigrid->comm );
                MPI_Wait(&Factor_Send_Requests[target],MPI_STATUS_IGNORE);
                MPI_Isend(&Achk(0,0),Achk.ByteSize(),MPI_BYTE,target,TAG_FACTOR, mpigrid->comm,&Factor_Send_Requests[target] );
                break;
              }
            }
          }
          TIMER_STOP(Send_Factor);

        }

        //MPI_Barrier(mpigrid->comm);

        //do the updates

        //DblNumMat RemoteFactorBuf2(n-j,jb);
        bool have_factor = false;

        //foreach processor, do I have an update from column j to any column ?
        for(int i=j+jb;i<n ;i+=jb){
          //compute the update
          if(iam==MAP(i,j)){
            DblNumMat * RemoteFactorPtr;
            if(iam!=MAP(j,j) ){
              RemoteFactorPtr = & RemoteFactorBuf;
              if(!have_factor){
#ifdef _DEBUG_
                logfileptr->OFS()<<"Receiving factor "<<j<<endl;
#endif
                //Receive factor
              TIMER_START(Factor_Recv);
                MPI_Recv(RemoteFactorPtr->Data(),(n-j)*jb*sizeof(double),MPI_BYTE,MAP(j,j),TAG_FACTOR,mpigrid->comm,MPI_STATUS_IGNORE);
              TIMER_STOP(Factor_Recv);
                have_factor=true;
              }
            }
            else{
            #ifdef _DEBUG_
              DblNumMat * LocalFactorPtr=AchunkLower[local_j/blksize];
              logfileptr->OFS()<<"Local factor "<<j<<" ("<<LocalFactorPtr->m()<<","<<LocalFactorPtr->n()<<")"<<endl;
            #endif
              RemoteFactorPtr = AchunkLower[local_j/blksize];
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

              if(target!=iam){
#ifdef _DEBUG_
              logfileptr->OFS()<<"Updating from "<<j<<", sending aggregate to P"<<target<<" for column "<<i<<endl;
#endif
//              MPI_Send(WLower[i/blksize]->Data(),WLower[i/blksize]->ByteSize(),MPI_BYTE,target,TAG_AGGREGATE,mpigrid->comm);
                MPI_Wait(&Aggregate_Send_Requests[target],MPI_STATUS_IGNORE);
                MPI_Isend(WLower[i/blksize]->Data(),WLower[i/blksize]->ByteSize(),MPI_BYTE,target,TAG_AGGREGATE_SEND,mpigrid->comm, &Aggregate_Send_Requests[target]);
              }
            }
            TIMER_STOP(Launch_Aggregates);



          }
        }

        //MPI_Barrier(mpigrid->comm);

      } 
      
      MPI_Barrier(mpigrid->comm);

  }

  //Check if j is the last column updating i
  bool FBMatrix_mpi::lastUpdate(Int j, Int i){
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

}

