/// @file run_sparse_fanboth.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2014-01-23
//#define _DEBUG_


#include <mpi.h>

#include <time.h>
#include <omp.h>

#include  "sympack.hpp"
#include  "sympack/SupernodalMatrix.hpp"

#include  "sympack/CommTypes.hpp"
#include  "sympack/Ordering.hpp"
#include  "sympack/DistSparseMatrixGraph.hpp"

#define SCALAR double
#define INSCALAR double

#undef DEBUG
#define USE_REDUCE

//#include "sympack/MemTrack.h"




using namespace SYMPACK;


//special reduction: first element contains the max local sum
void PtrSum( void *in, void *inout, int *len, MPI_Datatype *dptr ) 
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




int main(int argc, char **argv) 
{
  MPI_Init(&argc,&argv);
  upcxx::init(&argc, &argv);


  MPI_Comm worldcomm;
  MPI_Comm_dup(MPI_COMM_WORLD,&worldcomm);
  //  int np, iam;
  MPI_Comm_size(worldcomm,&np);
  MPI_Comm_rank(worldcomm,&iam);

  NGCholOptions optionsFact;

#if defined(PROFILE) || defined(PMPI)
  TAU_PROFILE_INIT(argc, argv);
#endif

  logfileptr = new LogFile(iam);
  logfileptr->OFS()<<"********* LOGFILE OF P"<<iam<<" *********"<<endl;
  logfileptr->OFS()<<"**********************************"<<endl;

  // *********************************************************************
  // Input parameter
  // *********************************************************************
  std::map<std::string,SYMPACK::vector<std::string> > options;

  OptionsCreate(argc, argv, options);

  std::string filename;
  if( options.find("-in") != options.end() ){
    filename= options["-in"].front();
  }
  std::string informatstr;
  if( options.find("-inf") != options.end() ){
    informatstr= options["-inf"].front();
  }

  if( options.find("-ordering") != options.end() ){
    if(options["-ordering"].front()=="AMD"){
      optionsFact.ordering = AMD;
    }
    else if(options["-ordering"].front()=="NDBOX"){
      optionsFact.ordering = NDBOX;
    }
    else if(options["-ordering"].front()=="NDGRID"){
      optionsFact.ordering = NDGRID;
    }
#ifdef USE_SCOTCH
    else if(options["-ordering"].front()=="SCOTCH"){
      optionsFact.ordering = SCOTCH;
    }
#endif
#ifdef USE_METIS
    else if(options["-ordering"].front()=="METIS"){
      optionsFact.ordering = METIS;
    }
#endif
#ifdef USE_PARMETIS
    else if(options["-ordering"].front()=="PARMETIS"){
      optionsFact.ordering = PARMETIS;
    }
#endif
#ifdef USE_PTSCOTCH
    else if(options["-ordering"].front()=="PTSCOTCH"){
      optionsFact.ordering = PTSCOTCH;
    }
#endif
    else{
      optionsFact.ordering = MMD;
    }
  }

  Int all_np = np;

  MPI_Comm workcomm;
  MPI_Comm_split(worldcomm,iam<np,iam,&workcomm);
  Int new_rank = (iam<np)?iam:iam-np;
  upcxx::team_all.split(iam<np,new_rank, workteam);
  
  if(iam<np){
    //  int np, iam;
    MPI_Comm_size(workcomm,&np);
    MPI_Comm_rank(workcomm,&iam);


    DistSparseMatrix<SCALAR> HMat(workcomm);
    ReadMatrix<SCALAR,INSCALAR>(filename , informatstr,  HMat);
    if(iam==0){ cout<<"Matrix order is "<<HMat.size<<endl; }

    {    

      optionsFact.commEnv = new CommEnvironment(workcomm);


      SparseMatrixStructure Local = HMat.GetLocalStructure();

      DistSparseMatrixGraph graph;
      graph.comm = HMat.comm;
      graph.baseval = 0;
      graph.keepDiag = 0;
      graph.sorted = 1;
      graph.FromStructure(Local);

      //for testing only
      //SparseMatrixStructure Global;
      //Local.ToGlobal(Global,optionsFact.commEnv->MPI_GetComm());
      //Global.ExpandSymmetric();

    Ordering Order_;
    Order_.SetCommEnvironment(optionsFact.commEnv);

#if 0
      SparseMatrixStructure Global;
      Local.ToGlobal(Global,optionsFact.commEnv->MPI_GetComm());
      Global.ExpandSymmetric();
    //Create an Ordering object to hold the permutation
    Ordering Order_;
    Order_.SetCommEnvironment(optionsFact.commEnv);
    Order_.SetStructure(Global);
    logfileptr->OFS()<<"Structure set"<<endl;

    TIMER_START(ORDERING);
    switch(optionsFact.ordering){
      case MMD:
        //Reoder the matrix with MMD
        Order_.MMD();
        break;
      case AMD:
        Order_.AMD();
        break;
      case NDBOX:
        Order_.NDBOX();
        break;
      case NDGRID:
        Order_.NDGRID();
        break;
#ifdef USE_SCOTCH
      case SCOTCH:
        Order_.SCOTCH();
        break;
#endif
#ifdef USE_METIS
      case METIS:
        Order_.METIS();
        break;
#endif
#ifdef USE_PARMETIS
      case PARMETIS:
        Order_.PARMETIS();
        break;
#endif
#ifdef USE_PTSCOTCH
      case PTSCOTCH:
        Order_.PTSCOTCH();
        break;
#endif
      default:
        //do nothing: either natural or user provided ordering
        break;
    }
    TIMER_STOP(ORDERING);
    if(iam==0){cout<<"Ordering done"<<endl;}
    logfileptr->OFS()<<"Ordering done"<<endl;

    //logfileptr->OFS()<<"Matrix reexpanding"<<endl;
    //Global_->ExpandSymmetric();
    logfileptr->OFS()<<"Matrix expanded"<<endl;








      delete optionsFact.commEnv;
#else

    double tstart = MPI_Wtime();
    graph.ExpandSymmetric();
    double tstop = MPI_Wtime();
    if(iam==0){
      cout<<"Time elapsed (Expanding): "<<tstop-tstart<<endl;
    }

    tstart = MPI_Wtime();
    Order_.PTSCOTCH(graph);
    tstop = MPI_Wtime();
    if(iam==0){
      cout<<"Time elapsed (Ordering): "<<tstop-tstart<<endl;
    }
    logfileptr->OFS()<<"perm: "<<Order_.perm<<endl;
    logfileptr->OFS()<<"invp: "<<Order_.invp<<endl;


    
    //logfileptr->OFS()<<"graph.colptr: "<<graph.colptr<<endl;
    //logfileptr->OFS()<<"graph.rowind: "<<graph.rowind<<endl;
    //graph.Permute(&Order_.perm[0],&Order_.invp[0]);

    //logfileptr->OFS()<<"graph.colptr (permuted): "<<graph.colptr<<endl;
    //logfileptr->OFS()<<"graph.rowind (permuted): "<<graph.rowind<<endl;


    {
      ETree etree2;

      SparseMatrixStructure Global;
      Local.ToGlobal(Global,optionsFact.commEnv->MPI_GetComm());
      Global.ExpandSymmetric();
      tstart = MPI_Wtime();
      etree2.ConstructETree(Global, Order_);
      tstop = MPI_Wtime();
      if(iam==0){
        cout<<"Time elapsed (ETree2): "<<tstop-tstart<<endl;
      }
      logfileptr->OFS()<<"ETREE2: "<<etree2<<endl;

    }

    tstart = MPI_Wtime();
    ETree etree;
    etree.ConstructETree(graph, Order_);
    
    tstop = MPI_Wtime();
    if(iam==0){
      cout<<"Time elapsed (ETree): "<<tstop-tstart<<endl;
    }
    logfileptr->OFS()<<"ETREE: "<<etree<<endl;

#endif
    }


  }

  MPI_Barrier(workcomm);
  MPI_Comm_free(&workcomm);

  MPI_Barrier(worldcomm);
  MPI_Comm_free(&worldcomm);


  //MemTrack::TrackDumpBlocks(logfileptr->OFS());
  //MemTrack::TrackListMemoryUsage(logfileptr->OFS());
  delete logfileptr;

  upcxx::finalize();
  //int finalized = 0;
  //MPI_Finalized(&finalized);
  //if(!finalized){ 
  //  MPI_Finalize();
  //}
  return 0;
}


