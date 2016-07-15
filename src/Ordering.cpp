#include "sympack/Ordering.hpp"
#include "sympack/mmd.hpp"
#include "sympack/utility.hpp"

#include <limits>

#ifdef USE_SCOTCH
#include <typeinfo>
#include <scotch.h>
#endif
#ifdef USE_PTSCOTCH
#include <typeinfo>
#include <ptscotch.h>
#endif

//#if defined(USE_METIS) && not defined(USE_SCOTCH)
#ifdef USE_METIS
#include <typeinfo>
//#define USE_METIS_INC
//#ifdef USE_METIS_INC
#include <metis.h>
//#define METIS_NOPTIONS 40
//#define METIS_NOPTIONS 40
#endif


namespace SYMPACK {

  extern "C" {
    int METIS_NodeND (int * N     , int* XADJ2 , int* ADJ2  , int * VWGT, int* OPTION, int* dback , int* dforw);
    int ParMETIS_V3_NodeND(int * vtxdist  , int* XADJ , int* ADJ  , int * numflag, int* OPTION, int* order , int* sizes, MPI_Comm * comm);
  }


}

namespace SYMPACK{


  void Ordering::SetStructure(SparseMatrixStructure & aGlobal){
    assert(pStructure==NULL);
    pStructure = &aGlobal;
    invp.resize(pStructure->size);
    perm.resize(pStructure->size);
    for(Int i =0;i<invp.size();++i){
      perm[i]=i+1;
      invp[i]=i+1;
    }
  }







  void Ordering::MMD(){
    //    assert(pStructure!=NULL);

    logfileptr->OFS()<<"MMD used"<<endl;
    if(iam==0){cout<<"MMD used"<<endl;}

    if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
      throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call MMD\n" );
    }

    invp.resize(pStructure->size);
    if(iam==0){

      SYMPACK::vector<MMDInt> perm(pStructure->size);
      SYMPACK::vector<MMDInt> perminv(pStructure->size);
      MMDInt iwsiz = 4*pStructure->size;
      SYMPACK::vector<MMDInt> iwork (iwsiz);

      MMDInt nadj = pStructure->expRowind.size();
      MMDInt nofsub =0;
      MMDInt iflag =0;

      SYMPACK::vector<MMDInt> tmpXadj(pStructure->expColptr.size());
      for(MMDInt col=0; col<tmpXadj.size();++col){
        tmpXadj[col] = pStructure->expColptr[col];
      }


      SYMPACK::vector<MMDInt> tmpAdj(pStructure->expRowind.size()-pStructure->size);
      MMDInt pos = 0;
      for(MMDInt col=0; col<tmpXadj.size()-1;++col){
        for(MMDInt j=tmpXadj[col]; j<tmpXadj[col+1];++j){
          if( pStructure->expRowind[j-1]-1 != col){
            tmpAdj[pos++] = pStructure->expRowind[j-1];
          }
        }
      }

      MMDInt rm = 0;
      for(MMDInt col=0; col<tmpXadj.size();++col){
        tmpXadj[col]-=rm;
        rm++;
      }

      MMDInt N = pStructure->size;
      FORTRAN(ordmmd)( &N , &nadj , &tmpXadj[0] , &tmpAdj[0], 
          &perminv[0] , &perm[0] , &iwsiz , &iwork[0] , &nofsub, &iflag ) ;


      assert(iflag == 0);

      for(Int i = 1; i <=N; ++i){
        invp[i-1] = perminv[i-1];
      }
      //logfileptr->OFS()<<perm<<endl;
    }
    // broadcast invp
    Int N = invp.size();
    MPI_Bcast(&invp[0],N*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
    perm.resize(pStructure->size);
    for(Int i = 1; i <=N; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }

    //logfileptr->OFS()<<perm<<endl;


  }

  void Ordering::NDBOX(){
    //    assert(pStructure!=NULL);

    if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
      throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call NDBOX\n" );
    }

    invp.resize(pStructure->size);
    if(iam==0){
      Int k = std::ceil(std::pow(pStructure->size,1.0/3.0));
      logfileptr->OFS()<<"BOX K = "<<k<<endl;
      Int iflag =1;

      {
        SYMPACK::vector<Int> stack(k*k>25?invp.size():2*invp.size());
        Int tmp = stack.size();
        FORTRAN(boxnd)( &k , &k, &k, &invp[0], &stack[0],&tmp, &iflag);
      }

    }

    // broadcast invp
    Int N = invp.size();
    MPI_Bcast(&invp[0],N*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
    perm.resize(pStructure->size);
    for(Int i = 1; i <=N; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }


  }

  void Ordering::NDGRID(){
    //    assert(pStructure!=NULL);

    if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
      throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call NDGRID\n" );
    }

    invp.resize(pStructure->size);
    if(iam==0){
      Int k = std::ceil(std::pow(pStructure->size,1.0/2.0));
      Int iflag =1;

      logfileptr->OFS()<<"GRID K = "<<k<<endl;

      {
        SYMPACK::vector<Int> stack(k*k>25?invp.size():2*invp.size());
        Int tmp = stack.size();
        FORTRAN(gridnd)( &k , &k, &invp[0], &stack[0],&tmp, &iflag);
        assert(iflag==0);
      }
    }

    // broadcast invp
    Int N = invp.size();
    MPI_Bcast(&invp[0],N*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
    perm.resize(pStructure->size);
    for(Int i = 1; i <=N; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }



    //logfileptr->OFS()<<"perm "<<perm<<endl;
    //logfileptr->OFS()<<"invperm "<<invp<<endl;
  }


  void Ordering::AMD(){
    //  assert(pStructure!=NULL);

    logfileptr->OFS()<<"AMD used"<<endl;
    if(iam==0){cout<<"AMD used"<<endl;}

    if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
      throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call AMD\n" );
    }

    invp.resize(pStructure->size);
    if(iam==0){
      SYMPACK::vector<AMDInt> perm(pStructure->size);
      SYMPACK::vector<AMDInt> perminv(pStructure->size);
      AMDInt iwsiz = 4*pStructure->size;
      SYMPACK::vector<AMDInt> iwork (iwsiz);

      AMDInt nadj = pStructure->expRowind.size();
      AMDInt nofsub =0;
      AMDInt iflag =0;

      AMDInt N = pStructure->size; 

      SYMPACK::vector<AMDInt> tmpXadj(pStructure->expColptr.size());
      for(AMDInt col=0; col<tmpXadj.size();++col){
        tmpXadj[col] = pStructure->expColptr[col];
      }
      //allocate extra N elements for AMD (elbow)
      SYMPACK::vector<AMDInt> tmpAdj(pStructure->expRowind.size()+N);
      for(AMDInt i=0; i<pStructure->expRowind.size();++i){
        tmpAdj[i] = pStructure->expRowind[i];
      }

      AMDInt IOVFLO = std::numeric_limits<AMDInt>::max();//2147483647;
      AMDInt NCMPA;

      SYMPACK::vector<AMDInt> VTXDEG(N);
      SYMPACK::vector<AMDInt> QSIZE(N);
      SYMPACK::vector<AMDInt> ECFORW(N);
      SYMPACK::vector<AMDInt> MARKER(N);
      SYMPACK::vector<AMDInt> NVTXS(N+1);
      for(AMDInt i=0;i<N;++i){
        NVTXS[i] = tmpXadj[i+1]-tmpXadj[i];
      }


      AMDInt IWLEN = tmpXadj[N-1] + NVTXS[N-1] + N - 1 ;
      AMDInt PFREE = tmpXadj[N-1] + NVTXS[N-1];


      FORTRAN(amdbar)( &N , &tmpXadj[0] , &tmpAdj[0] , &NVTXS[0] ,&IWLEN ,&PFREE ,  &QSIZE[0],&ECFORW[0] , &perm[0], &iwork[0], &perminv[0],&VTXDEG[0], &NCMPA , &MARKER[0] , &IOVFLO );

      for(AMDInt i = 0; i < perminv.size(); ++i){
        invp[i] = perminv[i];
      }
    }

    // broadcast invp
    Int N = invp.size();
    MPI_Bcast(&invp[0],N*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
    perm.resize(pStructure->size);
    for(Int i = 1; i <=N; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }


    logfileptr->OFS()<<"perm "<<perm<<endl;
    logfileptr->OFS()<<"invperm "<<invp<<endl;

    //    Permute(perm);

  }


#ifdef USE_METIS
  void Ordering::METIS(){
    //    assert(pStructure!=NULL);

    logfileptr->OFS()<<"METIS used"<<endl;
    if(iam==0){cout<<"METIS used"<<endl;}


    if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
      throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call METIS\n" );
    }


    invp.resize(pStructure->size);
    if(iam==0){
      perm.resize(pStructure->size);

#if 0
      int options[METIS_NOPTIONS];
      METIS_SetDefaultOptions(&options[0]);
      options[METIS_OPTION_NUMBERING] = 1;
      options[METIS_OPTION_DBGLVL] = 8;
#endif

      //options[METIS_OPTION_CCORDER] = 1;
      //options[METIS_OPTION_COMPRESS] = 1;
      // int OPTION[8];
      // OPTION[0] = 0;
      // OPTION[1] = 3;
      // OPTION[2] = 1;
      // OPTION[3] = 2;
      // OPTION[4] = 0;
      // OPTION[5] = 1;
      // OPTION[6] = 0;
      // OPTION[7] = 1;

      int N = pStructure->size; 

      SYMPACK::vector<int> tmpXadj(pStructure->expColptr.size());
      for(int i = 0; i<tmpXadj.size();++i){tmpXadj[i] = pStructure->expColptr[i];}
      SYMPACK::vector<int> tmpAdj;
      tmpAdj.reserve(pStructure->expRowind.size());

      for(int col=0; col<tmpXadj.size()-1;++col){
        for(int j=tmpXadj[col]; j<tmpXadj[col+1];++j){
          if( pStructure->expRowind[j-1]-1 != col){
            tmpAdj.push_back(pStructure->expRowind[j-1]);
          }
        }
      }

      int rm = 0;
      for(int col=0; col<tmpXadj.size();++col){
        tmpXadj[col]-=rm;
        rm++;
      }

      //    cout<<tmpXadj<<endl;
      //    cout<<tmpAdj<<endl;

      //switch everything to 0 based
      for(int col=0; col<tmpXadj.size();++col){ tmpXadj[col]--;}
      for(int col=0; col<tmpAdj.size();++col){ tmpAdj[col]--;}

      //    cout<<tmpXadj<<endl;
      //    cout<<tmpAdj<<endl;

#ifdef SCOTCH_METIS
      Int numflag = 0;
      METIS_NodeND( &N, &tmpXadj[0] , &tmpAdj[0]  , &numflag , NULL, &perm[0] , &invp[0] );
#else
      METIS_NodeND( &N, &tmpXadj[0] , &tmpAdj[0]  , NULL , NULL, &perm[0] , &invp[0] );
#endif

      //#ifdef USE_METIS_INC
      //    metis_nodend( &N, &tmpXadj[0] , &tmpAdj[0]  , NULL , NULL, &perm[0] , &invp[0] );
      //#else
      //    int numflag = 0;
      //    metis_nodend( &N, &tmpXadj[0] , &tmpAdj[0]  , &numflag , NULL, &perm[0] , &invp[0] );
      //#endif

      //switch everything to 1 based
      for(int col=0; col<invp.size();++col){ invp[col]++;}
      for(int col=0; col<perm.size();++col){ perm[col]++;}
    }


    // broadcast invp
    Int N = invp.size();
    MPI_Bcast(&invp[0],N*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
    perm.resize(pStructure->size);
    for(Int i = 1; i <=N; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }




  }



  void Ordering::METIS(const SparseMatrixGraph & g){
    logfileptr->OFS()<<"METIS used"<<endl;
    if(iam==0){cout<<"METIS used"<<endl;}

    if(iam == 0 && (!g.IsExpanded() || g.keepDiag==1) ){
      throw std::logic_error( "SparseMatrixGraph must be expanded and not including the diagonal in order to call METIS\n" );
    }


    idx_t baseval = g.baseval;
    idx_t N = (idx_t)g.size; 

    bool isSameInt = typeid(idx_t) == typeid(Int);
    bool isSameIdx = typeid(idx_t) == typeid(Idx);
    bool isSamePtr = typeid(idx_t) == typeid(Ptr);


    invp.resize(N);

    if(iam==0){
      idx_t * iperm;      
      if(!isSameInt){
        iperm = new idx_t[N];
      }
      else{
        iperm = (idx_t*)&invp[0];
      }

      idx_t * mperm;      
      if(!isSameInt){
        mperm = new idx_t[N];
      }
      else{
        mperm = (idx_t*)&perm[0];
      }


      idx_t * prowind = NULL;
      if(!isSameIdx){
        prowind = new idx_t[g.EdgeCount()];
        for(Ptr i = 0; i<g.rowind.size();i++){ prowind[i] = (idx_t)g.rowind[i];}
      }
      else{
        prowind = (idx_t*)&g.rowind[0];
      }

      idx_t * pcolptr = NULL;
      if(!isSamePtr){
        pcolptr = new idx_t[g.VertexCount()+1];
        for(Ptr i = 0; i<g.colptr.size();i++){ pcolptr[i] = (idx_t)g.colptr[i];}
      }
      else{
        pcolptr = (idx_t*)&g.colptr[0];
      }


      idx_t options[METIS_NOPTIONS];
      METIS_SetDefaultOptions(options);
      options[METIS_OPTION_NUMBERING] = g.baseval;

      METIS_NodeND(&N, pcolptr, prowind, NULL,
          options, mperm, iperm);


      if(!isSameInt){ 
        delete [] mperm;
      }

      if(!isSameInt ||g.baseval!=1){ 
        //switch everything to 1 based
        for(int col=0; col<N;++col){ invp[col] = iperm[col]+(1-g.baseval);}
      }

      if(!isSameInt){ 
        delete [] iperm;
      }

      if(!isSamePtr){
        delete [] pcolptr;
      }

      if(!isSameIdx){
        delete [] prowind;
      }
    }
    // broadcast invp
    MPI_Bcast(&invp[0],N*sizeof(Int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
    perm.resize(N);
    for(Int i = 1; i <=N; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }


  }







#endif


#ifdef USE_PARMETIS
  void Ordering::PARMETIS(){
    //  assert(pStructure!=NULL);

    logfileptr->OFS()<<"PARMETIS used"<<endl;
    if(iam==0){cout<<"PARMETIS used"<<endl;}


    if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
      throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call AMD\n" );
    }


    int baseval = 1;
    int N = pStructure->size; 
    invp.resize(N);
    perm.resize(N);
    int ndomains = (int)pow(2.0,std::floor(std::log2(np)));

    MPI_Comm ndcomm;
    Int mpirank;
    MPI_Comm_split(CommEnv_->MPI_GetComm(),iam<ndomains,iam,&ndcomm);
    if(iam<ndomains){
      SYMPACK::vector<int> sizes(2*ndomains);
      SYMPACK::vector<int> vtxdist(ndomains+1);

      int localN = N/ndomains;
      int fc = (iam)*localN+1;
      //build vtxdist SYMPACK::vector
      for(int i = 0; i<ndomains;++i){ vtxdist[i] = i*localN+baseval; } vtxdist[ndomains] = N+baseval;
      if(iam==ndomains-1){
        localN = N - (ndomains-1)*localN;
      }
      //logfileptr->OFS()<<"vtxdist: "<<vtxdist<<endl;
      //logfileptr->OFS()<<"ndomains = "<<ndomains<<endl;

      //build local colptr and count nnz in local rowind
      SYMPACK::vector<int> tmpXadj(localN+1);
      Ptr localNNZ = 0;
      int offset = pStructure->expColptr[fc-1]-1;
      tmpXadj[0] = 1;
      for(int col=1; col<=localN;++col){
        int column = col + fc -1;
        Ptr colbeg = pStructure->expColptr[column-1];
        Ptr colend = pStructure->expColptr[column]-1;
        //remove self from adjacency list
        tmpXadj[col] = tmpXadj[col-1] + colend - colbeg + 1 - 1;
        //logfileptr->OFS()<<col<<" <-> "<<column<<" | "<<colbeg<<"--"<<colend<<" <-> "<<tmpXadj[col-1]<<"--"<<tmpXadj[col]-1<<endl;
        localNNZ+=colend-colbeg+1 -1;
      }

      //build local rowind, discarding self
      SYMPACK::vector<int> tmpAdj;
      tmpAdj.resize(localNNZ);
      Ptr pos =0;
      for(int col=1; col<=localN;++col){
        int column = col + fc -1;
        Ptr colbeg = pStructure->expColptr[column-1];
        Ptr colend = pStructure->expColptr[column]-1;
        for(Ptr j=colbeg; j<=colend;++j){
          int row = pStructure->expRowind[j-1];
          if(row!=column){
            tmpAdj[pos++] = row;
          }
        }
      }


      //switch everything to 0 based
      //logfileptr->OFS()<<tmpXadj<<endl;
      for(int col=0; col<localN+1;++col){ tmpXadj[col]+=(baseval-1);}
      //logfileptr->OFS()<<tmpXadj<<endl;
      //logfileptr->OFS()<<tmpAdj<<endl;
      for(int col=0; col<localNNZ;++col){ tmpAdj[col]+=(baseval-1);}
      //logfileptr->OFS()<<tmpAdj<<endl;

      int options[3];
      options[0] = 0;
      int numflag = baseval;


      int npnd;
      MPI_Comm_size (ndcomm, &npnd);
      //logfileptr->OFS()<<"PROC ND: "<<npnd<<endl;
      ParMETIS_V3_NodeND( &vtxdist[0], &tmpXadj[0] , &tmpAdj[0], &numflag, &options[0], &perm[0], &sizes[0], &ndcomm );
      //    ParMETIS_V3_NodeND( &vtxdist[0], &tmpXadj[0] , &tmpAdj[0], &numflag, &options[0], &invp[0], &sizes[0], &comm );

      //logfileptr->OFS()<<"Sizes: "<<sizes<<endl;
      //logfileptr->OFS()<<"Order: "<<endl;
      //for(int i =0;i<localN;++i){
      //  //logfileptr->OFS()<<invp[i]<<" ";
      //  logfileptr->OFS()<<perm[i]<<" ";
      //}
      //logfileptr->OFS()<<endl;


      //compute displs
      SYMPACK::vector<int> mpidispls(ndomains,0);
      SYMPACK::vector<int> mpisizes(ndomains,0);
      for(int p = 1;p<=ndomains;++p){
        mpisizes[p-1] = (vtxdist[p] - vtxdist[p-1])*sizeof(int);
        mpidispls[p-1] = (vtxdist[p-1]-baseval)*sizeof(int);
      }

      //gather on the root
      MPI_Gatherv(&perm[0],mpisizes[iam],MPI_BYTE,&invp[0],&mpisizes[0],&mpidispls[0],MPI_BYTE,0,ndcomm);
      //    MPI_Gatherv(&invp[0],mpisizes[iam],MPI_BYTE,&perm[0],&mpisizes[0],&mpidispls[0],MPI_BYTE,0,comm);

      if(iam==0){
        //switch everything to 1 based
        for(int col=0; col<N;++col){ invp[col]+=(1-baseval);}

        //logfileptr->OFS()<<"Full Order: "<<endl;
        //for(int i =0;i<N;++i){
        //  logfileptr->OFS()<<invp[i]<<" ";
        //  //logfileptr->OFS()<<perm[i]<<" ";
        //}
        //logfileptr->OFS()<<endl;
      }


    }

    MPI_Comm_free(&ndcomm);

    // broadcast invp
    MPI_Bcast(&invp[0],N*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
    for(Int i = 1; i <=N; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }

    //  // broadcast perm
    //  MPI_Bcast(&perm[0],N*sizeof(int),MPI_BYTE,0,comm);
    //
    //  for(Int i = 1; i <= perm.size(); ++i){
    //    Int node = perm[i-1];
    //    invp[node-1] = i;
    //  }
  }


  void Ordering::PARMETIS(const DistSparseMatrixGraph & g){

    logfileptr->OFS()<<"PARMETIS used"<<endl;
    if(iam==0){cout<<"PARMETIS used"<<endl;}

    if(!g.IsExpanded() || g.keepDiag==1){
      throw std::logic_error( "DistSparseMatrixGraph must be expanded and not including the diagonal in order to call PARMETIS\n" );
    }


    idx_t baseval = g.baseval;
    idx_t N = g.size; 
    invp.resize(N);

    int mpirank;

    idx_t ndomains = (idx_t)pow(2.0,std::floor(std::log2(np)));
    //make sure every one have an element
    while(((double)N/(double)ndomains)<1.0){ndomains /= 2;}

    MPI_Comm ndcomm;
    MPI_Comm_split(g.comm,iam<ndomains,iam,&ndcomm);
    //      while(((double)N/(double)ndomains)<1.0){ndomains--;}

    MPI_Comm_rank(ndcomm,&mpirank);
    SYMPACK::vector<idx_t> vtxdist;
    idx_t localN;
    if(iam<ndomains){
      assert(mpirank==iam);
      //vtxdist.resize(ndomains+1);

      //localN = g.LocalVertexCount();
      ////build vtxdist SYMPACK::vector
      //for(idx_t i = 0; i<ndomains;++i){
      // vtxdist[i] = i*(N/ndomains)+baseval; 
      //} 
      //vtxdist[ndomains] = N+baseval;

      SYMPACK::vector<idx_t> sizes(2*ndomains);
      vtxdist.resize(g.vertexDist.size());
      for(int i = 0 ; i < g.vertexDist.size(); i++){
        vtxdist[i] = (idx_t)g.vertexDist[i];
      }

      logfileptr->OFS()<<vtxdist<<endl;

      //        if(iam==ndomains-1){
      //          localN = N - (ndomains-1)*localN;
      //        }

      idx_t * iperm = NULL;
      if(typeid(idx_t) != typeid(Int)){
        iperm = new idx_t[N];
      }
      else{
        iperm = (idx_t*)&invp[0];
      }

      idx_t * pperm = NULL;
      if(typeid(idx_t) != typeid(Int)){
        pperm = new idx_t[N];
      }
      else{
        perm.resize(N);
        pperm = (idx_t*)&perm[0];
      }

      idx_t * prowind = NULL;
      if(typeid(idx_t) != typeid(Idx)){
        prowind = new idx_t[g.LocalEdgeCount()];
        for(Ptr i = 0; i<g.rowind.size();i++){ prowind[i] = (idx_t)g.rowind[i];}
      }
      else{
        prowind = (idx_t*)&g.rowind[0];
      }

      idx_t * pcolptr = NULL;
      if(typeid(idx_t) != typeid(Ptr)){
        pcolptr = new idx_t[g.LocalVertexCount()+1];
        for(Ptr i = 0; i<g.colptr.size();i++){ pcolptr[i] = (idx_t)g.colptr[i];}
      }
      else{
        pcolptr = (idx_t*)&g.colptr[0];
      }








      idx_t options[3];
      options[0] = 0;
      idx_t numflag = g.baseval;

      int npnd;
      MPI_Comm_size (ndcomm, &npnd);
      ParMETIS_V3_NodeND( &vtxdist[0], pcolptr , prowind, &numflag, &options[0], pperm, &sizes[0], &ndcomm );

      //compute displs
      SYMPACK::vector<int> mpidispls(ndomains,0);
      SYMPACK::vector<int> mpisizes(ndomains,0);
      for(int p = 1;p<=ndomains;++p){
        mpisizes[p-1] = (vtxdist[p] - vtxdist[p-1])*sizeof(idx_t);
        mpidispls[p-1] = (vtxdist[p-1]-baseval)*sizeof(idx_t);
      }

      //gather on the root
      MPI_Gatherv(pperm,mpisizes[iam],MPI_BYTE,iperm,&mpisizes[0],&mpidispls[0],MPI_BYTE,0,ndcomm);
      //    MPI_Gatherv(&invp[0],mpisizes[iam],MPI_BYTE,&perm[0],&mpisizes[0],&mpidispls[0],MPI_BYTE,0,comm);

      if(iam==0 && (g.baseval!=1 || typeid(idx_t) != typeid(Int))){
        //switch everything to 1 based
        for(int col=0; col<N;++col){ invp[col] = iperm[col] + (1-baseval);}
      }

      if(typeid(idx_t) != typeid(Int)){
        delete [] pperm;
      }

      if(typeid(idx_t) != typeid(Int)){
        delete [] iperm;
      }

      if(typeid(idx_t) != typeid(Ptr)){
        delete [] pcolptr;
      }

      if(typeid(idx_t) != typeid(Idx)){
        delete [] prowind;
      }






    }

    MPI_Comm_free(&ndcomm);
    vtxdist.clear();

    // broadcast invp
    MPI_Bcast(&invp[0],N*sizeof(Int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
    //recompute perm
    perm.resize(N);
    for(Int i = 1; i <=N; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }

    logfileptr->OFS()<<perm<<endl;
    logfileptr->OFS()<<invp<<endl;


  }




#endif

#ifdef USE_SCOTCH
  void Ordering::SCOTCH(){
    //assert(pStructure!=NULL);

    logfileptr->OFS()<<"SCOTCH used"<<endl;
    if(iam==0){cout<<"SCOTCH used"<<endl;}


    if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
      throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call SCOTCH\n" );
    }


    invp.resize(pStructure->size);
    if(iam==0){
      perm.resize(pStructure->size);

      SYMPACK::vector<SCOTCH_Num> tmpInvp(pStructure->size);
      SYMPACK::vector<SCOTCH_Num> tmpPerm(pStructure->size);

      SCOTCH_Num N = pStructure->size; 

      SYMPACK::vector<SCOTCH_Num> tmpXadj(pStructure->expColptr.size());
      for(SCOTCH_Num i = 0; i<tmpXadj.size();++i){tmpXadj[i] = pStructure->expColptr[i];}

      SYMPACK::vector<SCOTCH_Num> tmpAdj;
      tmpAdj.reserve(pStructure->expRowind.size());

      for(SCOTCH_Num col=0; col<tmpXadj.size()-1;++col){
        for(SCOTCH_Num j=tmpXadj[col]; j<tmpXadj[col+1];++j){
          if( pStructure->expRowind[j-1]-1 != col){
            tmpAdj.push_back(pStructure->expRowind[j-1]);
          }
        }
      }

      SCOTCH_Num rm = 0;
      for(SCOTCH_Num col=0; col<tmpXadj.size();++col){
        tmpXadj[col]-=rm;
        rm++;
      }

      //switch everything to 0 based
      for(SCOTCH_Num col=0; col<tmpXadj.size();++col){ tmpXadj[col]--;}
      for(SCOTCH_Num col=0; col<tmpAdj.size();++col){ tmpAdj[col]--;}

      int numflag = 0;
      {
        SCOTCH_Graph        grafdat;                    /* Scotch graph object to interface with libScotch    */
        SCOTCH_Ordering     ordedat;                    /* Scotch ordering object to interface with libScotch */
        SCOTCH_Strat        stradat;

        SCOTCH_graphInit (&grafdat);

        if (SCOTCH_graphBuild (&grafdat,
              numflag, N, &tmpXadj[0], &tmpXadj[0] + 1, NULL, NULL,
              tmpXadj[N] - numflag, &tmpAdj[0], NULL) == 0) {

          //STRSTRING='n{sep=m{asc=b{width=3,strat=q{strat=f}},'//
          // 'low=q{strat=h},vert=1000,dvert=100,dlevl=0,'//
          // 'proc=1,seq=q{strat=m{type=h,vert=100,'//
          // 'low=h{pass=10},asc=b{width=3,bnd=f{bal=0.2},'//


          //switch (stratnum) {
          //case 1: // nested dissection with STRATPAR levels
          //strategy_string << "n{ole=s,ose=s,sep=(/levl>" << int(stratpar-1) << "?z:"
          //"(m{asc=b{bnd=f{move=200,pass=1000,bal=0.9},org=(|h{pass=10})f{move=200,pass=1000,bal=0.9},width=3},"
          //"low=h{pass=10},type=h,vert=100,rat=0.7});)}"; break;
          //case 2: // nested dissection with separators <= STRATPAR
          //strategy_string << "n{ole=s,ose=s,sep=(/vert<" << int(stratpar) << "?z:"
          //"(m{asc=b{bnd=f{move=200,pass=1000,bal=0.9},org=(|h{pass=10})"
          //"f{move=200,pass=1000,bal=0.9},width=3},"
          //"low=h{pass=10},type=h,vert=100,rat=0.7});)}"; break;
          //default: // Pure nested dissection
          //strategy_string << "n{ole=s,ose=s,sep=m{asc=b{bnd=f{move=200,pass=1000,bal=0.9},"
          //"org=(|h{pass=10})f{move=200,pass=1000,bal=0.9},width=3},low=h{pass=10},type=h,vert=100,rat=0.7}}";
          //}

          SCOTCH_stratInit (&stradat);

          stringstream strategy_string;
          //strategy_string << "n{ole=s,ose=s,sep=m{asc=b{bnd=f{move=200,pass=1000,bal=0.9},"
          //"org=(|h{pass=10})f{move=200,pass=1000,bal=0.9},width=3},low=h{pass=10},vert=100,rat=0.7}}";



          //int vert = 120;
          //int cmin = 1;//20;
          //int cmax = 100000;
          //int frat = 0.08;
          //
          //        strategy_string << "c{rat=0.7,"         
          //        "cpr=n{sep=/(vert>"<<vert<<")?m{vert=100,"         
          //                          "low=h{pass=10},"   
          //                          "asc=f{bal=0.2}}|"  
          //                        "m{vert=100,"         
          //                          "low=h{pass=10},"   
          //                          "asc=f{bal=0.2}};," 
          //        "ole=f{cmin="<<cmin<<",cmax="<<cmax<<",frat="<<cmax<<"},"   
          //        "ose=g},"                             
          //        "unc=n{sep=/(vert>"<<vert<<")?(m{vert=100,"        
          //                           "low=h{pass=10},"  
          //                           "asc=f{bal=0.2}})|"
          //                         "m{vert=100,"        
          //                           "low=h{pass=10},"  
          //                           "asc=f{bal=0.2}};,"
          //        "ole=f{cmin="<<cmin<<",cmax="<<cmax<<",frat="<<cmax<<"},"   
          //        "ose=g}}";
          //
          //
          //        int err = SCOTCH_stratGraphOrder(&stradat,strategy_string.str().c_str());
          ////"d{cmin=20,frat=0.08}");
          //        assert(err==0);

#ifdef SCOTCH_DEBUG_ALL
          if (SCOTCH_graphCheck (&grafdat) == 0)        /* TRICK: next instruction called only if graph is consistent */
#endif /* SCOTCH_DEBUG_ALL */
          {
            if (SCOTCH_graphOrderInit (&grafdat, &ordedat, &tmpInvp[0], &tmpPerm[0], /* MeTiS and Scotch have opposite definitions for (inverse) permutations */
                  /*nb de supernoeud*/NULL, /*ptr vers rank tab : xsuper*/ NULL, /*tree tab: parent structure*/ NULL) == 0) {
              SCOTCH_graphOrderCompute (&grafdat, &ordedat, &stradat);
              SCOTCH_graphOrderExit    (&grafdat, &ordedat);
            }
          }
          SCOTCH_stratExit (&stradat);
        }
        SCOTCH_graphExit (&grafdat);
        }


        //switch everything to 1 based
        for(int col=0; col<invp.size();++col){ invp[col] = tmpInvp[col]+1;}
        for(int col=0; col<perm.size();++col){ perm[col] = tmpPerm[col]+1;}
        }

        // broadcast invp
        Int N = invp.size();
        MPI_Bcast(&invp[0],N*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
        perm.resize(pStructure->size);
        for(Int i = 1; i <=N; ++i){
          Int node = invp[i-1];
          perm[node-1] = i;
        }


        }



        void Ordering::SCOTCH(const SparseMatrixGraph & g){
          logfileptr->OFS()<<"SCOTCH used"<<endl;
          if(iam==0){cout<<"SCOTCH used"<<endl;}

          if(iam == 0 && (!g.IsExpanded() || g.keepDiag==1) ){
            throw std::logic_error( "SparseMatrixGraph must be expanded and not including the diagonal in order to call SCOTCH\n" );
          }


          SCOTCH_Num baseval = g.baseval;
          SCOTCH_Num N = g.size; 

          bool isSameInt = typeid(SCOTCH_Num) == typeid(Int);
          bool isSameIdx = typeid(SCOTCH_Num) == typeid(Idx);
          bool isSamePtr = typeid(SCOTCH_Num) == typeid(Ptr);


          invp.resize(N);

          if(iam==0){
            SCOTCH_Num * permtab;      
            if(!isSameInt){
              permtab = new SCOTCH_Num[N];
            }
            else{
              permtab = (SCOTCH_Num*)&invp[0];
            }

            SCOTCH_Num * prowind = NULL;
            if(!isSameIdx){
              prowind = new SCOTCH_Num[g.EdgeCount()];
              for(Ptr i = 0; i<g.rowind.size();i++){ prowind[i] = (SCOTCH_Num)g.rowind[i];}
            }
            else{
              prowind = (SCOTCH_Num*)&g.rowind[0];
            }

            SCOTCH_Num * pcolptr = NULL;
            if(!isSamePtr){
              pcolptr = new SCOTCH_Num[g.VertexCount()+1];
              for(Ptr i = 0; i<g.colptr.size();i++){ pcolptr[i] = (SCOTCH_Num)g.colptr[i];}
            }
            else{
              pcolptr = (SCOTCH_Num*)&g.colptr[0];
            }


            SCOTCH_Num nnz = g.EdgeCount();
            SCOTCH_Num baseval = g.baseval;
            {
              SCOTCH_Graph        grafdat;                    /* Scotch graph object to interface with libScotch    */
              SCOTCH_Ordering     ordedat;                    /* Scotch ordering object to interface with libScotch */
              SCOTCH_Strat        stradat;

              SCOTCH_graphInit (&grafdat);

              if (SCOTCH_graphBuild (
                    &grafdat,
                    baseval, 
                    N, 
                    pcolptr, 
                    NULL, 
                    NULL,
                    NULL,
                    nnz, 
                    prowind, 
                    NULL) == 0) {

                SCOTCH_stratInit (&stradat);

                stringstream strategy_string;
#ifdef SCOTCH_DEBUG_ALL
                if (SCOTCH_graphCheck (&grafdat) == 0)        /* TRICK: next instruction called only if graph is consistent */
#endif /* SCOTCH_DEBUG_ALL */
                {
                  if (SCOTCH_graphOrderInit (
                        &grafdat,
                        &ordedat,
                        permtab,
                        NULL,
                        /*nb de supernoeud*/NULL,
                        /*ptr vers rank tab : xsuper*/ NULL,
                        /*tree tab: parent structure*/ NULL) == 0) {
                    SCOTCH_graphOrderCompute (&grafdat, &ordedat, &stradat);
                    SCOTCH_graphOrderExit    (&grafdat, &ordedat);
                  }
                }
                SCOTCH_stratExit (&stradat);
              }
              SCOTCH_graphExit (&grafdat);

              if(!isSameInt ||g.baseval!=1){ 
                //switch everything to 1 based
                for(int col=0; col<N;++col){ invp[col] = permtab[col]+(1-baseval);}
              }

              if(!isSameInt){ 
                delete [] permtab;
              }

              if(!isSamePtr){
                delete [] pcolptr;
              }

              if(!isSameIdx){
                delete [] prowind;
              }
            }
          }
          // broadcast invp
          MPI_Bcast(&invp[0],N*sizeof(Int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
          perm.resize(N);
          for(Int i = 1; i <=N; ++i){
            Int node = invp[i-1];
            perm[node-1] = i;
          }


        }







#endif



#ifdef USE_PTSCOTCH
        void Ordering::PTSCOTCH(const DistSparseMatrixGraph & g){

          logfileptr->OFS()<<"PTSCOTCH used"<<endl;
          if(iam==0){cout<<"PTSCOTCH used"<<endl;}

          if(!g.IsExpanded() || g.keepDiag==1){
            throw std::logic_error( "DistSparseMatrixGraph must be expanded and not including the diagonal in order to call PTSCOTCH\n" );
          }


          SCOTCH_Num baseval = g.baseval;
          SCOTCH_Num N = g.size; 
          invp.resize(N);

          int mpirank;

          //      SCOTCH_Num ndomains = (SCOTCH_Num)pow(2.0,std::floor(std::log2(np)));
          //      //make sure every one have an element
          //      while(((double)N/(double)ndomains)<1.0){ndomains /= 2;}
          //      MPI_Comm ndcomm;
          //      int color = iam < ndomains;
          //      MPI_Comm_split(g.comm,color,iam,&ndcomm);


          SCOTCH_Num ndomains = np;
          MPI_Comm ndcomm;
          MPI_Comm_dup(g.comm,&ndcomm);
          //      while(((double)N/(double)ndomains)<1.0){ndomains--;}

          MPI_Comm_rank(ndcomm,&mpirank);
          SYMPACK::vector<SCOTCH_Num> vtxdist;
          SCOTCH_Num localN;
          if(iam<ndomains){
            assert(mpirank==iam);
            //vtxdist.resize(ndomains+1);

            //localN = g.LocalVertexCount();
            ////build vtxdist SYMPACK::vector
            //for(SCOTCH_Num i = 0; i<ndomains;++i){
            // vtxdist[i] = i*(N/ndomains)+baseval; 
            //} 
            //vtxdist[ndomains] = N+baseval;

            vtxdist.resize(g.vertexDist.size());
            for(int i = 0 ; i < g.vertexDist.size(); i++){
              vtxdist[i] = (SCOTCH_Num)g.vertexDist[i];
            }

            logfileptr->OFS()<<vtxdist<<endl;

            //        if(iam==ndomains-1){
            //          localN = N - (ndomains-1)*localN;
            //        }


            SCOTCH_Num * prowind = NULL;
            if(typeid(SCOTCH_Num) != typeid(Idx)){
              prowind = new SCOTCH_Num[g.LocalEdgeCount()];
              for(Ptr i = 0; i<g.rowind.size();i++){ prowind[i] = (SCOTCH_Num)g.rowind[i];}
            }
            else{
              prowind = (SCOTCH_Num*)&g.rowind[0];
            }

            SCOTCH_Num * pcolptr = NULL;
            if(typeid(SCOTCH_Num) != typeid(Ptr)){
              pcolptr = new SCOTCH_Num[g.LocalVertexCount()+1];
              for(Ptr i = 0; i<g.colptr.size();i++){ pcolptr[i] = (SCOTCH_Num)g.colptr[i];}
            }
            else{
              pcolptr = (SCOTCH_Num*)&g.colptr[0];
            }








            SCOTCH_Num options[3];
            options[0] = 0;
            SCOTCH_Num numflag = baseval;

            int npnd;
            MPI_Comm_size (ndcomm, &npnd);

            SCOTCH_Dgraph       grafdat;                    /* Scotch distributed graph object to SCOTCH_Numerface with libScotch    */
            SCOTCH_Dordering    ordedat;                    /* Scotch distributed ordering object to SCOTCH_Numerface with libScotch */
            SCOTCH_Strat        stradat;
            SCOTCH_Num          vertlocnbr;
            SCOTCH_Num          edgelocnbr;

            vertlocnbr = g.LocalVertexCount();
            edgelocnbr = g.LocalEdgeCount();

            if(SCOTCH_dgraphInit (&grafdat, ndcomm) != 0){
              throw std::logic_error( "Error in SCOTCH_dgraphInit\n" );
            }

            //logfileptr->OFS()<<"vtxdist: "<<vtxdist<<endl;
            //logfileptr->OFS()<<"vertlocnbr: "<<vertlocnbr<<endl;
            //logfileptr->OFS()<<"colptr: "<<g.colptr<<endl;
            //logfileptr->OFS()<<"edgelocnbr: "<<edgelocnbr<<endl;
            //logfileptr->OFS()<<"rowind: "<<g.rowind<<endl;
            //assert(vertlocnbr == g.colptr.size()-1);

            //      logfileptr->OFS()<<"Entering SCOTCH_dgraphBuild"<<endl;
            if (SCOTCH_dgraphBuild (&grafdat, 
                  baseval,
                  vertlocnbr, 
                  vertlocnbr, 
                  pcolptr, 
                  NULL,
                  NULL, 
                  NULL,
                  edgelocnbr,
                  edgelocnbr,
                  prowind, 
                  NULL,
                  NULL
                  ) != 0) {
              throw std::logic_error( "Error in SCOTCH_dgraphBuild\n" );
            }

            if(SCOTCH_dgraphCheck (&grafdat) != 0){
              throw std::logic_error( "Error in SCOTCH_dgraphCheck\n" );
            }

            if(SCOTCH_stratInit (&stradat)!= 0){
              throw std::logic_error( "Error in SCOTCH_stratInit\n" );
            }

            if (SCOTCH_dgraphOrderInit (&grafdat, &ordedat) != 0) {
              throw std::logic_error( "Error in SCOTCH_dgraphOrderInit\n" );
            }


            if(SCOTCH_dgraphOrderCompute (&grafdat, &ordedat, &stradat)!=0){
              throw std::logic_error( "Error in SCOTCH_dgraphOrderCompute\n" );
            }










#if 1
            SYMPACK::vector<SCOTCH_Num> sc_permtab(N);
            SYMPACK::vector<SCOTCH_Num> sc_peritab(N);
            SCOTCH_stratExit (&stradat);
            SCOTCH_Ordering  ordering;

            SCOTCH_dgraphCorderInit (&grafdat,
                &ordering,
                &sc_permtab[0],
                &sc_peritab[0],
                NULL,
                NULL,
                NULL
                );


            logfileptr->OFS()<<sc_permtab<<endl;
            logfileptr->OFS()<<sc_peritab<<endl;



            //SCOTCH_dgraphOrderGather (&grafdat, &ordedat, &ordering);
            if (iam == 0) {
              SCOTCH_dgraphOrderGather (&grafdat, &ordedat, &ordering);
            }
            else {     
              SCOTCH_dgraphOrderGather (&grafdat, &ordedat, NULL);
            }




            SCOTCH_dgraphCorderExit( &grafdat, &ordering );
            SCOTCH_dgraphOrderExit (&grafdat, &ordedat);
            SCOTCH_dgraphExit (&grafdat);


            if(iam==0 && baseval!=1){
              //switch everything to 1 based
              for(int col=0; col<N;++col){ invp[col] = sc_permtab[col] + (1-baseval);}
            }


            if(typeid(SCOTCH_Num) != typeid(Ptr)){
              delete [] pcolptr;
            }

            if(typeid(SCOTCH_Num) != typeid(Idx)){
              delete [] prowind;
            }













#else
            SYMPACK::vector<SCOTCH_Num> sc_permtab(vertlocnbr);


            if(SCOTCH_dgraphOrderPerm(&grafdat, &ordedat, &sc_permtab[0])!=0){
              throw std::logic_error( "Error in SCOTCH_dgraphOrderPerm\n" );
            }
            SCOTCH_dgraphOrderExit (&grafdat, &ordedat);
            SCOTCH_stratExit (&stradat);
            SCOTCH_dgraphExit (&grafdat);


            if(typeid(SCOTCH_Num) != typeid(Ptr)){
              delete [] pcolptr;
            }

            if(typeid(SCOTCH_Num) != typeid(Idx)){
              delete [] prowind;
            }


            //logfileptr->OFS()<<"Order: "<<sc_perm<<endl;

            //compute displs
            SYMPACK::vector<int> mpidispls(ndomains,0);
            SYMPACK::vector<int> mpisizes(ndomains,0);
            for(int p = 1;p<=ndomains;++p){
              mpisizes[p-1] = (vtxdist[p] - vtxdist[p-1])*sizeof(SCOTCH_Num);
              mpidispls[p-1] = (vtxdist[p-1]-baseval)*sizeof(SCOTCH_Num);
            }


            SCOTCH_Num * rbuffer;
            if(iam==0){
              if(typeid(SCOTCH_Num) != typeid(Int)){
                rbuffer = new SCOTCH_Num[N];
              }
              else{
                rbuffer = (SCOTCH_Num*)&invp[0];
              }
            }
            MPI_Gatherv(&sc_perm[0],sc_perm.size()*sizeof(SCOTCH_Num),MPI_BYTE,&rbuffer[0],&mpisizes[0],&mpidispls[0],MPI_BYTE,0,ndcomm);

            if(iam==0){

              //switch everything to 1 based
              if((void*)rbuffer!=(void*)&invp[0]){
                for(int col=0; col<N;++col){ invp[col] = (Int)rbuffer[col] + (1-baseval);}
                delete [] rbuffer;
              }
              else if(baseval!=1){
                for(int col=0; col<N;++col){ invp[col]+=(1-baseval);}
              }

              //logfileptr->OFS()<<"Full Order: "<<endl;
              //for(SCOTCH_Num i =0;i<N;++i){
              //  logfileptr->OFS()<<invp[i]<<" ";
              //  //logfileptr->OFS()<<invp[i]<<" ";
              //}
              //logfileptr->OFS()<<endl;
            }
#endif

          }

          MPI_Comm_free(&ndcomm);
          vtxdist.clear();

          // broadcast invp
          MPI_Bcast(&invp[0],N*sizeof(Int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
          //recompute perm
          perm.resize(N);
          for(Int i = 1; i <=N; ++i){
            Int node = invp[i-1];
            perm[node-1] = i;
          }

          logfileptr->OFS()<<perm<<endl;
          logfileptr->OFS()<<invp<<endl;


        }





        void Ordering::PTSCOTCH(){
          //assert(pStructure!=NULL);

          logfileptr->OFS()<<"PTSCOTCH used"<<endl;
          if(iam==0){cout<<"PTSCOTCH used NOT WORKING"<<endl;}


          if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
            throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call PTSCOTCH\n" );
          }


          SCOTCH_Num baseval = 1;
          SCOTCH_Num N = pStructure->size; 
          invp.resize(N);
          //      SCOTCH_Num ndomains = np;
          SCOTCH_Num ndomains = (SCOTCH_Num)pow(2.0,std::floor(std::log2(np)));
          //make sure every one have an element
          while(((double)N/(double)ndomains)<1.0){ndomains /= 2;}
          //      while(((double)N/(double)ndomains)<1.0){ndomains--;}

          MPI_Comm ndcomm;
          int mpirank;
          int color = iam < ndomains;
          MPI_Comm_split(CommEnv_->MPI_GetComm(),color,iam,&ndcomm);
          MPI_Comm_rank(ndcomm,&mpirank);
          SYMPACK::vector<SCOTCH_Num> vtxdist;
          SYMPACK::vector<SCOTCH_Num> tmpXadj;
          SYMPACK::vector<SCOTCH_Num> tmpAdj;
          SCOTCH_Num localN;
          if(iam<ndomains){
            assert(mpirank==iam);
            vtxdist.resize(ndomains+1);

            localN = N/ndomains;
            Idx fc = (iam)*localN+1;
            //build vtxdist SYMPACK::vector
            for(SCOTCH_Num i = 0; i<ndomains;++i){
              vtxdist[i] = i*localN+baseval; 
            } 
            vtxdist[ndomains] = N+baseval;

            if(iam==ndomains-1){
              localN = N - (ndomains-1)*localN;
            }



            //logfileptr->OFS()<<"vtxdist: "<<vtxdist<<endl;
            //logfileptr->OFS()<<"ndomains = "<<ndomains<<endl;

            //build local colptr and count nnz in local rowind
            tmpXadj.resize(localN+1);
            Ptr localNNZ = 0;
            tmpXadj[0] = 1;
            for(Idx col=1; col<=localN;++col){
              SCOTCH_Num column = col + fc -1;
              Ptr colbeg = pStructure->expColptr[column-1];
              Ptr colend = pStructure->expColptr[column]-1;
              //remove self from adjacency list
              tmpXadj[col] = tmpXadj[col-1] + (colend - colbeg + 1) - 1;
              assert(tmpXadj[col]>=tmpXadj[col-1]);
              localNNZ+=(colend-colbeg+1)-1;
            }

            //build local rowind, discarding self
            tmpAdj.resize(localNNZ);
            SCOTCH_Num pos =0;
            for(Idx col=1; col<=localN;++col){
              SCOTCH_Num column = col + fc -1;
              Ptr colbeg = pStructure->expColptr[column-1];
              Ptr colend = pStructure->expColptr[column]-1;
              for(Ptr j=colbeg; j<=colend;++j){
                SCOTCH_Num row = pStructure->expRowind[j-1];
                if(row!=column){
                  tmpAdj[pos++] = row;
                }
              }
            }




            MPI_Barrier(ndcomm);

            assert(localNNZ==tmpAdj.size());

            //switch everything to 0 based
            if(baseval!=1){
              for(Idx col=0; col<localN+1;++col){ tmpXadj[col]+=(baseval-1);}
              for(Idx col=0; col<localNNZ;++col){ tmpAdj[col]+=(baseval-1);}
            }
            //        logfileptr->OFS()<<"XAdj: "<<tmpXadj<<endl;
            //        logfileptr->OFS()<<"Adj: "<<tmpAdj<<endl;

            if(0)
            {
              SCOTCH_Num vertlocnum=1628;
              SCOTCH_Num vertglbnum=vertlocnum+vtxdist[iam]-baseval;
              SYMPACK::vector<SCOTCH_Num> toPrint;
              logfileptr->OFS()<<"colptr of "<<vertlocnum<<"("<<vertglbnum<<") "<<" from "<<tmpXadj[vertlocnum-1]<<" to "<<tmpXadj[vertlocnum]<<endl;
              logfileptr->OFS()<<"rowind of "<<vertlocnum<<"("<<vertglbnum<<") "<<": ";
              for(SCOTCH_Num edgelocnum = tmpXadj[vertlocnum-1]; edgelocnum <tmpXadj[vertlocnum]; edgelocnum++){
                SCOTCH_Num vert =tmpAdj[edgelocnum-1];
                logfileptr->OFS()<<vert<<" ";

                if(vert<vtxdist[iam+1] && vert>=vtxdist[iam]){
                  toPrint.push_back(vert);
                }

              }
              logfileptr->OFS()<<endl;


              for(int i =0;i<toPrint.size();i++){
                vertlocnum=toPrint[i]-vtxdist[iam]+baseval;
                vertglbnum=vertlocnum+vtxdist[iam]-baseval;
                logfileptr->OFS()<<"colptr of "<<vertlocnum<<"("<<vertglbnum<<") "<<" from "<<tmpXadj[vertlocnum-1]<<" to "<<tmpXadj[vertlocnum]<<endl;
                logfileptr->OFS()<<"rowind of "<<vertlocnum<<"("<<vertglbnum<<") "<<": ";
                for(SCOTCH_Num edgelocnum = tmpXadj[vertlocnum-1]; edgelocnum <tmpXadj[vertlocnum]; edgelocnum++){
                  SCOTCH_Num vert =tmpAdj[edgelocnum-1];
                  logfileptr->OFS()<<vert<<" ";
                }
                logfileptr->OFS()<<endl;


              }



            }


            //}
            //
            //    //clear up some memory
            //    pStructure->ClearExpandedSymmetric();
            //
            //      if(iam<ndomains){
            SCOTCH_Num options[3];
            options[0] = 0;

            SCOTCH_Num numflag = baseval;

            int npnd;
            MPI_Comm_size (ndcomm, &npnd);
            logfileptr->OFS()<<"PROC ND: "<<npnd<<endl;

            SCOTCH_Dgraph       grafdat;                    /* Scotch distributed graph object to SCOTCH_Numerface with libScotch    */
            SCOTCH_Dordering    ordedat;                    /* Scotch distributed ordering object to SCOTCH_Numerface with libScotch */
            SCOTCH_Strat        stradat;
            SCOTCH_Num          vertlocnbr;
            SCOTCH_Num          edgelocnbr;

            if(SCOTCH_dgraphInit (&grafdat, ndcomm) != 0){
              abort;
            }

            vertlocnbr = vtxdist[iam + 1] - vtxdist[iam];
            edgelocnbr = tmpXadj[vertlocnbr] - baseval;

            //        SCOTCH_Num maxvertlocnbr = vertlocnbr;
            //        for(SCOTCH_Num i = 0; i<ndomains;++i){
            //          maxvertlocnbr = std::max(maxvertlocnbr,vtxdist[i+1]-vtxdist[i]);
            //        }


            //logfileptr->OFS()<<maxvertlocnbr<<endl;
            //logfileptr->OFS()<<vertlocnbr<<" vs "<<tmpXadj.size()<<endl;
            //logfileptr->OFS()<<edgelocnbr<<" vs "<<tmpAdj.size()<<endl;

            logfileptr->OFS()<<"Entering SCOTCH_dgraphBuild"<<endl;
            if (SCOTCH_dgraphBuild (&grafdat, baseval,
                  vertlocnbr, vertlocnbr, &tmpXadj[0], &tmpXadj[0]+1, NULL, NULL,
                  edgelocnbr, edgelocnbr, &tmpAdj[0], NULL, NULL) == 0) {
              SCOTCH_stratInit (&stradat);
              logfileptr->OFS()<<"Entering SCOTCH_dgraphCheck"<<endl;
              assert(SCOTCH_dgraphCheck (&grafdat) == 0);
              {
                logfileptr->OFS()<<"Entering SCOTCH_dgraphOrderInit"<<endl;
                if (SCOTCH_dgraphOrderInit (&grafdat, &ordedat) == 0) {
                  logfileptr->OFS()<<"Entering SCOTCH_dgraphOrderCompute"<<endl;
                  SCOTCH_dgraphOrderCompute (&grafdat, &ordedat, &stradat);
                  SYMPACK::vector<SCOTCH_Num> sc_perm(invp.size());
                  //for(SCOTCH_Num i=0;i<sc_perm.size();i++){ sc_perm[i] = (SCOTCH_Num)invp[i]; }
                  logfileptr->OFS()<<"Entering SCOTCH_dgraphOrderPerm"<<endl;
                  SCOTCH_dgraphOrderPerm(&grafdat, &ordedat, &sc_perm[0]);
                  for(SCOTCH_Num i=0;i<sc_perm.size();i++){ invp[i] = sc_perm[i]; }
                  SCOTCH_dgraphOrderExit (&grafdat, &ordedat);
                }
              }
              SCOTCH_stratExit (&stradat);
            }
            SCOTCH_dgraphExit (&grafdat);


            logfileptr->OFS()<<"Order: "<<endl;
            for(SCOTCH_Num i =0;i<localN;++i){
              //logfileptr->OFS()<<invp[i]<<" ";
              logfileptr->OFS()<<invp[i]<<" ";
            }
            logfileptr->OFS()<<endl;


            //compute displs
            SYMPACK::vector<int> mpidispls(ndomains,0);
            SYMPACK::vector<int> mpisizes(ndomains,0);
            for(int p = 1;p<=ndomains;++p){
              mpisizes[p-1] = (vtxdist[p] - vtxdist[p-1])*sizeof(Int);
              mpidispls[p-1] = (vtxdist[p-1]-baseval)*sizeof(Int);
            }

            //gather on the root
            if(iam==0){
              MPI_Gatherv(MPI_IN_PLACE,mpisizes[iam],MPI_BYTE,&invp[0],&mpisizes[0],&mpidispls[0],MPI_BYTE,0,ndcomm);
            }
            else{
              MPI_Gatherv(&invp[0],mpisizes[iam],MPI_BYTE,&invp[0],&mpisizes[0],&mpidispls[0],MPI_BYTE,0,ndcomm);
            }
            //    MPI_Gatherv(&invp[0],mpisizes[iam],MPI_BYTE,&invp[0],&mpisizes[0],&mpidispls[0],MPI_BYTE,0,comm);

            if(iam==0){
              //switch everything to 1 based
              if(baseval!=1){
                for(int col=0; col<N;++col){ invp[col]+=(1-baseval);}
              }

              logfileptr->OFS()<<"Full Order: "<<endl;
              for(SCOTCH_Num i =0;i<N;++i){
                logfileptr->OFS()<<invp[i]<<" ";
                //logfileptr->OFS()<<invp[i]<<" ";
              }
              logfileptr->OFS()<<endl;
            }


        }

        MPI_Comm_free(&ndcomm);

        vtxdist.clear();
        tmpXadj.clear();
        tmpAdj.clear();

        // broadcast invp
        MPI_Bcast(&invp[0],N*sizeof(Int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
        perm.resize(N);
        for(Int i = 1; i <=N; ++i){
          Int node = invp[i-1];
          perm[node-1] = i;
        }

        //  // broadcast perm
        //  MPI_Bcast(&perm[0],N*sizeof(SCOTCH_Num),MPI_BYTE,0,comm);
        //
        //  for(Int i = 1; i <= perm.size(); ++i){
        //    Int node = perm[i-1];
        //    invp[node-1] = i;
        //  }
        logfileptr->OFS()<<"Matrix clearing"<<endl;
        pStructure->ClearExpandedSymmetric();
        logfileptr->OFS()<<"Matrix reexpanding"<<endl;
        pStructure->ExpandSymmetric();
        }
#endif


        void Ordering::SetCommEnvironment(CommEnvironment * CommEnv){
          CommEnv_=CommEnv;
        }

        void Ordering::Compose(SYMPACK::vector<Int> & invp2){
          //Compose the two permutations
          for(Int i = 1; i <= invp.size(); ++i){
            Int interm = invp[i-1];
            invp[i-1] = invp2[interm-1];
          }
          for(Int i = 1; i <= invp.size(); ++i){
            Int node = invp[i-1];
            perm[node-1] = i;
          }
        }

        }
