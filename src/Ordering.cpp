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
#include <metis.h>
#endif
#ifdef USE_PARMETIS
#include <typeinfo>
#include <metis.h>
#endif

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#pragma message  "IDXTYPEWIDTH=" STR(IDXTYPEWIDTH)

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


  void Ordering::MMD(const SparseMatrixGraph & g){
    logfileptr->OFS()<<"MMD used"<<endl;
    if(iam==0){cout<<"MMD used"<<endl;}

    if(iam == 0 && (!g.IsExpanded() || g.keepDiag==1) ){
      throw std::logic_error( "SparseMatrixGraph must be expanded and not including the diagonal in order to call MMD\n" );
    }


    bool isSameInt = typeid(MMDInt) == typeid(Int);
    bool isSameIdx = typeid(MMDInt) == typeid(Idx);
    bool isSamePtr = typeid(MMDInt) == typeid(Ptr);


    MMDInt N = g.VertexCount();
    MPI_Bcast(&N,sizeof(MMDInt),MPI_BYTE,0,CommEnv_->MPI_GetComm());

    invp.resize(N);

    if(iam==0){
      MMDInt iwsiz = 4*N;
      SYMPACK::vector<MMDInt> iwork (iwsiz);

      MMDInt * MMDInvp;      
      if(!isSameInt){
        MMDInvp = new MMDInt[N];
      }
      else{
        MMDInvp = (MMDInt*)&invp[0];
      }

      MMDInt * MMDperm;      
      if(!isSameInt){
        MMDperm = new MMDInt[N];
      }
      else{
        perm.resize(N);
        MMDperm = (MMDInt*)&perm[0];
      }


      MMDInt * prowind = NULL;
      if(!isSameIdx || g.baseval!=1){
        prowind = new MMDInt[g.EdgeCount()];
        for(Ptr i = 0; i<g.rowind.size();i++){ prowind[i] = (MMDInt)g.rowind[i];}
      }
      else{
        prowind = (MMDInt*)&g.rowind[0];
      }

      MMDInt * pcolptr = NULL;
      if(!isSamePtr || g.baseval!=1){
        pcolptr = new MMDInt[g.VertexCount()+1];
        for(Ptr i = 0; i<g.colptr.size();i++){ pcolptr[i] = (MMDInt)g.colptr[i];}
      }
      else{
        pcolptr = (MMDInt*)&g.colptr[0];
      }


      MMDInt nadj = g.EdgeCount();
      MMDInt nofsub =0;
      MMDInt iflag =0;
      FORTRAN(ordmmd)( &N , &nadj , pcolptr, prowind, 
          MMDInvp , MMDperm , &iwsiz , &iwork[0] , &nofsub, &iflag ) ;


      assert(iflag == 0);


      if(!isSameInt){ 
        //switch everything to 1 based
        for(int col=0; col<N;++col){ invp[col] = MMDInvp[col];}
      }

      if(!isSameInt){ 
        delete [] MMDperm;
      }

      if(!isSameInt){ 
        delete [] MMDInvp;
      }
      if(!isSamePtr || g.baseval!=1){
        delete [] pcolptr;
      }

      if(!isSameIdx || g.baseval!=1){
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













  void Ordering::NDBOX(Int size){
    invp.resize(size);
    if(iam==0){
      Int k = std::ceil(std::pow(size,1.0/3.0));
      logfileptr->OFS()<<"BOX K = "<<k<<endl;
      Int iflag =1;

      {
        SYMPACK::vector<Int> stack(k*k>25?invp.size():2*invp.size());
        Int tmp = stack.size();
        FORTRAN(boxnd)( &k , &k, &k, &invp[0], &stack[0],&tmp, &iflag);
      }

    }

    // broadcast invp
    MPI_Bcast(&invp[0],size*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
    perm.resize(size);
    for(Int i = 1; i <=size; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }
  }

  void Ordering::NDGRID(Int size){
    invp.resize(size);
    if(iam==0){
      Int k = std::ceil(std::pow(size,1.0/2.0));
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
    MPI_Bcast(&invp[0],size*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
    perm.resize(size);
    for(Int i = 1; i <=size; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }

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

      AMDInt nadj = pStructure->expRowind.size();
      AMDInt nofsub =0;
      AMDInt iflag =0;
      AMDInt iwsiz = 4*pStructure->size;
      SYMPACK::vector<AMDInt> iwork (iwsiz);
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

  void Ordering::AMD(const SparseMatrixGraph & g){
    logfileptr->OFS()<<"AMD used"<<endl;
    if(iam==0){cout<<"AMD used"<<endl;}

    if(iam == 0 && (!g.IsExpanded() || g.keepDiag==1) ){
      throw std::logic_error( "SparseMatrixGraph must be expanded and not including the diagonal in order to call AMD\n" );
    }


    bool isSameInt = typeid(AMDInt) == typeid(Int);
    bool isSameIdx = typeid(AMDInt) == typeid(Idx);
    bool isSamePtr = typeid(AMDInt) == typeid(Ptr);


    AMDInt N = g.VertexCount();
    MPI_Bcast(&N,sizeof(AMDInt),MPI_BYTE,0,CommEnv_->MPI_GetComm());
    invp.resize(N);

    if(iam==0){
      AMDInt * AMDInvp;      
      if(!isSameInt){
        AMDInvp = new AMDInt[N];
      }
      else{
        AMDInvp = (AMDInt*)&invp[0];
      }

      AMDInt * AMDperm;      
      if(!isSameInt){
        AMDperm = new AMDInt[N];
      }
      else{
        perm.resize(N);
        AMDperm = (AMDInt*)&perm[0];
      }


      AMDInt * prowind = new AMDInt[g.EdgeCount() + N];
      for(Ptr i = 0; i<g.rowind.size();i++){ prowind[i] = (AMDInt)g.rowind[i];}

      AMDInt * pcolptr = NULL;
      if(!isSamePtr || g.baseval!=1){
        pcolptr = new AMDInt[g.VertexCount()+1];
        for(Ptr i = 0; i<g.colptr.size();i++){ pcolptr[i] = (AMDInt)g.colptr[i];}
      }
      else{
        pcolptr = (AMDInt*)&g.colptr[0];
      }

      AMDInt IOVFLO = std::numeric_limits<AMDInt>::max();//2147483647;
      AMDInt NCMPA;

      AMDInt nadj = g.EdgeCount();
      AMDInt nofsub =0;
      AMDInt iflag =0;
      AMDInt iwsiz = 4*N;
      SYMPACK::vector<AMDInt> iwork (iwsiz);
      SYMPACK::vector<AMDInt> VTXDEG(N);
      SYMPACK::vector<AMDInt> QSIZE(N);
      SYMPACK::vector<AMDInt> ECFORW(N);
      SYMPACK::vector<AMDInt> MARKER(N);
      SYMPACK::vector<AMDInt> NVTXS(N+1);
      for(AMDInt i=0;i<N;++i){
        NVTXS[i] = pcolptr[i+1]-pcolptr[i];
      }


      AMDInt IWLEN = pcolptr[N-1] + NVTXS[N-1] + N - 1 ;
      AMDInt PFREE = pcolptr[N-1] + NVTXS[N-1];


      FORTRAN(amdbar)( &N , pcolptr , prowind , &NVTXS[0] ,&IWLEN ,&PFREE ,  &QSIZE[0],&ECFORW[0] , &perm[0], &iwork[0], AMDInvp,&VTXDEG[0], &NCMPA , &MARKER[0] , &IOVFLO );

      assert(iflag == 0);

      if(!isSameInt){ 
        //switch everything to 1 based
        for(int col=0; col<N;++col){ invp[col] = AMDInvp[col];}
      }

      if(!isSameInt){ 
        delete [] AMDperm;
      }

      if(!isSameInt){ 
        delete [] AMDInvp;
      }
      if(!isSamePtr || g.baseval!=1){
        delete [] pcolptr;
      }

      delete [] prowind;
    }
    // broadcast invp
    MPI_Bcast(&invp[0],N*sizeof(Int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
    perm.resize(N);
    for(Int i = 1; i <=N; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }
  }

#ifdef USE_METIS
  void Ordering::METIS(const SparseMatrixGraph & g){
    logfileptr->OFS()<<"METIS used"<<endl;
    if(iam==0){cout<<"METIS used"<<endl;}

    if(iam == 0 && (!g.IsExpanded() || g.keepDiag==1) ){
      throw std::logic_error( "SparseMatrixGraph must be expanded and not including the diagonal in order to call METIS\n" );
    }


    idx_t baseval = g.baseval;
    idx_t N = (idx_t)g.size; 
    MPI_Bcast(&N,sizeof(idx_t),MPI_BYTE,0,CommEnv_->MPI_GetComm());

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
        void Ordering::SCOTCH(const SparseMatrixGraph & g){
          logfileptr->OFS()<<"SCOTCH used"<<endl;
          if(iam==0){cout<<"SCOTCH used"<<endl;}

          if(iam == 0 && (!g.IsExpanded() || g.keepDiag==1) ){
            throw std::logic_error( "SparseMatrixGraph must be expanded and not including the diagonal in order to call SCOTCH\n" );
          }


          SCOTCH_Num baseval = g.baseval;
          SCOTCH_Num N = g.size; 
          MPI_Bcast(&N,sizeof(SCOTCH_Num),MPI_BYTE,0,CommEnv_->MPI_GetComm());

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
