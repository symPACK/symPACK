#include "sympack/Ordering.hpp"
#include "sympack/utility.hpp"

#include <limits>
#include <typeinfo>

#ifdef USE_SCOTCH
#include <scotch.h>
#endif
#ifdef USE_PTSCOTCH
#include <ptscotch.h>
#endif

#ifdef USE_METIS
#include <metis.h>
#endif
#ifdef USE_PARMETIS
#include <metis.h>
#include <parmetis.h>
#endif

/* define pour l'affichage */
#define SCOTCH_STRAT_DIRECT                                             \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>120)?m{rat=0.8,"                                  \
  ""                        "vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}}|"                          \
  ""                      "m{rat=0.8,"                                  \
  ""                        "vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}};,"                         \
  ""      "ole=f{cmin=0,cmax=100000,frat=0.0},"                       \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>120)?(m{rat=0.8,"                                 \
  ""                         "vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}})|"                        \
  ""                        "m{rat=0.8,"                                \
  ""                          "vert=100,"                               \
  ""                          "low=h{pass=10},"                         \
  ""                          "asc=f{bal=0.2}};,"                       \
  ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \
  ""      "ose=g}}"

#define PTSCOTCH_STRAT_DIRECT                                           \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>120)?m{rat=0.8,"                                 \
  ""                        "vert=100,"                                \
  ""                        "low=h{pass=10},"                          \
  ""                        "asc=f{bal=0.2}}|"                         \
  ""                      "m{rat=0.8,"                                 \
  ""                        "vert=100,"                                \
  ""                        "low=h{pass=10},"                          \
  ""                        "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=0,cmax=100000,frat=0.0},"                         \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>120)?(m{type=h,"                                  \
  ""                         "rat=0.8,"                                 \
  ""                         "vert=100000,"                             \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=08.2}})|"                       \
  ""                       "m{type=h,"                                  \
  ""                         "rat=0.8,"                                 \
  ""                         "vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \
  ""      "ose=g}}"





namespace symPACK {

  extern "C" {
    void FORTRAN(ordmmd)( MMDInt * neqns , MMDInt * nadj  , MMDInt * xadj  ,
        MMDInt * adjncy, MMDInt * invp  , MMDInt * perm  , MMDInt * iwsiz ,
        MMDInt * iwork , MMDInt * nofsub, MMDInt * iflag);
  }


  extern "C" {
    void FORTRAN(amdbar) (AMDInt * N, AMDInt * PE, AMDInt * IW, AMDInt * LEN, AMDInt * IWLEN, AMDInt * PFREE, AMDInt * NV, AMDInt * NEXT, AMDInt *
        LAST, AMDInt * HEAD, AMDInt * ELEN, AMDInt * DEGREE, AMDInt * NCMPA, AMDInt * W, AMDInt * IOVFLO);
  }

  extern "C" {
    void FORTRAN(boxnd) (Int * P, Int * Q, Int * R, Int * IPERM, Int * WORK,Int * WORKSZ, Int * IERROR);
    void FORTRAN(gridnd) (Int * P, Int * Q, Int * IPERM, Int * WORK,Int * WORKSZ, Int * IERROR);
  }


  extern "C" {
    void FORTRAN(genrcm)(RCMInt *, RCMInt*, RCMInt *,RCMInt*,RCMInt*,RCMInt*);
  }

}

namespace symPACK{

  void Ordering::MMD(const SparseMatrixGraph & g,MPI_Comm comm){
    int iam =0;
    int np =1;
    MPI_Comm_rank(comm,&iam);
    MPI_Comm_size(comm,&np);

    logfileptr->OFS()<<"MMD used"<<std::endl;
    if(iam==0){symPACKOS<<"MMD used"<<std::endl;}

    if(iam == 0 && (!g.IsExpanded() || g.keepDiag==1) ){
      throw std::logic_error( "SparseMatrixGraph must be expanded and not including the diagonal in order to call MMD\n" );
    }


    bool isSameInt = typeid(MMDInt) == typeid(Int);
    bool isSameIdx = typeid(MMDInt) == typeid(Idx);
    bool isSamePtr = typeid(MMDInt) == typeid(Ptr);

    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(Int), MPI_BYTE, &type );
    MPI_Type_commit(&type);

    MMDInt N = g.VertexCount();


    if(iam==0){
      invp.resize(N);
      MMDInt iwsiz = 4*N;
      std::vector<MMDInt> iwork (iwsiz);

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


      MMDInt * prowind = nullptr;
      if(!isSameIdx || g.baseval!=1){
        prowind = new MMDInt[g.EdgeCount()];
        for(Ptr i = 0; i<g.rowind.size();i++){ prowind[i] = (MMDInt)g.rowind[i];}
      }
      else{
        prowind = (MMDInt*)&g.rowind[0];
      }

      MMDInt * pcolptr = nullptr;
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
      MPI_Bcast(&N,sizeof(MMDInt),MPI_BYTE,0,comm);
    }
    else{
      MPI_Bcast(&N,sizeof(MMDInt),MPI_BYTE,0,comm);
      invp.resize(N);
    }
    // broadcast invp
    MPI_Bcast(&invp[0],N,type,0,comm);
    perm.resize(N);
    for(Int i = 1; i <=N; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }

    MPI_Type_free(&type);




  }


  void Ordering::RCM(const SparseMatrixGraph & g,MPI_Comm comm){
    int iam =0;
    int np =1;
    MPI_Comm_rank(comm,&iam);
    MPI_Comm_size(comm,&np);

    logfileptr->OFS()<<"RCM used"<<std::endl;
    if(iam==0){symPACKOS<<"RCM used"<<std::endl;}

    if(iam == 0 && (!g.IsExpanded() || g.keepDiag==1) ){
      throw std::logic_error( "SparseMatrixGraph must be expanded and not including the diagonal in order to call RCM\n" );
    }


    bool isSameInt = typeid(RCMInt) == typeid(Int);
    bool isSameIdx = typeid(RCMInt) == typeid(Idx);
    bool isSamePtr = typeid(RCMInt) == typeid(Ptr);

    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(Int), MPI_BYTE, &type );
    MPI_Type_commit(&type);

    RCMInt N = g.VertexCount();


    if(iam==0){
      invp.resize(N);
      RCMInt iwsiz = 4*N;
      std::vector<RCMInt> iwork (iwsiz);

      RCMInt * RCMInvp;      
      if(!isSameInt){
        RCMInvp = new RCMInt[N];
      }
      else{
        RCMInvp = (RCMInt*)&invp[0];
      }

      RCMInt * RCMperm;      
      if(!isSameInt){
        RCMperm = new RCMInt[N];
      }
      else{
        perm.resize(N);
        RCMperm = (RCMInt*)&perm[0];
      }


      RCMInt * prowind = nullptr;
      if(!isSameIdx || g.baseval!=1){
        prowind = new RCMInt[g.EdgeCount()];
        for(Ptr i = 0; i<g.rowind.size();i++){ prowind[i] = (RCMInt)g.rowind[i];}
      }
      else{
        prowind = (RCMInt*)&g.rowind[0];
      }

      RCMInt * pcolptr = nullptr;
      if(!isSamePtr || g.baseval!=1){
        pcolptr = new RCMInt[g.VertexCount()+1];
        for(Ptr i = 0; i<g.colptr.size();i++){ pcolptr[i] = (RCMInt)g.colptr[i];}
      }
      else{
        pcolptr = (RCMInt*)&g.colptr[0];
      }

      std::vector<RCMInt>mark(2*N);
      std::vector<RCMInt>xls(2*N);
      FORTRAN(genrcm)(&N, pcolptr, prowind, RCMperm, mark.data(), xls.data() );
      for(Int i = 1; i <=N; ++i){
        Int node = RCMperm[i-1];
        RCMInvp[node-1] = i;
      }



      if(!isSameInt){ 
        //switch everything to 1 based
        for(int col=0; col<N;++col){ invp[col] = RCMInvp[col];}
      }

      if(!isSameInt){ 
        delete [] RCMperm;
      }

      if(!isSameInt){ 
        delete [] RCMInvp;
      }
      if(!isSamePtr || g.baseval!=1){
        delete [] pcolptr;
      }

      if(!isSameIdx || g.baseval!=1){
        delete [] prowind;
      }
      MPI_Bcast(&N,sizeof(RCMInt),MPI_BYTE,0,comm);
    }
    else{
      MPI_Bcast(&N,sizeof(RCMInt),MPI_BYTE,0,comm);
      invp.resize(N);
    }
    // broadcast invp
    MPI_Bcast(&invp[0],N,type,0,comm);
    perm.resize(N);
    for(Int i = 1; i <=N; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }
    MPI_Type_free(&type);
  }













  void Ordering::NDBOX(Int size, MPI_Comm comm){
    int iam =0;
    int np =1;
    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(Int), MPI_BYTE, &type );
    MPI_Type_commit(&type);
    MPI_Comm_rank(comm,&iam);
    MPI_Comm_size(comm,&np);
    invp.resize(size);
    if(iam==0){
      Int k = std::ceil(std::pow(size,1.0/3.0));
      logfileptr->OFS()<<"BOX K = "<<k<<std::endl;
      Int iflag =1;

      {
        std::vector<Int> stack(k*k>25?invp.size():2*invp.size());
        Int tmp = stack.size();
        FORTRAN(boxnd)( &k , &k, &k, &invp[0], &stack[0],&tmp, &iflag);
      }

    }

    // broadcast invp
    MPI_Bcast(&invp[0],size,type,0,comm);
    perm.resize(size);
    for(Int i = 1; i <=size; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }
    MPI_Type_free(&type);
  }

  void Ordering::NDGRID(Int size, MPI_Comm comm){
    int iam =0;
    int np =1;
    MPI_Comm_rank(comm,&iam);
    MPI_Comm_size(comm,&np);
    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(Int), MPI_BYTE, &type );
    MPI_Type_commit(&type);
    invp.resize(size);
    if(iam==0){
      Int k = std::ceil(std::pow(size,1.0/2.0));
      Int iflag =1;

      logfileptr->OFS()<<"GRID K = "<<k<<std::endl;

      {
        std::vector<Int> stack(k*k>25?invp.size():2*invp.size());
        Int tmp = stack.size();
        FORTRAN(gridnd)( &k , &k, &invp[0], &stack[0],&tmp, &iflag);
        assert(iflag==0);
      }
    }

    // broadcast invp
    MPI_Bcast(&invp[0],size,type,0,comm);
    perm.resize(size);
    for(Int i = 1; i <=size; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }

    MPI_Type_free(&type);
  }


  void Ordering::AMD(const SparseMatrixGraph & g, MPI_Comm comm){

    int iam =0;
    int np =1;
    MPI_Comm_rank(comm,&iam);
    MPI_Comm_size(comm,&np);

    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(Int), MPI_BYTE, &type );
    MPI_Type_commit(&type);
    logfileptr->OFS()<<"AMD used"<<std::endl;
    if(iam==0){symPACKOS<<"AMD used"<<std::endl;}

    if(iam == 0 && (!g.IsExpanded() || g.keepDiag==1) ){
      throw std::logic_error( "SparseMatrixGraph must be expanded and not including the diagonal in order to call AMD\n" );
    }


    bool isSameInt = typeid(AMDInt) == typeid(Int);
    bool isSameIdx = typeid(AMDInt) == typeid(Idx);
    bool isSamePtr = typeid(AMDInt) == typeid(Ptr);


    AMDInt N = g.VertexCount();

    if(iam==0){
      invp.resize(N);
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

      AMDInt * pcolptr = nullptr;
      if(!isSamePtr || g.baseval!=1){
        pcolptr = new AMDInt[g.VertexCount()+1];
        for(Ptr i = 0; i<g.colptr.size();i++){ pcolptr[i] = (AMDInt)g.colptr[i];}
      }
      else{
        pcolptr = (AMDInt*)&g.colptr[0];
      }

      AMDInt IOVFLO = std::numeric_limits<AMDInt>::max();
      AMDInt NCMPA;

      AMDInt nadj = g.EdgeCount();
      AMDInt nofsub =0;
      AMDInt iflag =0;
      AMDInt iwsiz = 4*N;
      std::vector<AMDInt> iwork (iwsiz);
      std::vector<AMDInt> VTXDEG(N);
      std::vector<AMDInt> QSIZE(N);
      std::vector<AMDInt> ECFORW(N);
      std::vector<AMDInt> MARKER(N);
      std::vector<AMDInt> NVTXS(N+1);
      for(AMDInt i=0;i<N;++i){
        NVTXS[i] = pcolptr[i+1]-pcolptr[i];
      }


      AMDInt IWLEN = pcolptr[N-1] + NVTXS[N-1] + N - 1 ;
      AMDInt PFREE = pcolptr[N-1] + NVTXS[N-1];


      FORTRAN(amdbar)( &N , pcolptr , prowind , &NVTXS[0] ,&IWLEN ,&PFREE ,  &QSIZE[0],&ECFORW[0] , AMDperm, &iwork[0], AMDInvp,&VTXDEG[0], &NCMPA , &MARKER[0] , &IOVFLO );

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
      MPI_Bcast(&N,sizeof(AMDInt),MPI_BYTE,0,comm);
    }
    else{
      MPI_Bcast(&N,sizeof(AMDInt),MPI_BYTE,0,comm);
      invp.resize(N);
    }
    // broadcast invp
    MPI_Bcast(&invp[0],N,type,0,comm);
    perm.resize(N);
    for(Int i = 1; i <=N; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }
    MPI_Type_free(&type);
  }

#ifdef USE_METIS
  void Ordering::METIS(const SparseMatrixGraph & g,MPI_Comm comm){

    int iam =0;
    int np =1;
    MPI_Comm_rank(comm,&iam);
    MPI_Comm_size(comm,&np);

    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(Int), MPI_BYTE, &type );
    MPI_Type_commit(&type);
    logfileptr->OFS()<<"METIS used"<<std::endl;
    if(iam==0){symPACKOS<<"METIS used"<<std::endl;}

    if(iam == 0 && (!g.IsExpanded() || g.keepDiag==1) ){
      throw std::logic_error( "SparseMatrixGraph must be expanded and not including the diagonal in order to call METIS\n" );
    }


    idx_t baseval = g.baseval;
    idx_t N = (idx_t)g.size; 

    bool isSameInt = typeid(idx_t) == typeid(Int);
    bool isSameIdx = typeid(idx_t) == typeid(Idx);
    bool isSamePtr = typeid(idx_t) == typeid(Ptr);



    if(iam==0){
      idx_t * iperm;      
      if(!isSameInt){
        iperm = new idx_t[N];
      }
      else{
        invp.resize(N);
        iperm = (idx_t*)&invp[0];
      }

      idx_t * mperm;      
      if(!isSameInt){
        mperm = new idx_t[N];
      }
      else{
        perm.resize(N);
        mperm = (idx_t*)&perm[0];
      }


      idx_t * prowind = nullptr;
      if(!isSameIdx){
        prowind = new idx_t[g.EdgeCount()];
        for(Ptr i = 0; i<g.rowind.size();i++){ prowind[i] = (idx_t)g.rowind[i];}
      }
      else{
        prowind = (idx_t*)&g.rowind[0];
      }

      idx_t * pcolptr = nullptr;
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

      METIS_NodeND(&N, pcolptr, prowind, nullptr,
          options, mperm, iperm);


      if(!isSameInt){ 
        delete [] mperm;
      }

      if(!isSameInt ||g.baseval!=1){ 
        if(invp.size()!=N){invp.resize(N);}
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
      MPI_Bcast(&N,sizeof(idx_t),MPI_BYTE,0,comm);
    }
    else{
      MPI_Bcast(&N,sizeof(idx_t),MPI_BYTE,0,comm);
      invp.resize(N);
    }
    // broadcast invp
    MPI_Bcast(&invp[0],N,type,0,comm);
    perm.resize(N);
    for(Int i = 1; i <=N; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }

    MPI_Type_free(&type);

  }







#endif


#ifdef USE_PARMETIS
  void Ordering::PARMETIS(const DistSparseMatrixGraph & g){

    int iam =0;
    int np =1;
    MPI_Comm_rank(g.comm,&iam);
    MPI_Comm_size(g.comm,&np);

    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(Int), MPI_BYTE, &type );
    MPI_Type_commit(&type);
    logfileptr->OFS()<<"PARMETIS used"<<std::endl;
    if(iam==0){symPACKOS<<"PARMETIS used"<<std::endl;}

    if(!g.IsExpanded() || g.keepDiag==1){
      throw std::logic_error( "DistSparseMatrixGraph must be expanded and not including the diagonal in order to call PARMETIS\n" );
    }


    idx_t baseval = g.baseval;
    idx_t N = g.size; 
    invp.resize(N);

    int mpirank;

    idx_t ndomains = np;
    NpOrdering = std::max(std::min(np,NpOrdering>0?NpOrdering:np),0);
    if(NpOrdering!=0){

      ndomains = NpOrdering;

      if(iam==0){symPACKOS<<"PARMETIS using "<<ndomains<< " / "<< np<<std::endl;}
      //make a copy
      DistSparseMatrixGraph ng = g;
      //create a new vertexDist, this is the important stuff
      Idx colPerProc = g.size / NpOrdering;
      std::vector<Idx> newVertexDist(np+1,colPerProc);
      newVertexDist[0]=1;
      std::partial_sum(newVertexDist.begin(),newVertexDist.end(),newVertexDist.begin());
      for(int i = NpOrdering;i<np+1;i++){ newVertexDist[i] = g.size+newVertexDist[0];}

      //redistribute ng
      ng.Redistribute(&newVertexDist[0]);

      MPI_Comm ndcomm;
      //TODO this needs to be experimented to distributed the load across multiple nodes
      int color = iam<ndomains;
      MPI_Comm_split(g.comm,color,iam,&ndcomm);

      MPI_Comm_rank(ndcomm,&mpirank);
      std::vector<idx_t> vtxdist;
      idx_t localN;
      if(iam<ndomains){
        assert(mpirank==iam);

        std::vector<idx_t> sizes(2*ndomains);
        vtxdist.resize(ng.vertexDist.size());
        for(int i = 0 ; i < ng.vertexDist.size(); i++){
          vtxdist[i] = (idx_t)ng.vertexDist[i];
        }

        idx_t * iperm = nullptr;
        if(typeid(idx_t) != typeid(Int)){
          iperm = new idx_t[N];
        }
        else{
          iperm = (idx_t*)&invp[0];
        }

        idx_t * pperm = nullptr;
        if(typeid(idx_t) != typeid(Int)){
          pperm = new idx_t[N];
        }
        else{
          perm.resize(N);
          pperm = (idx_t*)&perm[0];
        }

        idx_t * prowind = nullptr;
        if(typeid(idx_t) != typeid(Idx)){
          prowind = new idx_t[ng.LocalEdgeCount()];
          for(Ptr i = 0; i<ng.rowind.size();i++){ prowind[i] = (idx_t)ng.rowind[i];}
        }
        else{
          prowind = (idx_t*)&ng.rowind[0];
        }

        idx_t * pcolptr = nullptr;
        if(typeid(idx_t) != typeid(Ptr)){
          pcolptr = new idx_t[ng.LocalVertexCount()+1];
          for(Ptr i = 0; i<ng.colptr.size();i++){ pcolptr[i] = (idx_t)ng.colptr[i];}
        }
        else{
          pcolptr = (idx_t*)&ng.colptr[0];
        }


        idx_t options[3];
        options[0] = 0;
        idx_t numflag = ng.baseval;

        int npnd;
        MPI_Comm_size (ndcomm, &npnd);
        ParMETIS_V3_NodeND( &vtxdist[0], pcolptr , prowind, &numflag, &options[0], pperm, &sizes[0], &ndcomm );

        //compute displs
        std::vector<int> mpidispls(ndomains,0);
        std::vector<int> mpisizes(ndomains,0);
        for(int p = 1;p<=ndomains;++p){
          mpisizes[p-1] = (vtxdist[p] - vtxdist[p-1])*sizeof(idx_t);
          mpidispls[p-1] = (vtxdist[p-1]-baseval)*sizeof(idx_t);
        }

        //gather on the root
        MPI_Gatherv(pperm,mpisizes[iam],MPI_BYTE,iperm,&mpisizes[0],&mpidispls[0],MPI_BYTE,0,ndcomm);

        if(iam==0 && (ng.baseval!=1 || typeid(idx_t) != typeid(Int))){
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



    }
    else{
      MPI_Comm ndcomm;
      int color = iam<ndomains;
      MPI_Comm_split(g.comm,color,iam,&ndcomm);

      MPI_Comm_rank(ndcomm,&mpirank);
      std::vector<idx_t> vtxdist;
      idx_t localN;
      if(iam<ndomains){
        assert(mpirank==iam);

        std::vector<idx_t> sizes(2*ndomains);
        vtxdist.resize(g.vertexDist.size());
        for(int i = 0 ; i < g.vertexDist.size(); i++){
          vtxdist[i] = (idx_t)g.vertexDist[i];
        }


        idx_t * iperm = nullptr;
        if(typeid(idx_t) != typeid(Int)){
          iperm = new idx_t[N];
        }
        else{
          iperm = (idx_t*)&invp[0];
        }

        idx_t * pperm = nullptr;
        if(typeid(idx_t) != typeid(Int)){
          pperm = new idx_t[N];
        }
        else{
          perm.resize(N);
          pperm = (idx_t*)&perm[0];
        }

        idx_t * prowind = nullptr;
        if(typeid(idx_t) != typeid(Idx)){
          prowind = new idx_t[g.LocalEdgeCount()];
          for(Ptr i = 0; i<g.rowind.size();i++){ prowind[i] = (idx_t)g.rowind[i];}
        }
        else{
          prowind = (idx_t*)&g.rowind[0];
        }

        idx_t * pcolptr = nullptr;
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
        std::vector<int> mpidispls(ndomains,0);
        std::vector<int> mpisizes(ndomains,0);
        for(int p = 1;p<=ndomains;++p){
          mpisizes[p-1] = (vtxdist[p] - vtxdist[p-1])*sizeof(idx_t);
          mpidispls[p-1] = (vtxdist[p-1]-baseval)*sizeof(idx_t);
        }

        //gather on the root
        MPI_Gatherv(pperm,mpisizes[iam],MPI_BYTE,iperm,&mpisizes[0],&mpidispls[0],MPI_BYTE,0,ndcomm);

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

    }


    // broadcast invp
    MPI_Bcast(&invp[0],N,type,0,g.comm);
    //recompute perm
    perm.resize(N);
    for(Int i = 1; i <=N; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }

    MPI_Type_free(&type);

  }




#endif

#ifdef USE_SCOTCH
  void Ordering::SCOTCH(const SparseMatrixGraph & g,MPI_Comm comm){

    int iam =0;
    int np =1;
    MPI_Comm_rank(comm,&iam);
    MPI_Comm_size(comm,&np);

    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(Int), MPI_BYTE, &type );
    MPI_Type_commit(&type);
    logfileptr->OFS()<<"SCOTCH used"<<std::endl;
    if(iam==0){symPACKOS<<"SCOTCH used"<<std::endl;}

    if(iam == 0 && (!g.IsExpanded() || g.keepDiag==1) ){
      throw std::logic_error( "SparseMatrixGraph must be expanded and not including the diagonal in order to call SCOTCH\n" );
    }


    SCOTCH_Num baseval = g.baseval;
    SCOTCH_Num N = g.size; 

    bool isSameInt = typeid(SCOTCH_Num) == typeid(Int);
    bool isSameIdx = typeid(SCOTCH_Num) == typeid(Idx);
    bool isSamePtr = typeid(SCOTCH_Num) == typeid(Ptr);



    if(iam==0){
      invp.resize(N);
      SCOTCH_Num * permtab;      
      if(!isSameInt){
        permtab = new SCOTCH_Num[N];
      }
      else{
        permtab = (SCOTCH_Num*)&invp[0];
      }

      SCOTCH_Num * prowind = nullptr;
      if(!isSameIdx){
        prowind = new SCOTCH_Num[g.EdgeCount()];
        for(Ptr i = 0; i<g.rowind.size();i++){ prowind[i] = (SCOTCH_Num)g.rowind[i];}
      }
      else{
        prowind = (SCOTCH_Num*)&g.rowind[0];
      }

      SCOTCH_Num * pcolptr = nullptr;
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
              nullptr, 
              nullptr,
              nullptr,
              nnz, 
              prowind, 
              nullptr) == 0) {

          SCOTCH_stratInit (&stradat);

          char strat [550];
          sprintf(strat, SCOTCH_STRAT_DIRECT);
          SCOTCH_stratGraphOrder (&stradat, strat);

#ifdef SCOTCH_DEBUG_ALL
          if (SCOTCH_graphCheck (&grafdat) == 0)        /* TRICK: next instruction called only if graph is consistent */
#endif /* SCOTCH_DEBUG_ALL */
          {
            if (SCOTCH_graphOrderInit (
                  &grafdat,
                  &ordedat,
                  permtab,
                  nullptr,
                  /*nb de supernoeud*/nullptr,
                  /*ptr vers rank tab : xsuper*/ nullptr,
                  /*tree tab: parent structure*/ nullptr) == 0) {
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
      MPI_Bcast(&N,sizeof(SCOTCH_Num),MPI_BYTE,0,comm);
    }
    else{
      MPI_Bcast(&N,sizeof(SCOTCH_Num),MPI_BYTE,0,comm);
      invp.resize(N);
    }
    // broadcast invp
    MPI_Bcast(&invp[0],N,type,0,comm);
    perm.resize(N);
    for(Int i = 1; i <=N; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }

    MPI_Type_free(&type);

  }







#endif



#ifdef USE_PTSCOTCH
  void Ordering::PTSCOTCH(const DistSparseMatrixGraph & g){
    int iam =0;
    int np =1;
    MPI_Comm_rank(g.comm,&iam);
    MPI_Comm_size(g.comm,&np);

    MPI_Datatype type;
    MPI_Type_contiguous( sizeof(Int), MPI_BYTE, &type );
    MPI_Type_commit(&type);
    logfileptr->OFS()<<"PTSCOTCH used"<<std::endl;
    if(iam==0){symPACKOS<<"PTSCOTCH used"<<std::endl;}

    if(!g.IsExpanded() || g.keepDiag==1){
      throw std::logic_error( "DistSparseMatrixGraph must be expanded and not including the diagonal in order to call PTSCOTCH\n" );
    }


    SCOTCH_Num baseval = g.baseval;
    SCOTCH_Num N = g.size; 
    invp.resize(N);

    int mpirank;

    NpOrdering = std::max(std::min(np,NpOrdering>0?NpOrdering:np),0);
    SCOTCH_Num ndomains = NpOrdering;
    MPI_Comm ndcomm;
    int color = iam < ndomains;
    MPI_Comm_split(g.comm,color,iam,&ndcomm);

    MPI_Comm_rank(ndcomm,&mpirank);
    std::vector<SCOTCH_Num> vtxdist;
    SCOTCH_Num localN;

    SCOTCH_Num * prowind = nullptr;
    SCOTCH_Num * pcolptr = nullptr;
    SCOTCH_Num          vertlocnbr;
    SCOTCH_Num          edgelocnbr;
    DistSparseMatrixGraph ng;
    if(ndomains!=np){
      //make a copy
      ng = g;
      //create a new vertexDist, this is the important stuff
      Idx colPerProc = g.size / ndomains;
      std::vector<Idx> newVertexDist(np+1,colPerProc);
      newVertexDist[0]=1;
      std::partial_sum(newVertexDist.begin(),newVertexDist.end(),newVertexDist.begin());
      for(int i = ndomains;i<np+1;i++){ newVertexDist[i] = g.size+newVertexDist[0];}

      //redistribute ng
      ng.Redistribute(&newVertexDist[0]);

      vtxdist.resize(ng.vertexDist.size());
      for(int i = 0 ; i < ng.vertexDist.size(); i++){
        vtxdist[i] = (SCOTCH_Num)ng.vertexDist[i];
      }

      if(typeid(SCOTCH_Num) != typeid(Idx)){
        prowind = new SCOTCH_Num[ng.LocalEdgeCount()];
        for(Ptr i = 0; i<ng.rowind.size();i++){ prowind[i] = (SCOTCH_Num)ng.rowind[i];}
      }
      else{
        prowind = (SCOTCH_Num*)&ng.rowind[0];
      }

      if(typeid(SCOTCH_Num) != typeid(Ptr)){
        pcolptr = new SCOTCH_Num[ng.LocalVertexCount()+1];
        for(Ptr i = 0; i<ng.colptr.size();i++){ pcolptr[i] = (SCOTCH_Num)ng.colptr[i];}
      }
      else{
        pcolptr = (SCOTCH_Num*)&ng.colptr[0];
      }

      vertlocnbr = ng.LocalVertexCount();
      edgelocnbr = ng.LocalEdgeCount();

    }
    else{
      vtxdist.resize(g.vertexDist.size());
      for(int i = 0 ; i < g.vertexDist.size(); i++){
        vtxdist[i] = (SCOTCH_Num)g.vertexDist[i];
      }

      if(typeid(SCOTCH_Num) != typeid(Idx)){
        prowind = new SCOTCH_Num[g.LocalEdgeCount()];
        for(Ptr i = 0; i<g.rowind.size();i++){ prowind[i] = (SCOTCH_Num)g.rowind[i];}
      }
      else{
        prowind = (SCOTCH_Num*)&g.rowind[0];
      }

      if(typeid(SCOTCH_Num) != typeid(Ptr)){
        pcolptr = new SCOTCH_Num[g.LocalVertexCount()+1];
        for(Ptr i = 0; i<g.colptr.size();i++){ pcolptr[i] = (SCOTCH_Num)g.colptr[i];}
      }
      else{
        pcolptr = (SCOTCH_Num*)&g.colptr[0];
      }

      vertlocnbr = g.LocalVertexCount();
      edgelocnbr = g.LocalEdgeCount();
    }

    if(iam<ndomains){
      assert(mpirank==iam);
      SCOTCH_Num options[3];
      options[0] = 0;
      SCOTCH_Num numflag = baseval;

      int npnd;
      MPI_Comm_size (ndcomm, &npnd);

      SCOTCH_Dgraph       grafdat;                    /* Scotch distributed graph object to SCOTCH_Numerface with libScotch    */
      SCOTCH_Dordering    ordedat;                    /* Scotch distributed ordering object to SCOTCH_Numerface with libScotch */
      SCOTCH_Strat        stradat;


      logfileptr->OFS()<<"SCOTCH_dgraphInit"<<std::endl;
      if(SCOTCH_dgraphInit (&grafdat, ndcomm) != 0){
        throw std::logic_error( "Error in SCOTCH_dgraphInit\n" );
      }

      logfileptr->OFS()<<"SCOTCH_dgraphBuild"<<std::endl;
      if (SCOTCH_dgraphBuild (&grafdat, 
            baseval,
            vertlocnbr, 
            vertlocnbr, 
            pcolptr, 
            nullptr,
            nullptr, 
            nullptr,
            edgelocnbr,
            edgelocnbr,
            prowind, 
            nullptr,
            nullptr
            ) != 0) {
        throw std::logic_error( "Error in SCOTCH_dgraphBuild\n" );
      }

#ifdef _DEBUG_
      if(SCOTCH_dgraphCheck (&grafdat) != 0){
        throw std::logic_error( "Error in SCOTCH_dgraphCheck\n" );
      }
#endif

      logfileptr->OFS()<<"SCOTCH_stratInit"<<std::endl;
      if(SCOTCH_stratInit (&stradat)!= 0){
        throw std::logic_error( "Error in SCOTCH_stratInit\n" );
      }
      logfileptr->OFS()<<"SCOTCH_orderInit"<<std::endl;
      if (SCOTCH_dgraphOrderInit (&grafdat, &ordedat) != 0) {
        throw std::logic_error( "Error in SCOTCH_dgraphOrderInit\n" );
      }

      logfileptr->OFS()<<"SCOTCH_orderCompute"<<std::endl;
      if(SCOTCH_dgraphOrderCompute (&grafdat, &ordedat, &stradat)!=0){
        throw std::logic_error( "Error in SCOTCH_dgraphOrderCompute\n" );
      }




      std::vector<SCOTCH_Num> sc_permtab(N);
      std::vector<SCOTCH_Num> sc_peritab(N);
      SCOTCH_stratExit (&stradat);
      SCOTCH_Ordering  ordering;

      logfileptr->OFS()<<"SCOTCH_dgraphCorderInit"<<std::endl;
      SCOTCH_dgraphCorderInit (&grafdat,
          &ordering,
          &sc_permtab[0],
          &sc_peritab[0],
          nullptr,
          nullptr,
          nullptr
          );

      logfileptr->OFS()<<"SCOTCH_dgraphOrderGather"<<std::endl;
      SCOTCH_dgraphOrderGather (&grafdat, &ordedat, (iam==0?&ordering:nullptr));

      logfileptr->OFS()<<"SCOTCH_dgraphCorderExit"<<std::endl;
      SCOTCH_dgraphCorderExit( &grafdat, &ordering );
      logfileptr->OFS()<<"SCOTCH_dgraphOrderExit"<<std::endl;
      SCOTCH_dgraphOrderExit (&grafdat, &ordedat);
      logfileptr->OFS()<<"SCOTCH_dgraphExit"<<std::endl;
      SCOTCH_dgraphExit (&grafdat);


      if(iam==0){
        //switch everything to 1 based
        for(int col=0; col<N;++col){ invp[col] = sc_permtab[col] + (1-baseval);}
      }


      if(typeid(SCOTCH_Num) != typeid(Ptr)){
        delete [] pcolptr;
      }

      if(typeid(SCOTCH_Num) != typeid(Idx)){
        delete [] prowind;
      }

    }

    MPI_Comm_free(&ndcomm);
    vtxdist.clear();

    // broadcast invp
    MPI_Bcast(&invp[0],N,type,0,g.comm);
    //recompute perm
    perm.resize(N);
    for(Int i = 1; i <=N; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }

    MPI_Type_free(&type);

  }

#endif

  void Ordering::GetRelativeInvp(const std::vector<Int> & frominvp, std::vector<Int> & relinvp){
    assert(frominvp.size()==invp.size());
    relinvp.resize(invp.size());
    Int n = invp.size();
#pragma ivdep
    for(Int node = 1; node <= n; ++node){
      Int i = perm[node-1];
      Int interm = frominvp[i-1];
      relinvp[interm-1] = node;
    }
  }

  void Ordering::Compose(std::vector<Int> & invp2){
    //Compose the two permutations
    Int n = invp.size();
#pragma ivdep
    for(Int i = 1; i <= n; ++i){
      Int interm = invp[i-1];
      invp[i-1] = invp2[interm-1];
    }
#pragma ivdep
    for(Int i = 1; i <= n; ++i){
      Int node = invp[i-1];
      perm[node-1] = i;
    }
  }


}
