#include "sympack/Ordering.hpp"
#include "sympack/mmd.hpp"
#include "sympack/utility.hpp"


#ifdef USE_SCOTCH
#include <scotch.h>
#endif
#ifdef USE_PTSCOTCH
#include <ptscotch.h>
#endif

//#define USE_METIS_INC
//#ifdef USE_METIS_INC
//#include <metis.h>
//#define METIS_NOPTIONS 40
//#define METIS_NOPTIONS 40
//#endif


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
    invp.Resize(pStructure->size);
    perm.Resize(pStructure->size);
    for(Int i =0;i<invp.m();++i){
      perm[i]=i+1;
      invp[i]=i+1;
    }
}







void Ordering::MMD(){
    assert(pStructure!=NULL);

  logfileptr->OFS()<<"MMD used"<<endl;
  if(iam==0){cout<<"MMD used"<<endl;}

    if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
			throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call MMD\n" );
    }

    invp.Resize(pStructure->size);
    if(iam==0){
      Int iwsiz = 4*pStructure->size;
      std::vector<Int> iwork (iwsiz);

      Int nadj = pStructure->expRowind.m();
      Int nofsub =0;
      Int iflag =0;

      IntNumVec tmpXadj = pStructure->expColptr;
      IntNumVec tmpAdj(pStructure->expRowind.m()-pStructure->size);

      int pos = 0;
      for(int col=0; col<tmpXadj.m()-1;++col){
        for(int j=tmpXadj[col]; j<tmpXadj[col+1];++j){
          if( pStructure->expRowind[j-1]-1 != col){
            tmpAdj[pos++] = pStructure->expRowind[j-1];
          }
        }
      }

      int rm = 0;
      for(int col=0; col<tmpXadj.m();++col){
        tmpXadj[col]-=rm;
        rm++;
      }



      perm.Resize(pStructure->size);





      FORTRAN(ordmmd)( &pStructure->size , &nadj , tmpXadj.Data() , tmpAdj.Data(), 
          &invp[0] , &perm[0] , &iwsiz , &iwork[0] , &nofsub, &iflag ) ;


      assert(iflag == 0);

      //logfileptr->OFS()<<perm<<endl;
    }
      // broadcast invp
      Int N = invp.m();
      MPI_Bcast(&invp[0],N*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
      perm.Resize(pStructure->size);
      for(Int i = 1; i <=N; ++i){
        Int node = invp[i-1];
        perm[node-1] = i;
      }

      //logfileptr->OFS()<<perm<<endl;

    
}

void Ordering::NDBOX(){
    assert(pStructure!=NULL);

    if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
			throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call NDBOX\n" );
    }

    invp.Resize(pStructure->size);
    if(iam==0){
    Int k = std::ceil(std::pow(pStructure->size,1.0/3.0));
    logfileptr->OFS()<<"BOX K = "<<k<<endl;
    Int iflag =1;

    {
      vector<Int> stack(k*k>25?invp.m():2*invp.m());
      Int tmp = stack.size();
    FORTRAN(boxnd)( &k , &k, &k, &invp[0], &stack[0],&tmp, &iflag);
    }

    }

      // broadcast invp
      Int N = invp.m();
      MPI_Bcast(&invp[0],N*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
      perm.Resize(pStructure->size);
      for(Int i = 1; i <=N; ++i){
        Int node = invp[i-1];
        perm[node-1] = i;
      }


}

void Ordering::NDGRID(){
    assert(pStructure!=NULL);

    if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
			throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call NDGRID\n" );
    }

    invp.Resize(pStructure->size);
    if(iam==0){
    Int k = std::ceil(std::pow(pStructure->size,1.0/2.0));
    Int iflag =1;

    logfileptr->OFS()<<"GRID K = "<<k<<endl;

    {
      vector<Int> stack(k*k>25?invp.m():2*invp.m());
      Int tmp = stack.size();
    FORTRAN(gridnd)( &k , &k, &invp[0], &stack[0],&tmp, &iflag);
assert(iflag==0);
    }
    }

      // broadcast invp
      Int N = invp.m();
      MPI_Bcast(&invp[0],N*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
      perm.Resize(pStructure->size);
      for(Int i = 1; i <=N; ++i){
        Int node = invp[i-1];
        perm[node-1] = i;
      }


  
  //logfileptr->OFS()<<"perm "<<perm<<endl;
  //logfileptr->OFS()<<"invperm "<<invp<<endl;
}


void Ordering::AMD(){
    assert(pStructure!=NULL);

  logfileptr->OFS()<<"AMD used"<<endl;
  if(iam==0){cout<<"AMD used"<<endl;}

    if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
			throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call AMD\n" );
    }

    invp.Resize(pStructure->size);
    if(iam==0){
    Int iwsiz = 4*pStructure->size;
    std::vector<Int> iwork (iwsiz);
    perm.Resize(pStructure->size);

    Int nadj = pStructure->expRowind.m();
    Int nofsub =0;
    Int iflag =0;

    Int N = pStructure->size; 
    IntNumVec tmpXadj = pStructure->expColptr;
    //allocate extra N elements for AMD (elbow)
    IntNumVec tmpAdj(pStructure->expRowind.m()+N);
    std::copy(&pStructure->expRowind[0],&pStructure->expRowind[pStructure->expRowind.m()-1]+1,&tmpAdj[0]);

    Int IOVFLO = 2147483647;
    Int NCMPA;

    IntNumVec VTXDEG(N);
    IntNumVec QSIZE(N);
    IntNumVec ECFORW(N);
    IntNumVec MARKER(N);
    IntNumVec NVTXS(N+1);
    for(Int i=0;i<N;++i){
      NVTXS[i] = tmpXadj[i+1]-tmpXadj[i];
    }


    Int IWLEN = tmpXadj[N-1] + NVTXS[N-1] + N - 1 ;
    Int PFREE = tmpXadj[N-1] + NVTXS[N-1];


    FORTRAN(amdbar)( &N , &tmpXadj[0] , &tmpAdj[0] , &NVTXS[0] ,&IWLEN ,&PFREE ,  &QSIZE[0],&ECFORW[0] , &perm[0], &iwork[0], &invp[0],&VTXDEG[0], &NCMPA , &MARKER[0] , &IOVFLO );
    }

      // broadcast invp
      Int N = invp.m();
      MPI_Bcast(&invp[0],N*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
      perm.Resize(pStructure->size);
      for(Int i = 1; i <=N; ++i){
        Int node = invp[i-1];
        perm[node-1] = i;
      }


//  logfileptr->OFS()<<"perm "<<perm<<endl;
//  logfileptr->OFS()<<"invperm "<<invp<<endl;

//    Permute(perm);

}


#ifdef USE_METIS
  void Ordering::METIS(){

  logfileptr->OFS()<<"METIS used"<<endl;
  if(iam==0){cout<<"METIS used"<<endl;}

    assert(pStructure!=NULL);

    if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
      throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call AMD\n" );
    }


    invp.Resize(pStructure->size);
    if(iam==0){
      perm.Resize(pStructure->size);

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

    vector<int> tmpXadj(pStructure->expColptr.m());
    for(int i = 0; i<tmpXadj.size();++i){tmpXadj[i] = pStructure->expColptr[i];}
    vector<int> tmpAdj;
    tmpAdj.reserve(pStructure->expRowind.m());

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
    for(int col=0; col<invp.m();++col){ invp[col]++;}
    for(int col=0; col<perm.m();++col){ perm[col]++;}
    }


      // broadcast invp
      Int N = invp.m();
      MPI_Bcast(&invp[0],N*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
      perm.Resize(pStructure->size);
      for(Int i = 1; i <=N; ++i){
        Int node = invp[i-1];
        perm[node-1] = i;
      }




  }
#endif


#ifdef USE_PARMETIS
void Ordering::PARMETIS(){

  logfileptr->OFS()<<"PARMETIS used"<<endl;
  if(iam==0){cout<<"PARMETIS used"<<endl;}

  assert(pStructure!=NULL);

  if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
    throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call AMD\n" );
  }


  int baseval = 1;
  int N = pStructure->size; 
  invp.Resize(N);
  perm.Resize(N);
  int ndomains = (int)pow(2.0,std::floor(std::log2(np)));

  MPI_Comm ndcomm;
  Int mpirank;
  MPI_Comm_split(CommEnv_->MPI_GetComm(),iam<ndomains,mpirank,&ndcomm);
  if(iam<ndomains){
    vector<int> sizes(2*ndomains);
    vector<int> vtxdist(ndomains+1);

    int localN = N/ndomains;
    int fc = (iam)*localN+1;
    //build vtxdist vector
    for(int i = 0; i<ndomains;++i){ vtxdist[i] = i*localN+baseval; } vtxdist[ndomains] = N+baseval;
    if(iam==ndomains-1){
      localN = N - (ndomains-1)*localN;
    }
    //logfileptr->OFS()<<"vtxdist: "<<vtxdist<<endl;
    //logfileptr->OFS()<<"ndomains = "<<ndomains<<endl;

    //build local colptr and count nnz in local rowind
    vector<int> tmpXadj(localN+1);
    Ptr localNNZ = 0;
    int offset = pStructure->expColptr[fc-1]-1;
    tmpXadj[1-1] = 1;
    for(int col=1; col<=localN;++col){
      int column = col + fc -1;
      Ptr colbeg = pStructure->expColptr[column-1];
      Ptr colend = pStructure->expColptr[column]-1;
      //remove self from adjacency list
      tmpXadj[col] = tmpXadj[col-1] + colend - colbeg + 1 - 1;
      //logfileptr->OFS()<<col<<" <-> "<<column<<" | "<<colbeg<<"--"<<colend<<" <-> "<<tmpXadj[col-1]<<"--"<<tmpXadj[col]-1<<endl;
      localNNZ+=colend-colbeg+1;
    }

    //build local rowind, discarding self
    vector<int> tmpAdj;
    tmpAdj.reserve(localNNZ);
    for(int col=1; col<=localN;++col){
      int column = col + fc -1;
      Ptr colbeg = pStructure->expColptr[column-1];
      Ptr colend = pStructure->expColptr[column]-1;
      for(Ptr j=colbeg; j<=colend;++j){
        int row = pStructure->expRowind[j-1];
        if(row!=column){
          tmpAdj.push_back(row);
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
    vector<int> mpidispls(ndomains,0);
    vector<int> mpisizes(ndomains,0);
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
//  for(Int i = 1; i <= perm.m(); ++i){
//    Int node = perm[i-1];
//    invp[node-1] = i;
//  }
}
#endif

#ifdef USE_SCOTCH
void Ordering::SCOTCH(){

  logfileptr->OFS()<<"SCOTCH used"<<endl;
  if(iam==0){cout<<"SCOTCH used"<<endl;}

  assert(pStructure!=NULL);

  if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
    throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call SCOTCH\n" );
  }


  invp.Resize(pStructure->size);
  if(iam==0){
  perm.Resize(pStructure->size);

  vector<int64_t> tmpInvp(pStructure->size);
  vector<int64_t> tmpPerm(pStructure->size);

  int N = pStructure->size; 

  vector<int64_t> tmpXadj(pStructure->expColptr.m());
  for(int i = 0; i<tmpXadj.size();++i){tmpXadj[i] = pStructure->expColptr[i];}

  vector<int64_t> tmpAdj;
  tmpAdj.reserve(pStructure->expRowind.m());

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

  //switch everything to 0 based
  for(int col=0; col<tmpXadj.size();++col){ tmpXadj[col]--;}
  for(int col=0; col<tmpAdj.size();++col){ tmpAdj[col]--;}

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
    for(int col=0; col<invp.m();++col){ invp[col] = tmpInvp[col]+1;}
    for(int col=0; col<perm.m();++col){ perm[col] = tmpPerm[col]+1;}
    }

      // broadcast invp
      Int N = invp.m();
      MPI_Bcast(&invp[0],N*sizeof(int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
      perm.Resize(pStructure->size);
      for(Int i = 1; i <=N; ++i){
        Int node = invp[i-1];
        perm[node-1] = i;
      }


    }
#endif



#ifdef USE_PTSCOTCH
    void Ordering::PTSCOTCH(){

      logfileptr->OFS()<<"PTSCOTCH used"<<endl;
      if(iam==0){cout<<"PTSCOTCH used NOT WORKING"<<endl;}

      assert(pStructure!=NULL);

      if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
        throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call PTSCOTCH\n" );
      }


      SCOTCH_Num baseval = 0;
      SCOTCH_Num N = pStructure->size; 
      invp.Resize(N);
      perm.Resize(N);
      SCOTCH_Num ndomains = (SCOTCH_Num)pow(2.0,std::floor(std::log2(np)));

      MPI_Comm ndcomm;
      Int mpirank;
      MPI_Comm_split(CommEnv_->MPI_GetComm(),iam<ndomains,mpirank,&ndcomm);
      if(iam<ndomains){
        vector<SCOTCH_Num> vtxdist(ndomains+1);

        SCOTCH_Num localN = N/ndomains;
        SCOTCH_Num fc = (iam)*localN+1;
        //build vtxdist vector
        for(SCOTCH_Num i = 0; i<ndomains;++i){ vtxdist[i] = i*localN+baseval; } vtxdist[ndomains] = N+baseval;
        if(iam==ndomains-1){
          localN = N - (ndomains-1)*localN;
        }
        //logfileptr->OFS()<<"vtxdist: "<<vtxdist<<endl;
        //logfileptr->OFS()<<"ndomains = "<<ndomains<<endl;

        //build local colptr and count nnz in local rowind
        vector<SCOTCH_Num> tmpXadj(localN+1);
        Ptr localNNZ = 0;
        SCOTCH_Num offset = pStructure->expColptr[fc-1]-1;
        tmpXadj[1-1] = 1;
        for(SCOTCH_Num col=1; col<=localN;++col){
          SCOTCH_Num column = col + fc -1;
          Ptr colbeg = pStructure->expColptr[column-1];
          Ptr colend = pStructure->expColptr[column]-1;
          //remove self from adjacency list
          tmpXadj[col] = tmpXadj[col-1] + colend - colbeg + 1 - 1;
          //logfileptr->OFS()<<col<<" <-> "<<column<<" | "<<colbeg<<"--"<<colend<<" <-> "<<tmpXadj[col-1]<<"--"<<tmpXadj[col]-1<<endl;
          localNNZ+=colend-colbeg+1;
        }

        //build local rowind, discarding self
        vector<SCOTCH_Num> tmpAdj;
        tmpAdj.reserve(localNNZ);
        for(SCOTCH_Num col=1; col<=localN;++col){
          SCOTCH_Num column = col + fc -1;
          Ptr colbeg = pStructure->expColptr[column-1];
          Ptr colend = pStructure->expColptr[column]-1;
          for(Ptr j=colbeg; j<=colend;++j){
            SCOTCH_Num row = pStructure->expRowind[j-1];
            if(row!=column){
              tmpAdj.push_back(row);
            }
          }
        }


        //switch everything to 0 based
        //logfileptr->OFS()<<tmpXadj<<endl;
        for(SCOTCH_Num col=0; col<localN+1;++col){ tmpXadj[col]+=(baseval-1);}
        //logfileptr->OFS()<<tmpXadj<<endl;
        //logfileptr->OFS()<<tmpAdj<<endl;
        for(SCOTCH_Num col=0; col<localNNZ;++col){ tmpAdj[col]+=(baseval-1);}
        //logfileptr->OFS()<<tmpAdj<<endl;

        SCOTCH_Num options[3];
        options[0] = 0;
        SCOTCH_Num numflag = baseval;


        int npnd;
        MPI_Comm_size (ndcomm, &npnd);
        //logfileptr->OFS()<<"PROC ND: "<<npnd<<endl;



        SCOTCH_Dgraph       grafdat;                    /* Scotch distributed graph object to SCOTCH_Numerface with libScotch    */
        SCOTCH_Dordering    ordedat;                    /* Scotch distributed ordering object to SCOTCH_Numerface with libScotch */
        SCOTCH_Strat        stradat;
        SCOTCH_Num          vertlocnbr;
        SCOTCH_Num          edgelocnbr;

        assert(SCOTCH_dgraphInit (&grafdat, ndcomm) == 0);

        vertlocnbr = vtxdist[iam + 1] - vtxdist[iam];
        edgelocnbr = tmpXadj[vertlocnbr] - baseval;

        if (SCOTCH_dgraphBuild (&grafdat, baseval,
              vertlocnbr, vertlocnbr, &tmpXadj[0], NULL, NULL, NULL,
              edgelocnbr, edgelocnbr, &tmpAdj[0], NULL, NULL) == 0) {
          SCOTCH_stratInit (&stradat);
          {
            if (SCOTCH_dgraphOrderInit (&grafdat, &ordedat) == 0) {
              SCOTCH_dgraphOrderCompute (&grafdat, &ordedat, &stradat);
              vector<SCOTCH_Num> sc_perm(perm.m());
              for(SCOTCH_Num i=0;i<sc_perm.size();i++){ sc_perm[i] = (SCOTCH_Num)perm[i]; }
              SCOTCH_dgraphOrderPerm(&grafdat, &ordedat, &sc_perm[0]);
              for(SCOTCH_Num i=0;i<sc_perm.size();i++){ perm[i] = sc_perm[i]; }
              SCOTCH_dgraphOrderExit (&grafdat, &ordedat);
            }
          }
          SCOTCH_stratExit (&stradat);
        }
        SCOTCH_dgraphExit (&grafdat);


        //logfileptr->OFS()<<"Order: "<<endl;
        //for(SCOTCH_Num i =0;i<localN;++i){
        //  //logfileptr->OFS()<<invp[i]<<" ";
        //  logfileptr->OFS()<<perm[i]<<" ";
        //}
        //logfileptr->OFS()<<endl;


        //compute displs
        vector<int> mpidispls(ndomains,0);
        vector<int> mpisizes(ndomains,0);
        for(int p = 1;p<=ndomains;++p){
          mpisizes[p-1] = (vtxdist[p] - vtxdist[p-1])*sizeof(Int);
          mpidispls[p-1] = (vtxdist[p-1]-baseval)*sizeof(Int);
        }

        //gather on the root
        MPI_Gatherv(&perm[0],mpisizes[iam],MPI_BYTE,&invp[0],&mpisizes[0],&mpidispls[0],MPI_BYTE,0,ndcomm);
        //    MPI_Gatherv(&invp[0],mpisizes[iam],MPI_BYTE,&perm[0],&mpisizes[0],&mpidispls[0],MPI_BYTE,0,comm);

        if(iam==0){
          //switch everything to 1 based
          for(int col=0; col<N;++col){ invp[col]+=(1-baseval);}

          //logfileptr->OFS()<<"Full Order: "<<endl;
          //for(SCOTCH_Num i =0;i<N;++i){
          //  logfileptr->OFS()<<invp[i]<<" ";
          //  //logfileptr->OFS()<<perm[i]<<" ";
          //}
          //logfileptr->OFS()<<endl;
        }


      }

      MPI_Comm_free(&ndcomm);

      // broadcast invp
      MPI_Bcast(&invp[0],N*sizeof(Int),MPI_BYTE,0,CommEnv_->MPI_GetComm());
      for(Int i = 1; i <=N; ++i){
        Int node = invp[i-1];
        perm[node-1] = i;
      }

      //  // broadcast perm
      //  MPI_Bcast(&perm[0],N*sizeof(SCOTCH_Num),MPI_BYTE,0,comm);
      //
      //  for(Int i = 1; i <= perm.m(); ++i){
      //    Int node = perm[i-1];
      //    invp[node-1] = i;
      //  }
    }
#endif


    void Ordering::SetCommEnvironment(CommEnvironment * CommEnv){
      CommEnv_=CommEnv;
    }

void Ordering::Compose(IntNumVec & invp2){
    assert(pStructure!=NULL);
      //Compose the two permutations
      for(Int i = 1; i <= invp.m(); ++i){
            Int interm = invp[i-1];
            invp[i-1] = invp2[interm-1];
      }
      for(Int i = 1; i <= invp.m(); ++i){
        Int node = invp[i-1];
        perm[node-1] = i;
      }
}

}
