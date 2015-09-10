#include "ngchol/Ordering.hpp"
#include "ngchol/Utility.hpp"
#include <cassert>
#include <stdexcept>
#include <iostream>
using namespace std;

#ifndef Add_
#define FORTRAN(name) name
#else
#define FORTRAN(name) name##_
#endif

#if 0
#include <metis.h>
#endif

namespace LIBCHOLESKY {

extern "C" {

void FORTRAN(ordmmd)( int * neqns , int * nadj  , int * xadj  ,
         int * adjncy, int * invp  , int * perm  , int * iwsiz ,
                         int * iwork , int * nofsub, int * iflag);


void FORTRAN(amdbar) (int * N, int * PE, int * IW, int * LEN, int * IWLEN, int * PFREE, int * NV, int * NEXT, int *
         LAST, int * HEAD, int * ELEN, int * DEGREE, int * NCMPA, int * W, int * IOVFLO);

int metis_nodend (int * N     , int* XADJ2 , int* ADJ2  , int * VWGT, int* OPTION, int* dback , int* dforw);
 }


}

namespace LIBCHOLESKY{


void Ordering::SetStructure(SparseMatrixStructure & aGlobal){
    assert(pStructure==NULL);
    pStructure = &aGlobal;
    invp.resize(pStructure->size);
    perm.resize(pStructure->size);
    for(int i =0;i<invp.size();++i){
      perm[i]=i+1;
      invp[i]=i+1;
    }
}


void Ordering::MMD(){
    assert(pStructure!=NULL);

    if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
			throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call MMD\n" );
    }

    int iwsiz = 4*pStructure->size;
    std::vector<int> iwork (iwsiz);
    invp.resize(pStructure->size);
    perm.resize(pStructure->size);

    int nadj = pStructure->expRowind.size();
    int nofsub =0;
    int iflag =0;

    vector<int> tmpXadj = pStructure->expColptr;
    vector<int> tmpAdj;
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


    FORTRAN(ordmmd)( &pStructure->size , &nadj , &tmpXadj[0] , &tmpAdj[0], 
            &invp[0] , &perm[0] , &iwsiz , &iwork[0] , &nofsub, &iflag ) ;


//    Permute(perm);

    assert(iflag == 0);
}

void Ordering::METIS(){
  assert(pStructure!=NULL);

  if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
    throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call AMD\n" );
  }


  invp.resize(pStructure->size);
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

    vector<int> tmpXadj = pStructure->expColptr;
    vector<int> tmpAdj;
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

#ifdef USE_METIS
  metis_nodend( &N, &tmpXadj[0] , &tmpAdj[0]  , NULL , NULL, &perm[0] , &invp[0] );
#else
  int numflag = 0;
  metis_nodend( &N, &tmpXadj[0] , &tmpAdj[0]  , &numflag , NULL, &perm[0] , &invp[0] );
#endif

    //switch everything to 1 based
    for(int col=0; col<invp.size();++col){ invp[col]++;}
    for(int col=0; col<perm.size();++col){ perm[col]++;}
//    for(int col=0; col<tmpXadj.size();++col){ tmpXadj[col]++;}
//    for(int col=0; col<tmpAdj.size();++col){ tmpAdj[col]++;}
}

void Ordering::AMD(){
    assert(pStructure!=NULL);

    if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
			throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call AMD\n" );
    }

    int iwsiz = 4*pStructure->size;
    std::vector<int> iwork (iwsiz);
    invp.resize(pStructure->size);
    perm.resize(pStructure->size);


    int nadj = pStructure->expRowind.size();
    int nofsub =0;
    int iflag =0;

    int N = pStructure->size; 
    vector<int> tmpXadj = pStructure->expColptr;
    //allocate extra N elements for AMD (elbow)
    vector<int> tmpAdj(pStructure->expRowind.size()+N);
    std::copy(&pStructure->expRowind[0],&pStructure->expRowind[pStructure->expRowind.size()-1]+1,&tmpAdj[0]);

    int IOVFLO = 2147483647;
    int NCMPA;

    vector<int> VTXDEG(N);
    vector<int> QSIZE(N);
    vector<int> ECFORW(N);
    vector<int> MARKER(N);
    vector<int> NVTXS(N+1);
    for(int i=0;i<N;++i){
      NVTXS[i] = tmpXadj[i+1]-tmpXadj[i];
    }


    int IWLEN = tmpXadj[N-1] + NVTXS[N-1] + N - 1 ;
    int PFREE = tmpXadj[N-1] + NVTXS[N-1];


    FORTRAN(amdbar)( &N , &tmpXadj[0] , &tmpAdj[0] , &NVTXS[0] ,&IWLEN ,&PFREE ,  &QSIZE[0],&ECFORW[0] , &perm[0], &iwork[0], &invp[0],&VTXDEG[0], &NCMPA , &MARKER[0] , &IOVFLO );

//  logfileptr->OFS()<<"perm "<<perm<<endl;
//  logfileptr->OFS()<<"invperm "<<invp<<endl;

//    Permute(perm);

}

void Ordering::Compose(vector<int> & invp2){
    assert(pStructure!=NULL);

//    cout<<"composing: ";
//    for(int i=0;i<invp2.size();++i){
//      cout<<invp2[i]<<" ";
//    }
//    cout<<endl;


      //Compose the two permutations
      for(int i = 1; i <= invp.size(); ++i){
            int interm = invp[i-1];
            invp[i-1] = invp2[interm-1];
      }
      for(int i = 1; i <= invp.size(); ++i){
        int node = invp[i-1];
        perm[node-1] = i;
      }
}

}
