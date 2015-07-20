#include "ngchol/Ordering.hpp"
#include <cassert>
#include <stdexcept>

#ifndef Add_
#define FORTRAN(name) name
#else
#define FORTRAN(name) name##_
#endif


namespace LIBCHOLESKY {

extern "C" {

void FORTRAN(ordmmd)( int * neqns , int * nadj  , int * xadj  ,
         int * adjncy, int * invp  , int * perm  , int * iwsiz ,
                         int * iwork , int * nofsub, int * iflag);


void FORTRAN(amdbar) (int * N, int * PE, int * IW, int * LEN, int * IWLEN, int * PFREE, int * NV, int * NEXT, int *
         LAST, int * HEAD, int * ELEN, int * DEGREE, int * NCMPA, int * W, int * IOVFLO);


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
    vector<int> tmpAdj = pStructure->expRowind;
    FORTRAN(ordmmd)( &pStructure->size , &nadj , &tmpXadj[0] , &tmpAdj[0], 
            &invp[0] , &perm[0] , &iwsiz , &iwork[0] , &nofsub, &iflag ) ;


//  logfileptr->OFS()<<"perm "<<perm<<endl;
//  logfileptr->OFS()<<"invperm "<<invp<<endl;

//    Permute(perm);

    assert(iflag == 0);
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
