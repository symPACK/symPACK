#include "ngchol/Ordering.hpp"
#include "ngchol/mmd.hpp"

namespace LIBCHOLESKY{


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

    if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
			throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call MMD\n" );
    }

    Int iwsiz = 4*pStructure->size;
    std::vector<Int> iwork (iwsiz);
    invp.Resize(pStructure->size);
    perm.Resize(pStructure->size);

    Int nadj = pStructure->expRowind.m();
    Int nofsub =0;
    Int iflag =0;

    IntNumVec tmpXadj = pStructure->expColptr;
    IntNumVec tmpAdj = pStructure->expRowind;
    FORTRAN(ordmmd)( &pStructure->size , &nadj , tmpXadj.Data() , tmpAdj.Data(), 
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

    Int iwsiz = 4*pStructure->size;
    std::vector<Int> iwork (iwsiz);
    invp.Resize(pStructure->size);
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

//  logfileptr->OFS()<<"perm "<<perm<<endl;
//  logfileptr->OFS()<<"invperm "<<invp<<endl;

//    Permute(perm);

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
