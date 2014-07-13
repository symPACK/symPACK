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
