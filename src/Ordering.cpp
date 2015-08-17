#include "ngchol/Ordering.hpp"
#include "ngchol/mmd.hpp"
#include "ngchol/utility.hpp"

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








    FORTRAN(ordmmd)( &pStructure->size , &nadj , tmpXadj.Data() , tmpAdj.Data(), 
            &invp[0] , &perm[0] , &iwsiz , &iwork[0] , &nofsub, &iflag ) ;


  logfileptr->OFS()<<"perm "<<perm<<endl;
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



void Ordering::METIS(){
/////  assert(pStructure!=NULL);
/////
/////  if(!pStructure->bIsGlobal || !pStructure->bIsExpanded){
/////    throw std::logic_error( "SparseMatrixpStructure->must be global and expanded in order to call AMD\n" );
/////  }
/////
/////
/////  invp.resize(pStructure->size);
/////  perm.resize(pStructure->size);
/////  int options[METIS_NOPTIONS];
/////  METIS_SetDefaultOptions(&options[0]);
/////options[METIS_OPTION_NUMBERING] = 1;
///////options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM;
/////////options[METIS_OPTION_OFLAGS] = 0;
///////options[METIS_OPTION_PFACTOR] = 0;
///////options[METIS_OPTION_NSEPS] = 0;
/////options[METIS_OPTION_DBGLVL] = 8;
/////
/////
/////
///////    Opt [0] = 0 ; /* use defaults */
/////////    Opt [1] = 3 ; /* matching type */
///////options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
///////    Opt [2] = 1 ; /* init. partitioning algo*/
///////    Opt [3] = 2 ; /* refinement algorithm */
///////    Opt [4] = 0 ; /* no debug */
///////    Opt [5] = 1 ; /* initial compression */
///////    Opt [6] = 0 ; /* no dense node removal */
///////    Opt [7] = 1 ; /* number of separators @ each step */
/////
/////
/////
/////
/////
/////
/////
/////
///////METIS_OPTION_CTYPE, METIS_OPTION_RTYPE, METIS_OPTION_NO2HOP,
///////METIS_OPTION_NSEPS, METIS_OPTION_NITER, METIS_OPTION_UFACTOR,
///////METIS_OPTION_COMPRESS, METIS_OPTION_CCORDER, METIS_OPTION_SEED,
///////METIS_OPTION_PFACTOR, METIS_OPTION_NUMBERING, METIS_OPTION_DBGLVL
/////
/////
/////
///////options[METIS_OPTION_CCORDER] = 1;
///////options[METIS_OPTION_COMPRESS] = 1;
///////  OPTION[0] = 0;
///////    OPTION[1] = 3;
///////    OPTION[2] = 1;
///////    OPTION[3] = 2;
///////    OPTION[4] = 0;
///////    OPTION[5] = 1;
///////    OPTION[6] = 0;
///////    OPTION[7] = 1;
/////    int N = pStructure->size; 
/////
/////    vector<int> tmpXadj = pStructure->expColptr;
/////    vector<int> tmpAdj;
/////    tmpAdj.reserve(pStructure->expRowind.size());
/////
/////    for(int col=0; col<tmpXadj.size()-1;++col){
/////      for(int j=tmpXadj[col]; j<tmpXadj[col+1];++j){
/////        if( pStructure->expRowind[j-1]-1 != col){
/////          tmpAdj.push_back(pStructure->expRowind[j-1]);
/////        }
/////      }
/////    }
/////
/////    int rm = 0;
/////    for(int col=0; col<tmpXadj.size();++col){
/////      tmpXadj[col]-=rm;
/////      rm++;
/////    }
/////
/////
///////    vector<int> tmpXadj = pStructure->expColptr;
///////    vector<int> tmpAdj = pStructure->expRowind;
/////
/////  //convert to 0 based
///////  for(int i =0; i<tmpXadj.size();++i){tmpXadj[i]--;}
///////  for(int i =0; i<tmpAdj.size();++i){tmpAdj[i]--;}
/////  METIS_NodeND( &N, &tmpXadj[0] , &tmpAdj[0]  , NULL, &options[0], &perm[0] , &invp[0] );
///////  METIS_NodeND( &N, &tmpXadj[0] , &tmpAdj[0]  , NULL, &options[0], &invp[0] , &perm[0] );
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
