/// @file FBMatrix.cpp
/// @brief Matrix structure and methods for the FAN-BOTH method.
/// @author Mathias Jacquelin
/// @date 2014-01-21

#include "ngchol/FBMatrix.hpp"

#include  "ngchol/utility.hpp"
#include  "ngchol/blas.hpp"
#include  "ngchol/lapack.hpp"

#include  "ngchol/LogFile.hpp"
#include  "ngchol/upcxx_additions.hpp"




namespace LIBCHOLESKY{



  FBMatrix::FBMatrix():blksize(1),n(0){
  }

  FBMatrix::~FBMatrix(){
    AchunkLower.clear();
    WLower.clear();
  }

  void FBMatrix::ClearTmp(){
    WLower.clear();
  }


  void FBMatrix::Allocate(Int & np, Int pn, Int pblksize){
    n=pn;
    blksize=pblksize;


    //for 2D maps
    prow = (Int)sqrt((double)np);


#ifdef MAP2D  
    np = (Int)(np/prow) * prow;
#endif

    this->np = np;
    //determine pcol and prow
    Int firstproc = MAP(0,0);
    Int j=0;
    Int curproc = -1;
    pcol =0;
    do{
      j+=blksize;
      if(j<pn){
        curproc=MAP(j,j);
#ifdef _DEBUG_
        logfileptr->OFS()<<"firstproc = "<<firstproc<<" vs "<<curproc<<endl;
#endif
        ++pcol;
      }
    }while(firstproc!=curproc && j<pn);

    assert(pcol <= np);

    Int totBlk = n/blksize;
    Int remaining = n-totBlk*blksize;
    if(remaining>0) totBlk++;

    Int localBlk = totBlk/pcol;
    Int additionalBlock = totBlk - pcol*localBlk;
    if(iam<additionalBlock){ localBlk++; }


    Int prevBlk = (iam) +pcol*(localBlk-1); 
    Int chksize = (localBlk-1)*blksize + min(n-(prevBlk)*blksize,blksize);

    Int lowerChkSize = (localBlk-1)*(n-iam*blksize)*blksize - blksize*blksize*pcol*(localBlk-1)*localBlk/2 + (n-(prevBlk)*blksize)*min(n-(prevBlk)*blksize,blksize); 

    Int numStep = ceil((double)n/(double)blksize);
    WLower.resize(numStep,NULL); //, DblNumMat(0,0) );
    AchunkLower.resize(localBlk,NULL);//, DblNumMat(0,0) );

    //WARNING this function needs to be completed by specialized classes
  } 

  void FBMatrix::Aggregate(Int j, DblNumMat &DistW){
    TIMER_START(Aggregate);
    //j is local
    Int local_j = this->global_col_to_local(j);
    DblNumMat & LocalChunk = *AchunkLower[local_j/blksize];

    Int jb = min(blksize, n-j);
    //aggregate previous updates

    //assert(DistW.Size()==LocalChunk.Size());
    blas::Axpy(LocalChunk.Size(), -1.0, &DistW(0,0),1,&LocalChunk(0,0),1);

    TIMER_STOP(Aggregate);
  }

  void FBMatrix::Factor(Int j){
    TIMER_START(Factor);
    Int jb = min(blksize, n-j);
    //j is local
    Int local_j = global_col_to_local(j);
    DblNumMat * LocalChunk = AchunkLower[local_j/blksize];

    lapack::Potrf( 'L', jb, &LocalChunk->at(0,0 ), n-j);
    if(n-j-jb>0){
      blas::Trsm('R','L','T','N',n-j-jb,jb, 1.0, &LocalChunk->at(0,0), n-j, &LocalChunk->at(jb,0), n-j);
    }
    TIMER_STOP(Factor);
  }

  void FBMatrix::Factor_ref(Int j){
    TIMER_START(Factor);
    Int jb = min(blksize, n-j);
    //j is local
    Int local_j = global_col_to_local(j);
    DblNumMat & LocalChunk = *AchunkLower[local_j/blksize];

    lapack::Potrf( 'L', jb, &LocalChunk.at(0,0 ), n-j);
    if(n-j-jb>0){
      blas::Trsm('R','L','T','N',n-j-jb,jb, 1.0, &LocalChunk.at(0,0), n-j, &LocalChunk.at(jb,0), n-j);
    }
    TIMER_STOP(Factor);
  }

  void FBMatrix::Update(Int j, Int i, DblNumMat & Factor){
    TIMER_START(Update);
    DblNumMat & WChunk = *WLower[i/blksize];
    //i is local, j is global
    Int jb = min(blksize, n-j);
    Int jbi = min(blksize, n-i);

    Int mf = Factor.m();
    Int nf = Factor.n();

    //call dgemm
    blas::Gemm('N','T',n-i,jbi,jb,1.0,&Factor(i-j,0),n-j,&Factor(i-j,0),n-j,1.0,&WChunk(0,0),n-i);
    TIMER_STOP(Update);
  }


  void FBMatrix::WaitFactorization(){     
  }

}


