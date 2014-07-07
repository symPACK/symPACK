#include "ngchol/Cholesky.hpp"

#include "ngchol/NumMat.hpp"
#include <vector>

namespace LIBCHOLESKY{


void cmod(DistSparseMatrix<Scalar> & Amat, IntNumVec & xsuper,Int J, Int K, NumMat<Scalar> & W){
    //determine L_JK
    Int firstRow = 0;
    Int lastRow = 0;
    Int firstRowIdx = 0;
    Int lastRowIdx = 0;
    for(Int rowIdx=Amat.Global.colptr(K.lastCol-1); rowIdx<Amat.Global.colptr(K.lastCol); rowIdx++){
      Int row = Amat.Global.rowind(rowIdx-1);
      if(row< xsuper(J-1)) continue;
      if(row> xsuper(J)-1) break;
      if(firstRow==0) { firstRowIdx=rowIdx; firstRow = row; }
      lastRow = row;
      lastRowIdx = rowIdx;
    }

    logfileptr->OFS()<<"L_JK would be L_"<<firstRow<<".."<<lastRow<<","<<xsuper(K-1)<<".."<<xsuper(K)-1<<std::endl;
 
     

}


void LL_SS_Fact(DistSparseMatrix<Scalar> & Amat, IntNumVec & xsuper, IntNumVec & cc, IntNumVec & rc ){
  NumMat<Scalar> workspace;
  //What is the max workspace required ? UPPER BOUND : max cc
  Int size = 0;
  for(Int i = 0; i<cc.m();++i){ size = std::max(size, cc(i) );}
  workspace.Resize(size);
  SetValue(workspace, SCALAR_ZERO);


  IntNumVec S(xsuper.());
  //zero means empty set or end of set
  SetValue(S,I_ZERO);

  for(Int J=1; J<=xsuper.m(); J++){
    Int K = S(J-1);
    Int Kidx = J-1;
    while(S(Kidx)!=0){
      //Determine J', the subset of columns in J updated by columns in K
      
      cmod(Amat, xsuper, J, K, workspace);
      
      Kidx = S(Kidx); 
    }

  }
  
}

}
