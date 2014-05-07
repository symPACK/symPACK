#ifndef _MAPPING_HPP_
#define _MAPPING_HPP_

#include "Environment.hpp"


#ifndef MAPCLASS
#define MAPCLASS Modwrap2D
#endif

namespace LIBCHOLESKY{


class Mapping{
  protected:
    Int iNumProc_;
    Int iPRows_;
    Int iPCols_;
    Int iBlockSize_;
//      inline Int row2D(Int i, Int j) {return (i/iBlockSize_)%iNumProc_;}
//      inline Int col2D(Int i, Int j) {return (j/iBlockSize_)%iNumProc_;}
//
//      inline Int chevron2D(Int i, Int j) {return (min(i,j)/iBlockSize_)%iNumProc_;}
//      inline Int antichevron2D(Int i, Int j) {return (max(i,j)/iBlockSize_)%iNumProc_;}
//
//      inline Int diag2D(Int i, Int j) {return (abs(i-j)/iBlockSize_)%iNumProc_;}
//      inline Int antidiag2D(Int i, Int j) {return ((i+j)/iBlockSize_)%iNumProc_;}
//
//
//      inline Int modwrap2D(Int i, Int j) {return min(i/iBlockSize_,j/iBlockSize_)%iPRows_ + iPRows_*floor((double)(max(i/iBlockSize_,j/iBlockSize_)%iNumProc_)/(double)iPRows_);}
//      inline Int modwrap2Dns(Int i, Int j) {return(i/iBlockSize_)%iPRows_ + iPRows_*floor((double)((j/iBlockSize_)%iNumProc_)/(double)iPRows_);}

  public:

      Mapping(Int aiNumProc, Int aiPRows, Int aiPCols, Int aiBlockSize = 1):iNumProc_(aiNumProc),iPRows_(aiPRows),iPCols_(aiPCols),iBlockSize_(aiBlockSize){};
      Mapping(Mapping & C):Mapping(C.iNumProc_,C.iPRows_,C.iPCols_,C.iBlockSize_){};
      Mapping():Mapping(0,0,0,0){};

      Int Map(Int i, Int j){ abort(); return 0;}
};



class Modwrap2D: public Mapping{
  protected:
    inline Int modwrap2D_(Int i, Int j) {return std::min(i/iBlockSize_,j/iBlockSize_)%iPRows_ + iPRows_*floor((double)(std::max(i/iBlockSize_,j/iBlockSize_)%iNumProc_)/(double)iPRows_);}
  public:
      Modwrap2D(Int aiNumProc, Int aiPRows, Int aiPCols, Int aiBlockSize = 1):Mapping(aiNumProc,aiPRows,aiPCols,aiBlockSize){};
      Modwrap2D(Modwrap2D & C):Mapping(C.iNumProc_,C.iPRows_,C.iPCols_,C.iBlockSize_){};
      Modwrap2D():Mapping(0,0,0,0){};
      Int Map(Int i, Int j){ return modwrap2D_(i,j);}
};

}

#endif //_MAPPING_HPP_
