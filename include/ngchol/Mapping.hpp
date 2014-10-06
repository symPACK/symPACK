#ifndef _MAPPING_HPP_
#define _MAPPING_HPP_

#include "ngchol/Environment.hpp"


#ifndef MAPCLASS
#define MAPCLASS Modwrap2D
#endif

namespace LIBCHOLESKY{


class Mapping{
  protected:
    std::vector<Int> * pProcMap_;
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
      void Update(std::vector<Int> & aProcMap){
        if(pProcMap_==NULL){
          pProcMap_ = new std::vector<Int>(aProcMap);
        }
        else{
          pProcMap_->clear();
          pProcMap_->insert(aProcMap.begin(),aProcMap.end(),pProcMap_->begin());
        }
      }
  
      Mapping(Int aiNumProc, Int aiPRows, Int aiPCols, Int aiBlockSize = 1):iNumProc_(aiNumProc),iPRows_(aiPRows),iPCols_(aiPCols),iBlockSize_(aiBlockSize),pProcMap_(NULL){};
      Mapping(Int aiNumProc, Int aiPRows, Int aiPCols, std::vector<Int> & aProcMap, Int aiBlockSize = 1):
          iNumProc_(aiNumProc),iPRows_(aiPRows),iPCols_(aiPCols),iBlockSize_(aiBlockSize){
          pProcMap_ = new std::vector<Int>(aProcMap);
      };

      Mapping(Mapping & C){
        iNumProc_   = C.iNumProc_;
        iPRows_     = C.iPRows_;
        iPCols_     = C.iPCols_;
        iBlockSize_ = C.iBlockSize_;
        if(C.pProcMap_!=NULL){
          pProcMap_ = new std::vector<Int>(*C.pProcMap_);
        }
      }
      Mapping(){
        pProcMap_=NULL;
      };
      virtual ~Mapping(){
        if(pProcMap_!=NULL){
          delete pProcMap_;
        }
      };
      virtual inline Int Map(Int i, Int j)=0;
};



class Modwrap2D: public Mapping{
  protected:
    inline Int modwrap2D_(Int i, Int j) {
      Int p = std::min(i/iBlockSize_,j/iBlockSize_)%iPRows_ +
               iPRows_*floor((double)(std::max(i/iBlockSize_,j/iBlockSize_)%iNumProc_)/(double)iPRows_);
      return p;
    }
  public:
      Modwrap2D(Int aiNumProc, Int aiPRows, Int aiPCols, Int aiBlockSize = 1):Mapping(aiNumProc,aiPRows,aiPCols,aiBlockSize){};
      Modwrap2D(Int aiNumProc, Int aiPRows, Int aiPCols, std::vector<Int> & aProcMap, Int aiBlockSize = 1):Mapping(aiNumProc,aiPRows,aiPCols,aProcMap,aiBlockSize){};
//      Modwrap2D(Modwrap2D & C):Mapping(C.iNumProc_,C.iPRows_,C.iPCols_,C.iBlockSize_){};
      Modwrap2D(Modwrap2D & C):Mapping(C){};
      Modwrap2D():Mapping(0,0,0,0){};
      inline Int Map(Int i, Int j){ 
        Int proc =modwrap2D_(i,j);
        if(pProcMap_!=NULL){
          proc = modwrap2D_((*pProcMap_)[i],(*pProcMap_)[j]);
        }
        return proc;
      }
      ~Modwrap2D(){};
};


class Modwrap2DNS: public Mapping{
  protected:
    inline Int modwrap2DNS_(Int i, Int j) { Int p = i/iBlockSize_%iPRows_ + iPRows_*floor((double)((j/iBlockSize_)%iNumProc_)/(double)iPRows_); return p;}
  public:
      Modwrap2DNS(Int aiNumProc, Int aiPRows, Int aiPCols, Int aiBlockSize = 1):Mapping(aiNumProc,aiPRows,aiPCols,aiBlockSize){};
//      Modwrap2DNS(Modwrap2DNS & C):Mapping(C.iNumProc_,C.iPRows_,C.iPCols_,C.iBlockSize_){};
      Modwrap2DNS(Modwrap2DNS & C):Mapping(C){};
      Modwrap2DNS():Mapping(0,0,0,0){};
      inline Int Map(Int i, Int j){ return modwrap2DNS_(i,j);}
      ~Modwrap2DNS(){};
};


class Row2D: public Mapping{
  protected:
    inline Int row2D_(Int i, Int j) { Int p = (i/iBlockSize_)%iNumProc_; return p;}
  public:
      Row2D(Int aiNumProc, Int aiPRows, Int aiPCols, Int aiBlockSize = 1):Mapping(aiNumProc,aiNumProc,aiPCols,aiBlockSize){};
      //Row2D(Row2D & C):Mapping(C.iNumProc_,C.iPRows_,C.iPCols_,C.iBlockSize_){};
      Row2D(Row2D & C):Mapping(C){};
      Row2D():Mapping(0,0,0,0){};
      inline Int Map(Int i, Int j){ return row2D_(i,j);}
      ~Row2D(){};
};

class Col2D: public Mapping{
  protected:
    inline Int col2D_(Int i, Int j) { Int p = (j/iBlockSize_)%iNumProc_; return p;}
  public:
      Col2D(Int aiNumProc, Int aiPRows, Int aiPCols, Int aiBlockSize = 1):Mapping(aiNumProc,aiNumProc,aiPCols,aiBlockSize){};
      //Col2D(Row2D & C):Mapping(C){}
      Col2D(Col2D & C):Mapping(C){};
      Col2D():Mapping(0,0,0,0){};
      inline Int Map(Int i, Int j){ return col2D_(i,j);}
      ~Col2D(){};
};


}

#endif //_MAPPING_HPP_
