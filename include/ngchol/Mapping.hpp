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

      void Dump(Int n){
        logfileptr->OFS()<<"Resulting mapping: "<<endl;
        logfileptr->OFS()<<"np = "<<iNumProc_<<endl;
        logfileptr->OFS()<<"prows = "<<iPRows_<<endl;
        logfileptr->OFS()<<"pcols = "<<iPCols_<<endl;
        for(Int i =0;i<n;++i){
          for(Int j = 0;j<=i;++j){
            logfileptr->OFS()<<Map(i,j)<<" ";
          }
          logfileptr->OFS()<<endl;
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
        Int proc = 0;
        if(pProcMap_!=NULL){
          proc = modwrap2D_((*pProcMap_)[i],(*pProcMap_)[j]);
        }
        else{
          proc = modwrap2D_(i,j);
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
      Modwrap2DNS(Int aiNumProc, Int aiPRows, Int aiPCols, std::vector<Int> & aProcMap, Int aiBlockSize = 1):Mapping(aiNumProc,aiPRows,aiPCols,aProcMap,aiBlockSize){};
//      Modwrap2DNS(Modwrap2DNS & C):Mapping(C.iNumProc_,C.iPRows_,C.iPCols_,C.iBlockSize_){};
      Modwrap2DNS(Modwrap2DNS & C):Mapping(C){};
      Modwrap2DNS():Mapping(0,0,0,0){};
      inline Int Map(Int i, Int j){
       Int proc = 0;
        if(pProcMap_!=NULL){
          proc = modwrap2DNS_((*pProcMap_)[i],(*pProcMap_)[j]);
        }
        else{
          proc = modwrap2DNS_(i,j);
        }
        return proc;
      }
      ~Modwrap2DNS(){};
};


class Row2D: public Mapping{
  protected:
    inline Int row2D_(Int i, Int j) { Int p = (i/iBlockSize_)%iNumProc_; return p;}
  public:
      Row2D(Int aiNumProc, Int aiPRows, Int aiPCols, Int aiBlockSize = 1):Mapping(aiNumProc,aiPRows,aiPCols,aiBlockSize){};
      Row2D(Int aiNumProc, Int aiPRows, Int aiPCols, std::vector<Int> & aProcMap, Int aiBlockSize = 1):Mapping(aiNumProc,aiPRows,aiPCols,aProcMap,aiBlockSize){};
      //Row2D(Row2D & C):Mapping(C.iNumProc_,C.iPRows_,C.iPCols_,C.iBlockSize_){};
      Row2D(Row2D & C):Mapping(C){};
      Row2D():Mapping(0,0,0,0){};
      inline Int Map(Int i, Int j){
       Int proc = 0;
        if(pProcMap_!=NULL){
          proc = row2D_((*pProcMap_)[i],(*pProcMap_)[j]);
        }
        else{
          proc = row2D_(i,j);
        }
        return proc;
      }

      ~Row2D(){};
};

class Wrap2D: public Mapping{
  protected:
    inline Int wrap2D_(Int i, Int j) { Int p = (i/iBlockSize_)%iPRows_ + iPRows_*((j/iBlockSize_)%iPCols_); return p;}
  public:
      Wrap2D(Int aiNumProc, Int aiPRows, Int aiPCols, Int aiBlockSize = 1):Mapping(aiNumProc,aiPRows,aiPCols,aiBlockSize){};
      Wrap2D(Int aiNumProc, Int aiPRows, Int aiPCols, std::vector<Int> & aProcMap, Int aiBlockSize = 1):Mapping(aiNumProc,aiPRows,aiPCols,aProcMap,aiBlockSize){};
      //Wrap2D(Row2D & C):Mapping(C){}
      Wrap2D(Wrap2D & C):Mapping(C){};
      Wrap2D():Mapping(0,0,0,0){};
      inline Int Map(Int i, Int j){
       Int proc = 0;
        if(pProcMap_!=NULL){
          proc = wrap2D_((*pProcMap_)[i],(*pProcMap_)[j]);
        }
        else{
          proc = wrap2D_(i,j);
        }
        return proc;
      }
      ~Wrap2D(){};
};

class Wrap2DForced: public Mapping{
  protected:
    inline Int wrap2D_(Int i, Int j) { Int p = (i!=j)?(i/iBlockSize_)%iPRows_ + iPRows_*((j/iBlockSize_)%iPCols_):j%iNumProc_; return p;}
  public:
      Wrap2DForced(Int aiNumProc, Int aiPRows, Int aiPCols, Int aiBlockSize = 1):Mapping(aiNumProc,aiPRows,aiPCols,aiBlockSize){};
      Wrap2DForced(Int aiNumProc, Int aiPRows, Int aiPCols, std::vector<Int> & aProcMap, Int aiBlockSize = 1):Mapping(aiNumProc,aiPRows,aiPCols,aProcMap,aiBlockSize){};
      //Wrap2D(Row2D & C):Mapping(C){}
      Wrap2DForced(Wrap2DForced & C):Mapping(C){};
      Wrap2DForced():Mapping(0,0,0,0){};
      inline Int Map(Int i, Int j){
       Int proc = 0;
        if(pProcMap_!=NULL){
          proc = wrap2D_((*pProcMap_)[i],(*pProcMap_)[j]);
        }
        else{
          proc = wrap2D_(i,j);
        }
        return proc;
      }
      ~Wrap2DForced(){};
};





class AntiDiag2D: public Mapping{
  protected:
    inline Int antidiag2D_(Int i, Int j) { Int p = (i/iBlockSize_+j/iBlockSize_)%iNumProc_; return p;}
  public:
      AntiDiag2D(Int aiNumProc, Int aiPRows, Int aiPCols, Int aiBlockSize = 1):Mapping(aiNumProc,aiPRows,aiPCols,aiBlockSize){};
      AntiDiag2D(Int aiNumProc, Int aiPRows, Int aiPCols, std::vector<Int> & aProcMap, Int aiBlockSize = 1):Mapping(aiNumProc,aiPRows,aiPCols,aProcMap,aiBlockSize){};
      //AntiDiag2D(Row2D & C):Mapping(C){}
      AntiDiag2D(AntiDiag2D & C):Mapping(C){};
      AntiDiag2D():Mapping(0,0,0,0){};
      inline Int Map(Int i, Int j){
       Int proc = 0;
        if(pProcMap_!=NULL){
          proc = antidiag2D_((*pProcMap_)[i],(*pProcMap_)[j]);
        }
        else{
          proc = antidiag2D_(i,j);
        }
        return proc;
      }
      ~AntiDiag2D(){};
};



class Col2D: public Mapping{
  protected:
    inline Int col2D_(Int i, Int j) { Int p = (j/iBlockSize_)%iNumProc_; return p;}
  public:
      Col2D(Int aiNumProc, Int aiPRows, Int aiPCols, Int aiBlockSize = 1):Mapping(aiNumProc,aiPRows,aiPCols,aiBlockSize){};
      Col2D(Int aiNumProc, Int aiPRows, Int aiPCols, std::vector<Int> & aProcMap, Int aiBlockSize = 1):Mapping(aiNumProc,aiPRows,aiPCols,aProcMap,aiBlockSize){};
      //Col2D(Col2D & C):Mapping(C.iNumProc_,C.iPRows_,C.iPCols_,C.iBlockSize_){};
      Col2D(Col2D & C):Mapping(C){};
      Col2D():Mapping(0,0,0,0){};
      inline Int Map(Int i, Int j){
       Int proc = 0;
        if(pProcMap_!=NULL){
          proc = col2D_((*pProcMap_)[i],(*pProcMap_)[j]);
        }
        else{
          proc = col2D_(i,j);
        }
        return proc;
      }

      ~Col2D(){};
};





}

#endif //_MAPPING_HPP_
