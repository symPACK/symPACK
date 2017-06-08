/*
	 Copyright (c) 2016 The Regents of the University of California,
	 through Lawrence Berkeley National Laboratory.  

   Author: Mathias Jacquelin
	 
   This file is part of symPACK. All rights reserved.

	 Redistribution and use in source and binary forms, with or without
	 modification, are permitted provided that the following conditions are met:

	 (1) Redistributions of source code must retain the above copyright notice, this
	 list of conditions and the following disclaimer.
	 (2) Redistributions in binary form must reproduce the above copyright notice,
	 this list of conditions and the following disclaimer in the documentation
	 and/or other materials provided with the distribution.
	 (3) Neither the name of the University of California, Lawrence Berkeley
	 National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
	 be used to endorse or promote products derived from this software without
	 specific prior written permission.

	 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
	 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
	 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
	 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
	 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

	 You are under no obligation whatsoever to provide any bug fixes, patches, or
	 upgrades to the features, functionality or performance of the source code
	 ("Enhancements") to anyone; however, if you choose to make your Enhancements
	 available either publicly, or directly to Lawrence Berkeley National
	 Laboratory, without imposing a separate written license agreement for such
	 Enhancements, then you hereby grant the following license: a non-exclusive,
	 royalty-free perpetual license to install, use, modify, prepare derivative
	 works, incorporate into other computer software, distribute, and sublicense
	 such enhancements or derivative works thereof, in binary and source code form.
*/
#ifndef _MAPPING_HPP_
#define _MAPPING_HPP_

#include "sympack/Environment.hpp"
#include "sympack/LoadBalancer.hpp"


#ifndef MAPCLASS
#define MAPCLASS Modwrap2D
#endif

namespace symPACK{


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
      virtual void Update(std::vector<Int> & aProcMap){
        if(pProcMap_==NULL){
          pProcMap_ = new std::vector<Int>(aProcMap);
        }
        else{
          pProcMap_->clear();
          pProcMap_->insert(aProcMap.begin(),aProcMap.end(),pProcMap_->begin());
        }
      }

      void Dump(Int n){
        logfileptr->OFS()<<"Resulting mapping: "<<std::endl;
        logfileptr->OFS()<<"np = "<<iNumProc_<<std::endl;
        logfileptr->OFS()<<"prows = "<<iPRows_<<std::endl;
        logfileptr->OFS()<<"pcols = "<<iPCols_<<std::endl;
        for(Int i =0;i<n;++i){
          for(Int j = 0;j<=i;++j){
            logfileptr->OFS()<<Map(i,j)<<" ";
          }
          logfileptr->OFS()<<std::endl;
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
      virtual Int Map(Int i, Int j)=0;
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













class TreeMapping: public Mapping{
  protected:
    TreeLoadBalancer * pBalancer_;

    TreeLoadBalancer::ProcGroup * active_group_;



    Int iNumProc_;
    Int iPRows_;
    Int iPCols_;
    Int iBlockSize_;
  public:
      virtual void Update(TreeLoadBalancer * apBalancer){
        pBalancer_ = apBalancer;
      }

 
      TreeMapping(Int aiNumProc, Int aiPRows, Int aiPCols, Int aiBlockSize = 1):iNumProc_(aiNumProc),iPRows_(aiPRows),iPCols_(aiPCols),iBlockSize_(aiBlockSize),pBalancer_(NULL){};
//      TreeMapping(Int aiNumProc, Int aiPRows, Int aiPCols, std::vector<Int> & aProcMap, Int aiBlockSize = 1):
//          iNumProc_(aiNumProc),iPRows_(aiPRows),iPCols_(aiPCols),iBlockSize_(aiBlockSize){
//          pProcMap_ = new std::vector<Int>(aProcMap);
//      };
//
//      TreeMapping(TreeMapping & C){
//        iNumProc_   = C.iNumProc_;
//        iPRows_     = C.iPRows_;
//        iPCols_     = C.iPCols_;
//        iBlockSize_ = C.iBlockSize_;
//        if(C.pProcMap_!=NULL){
//          pProcMap_ = new std::vector<Int>(*C.pProcMap_);
//        }
//      }
      TreeMapping(){
        pBalancer_=NULL;
      };

      

      virtual Int Map(Int i, Int j){
i = i+1;
j = j+1;
        assert(pBalancer_!=NULL);
        //call GroupMap on the correct group, which is the one corresponding to j
        std::vector<Int> & groupIdx = pBalancer_->GroupIdx();
        std::vector<Int> & groupWorker = pBalancer_->GroupWorker();
        std::vector<TreeLoadBalancer::ProcGroup> & levelGroups = pBalancer_->LevelGroups();

#ifdef VERBOSE
      logfileptr->OFS()<<"MAP("<<i<<","<<j<<")"<<std::endl;
      logfileptr->OFS()<<"ETREE:"<<std::endl<<pBalancer_->SupETree()<<std::endl;
#endif

        //We should build a procmap local to the subtree rooted in j
        Int parentJ = pBalancer_->SupETree().PostParent(j-1);

        Int node = 0;
        for(Int I = 1; I<j;I++){
          Int parent = pBalancer_->SupETree().PostParent(I-1);
          if(parent==parentJ){
            node=I;
          }
        }
        //node is now the first node of subtree rooted in j
        node = node+1;
assert(i>=node);
        Int treesize = j - node + 1;

        std::vector<Int> locProcMap(treesize);
        for(Int I = node;I<=j;++I){
          Int idx = groupIdx[I];
          Int worker = groupWorker[I];
          locProcMap[I-node] = worker;
        }


      Int iPRows = ceil(sqrt(treesize));
#ifdef VERBOSE
        logfileptr->OFS()<<"local procMap"<<locProcMap<<std::endl;

      for(Int myI=0; myI<treesize;myI++){
      for(Int myJ=0; myJ<=myI;myJ++){
      Int p =GroupMap(myI,myJ,iPRows,treesize);
        logfileptr->OFS()<<p<<" ";
      }
        logfileptr->OFS()<<std::endl;
      }
        logfileptr->OFS()<<std::endl;
#endif

      //Int myI =locProcMap[ i-node ];
      //Int myJ =locProcMap[ j-node ];
      Int myI = i-node ;
      Int myJ = j-node ;
      Int p =GroupMap(myI-1,myJ-1,iPRows,treesize);
 //std::min(myI,myJ)%iPRows + iPRows*floor((double)(std::max(myI,myJ)%treesize)/(double)iPRows);

#ifdef VERBOSE
        logfileptr->OFS()<<"local i "<<myI<<std::endl;
        logfileptr->OFS()<<"local j "<<myJ<<std::endl;
        logfileptr->OFS()<<"local p "<<p<<std::endl;
assert(p<treesize);
        logfileptr->OFS()<<"global p "<<locProcMap[p]<<std::endl;
#endif
        //convert i and j to "local indices" within the group.
        //Subtrees are postordered to we just need to remove the smallest
        //value, which is j - active_group_->Ranks().size() 

//        Int local
//        Int proc = GroupMap();
//        if(pProcMap_!=NULL){
//          proc = antidiag2D_((*pProcMap_)[i],(*pProcMap_)[j]);
//        }
//        else{
//          proc = antidiag2D_(i,j);
//        }

        return locProcMap[p];
      }

  protected:
      virtual Int GroupMap(Int i, Int j,Int iPRows, Int iNp) =0;
};




class ModWrap2DTreeMapping: public TreeMapping{
  public:
     ModWrap2DTreeMapping(Int aiNumProc, Int aiPRows, Int aiPCols, Int aiBlockSize = 1):TreeMapping(aiNumProc, aiPRows, aiPCols, aiBlockSize ){};
      ModWrap2DTreeMapping():TreeMapping(){
      };

  protected:
      virtual Int GroupMap(Int i, Int j,Int iPRows, Int iNp){
      return std::min(iNp-1,(Int)(std::min(i,j)%iPRows + iPRows*floor((double)(std::max(i,j)%iNp)/(double)iPRows)));
}
};












}

#endif //_MAPPING_HPP_
