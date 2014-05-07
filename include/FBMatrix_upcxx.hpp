/// @file FBMatrix.hpp
/// @brief Matrix structure and methods for the FAN-BOTH method.
/// @author Mathias Jacquelin
/// @date 2014-01-21
#ifndef _FBMATRIX_UPCXX_HPP_ 
#define _FBMATRIX_UPCXX_HPP_

#ifndef UPCXX
#define UPCXX
#define __SAVE_UPCXX__
#endif

#include <upcxx.h>

#include "Environment.hpp"
#include "NumMat_upcxx.hpp"
#include "FBMatrix.hpp"


using namespace std;

namespace LIBCHOLESKY{




  class FBMatrix_upcxx: public FBMatrix{
    public:

      //Parameters
      double aggregate_comm_time;
      double factor_comm_time;
      Int outstdAggreg=0;
      Int outstdUpdate=0;
      Int prefetch;
      std::vector< upcxx::global_ptr<FBMatrix_upcxx> > * pRemoteObjPtrs;

      FBMatrix_upcxx();
      virtual ~FBMatrix_upcxx();

      void Initialize();

      virtual void Allocate(Int np, Int pn, Int pblksize);

      virtual void Distribute( DblNumMat & Aorig);
//      virtual void DistributeSparse( DistSparseMatrix<double> & Aorig);

      virtual void Gather( DblNumMat_upcxx & Adest);

      virtual void WaitFactorization();

      virtual void NumericalFactorization();
      virtual void NumericalFactorizationLoop();

      bool lastUpdate(Int j, Int i);
  };




}

#ifdef __SAVE_UPCXX__
#undef UPCXX
#endif

#endif
