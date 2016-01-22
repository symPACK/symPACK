/// @file FBMatrix.hpp
/// @brief Matrix structure and methods for the FAN-BOTH method.
/// @author Mathias Jacquelin
/// @date 2014-01-21
#ifndef _FBMATRIX_HPP_ 
#define _FBMATRIX_HPP_

#include "sympack/Environment.hpp"
#include "sympack/NumVec.hpp"
#include "sympack/NumMat.hpp"
#include "sympack/DistSparseMatrix.hpp"

#include <vector>

#ifndef MAP
#define MAP modwrap2D
#endif

using namespace std;

namespace SYMPACK{



  class FBMatrix{
    public:

      //Parameters
      Int blksize;

      std::vector<DblNumMat *> AchunkLower;
      std::vector<DblNumMat *> WLower;
  
      //lock should be initialized with the number of contribution a block column is receiving
      IntNumVec AggLock;

      Int n;
      Int pcol, prow;
      Int np, iam;

      FBMatrix();
      virtual ~FBMatrix();

      void ClearTmp();

      inline Int row2D(Int i, Int j) {return (i/blksize)%np;}
      inline Int col2D(Int i, Int j) {return (j/blksize)%np;}

      inline Int chevron2D(Int i, Int j) {return (min(i,j)/blksize)%np;}
      inline Int antichevron2D(Int i, Int j) {return (max(i,j)/blksize)%np;}

      inline Int diag2D(Int i, Int j) {return (abs(i-j)/blksize)%np;}
      inline Int antidiag2D(Int i, Int j) {return ((i+j)/blksize)%np;}


      inline Int modwrap2D(Int i, Int j) {return min(i/blksize,j/blksize)%prow + prow*floor((double)(max(i/blksize,j/blksize)%np)/(double)prow);}
      inline Int modwrap2Dns(Int i, Int j) {return(i/blksize)%prow + prow*floor((double)((j/blksize)%np)/(double)prow);}

      inline Int global_col_to_local(Int j){ return ((j)/(pcol*blksize))*blksize; }

      virtual void Allocate(Int & np, Int pn, Int pblksize);

      virtual void Distribute( DblNumMat & Aorig) = 0;

      virtual void Gather( DblNumMat & Adest){};

      virtual void WaitFactorization();

      virtual void NumericalFactorization() = 0;



      void Aggregate(Int j, DblNumMat &DistW);

      void Factor(Int j);
      void Factor_ref(Int j);

      void Update(Int j, Int i, DblNumMat & Factor);



  };

}





#endif
