/// @file FBMatrix.hpp
/// @brief Matrix structure and methods for the FAN-BOTH method.
/// @author Mathias Jacquelin
/// @date 2014-01-21
#ifndef _FBMATRIX_MPI_HPP_ 
#define _FBMATRIX_MPI_HPP_


#ifdef UPCXX
#define __SAVE_UPCXX__
#undef UPCXX
#endif


#include "Environment.hpp"
#include "FBMatrix.hpp"


using namespace std;

namespace LIBCHOLESKY{


  class MPIGrid{
    public:
    // Data
    MPI_Comm    comm;
    MPI_Comm    rowComm;
    MPI_Comm    colComm;
    Int         mpirank;
    Int         mpisize; 
    Int         numProcRow;
    Int         numProcCol;

    // Member function
    MPIGrid( MPI_Comm Bcomm, int nprow, int npcol );
    ~MPIGrid();
  };

  class FBMatrix_mpi : public FBMatrix{
    protected:

      virtual bool lastUpdate(Int j, Int i);

    public:

      //Parameters
      MPIGrid * mpigrid;

      void Initialize(MPIGrid & grid);

      virtual void Allocate(Int np,Int pn, Int pblksize);

      virtual void Distribute( DblNumMat & Aorig);

      virtual void Gather( DblNumMat & Adest);

      virtual void NumericalFactorization();

#ifdef DRAW_GRAPH
      LogFile * graphfileptr;
      virtual ~FBMatrix_mpi();
      void Draw_Graph();
#endif
  };



}


#ifdef __SAVE_UPCXX__
#define UPCXX
#endif


#endif
