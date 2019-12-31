#ifndef _SYMPACK_HEADER_
#define _SYMPACK_HEADER_

//#define SP_THREADS

#ifdef __cplusplus
#include  "sympack/Environment.hpp"
#include  "sympack/utility.hpp"
#include  "sympack/DistSparseMatrix.hpp"
#include  "sympack/symPACKMatrix.hpp"
#include  "sympack/Mapping.hpp"
#include  "sympack/CommTypes.hpp"
#include  "sympack/Ordering.hpp"
#include  "sympack/LogFile.hpp"
#endif

#include <mpi.h>






#ifdef __cplusplus
extern "C" {
#endif

extern bool libMPIInit;
extern bool libUPCXXInit;

  int symPACK_Init(int *argc=nullptr, char ***argv=nullptr);
  int symPACK_Finalize();
  int symPACK_Rank(int * rank);



  //returns an integer corresponding to a symPACK handle

  int symPACK_C_InitInstanceFloat(MPI_Comm ccomm, bool is2D = false);
  int symPACK_C_InitInstanceDouble(MPI_Comm ccomm, bool is2D = false);
  int symPACK_C_InitInstanceComplex(MPI_Comm ccomm, bool is2D = false);
  int symPACK_C_InitInstanceDoubleComplex(MPI_Comm ccomm, bool is2D = false);
  
  int symPACK_InitInstanceFloat(MPI_Fint * Fcomm, int is2D = 0);
  int symPACK_InitInstanceDouble(MPI_Fint * Fcomm, int is2D = 0);
  int symPACK_InitInstanceComplex(MPI_Fint * Fcomm, int is2D = 0);
  int symPACK_InitInstanceDoubleComplex(MPI_Fint * Fcomm, int is2D = 0);
  
  void symPACK_SymbolicFactorize(int * sp_handle, int * n, int * colptr , int * rowind);
  
  void symPACK_DistributeFloat(int * sp_handle, float * nzvals);
  void symPACK_DistributeDouble(int * sp_handle, double * nzvals);
  void symPACK_DistributeComplex(int * sp_handle, float * nzvals);
  void symPACK_DistributeDoubleComplex(int * sp_handle, double * nzvals);
  
  void symPACK_NumericalFactorize(int * sp_handle);
  
  void symPACK_NumericalSolveFloat(int * sp_handle, int * nrhs, float * rhs);
  void symPACK_NumericalSolveDouble(int * sp_handle, int * nrhs, double * rhs);
  void symPACK_NumericalSolveComplex(int * sp_handle, int * nrhs, float * rhs);
  void symPACK_NumericalSolveDoubleComplex(int * sp_handle, int * nrhs, double * rhs);
  
  void symPACK_FinalizeInstance(int * sp_handle);




#ifdef __cplusplus
}



#endif

#endif // _SYMPACK_HEADER_
