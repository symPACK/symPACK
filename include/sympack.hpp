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
#ifdef CUDA_MODE
  int symPACK_cuda_setup(symPACK::symPACKOptions optionsFact);
#endif
  int symPACK_Finalize();
  int symPACK_Rank(int * rank);

#ifdef __cplusplus
}



#endif

#endif // _SYMPACK_HEADER_
