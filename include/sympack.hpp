#ifndef _SYMPACK_HEADER_
#define _SYMPACK_HEADER_


#include  "sympack/Environment.hpp"
#include  "sympack/utility.hpp"
#include  "sympack/DistSparseMatrix.hpp"
#include  "sympack/symPACKMatrix.hpp"
#include  "sympack/Mapping.hpp"
#include  "sympack/CommTypes.hpp"
#include  "sympack/Ordering.hpp"
#include  "sympack/LogFile.hpp"

#ifdef __cplusplus
extern "C" {
#endif
  int symPACK_Init(int *argc=NULL, char ***argv=NULL);
  int symPACK_Finalize();
#ifdef __cplusplus
}
#endif

#endif // _SYMPACK_HEADER_
