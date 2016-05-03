#ifndef _MMD_HEADER_
#define _MMD_HEADER_

#include "sympack/Environment.hpp"

namespace SYMPACK {



extern "C" {

#if MMDIDXSIZE==64
  typedef    int64_t   MMDInt;
#else
  typedef    int32_t   MMDInt;
#endif

#if AMDIDXSIZE==64
  typedef    int64_t   AMDInt;
#else
  typedef    int32_t   AMDInt;
#endif

void FORTRAN(ordmmd)( MMDInt * neqns , MMDInt * nadj  , MMDInt * xadj  ,
         MMDInt * adjncy, MMDInt * invp  , MMDInt * perm  , MMDInt * iwsiz ,
                         MMDInt * iwork , MMDInt * nofsub, MMDInt * iflag);


void FORTRAN(amdbar) (AMDInt * N, AMDInt * PE, AMDInt * IW, AMDInt * LEN, AMDInt * IWLEN, AMDInt * PFREE, AMDInt * NV, AMDInt * NEXT, AMDInt *
         LAST, AMDInt * HEAD, AMDInt * ELEN, AMDInt * DEGREE, AMDInt * NCMPA, AMDInt * W, AMDInt * IOVFLO);

void FORTRAN(boxnd) (Int * P, Int * Q, Int * R, Int * IPERM, Int * WORK,Int * WORKSZ, Int * IERROR);
void FORTRAN(gridnd) (Int * P, Int * Q, Int * IPERM, Int * WORK,Int * WORKSZ, Int * IERROR);


 }


}
#endif //_MMD_HEADER_

