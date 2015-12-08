#ifndef _MMD_HEADER_
#define _MMD_HEADER_

#include "ngchol/Environment.hpp"

namespace LIBCHOLESKY {

extern "C" {

void FORTRAN(ordmmd)( Int * neqns , Int * nadj  , Int * xadj  ,
         Int * adjncy, Int * invp  , Int * perm  , Int * iwsiz ,
                         Int * iwork , Int * nofsub, Int * iflag);


void FORTRAN(amdbar) (Int * N, Int * PE, Int * IW, Int * LEN, Int * IWLEN, Int * PFREE, Int * NV, Int * NEXT, Int *
         LAST, Int * HEAD, Int * ELEN, Int * DEGREE, Int * NCMPA, Int * W, Int * IOVFLO);

void FORTRAN(boxnd) (Int * P, Int * Q, Int * R, Int * IPERM, Int * WORK, Int * IERROR);

 }


}
#endif //_MMD_HEADER_

