#ifndef _MMD_HEADER_
#define _MMD_HEADER_

#include "ngchol/Environment.hpp"

namespace LIBCHOLESKY {

extern "C" {

void FORTRAN(ordmmd)( Int * neqns , Int * nadj  , Int * xadj  ,
         Int * adjncy, Int * invp  , Int * perm  , Int * iwsiz ,
                         Int * iwork , Int * nofsub, Int * iflag);

 }




}
#endif //_MMD_HEADER_

