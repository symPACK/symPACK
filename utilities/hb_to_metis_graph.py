# !/bin/python

import SparseUtil.HB as hb
import math
import sys



file = sys.argv[1]
(m,n,nnz,struct,colptr,rowind,nzval) = hb.read(file)

assert(struct[1]=='s')

(m,n,nnz2,struct,colptr,rowind,_) = hb.extend(m,n,nnz,struct,colptr,rowind,nzval)


print n," ",nnz-n
for col in range(1, n+1):
  colbeg = colptr[col-1]
  colend = colptr[col]-1

  for i in range(colbeg,colend+1):
    row = rowind[i-1]
    if(row!=col):
      print row," ",
  print ''



