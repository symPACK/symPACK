# !/bin/python

import SparseUtil.HB as hb
import math
import sys



file = sys.argv[1]
(m,n,nnz,struct,colptr,rowind,nzval) = hb.read(file)
(m,n,nnz,I,J,NZs) = hb.CSCtoCOO(m,n,nnz,colptr,rowind,nzval)


type = "real"
if struct[0] == 'c':
  type = "complex"

symm = "symmetric"
if struct[1] == 'u':
  symm = "unsymmetric"

print "%%%%MatrixMarket matrix coordinate %s %s" % (type,symm)
print "%d %d %d" % (m,n,nnz)

for idx,row in enumerate(I):
    col = J[idx]
    nz = NZs[idx]
    print "%d %d %.13e" % (row,col,nz)



