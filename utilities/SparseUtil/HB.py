# !/bin/python

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import coo_matrix

#HB file driver
def read(file):
  n = 0
  m = 0
  nnz = 0
  colptr = []
  rowind = []
  nzval  = []
  
  
  with open(file,'r') as f:
    for lineno,line in enumerate(f):
      if lineno==0:
        #skip title line
        continue
      elif lineno==1:
        #read info on numRows, numRowColPtr, numRowRowind, numRowNzval
        (numRows, numRowColptr, numRowRowind, numRowNzVal) = [t(s) for t,s in zip((int,int,int,int),line.split())]
      elif lineno==2:
        #read info on structure, m, n, nnz
        (struct, m, n, nnz, dummy) = [t(s) for t,s in zip((str,int,int,int,int),line.split())]
        #analyze struct
        struct = list(struct)
      elif lineno==3:
        #skip format line
        continue
      elif lineno>=4 and lineno<4+numRowColptr:
        colptr.extend(map(int, line.split()))
      elif lineno>=4+numRowColptr and lineno<4+numRowColptr+numRowRowind:
        rowind.extend(map(int, line.split()))
      elif lineno>=4+numRowColptr+numRowRowind and lineno<4+numRows:
        #read nzvals
        if struct[0]=='r' :
          nzval.extend(map(float, line.split()))
        else:
          nzval.extend(map(complex, line.split()))
  
  return(m,n,nnz,struct,colptr,rowind,nzval)

def CSCtoCOO(m,n,nnz,colptr,rowind,nzval):
  I=[0]*nnz
  J=[0]*nnz
  NZ = [0]*nnz

  tail = 0
  for col in range(1, n+1):
    colbeg = colptr[col-1]
    colend = colptr[col]-1
    for i in range(colbeg,colend+1):
      row = rowind[i-1]
      I[tail] = row
      J[tail] = col
      NZ[tail] = nzval[i-1]
      tail+=1
  return (m,n,nnz,I,J,NZ)

def COOtoCSC(m,n,nnz,I,J,NZ):
  
  rowinds= [ [] for i in range(0,n) ]
  nzvals= [ [] for i in range(0,n) ]

  for idx,row in enumerate(I):
    col = J[idx]
    nz = NZ[idx]
    rowinds[col-1].append(row)
    nzvals[col-1].append(nz)
     
  
  colptr= [1]
  colptr.extend([ len(x) for x in rowinds ])
  for i in range(1,len(colptr)):
    colptr[i] = colptr[i-1]+colptr[i]

  rowind = [item for sublist in rowinds for item in sublist]
  nzval = [item for sublist in nzvals for item in sublist]
  
  
  return (m,n,nnz,colptr,rowind,nzval)


def extend(m,n,nnz,struct,colptr,rowind,nzval,issorted=True):
  if(struct[1]=='u'):
    return (m,n,nnz,struct,colptr,rowind,nzval)
  
  #convert to COO as it is easier to extend
  (m2,n2,nnz2,I,J,NZ) = CSCtoCOO(m,n,nnz,colptr,rowind,nzval)

  
  newnnz = nnz-n
  eI  = [0]*(newnnz)
  eJ  = [0]*(newnnz)
  eNZ = [0]*(newnnz)
  tail = 0
  for idx,row in enumerate(I):
    col = J[idx]
    nz = NZ[idx]
    if(row>col):
      eI[tail] = col
      eJ[tail] = row
      eNZ[tail] = nz
      tail+=1

  I.extend(eI)
  J.extend(eJ)
  NZ.extend(eNZ)

  nnz=nnz+newnnz

  #convert back to CSC
  (m2,n2,nnz2,ecolptr,erowind,enzval) = COOtoCSC(m,n,nnz,I,J,NZ)

  if issorted:
    for col in range(1, n+1):
      colbeg = ecolptr[col-1]
      colend = ecolptr[col]-1

      #sort rowind and nzval the same way     
      erowind[colbeg-1:colend],enzval[colbeg-1:colend] = (list(t) for t in zip(*sorted(zip(erowind[colbeg-1:colend],enzval[colbeg-1:colend]))) )


  struct2 = struct
  struct2[1] = 'u'
  return (m2,n2,nnz2,struct2,ecolptr,erowind,enzval)
  
  


#def readToCSC(file):
#  (m,n,nnz,struct,colptr,rowind,nzval) = read(file)
#  A = csc_matrix((nzval,colptr,rowind),shape=(m,n),dtype=np.float)
#  return A
#
#  
#def extend(A):
#  B = A.tocsc()
#  print B
