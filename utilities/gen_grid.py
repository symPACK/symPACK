# !/bin/python

import math
import sys

k = int(sys.argv[1])
p = int(sys.argv[2])
stencil = int(sys.argv[3])


colptr = [0]*(k*p+1)
rowind = []

def toidx(i,j):
  #row major
  return (i-1)*p + j
  #col major
  #return (j-1)*k + i

#row major
#loop over rows
for i in range(1,k+1):
  #loop over columns
  for j in range(1,p+1):
    idx = toidx(i,j) 
    colptr[idx-1] = len(rowind)+1;
    
    #connect to node below
    if(i>1):
      rowind.append(toidx(i-1,j))
  
    #connect to node above
    if(i<k):
      rowind.append(toidx(i+1,j))

    #connect to node left
    if(j>1):
      rowind.append(toidx(i,j-1))

    #connect to node right
    if(j<p):
      rowind.append(toidx(i,j+1))

    if(stencil >=7):
      #connect to node bottom left
      if(j>1 and i>1):
        rowind.append(toidx(i-1,j-1))

      #connect to node top right
      if(j<p and i<k):
        rowind.append(toidx(i+1,j+1))

    if(stencil ==9):
      #connect to node bottom right
      if(j<p and i>1):
        rowind.append(toidx(i-1,j+1))

      #connect to node top left
      if(j>1 and i<k):
        rowind.append(toidx(i+1,j-1))
 
colptr[p*k] = len(rowind)+1;

#n = len(colptr)-1
#for col in range(1, n+1):
#  colbeg = colptr[col-1]
#  colend = colptr[col]-1
#  print "Column ",col,": ", rowind[colbeg-1:colend]

for k in colptr:
  print k,' ',
print ''

for k in rowind:
  print k,' ',
print ''

#print(' '.join('{}'.format(*k) for k in rowind))
#print colptr
#print rowind 



