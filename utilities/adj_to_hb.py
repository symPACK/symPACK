# !/bin/python

import fileinput
import math
import sys

n = -1
m = -1
nnz = -1
lennnz = -1
colptr = []
rowind = []
for line in fileinput.input():
    if fileinput.filelineno() == 1:
        if n!=-1 and m!=-1 and nnz!=-1:
            #print prev HB file
            print m,n,nnz,0
            print colptr,
            print rowind,
            print rowind,
        colptr = [ int(x) for x in line.split() ]
        nnz = colptr[-1]-1
        n = len(colptr) - 1
    else:
        rowind = [ int(x) for x in line.split()]
        
        m = max( rowind )
      


if n!=-1 and m!=-1 and nnz!=-1:
    #add diagonal entries to colptr and rowind
    nnz = nnz + n
    old_colptr = colptr
    colptr = [ x+i for i,x in enumerate(colptr) ]
    old_rowind = rowind
    rowind = []
    for col in range(1, n+1):
      colbeg = old_colptr[col-1]
      colend = old_colptr[col]-1
      rowind.append(col)
      for i in range(colbeg,colend+1):
        row = old_rowind[i-1]
        rowind.append(row)

      colbeg = colptr[col-1]
      colend = colptr[col]-1
      rowind[colbeg-1:colend] = sorted(rowind[colbeg-1:colend])

    old_colptr = colptr
    old_rowind = rowind
    colptr = []
    rowind = []
    colptr.append(1)
    for col in range(1, n+1):
      colbeg = old_colptr[col-1]
      colend = old_colptr[col]-1
      for i in range(colbeg,colend+1):
        row = old_rowind[i-1]
        if(row>=col):
          rowind.append(row)
      colptr.append(len(rowind))
    nnz = colptr[-1]-1

    #print prev HB file
    lennnz = len("%d" % nnz)+1
    lenm = len("%d" % m)+1
    lenmnz = len("%0.12E" % float(m))+1
    numElemColptr = math.floor(80 / (lennnz))
    numElemRowind = math.floor(80 / (lenm))
    numElemNzval = math.floor(80 / (lenmnz))
    numRowColptr = math.ceil(len(colptr) / numElemColptr)
    numRowRowind = math.ceil(len(rowind) / numElemRowind)
    numRowNzval = math.ceil(len(rowind) / numElemNzval)

    print "Adjacency to HB converted file                                         |23      "
    print "%14d%14d%14d%14d" % (numRowColptr+numRowRowind+numRowNzval,numRowColptr,numRowRowind,numRowNzval)
    print "rsa           %14d%14d%14d%14d" % (m,n,nnz,0)
    print "%s%s%s" % ( '{:<16}'.format("(%dI%d)" %(numElemColptr,lennnz)), '{:<16}'.format("(%dI%d)" %(numElemRowind,lenm)) , '{:<20}'.format("(%dE%d.12)" %(numElemNzval,lenmnz)))
    for idx,val in enumerate(colptr):
       if(idx!=0 and idx % numElemColptr == 0):
          sys.stdout.write('\n')
       sys.stdout.write(( "%%%dd" % lennnz ) % val)
    sys.stdout.write('\n')

    for idx,val in enumerate(rowind):
       if(idx!=0 and idx % numElemRowind == 0):
          sys.stdout.write('\n')
       sys.stdout.write(( "%%%dd" % lenm ) % val)
    sys.stdout.write('\n')

    for idx,val in enumerate(rowind):
       if(idx!=0 and idx % numElemNzval == 0):
          sys.stdout.write('\n')
       sys.stdout.write(( "%%%d.12E" % lenmnz ) % float(val))
    sys.stdout.write('\n')

    sys.stdout.flush()

#    print rowind,
#    print rowind,
