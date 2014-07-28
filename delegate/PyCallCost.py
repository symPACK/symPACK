#Kevin Carlson
#Python Call Cost Interface
#LBL Summer 2014

import sys
import ctypes

#Read the Adjacency File from the command line
adjacencyFile = bytes(sys.argv[1], 'utf-8')
permutation = sys.argv[2]

#Import the C Functions from util.cpp and set the arguments and return types.
utilLib = ctypes.CDLL('/Users/kevin/Desktop/LBLrepo/ngchol/delegate/sharedLib.so')

GetCost = utilLib.GetCost
GetCost.restype = ctypes.c_double
GetCost.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)] 

GetCostPerCol  = utilLib.GetCostPerCol
GetCostPerCol.restype = ctypes.c_double
GetCostPerCol.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]

GetPrefixSum = utilLib.GetPrefixSum
GetPrefixSum.restype = None
GetPrefixSum.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]

ReadAdjacency = utilLib.ReadAdjacency
ReadAdjacency.restype = ctypes.c_int
ReadAdjacency.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.POINTER(ctypes.POINTER(ctypes.c_int)), ctypes.POINTER(ctypes.POINTER(ctypes.c_int)), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]


#Set up variables and call ReadAdjacency.
cAdjacencyFile = ctypes.create_string_buffer(adjacencyFile)
XAdjint = ctypes.c_int()
Adjint = ctypes.c_int()
cXAdj = ctypes.pointer(ctypes.pointer(XAdjint))
cAdj = ctypes.pointer(ctypes.pointer(Adjint))
cN = ctypes.pointer(ctypes.c_int(0))
cNnz = ctypes.pointer(ctypes.c_int(0))
ReadAdjacency(cAdjacencyFile, cXAdj, cAdj, cN, cNnz)

#Set up variables and call GetCost.
XAdj = cXAdj.contents
Adj = cAdj.contents
N = cN.contents
Nnz = cNnz.contents
perm = (ctypes.c_int * N.value)()
for i in range(len(permutation.split())):
    perm[i] = int(permutation.split()[i])
cost = GetCost(N, Nnz, XAdj, Adj, perm)
print(cost)

#Free the Adj and XAdj Arrays and exit.
utilLib.freeIntPointer(XAdj)
utilLib.freeIntPointer(Adj)
