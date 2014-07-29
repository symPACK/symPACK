#Kevin Carlson
#Python Delegate Ordering Algorithm
#LBL Summer 2014

import sys
import ctypes

#Read the Adjacency File from the command line
adjacencyFile = bytes(sys.argv[1], 'utf-8')
inputPermutation = sys.argv[2]

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

#Set up variables and call GetCost for the input permutation.
xAdj = cXAdj.contents
Adj = cAdj.contents
N = cN.contents
Nnz = cNnz.contents
bestPerm = (ctypes.c_int * N.value)()
currentPerm = (ctypes.c_int * N.value)()
for i in range(len(inputPermutation.split())):
    bestPerm[i] = int(inputPermutation.split()[i])
    currentPerm[i] = int(inputPermutation.split()[i])
inputCost = GetCost(N, Nnz, xAdj, Adj, bestPerm)

#Set up Python Lists for the Adj and xAdj arays
pyAdj = []
for i in range(Nnz.value):
    pyAdj.append(Adj[i]) 

pyXAdj = []
for i in range(N.value + 1):
    pyXAdj.append(xAdj[i])


#Returns a list of all of the neightbors of 
#Value within depth range.
def getNeighbors(value, depth=1):
    start = pyXAdj[value - 1] - 1
    end = pyXAdj[value]-1 
    return pyAdj[start:end] 


bestCost = float("inf")


def permute():
    global bestCost
    global bestPerm
    global CurrentPerm
    for start in range(N.value):
        bestIndex = start
        currentIndex = start
        while currentIndex < N.value:
            temp = currentPerm[currentIndex]
            currentPerm[currentIndex] = currentPerm[start]
            currentPerm[start] = temp
            nextCost = GetCost(N, Nnz, xAdj, Adj, currentPerm)
            if nextCost < bestCost:
                bestIndex = currentIndex
                bestCost = nextCost
            temp = currentPerm[currentIndex]
            currentPerm[currentIndex] = currentPerm[start]
            currentPerm[start] = temp
            currentIndex += 1
        swap = currentPerm[bestIndex]
        currentPerm[bestIndex] = currentPerm[start]
        bestPerm[bestIndex] = bestPerm[start]
        currentPerm[start] = swap
        bestPerm[start] = swap
    return bestCost




def legalPerm(permutation):
    testSet = set()
    for i in range(N.value):
        testSet.add(permutation[i])
    if len(testSet) == N.value:
        return True
    return False
    


