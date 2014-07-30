#Kevin Carlson
#Python Delegate Ordering Algorithm
#LBL Summer 2014

import sys
import ctypes

#Read the arguments supplied from the command line
if len(sys.argv) != 3:
    print('Usage: Adjacency file, "ordering (space seperated)"', file=sys.stderr)
    sys.exit()
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
bestCost = float("inf")

#Set up Python Lists for the Adj and xAdj arrays
pyAdj = []
for i in range(Nnz.value):
    pyAdj.append(Adj[i]) 

pyXAdj = []
for i in range(N.value + 1):
    pyXAdj.append(xAdj[i])


# Function to complete delegate algorithm on
# the global variable bestPerm.
def permute(depth):
    global bestCost
    global bestPerm
    global currentPerm
    while depth >= 1:
        for start in range(N.value):
            bestIndex = start
            currentIndex = start
            while currentIndex < N.value:
               
                #Swap Neighborhood
                oldPerm = makePermList(currentPerm)
                modifiedList = makePermList(currentPerm)
                neighborList = getNeighborhood(currentPerm[currentIndex],depth)
                for el in oldPerm:
                    if el in neighborList:
                        modifiedList.append(modifiedList.pop(modifiedList.index(el)))
                fillPermArray(currentPerm, modifiedList)
                
                nextCost = GetCost(N, Nnz, xAdj, Adj, currentPerm)
                if nextCost < bestCost:
                    bestIndex = currentIndex
                    bestCost = nextCost

                #Restore Neighborhood
                fillPermArray(currentPerm,oldPerm)
                currentIndex += 1

            #Swap best Neighborhood in bestPerm and currentPerm
            tempList = makePermList(currentPerm)
            modifiedList = makePermList(currentPerm)
            neighborList = getNeighborhood(currentPerm[bestIndex], depth)
            for el in tempList:
                if el in neighborList:
                    modifiedList.append(modifiedList.pop(modifiedList.index(el)))
            fillPermArray(currentPerm, modifiedList)
            fillPermArray(bestPerm, modifiedList)
            
        depth -= 1
    return bestCost        
#    global bestCost
#    global bestPerm
#    global CurrentPerm
#    for start in range(N.value):
#        bestIndex = start
#        currentIndex = start
#        while currentIndex < N.value:
#            temp = currentPerm[currentIndex]
#            currentPerm[currentIndex] = currentPerm[start]
#            currentPerm[start] = temp
#            nextCost = GetCost(N, Nnz, xAdj, Adj, currentPerm)
#            if nextCost < bestCost:
#                bestIndex = currentIndex
#                bestCost = nextCost
#            temp = currentPerm[currentIndex]
#            currentPerm[currentIndex] = currentPerm[start]
#            currentPerm[start] = temp
#            currentIndex += 1
#        swap = currentPerm[bestIndex]
#        currentPerm[bestIndex] = currentPerm[start]
#        bestPerm[bestIndex] = bestPerm[start]
#        currentPerm[start] = swap
#        bestPerm[start] = swap
#    return bestCost



def makePermList(arg):
    pylist = []
    for i in range(N.value):
        pylist.append(arg[i])
    return pylist

def fillPermArray(array, inputList):
    for i in range(len(inputList)):
        array[i] = inputList[i]
    return array


#Returns a list of all of the immediate neighbors.
def getNeighbors(value):
    start = pyXAdj[value - 1] - 1
    end = pyXAdj[value]-1 
    return pyAdj[start:end]

#Returns a list of all of the neighbors up to depth.
def getNeighborhood(value, depth):
    if depth == 0: 
        return set([value])
    if depth == 1:
        return set(getNeighbors(value))
    else:
        neighborlist = getNeighbors(value)  #List of all of the neighbors
        neighborhood = set(neighborlist)
        for el in neighborlist:
            neighborhood |= getNeighborhood(el, depth - 1)  #| is union of sets
        return neighborhood
        
        



# Utility function to determine legality of
# an ordering permutation.
def legalPerm(permutation):
    testSet = set()
    for i in range(N.value):
        if i in range(1, N.value+1):
            testSet.add(permutation[i])
    if len(testSet) == N.value:
        return True
    return False

# Utility function to print a permutation.
def printPerm(permutation):
    for i in range(N.value):
        print(permutation[i])

