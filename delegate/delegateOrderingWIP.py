#kevin carlson
#Python Delegate Ordering Algorithm
#LBL Summer 2014

#Correctly swaps and checks all permutations. Most up to date version as of 8/5/2014.

import sys
import ctypes

#Read the arguments supplied from the command line
if len(sys.argv) != 3 and len(sys.argv) != 4:
    print('Usage: Adjacency file, "ordering (space seperated)", optional update flag', file=sys.stderr)
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
currentPerm = (ctypes.c_int * N.value)()
for i in range(len(inputPermutation.split())):
    currentPerm[i] = int(inputPermutation.split()[i])
inputCost = GetCost(N, Nnz, xAdj, Adj, currentPerm)
bestMove = [[], inputCost]                  #permutation cost 
bestCost = inputCost

#Set up Python Lists for the Adj and xAdj arrays
pyAdj = []
for i in range(Nnz.value):
    pyAdj.append(Adj[i]) 

pyXAdj = []
for i in range(N.value + 1):
    pyXAdj.append(xAdj[i])

def get_reach (xAdj, Adj, elimNode, step, inputPerm):
    reachableSet = set()
    exploreSet = set()
    exploreSet |= set(getNeighbors(elimNode)) #Initialize exploredSet to neighbors of elimNode
    explored = set(exploreSet)

    permList = []
    splitPerm = inputPerm
    if type(inputPerm) == str:
        splitPerm = inputPerm.split()
    for i in range(len(splitPerm)):
        permList.append(int(splitPerm[i]))
        
    while (len(exploreSet) != 0):
        node = exploreSet.pop()
        if node == elimNode:
            continue
        elimStep = permList.index(node) + 1
        if elimStep > step-1:
            reachableSet.add(node)
        else:
            neighbors =set(getNeighbors(node))
            exploreSet |= ( neighbors - explored)
            explored |= neighbors
    return reachableSet

def updatePermAdj(perm):
    global pyAdj
    global pyXAdj
    newPyAdj = []
    newPyXAdj = []
    newAdjStruct(pyXAdj, Adj, perm, newPyXAdj, newPyAdj)
    pyAdj = newPyAdj
    pyXAdj = newPyXAdj


def newAdjStruct (xAdj, Adj, inputPerm, newXAdj, newAdj):
    permList = []
    if type(inputPerm) == str:
        inputPerm = inputPerm.split()
    for i in range(len(inputPerm)):
        permList.append(int(inputPerm[i]))
    newXAdj += [1]
    for el in range(len(permList)):
        newAdj.extend([permList[el]])
        sortedList = list( get_reach(xAdj, Adj, permList[el], el + 1, inputPerm ))
        sortedList.sort()
        newAdj.extend(sortedList)
        newXAdj += [ len(newAdj) + 1 ]
    #newXAdj += [len(newAdj) + 1]
         
# Function to complete delegate algorithm on
# the global variable bestPerm.
def permute(depth):
    if len(sys.argv) == 4 and sys.argv[3] == 'update':
       updatePermAdj(inputPermutation) 
    global bestMove
    global bestCost
    global currentPerm
    originalPerm = makePermList(currentPerm)
    k = 0
    #print(originalPerm)
    while k <= depth:
        start = 0
        while start < N.value:
            permList = makePermList(originalPerm)
            startNeighborhood = getNeighborhood(permList[start], k)
            #print(startNeighborhood)
            orderedNeighborhood = []
            for i in range(N.value):
                if permList[i] in startNeighborhood:
                    orderedNeighborhood.append(permList[i])
            if len(startNeighborhood) != len(orderedNeighborhood):
                print("Neighborhoodlengths are wrong")

            movingAfter = 0
            while movingAfter < (len(orderedNeighborhood) - 1):
                while (permList.index(orderedNeighborhood[movingAfter]) + 1) != permList.index(orderedNeighborhood[movingAfter+1]):
                    indexFirstNeighbor = permList.index(orderedNeighborhood[0])
                    permList.insert(indexFirstNeighbor, permList.pop(permList.index(orderedNeighborhood[movingAfter]) + 1)) 
                    #If this position is better keep it.
                    fillPermArray(currentPerm, permList)
                    #print(permList)
                    currentCost = GetCost(N, Nnz, xAdj, Adj, currentPerm)
                    if currentCost < bestMove[1]:
                        bestCost = currentCost
                        bestMove [1]  = currentCost
                        bestMove[0] = list(permList)
                movingAfter += 1
            lastNeighbor = permList.index(orderedNeighborhood[-1])
            firstNeighbor = permList.index(orderedNeighborhood[0])
            while lastNeighbor != N.value-1:
                permList.insert(firstNeighbor, permList.pop(lastNeighbor + 1))
                fillPermArray(currentPerm, permList)
                #print(permList)
                currentCost = GetCost(N, Nnz, xAdj, Adj, currentPerm)
                if currentCost < bestMove[1]:
                    bestCost = currentCost
                    bestMove[1] = currentCost
                    bestMove[0] = list(permList)
                lastNeighbor = permList.index(orderedNeighborhood[-1])
                firstNeighbor = permList.index(orderedNeighborhood[0])
            start += 1
        k += 1
    fillPermArray(currentPerm, bestMove[0])     
    return bestCost        

def fullPermute(depth):
    lastCost = inputCost
    while True:
        nextCost = permute(depth)
        print(nextCost)
        if nextCost == lastCost:
            return nextCost
        else: lastCost = nextCost


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
    correctSum = 0
    permSum = 0
    for i in range (1, N.value + 1):
        correctSum += i
    for el in range (N.value):
        permSum += permutation[el]
    if permSum == correctSum:
        return True
    return False


# Utility function to print a permutation.
def printPerm(permutation):
    for i in range(N.value):
        print(permutation[i])

