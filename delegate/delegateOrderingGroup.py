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
    newPerm = []
    k = 0
    while k <= depth:
        start = 0
        moved = False
        while start < N.value:
            #originalPerm maintains the original permutation for this iteration
            #permList is the one to play with and permute
 #           print("Start is " + str(start))
            originalPerm = makePermList(bestPerm)
            permList = makePermList(bestPerm)


            #StartNeighborhood is a set of the neighbors we are looking at.
            #orderedNeighborhood gets them in the order they appeared in the originalPerm
            startNeighborhood = getNeighborhood(bestPerm[start], k)
            orderedNeighborhood = []
            for i in range(N.value):
                if bestPerm[i] in startNeighborhood:
                    orderedNeighborhood.append(bestPerm[i])
            if len(startNeighborhood) != len(orderedNeighborhood):
                print("Neighborhoodlengths are wrong")
#            print(orderedNeighborhood)

            #lastNeighborIndex is the index of the last neighbor. Where we are going to start putting the neighborhood.
            lastNeighborIndex = originalPerm.index(orderedNeighborhood[-1])
            #print(orderedNeighborhood)
            movingAfter = 0
            while movingAfter < (len(orderedNeighborhood) - 1):
                while (permList.index(orderedNeighborhood[movingAfter]) + 1) != permList.index(orderedNeighborhood[movingAfter+1]):
                    #print("Moving Together: " + str(permList))
                    indexFirstNeighbor = permList.index(orderedNeighborhood[0])
                    permList.insert(indexFirstNeighbor, permList.pop(permList.index(orderedNeighborhood[movingAfter]) + 1)) 

                    #If this position is better keep it.
                    fillPermArray(currentPerm, permList)
                    currentCost = GetCost(N, Nnz, xAdj, Adj, currentPerm)
                    if currentCost < bestCost:
                        moved = True
                        newPerm = permList
                        bestCost = currentCost
                        bestPerm = currentPerm
                movingAfter += 1 

            #Try moving elements from after the neighborhood forward and check them.
            lastNeighbor = permList.index(orderedNeighborhood[-1])
            firstNeighbor = permList.index(orderedNeighborhood[0])
            #print("Before: " + str(permList))
            while lastNeighbor != N.value-1:
                permList.insert(firstNeighbor, permList.pop(lastNeighbor + 1))
                #print("After: " + str(permList))
                fillPermArray(currentPerm, permList)
                currentCost = GetCost(N, Nnz, xAdj, Adj, currentPerm)
                if currentCost < bestCost:
                    moved = True
                    newPerm = permList
                    bestCost = currentCost
                lastNeighbor = permList.index(orderedNeighborhood[-1])
                firstNeighbor = permList.index(orderedNeighborhood[0])

            #Update the permutations with the best one found this cycle and increase start.
            fillPermArray(currentPerm, newPerm)
            fillPermArray(bestPerm, newPerm)
            if moved == False:
                start += 1
            moved = False
        k += 1
    return bestCost        


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

