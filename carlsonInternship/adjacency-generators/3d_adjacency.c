#include <stdlib.h>
#include <stdio.h>

/*3D representation for adjacency arrays. Takes 4 integer arguments
  from the command line LENGTH WIDTH DEPTH and OP. If OP is 27, a 27
  point operation will be used, if it is anything else it will default
  to a 7 point operation.
  Author: Kevin Carlson */

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Invalid arguments, please give height, width, depth and operator.\n");
        return 1;
    }

    //Declare the ints needed for generating the arrays.
    unsigned int height;
    unsigned int width;
    unsigned int depth;
    unsigned int op;
    unsigned int* nodeArray;
    unsigned int* neighborArray;

    //Read the input to the program.
    height = atoi(argv[1]);
    width = atoi(argv[2]);
    depth = atoi (argv[3]);
    op = atoi(argv[4]);

    //Create and Malloc necessary arrays.
    int numNodes;
    int numNeighbors;
    int singleNeighbors;
    int numPerPlane;
    singleNeighbors = op - 1;
    numNodes = height * width * depth;
    numNeighbors = numNodes * singleNeighbors;
    numPerPlane = height * width;
    nodeArray = (unsigned int*) malloc((numNodes+1)*sizeof(unsigned int));
    neighborArray = (unsigned int*) malloc(numNeighbors*sizeof(unsigned int));

    //Fill the Arrays
    int nodeCounter = 0;
    int neighborCounter = 0;
    for (int d = 0; d < depth; d++) {
    	for (int h = 0; h < height; h++) {
			for (int w = 0; w < width; w++) {
				
				//Do Previous Depth
				if (d != 0) {
					if (op == 27) {
						// DO Everything but center prev here
						if (w != 0 && h != 0) {
							neighborArray[neighborCounter] = nodeCounter - numPerPlane - width - 1;
							neighborCounter += 1;
						}
						if (h != 0) {
							neighborArray[neighborCounter] = nodeCounter - numPerPlane - width;
							neighborCounter +=1;
						}
						if (w != (width - 1) && h != 0) {
							neighborArray[neighborCounter] = nodeCounter - numPerPlane - width + 1;
							neighborCounter += 1;
						}
						if (w!= 0) {
							neighborArray[neighborCounter] = nodeCounter - numPerPlane - 1;
							neighborCounter += 1;
						}
					}
					neighborArray[neighborCounter] = nodeCounter - numPerPlane;
					neighborCounter += 1;
					if (op == 27) {
						if (w!= width - 1) {
							neighborArray[neighborCounter] = nodeCounter - numPerPlane + 1;
							neighborCounter += 1;
						}
						if (w!= 0 && h != (height - 1)) {
							neighborArray[neighborCounter] = nodeCounter - numPerPlane + width - 1;
							neighborCounter += 1;
						}
						if (h != (height - 1)) {
							neighborArray[neighborCounter] = nodeCounter - numPerPlane + width;
							neighborCounter +=1;
						}
						if (w != (width - 1) && h != (height - 1)) {
							neighborArray[neighborCounter] = nodeCounter - numPerPlane + width + 1;
							neighborCounter += 1;
						}
					}
				}

				//Do Current Depth
				if (h != 0 && w != 0 && op == 27) {
					neighborArray[neighborCounter] = nodeCounter - width - 1;
					neighborCounter += 1;
				}
				if (h != 0) {
					neighborArray[neighborCounter] = nodeCounter - width;
					neighborCounter += 1;
				}
				if (h != 0 && w != (width - 1) && op == 27) {
					neighborArray[neighborCounter] = nodeCounter - width + 1;
					neighborCounter += 1;
				}
				if (w != 0) {
					neighborArray[neighborCounter] = nodeCounter - 1;
					neighborCounter += 1;
				}
				if (w != width - 1) {
					neighborArray[neighborCounter] = nodeCounter + 1;
					neighborCounter += 1;
				}
				if (h != (height - 1) && w != 0 && op == 27) {
					neighborArray[neighborCounter] = nodeCounter + width - 1;
					neighborCounter += 1;
				}
				if (h != height - 1) {
					neighborArray[neighborCounter] = nodeCounter + width;
					neighborCounter += 1;
				}
				if (h != (height - 1) && w != (width - 1) && op == 27) {
					neighborArray[neighborCounter] = nodeCounter + width + 1;
					neighborCounter += 1;
				}

				//Do Next Depth
				if (d != depth - 1) {
					if (op == 27) {
						// DO Everything but center prev here
						if (w != 0 && h != 0) {
							neighborArray[neighborCounter] = nodeCounter + numPerPlane - width - 1;
							neighborCounter += 1;
						}
						if (h != 0) {
							neighborArray[neighborCounter] = nodeCounter + numPerPlane - width;
							neighborCounter +=1;
						}
						if (w != (width - 1) && h != 0) {
							neighborArray[neighborCounter] = nodeCounter + numPerPlane - width + 1;
							neighborCounter += 1;
						}
						if (w!= 0) {
							neighborArray[neighborCounter] = nodeCounter + numPerPlane - 1;
							neighborCounter += 1;
						}
					}

					neighborArray[neighborCounter] = nodeCounter + numPerPlane;
					neighborCounter +=1;

					if (op == 27) {
						if (w!= width - 1) {
							neighborArray[neighborCounter] = nodeCounter + numPerPlane + 1;
							neighborCounter += 1;
						}
						if (w!= 0 && h != (height - 1)) {
							neighborArray[neighborCounter] = nodeCounter + numPerPlane + width - 1;
							neighborCounter += 1;
						}
						if (h != (height - 1)) {
							neighborArray[neighborCounter] = nodeCounter + numPerPlane + width;
							neighborCounter +=1;
						}
						if (w != (width - 1) && h != (height - 1)) {
							neighborArray[neighborCounter] = nodeCounter + numPerPlane + width + 1;
							neighborCounter += 1;
						}
					}
				}
				nodeArray[nodeCounter + 1] = neighborCounter;
				nodeCounter += 1;
			}
 		}
    }

    //Reallocate the second Array to correct size.
    neighborArray = realloc(neighborArray, neighborCounter*sizeof(int));

    //Print created Arrays.
    for (int j = 0; j <= numNodes; j++) {
    	int temp = nodeArray[j] + 1;
        printf("%d ", temp);
    }
    printf("\n");
    for (int k = 0; k < neighborCounter; k++) {
    	int temp2 = neighborArray[k] + 1;
    	printf("%d ", temp2);
	}
    printf("\n");


    //Free both Arrays once printed.
    free(nodeArray);
    free(neighborArray);
    return 0;
}

