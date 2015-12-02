#include <stdlib.h>
#include <stdio.h>

/*Torus representation for adjacency arrays.
 Author: Kevin Carlson */

void fill_adjacency_arrays(unsigned int height, unsigned int width,
                             unsigned int op, unsigned int* nodeArray,
                             unsigned int* neighborArray);

int comp (const void * elem1, const void * elem2) {
    int f = *((int*)elem1);
    int s = *((int*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Invalid arguments, please give length, width, and operator.");
    }

    //Declare the ints needed for generating the arrays.
    unsigned int height;
    unsigned int width;
    unsigned int op;
    unsigned int* nodeArray;
    unsigned int* neighborArray;

    //Read the input to the program.
    height = atoi(argv[1]);
    width = atoi(argv[2]);
    op = atoi(argv[3]);

    //Calcualte the size of arrays and malloc.
    int numNodes;
    int numAdjacents;
    int singleAdjacents;
    numNodes = height * width;
    singleAdjacents = op - 1;
    numAdjacents = numNodes * singleAdjacents;
    nodeArray = (unsigned int*) malloc((numNodes+1)*sizeof(unsigned int));
    neighborArray = (unsigned int*) malloc(numAdjacents*sizeof(unsigned int));

    //Call the function to create the arrays.
    fill_adjacency_arrays(height, width, op, nodeArray, neighborArray);

    //Print created Arrays.
    for (int j = 0; j <= height * width; j++) {
        printf("%d ", nodeArray[j]);
    }
    printf("\n");
    for (int k = 0; k < numAdjacents; k++) {
    	printf("%d ", neighborArray[k]);
    }
    printf("\n");


    //Free both Arrays once printed.
    free(nodeArray);
    free(neighborArray);
    return 0;
}

void fill_adjacency_arrays(unsigned int height, unsigned int width,
                             unsigned int op, unsigned int* nodeArray,
                             unsigned int* neighborArray) {
	int bottomRowStart = 0;
	int bottomRowEnd = width - 1;
	int topRowStart = width * (height - 1);
	int topRowEnd = (height * width) - 1;
	int leftEdge = 0;
	int rightEdge = width - 1;
	int neighborArrayIndex = 0;
    int tempArray[100];
    int tempArrayIndex = 0;
    for (int i = 0; i < (height * width); i ++) {
        nodeArray[i] = i * (op - 1);
        
        //Neighbor Below Left
        if (op >= 7) {
        	if (i == bottomRowStart) {
        		//Do bottomCorner stuff
        		tempArray[tempArrayIndex] = (height*width) - 1;
//        		leftEdge += width;
        	} else if ((bottomRowStart < i) && (bottomRowEnd >= i)) {
        		//Do BottomRow Stuff
        		tempArray[tempArrayIndex] = i + (width * (height - 1) - 1);
        	} else if (i == leftEdge) {
        		//Do LeftEdge Stuff
        		tempArray[tempArrayIndex] = i - 1;
//        		leftEdge += width;
        	} else {
        		//Do normal stuff
        		tempArray[tempArrayIndex] = i - width - 1;
        	}
        	tempArrayIndex += 1;

        }

        //Neighbor Below Center
        if (i >= bottomRowStart && i <= bottomRowEnd) {
        	//Do Bottom row
        	tempArray[tempArrayIndex] = i + (width * (height - 1));
        } else {
        	tempArray[tempArrayIndex] = i - width;
        }
        tempArrayIndex += 1;

        //Neighbor Below Right
        if (op == 9) {
        	if (i == bottomRowEnd) {
        		tempArray[tempArrayIndex] = width * (height - 1);
        	} else if ((bottomRowStart <= i) && (bottomRowEnd > i)) {
        		tempArray[tempArrayIndex] = i + (width * (height - 1)) + 1;
        	} else if (i == rightEdge) {
        		tempArray[tempArrayIndex] = i - (2*width) + 1;
        	} else {
        		tempArray[tempArrayIndex] = i - width + 1;
        	}
        	tempArrayIndex += 1;
        }

        //Neighbor Left
        if (i == leftEdge) {
        	tempArray[tempArrayIndex] = i + width - 1;
        } else {
        	tempArray[tempArrayIndex] = i - 1;
        }
        tempArrayIndex += 1;

        //Neighbor Right
        if (i == rightEdge) {
        	tempArray[tempArrayIndex] = i - width + 1;
        } else {
        	tempArray[tempArrayIndex] = i + 1;
        }
        tempArrayIndex += 1;

        //Neighbor Above Left
        if (op == 9) {
        	if (i == topRowStart) {
        		tempArray[tempArrayIndex] = width - 1;
        	} else if ((i > topRowStart) && (i <= topRowEnd)) {
        		tempArray[tempArrayIndex] = i - (width * (height - 1)) - 1;
        	} else if (i == leftEdge) {
        		tempArray[tempArrayIndex] = i + (2*width) - 1;
        	} else {
        		tempArray[tempArrayIndex] = i + width - 1;
        	}
        	tempArrayIndex += 1;
        }

        //Neighbor Above Center
		if (i >= topRowStart && i <= topRowEnd) {
			tempArray[tempArrayIndex] = i - (width * (height - 1));
		} else {
			tempArray[tempArrayIndex] = i + width;
		}
		tempArrayIndex += 1;		

        //Neighbor above Right
        if (op >= 7) {
        	if (i == topRowEnd) {
        		tempArray[tempArrayIndex] = 0;
        	} else if (i >= topRowStart && i < topRowEnd) {
        		tempArray[tempArrayIndex] = i - (width*(height - 1)) + 1;
        	} else if (i ==rightEdge) {
        		tempArray[tempArrayIndex] = i + 1;
        	} else {
        		tempArray[tempArrayIndex] = i + width + 1;
        	}
        	tempArrayIndex += 1;
        }

        //Update left edge
        if (i == leftEdge) {
        	leftEdge += width;
        }

        if (i == rightEdge) {
        	rightEdge += width;
        }

        qsort(tempArray, (op - 1), sizeof(int), comp);
        for (int z = 0; z < op - 1; z++) {
            neighborArray[neighborArrayIndex] = tempArray[z];
            neighborArrayIndex += 1;
        }
        tempArrayIndex = 0;
    }
    nodeArray[height*width] = (((height*width)) *(op-1));
}