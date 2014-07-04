#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "util.h"

//#define STOPCOUNT 500000000
#define STENCIL 4
#define POPSIZE 20
//#define MAXGENS 50000
//#define PMUTATION  .5
#define NUM_GENES 100
//#define INPUT_FILE "population10x10.dat"
//#define ADJ_FILE "input10x10_5.adj"



char* INPUT_FILE;
char* ADJ_FILE;
int MAXGENS;
int STOPCOUNT;
//int STENCIL;
double PMUTATION;
//int POPSIZE;
int swapPercent;
int swapLength;


struct averageStruct {
    int average;
    int value;
};

struct individual {
    double fitness;
    double rfitness;
    double cfitness;
    int ordering[NUM_GENES];
};

struct individual currentPop[POPSIZE];
struct individual nextPop[2 * POPSIZE];

int adjArray1[NUM_GENES + 1];
int adjArray2[NUM_GENES * STENCIL];

int main(int arc, char *argv[]);
struct individual cross (int firstParent, int secondParent);
void evaluateOrdering(struct individual indivs[], int size);
void init();
void keepBest();
void mutatePop(struct individual mutated[]);
int pickParent(struct individual possibleParents[]);
void printPop();
int costComp(const void * a, const void * b);
int averageStructComp(const void * a, const void * b);
int isPermutation(int possiblePermutation[]);

double totalFitness;
int currentGen = 0;

double costStopCounter = 0;
double costStopScore = 0;
int breakThreshold;

char* mutType;
char* crossType;




/* The main function for a program to generate Matrix Orderings which
   require the least possible fill during factorization and other
   matrix operations. Will implement a genetic algorithm to search
   for the optimal matrix ordering. Does not take any commandline
   arguments as of now.*/
int main (int argc, char *argv[]) {

    if (argc != 11) {
        fprintf(stderr,"Error: Wrong number of Arugments. Should be the following:\nPopulation File\nAdjaceny File \nMax # of Generations\nChance of mutating individual\nPercentage of the indiv to be mutated\nLength of genes to be mutated at once\n# of generations to stop after no improvement\nThreshold to stop the generation\nType of Mutation\nType of Crossover\n");
        exit(5);
    } else {
        INPUT_FILE = argv[1];
        ADJ_FILE = argv[2];
        MAXGENS = atoi(argv[3]);
        PMUTATION = atof(argv[4]);
        swapPercent = atoi(argv[5]);
        swapLength = atoi(argv[6]);
        STOPCOUNT = atoi(argv[7]);
        breakThreshold = atoi (argv[8]);
        mutType = argv[9];
        crossType = argv[10];
    }
    
    setbuf(stdout, NULL);
    init();

    printf("Initial Population \n");
    printPop();
    struct individual child;
    int parent1;
    int parent2;
    while (currentGen < MAXGENS) {

        evaluateOrdering(currentPop, POPSIZE);
        memcpy(nextPop, currentPop, POPSIZE * sizeof(struct individual));
        for (int j = POPSIZE; j < ((2 * POPSIZE)); j++) {

            parent1 = pickParent(currentPop);
            parent2 = pickParent(currentPop);
            child = cross(parent1, parent2);
            nextPop[j] = child;
        }
        mutatePop(nextPop);
        evaluateOrdering(nextPop, 2*POPSIZE);
        qsort(nextPop, POPSIZE * 2, sizeof(struct individual), costComp);

        if (nextPop[0].fitness > costStopScore) {
            costStopScore = nextPop[0].fitness;
            costStopCounter = 0;
        }
        costStopCounter += 1;

        if ((1/nextPop[0].fitness) < breakThreshold) {
            break;
        }

        memcpy(currentPop, nextPop, POPSIZE * sizeof(struct individual));
        currentGen += 1;
        if (costStopCounter == STOPCOUNT) {
            break;
        }
    } 
    printf("Final Population with %d generations.\n", currentGen);
    printPop();


    if (isPermutation(currentPop[0].ordering)) {
        printf("The best individual in this population is a real permutation\n");
    }


}
/*Comparison function used to order individuals by fitness*/
int costComp(const void * a, const void * b) {
    struct individual *indiv1 = (struct individual *) a;
    struct individual *indiv2 = (struct individual *) b;
    if (indiv2->fitness == indiv1->fitness) {
        return 0;
    } else if (indiv2->fitness > indiv1->fitness) {
        return 1;
    } else {
        return -1;
    }
}

/*Comparison function used for the average crossover algorithm*/
int averageStructComp(const void * a, const void * b) {
    struct averageStruct *pair1 = (struct averageStruct *) a;
    struct averageStruct *pair2 = (struct averageStruct *) b;
    if (pair1->average == pair2->average) {
        return 0;
    } else if (pair1->average > pair2->average) {
        return 1;
    } else {
        return -1;
    }
}

/*Prints the genes of the current population. */
void printPop() {
    double average = 0;
    for (int i = 0; i < POPSIZE; i++) {
        for (int j = 0; j < NUM_GENES; j++) {
            printf("%d ", currentPop[i].ordering[j]);
        }
        double thisFit;
        if (currentGen != 0) {
            thisFit = (1 / currentPop[i].fitness);
        } else {
            thisFit = ((GetCost(NUM_GENES, NUM_GENES * STENCIL, adjArray1, adjArray2, currentPop[i].ordering)));
        }
        average += thisFit;
        printf("Fitness: %f\n", thisFit);
    }
    average = average / POPSIZE;
    printf("Average is %f", average);
    printf("\n");
}

/* Implements the crossing of two selected matrix orderings in order
   to complete the descendant stage of the genetic algorithm. */
struct individual cross (int firstParent, int secondParent) {
    struct individual child;
    if (!strcmp("order", crossType)) {
        int crossPoint1 = rand() % NUM_GENES;
        int crossPoint2 = 0;
        while (crossPoint2 < crossPoint1) {
            crossPoint2 = rand() % NUM_GENES;
        }
        int takenGenes[NUM_GENES + 1];
        for (int zeroing = 0; zeroing < NUM_GENES + 1; zeroing++){
            takenGenes[zeroing] = 0; 
        }
        for (int genes = crossPoint1; genes < crossPoint2; genes++) {
            child.ordering[genes] = currentPop[firstParent].ordering[genes];
            takenGenes[child.ordering[genes]] = 1;
        }

        int counter = 0;
        for (int secondPass = 0; secondPass < NUM_GENES; secondPass++) {
            if (!takenGenes[currentPop[secondParent].ordering[secondPass]]) {
                child.ordering[counter] = currentPop[secondParent].ordering[secondPass];
                takenGenes[child.ordering[counter]] = 1;
                counter++;
                if (counter == crossPoint1) {break;}
            }
        }
        counter = crossPoint2;
        for (int thirdPass = 0; thirdPass < NUM_GENES; thirdPass++) {
            if (!takenGenes[currentPop[secondParent].ordering[thirdPass]]) {
                child.ordering[counter] = currentPop[secondParent].ordering[thirdPass];
                counter++;
            }
        }
        if (!isPermutation(child.ordering)) {
            printf("Child is illegal permutation.\n");
            printf("crosspoint1 is %d and crosspoint2 is %d\n", crossPoint1, crossPoint2);
            printf("Dumping Child ordering\n");
            for (int z = 0; z < NUM_GENES; z++) {
                printf("%d %d \n", z,child.ordering[z] );
            }
            printf("\n");
            exit(1);
        }
    } else if (!strcmp("average", crossType)) {
        struct averageStruct newOrder[NUM_GENES +1];
        for (int init = 0; init < NUM_GENES + 1; init++) {
            newOrder[init].value = init;
        }
        for (int location = 0; location < NUM_GENES; location++) {
            newOrder[currentPop[firstParent].ordering[location]].average = location;
        }
        for (int location2 = 0; location2 < NUM_GENES; location2++) {
            newOrder[currentPop[secondParent].ordering[location2]].average += location2;
        }
        for (int division = 0; division < NUM_GENES; division++) {
            newOrder[division].average /= 2;
        }
        qsort(newOrder, NUM_GENES + 1, sizeof(struct averageStruct), averageStructComp);
        int childCounter = 0;
        for (int i = 0; i < (NUM_GENES + 1); i++) {
            if (newOrder[i].value != 0) {
                child.ordering[childCounter] = newOrder[i].value;
                childCounter += 1;
            }
        }
        if (!isPermutation(child.ordering)) {
            printf("Child is illegal permutation.\n");
            exit(1);
        }
    } else {
        fprintf(stderr,"Not a legal type of crossover.");
        exit(1);
    }
    return child;
}

/* Picks an individual in the current generation to serve as a parent
   for the next generation. Will use the fraction of total population
   fitness in order to select (on average) those parents which are
   more fit to produce the next generation. Algorithm for this borrowed
   from FSU tutorial, perhaps can be improved. */
int pickParent(struct individual possibleParents[]) {
    totalFitness=0;
    //Find the total fitness of the population.
    for (int totalSum = 0; totalSum < POPSIZE; totalSum++) {
        totalFitness += currentPop[totalSum].fitness;
    }
    for (int fitCounter = 0; fitCounter < POPSIZE; fitCounter++) {
        currentPop[fitCounter].rfitness = currentPop[fitCounter].fitness / totalFitness;
    }
    //Find the cumulative fitness of the individuals.
    currentPop[0].cfitness = currentPop[0].rfitness;
    for (int cfitCounter = 1; cfitCounter < POPSIZE; cfitCounter ++) {
        currentPop[cfitCounter].cfitness = currentPop[cfitCounter - 1].cfitness + currentPop[cfitCounter].rfitness;
    }
    
    //Select and return 1 parent, based off of cfitness.
    double random;
    random =((double)rand()/(double)RAND_MAX);
    if (random < currentPop[0].cfitness) {
        return 0;
    } else {
        for (int picker = 0; picker < POPSIZE; picker++) {
            if (currentPop[picker].cfitness <= random && random < currentPop[picker + 1].cfitness) {
                return picker + 1;
            }
        }
    }
}

/* Implements the evaluation function for a matrix ordering. Will store
   the value of this ordering into the structure that contains it. */
void evaluateOrdering(struct individual indivs[], int size) {
    for (int in = 0; in < size; in++) {
        indivs[in].fitness = 1 / (GetCost(NUM_GENES, NUM_GENES * STENCIL, adjArray1, adjArray2 ,indivs[in].ordering));
    }
}

/* Initializze the original generation to begin the simulation.
   probably involves reading in my other program and using that
   adjacency list to somehow generate a population. Also seeds
   the random number generator used for ohter functions of this
   program. Also reads in the adjacency input file.*/
void init() {
    srand(time(NULL));
    FILE *inFile;
    FILE *adjFile;
    int number;
    inFile = fopen(INPUT_FILE, "r");
    adjFile = fopen(ADJ_FILE, "r");
    if (inFile == NULL) {
        fprintf(stderr, "Unable to open file\n");
        exit(1);
    }
    for (int indiv = 0; indiv < POPSIZE; indiv ++) {
        for (int gene = 0; gene < NUM_GENES; gene++) {
            fscanf(inFile, "%d", &number);
            if (feof(inFile)) {
                fprintf(stderr, "Reached EOF Early, check parameters\n");
                exit(1);
            }
            currentPop[indiv].ordering[gene] = number;
        }
    }
    fclose(inFile);
    if (adjFile == NULL) {
        fprintf(stderr, "Unable to open adjacency file\n");
        exit(1);
    }
    for (int i=0; i < NUM_GENES + 1; i++) {
        fscanf(adjFile, "%d", &number);
        adjArray1[i] = number;
    }
    for (int k = 0; k < (NUM_GENES * STENCIL); k++) {
        fscanf(adjFile, "%d", &number);
        adjArray2[k] = number;
    }
}

/* Mutates the entire current population. */
void mutatePop(struct individual mutated[]) {
    for (int indiv = 0; indiv < POPSIZE; indiv++) { 
        int randomMutate;
        int mutateBarrier;
        mutateBarrier = PMUTATION * RAND_MAX;
        randomMutate=rand();
        if (randomMutate < mutateBarrier) {
            int numSwaps = swapPercent * POPSIZE * .01;
            for (int x = 0; x < numSwaps; x++) {
                if (!strcmp("swap", mutType)) {
                    int temp[NUM_GENES]; //remove * to revert also is temp1
                    int gene1 = rand() % (NUM_GENES - (swapLength + 1));
                    int gene2 = gene1;
                    while (((gene1 - swapLength) <= gene2) && ((gene1 + swapLength) >= gene2)) {
                        gene2 = rand() % (NUM_GENES - (swapLength + 1));
                    }
                    memcpy(&temp[0], &mutated[indiv].ordering[gene1], (swapLength * sizeof(int)));//temp = mutated[indiv].ordering[gene1];
                    memcpy(&mutated[indiv].ordering[gene1], &mutated[indiv].ordering[gene2], (swapLength * sizeof(int))); //mutated[indiv].ordering[gene1] = mutated[indiv].ordering[gene2];
                    memcpy(&mutated[indiv].ordering[gene2], &temp[0], (swapLength * sizeof(int)));//mutated[indiv].ordering[gene2] = temp;
                    if (! isPermutation(mutated[indiv].ordering)) {
                        printf(" Gene 1 is %d and Gene 2 is %d\n", gene1, gene2);
                        printf("Dumping temp");
                        for (int z = 0; z < NUM_GENES; z++) {
                            printf("%d ", temp[z] );
                        }
                        printf("\n");
                        exit(1);
                    }
                } else if (!strcmp("shift", mutType)) {
                    int temp[NUM_GENES];
                    int shiftPoint = rand() % (NUM_GENES - (swapLength + 1));
                    memcpy(&temp[0], &mutated[indiv].ordering[shiftPoint], (swapLength * sizeof(int)));
                    memcpy(&mutated[indiv].ordering[shiftPoint], &mutated[indiv].ordering[shiftPoint + swapLength], ((NUM_GENES - (shiftPoint + swapLength)) * sizeof(int)));
                    memcpy(&mutated[indiv].ordering[NUM_GENES - swapLength], &temp[0], swapLength * sizeof(int));
                    if (! isPermutation(mutated[indiv].ordering)) {
                        printf("shiftPoint is %d\n", shiftPoint);
                        printf("Dumping temp");
                        for (int z = 0; z < NUM_GENES; z++) {
                            printf("%d ", temp[z] );
                        }
                        printf("\n");
                        exit(1);
                    }


                } else if (!strcmp("invert", mutType)) {
                    int temp[NUM_GENES];
                    int flipPoint = rand() % (NUM_GENES - (swapLength + 1));
                    memcpy(&temp[0], &mutated[indiv].ordering[flipPoint], (swapLength * sizeof(int)));
                    for (int swapper = swapLength - 1; swapper >= 0; swapper--) {
                        int counter = 0;
                        mutated[indiv].ordering[flipPoint + counter] = temp[swapper];
                        counter++;
                    }
                    if (! isPermutation(mutated[indiv].ordering)) {
                        printf("flipPoint is %d\n", flipPoint);
                        printf("Dumping temp");
                        for (int z = 0; z < NUM_GENES; z++) {
                            printf("%d ", temp[z] );
                        }
                        printf("\n");
                        exit(1);
                    }

                } else {
                    fprintf(stderr, "Invalid mutation type. Mutation must be either swap, shift, or invert.\n");
                    exit(1);
                }
            }
        }
    }
}

/*Utility function to determine if argument is a legal permutation*/
int isPermutation(int possiblePermutation[]) {
    int boolArray[NUM_GENES];
    int testNum;
    for (int i = 0; i < NUM_GENES; i++) {
        boolArray[i] = 0;
    }
    for (int j = 0; j < NUM_GENES; j++) {
        testNum = possiblePermutation[j];
        if (boolArray[(testNum - 1)] == 0) {
            boolArray[(testNum - 1)] = 1;
        } else {
            printf("This is not a legal permutation. Printing and exiting.\n");
            printf("Printing Possible Permutation ");
            for (int k = 0; k < NUM_GENES; k++) {
                printf("%d ", possiblePermutation[k]);
            }
            printf("Printing double number: %d \n", testNum);
            return 0;
        }
    }
    return 1;
}


