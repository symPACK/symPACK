#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include "util.h"
#include <omp.h>

#define POPSIZE 20
#define NUM_GENES 100


char* INPUT_FILE;
char* ADJ_FILE;
int MAXGENS;
int STOPCOUNT;
double PMUTATION;
//int POPSIZE = 20;
double swapPercent;
int swapLength;


struct averageStruct {
    int average;
    int value;
};

struct individual {
    double fitness{-1};
    double rfitness;
    double cfitness;
    int ordering[NUM_GENES];
};

struct individual currentPop[POPSIZE];
struct individual nextPop[2*POPSIZE];
int* adjArray1;
int* adjArray2;

int nnz;
int n;

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
char* selectionType;

int costCounter = 0;




/* The main function for a program to generate Matrix Orderings which
   require the least possible fill during factorization and other
   matrix operations. Will implement a genetic algorithm to search
   for the optimal matrix ordering. Does not take any commandline
   arguments as of now.*/
int main (int argc, char *argv[]) {

    if (argc != 12) {
        fprintf(stderr,"Error: Wrong number of Arugments. Should be the following:\nPopulation File\nAdjaceny File \nMax # of Generations\nChance of mutating individual\nPercentage of the indiv to be mutated\nLength of genes to be mutated at once\n# of generations to stop after no improvement\nThreshold to stop the generation\nType of Mutation\nType of Crossover\n");
        exit(5);
    } else {
        INPUT_FILE = argv[1];
        ADJ_FILE = argv[2];
        MAXGENS = atoi(argv[3]);
        PMUTATION = atof(argv[4]);
        swapPercent = atof(argv[5]);
        swapLength = atoi(argv[6]);
        STOPCOUNT = atoi(argv[7]);
        breakThreshold = atoi(argv[8]);
        mutType = argv[9];
        crossType = argv[10];
        selectionType = argv[11];
    }

    setbuf(stdout, NULL);
    init();

    printf("Initial Population \n");
    printPop();
    evaluateOrdering(currentPop, POPSIZE);
    while (currentGen < MAXGENS) {

        memcpy(nextPop, currentPop, POPSIZE * sizeof(struct individual));

        //#pragma omp parallel for
        for (int j = POPSIZE; j < ((2 * POPSIZE)); j++) {
            struct individual child;
            int parent1;
            int parent2 = -1;
            parent1 = pickParent(currentPop);
            while (parent2 < 0 || parent2 == parent1) {
                parent2 = pickParent(currentPop);
            }
            //printf("Parent 1 is %d and Parent 2 is %d\n", parent1, parent2);
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

        if ((1/nextPop[POPSIZE - 1].fitness) < breakThreshold) {
            break;
        }

        memcpy(currentPop, nextPop, POPSIZE * sizeof(struct individual));
        currentGen += 1;
        if (costStopCounter == STOPCOUNT) {
            break;
        }
        //printf("Best Indiv is fitness %f\n", 1/currentPop[0].fitness);
        //printf("Current Generation is %d\n", currentGen);
    } 
    printf("Final Population with %d generations.\n", currentGen);
    printPop();


    if (isPermutation(currentPop[0].ordering)) {
        printf("The best individual in this population is a real permutation.\n");
    }
    printf("Cost was called %d times this run.\n", costCounter);
    printf("This run took %d generations.\n", currentGen);
    free(adjArray1);
    free(adjArray2);

    for (int tired = 1; tired < 1; tired++) {
        printf("Does this run once?");
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
        for (int j = 0; j < n; j++) {
            printf("%d ", currentPop[i].ordering[j]);
        }
        double thisFit;
        if (currentGen != 0) {
            thisFit = (1 / currentPop[i].fitness);
        } else {
            thisFit = ((GetCost(n, nnz, adjArray1, adjArray2, currentPop[i].ordering)));
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
    //printf("Using Parent1: %d and Parent2: %d\n", firstParent, secondParent);
    if (!strcmp("order", crossType)) {
        int crossOver1 = (rand() % (n-2));                   //CrossOver1 should be between 0 and n-2
        int crossOver2 = -1;
        while (crossOver2 <= (crossOver1) + 1)   {
            crossOver2 = rand() % n;                         //CrossOver2 should be > crossOver1 and < n
        }
        //printf("Crossover1 is %d and crossover2 is %d\n", crossOver1, crossOver2);

        int takenGenes[n + 1];
        for (int zeroing = 0; zeroing < n + 1; zeroing++){
            takenGenes[zeroing] = 0; 
        }

        for (int copySeg = crossOver1; copySeg < crossOver2; copySeg++) {
            child.ordering[copySeg] = currentPop[firstParent].ordering[copySeg];
            takenGenes[currentPop[firstParent].ordering[copySeg]] = 1;
        }

        int secondCrawler = 0;
        for (int i = 0; i < crossOver1; i++) {
            while (takenGenes[currentPop[secondParent].ordering[secondCrawler]] == 1) {
                secondCrawler += 1;
            }
            child.ordering[i] = currentPop[secondParent].ordering[secondCrawler];
            takenGenes[currentPop[secondParent].ordering[secondCrawler]] = 1;
        }

        for (int j = crossOver2; j < n; j++) {
            while (takenGenes[currentPop[secondParent].ordering[secondCrawler]] == 1) {
                secondCrawler += 1;
            }
            child.ordering[j] = currentPop[secondParent].ordering[secondCrawler];
            takenGenes[currentPop[secondParent].ordering[secondCrawler]] = 1;
        }

        if (!isPermutation(child.ordering)) {
            printf("Child is illegal permutation.\n");
            printf("Dumping orderings (child p1 p2)\n");
            for (int z = 0; z < n; z++) {
                printf("%d %d %d %d\n", z, child.ordering[z], currentPop[firstParent].ordering[z], currentPop[secondParent].ordering[z]);
            }
            printf("\n");
            exit(1);
        }
    } else if (!strcmp("average", crossType)) {
        struct averageStruct newOrder[n +1];
        for (int starter = 0; starter < n + 1; starter++) {
            newOrder[starter].value = starter;
        }
        for (int location = 0; location < n; location++) {
            newOrder[currentPop[firstParent].ordering[location]].average = location;
        }
        for (int location2 = 0; location2 < n; location2++) {
            newOrder[currentPop[secondParent].ordering[location2]].average += location2;
        }
        for (int division = 0; division < n; division++) {
            newOrder[division].average /= 2;
        }
        qsort(newOrder, n + 1, sizeof(struct averageStruct), averageStructComp);
        int childCounter = 0;
        for (int i = 0; i < (n + 1); i++) {
            if (newOrder[i].value != 0) {
                child.ordering[childCounter] = newOrder[i].value;
                childCounter += 1;
            }
        }
        if (!isPermutation(child.ordering)) {
            printf("Child is illegal permutation.\n");
            exit(1);
        }
    } else if (!strcmp("none", crossType)) {
        for (int sillyCopy = 0; sillyCopy < n; sillyCopy++) {
            child.ordering[sillyCopy] = currentPop[firstParent].ordering[sillyCopy];
        } 
        if (!isPermutation(child.ordering)) {
            printf("Child is illegal permutation.\n");
            exit(1);
        }
    } else {
        fprintf(stderr,"Not a legal type of crossover. Should be Average, order, or none.");
        exit(1);
    }

    int parentATesting = 0;
    int parentBTesting = 0;
    for (int testingVar = 0; testingVar < n; testingVar++) {
        if (child.ordering[testingVar] == currentPop[firstParent].ordering[testingVar]) {
            parentATesting += 1;
        }
        if (child.ordering[testingVar] == currentPop[secondParent].ordering[testingVar]) {
            parentBTesting += 1;
        }
    }
    if (parentATesting == n) {
        int newParent1 = pickParent(currentPop);
        int newParent2 = pickParent(currentPop);
        return cross(newParent1, newParent2);
    }
    if (parentBTesting == n) {
        int newParent1 = pickParent(currentPop);
        int newParent2 = pickParent(currentPop);
        return cross(newParent1, newParent2);
    }
    printf("Child generated with cost %f\n",((GetCost(n, nnz, adjArray1, adjArray2, child.ordering))));


    return child;
}

/* Picks an individual in the current generation to serve as a parent
   for the next generation. Will use the fraction of total population
   fitness in order to select (on average) those parents which are
   more fit to produce the next generation. Algorithm for this borrowed
   from FSU tutorial, perhaps can be improved. */
int pickParent(struct individual possibleParents[]) {
    if (!strcmp("rw", selectionType)) {
        totalFitness=0;
        //Find the total fitness of the population.
        for (int totalSum = 0; totalSum < POPSIZE; totalSum++) {
            totalFitness += currentPop[totalSum].fitness;
        }
        //Find the relative fitness of the population.
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
    } else if (!strcmp("ro", selectionType)) {
        totalFitness=0;
        //Find the total fitness of the population.
        for (int totalSum = 0; totalSum < POPSIZE; totalSum++) {
            totalFitness += totalSum + 1;
            //printf("Total Fitness is %f\n", totalFitness);
        }
        //Find the relative fitness of the population.
        int ranking = 1;
        for (int fitCounter = POPSIZE - 1; fitCounter >= 0; fitCounter--) {
            currentPop[fitCounter].rfitness = ranking / totalFitness;
            ranking++;
            //printf("Rfitnesses are %f\n", currentPop[fitCounter].rfitness);
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
    } else {
        fprintf(stderr,"Not a legal type of selection. Should be ro for Rank Ordering or RW for Roulette Wheel");
        exit(1);
    }
}

/* Implements the evaluation function for a matrix ordering. Will store
   the value of this ordering into the structure that contains it. */
void evaluateOrdering(struct individual indivs[], int size) {
    //#pragma omp parallel for
    for (int in = 0; in < size; in++) {
        if (indivs[in].fitness == -1) {
            indivs[in].fitness = 1.0 / (GetCost(n, nnz, adjArray1, adjArray2 ,indivs[in].ordering));
            costCounter += 1;
        }
    }
}

/* Initializze the original generation to begin the simulation.
   probably involves reading in my other program and using that
   adjacency list to somehow generate a population. Also seeds
   the random number generator used for ohter functions of this
   program. Also reads in the adjacency input file.*/
void init() {
    srand(time(NULL));
    ReadAdjacency(ADJ_FILE, &adjArray1, &adjArray2, &n, &nnz);    
    FILE *inFile;
    int number;
    inFile = fopen(INPUT_FILE, "r");
    if (inFile == NULL) {
        fprintf(stderr, "Unable to open file\n");
        exit(1);
    }
    for (int indiv = 0; indiv < POPSIZE; indiv ++) {
        for (int gene = 0; gene < n; gene++) {
            fscanf(inFile, "%d", &number);
            if (feof(inFile)) {
                fprintf(stderr, "Reached EOF Early, check parameters\n");
                exit(1);
            }
            currentPop[indiv].ordering[gene] = number;
        }
    }
    fclose(inFile);
}

/* Mutates the entire current population. */
void mutatePop(struct individual mutated[]) {
    //#pragma omp parallel for
    for (int indiv = POPSIZE; indiv < POPSIZE * 2; indiv++) {
        double randomMutate;
        double mutateBarrier;
        mutateBarrier = PMUTATION * RAND_MAX;
        randomMutate=rand();
        if (randomMutate < mutateBarrier) {
            mutated[indiv].fitness = - 1;
            int numSwaps = swapPercent * NUM_GENES * .01;//playing with this variable
            // printf("numSwaps is: %d\n", numSwaps);
            for (int x = 0; x < numSwaps; x++) {
                if (!strcmp("swap", mutType)) {
                    int temp[n]; 
                    int gene1 = rand() % (n - (swapLength + 1));
                    int gene2 = gene1;
                    while (((gene1 - swapLength) <= gene2) && ((gene1 + swapLength) >= gene2)) {
                        gene2 = rand() % (n - (swapLength + 1));
                    }
                    memcpy(&temp[0], &mutated[indiv].ordering[gene1], (swapLength * sizeof(int)));
                    memcpy(&mutated[indiv].ordering[gene1], &mutated[indiv].ordering[gene2], (swapLength * sizeof(int))); 
                    memcpy(&mutated[indiv].ordering[gene2], &temp[0], (swapLength * sizeof(int)));
                    if (! isPermutation(mutated[indiv].ordering)) {
                        printf(" Gene 1 is %d and Gene 2 is %d\n", gene1, gene2);
                        printf("Dumping temp");
                        for (int z = 0; z < n; z++) {
                            printf("%d ", temp[z] );
                        }
                        printf("\n");
                        exit(1);
                    }
                } else if (!strcmp("shift", mutType)) {
                    int temp[n];
                    int shiftPoint = rand() % (n - (swapLength + 1));
                    memcpy(&temp[0], &mutated[indiv].ordering[shiftPoint], (swapLength * sizeof(int)));
                    memcpy(&mutated[indiv].ordering[shiftPoint], &mutated[indiv].ordering[shiftPoint + swapLength], ((n - (shiftPoint + swapLength)) * sizeof(int)));
                    memcpy(&mutated[indiv].ordering[n - swapLength], &temp[0], swapLength * sizeof(int));
                    if (! isPermutation(mutated[indiv].ordering)) {
                        printf("shiftPoint is %d\n", shiftPoint);
                        printf("Dumping temp");
                        for (int z = 0; z < n; z++) {
                            printf("%d ", temp[z] );
                        }
                        printf("\n");
                        exit(1);
                    }


                } else if (!strcmp("invert", mutType)) {
                    int temp[n];
                    int flipPoint = rand() % (n - (swapLength + 1));
                    memcpy(&temp[0], &mutated[indiv].ordering[flipPoint], (swapLength * sizeof(int)));
                    for (int swapper = swapLength - 1; swapper >= 0; swapper--) {
                        int counter = 0;
                        mutated[indiv].ordering[flipPoint + counter] = temp[swapper];
                        counter++;
                    }
                    if (! isPermutation(mutated[indiv].ordering)) {
                        printf("flipPoint is %d\n", flipPoint);
                        printf("Dumping temp");
                        for (int z = 0; z < n; z++) {
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
    int boolArray[n];
    int testNum;
    for (int i = 0; i < n; i++) {
        boolArray[i] = 0;
    }
    for (int j = 0; j < n; j++) {
        testNum = possiblePermutation[j];
        if (boolArray[(testNum - 1)] == 0) {
            boolArray[(testNum - 1)] = 1;
        } else {
            printf("This is not a legal permutation. Printing and exiting.\n");
            printf("Printing Possible Permutation ");
            for (int k = 0; k < n; k++) {
                printf("%d ", possiblePermutation[k]);
            }
            printf("Printing double number: %d \n", testNum);
            return 0;
        }
    }
    return 1;
}

