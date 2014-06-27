#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "cost.h"

#define STOPCOUNT 500000000
#define PARAMETER 4
#define POPSIZE 20
#define MAXGENS 100000
#define PMUTATION  .5
#define NUM_GENES 9
#define INPUT_FILE "population.dat"
#define ADJ_FILE "input.adj"

struct individual {
    double fitness;
    double rfitness;
    double cfitness;
    int ordering[NUM_GENES];
};

struct individual currentPop[POPSIZE];
struct individual nextPop[2 * POPSIZE];

int adjArray1[NUM_GENES + 1];
int adjArray2[NUM_GENES * PARAMETER];

int main();
struct individual cross (int firstParent, int secondParent);
void evaluateOrdering(struct individual indivs[]);
void init();
void keepBest();
void mutatePop(struct individual mutated[]);
int pickParent(struct individual possibleParents[]);
void printPop();
int costComp(const void * a, const void * b);
int isPermutation(int possiblePermutation[]);

double totalFitness;
int currentGen = 0;
double bestEverScore = 10000000;
double bestCost;

double costStopCounter = 0;
double costStopScore = 0;



/* The main function for a program to generate Matrix Orderings which
   require the least possible fill during factorization and other
   matrix operations. Will implement a genetic algorithm to search
   for the optimal matrix ordering. Does not take any commandline
   arguments as of now.*/
int main (void) {
    
    setbuf(stdout, NULL);
    init();

    /*Adding stuff to test total perms here.*/
    int testordering[NUM_GENES];
    int testNum;
    int bestOrdering[NUM_GENES];
    double testScore;
    double bestScore = 10000;
    /*FILE *permsFile;
    permsFile = fopen("allperms.txt", "r");
    if (permsFile == NULL) {
        fprintf(stderr, "Unable to open perms file\n");
        exit(1);
    }
    for (unsigned int i = 0; i < 362880; i++) {
        for (int j = 0; j < NUM_GENES; j++) {
            fscanf(permsFile, "%d", &testNum);
            testordering[j] = testNum;
        }
        if (! isPermutation(testordering)) {
            printf("i is %d \n", i);
        }
        if (i == 0) {
            printf(" The first ordering read in is: ");
            for (int printing= 0; printing< NUM_GENES; printing++) {
                printf("%d ", testordering[printing]);
            }
            printf("\n");
        }
        testScore = GetCost(NUM_GENES, NUM_GENES * PARAMETER, adjArray1, adjArray2, testordering);
        //printf("%f\n", testScore/*GetCost(NUM_GENES, NUM_GENES * PARAMETER, adjArray1, adjArray2, testordering));
        if (testScore < bestScore) {
            bestScore = testScore;
            memcpy(bestOrdering, testordering, NUM_GENES * sizeof(int));
        }
    }
    printf("The best ordering has fill %f\n", bestScore);
    printf("The best ordering is ");
    for (int k = 0; k < NUM_GENES; k++) {
        printf("%d ",bestOrdering[k]);
    }
    printf("\n");

    printf("HardCoding a testCase \n");
    int hardCode[9] = {9, 2, 1, 4, 5, 6, 7, 8, 9};
    printf("Fitness of HardCode: %f\n", ((GetCost(NUM_GENES, NUM_GENES * PARAMETER, adjArray1, adjArray2, hardCode))));
    printf("HardCoded Ordering is:");
    for (int k = 0; k < NUM_GENES; k++) {
        printf("%d ",hardCode[k]);
    }
    printf("\n");*/



    printf("Initial Population \n");
    printPop();
    struct individual child;
    int parent1;
    int parent2;
    while (currentGen < MAXGENS) {


        //printf("Entered the while loop before crashing\n");
        memcpy(nextPop, currentPop, POPSIZE * sizeof(struct individual));
        evaluateOrdering(nextPop);
        //printf("Copy one works\n");
        //for (int test2 = 0; test2 < POPSIZE; test2++) {
        //        printf(" %.2f", (1/nextPop[test2].fitness));
        //}
        //printf("\n");
        for (int j = POPSIZE; j < ((2 * POPSIZE)); j++) {

            parent1 = pickParent(currentPop);
            //printf("Can I pick a parent?%d\n", parent1);
            parent2 = pickParent(currentPop);
            //printf("Can I pick a second parent?%d\n", parent2);
            child = cross(parent1, parent2);
            //printf("Cross is happening!");
            nextPop[j] = child;
        }
        //printf("Children made\n");
        mutatePop(nextPop);
        //printf("Mutated\n");
        evaluateOrdering(nextPop);
        //printf("Evaled\n");
        qsort(nextPop, POPSIZE * 2, sizeof(struct individual), costComp);

        //for (int test = 0; test < 2*POPSIZE; test++) {
        //    printf(" %.2f", 1/(nextPop[test].fitness));
        //}
        //printf("\n");

        if (nextPop[0].fitness > costStopScore) {
            costStopScore = nextPop[0].fitness;
            costStopCounter = 0;
        }
        costStopCounter += 1;

        memcpy(currentPop, nextPop, POPSIZE * sizeof(struct individual));
        currentGen += 1;
        if (costStopCounter == STOPCOUNT) {
            break;
        }
    } 
    printf("Final Population with %d generations.\n", currentGen);
    printPop();

    //printf("Best Pop fitness seen this whole run is %f\n", bestEverScore);
}

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


/*Prints the genes of the current population. */
void printPop() {
    for (int i = 0; i < POPSIZE; i++) {
        for (int j = 0; j < NUM_GENES; j++) {
            printf("%d ", currentPop[i].ordering[j]);
        }
        //printf("\n");
        printf("Fitness: %f\n", ((GetCost(NUM_GENES, NUM_GENES * PARAMETER, adjArray1, adjArray2, currentPop[i].ordering))));
 //       printf("The fitness of this individual is %f\n", currentPop[i].fitness);
    }
    printf("\n");
}



/* Implements the crossing of two selected matrix orderings in order
   to complete the descendant stage of the genetic algorithm.
   TODO: Determine what arguments this function will need. */
struct individual cross (int firstParent, int secondParent) {
    int crossPoint = rand() % NUM_GENES;
    struct individual child;


    int takenGenes[NUM_GENES + 1];
    for (int zeroing = 0; zeroing < NUM_GENES + 1; zeroing++){
        takenGenes[zeroing] = 0; 
    }

    int counter = 0;
    for (int genes = 0; genes < crossPoint; genes++) {
        child.ordering[counter] = currentPop[firstParent].ordering[genes];
        takenGenes[child.ordering[counter]] = 1;
        counter++;
    }
    for (int secondPass = 0; secondPass < NUM_GENES; secondPass++) {
        if (!takenGenes[currentPop[secondParent].ordering[secondPass]]) {
            child.ordering[counter] = currentPop[secondParent].ordering[secondPass];
            counter++;
        }
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
    //Find the relative fitness of each individual.
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
    int test;
    random =((double)rand()/(double)RAND_MAX);
    if (random < currentPop[0].cfitness) {
        //printf("Selected 0 with fitness %f\n", currentPop[0].fitness);
        return 0;
    } else {
        for (int picker = 0; picker < POPSIZE; picker++) {
            if (currentPop[picker].cfitness <= random && random < currentPop[picker + 1].cfitness) {
                //printf("Selected %d with fitness %f\n", picker + 1, currentPop[picker + 1].fitness);
                return picker + 1;
            } else {return 1;}
        }
    }
}


/* Implements the evaluation function for a matrix ordering. Will store
   the value of this ordering into the structure that contains it.
   TODO: Determine what arguments this function will need. */
void evaluateOrdering(struct individual indivs[]) {
    for (int i = 0; i < POPSIZE * 2; i++) {
        indivs[i].fitness = 1 / (GetCost(NUM_GENES, NUM_GENES * PARAMETER, adjArray1, adjArray2 ,indivs[i].ordering));
        double tempcost = GetCost(NUM_GENES, NUM_GENES * PARAMETER, adjArray1, adjArray2 ,indivs[i].ordering);
        if (tempcost < bestEverScore) {
            bestEverScore = tempcost;
        }
        double costTest = 1/tempcost;
        if (costTest > bestCost) {
            bestCost = costTest;
        }
     //   printf("Fitness is %f\n",indivs[i].fitness);
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
    adjFile = fopen(ADJ_FILE, "r");
    if (adjFile == NULL) {
        fprintf(stderr, "Unable to open adjacency file\n");
        exit(1);
    }
    for (int i=0; i < NUM_GENES + 1; i++) {
        fscanf(adjFile, "%d", &number);
        adjArray1[i] = number;
    }
    for (int k = 0; k < (NUM_GENES * PARAMETER); k++) {
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
            /*DO MUTATE STUFF*/
            int temp1;
            int temp2;
            int temp3;
            int gene1 = rand() % NUM_GENES;
            int gene2 = rand() % NUM_GENES;
            int gene3 = rand() % NUM_GENES;
            int gene4 = rand() % NUM_GENES;
            int gene5 = rand() % NUM_GENES;
            int gene6 = rand() % NUM_GENES;
            temp1 = mutated[indiv].ordering[gene1];
            mutated[indiv].ordering[gene1] = mutated[indiv].ordering[gene2];
            mutated[indiv].ordering[gene2] = temp1;
            temp2 = mutated[indiv].ordering[gene3];
            mutated[indiv].ordering[gene3] = mutated[indiv].ordering[gene4];
            mutated[indiv].ordering[gene4] = temp2;
            temp3 = mutated[indiv].ordering[gene5];
            mutated[indiv].ordering[gene5] = mutated[indiv].ordering[gene6];
            mutated[indiv].ordering[gene6] = temp3;
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
            printf("\nPrinting double number: %d \n", testNum);
            return 0;
        }
    }
    return 1;
}


