#ifndef GENETICORDERING_H
#define GENETICORDERING_H

struct averageStruct {
    int average;
    int value;
};

struct individual {
    double fitness;
    double rfitness;
    double cfitness;
    int* ordering;
};

struct individual** population;

int* adjArray1;
int* adjArray2;

char* inputFile;
char* adjacencyFile;
char* mutType;
char* crossType;
char* selectionType;
char* generationType;

int runningProgram = 1;
int nnz;
int n;
int maxGens;
int stopCount;
int popSize;
int maxPopSize;
int swapLength;
int growthNumber;
int currentGen;
int breakThreshold;
int costCounter;

double pMutation;
double maxPopScale;
double swapPercent;
double totalFitness;
double costStopCounter;
double costStopScore;
double fillPercent;
double childPercent;

time_t seed;

int main (int argc, char *argv[]);

struct individual* cross (int firstParent, int secondParent);

void evaluateOrdering(struct individual* indivs[], int size);
void init(int argc, char *argv[]);
void printPop();
void mutatePop(struct individual* mutated[]);
void swapMutate(struct individual* mutated[], int numMutations, int indiv);
void shiftMutate(struct individual* mutated[], int numMutations, int indiv);
void invertMutate(struct individual* mutated[], int numMutations, int indiv);
void orderCrossover(struct individual* child, int firstParent, int secondParent);
void averageCrossover(struct individual* child, int firstParent, int secondParent);
void noneCrossover(struct individual* child, int firstParent, int secondParent);
void prefixCrossover(struct individual* child, int firstParent, int secondParent);
void prefix2Crossover(struct individual* child, int firstParent, int secondparent);
void uxCrossover(struct individual* child, int firstParent, int secondParent);
void fractionGeneration(individual** population);
void sortGeneration(individual** population);
void signalStopper(int signum);

int costComp(const void * a, const void * b);
int averageStructComp(const void * a, const void * b);
int isPermutation(int possiblePermutation[]);
int pickParent(struct individual* possibleParents[]);
int rankPicking(struct individual* possibleParents[]);
int fitPicking(struct individual* possibleParents[]);


#endif /* GENETICORDERING_H */
