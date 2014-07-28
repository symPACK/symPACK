#ifndef COST_HEADER
#define COST_HEADER


#ifdef __cplusplus
#include <vector>
#include "ETree.hpp"

using namespace std;

extern "C" {
#endif
extern double GetCost(int n, int nnz, int * xadj, int * adj,int * perm);
extern double GetCostPerCol(int n, int nnz, int * xadj, int * adj,int * perm, int * costc);

extern void GetPrefixSum(int n, int * arr, int * arrout);

extern int ReadAdjacency(const char * pfilename, int ** pxadj, int ** padj, int * n , int * nnz);

extern void freeIntPointer(int* pointer);

extern void freeCharPointer(char* pointer);


#ifdef __cplusplus
}


void ExpandSymmetric(int size,const int * colptr,const int * rowind, vector<int> & expColptr, vector<int> & expRowind);
void GetPermutedGraph(int n, int nnz, int * xadj, int * adj, int * perm,int * invp, int * newxadj, int * newadj);
void GetPermutedGraph(int n, int nnz, int * xadj, int * adj, int * perm, int * newxadj, int * newadj);
void SymbolicFactorization(ETree& tree,const vector<int> & colptr,const vector<int> & rowind,const vector<int> & cc, vector<int> & xlindx, vector<int> & lindx);
void GetLColRowCount(ETree & tree,const int * xadj, const int * adj, vector<int> & cc, vector<int> & rc);
void displayMatrix(vector<int> & xadj, vector<int> & adj);
int ReadAdjacency(const char * pfilename, vector<int> & xadj, vector<int> & adj);

#endif

#endif
