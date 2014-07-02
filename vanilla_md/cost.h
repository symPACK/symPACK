#ifndef COST_HEADER
#define COST_HEADER


#ifdef __cplusplus
#include <vector>
#include "ETree.hpp"

using namespace std;

extern "C" {
#endif
extern double GetCost(int n, int nnz, int * xadj, int * adj,int * perm);
#ifdef __cplusplus
}


void GetPermutedGraph(int n, int nnz, int * xadj, int * adj, int * perm, int * newxadj, int * newadj);
void SymbolicFactorization(ETree& tree,const vector<int> & colptr,const vector<int> & rowind,const vector<int> & cc, vector<int> & xlindx, vector<int> & lindx);
void GetLColRowCount(ETree & tree,const int * xadj, const int * adj, vector<int> & cc, vector<int> & rc);

#endif

#endif
