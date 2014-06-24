#ifndef COST_HEADER
#define COST_HEADER

#ifdef __cplusplus
extern "C" {
#endif
extern double GetCost(int n, int nnz, int * xadj, int * adj,int * perm);
#ifdef __cplusplus
}
#endif

#endif
