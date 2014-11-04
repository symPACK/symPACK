#include <vector>
#include <list>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <iostream>
#include <algorithm>

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#include "timer.h"
#include "util.h"

using namespace std;


#ifndef Add_
#define FORTRAN(name) name
#else
#define FORTRAN(name) name##_
#endif




#define STAT
//#define ALLOC_BOOL
//#define MY_LIST

extern "C" {

void FORTRAN(ordmmd)( int * neqns , int * nadj  , int * xadj  ,
         int * adjncy, int * invp  , int * perm  , int * iwsiz ,
                         int * iwork , int * nofsub, int * iflag);

}




using namespace LIBCHOLESKY;

double mysecond()
{
  struct timeval tv; 
  gettimeofday(&tv, 0); 
  return tv.tv_sec + ((double) tv.tv_usec / 1000000);
}

int main(int argc, char *argv[]) {
  /* initialize random seed: */
  int seed =time(NULL); 
  srand (seed);

  vector<int> ixadj;
  vector<int> iadj;

  if(argc<2){
    cerr<<"Usage is: "<<argv[0]<<" input.file"<<endl;
    return -1;
  }


  ReadAdjacencyHB(argv[1], ixadj, iadj);

  int n = ixadj.size()-1;


  //expand to asymmetric storage
  vector<int> xadj;
  vector<int> adj;
  ExpandSymmetric(n,&ixadj[0],&iadj[0], xadj, adj);

  double mdtime = mysecond();

    int iwsiz = 4*n;
    std::vector<int> iwork (iwsiz);
  vector<int> iperm(n);
  vector<int> perm(n);

    int nadj = adj.size();
    int nofsub =0;
    int iflag =0;

    FORTRAN(ordmmd)( &n , &nadj , &xadj[0], &adj[0], 
            &iperm[0] , &perm[0] , &iwsiz , &iwork[0] , &nofsub, &iflag ) ;


  mdtime = mysecond() - mdtime;


  
  cout<<"Schedule: ";
  for(int step = 1; step<=n;++step){
    cout<<" "<<iperm[step-1];
  }
  cout<<endl;

 

 printf("  main algorithm time (compute the minimum degree ordering): %g s\n", mdtime);

}




