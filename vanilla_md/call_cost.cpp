#include <vector>
#include <list>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "util.h"
#include "ETree.hpp"

//#define verbose

using namespace std;



int main(int argc, char *argv[]) {
  vector<int> ixadj;
  vector<int> iadj;

  if(argc<3){
    cerr<<"Usage is: "<<argv[0]<<" input.file \"ordering\""<<endl;
    return -1;
  }


  ReadAdjacency(argv[1], ixadj, iadj);
  int n = ixadj.size()-1;


  //expand to asymmetric storage
  vector<int> xadj;
  vector<int> adj;
  ExpandSymmetric(n,&ixadj[0],&iadj[0], xadj, adj);

  vector<int> perm;
  {
    string orderstr(argv[2]);
    istringstream iss(orderstr);
    perm.reserve(n);

    int i;
    while(iss>> i){ perm.push_back(i);}
    if(perm.size()!=n){
      cerr<<"Wrong ordering, not the same size as graph"<<endl;
      return -2;
    }  
  }

  vector<int> costc(n);
  double cost2  = GetCost(n,adj.size(),&xadj[0],&adj[0],&perm[0]);
  double cost  = GetCostPerCol(n,adj.size(),&xadj[0],&adj[0],&perm[0],&costc[0]);


  cout<<"Cost(fast) is "<<cost2<<endl;
  cout<<"Cost is "<<cost<<endl;
  cout<<"Cost per col is: ";
  double sum = 0;
  for(int i=0;i<n;++i){
    cout<<" "<<costc[i];
    sum += costc[i];
  }
  cout<<endl;

  assert(sum==cost);

  vector<int> psum(n);
  GetPrefixSum(n,&costc[0],&psum[0]);

  cout<<"Prefix sum of cost is: ";
  for(int i=0;i<n;++i){
    cout<<" "<<psum[i];
  }
  cout<<endl;
  //  assert(psum[n-1]==cost);



}
