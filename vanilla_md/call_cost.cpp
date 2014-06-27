#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <sstream>

#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "cost.h"
#include "ETree.hpp"

//#define verbose

using namespace std;
 



int main(int argc, char *argv[]) {
  /* initialize random seed: */
  int seed =time(NULL); 
  srand (seed);

  vector<int> xadj;
  vector<int> adj;

  if(argc<3){
    cerr<<"Usage is: "<<argv[0]<<" input.file \"ordering\""<<endl;
    return -1;
  }

  string filename(argv[1]);
  ifstream infile;
  infile.open(filename.c_str());

  string line;
  //read xadj on the first line of the input file
  if(getline(infile, line))
  {
    istringstream iss(line);
    int i;
    while(iss>> i){ xadj.push_back(i);}
  }    
  else{
    return -2;
  }

  //read adj on the second line of the input file
  if(getline(infile, line))
  {
    istringstream iss(line);
    int i;
    while(iss>> i){ adj.push_back(i);}
  }    
  else{
    return -2;
  }


  infile.close();

  int n = xadj.size()-1;

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
  


  double cost  = GetCost(n,adj.size(),&xadj[0],&adj[0],&perm[0]);

    cout<<cost<<endl;
}
