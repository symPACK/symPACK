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

  if(argc<2){
    cerr<<"Usage is: "<<argv[0]<<" input.file"<<endl;
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
  int  perm[16] = { 3, 1, 6, 2, 4, 5, 7, 13, 10, 15, 11, 12, 8, 16, 9, 14};

  double cost  = GetCost(n,adj.size(),&xadj[0],&adj[0],&perm[0]);

    cout<<cost<<endl;
}
