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

  if(argc<2){
    cerr<<"Usage is: "<<argv[0]<<" input.file"<<endl;
    return -1;
  }


  ReadAdjacency(argv[1], ixadj, iadj);
  int n = ixadj.size()-1;


  //expand to asymmetric storage
  vector<int> xadj;
  vector<int> adj;
  ExpandSymmetric(n,&ixadj[0],&iadj[0], xadj, adj);

  //initialize perm with natural ordering
  vector<int> perm;
  for(int i =1; i<=n;++i){perm.push_back(i);}

  random_shuffle(perm.begin(),perm.end());

  vector<int> costc(n);

  int oldpos = -1;
  int begin = 0;
  double cost = 0.0;
  for(int step = 1; step<=n;++step){

    cout<<"        Perm is: ";
    for(int i=0;i<n;++i){
      cout<<" "<<perm[i];
    }
    cout<<endl;

    //compute the cost
    cost = GetCostPerCol(n,adj.size(),&xadj[0],&adj[0],&perm[0],&costc[0]);
    cout<<"        Cost is "<<cost<<endl;
    cout<<"        Cost per col is: ";
    for(int i=0;i<n;++i){
      cout<<" "<<costc[i];
    }
    cout<<endl;


    //find the largest jump in fill
    int * it = max_element(&costc[begin],&costc[0]+n);
    //find the first big jump in fill
    //    int * it = upper_bound(&costc[begin],&costc[0]+n,0);
    int pos = distance(&costc[0],it);
    if(oldpos!=pos){   

      cout<<"        Max at pos: "<<pos<<endl; 
      int col = perm[pos];

      //place it at the end
      perm.erase(perm.begin()+pos);
      perm.push_back(col);
      //    begin++;

      oldpos=pos;
    }
    else{
//      break;
    }

///    int pos=0;
///    while(costc[pos++]==0);
/////    int * it = upper_bound(&costc[begin],&costc[0]+n,0);
/////    int pos = distance(&costc[0],it);
///    cout<<"        First at pos: "<<pos<<endl; 
///    int col = perm[pos];
///
///    //try to place it either at the end or the beginning (the one which pushes the fill to the right)
///    
///    //place it at the end
///    perm.erase(perm.begin()+pos);
///    perm.insert(perm.end(),col);
///
///
/////    //get cost
/////    GetCostPerCol(n,adj.size(),&xadj[0],&adj[0],&perm[0],&costc[0]);
/////    int pose=0;
/////    while(costc[pose++]==0);
/////
/////    //place it at the beginning
/////    perm.erase(perm.begin()+n-1);
/////    perm.insert(perm.begin(),col);
/////    //get cost
/////    GetCostPerCol(n,adj.size(),&xadj[0],&adj[0],&perm[0],&costc[0]);
/////    int posb=0;
/////    while(costc[posb++]==0);
/////
/////    if(posb>pose){
/////      //put it at the end
/////      perm.erase(perm.begin());
/////      perm.insert(perm.end(),col);
/////    }

  }


  cost  = GetCost(n,adj.size(),&xadj[0],&adj[0],&perm[0]);

  cout<<"Cost is "<<cost<<endl;

  cout<<"Schedule is: ";
  for(int i=0;i<n;++i){
    cout<<" "<<perm[i];
  }
  cout<<endl;
}
