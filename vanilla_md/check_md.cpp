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

#include <assert.h> 

#include "cost.h"
#include "ETree.hpp"

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

  //read the ordering
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
 
  //recompute adj and xadj
  int nnz = adj.size();


  vector<int> invperm(n);

  for(int step = 1; step<=n;++step){
    invperm[perm[step-1]-1] = step;
  }




  vector<int> newxadj(n+1);
  vector<int> newadj(nnz);
  GetPermutedGraph(n,nnz,&xadj[0], &adj[0], &perm[0], &newxadj[0], &newadj[0]);


////  //Perform the symbolic factorization
////  ETree  tree;
////  tree.ConstructETree(n,&newxadj[0],&newadj[0]);
////
////  vector<int> cc,rc;
////  GetLColRowCount(tree,&newxadj[0],&newadj[0],cc,rc);
////
////  vector<int> xlindx,lindx;
////
////  
////  vector<int> lxadj,ladj;
////
////  lxadj.resize(n+1);
////  //ladj.resize(newadj.size()+n);
////  lxadj[0] = 1;
////  for(int i =2; i<=n+1;++i){
////    ladj.insert(ladj.end(),&newadj[newxadj[i-2]-1],&newadj[newxadj[i-1]-1]);
////    ladj.push_back(i-1);
////    lxadj[i-1] = ladj.size()+1; 
////  }
////
////
////  cout<<"newxadj: ";
////  for(int step = 1; step<=newxadj.size();++step){
////    cout<<" "<<newxadj[step-1];
////  }
////  cout<<endl;
////  cout<<"newadj: ";
////  for(int step = 1; step<=newadj.size();++step){
////    cout<<" "<<newadj[step-1];
////  }
////  cout<<endl;
////
////
////
////
////
////  cout<<"lxadj: ";
////  for(int step = 1; step<=lxadj.size();++step){
////    cout<<" "<<lxadj[step-1];
////  }
////  cout<<endl;
////  cout<<"ladj: ";
////  for(int step = 1; step<=ladj.size();++step){
////    cout<<" "<<ladj[step-1];
////  }
////  cout<<endl;
////
////
////
////  SymbolicFactorization(tree,lxadj,ladj,cc, xlindx, lindx);
////
////  //Compute the boolean structure matrix
  vector< vector<bool> > matrix(n, vector<bool>(n,false));

  //Fill the matrix with the content of newadj
  for(int col=1;col<=n;++col){
    int fi = newxadj[col-1];
    int li = newxadj[col]-1;
      
    vector<bool> & colstruct = matrix[col-1]; 
   
    colstruct[col-1] = true;
 
    for(int i=fi; i<=li;++i){
      int row = newadj[i-1];
      
      vector<bool> & rowstruct = matrix[row-1]; 

      colstruct[row-1] = true;
      rowstruct[col-1] = true;

    }

//    assert(std::count(colstruct.begin(),colstruct.end(),true)==cc[col-1]);

  }
////
////  for(int row = 1; row <=n; ++row){
////    for(int col = 1; col <=n; ++col){
////      vector<bool> & colstruct = matrix[col-1]; 
////      cerr<<" "<<colstruct[row-1]?1:0;
////    }
////    cerr<<endl;
////  }
////
////
////
////if(0){
////
//// //Compute the boolean structure matrix
////  vector< vector<bool> > matrix2(n, vector<bool>(n,false));
////
////  //Fill the matrix with the content of lindx
////  for(int col=1;col<=n;++col){
////    int fi = xlindx[col-1];
////    int li = xlindx[col]-1;
////      
////    vector<bool> & colstruct = matrix2[col-1]; 
////    
////    for(int i=fi; i<=li;++i){
////      int row = lindx[i-1];
////      
//////      vector<bool> & rowstruct = matrix2[row-1]; 
////
////      colstruct[row-1] = true;
//////      rowstruct[col-1] = true;
////
////    }
////
//////    assert(std::count(colstruct.begin(),colstruct.end(),true)==cc[col-1]);
////
////  }
////
////    cerr<<endl;
////    cerr<<endl;
////    cerr<<endl;
////  for(int row = 1; row <=n; ++row){
////    for(int col = 1; col <=n; ++col){
////      vector<bool> & colstruct = matrix2[col-1]; 
////      cerr<<" "<<colstruct[row-1]?1:0;
////    }
////    cerr<<endl;
////  }
////
////}







  //Simulate the factorization and check the degrees
  for(int col = 1; col <=n; ++col){
    vector<bool> & colstruct = matrix[col-1]; 

    vector<bool>::iterator it = colstruct.begin();
    advance(it,col-1);
    int deg = std::count(it,colstruct.end(),true)-1;

    cout<<"Degree of col "<<perm[col-1]<<" is "<<deg<<endl;
    //assert(deg +1  == cc[col-1]);

//    cout<<"Degree of col "<<perm[col-1]<<" is "<<deg<<endl;
    //Check that it is indeed the minimum degree
    for(int ncol = col+1;ncol<=n;++ncol){
      //if that column is updated
        vector<bool> & ncolstruct = matrix[ncol-1]; 

        vector<bool>::iterator nit = ncolstruct.begin();
        advance(nit,col-1);
        int ndeg = std::count(nit,ncolstruct.end(),true)-1;

      if( colstruct[ncol-1]){
        cout<<"      Degree of col "<<perm[ncol-1]<<" is "<<ndeg<<endl;
//        for(int nrow = col;nrow<=n;++nrow){
//          cout<<" "<<perm[nrow-1];
//        }
//        cout<<endl;
//        for(int nrow = col;nrow<=n;++nrow){
//          if(perm[nrow-1]>9){
//            cout<<"  "<<ncolstruct[nrow-1];
//          }
//          else{
//            cout<<" "<<ncolstruct[nrow-1];
//          }
//        }
//        cout<<endl;
//        cout<<endl;


        for(int nrow = col;nrow<=n;++nrow){
          ncolstruct[nrow-1] = ncolstruct[nrow-1] || colstruct[nrow-1];
        }

      }

        if(ndeg<deg){
          cerr<<"This ordering is not a minimum degree ordering"<<endl;
          cerr<<"Column "<<ncol<<" has a lower degree than column "<<col<<endl;
          return -1;
        }

    }

  }





 






  cerr<<"This ordering is a minimum degree ordering"<<endl;
  //This was indeed a md ordering, return without error
  return 0;


}




