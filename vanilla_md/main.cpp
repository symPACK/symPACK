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


struct node_t{
  int id;
  int degree;
  int elim_step;
};

bool node_comp(node_t * & a, node_t * & b){
  //if (a->degree>b->degree)
  if (a->degree>b->degree)
    return true;
  else{
    if(a->degree == b->degree){
      double tmp = fmod(rand(),10.0);
      return (tmp>5.0);
    }
    else{
      return false;
    }
  }
}



void get_reach(const vector<int> & xadj,const vector<int> & adj,vector<node_t*> & nodes, const node_t & min_node, const int elim_step, list<node_t*> & reach_set){
  //this array is used to mark the node so that 
  //they are not added twice into the reachable set
  vector<int> explored(nodes.size(),0);
  //this list contains the nodes to explore
  list<node_t *> explore_set;

  //initialize explore_set with the neighborhood
  // of min_node in the original graph 
  int beg = xadj[min_node.id-1];  
  int end = xadj[min_node.id]-1;
  for(int i = beg;i<=end;++i){
    if(adj[i-1]!=0){
      node_t * next_node = nodes[adj[i-1]-1];
      explore_set.push_back(next_node);
      explored[adj[i-1]-1]=1;
    }
    else{
      abort();
    }
  }  

  //now find path between min_nodes and other nodes
  while(explore_set.size()>0){
    //pop a node
    node_t * cur_node = explore_set.back();

    explore_set.pop_back();

    if(cur_node->id == min_node.id ){
      continue;
    }

    if(cur_node->elim_step == -1 ){
      reach_set.push_back(cur_node);
    }
    else{
      int beg = xadj[cur_node->id-1];  
      int end = xadj[cur_node->id]-1;
      for(int i = beg;i<=end;++i){
        if(adj[i-1]!=0){
          if(!explored[adj[i-1]-1]){
            node_t * next_node = nodes[adj[i-1]-1];
            explore_set.push_back(next_node);
            explored[adj[i-1]-1]=1;
          }
        }
        else{
          abort();
        }
      }
    }
  }
}



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
  vector<node_t *> nodes(n);
  vector<node_t *> node_heap(n);

  int initEdgeCnt = 0;
  int edgeCnt = 0;

  //initialize nodes
  for(int i=0;i<n;++i){
    nodes[i] = new node_t;
    node_t & cur_node = *nodes[i];
    cur_node.id = i+1;
    cur_node.degree = xadj[i+1] - xadj[i];
    cur_node.elim_step = -1;

    for(int idx = xadj[i]; idx <= xadj[i+1]-1;++idx){
      if(adj[idx-1]>i+1){
        initEdgeCnt++;
      }
    }
    node_heap[i] = &cur_node;
  }

  //make heap
  std::make_heap(node_heap.begin(),node_heap.end(),node_comp);

  vector<node_t*> schedule; 
  vector<int> perm; 
  vector<int> invperm(n); 
  //process with MD algorithm
  for(int step = 1; step<=n;++step){
    //sort heap
    std::make_heap(node_heap.begin(),node_heap.end(),node_comp);
    //std::sort_heap (node_heap.begin(),node_heap.end(),node_comp);
    //get the min degree node
    std::pop_heap (node_heap.begin(),node_heap.end(),node_comp);
    node_t & min_node = *node_heap.back();
    node_heap.pop_back();
    schedule.push_back(&min_node);

    perm.push_back(min_node.id);

    min_node.elim_step = step;
#ifdef verbose
    cout<<"Node "<<min_node.id<<" scheduled at step "<<step<<" degree is "<<min_node.degree<<endl;
#endif
    //update the degree of its reachable set
    list<node_t *> reach;
    get_reach(xadj,adj,nodes,min_node,step,reach);

    for(list<node_t *>::iterator it = reach.begin();it!=reach.end();++it){
      node_t * cur_neighbor = *it;
      list<node_t *> nghb_reach;
      get_reach(xadj,adj,nodes,*cur_neighbor,step+1,nghb_reach);
#ifdef verbose
      cout<<"     Node "<<cur_neighbor->id<<" degree "<<cur_neighbor->degree<<" -> "<<nghb_reach.size()<<" {";
      for(list<node_t *>::iterator it2 = nghb_reach.begin();it2!=nghb_reach.end();++it2){
        cout<<" "<<(*it2)->id;
      }
      cout<<" }"<<endl;
#endif
      cur_neighbor->degree = nghb_reach.size();

    }
  }





  cout<<"Schedule: ";
  for(int step = 1; step<=n;++step){
    cout<<" "<<schedule[step-1]->id;
  }
  cout<<endl;



  for(int i=0;i<n;++i){
    delete nodes[i];
  }

  double cost  = GetCost(n,adj.size(),&xadj[0],&adj[0],&perm[0]);

  cout<<"Cost is "<<cost<<endl;
}




