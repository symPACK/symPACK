#include <vector>
#include <list>
#include <iostream>
#include <algorithm>

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#include "util.h"

using namespace std;

double mysecond()
{
  struct timeval tv; 
  gettimeofday(&tv, 0); 
  return tv.tv_sec + ((double) tv.tv_usec / 1000000);
}

struct node_t{
  int id;
  int degree;
  int elim_step;
};


bool node_comp(node_t * & a, node_t * & b){
  if (a->degree>b->degree)
    return false;
  else{
    if(a->degree == b->degree){
//      double tmp = fmod(rand(),10.0);
//      return (tmp>5.0);
      return false;
    }
    else{
      return true;
    }
  }
}



void get_reach(const vector<int> & xadj,const vector<int> & adj,vector<node_t*> & nodes, const node_t & min_node, const int elim_step, vector<node_t*> & reach_set){
  //this array is used to mark the node so that 
  //they are not added twice into the reachable set
  vector<bool> explored(nodes.size(),false);
  //this list contains the nodes to explore
  //list<node_t *> explore_set;
  vector<node_t *> explore_set;

  //initialize explore_set with the neighborhood
  // of min_node in the original graph 
  int beg = xadj[min_node.id-1];  
  int end = xadj[min_node.id]-1;
  for(int i = beg;i<=end;++i){
    if(adj[i-1]!=0){
      node_t * next_node = nodes[adj[i-1]-1];
      explore_set.push_back(next_node);
      explored[adj[i-1]-1]=true;
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
            explored[adj[i-1]-1]=true;
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

  vector<node_t *> nodes(n);
  vector<node_t *> node_heap(n);

  int initEdgeCnt = 0;
  int edgeCnt = 0;

  //initialize nodes
  for(int i=0;i<n;++i){
    nodes[i] = new node_t;
    node_t & cur_node = *nodes[i];
    cur_node.id = i+1;
    cur_node.degree = xadj[i+1] - xadj[i] -1;
    cur_node.elim_step = -1;

    for(int idx = xadj[i]; idx <= xadj[i+1]-1;++idx){
      if(adj[idx-1]>i+1){
        initEdgeCnt++;
      }
    }
    node_heap[i] = &cur_node;
    //node_heap.push_back(&cur_node);


    //cout<<"deg("<<cur_node.id<<")="<<cur_node.degree<<endl;
  }

  double mdtime = mysecond();

  vector<node_t*> schedule; 
  vector<int> perm; 
  vector<int> invperm(n); 
  //process with MD algorithm
  auto node_heap_end = node_heap.end();
  for(int step = 1; step<=n;++step){
    //int min_degree = -1;
    //int min_id = 0;
    //for(int i=0;i<node_heap_end-node_heap.begin();++i){
    //  if(min_degree == -1 || min_degree >= node_heap[i]->degree){
    //    min_id = i;
    //    min_degree = node_heap[i]->degree;
    //  }
    //}

    auto it = std::min_element(node_heap.begin(),node_heap_end,node_comp);
    node_t & min_node = *(*it);

    //assert(min_node.degree==min_degree);

    *it = *(node_heap_end-1);
    //remove one node by moving the end iterator down one by one element
    node_heap_end--;
    schedule.push_back(&min_node);

    perm.push_back(min_node.id);

    min_node.elim_step = step;
    //update the degree of its reachable set
    vector<node_t *> reach;
    get_reach(xadj,adj,nodes,min_node,step,reach);

    //cout<<"Reach of node "<<min_node.id<<endl;
    //for(auto it2 = reach.begin(); it2!=reach.end();++it2){
    //  cout<<(*it2)->id<<" ";
    //}
    //cout<<endl;

    for(auto it = reach.begin();it!=reach.end();++it){
      node_t * cur_neighbor = *it;
      if(cur_neighbor->id != min_node.id){
        vector<node_t *> nghb_reach;
        get_reach(xadj,adj,nodes,*cur_neighbor,step+1,nghb_reach);
        cur_neighbor->degree = nghb_reach.size();

    //cout<<"Nghb reach of node "<<cur_neighbor->id<<endl;
    //for(auto it2 = nghb_reach.begin(); it2!=nghb_reach.end();++it2){
    //  cout<<(*it2)->id<<" ";
    //}
    //cout<<endl;
      }
    }
  }


  mdtime = mysecond() - mdtime;



  cout<<"Schedule: ";
  for(int step = 1; step<=n;++step){
    cout<<" "<<schedule[step-1]->id;
  }
  cout<<endl;

 printf("  main algorithm time (compute the minimum degree ordering): %g s\n", mdtime);

}




