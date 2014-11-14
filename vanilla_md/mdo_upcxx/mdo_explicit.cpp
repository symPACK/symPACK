/**
 * Minimum degree ordering example -- compute a permutation of the rows
 * or columns of a sparse matrix that minimizes the fill-ins (non-zeros)
 * in Cholesky factorization
 */

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <sstream>
#include <climits>

#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#define USE_UPCXX

#ifdef USE_UPCXX
#include <upcxx.h>
#endif

#include "util.h"
#include "timer.h"

using namespace std;

double mysecond()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + ((double) tv.tv_usec / 1000000);
}

struct my_node_t {
  int id;
  int degree;
  int label;
  int adj_sz; // # of edges of this node, which is the size of adj
  upcxx::global_ptr<int> adj;
};

void find_indist(upcxx::shared_array<node_t> & nodes, upcxx::global_ptr<node_t> & ref_node, std::vector<int> & marker, int & tag){
  int ref_adj_sz = memberof(ref_node,adj_sz);
  int ref_id = memberof(ref_node,id); 
  vector<int> local_adj(ref_adj_sz);
  upcxx::copy(memberof(ref_node,adj).get(),upcxx::global_ptr<int>(&local_adj[0]),local_adj.size());

    //mark nodes
    tag++;
    if(tag==INT_MAX){tag=0;}
    marker[ref_id] = tag;
    for(int i=0;i<local_adj.size();++i){ marker[local_adj[i]] = tag; }



  for(int i=upcxx::myrank();i<nodes.size();i+=upcxx::ranks()){
    upcxx::global_ptr<node_t> cur_ptr = &nodes[i];
    int adj_sz = memberof(cur_ptr,adj_sz); 
    if(i == ref_id || adj_sz!=ref_adj_sz){
      continue;
    }




    int * loc_adj = memberof(cur_ptr,adj).get(); 
    bool indist = true;
    for(int j=0;j<adj_sz;++j){
      if(marker[loc_adj[j]]!=tag && marker[loc_adj[j]]!=INT_MAX){
        indist = false;
        break;
      }
    }

    if(indist){
      marker[i+1] = INT_MAX;
    }
  }

  //all local indistinguishable nodes are marked
  //we need to do a reduce 
  upcxx::upcxx_reduce(&marker[0],&marker[0],marker.size(),ref_node.where(),UPCXX_MAX,UPCXX_INT);
  upcxx::upcxx_bcast(&marker[0], &marker[0], marker.size()*sizeof(int), ref_node.where());


////  logfile<<"Nodes "<<ref_id<<" and [";
////  for(int i=0;i<local_adj.size();++i){
////    if(marker[local_adj[i]]==INT_MAX){
////      logfile<<local_adj[i]<<" ";
////    }
//////    cerr<<marker[local_adj[i]]<<" ";
////  }
//////  cerr<<endl;
////  logfile<<"] are indist."<<endl;
    

}


int main(int argc, char *argv[]) 
{
  upcxx::init(&argc, &argv);
//  double io_time = mysecond();
//
//  /* initialize random seed: */
//  int seed =time(NULL); 
//  srand (seed);
//
//
//
//  // Don't need xadj any more as its info is stored in node_t
//  vector<int> xadj;
//  // upcxx::shared_array<int> xadj_shared;
//  // xadj doesn't need to be shared because it's static and relatively small.
//  
//  vector<int> local_adj;
//  vector<int> ixadj;
//  vector<int> iadj;
//
//
//  int n;
//  ReadAdjacencyHB(argv[1], ixadj, iadj);
//  n = ixadj.size()-1;
//
//  //expand to asymmetric storage
//  ExpandSymmetric(n,&ixadj[0],&iadj[0], xadj, local_adj);
//
//
//  io_time = mysecond() - io_time;
//  double init_time = mysecond();
//
//  // vector<node_t> nodes(n);
//
//  upcxx::shared_array<node_t> nodes(n);
//  nodes.init();



  stringstream filename;
  filename<<"logTest"<<upcxx::myrank();
  logfile.open(filename.str().c_str());

  bool doMassElim = false;
  if(argc>2){
    if(atoi(argv[2])==1){
      doMassElim = true;
    }
  }


TIMER_START(io);



  int n,nnz;
  upcxx::shared_array<node_t> * nodes_ptr;
  ReadAdjacencyHB_PARA(argv[1], nodes_ptr, n, nnz);

  //io_time = mysecond() - io_time;
  upcxx::shared_array<node_t> & nodes = *nodes_ptr;

  logfile<<"Input file read"<<endl;

    upcxx::barrier(); 
TIMER_STOP(io);

  TIMER_START(init);
  ExpandSymmetric_PARA(n, nodes);

  logfile<<"Matrix expanded"<<endl;
  double init_time = mysecond();











  // YZ: comment out unused variables
  // int initEdgeCnt = 0;
  // int edgeCnt = 0;
  //vector<unsigned int> marker(n+1,0);
  //unsigned int tag =0;

  // Finish reading the data from file and initializing the data structures
  // dump_local_nodes((node_t *)&nodes[upcxx::myrank()], nodes.size());

  upcxx::shared_array<int> all_min_degrees(upcxx::ranks());
  upcxx::shared_array<int> all_min_ids(upcxx::ranks());
  all_min_degrees.init();
  all_min_ids.init();
  
  init_time = mysecond() - init_time;
  TIMER_STOP(init);
  TIMER_START(mdo_time);
  double mdo_time = mysecond();

  vector<int> schedule;
  vector<int> marker(n+1,0);
  int tag=0;
  // vector< upcxx::global_ptr<node_t> > schedule_shared;

  // YZ: loop-carried dependencies between steps
  //process with MD algorithm
  //for (int step=1; step<=n; ++step) {
  int label=1;
  int step = 1;
  while(label<=n) {
    node_t *local_nodes = (node_t *)&nodes[upcxx::myrank()];
    node_t *my_min = NULL;
    // YZ: need to handle the case when n is not a multiple of ranks()!!
    int local_size = n / upcxx::ranks();
    if (upcxx::myrank() < (n - local_size * upcxx::ranks())) {
      local_size++;
    }
    // printf("step %d, local_size %d\n", step, local_size);
    int cur_min_degree = -1;
    for (int i = 0; i < local_size; i++) {
      node_t *cur_node = &local_nodes[i];
      if (cur_node->label == 0) {
        if (cur_min_degree == -1 || cur_node->degree < cur_min_degree) {
          cur_min_degree = cur_node->degree;
          my_min = cur_node;
        }
      }
    }
    
    if (my_min != NULL) {
      all_min_degrees[upcxx::myrank()] = my_min->degree;
      all_min_ids[upcxx::myrank()] = my_min->id;
    } else {
      all_min_degrees[upcxx::myrank()] = INT_MAX;
      all_min_ids[upcxx::myrank()] = -1;
    }

    upcxx::barrier();
    

    int global_min_id;
    int global_min_degree;
    
    if (upcxx::myrank() == 0) {
      int cur_min_degree = INT_MAX;
      int min_rank = 0;
      // find the global min node
      for (int i=0; i<upcxx::ranks(); i++) {
        if (cur_min_degree > all_min_degrees[i]) {
          cur_min_degree = all_min_degrees[i];
          min_rank = i;
        }
      }

      upcxx::global_ptr<node_t> min_node = &nodes[all_min_ids[min_rank]-1];
      global_min_id = all_min_ids[min_rank]; // memberof(min_node, id);
      memberof(min_node, label) = label;
      global_min_degree = cur_min_degree;
    }

    upcxx::barrier();
    upcxx::upcxx_bcast(&global_min_id, &global_min_id, sizeof(int), 0);
    upcxx::upcxx_bcast(&global_min_degree, &global_min_degree, sizeof(int), 0);
    assert(global_min_id != -1);



    if(global_min_degree == n-label){
      //add remaining nodes to the schedule
      //TODO this is very inefficient
      for(int i = 0;i<n;i++){
        upcxx::global_ptr<node_t> cur_ptr = &nodes[i];
        if(memberof(cur_ptr,label)==0){
          schedule.push_back(i+1);
          label++;
        }
      }
      break;
    }   
    schedule.push_back(global_min_id);
    label++;

    //update the degree of its reachable set
    upcxx::global_ptr<node_t> min_ptr = &nodes[global_min_id-1];
    vector< int > reach(memberof(min_ptr,adj_sz));
    //if(upcxx::myrank()==min_ptr.where()){
      upcxx::copy(memberof(min_ptr,adj).get(),upcxx::global_ptr<int>(&reach[0]),memberof(min_ptr,adj_sz).get());
    //}
    //upcxx::upcxx_bcast(&reach[0], &reach[0], reach.size()*sizeof(int), min_ptr.where());




    if(doMassElim){
      find_indist(nodes, min_ptr, marker, tag);
    }






    // YZ: Use UPC++ to parallel this loop.  There is no data dependencies
    // inside the for loop because the get_reach function does not change
    // the original graph (and the adjacency list) though some nodes' degree
    // might be changed.
#ifdef USE_UPCXX
    // vector< upcxx::global_ptr<node_t> >::iterator
    for (auto it = reach.begin();
        it < reach.end();
        it++) {
#else
#pragma omp parallel for
    for (vector<node_t *>::iterator it = reach.begin();
         it != reach.end();
         ++it) {
#endif


      upcxx::global_ptr<node_t> cur_neighbor = &nodes[*it-1];
      //if this node is indist, label it now
      if(*it != global_min_id && marker[*it]==INT_MAX){
        if(cur_neighbor.where()==upcxx::myrank()){
          memberof(cur_neighbor,label) = label;
        }
        schedule.push_back(*it);
        label++;
        continue;
      }


      if(*it == global_min_id || (*it-1)%(upcxx::ranks())!=upcxx::myrank() ){
        continue;
      }


      node_t * node_ptr = cur_neighbor;



      int old_adj_sz = node_ptr->adj_sz;
      int new_adj_sz = 0;
      int *new_local_adj = new int[reach.size()+old_adj_sz];
      tag++;
      if(tag==INT_MAX){tag=0;}
      marker[global_min_id]=tag;
      for(int i =0;i<old_adj_sz;++i){
        int mark = marker[node_ptr->adj[i]];
        if(mark!=tag && mark!=INT_MAX){
          new_local_adj[new_adj_sz++]=node_ptr->adj[i];
          marker[node_ptr->adj[i]]=tag;
        }
      }

      for(int i =0; i<reach.size();++i){
        int mark = marker[reach[i]];
        if(mark!=tag && mark!=INT_MAX){
          new_local_adj[new_adj_sz++]=reach[i];
          marker[reach[i]]=tag;
        }
      }

      node_ptr->adj_sz = new_adj_sz; 
      node_ptr->degree = new_adj_sz;

      if(old_adj_sz < new_adj_sz){
        upcxx::deallocate(node_ptr->adj);
        upcxx::global_ptr<int> adj = upcxx::allocate<int>(cur_neighbor.where(), new_adj_sz);
        node_ptr->adj = adj; 
      }

      upcxx::copy(upcxx::global_ptr<int>(new_local_adj), node_ptr->adj, new_adj_sz);

      delete new_local_adj;

    }
    step++;
    upcxx::barrier();
  } // close of  for (int step=1; step<=n; ++step)

  mdo_time = mysecond() - mdo_time;
  TIMER_STOP(mdo_time);

//  for (int i = 0; i < upcxx::ranks(); i++) {
    if (upcxx::myrank() == 0) {
      cout << "\n";
      cout<<"Rank " << upcxx::myrank() << " Schedule: ";
      for (int step=1; step<=n; ++step) {
        cout << " " << schedule[step-1];
      }
      cout << "\n";
    }
//    upcxx::barrier();
//  }

  if (upcxx::myrank() == 0) {
    printf("\nMinimum degree ordering done in %d vs %d steps\n",step,n);
    printf("\nMinimum degree ordering algorithm time breakdown on rank 0:\n");
//    printf("  io time (read graph from file into memory): %g s\n", io_time);
    printf("  setup time (initialize data structures): %g s\n", init_time);
    printf("  main algorithm time (compute the minimum degree ordering): %g s\n", mdo_time);
    printf("\n");
  }

  upcxx::barrier();
  upcxx::finalize();

  return 0;
}
  


