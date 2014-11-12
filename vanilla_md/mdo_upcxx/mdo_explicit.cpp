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
  int elim_step;
  int adj_sz; // # of edges of this node, which is the size of adj
  upcxx::global_ptr<int> adj;
};

void dump_local_nodes(my_node_t *local_nodes, int n);

bool node_comp(my_node_t * & a, my_node_t * & b)
{
  if (a->degree < b->degree) {
    return true;
  } else {
    if(a->degree == b->degree) {
#ifdef USE_RANDOM_TIE_BREAK
      double tmp = fmod(rand(),10.0);
      return (tmp>5.0);
#else
      return a->id < b->id; // use deterministic tie break
#endif
    } else {
      return false;
    }
  }
}



int main(int argc, char *argv[]) 
{
  upcxx::init(&argc, &argv);
  double io_time = mysecond();

  /* initialize random seed: */
  int seed =time(NULL); 
  srand (seed);



  // Don't need xadj any more as its info is stored in my_node_t
  vector<int> xadj;
  // upcxx::shared_array<int> xadj_shared;
  // xadj doesn't need to be shared because it's static and relatively small.
  
  vector<int> local_adj;
  vector<int> ixadj;
  vector<int> iadj;


  int n;
  ReadAdjacencyHB(argv[1], ixadj, iadj);
  n = ixadj.size()-1;

  //expand to asymmetric storage
  ExpandSymmetric(n,&ixadj[0],&iadj[0], xadj, local_adj);


  io_time = mysecond() - io_time;
  double init_time = mysecond();

  // vector<my_node_t> nodes(n);








  upcxx::shared_array<my_node_t> nodes(n);
  nodes.init();



  // YZ: comment out unused variables
  // int initEdgeCnt = 0;
  // int edgeCnt = 0;
  vector<unsigned int> marker(n+1,0);
  unsigned int tag =0;

  //initialize nodes
  for (int i = upcxx::myrank(); i < n; i += upcxx::ranks()) {
    upcxx::global_ptr<my_node_t> cur_node = &nodes[i];
    // cur_node.id = i+1;
    memberof(cur_node, id) = i+1;
    // cur_node.degree = xadj[i+1] - xadj[i];
    memberof(cur_node, degree) = xadj[i+1] - xadj[i];
    // cur_node.elim_step = -1;
    memberof(cur_node, elim_step) = -1;

    int adj_sz = xadj[i+1] - xadj[i];
    upcxx::global_ptr<int> adj = upcxx::allocate<int>(upcxx::myrank(), adj_sz);
    assert(adj != NULL);
    memberof(cur_node, adj) = adj;
    memberof(cur_node, adj_sz) = adj_sz;
    for (int j=0; j<adj_sz; j++) {
      adj[j] = local_adj[xadj[i] - 1 + j];
      // printf("adj[%lu] %d ", i+j, (int)adj[j]);
    }
    /*
    for(int idx = xadj[i]; idx <= xadj[i+1]-1; ++idx) {
      if(local_adj[idx-1]>i+1){
        initEdgeCnt++;
      }
    }
    */
  }

  // Finish reading the data from file and initializing the data structures
  // dump_local_nodes((my_node_t *)&nodes[upcxx::myrank()], nodes.size());

  upcxx::shared_array<int> all_min_degrees(upcxx::ranks());
  upcxx::shared_array<int> all_min_ids(upcxx::ranks());
  all_min_degrees.init();
  all_min_ids.init();
  
  init_time = mysecond() - init_time;
  double mdo_time = mysecond();

  vector<int> schedule;
  // vector< upcxx::global_ptr<my_node_t> > schedule_shared;

  // YZ: loop-carried dependencies between steps
  //process with MD algorithm
  for (int step=1; step<=n; ++step) {
    my_node_t *local_nodes = (my_node_t *)&nodes[upcxx::myrank()];
    my_node_t *my_min = NULL;
    // YZ: need to handle the case when n is not a multiple of ranks()!!
    int local_size = n / upcxx::ranks();
    if (upcxx::myrank() < (n - local_size * upcxx::ranks())) {
      local_size++;
    }
    // printf("step %d, local_size %d\n", step, local_size);
    int cur_min_degree = -1;
    for (int i = 0; i < local_size; i++) {
      my_node_t *cur_node = &local_nodes[i];
      if (cur_node->elim_step == -1) {
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

      upcxx::global_ptr<my_node_t> min_node = &nodes[all_min_ids[min_rank]-1];
      global_min_id = all_min_ids[min_rank]; // memberof(min_node, id);
      memberof(min_node, elim_step) = step;
    }

    upcxx::barrier();
    upcxx::upcxx_bcast(&global_min_id, &global_min_id, sizeof(int), 0);
    assert(global_min_id != -1);
    
    schedule.push_back(global_min_id);

    //update the degree of its reachable set
    upcxx::global_ptr<my_node_t> min_ptr = &nodes[global_min_id-1];
    vector< int > reach(memberof(min_ptr,adj_sz));
    upcxx::copy(memberof(min_ptr,adj).get(),upcxx::global_ptr<int>(&reach[0]),memberof(min_ptr,adj_sz).get());

    // YZ: Use UPC++ to parallel this loop.  There is no data dependencies
    // inside the for loop because the get_reach function does not change
    // the original graph (and the adjacency list) though some nodes' degree
    // might be changed.
#ifdef USE_UPCXX
    // vector< upcxx::global_ptr<my_node_t> >::iterator
    for (auto it = reach.begin();
        it < reach.end();
        it++) {
#else
#pragma omp parallel for
    for (vector<my_node_t *>::iterator it = reach.begin();
         it != reach.end();
         ++it) {
#endif

      if(*it == global_min_id || (*it-1)%(upcxx::ranks())!=upcxx::myrank() ){
        continue;
      }

      // my_node_t *cur_neighbor = *it;
      upcxx::global_ptr<my_node_t> cur_neighbor = &nodes[*it-1];
      assert(cur_neighbor.where()==upcxx::myrank());

      my_node_t * node_ptr = cur_neighbor;

      int old_adj_sz = node_ptr->adj_sz;
      int new_adj_sz = 0;
      int *new_local_adj = new int[reach.size()+old_adj_sz];
      tag++;
      marker[global_min_id]=tag;
      for(int i =0;i<old_adj_sz;++i){
        if(marker[node_ptr->adj[i]]!=tag){
          new_local_adj[new_adj_sz++]=node_ptr->adj[i];
          marker[node_ptr->adj[i]]=tag;
        }
      }

      for(int i =0; i<reach.size();++i){
        if(marker[reach[i]]!=tag){
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
      //delete local_adj;

//
//
//
//      //update the adjacency of that node
//
//      int old_adj_sz = memberof(cur_neighbor,adj_sz);
//      int new_adj_sz = 0;
//      int *new_local_adj = new int[reach.size()+old_adj_sz];
//      upcxx::global_ptr<int> cur_adj = memberof(cur_neighbor,adj);
//      tag++;
//      marker[global_min_id]=tag;
//      for(int i =0;i<old_adj_sz;++i){
//        if(marker[cur_adj[i]]!=tag){
//          new_local_adj[new_adj_sz++]=cur_adj[i];
//          marker[cur_adj[i]]=tag;
//        }
//      }
//
//      for(int i =0; i<reach.size();++i){
//        if(marker[reach[i]]!=tag){
//          new_local_adj[new_adj_sz++]=reach[i];
//          marker[reach[i]]=tag;
//        }
//      }
//
//      memberof(cur_neighbor,adj_sz) = new_adj_sz; 
//      memberof(cur_neighbor, degree) = new_adj_sz-1;
//
//      if(old_adj_sz < new_adj_sz){
//        upcxx::deallocate(cur_adj);
//        upcxx::global_ptr<int> adj = upcxx::allocate<int>(cur_neighbor.where(), new_adj_sz);
//        memberof(cur_neighbor,adj) = adj; 
//      }
//
//      upcxx::copy(upcxx::global_ptr<int>(new_local_adj), memberof(cur_neighbor,adj).get(), new_adj_sz);
//
//      delete new_local_adj;
//      //delete local_adj;

    }
    upcxx::barrier();
  } // close of  for (int step=1; step<=n; ++step)

  mdo_time = mysecond() - mdo_time;

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
    printf("\nMinimum degree ordering algorithm time breakdown on rank 0:\n");
    printf("  io time (read graph from file into memory): %g s\n", io_time);
    printf("  setup time (initialize data structures): %g s\n", init_time);
    printf("  main algorithm time (compute the minimum degree ordering): %g s\n", mdo_time);
    printf("\n");
  }

  upcxx::barrier();
  upcxx::finalize();

  return 0;
}
  
void dump_local_nodes(my_node_t *local_nodes, int n)
{
  int local_size = n / upcxx::ranks();
  if (upcxx::myrank() < (n - local_size * upcxx::ranks())) {
    local_size++;
  }
  for (int i = 0; i < local_size; i++) {
    my_node_t *cur_node = &local_nodes[i];
    fprintf(stdout, "local_nodes[%d], id %d, degree %d, elim_step %d, adj_sz %d, adj: ",
            i, cur_node->id, cur_node->degree, cur_node->elim_step, cur_node->adj_sz);
    for (int j = 0; j < cur_node->adj_sz; j++) {
      fprintf(stdout, "%d ", (int)cur_node->adj[j]);
    }
    fprintf(stdout, "\n");
  }
  printf("\n");
}

