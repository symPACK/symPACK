/**
 * Minimum degree ordering example -- compute a permutation of the rows
 * or columns of a sparse matrix that minimizes the fill-ins (non-zeros)
 * in Cholesky factorization
 */

#include <vector>
#include <stack>
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

void dump_local_nodes(node_t *local_nodes, int n);


void get_reach(upcxx::shared_array<node_t> &nodes,
               const int root_node_id,
               const unsigned int tag, 
               vector< upcxx::global_ptr<node_t> > &reach_set)
{
  //this array is used to mark the node so that 
  //they are not added twice into the reachable set
  vector<bool> explored(nodes.size(), false);

  //this list contains the nodes to explore
  stack< upcxx::global_ptr<node_t> > explore_set;
  
  //initialize explore_set with the neighborhood
  // of min_node in the original graph 
//  node_t min_node = nodes[min_node_id-1];
//  // copy adj to local
//  int *local_adj = new int[min_node.adj_sz];
//  upcxx::copy(min_node.adj, upcxx::global_ptr<int>(local_adj), min_node.adj_sz);
//
//  for(int i = 0; i < min_node.adj_sz; ++i){
//    int curr_adj = local_adj[i];
//    if(curr_adj != 0) {
//      // node_t * next_node = &nodes[curr_adj-1];
//      upcxx::global_ptr<node_t> next_node = &nodes[curr_adj-1];
//      explore_set.push_back(next_node);
//      explored[curr_adj-1] = true;
//    } else {
//      fprintf(stderr, "get_reach() Fatal error 1: adj[%d]=0 \n", i-1);
//      exit(1);
//    }
//  }  
//  delete local_adj;

  upcxx::global_ptr<node_t> root_node_p = &nodes[root_node_id-1];
  explore_set.push( root_node_p );
//  memberof(root_node_p,marker)=tag;
  explored[root_node_id-1]=true;

  //now find path between min_nodes and other nodes
  while (explore_set.size()>0) {
    //pop a node
    upcxx::global_ptr<node_t> cur_node_p = explore_set.top();
    explore_set.pop();
    assert(cur_node_p.raw_ptr() != NULL);
    node_t cur_node = *cur_node_p;

//    if (cur_node.id == min_node_id) {
//      continue;
//    }
//
//    if (cur_node.elim_step == -1) {
//      reach_set.push_back(cur_node_p);
//    } else {
//      int *local_adj = new int[cur_node.adj_sz];
//      upcxx::copy(cur_node.adj, upcxx::global_ptr<int>(local_adj), cur_node.adj_sz);
//      for(int i = 0; i< cur_node.adj_sz; ++i){
//        int curr_adj = local_adj[i];
//        if (curr_adj != 0) {
//          if (!explored[curr_adj-1]) {
//            upcxx::global_ptr<node_t> next_node = &nodes[curr_adj-1];
//            explore_set.push(next_node);
//            explored[curr_adj-1] = true;
//          }
//        } else {
//          fprintf(stderr, "get_reach() Fatal error 2: adj[%d]=0 \n", i-1);
//          exit(1);
//        }
//      }
//      delete local_adj;
//    }

      int *local_adj = new int[cur_node.adj_sz];
      upcxx::copy(cur_node.adj, upcxx::global_ptr<int>(local_adj), cur_node.adj_sz);
      for(int i = 0; i< cur_node.adj_sz; ++i){
        int curr_adj = local_adj[i];

        if (curr_adj != 0) {
//          if (next_node.marker != tag) {
          if (!explored[curr_adj-1]) {
            upcxx::global_ptr<node_t> next_node_p = &nodes[curr_adj-1];
            if(memberof(next_node_p,label).get()==0){
              reach_set.push_back(next_node_p);
            }
            else{
              explore_set.push(next_node_p);
            }
//            memberof(next_node_p,marker)=tag;
            explored[curr_adj-1] = true;
          }
        } else {
          fprintf(stderr, "get_reach() Fatal error 2: adj[%d]=0 \n", i-1);
          exit(1);
        }
      }
      delete local_adj;



  }
}

int main(int argc, char *argv[]) 
{



  upcxx::init(&argc, &argv);


  if(argc<2){
    cerr<<"Usage is: "<<argv[0]<<" input.file"<<endl;
    return -1;
  }
  stringstream filename;
  filename<<"logTest"<<upcxx::myrank();
  logfile.open(filename.str().c_str());


  double io_time = mysecond();

  /* initialize random seed: */
  int seed =time(NULL); 
  srand (seed);

  // Don't need xadj any more as its info is stored in node_t
  vector<int> xadj;
  // upcxx::shared_array<int> xadj_shared;
  // xadj doesn't need to be shared because it's static and relatively small.
  
  vector<int> local_adj;
  vector<int> ixadj;
  vector<int> iadj;



  ReadAdjacencyHB(argv[1], ixadj, iadj);
  //expand to asymmetric storage
  int n = ixadj.size()-1;
  ExpandSymmetric(n,&ixadj[0],&iadj[0], xadj, local_adj);



  io_time = mysecond() - io_time;
  double init_time = mysecond();

  // vector<node_t> nodes(n);








  upcxx::shared_array<node_t> nodes(n);
  nodes.init();



  //initialize nodes
  for (int i = upcxx::myrank(); i < n; i += upcxx::ranks()) {
    upcxx::global_ptr<node_t> cur_node = &nodes[i];
    // cur_node.id = i+1;
    memberof(cur_node, id) = i+1;
    // cur_node.degree = xadj[i+1] - xadj[i] -1;
    memberof(cur_node, degree) = xadj[i+1] - xadj[i] ;
    // cur_node.label = 0;
    memberof(cur_node, label) = 0;

    //memberof(cur_node, marker) = 0;

    int adj_sz = xadj[i+1] - xadj[i];
    upcxx::global_ptr<int> adj = upcxx::allocate<int>(upcxx::myrank(), adj_sz);
    assert(adj != NULL);
    memberof(cur_node, adj) = adj;
    memberof(cur_node, adj_sz) = adj_sz;
    for (int j=0; j<adj_sz; j++) {
      adj[j] = local_adj[xadj[i] - 1 + j];
    }
  }












  // Finish reading the data from file and initializing the data structures

  upcxx::shared_array<int> all_min_degrees(upcxx::ranks());
  upcxx::shared_array<int> all_min_ids(upcxx::ranks());
  all_min_degrees.init();
  all_min_ids.init();
  
  init_time = mysecond() - init_time;
  double mdo_time = mysecond();

  vector<int> schedule;
  schedule.reserve(n);


  unsigned int tag = 0;

  // vector< upcxx::global_ptr<node_t> > schedule_shared;

  // YZ: loop-carried dependencies between steps
  //process with MD algorithm
  for (int step=1; step<=n; ++step) {
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
        if (cur_min_degree == -1 || cur_node->degree <= cur_min_degree) {
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

      upcxx::global_ptr<node_t> min_node = &nodes[all_min_ids[min_rank]-1];
      global_min_id = all_min_ids[min_rank]; // memberof(min_node, id);
      memberof(min_node, label) = step;
    }

    upcxx::barrier();
    upcxx::upcxx_bcast(&global_min_id, &global_min_id, sizeof(int), 0);
    assert(global_min_id != -1);



    if(global_min_id == n-step){
      //add remaining nodes to the schedule
    }




    //cout<<"Global min id is "<<global_min_id<<endl;

    
    schedule.push_back(global_min_id);
    //cerr << " " << schedule.back();

    //update the degree of its reachable set
    tag++;
    vector< upcxx::global_ptr<node_t> > reach;
    get_reach(nodes, global_min_id, tag, reach);

    // YZ: Use UPC++ to parallel this loop.  There is no data dependencies
    // inside the for loop because the get_reach function does not change
    // the original graph (and the adjacency list) though some nodes' degree
    // might be changed.
#ifdef USE_UPCXX
    // vector< upcxx::global_ptr<node_t> >::iterator
    for (auto it = reach.begin() + upcxx::myrank();
        it < reach.end();
        it += upcxx::ranks()) {

#else
#pragma omp parallel for
    for (vector<node_t *>::iterator it = reach.begin();
         it != reach.end();
         ++it) {
#endif
      upcxx::global_ptr<node_t> cur_neighbor = *it;
      
      vector<upcxx::global_ptr<node_t> > nghb_reach;
      
      tag++;
      get_reach(nodes, memberof(cur_neighbor, id).get(), tag, nghb_reach);

     memberof(cur_neighbor, degree) = nghb_reach.size()+1;
    }
    upcxx::barrier();


    if(upcxx::myrank()==0 && step%(n/10)==0){
      cerr<<ceil(step*100/n)<<"  ";
    }

  } // close of  for (int step=1; step<=n; ++step)

  mdo_time = mysecond() - mdo_time;

  for (int i = 0; i < upcxx::ranks(); i++) {
    if (upcxx::myrank() == i) {
      cout << "\n";
      cout<<"Rank " << upcxx::myrank() << " Schedule: ";
      for (int step=1; step<=n; ++step) {
        cout << " " << schedule[step-1];
      }
      cout << "\n";
    }
    upcxx::barrier();
  }
  if (upcxx::myrank() == 0) {
    printf("\nMinimum degree ordering algorithm time breakdown on rank 0:\n");
    printf("  io time (read graph from file into memory): %g s\n", io_time);
    printf("  setup time (initialize data structures): %g s\n", init_time);
    printf("  main algorithm time (compute the minimum degree ordering): %g s\n", mdo_time);
    printf("\n");
  }

  upcxx::barrier();
  upcxx::finalize();

  logfile.close();
  return 0;
}
  
void dump_local_nodes(node_t *local_nodes, int n)
{
  int local_size = n / upcxx::ranks();
  if (upcxx::myrank() < (n - local_size * upcxx::ranks())) {
    local_size++;
  }
  for (int i = 0; i < local_size; i++) {
    node_t *cur_node = &local_nodes[i];
    fprintf(stdout, "local_nodes[%d], id %d, degree %d, label %d, adj_sz %d, adj: ",
            i, cur_node->id, cur_node->degree, cur_node->label, cur_node->adj_sz);
    for (int j = 0; j < cur_node->adj_sz; j++) {
      fprintf(stdout, "%d ", (int)cur_node->adj[j]);
    }
    fprintf(stdout, "\n");
  }
  printf("\n");
}

