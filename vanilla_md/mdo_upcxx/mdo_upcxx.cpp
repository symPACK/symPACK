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

#define advance_tag(tag) do{tag++; if(tag==INT_MAX){tag=0;}}while(0)

using namespace std;

double mysecond()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + ((double) tv.tv_usec / 1000000);
}

void dump_local_nodes(node_t *local_nodes, int n);


void get_reach2(upcxx::shared_array<node_t> &nodes,
               const int root_node_id,
               const int tag, 
               vector< upcxx::global_ptr<node_t> > &reach_set,
               vector< int> &marker
              ) 
{
  //this array is used to mark the node so that 
  //they are not added twice into the reachable set
//  vector<bool> explored(nodes.size(), false);

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
  
  marker[root_node_id-1]=tag;

  //now find path between min_nodes and other nodes
  while (explore_set.size()>0) {
    //pop a node
    upcxx::global_ptr<node_t> cur_node_p = explore_set.top();
    explore_set.pop();
    assert(cur_node_p.raw_ptr() != NULL);
    node_t cur_node = *cur_node_p;

      int *local_adj = new int[cur_node.adj_sz];
      upcxx::copy(cur_node.adj, upcxx::global_ptr<int>(local_adj), cur_node.adj_sz);
      for(int i = 0; i< cur_node.adj_sz; ++i){
        int curr_adj = local_adj[i];
        if(curr_adj==cur_node.id){ continue;}

        if (curr_adj != 0) {
//          if (next_node.marker != tag) {
          if (marker[curr_adj-1]!=tag) {
            upcxx::global_ptr<node_t> next_node_p = &nodes[curr_adj-1];
            if(memberof(next_node_p,label).get()==0){
              reach_set.push_back(next_node_p);
            }
            else{
              explore_set.push(next_node_p);
            }
//            memberof(next_node_p,marker)=tag;
            marker[curr_adj-1] = tag;
          }
        } else {
          fprintf(stderr, "get_reach() Fatal error 2: adj[%d]=0 \n", i);
          exit(1);
        }
      }
      delete local_adj;



  }
}


void get_reach3(upcxx::shared_array<node_t> &nodes,
               const int root_node_id,
               const int tag, 
               vector< int > &reach_set,
               vector< int> &marker,
               vector< int> &perm
              ) 
{
  //this array is used to mark the node so that 
  //they are not added twice into the reachable set

  TIMER_START(GET_REACH);
  //this list contains the nodes to explore
  stack< int > explore_set;
  reach_set.resize(0);

  explore_set.push( root_node_id );
  marker[root_node_id-1]=tag;
  //  memberof(root_node_p,marker)=tag;


  //now find path between min_nodes and other nodes
  while (explore_set.size()>0) {
    //pop a node
//  TIMER_START(GET_EXPLORE_NODE);
    upcxx::global_ptr<node_t> cur_node_p = &nodes[explore_set.top()-1];
//  TIMER_STOP(GET_EXPLORE_NODE);
  
    explore_set.pop();
    assert(cur_node_p.raw_ptr() != NULL);
    node_t cur_node = *cur_node_p;

    int *local_adj = new int[cur_node.adj_sz];
//  TIMER_START(GET_ADJ_SET);
    upcxx::copy(cur_node.adj, upcxx::global_ptr<int>(local_adj), cur_node.adj_sz);
//  TIMER_STOP(GET_ADJ_SET);
    for(int i = 0; i< cur_node.adj_sz; ++i){
      int curr_adj = local_adj[i];
      if(curr_adj==cur_node.id){ continue;}

      if (curr_adj != 0) {
        if (marker[curr_adj-1]!=tag) {
  //TIMER_START(GET_ADJ_NODE);
//          upcxx::global_ptr<node_t> next_node_p = &nodes[curr_adj-1];
  //TIMER_STOP(GET_ADJ_NODE);
          if(perm[curr_adj-1]==0){
            reach_set.push_back(curr_adj);
          }
          else{
            explore_set.push(curr_adj);
          }
          marker[curr_adj-1] = tag;
        }
      } else {
        fprintf(stderr, "get_reach() Fatal error 2: adj[%d]=0 \n", i);
        exit(1);
      }
    }
    delete local_adj;
  }

  TIMER_STOP(GET_REACH);
}



void get_reach(upcxx::shared_array<node_t> &nodes,
               const int root_node_id,
               const int tag, 
               vector< int > &reach_set,
               vector< int> &marker
              ) 
{
  //this array is used to mark the node so that 
  //they are not added twice into the reachable set

  TIMER_START(GET_REACH);
  //this list contains the nodes to explore
  stack< int > explore_set;
  reach_set.resize(0);

  explore_set.push( root_node_id );
  marker[root_node_id-1]=tag;
  //  memberof(root_node_p,marker)=tag;


  //now find path between min_nodes and other nodes
  while (explore_set.size()>0) {
    //pop a node
//  TIMER_START(GET_EXPLORE_NODE);
    upcxx::global_ptr<node_t> cur_node_p = &nodes[explore_set.top()-1];
//  TIMER_STOP(GET_EXPLORE_NODE);
  
    explore_set.pop();
    assert(cur_node_p.raw_ptr() != NULL);
    node_t cur_node = *cur_node_p;

    int *local_adj = new int[cur_node.adj_sz];
//  TIMER_START(GET_ADJ_SET);
    upcxx::copy(cur_node.adj, upcxx::global_ptr<int>(local_adj), cur_node.adj_sz);
//  TIMER_STOP(GET_ADJ_SET);
    for(int i = 0; i< cur_node.adj_sz; ++i){
      int curr_adj = local_adj[i];
      if(curr_adj==cur_node.id){ continue;}

      if (curr_adj != 0) {
        if (marker[curr_adj-1]!=tag) {
  //TIMER_START(GET_ADJ_NODE);
          upcxx::global_ptr<node_t> next_node_p = &nodes[curr_adj-1];
  //TIMER_STOP(GET_ADJ_NODE);
          if(memberof(next_node_p,label).get()==0){
            reach_set.push_back(curr_adj);
          }
          else{
            explore_set.push(curr_adj);
          }
          marker[curr_adj-1] = tag;
        }
      } else {
        fprintf(stderr, "get_reach() Fatal error 2: adj[%d]=0 \n", i);
        exit(1);
      }
    }
    delete local_adj;
  }

  TIMER_STOP(GET_REACH);
}


void merge_path(upcxx::shared_array<node_t> &nodes,
                const int root_node_id, 
                const int tag,  
                vector< int > &reach_set,
                vector< int > &marker,
                vector<int> & perm){
  
  //now the reachable set is complete, and the visited node list as well
    TIMER_START(COMPRESS_PATH);
    upcxx::global_ptr<node_t> root_node_ptr = &nodes[root_node_id-1];
   if(upcxx::myrank()==root_node_ptr.where()){
      //node_t root_node = *root_node_ptr; 
      int new_label = -1;
      if(perm[root_node_id-1]!=0){
        new_label = root_node_id;
      }


      vector<int> local_adj;
//      for(auto it = reach_set.begin()+upcxx::myrank();it<reach_set.end();it+=upcxx::ranks())
      for(auto it = reach_set.begin();it!=reach_set.end();++it)
      {

        if(perm[*it-1]!=0){
          continue;
        }

        upcxx::global_ptr<node_t> cur_node_ptr = &nodes[*it-1];
        node_t cur_node = *cur_node_ptr;
        
        //node_t * cur_node = *it;
        int old_sz = cur_node.adj_sz;

        //int * local_adj = new int[old_sz];
        local_adj.resize(old_sz);
        upcxx::copy(cur_node.adj,upcxx::global_ptr<int>(&local_adj[0]),old_sz);

        int labeled = 0;
        //find last
        int last = old_sz;
        for(int i = 0;i<last;++i){ 
          int node_id = local_adj[i];
          //upcxx::global_ptr<node_t> adj_node = &nodes[node_id-1];
          if(marker[node_id-1]==tag && perm[node_id-1]!=0 ){
            if(!labeled){
              if(new_label==-1){
                new_label = node_id;
              }
              //nodes[node_id-1]->merged = -new_label;
              local_adj[i] = new_label;
              labeled=1;
            }
            else{
              //nodes[cur_node->adj[i]-1]->merged = -new_label;
              local_adj[i]=local_adj[last-1];
              --i;
              --last;
            }
          }
        }
        memberof(cur_node_ptr,adj_sz)=last;
        upcxx::copy(upcxx::global_ptr<int>(&local_adj[0]),cur_node.adj,last);
//        delete [] local_adj;
      }

      if(new_label!=-1){
        upcxx::global_ptr<node_t> merged_node_ptr = &nodes[new_label-1];
        node_t merged_node = *merged_node_ptr;
        //go through their adjacency list
        int new_adj_sz = reach_set.size()+1;
        upcxx::global_ptr<int> new_adj = merged_node.adj;
        if(new_adj_sz>merged_node.adj_sz){
          new_adj = upcxx::allocate<int>(merged_node_ptr.where(),new_adj_sz);
          upcxx::deallocate(merged_node.adj);
          memberof(merged_node_ptr,adj) = new_adj;
        }

        memberof(merged_node_ptr,adj_sz)=new_adj_sz;

        int * new_local_adj = new int[new_adj_sz];
        new_local_adj[0] = merged_node.id;
        for(int i=1;i<new_adj_sz;++i){
          //new_local_adj[i] = memberof(reach_set[i-1],id);
          new_local_adj[i] = reach_set[i-1];
        }
        upcxx::copy(upcxx::global_ptr<int>(new_local_adj),new_adj,new_adj_sz);
  
        delete [] new_local_adj;
      }
    }
    upcxx::barrier();
    TIMER_STOP(COMPRESS_PATH);
}

   void merge_path3(upcxx::shared_array<node_t> &nodes,
                const int root_node_id, 
                const int tag,  
                vector< int > &reach_set,
                vector< int > &marker,
                vector<int> & perm){

     //now the reachable set is complete, and the visited node list as well
     TIMER_START(COMPRESS_PATH);
     upcxx::global_ptr<node_t> root_node_ptr = &nodes[root_node_id-1];
     //node_t root_node = *root_node_ptr; 
     int new_label = -1;
     if(perm[root_node_id-1]!=0){
       new_label = root_node_id;
     }


     vector<int> local_adj;
     for(auto it = reach_set.begin()+upcxx::myrank();it<reach_set.end();it+=upcxx::ranks())
       //      for(auto it = reach_set.begin();it!=reach_set.end();++it)
     {

       if(perm[*it-1]!=0){
         continue;
       }

       upcxx::global_ptr<node_t> cur_node_ptr = &nodes[*it-1];
       node_t cur_node = *cur_node_ptr;

       //node_t * cur_node = *it;
       int old_sz = cur_node.adj_sz;

       //int * local_adj = new int[old_sz];
       local_adj.resize(old_sz);
       upcxx::copy(cur_node.adj,upcxx::global_ptr<int>(&local_adj[0]),old_sz);

       int labeled = 0;
       //find last
       int last = old_sz;
       for(int i = 0;i<last;++i){ 
         int node_id = local_adj[i];
         //upcxx::global_ptr<node_t> adj_node = &nodes[node_id-1];
         if(marker[node_id-1]==tag && perm[node_id-1]!=0 ){
           if(!labeled){
             if(new_label==-1){
               new_label = node_id;
             }
             //nodes[node_id-1]->merged = -new_label;
             local_adj[i] = new_label;
             labeled=1;
           }
           else{
             //nodes[cur_node->adj[i]-1]->merged = -new_label;
             local_adj[i]=local_adj[last-1];
             --i;
             --last;
           }
         }
       }
       memberof(cur_node_ptr,adj_sz)=last;

       //        logfile<<"New Adj["<<*it<<"]: ";
       //        for(int i=0;i<last;++i){ logfile<<local_adj[i]<<" "; }
       //        logfile<<endl;

       upcxx::copy(upcxx::global_ptr<int>(&local_adj[0]),cur_node.adj,last);
       //        delete [] local_adj;
     }
     upcxx::barrier();

     if(upcxx::myrank()==root_node_ptr.where()){
       if(new_label!=-1){
         upcxx::global_ptr<node_t> merged_node_ptr = &nodes[new_label-1];
         node_t merged_node = *merged_node_ptr;
         //go through their adjacency list
         int new_adj_sz = reach_set.size()+1;
         upcxx::global_ptr<int> new_adj = merged_node.adj;
         if(new_adj_sz>merged_node.adj_sz){
           new_adj = upcxx::allocate<int>(merged_node_ptr.where(),new_adj_sz);
           upcxx::deallocate(merged_node.adj);
           memberof(merged_node_ptr,adj) = new_adj;
         }

         memberof(merged_node_ptr,adj_sz)=new_adj_sz;

         int * new_local_adj = new int[new_adj_sz];
         new_local_adj[0] = merged_node.id;
         for(int i=1;i<new_adj_sz;++i){
           //new_local_adj[i] = memberof(reach_set[i-1],id);
           new_local_adj[i] = reach_set[i-1];
         }
         upcxx::copy(upcxx::global_ptr<int>(new_local_adj),new_adj,new_adj_sz);

         delete [] new_local_adj;
       }
     }
     upcxx::barrier();
     TIMER_STOP(COMPRESS_PATH);
   }



void merge_path2(upcxx::shared_array<node_t> &nodes,
                const int root_node_id, 
                const int tag,  
                vector< upcxx::global_ptr<node_t> > &reach_set,
                vector<int> & marker){
  
  //now the reachable set is complete, and the visited node list as well
    TIMER_START(COMPRESS_PATH);
    upcxx::global_ptr<node_t> root_node_ptr = &nodes[root_node_id-1];
    if(upcxx::myrank()==root_node_ptr.where()){
      node_t root_node = *root_node_ptr; 
      int new_label = -1;
      if(root_node.label!=0){
        new_label = root_node.id;
      }


      for(auto it = reach_set.begin();it!=reach_set.end();++it)
      {
        upcxx::global_ptr<node_t> cur_node_ptr = *it;

        node_t cur_node = *cur_node_ptr;
        if(cur_node.label!=0){
          continue;
        }

        
        //node_t * cur_node = *it;
        int old_sz = cur_node.adj_sz;

        int * local_adj = new int[old_sz];
        upcxx::copy(cur_node.adj,upcxx::global_ptr<int>(local_adj),old_sz);

        int labeled = 0;
        //find last
        int last = old_sz;
        for(int i = 0;i<last;++i){ 
          int node_id = local_adj[i];
          upcxx::global_ptr<node_t> adj_node = &nodes[node_id-1];
          if(marker[node_id-1]==tag && memberof(adj_node,label)!=0 ){
            if(!labeled){
              if(new_label==-1){
                new_label = node_id;
              }
              //nodes[node_id-1]->merged = -new_label;
              local_adj[i] = new_label;
              labeled=1;
            }
            else{
              //nodes[cur_node->adj[i]-1]->merged = -new_label;
              local_adj[i]=local_adj[last-1];
              --i;
              --last;
            }
          }
        }
        memberof(cur_node_ptr,adj_sz)=last;
        upcxx::copy(upcxx::global_ptr<int>(local_adj),cur_node.adj,last);
        delete [] local_adj;
      }

      if(new_label!=-1){
        upcxx::global_ptr<node_t> merged_node_ptr = &nodes[new_label-1];
        node_t merged_node = *merged_node_ptr;
        //go through their adjacency list
        int new_adj_sz = reach_set.size()+1;
        upcxx::global_ptr<int> new_adj = merged_node.adj;
        if(new_adj_sz>merged_node.adj_sz){
          new_adj = upcxx::allocate<int>(merged_node_ptr.where(),new_adj_sz);
          upcxx::deallocate(merged_node.adj);
          memberof(merged_node_ptr,adj) = new_adj;
        }

        memberof(merged_node_ptr,adj_sz)=new_adj_sz;

        int * new_local_adj = new int[new_adj_sz];
        new_local_adj[0] = merged_node.id;
        for(int i=1;i<new_adj_sz;++i){
          new_local_adj[i] = memberof(reach_set[i-1],id);
        }
        upcxx::copy(upcxx::global_ptr<int>(new_local_adj),new_adj,new_adj_sz);
  
        delete [] new_local_adj;
      }
    }
    upcxx::barrier();
    TIMER_STOP(COMPRESS_PATH);
}




int main(int argc, char *argv[]) 
{



  upcxx::init(&argc, &argv);


  if(argc<2){
    cerr<<"Usage is: "<<argv[0]<<" input.file"<<endl;
    return -1;
  }


  bool always_merge = false;
  bool never_merge = false;
  int period = -1;
  int threshold = -1;
  int mergeAt = -1;
  if(argc>2){
    //process the doMerge parameter
    string strMerge(argv[2]);
    if(strMerge == "YES"){
      always_merge = true;
    }
    else if(strMerge=="NO"){
      never_merge = true;
    }
    else if(strMerge=="PERIODIC"){
      if(argc<4){
        cerr<<"You must specify a period"<<endl;
        return -1;
      }
      period = atoi(argv[3]);
    }
    else if(strMerge == "AFTER"){
      if(argc<4){
        cerr<<"You must specify a starting step"<<endl;
        return -1;
      }
      threshold = atoi(argv[3]);
    }
    else if(strMerge == "AT"){
      if(argc<4){
        cerr<<"You must specify a step"<<endl;
        return -1;
      }
      mergeAt = atoi(argv[3]);
    }
  }


  double io_time = mysecond();


  stringstream filename;
  filename<<"logTest"<<upcxx::myrank();
  logfile.open(filename.str().c_str());


  TIMER_START(io);

  int n,nnz;
  upcxx::shared_array<node_t> * nodes_ptr;
  ReadAdjacencyHB_PARA(argv[1], nodes_ptr, n, nnz);

  io_time = mysecond() - io_time;
  upcxx::shared_array<node_t> & nodes = *nodes_ptr;

  logfile<<"Input file read"<<endl;

  TIMER_STOP(io);

  TIMER_START(init);
    upcxx::barrier(); 
  ExpandSymmetric_PARA(n, nodes);

  logfile<<"Matrix expanded"<<endl;
  double init_time = mysecond();

  if(threshold==-1){threshold=n+1;}
  if(period==-1){period=n+1;}




  // Finish reading the data from file and initializing the data structures
  // dump_local_nodes((node_t *)&nodes[upcxx::myrank()], nodes.size());

  upcxx::shared_array<int> all_min_degrees(upcxx::ranks());
  upcxx::shared_array<int> all_min_ids(upcxx::ranks());
  all_min_degrees.init();
  all_min_ids.init();
  
  init_time = mysecond() - init_time;

    vector< int > reach;
    vector< int > nghb_reach;
    reach.reserve(n);
    nghb_reach.reserve(n);
  TIMER_STOP(init);




TIMER_START(mdo_time);
  double mdo_time = mysecond();

  vector<int> schedule;
  schedule.reserve(n);

  vector<int> perm(n,0);
//  schedule.reserve(n);

  int tag = 0;
  vector<int> marker(n+1);

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
      memberof(min_node, label) = step;
      global_min_degree = cur_min_degree;
    }

    //upcxx::barrier();
    upcxx::upcxx_bcast(&global_min_id, &global_min_id, sizeof(int), 0);
    upcxx::upcxx_bcast(&global_min_degree, &global_min_degree, sizeof(int), 0);
    assert(global_min_id != -1);

//    upcxx::global_ptr<node_t> min_ptr = &nodes[global_min_id-1];

    if(global_min_degree == n-step){
//    if(memberof(min_ptr,degree) == n-step){
      //add remaining nodes to the schedule
      //TODO this is very inefficient
      for(int i = 0;i<n;i++){
        upcxx::global_ptr<node_t> cur_ptr = &nodes[i];
        if(memberof(cur_ptr,label)==0){
          schedule.push_back(i+1);
        }
      }
      break;
    }


    //cout<<"Global min id is "<<global_min_id<<endl;

    
    schedule.push_back(global_min_id);
    perm[global_min_id-1]=schedule.size();

    //cerr << " " << schedule.back();

    //update the degree of its reachable set
    //logfile<<"Reach of root "<<global_min_id<<endl;
    advance_tag(tag);
    get_reach3(nodes, global_min_id, tag, reach,marker,perm);
//    vector<upcxx::global_ptr<node_t>> reach2;
//    get_reach2(nodes, global_min_id, tag, reach2,marker);
//    reach.resize(0);
//    for(int i =0;i<reach2.size();++i){ reach.push_back(memberof(reach2[i],id)); }


    bool doMerge = (always_merge || (step%period==0) || (step>=threshold) || (step==mergeAt)) && !never_merge;
    if(doMerge){
      //merge_path2(nodes,global_min_id,tag,reach2,marker);
      //merge_path(nodes,global_min_id,tag,reach,marker);
      merge_path3(nodes,global_min_id,tag,reach,marker,perm);
//      merge_path(nodes,global_min_id,tag,reach,marker,perm);
    }





    // YZ: Use UPC++ to parallel this loop.  There is no data dependencies
    // inside the for loop because the get_reach function does not change
    // the original graph (and the adjacency list) though some nodes' degree
    // might be changed.
#ifdef USE_UPCXX
    // vector< upcxx::global_ptr<node_t> >::iterator
//    for (auto it = reach2.begin() + upcxx::myrank();
//        it < reach2.end();
//        it += upcxx::ranks()) {

    for (auto it = reach.begin() + upcxx::myrank();
        it < reach.end();
        it += upcxx::ranks()) {
//    for (auto it = reach.begin();
//        it != reach.end();
//        it++) {
//        if((*it-1)%upcxx::ranks()!=upcxx::myrank()){continue;}
#else
#pragma omp parallel for
    for (vector<node_t *>::iterator it = reach.begin();
         it != reach.end();
         ++it) {
#endif
      //upcxx::global_ptr<node_t> cur_neighbor = *it;
      upcxx::global_ptr<node_t> cur_neighbor = &nodes[*it-1];
    
      //vector<int> nghb_reach;  
      advance_tag(tag);
      get_reach3(nodes, *it, tag, nghb_reach,marker,perm);
//    vector<upcxx::global_ptr<node_t>> nghb_reach2;
//    get_reach2(nodes, memberof(cur_neighbor,id).get(), tag, nghb_reach2,marker);
    //nghb_reach.resize(0);
    //for(int i =0;i<nghb_reach2.size();++i){ nghb_reach.push_back(memberof(nghb_reach2[i],id));} 

     //memberof(cur_neighbor, degree) = nghb_reach2.size()+1;
     memberof(cur_neighbor, degree) = nghb_reach.size()+1;
    }
    //updated nodes are now local so barrier is probably not necessary
    upcxx::barrier();


    if(upcxx::myrank()==0 && step%(n/10)==0){
      cerr<<ceil(step*100/n)<<"  ";
    }

  } // close of  for (int step=1; step<=n; ++step)
  if(upcxx::myrank()==0){
  cerr<<endl;
  }

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

