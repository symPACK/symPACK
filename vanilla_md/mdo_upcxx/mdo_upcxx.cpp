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
#include <numeric>

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




void update_degree(upcxx::shared_array<node_t> &nodes,
               const int root_node_id,
               const int tag,
               const int elim_node_id,
               const int deg_v, 
               const int tag_v, 
               vector< int> &marker,
               vector< int> &mask,
               vector< int> &perm,
               bool & indist,
               int & newDeg,
               vector<int> &work_adj
              ) 
{
  TIMER_START(UPDATE_DEGREE);
  //this list contains the nodes to explore
  stack< int > explore_set;
  explore_set.push( root_node_id );
  marker[root_node_id]=tag;

  indist = true;
  //newDeg = 1 + | reach(v) \ {u} |
  newDeg = deg_v-1;//memberof(&nodes[elim_node_id-1],degree);
  //u is in reach(v) so take it into account
  int count = 1;

#ifdef VERBOSE
  logfile<<"[ "<<root_node_id<<" ";
#endif
  //now find path between min_nodes and other nodes
  while (explore_set.size()>0) {
    //pop a node
    upcxx::global_ptr<node_t> cur_node_p = &nodes[explore_set.top()-1];
    explore_set.pop();
    node_t cur_node = *cur_node_p;

    work_adj.resize(cur_node.adj_sz);
    int * local_adj = &work_adj[0];

    //int *local_adj = new int[cur_node.adj_sz];
    upcxx::copy(cur_node.adj, upcxx::global_ptr<int>(local_adj), cur_node.adj_sz);
    for(int i = 0; i< cur_node.adj_sz; ++i){
      int curr_adj = local_adj[i];
      if(curr_adj==cur_node.id){ continue;}

      if (curr_adj != 0) {
        if (marker[curr_adj]!=tag) {
          if(perm[curr_adj-1]!=0 ){
            if(elim_node_id != curr_adj){
              explore_set.push(curr_adj);
            }
          }
          else{
            if(mask[curr_adj] == tag_v || mask[curr_adj]==INT_MAX){
#ifdef VERBOSE
//logfile<<""<<curr_adj<<","<<mask[curr_adj]<<" ";
//logfile<<""<<curr_adj<<" ";
#endif
              count++;
            }
            else{
#ifdef VERBOSE
logfile<<"["<<curr_adj<<","<<mask[curr_adj]<<"] ";
//logfile<<""<<curr_adj<<" ";
#endif
              newDeg++;
              indist=false;
            }
          }
          marker[curr_adj] = tag;
        }
      }
      else {
        fprintf(stderr, "get_reach() Fatal error 2: adj[%d]=0 \n", i);
        exit(1);
      }
    }
    //delete local_adj;
  }
#ifdef VERBOSE
  logfile<<"]"<<endl;
#endif

#ifdef VERBOSE
  if(indist){
    logfile<<"count = "<<count+1<<" vs "<<deg_v<<endl;
  }
#endif
  if(indist && count+1 != deg_v){
    indist = false;
  }

  if(indist){
#ifdef VERBOSE
    logfile<<"Node "<<root_node_id<<" is indist"<<endl;
#endif
    mask[root_node_id]=INT_MAX;
  }
  TIMER_STOP(UPDATE_DEGREE);
}





void get_reach(upcxx::shared_array<node_t> &nodes,
               const int root_node_id,
               const int tag, 
               vector< int > &reach_set,
               vector< int> &marker,
               vector< int> &perm,
               vector<int> &work_adj
              ) 
{
  //this array is used to mark the node so that 
  //they are not added twice into the reachable set

  TIMER_START(GET_REACH);
  //this list contains the nodes to explore
  stack< int > explore_set;
  reach_set.resize(0);

  explore_set.push( root_node_id );
  marker[root_node_id]=tag;
  //  memberof(root_node_p,marker)=tag;


  //now find path between min_nodes and other nodes
  while (explore_set.size()>0) {
    //pop a node
    upcxx::global_ptr<node_t> cur_node_p = &nodes[explore_set.top()-1];
  
    explore_set.pop();
    assert(cur_node_p.raw_ptr() != NULL);
    node_t cur_node = *cur_node_p;

    work_adj.resize(cur_node.adj_sz);
    int * local_adj = &work_adj[0];
    //int *local_adj = new int[cur_node.adj_sz];
    upcxx::copy(cur_node.adj, upcxx::global_ptr<int>(local_adj), cur_node.adj_sz);
    for(int i = 0; i< cur_node.adj_sz; ++i){
      int curr_adj = local_adj[i];
      if(curr_adj==cur_node.id){ continue;}

      if (curr_adj != 0) {
        if (marker[curr_adj]!=tag) {
          if(perm[curr_adj-1]==0){
            reach_set.push_back(curr_adj);
          }
          else{
            explore_set.push(curr_adj);
          }
          marker[curr_adj] = tag;
        }
      } else {
        fprintf(stderr, "get_reach() Fatal error 2: adj[%d]=0 \n", i);
        exit(1);
      }
    }
    //delete local_adj;
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
         if(marker[node_id]==tag && perm[node_id-1]!=0 ){
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




int main(int argc, char *argv[]) 
{



  upcxx::init(&argc, &argv);


  if(argc<2){
    cerr<<"Usage is: "<<argv[0]<<" input.file"<<endl;
    return -1;
  }

bool doMassElim = false;

  if(argc>2){
    //process the doMerge parameter
    string strMass(argv[2]);
    if(strMass == "YES"){
      doMassElim = true;
    }
    if(strMass == "NO"){
      doMassElim = false;
    }
   }



bool doEarlyExit = true;

  if(argc>3){
    //process the doMerge parameter
    string str(argv[3]);
    if(str == "YES"){
      doEarlyExit = true;
    }
    if(str == "NO"){
      doEarlyExit = false;
    }
   }
  bool always_merge = false;
  bool never_merge = false;
  int period = -1;
  int threshold = -1;
  int mergeAt = -1;
  if(argc>4){
    //process the doMerge parameter
    string strMerge(argv[4]);
    if(strMerge == "YES"){
      always_merge = true;
    }
    else if(strMerge=="NO"){
      never_merge = true;
    }
    else if(strMerge=="PERIODIC"){
      if(argc<6){
        cerr<<"You must specify a period"<<endl;
        return -1;
      }
      period = atoi(argv[5]);
    }
    else if(strMerge == "AFTER"){
      if(argc<6){
        cerr<<"You must specify a starting step"<<endl;
        return -1;
      }
      threshold = atoi(argv[5]);
    }
    else if(strMerge == "AT"){
      if(argc<6){
        cerr<<"You must specify a step"<<endl;
        return -1;
      }
      mergeAt = atoi(argv[5]);
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

  upcxx::shared_array<int> all_indist_count(upcxx::ranks());
  upcxx::shared_array<int> all_min_degrees(upcxx::ranks());
  upcxx::shared_array<int> all_min_ids(upcxx::ranks());
  all_indist_count.init();
  all_min_degrees.init();
  all_min_ids.init();
  
  init_time = mysecond() - init_time;

    vector< int > reach;
    vector< int > nghb_reach;
    vector< int > work_adj;
    reach.reserve(n);
    nghb_reach.reserve(n);
    work_adj.reserve(n);
  TIMER_STOP(init);


TIMER_START(mdo_time);
  double mdo_time = mysecond();

  vector<int> schedule;
  schedule.reserve(n);

  vector<int> perm(n,0);
//  schedule.reserve(n);

  int tag = 0;
  vector<int> marker(n+1,0);
  vector<int> mask(n+1,0);

  // vector< upcxx::global_ptr<node_t> > schedule_shared;

  // YZ: loop-carried dependencies between steps
  //process with MD algorithm
  int step=0;
  int label =1;

  node_t *local_nodes = (node_t *)&nodes[upcxx::myrank()];
  // YZ: need to handle the case when n is not a multiple of ranks()!!
  int other_size = n / upcxx::ranks();
  int local_size = n / upcxx::ranks();
  if (upcxx::myrank() < (n - local_size * upcxx::ranks())) {
    local_size++;
  }

  int * work_prefix = new int[upcxx::ranks()+1];

  while(label<=n){
    node_t *my_min = NULL;
    // printf("step %d, local_size %d\n", step, local_size);
    int cur_min_degree = -1;
    for (int i = 0; i < local_size; i++) {
      node_t *cur_node = &local_nodes[i];
      //if (cur_node->label == 0) {
      if (perm[cur_node->id-1] == 0) {
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
        //logfile<<all_min_degrees[i]<<endl;
        if (cur_min_degree > all_min_degrees[i]) {
          cur_min_degree = all_min_degrees[i];
          min_rank = i;
        }
      }

      upcxx::global_ptr<node_t> min_node = &nodes[all_min_ids[min_rank]-1];
      global_min_id = all_min_ids[min_rank]; // memberof(min_node, id);
      //memberof(min_node, label) = step;
      global_min_degree = all_min_degrees[min_rank];
      //logfile<<memberof(min_node,degree)<<" vs "<<global_min_degree<<endl;
      assert(memberof(min_node,id).get() == global_min_id);
      assert(memberof(min_node,degree).get() == global_min_degree);
    }

    //upcxx::barrier();
    upcxx::upcxx_bcast(&global_min_id, &global_min_id, sizeof(int), 0);
    upcxx::upcxx_bcast(&global_min_degree, &global_min_degree, sizeof(int), 0);
    assert(global_min_id != -1);

    upcxx::global_ptr<node_t> min_ptr = &nodes[global_min_id-1];
//    int global_min_degree = memberof(min_ptr,degree);
    assert(global_min_degree == memberof(min_ptr,degree));

    if(doEarlyExit && global_min_degree == n-step){
      //add remaining nodes to the schedule
      //TODO this is very inefficient

      int local_label_count = 0;
      for (int i = 0; i < local_size; i++) {
        node_t *cur_node = &local_nodes[i];
        if(perm[cur_node->id-1] == 0){
          local_label_count++;
        }
      }
      all_indist_count[upcxx::myrank()]=local_label_count;
      upcxx::barrier();
      //compute sequentially for now
      int total_label_count =0;
      if(upcxx::myrank()==0){
        for(int i=0;i<upcxx::ranks();++i){ work_prefix[i] = all_indist_count[i]; }

        for(int i=0;i<upcxx::ranks();++i){ logfile<<work_prefix[i]<<" "; } logfile<<" >> ";

        //Get prefix sum
        partial_sum(work_prefix,work_prefix+upcxx::ranks(),work_prefix);

        for(int i=0;i<upcxx::ranks();++i){ logfile<<work_prefix[i]<< " "; } logfile<<endl;

        //Shift by one
        for(int i=1;i<upcxx::ranks();++i){ all_indist_count[i] = work_prefix[i-1]; }
        all_indist_count[0] = 0 ;
        total_label_count = work_prefix[upcxx::ranks()-1];
      }
      upcxx::upcxx_bcast(&total_label_count,&total_label_count, sizeof(int), 0);

      label+=all_indist_count[upcxx::myrank()];
      //logfile<<"Labeling from "<<label<<endl;
      for (int i = 0; i < local_size; i++) {
        node_t *cur_node = &local_nodes[i];
        if(perm[cur_node->id-1] == 0){
          //logfile<<cur_node->id<<" is labeled "<<label<<endl;
          perm[cur_node->id-1]=label++;
        }
      }

#ifdef VERBOSE
      logfile<<"PERM     : "; for(int i=0;i<perm.size();++i){logfile<<perm[i]<<" ";}; logfile<<endl;
#endif
      upcxx::upcxx_reduce(&perm[0],&perm[0],perm.size(),0,UPCXX_MAX,UPCXX_INT);
      upcxx::upcxx_bcast(&perm[0], &perm[0],perm.size()*sizeof(int),0);
#ifdef VERBOSE
      logfile<<"RED PERM : "; for(int i=0;i<perm.size();++i){logfile<<perm[i]<<" ";}; logfile<<endl;
#endif
      break;
    }


#ifdef VERBOSE
    logfile<<"Node "<<global_min_id<<" is labeled "<<label<<endl;
#endif
    perm[global_min_id-1]=label++;

    //cerr << " " << schedule.back();

    //update the degree of its reachable set
    //logfile<<"Reach of root "<<global_min_id<<endl;
    advance_tag(tag);
    get_reach(nodes, global_min_id, tag, reach,marker,perm,work_adj);
    assert(reach.size() == global_min_degree-1);

    bool doMerge = (always_merge || (step%period==0) || (step>=threshold) || (step==mergeAt)) && !never_merge;
    if(doMerge){
      merge_path(nodes,global_min_id,tag,reach,marker,perm);
    }

    //all_indist_count[upcxx::myrank()] = 0;        
    int local_indist_count = 0;


    //mark the nodes in reach(v) with tag_v
    advance_tag(tag);
    int tag_v = tag;
    for (auto it = reach.begin(); it < reach.end(); ++it) { mask[*it] = tag_v; }

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
      //upcxx::global_ptr<node_t> cur_neighbor = *it;
      upcxx::global_ptr<node_t> cur_neighbor = &nodes[*it-1];
  
#ifdef VERBOSE 
      //vector<int> nghb_reach;  
      advance_tag(tag);
      get_reach(nodes, *it, tag, nghb_reach,marker,perm,work_adj);
      memberof(cur_neighbor, degree) = nghb_reach.size()+1;

    logfile<<"{ "<<global_min_id<<"| ";
    for (auto it2 = reach.begin(); it2 < reach.end(); ++it2) {
      logfile<<*it2<<" ";
    }
    logfile<<"}"<<endl;

    logfile<<"[[ "<<*it<<"| ";
    for (auto it2 = nghb_reach.begin(); it2 < nghb_reach.end(); ++it2) {
      logfile<<*it2<<" ";
    }
    logfile<<"]]"<<endl;

    for (auto it2 = reach.begin(); it2 < reach.end(); ++it2) {
      if(mask[*it2]!=INT_MAX){ mask[*it2] = tag_v; }
    }
#endif

      advance_tag(tag);
      bool indist = false;
      int new_degree =0;
      update_degree(nodes, *it, tag, global_min_id,global_min_degree, tag_v, marker,mask, perm, indist, new_degree, work_adj);
      memberof(cur_neighbor, degree) = new_degree;
//      memberof(cur_neighbor, degree) = nghb_reach.size()+1;

      if(doMassElim && indist){
//        //mark as indist
        local_indist_count++;
      }


#ifdef VERBOSE
//      logfile<<global_min_degree<<endl;
//      logfile<<new_degree<<" vs "<<nghb_reach.size()+1<<endl;
      assert(new_degree == nghb_reach.size()+1);
#endif
    }

    if(doMassElim){
      //compute prefix sum of local indist nodes counts
      all_indist_count[upcxx::myrank()]=local_indist_count;        
      upcxx::barrier();
      //compute sequentially for now
      int total_indist_count =0;
      if(upcxx::myrank()==0){


        for(int i=0;i<upcxx::ranks();++i){
          work_prefix[i] = all_indist_count[i];
        }
#ifdef VERBOSE
        for(int i=0;i<upcxx::ranks();++i){ logfile<<work_prefix[i]<<" "; } logfile<<" >> ";
#endif
        //Get prefix sum
        partial_sum(work_prefix,work_prefix+upcxx::ranks(),work_prefix);

#ifdef VERBOSE
        for(int i=0;i<upcxx::ranks();++i){ logfile<<work_prefix[i]<< " "; } logfile<<endl;
#endif

        for(int i=1;i<upcxx::ranks();++i){
          all_indist_count[i] = work_prefix[i-1];
        }
        all_indist_count[0] = 0 ;
        total_indist_count = work_prefix[upcxx::ranks()-1];


//        total_indist_count = all_indist_count[0];
//        for(int i=1;i<upcxx::ranks();++i){
//          int cur = all_indist_count[i];
//          total_indist_count+= cur;
//          all_indist_count[i] = cur + all_indist_count[i-1] - all_indist_count[0];
//        }
//        all_indist_count[0] = 0 ;
      }
      upcxx::upcxx_bcast(&total_indist_count,&total_indist_count, sizeof(int), 0);

      //now assign labels and update degrees of other nodes 
      int newlabel = label + all_indist_count[upcxx::myrank()];     
      for (auto it = reach.begin() + upcxx::myrank();
          it < reach.end();
          it += upcxx::ranks()) {
        if(mask[*it]==INT_MAX){
#if 1
          memberof(&nodes[*it-1],label) = newlabel;
#endif
          perm[*it-1] = newlabel++;
        }
        else{
          memberof(&nodes[*it-1],degree) = memberof(&nodes[*it-1],degree) - total_indist_count;
        }
      }
      label += total_indist_count;
//      logfile<<"Label is now "<<label<<endl;
    }

    if(doMassElim){
TIMER_START(REDUCE_DEG_UPDATES);
#if 1
      upcxx::barrier();
      //get global node updates
      for (auto it = reach.begin();
          it < reach.end();
          it++) {
        perm[*it-1] = memberof(&nodes[*it-1],label);
      }
      upcxx::barrier();
#else
    //reduce perm to get global nodes updates
#ifdef VERBOSE
    logfile<<"PERM     : "; for(int i=0;i<perm.size();++i){logfile<<perm[i]<<" ";}; logfile<<endl;
#endif

    upcxx::upcxx_reduce(&perm[0],&perm[0],perm.size(),0,UPCXX_MAX,UPCXX_INT);
    upcxx::upcxx_bcast(&perm[0], &perm[0], perm.size()*sizeof(int),0);

#ifdef VERBOSE
    logfile<<"RED PERM : "; for(int i=0;i<perm.size();++i){logfile<<perm[i]<<" ";}; logfile<<endl;
#endif
#endif
TIMER_STOP(REDUCE_DEG_UPDATES);
    }
    else{
      upcxx::barrier();
    }

    if(n>10){
      if(upcxx::myrank()==0 && step%(n/10)==0){
        cerr<<ceil(step*100/n)<<"  ";
      }
    }
    else{
      if(upcxx::myrank()==0){
        cerr<<step<<"  ";
      }
    }

    step++;
  } // close of  for (int step=1; step<=n; ++step)
  if(upcxx::myrank()==0){
  cerr<<endl;
  }

  //reduce the perm array


  mdo_time = mysecond() - mdo_time;
TIMER_STOP(mdo_time);

//  for (int i = 0; i < upcxx::ranks(); i++) {
    if (upcxx::myrank() == 0) {
      cout << "\n";
      cout<<"Rank " << upcxx::myrank() << " Schedule: ";
      schedule.resize(n);
      for (int node=1; node<=n; ++node) {
        int label = perm[node-1];
        schedule[label-1] = node;
      }

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
    printf("  io time (read graph from file into memory): %g s\n", io_time);
    printf("  setup time (initialize data structures): %g s\n", init_time);
    printf("  main algorithm time (compute the minimum degree ordering): %g s\n", mdo_time);
    printf("\n");
  }

  delete [] work_prefix;

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

