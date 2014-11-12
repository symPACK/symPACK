#include <vector>
#include <list>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <iostream>
#include <algorithm>

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#include "timer.h"
#include "util.h"
  

using namespace std;

#define BFS
//#define STAT
//#define ALLOC_BOOL

using namespace LIBCHOLESKY;

double mysecond()
{
  struct timeval tv; 
  gettimeofday(&tv, 0); 
  return tv.tv_sec + ((double) tv.tv_usec / 1000000);
}


struct node_stat_t{
  //stats on paths starting from THIS node
  int max_path;
  int total_path;
  int count_path;
  int max_path_step;

  //stats on paths explored when THIS node is eliminated
  int max_path_elim;
  int total_path_elim;
  int count_path_elim;
  int max_path_step_elim;

  node_stat_t(){
    max_path=0;
    total_path=0;
    count_path=0;
    max_path_step=0;
    max_path_elim=0;
    total_path_elim=0;
    count_path_elim=0;
    max_path_step_elim=0;
  }
  void fillStat(int pathLength, int step){
    if(max_path<pathLength){
      max_path = pathLength;
      max_path_step = step;
    }
    total_path += pathLength;
    count_path++;
  }
  void fillStatElim(int pathLength, int step){
    if(max_path_elim<pathLength){
      max_path_elim = pathLength;
      max_path_step_elim = step;
    }
    total_path_elim += pathLength;
    count_path_elim++;
  }
};

struct node_t{
  int id;
  int degree;
  int label;
//  int adj_sz;
  vector<int> adj;
//  int merged;
#ifndef ALLOC_BOOL
  unsigned int marker;
#endif
};


//typedef set<node_t *> deglist;
typedef list<node_t *> deglist;


#ifdef STAT
inline void get_reach(const vector<int> & xadj,vector<int> & adj,vector<node_t*> & nodes, const int root_node_id, const unsigned int tag,  vector<int> & reach_set, bool mergePath, vector<node_stat_t *> & stat, int step, const int elim_node_id )
#else
inline void get_reach(const vector<int> & xadj,vector<int> & adj,vector<node_t*> & nodes, const int root_node_id, const unsigned int tag,  vector<int> & reach_set,bool mergePath )
#endif
{
  TIMER_START(GET_REACH);
  //this array is used to mark the node so that 
  //they are not added twice into the reachable set
#ifdef ALLOC_BOOL
  TIMER_START(ALLOC_EXPLORED);
//  vector<bool> explored(nodes.size(),false);
  vector<char> explored(nodes.size(),false);
  TIMER_STOP(ALLOC_EXPLORED);
#endif

#ifdef BFS
 queue<int> explore_set;
#else
  //using a stack behaves like a Depth-First-Search
 stack<int> explore_set;
#endif

  reach_set.resize(0);


  node_t * root_node = nodes[root_node_id-1];
#ifndef ALLOC_BOOL
  root_node->marker = tag;
#endif
  explore_set.push(root_node->id);

#ifdef STAT
stack<int> path_lengths;
  path_lengths.push(0);
#endif


//  vector<int> visited_nodes;
//  visited_nodes.reserve(nodes.size());
  int visited = 0;






  //now find path between min_nodes and other nodes
  while(explore_set.size()>0){
    //pop a node
#ifdef BFS
    int cur_id = explore_set.front();
    node_t * cur_node = nodes[cur_id-1];
    //node_t * cur_node = explore_set.front();
#else
    
    int cur_id = explore_set.top();
    node_t * cur_node = nodes[cur_id-1];
    //node_t * cur_node = explore_set.top();
#endif
    explore_set.pop();

#ifdef STAT
      int cur_length = path_lengths.top();
#endif


      for(auto it = cur_node->adj.begin(); it!=cur_node->adj.end();++it){
        int node_id = *it;
          node_t * next_node = nodes[node_id-1];
          if(next_node->marker!=tag){
            if(next_node->label==0){
              //reach_set.push_back(next_node);
              reach_set.push_back(node_id);

#ifdef STAT
              stat[root_node_id-1]->fillStat(cur_length+1,step);
              stat[elim_node_id-1]->fillStatElim(cur_length+1,step);
#endif

            }
            else
            {
              explore_set.push(node_id);
#ifdef STAT
              path_lengths.push(cur_length+1);
#endif
              //visited_nodes.push_back(node_id);
              visited++;

            }
            next_node->marker=tag;
          }
      }

  }


  //now the reachable set is complete, and the visited node list as well
  if(mergePath && visited>1 && root_node->label!=0/**/){
    TIMER_START(COMPRESS_PATH);
    int new_label = -1;
    if(root_node->label!=0){
      new_label = root_node->id;
    }


    for(auto it = reach_set.begin();it!=reach_set.end();++it)
    {
      node_t * cur_node = nodes[*it-1];
      //node_t * cur_node = *it;
        int old_sz = cur_node->adj.size();
        int labeled = 0;
        //find last
        int last = cur_node->adj.size();
        for(int i = 0;i<last;++i){ 
          int node_id = cur_node->adj[i];
            node_t * adj_node = nodes[node_id-1];
            if(adj_node->marker==tag && adj_node->label!=0 ){
              if(!labeled){
                if(new_label==-1){
                  new_label = node_id;
                }
                //nodes[node_id-1]->merged = -new_label;
                cur_node->adj[i] = new_label;
                labeled=1;
              }
              else{
                //nodes[cur_node->adj[i]-1]->merged = -new_label;
                cur_node->adj[i]=cur_node->adj[last-1];
//                cur_node->adj_sz--;
                --i;
                --last;
                
              }
            }          
        }

//        assert(cur_node->adj_sz>0);
//        cout<<"adj_sz of node "<<cur_node->id<<" is now "<<cur_node->adj_sz<<" vs "<<old_sz<<endl;
        cur_node->adj.resize(last);


//      cout<<"B Node -----"<<cur_node->id<<endl;
//        for(auto it2 = cur_node->adj.begin(); it2!=cur_node->adj.end();++it2){
//          cout<<*it2<<" "; 
//        }
//      cout<<endl;
}

if(new_label!=-1){
    node_t * merged_node = nodes[new_label-1];
        //go through their adjacency list
    merged_node->adj.resize(1+reach_set.size());
    merged_node->adj[0]=merged_node->id;
    std::copy(reach_set.begin(),reach_set.end(),&(merged_node->adj[1])); 
}

    TIMER_STOP(COMPRESS_PATH);
  }


  TIMER_STOP(GET_REACH);
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


  ReadAdjacencyHB(argv[1], ixadj, iadj);

  int n = ixadj.size()-1;


  //expand to asymmetric storage
  vector<int> xadj;
  vector<int> adj;
  ExpandSymmetric(n,&ixadj[0],&iadj[0], xadj, adj);

  vector<node_t *> nodes(n);


#ifdef STAT
  vector<node_stat_t *> stat(n);
  for(auto it=stat.begin();it!=stat.end();it++){
    *it = new node_stat_t();
  }
#endif

  unsigned int tag = 0;


  vector<deglist * > degrees(n+1,NULL);
  int min_degree = n+1;


  //initialize nodes
  for(int i=0;i<n;++i){
    nodes[i] = new node_t;
    node_t & cur_node = *nodes[i];
    cur_node.id = i+1;
    cur_node.degree = xadj[i+1] - xadj[i] -1;

//    cur_node.adj_sz = xadj[i+1] - xadj[i];
    //coppy the adj structure
    cur_node.adj.insert(cur_node.adj.begin(),&adj[xadj[i]-1],&adj[xadj[i+1]-1]);

    cur_node.label = 0;
#ifndef ALLOC_BOOL
    cur_node.marker = 0;
#endif
    //cur_node.merged = cur_node.id;

     min_degree = min(min_degree,cur_node.degree);
    if(degrees[cur_node.degree]==NULL){
      degrees[cur_node.degree] = new deglist;
    }
    //degrees[cur_node.degree]->insert(&cur_node);
    degrees[cur_node.degree]->push_front(&cur_node);

  }


adj.swap(std::vector<int>());
xadj.swap(std::vector<int>());

vector<int> nghb_reach;
vector<int> reach;
//nghb_reach.reserve(n);
//reach.reserve(n);


  double mdtime = mysecond();

  TIMER_START(MD_ORDERING);
  //process with MD algorithm
  for(int step = 1; step<=n;++step){
   
  TIMER_START(FIND_MIN); 

    for(min_degree;min_degree<degrees.size();++min_degree){
      if(degrees[min_degree]!=NULL){
        if(!degrees[min_degree]->empty()){
          break;
        }
      }
    }
//    auto pos = std::find_if(degrees.begin()+min_degree,degrees.end(),[](set<node_t *> * i) -> bool{ bool retval = false; if(i!=NULL){retval= !i->empty();} return retval;});
//    min_degree = pos-degrees.begin();
  TIMER_STOP(FIND_MIN);
 
    //check for elimination if all degrees are the same
    if(min_degree==n-step){
      for(auto it = degrees[min_degree]->begin();it!=degrees[min_degree]->end();++it){
        (*it)->label=step++;
      }
      break;
    }


    //node_t & min_node = *(*degrees[min_degree]->begin());
    node_t & min_node = *(degrees[min_degree]->front());
    min_node.label = step;
//    assert(min_node.degree==min_degree);
    

TIMER_START(UPDATE_DEGREE_LIST);
    //degrees[min_node.degree]->erase(&min_node);
    degrees[min_node.degree]->pop_front();
    //auto pos = std::find(degrees[min_node.degree]->begin(),degrees[min_node.degree]->end(),&min_node);
    //degrees[min_node.degree]->erase(pos);
    
TIMER_STOP(UPDATE_DEGREE_LIST);

    //update the degree of its reachable set
    ++tag;
//    list<node_t *> reach;

    //bool doMerge = 100*step/n >=60;
//    bool doMerge = ((100*step/n)%10 == 0 ) || (100*step/n >=60);
    bool doMerge = true;

#ifdef STAT
    get_reach(xadj,adj,nodes,min_node.id,tag,reach,doMerge,stat,step,min_node.id);
#else
    get_reach(xadj,adj,nodes,min_node.id,tag,reach,doMerge);
#endif
    
    min_node.degree = reach.size();
//    assert(min_node.degree>0);



TIMER_START(REACH_UPDATE);
    for(auto it = reach.begin();it!=reach.end();++it){
      node_t * cur_neighbor = nodes[*it-1];
        ++tag;
        bool doMerge = false;
#ifdef STAT
        get_reach(xadj,adj,nodes,cur_neighbor->id,tag,nghb_reach,doMerge,stat,step,min_node.id);
#else
        get_reach(xadj,adj,nodes,cur_neighbor->id,tag,nghb_reach,doMerge);
#endif

        int old_degree = cur_neighbor->degree;
        cur_neighbor->degree = nghb_reach.size();
TIMER_START(UPDATE_DEGREE_LIST);
   //update the linked lists 
//  degrees[old_degree]->erase(cur_neighbor);
    auto pos = std::find(degrees[old_degree]->begin(),degrees[old_degree]->end(),cur_neighbor);
    degrees[old_degree]->erase(pos);

   min_degree = min(min_degree,cur_neighbor->degree);
  
    if(degrees[cur_neighbor->degree]==NULL){
      degrees[cur_neighbor->degree] = new deglist;
    }
    //degrees[cur_neighbor->degree]->insert(cur_neighbor);
    degrees[cur_neighbor->degree]->push_front(cur_neighbor);

TIMER_STOP(UPDATE_DEGREE_LIST);
    }
TIMER_STOP(REACH_UPDATE);

    if(step%(n/10)==0){
      cerr<<ceil(step*100/n)<<"  ";  
    }

  }
  cerr<<endl;
  TIMER_STOP(MD_ORDERING);

  mdtime = mysecond() - mdtime;

  //compute invperm
  vector<int> iperm(n);
  for(int i =1; i<=n;++i){
    int label = nodes[i-1]->label;
 //   assert(label>0 && label<=n);
    iperm[label-1]=i;
  }
  
  cout<<"Schedule: ";
  for(int step = 1; step<=n;++step){
    cout<<" "<<iperm[step-1];
  }
  cout<<endl;

for(int i =1; i<=n;++i){
    int label = nodes[i-1]->label;
    if(!(label>0 && label<=n)){
      for(auto it = degrees.begin();it!=degrees.end();++it){
        if(*it!=NULL){
          cout<<it-degrees.begin()<<": ";
          for(auto it2 = (*it)->begin();it2!=(*it)->end();++it2){
            cout<<(*it2)->id<<" ";
          }
          cout<<endl;
        }
      }

      abort();
    }
  }
  



 printf("  main algorithm time (compute the minimum degree ordering): %g s\n", mdtime);

#ifdef STAT
    cout<<"label"<<"\t"<<"node"<<"\t"<<"max path length"<<"\t"<<"max path step"<<"\t"<<"path count"<<"\t"<<"avg"<<endl;
  for(int step = 1; step<=n;++step){
    int node_id = iperm[step-1];
    node_stat_t * statptr = stat[node_id-1];

    double avg = 0;
    if((statptr)->count_path>0){
      avg = (statptr)->total_path/(double)(statptr)->count_path;
    }

    cout<<step<<"\t"<<node_id<<"\t"<<(statptr)->max_path<<"\t"<<(statptr)->max_path_step<<"\t"<<(statptr)->count_path<<"\t"<<avg<<endl;

//    cout<<"----- node "<<node_id<<"-----"<<endl;
//    cout<<"max path length: "<<(statptr)->max_path<<endl;
//    cout<<"max path step: "<<(statptr)->max_path_step<<endl;
//    cout<<"path count: "<<(statptr)->count_path<<endl;
//    cout<<"avg path length: "<<avg<<endl;
  }

  cout<<endl<<endl;

  cout<<"ELIMlabel"<<"\t"<<"node"<<"\t"<<"max path length"<<"\t"<<"max path step"<<"\t"<<"path count"<<"\t"<<"avg"<<endl;
  for(int step = 1; step<=n;++step){
    int node_id = iperm[step-1];
    node_stat_t * statptr = stat[node_id-1];

    double avg = 0;
    if((statptr)->count_path_elim>0){
      avg = (statptr)->total_path_elim/(double)(statptr)->count_path_elim;
    }

    cout<<step<<"\t"<<node_id<<"\t"<<(statptr)->max_path_elim<<"\t"<<(statptr)->max_path_step_elim<<"\t"<<(statptr)->count_path_elim<<"\t"<<avg<<endl;

//    cout<<"----- node "<<node_id<<"-----"<<endl;
//    cout<<"max path length: "<<(statptr)->max_path<<endl;
//    cout<<"max path step: "<<(statptr)->max_path_step<<endl;
//    cout<<"path count: "<<(statptr)->count_path<<endl;
//    cout<<"avg path length: "<<avg<<endl;
  }


    for(auto it = degrees.begin();it!=degrees.end();++it){
      if(*it!=NULL){
        delete *it;
      }
    }



  for(auto it=stat.begin();it!=stat.end();it++){
    delete *it;
  }
#endif

  for(auto it=nodes.begin();it!=nodes.end();it++){
    delete *it;
  }
}




