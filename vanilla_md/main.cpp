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

#include "ETree.hpp"

//#define verbose

using namespace std;

void displayMatrix(vector<int> & xadj, vector<int> & adj);
  void GetLColRowCount(ETree & tree,const int * xadj, const int * adj, vector<int> & cc, vector<int> & rc);

struct node_t{
  int id;
  int degree;
  int elim_step;
};

bool node_comp(node_t * & a, node_t * & b){
      if (a->degree<b->degree)
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



void get_reach(const vector<int> & xadj,const vector<int> & adj,const vector<node_t> & nodes, const node_t & min_node, const int elim_step, list<node_t*> & reach_set){
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
    explore_set.push_back(const_cast<node_t*>(&nodes[adj[i-1]-1]));
    explored[adj[i-1]-1]=1;
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
        if(!explored[adj[i-1]-1]){
          explore_set.push_back(const_cast<node_t*>(&nodes[adj[i-1]-1]));
          explored[adj[i-1]-1]=1;
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
  vector<node_t> nodes(n);
  vector<node_t *> node_heap(n);

  int initEdgeCnt = 0;
  int edgeCnt = 0;

  //initialize nodes
  for(int i=0;i<n;++i){
    node_t & cur_node = nodes[i];
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
    std::sort_heap (node_heap.begin(),node_heap.end(),node_comp);
    //get the min degree node
    std::pop_heap (node_heap.begin(),node_heap.end(),node_comp);
    node_t & min_node = *node_heap.back();
    node_heap.pop_back();
    schedule.push_back(&min_node);

    perm.push_back(min_node.id);
    invperm[min_node.id-1] = perm.size();

    min_node.elim_step = step;
#ifdef verbose
    cout<<"Node "<<min_node.id<<" scheduled at step "<<step<<endl;
#endif
    //update the degree of its reachable set
    list<node_t *> reach;
    get_reach(xadj,adj,nodes,min_node,step,reach);

    edgeCnt += reach.size();

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

  cout<<"perm: ";
  for(int step = 1; step<=n;++step){
    cout<<" "<<perm[step-1];
  }
  cout<<endl;

  cout<<"invperm: ";
  for(int step = 1; step<=n;++step){
    cout<<" "<<invperm[step-1];
  }
  cout<<endl;

  cout<<"Initial edge count: "<<initEdgeCnt<<endl;
  cout<<"Final edge count: "<<edgeCnt<<endl;

  //sort the adjacency structure according to the new labelling
  vector<int> newxadj(n+1);
  vector<int> newadj(adj.size());

  newxadj[0] = 1;
  for(int step = 1; step<=n;++step){
    newxadj[step] = newxadj[step-1] + xadj[perm[step-1]]-xadj[perm[step-1]-1];
  }

  cout<<"xadj: ";
  for(int step = 1; step<=n+1;++step){
    cout<<" "<<xadj[step-1];
  }
  cout<<endl;

  cout<<"newxadj: ";
  for(int step = 1; step<=n+1;++step){
    cout<<" "<<newxadj[step-1];
  }
  cout<<endl;

  for(int step = 1; step<=n;++step){
    int fn = newxadj[step-1];
    int ln = newxadj[step]-1;

    int oldlabel = perm[step-1];
    int ofn = xadj[oldlabel-1];
    int oln = xadj[oldlabel]-1;
    for(int i =ofn;i<=oln;++i){
      newadj[fn+i-ofn-1] = invperm[adj[i-1]-1];
    }

    //now sort them
    sort(&newadj[fn-1],&newadj[ln]);

  }
  
  cout<<"adj: ";
  for(int step = 1; step<=adj.size();++step){
    cout<<" "<<adj[step-1];
  }
  cout<<endl;

  cout<<"newadj: ";
  for(int step = 1; step<=newadj.size();++step){
    cout<<" "<<newadj[step-1];
  }
  cout<<endl;


  //Get the elimination tree
  ETree  tree;
//  tree.ConstructETree(n,&newxadj[0],&newadj[0]);
  tree.ConstructETree(n,&xadj[0],&adj[0]);

  vector<int> cc,rc;
//  GetLColRowCount(tree,&newxadj[0],&newadj[0],cc,rc);
  GetLColRowCount(tree,&xadj[0],&adj[0],cc,rc);

  int sum = 0;
  cout<<"Column count: ";
  for(int i =0; i<cc.size(); ++i){
    sum+= cc[i];
    cout<<" "<<cc[i];
  }
  cout<<endl;

  cout<<"Sum is "<<sum<<endl;

  displayMatrix(xadj, adj);
cout<<endl<<endl;
  displayMatrix(newxadj, newadj);
}


  void GetLColRowCount(ETree & tree,const int * xadj, const int * adj, vector<int> & cc, vector<int> & rc){
     //The tree need to be postordered
    if(!tree.IsPostOrdered()){
      tree.PostOrderTree();
    }

    int size = tree.Size();
 
    cc.resize(size);
    rc.resize(size);

    vector<int> level(size+1);
    vector<int> weight(size+1);
    vector<int> fdesc(size+1);
    vector<int> nchild(size+1);
    vector<int> set(size);
    vector<int> prvlf(size);
    vector<int> prvnbr(size);


        int xsup = 1;
        level[0] = 0;
      for(int k = size; k>=1; --k){
            rc[k-1] = 1;
            cc[k-1] = 0;
            set[k-1] = k;
            prvlf[k-1] = 0;
            prvnbr[k-1] = 0;
            level[k] = level[tree.PostParent(k-1)] + 1;
            weight[k] = 1;
            fdesc[k] = k;
            nchild[k] = 0;
      }

      nchild[0] = 0;
      fdesc[0] = 0;
      for(int k =1; k<size; ++k){
            int parent = tree.PostParent(k-1);
            weight[parent] = 0;
            ++nchild[parent];
            int ifdesc = fdesc[k];
            if  ( ifdesc < fdesc[parent] ) {
                fdesc[parent] = ifdesc;
            }
      }







      for(int lownbr = 1; lownbr<=size; ++lownbr){
        int lflag = 0;
        int ifdesc = fdesc[lownbr];
        int oldnbr = tree.FromPostOrder(lownbr);
        int jstrt = xadj[oldnbr-1];
        int jstop = xadj[oldnbr] - 1;


        //           -----------------------------------------------
        //           for each ``high neighbor'', hinbr of lownbr ...
        //           -----------------------------------------------
        for(int j = jstrt; j<=jstop;++j){
          int hinbr = tree.ToPostOrder(adj[j-1]);
          if  ( hinbr > lownbr )  {
            if  ( ifdesc > prvnbr[hinbr-1] ) {
              //                       -------------------------
              //                       increment weight(lownbr).
              //                       -------------------------
              ++weight[lownbr];
              int pleaf = prvlf[hinbr-1];
              //                       -----------------------------------------
              //                       if hinbr has no previous ``low neighbor'' 
              //                       then ...
              //                       -----------------------------------------
              if  ( pleaf == 0 ) {
                //                           -----------------------------------------
                //                           ... accumulate lownbr-->hinbr path length 
                //                               in rowcnt(hinbr).
                //                           -----------------------------------------
                rc[hinbr-1] += level[lownbr] - level[hinbr];
              }
              else{
                //                           -----------------------------------------
                //                           ... otherwise, lca <-- find(pleaf), which 
                //                               is the least common ancestor of pleaf 
                //                               and lownbr.
                //                               (path halving.)
                //                           -----------------------------------------
                int last1 = pleaf;
                int last2 = set[last1-1];
                int lca = set[last2-1];
                while(lca != last2){
                  set[last1-1] = lca;
                  last1 = lca;
                  last2 = set[last1-1];
                  lca = set[last2-1];
                }
                //                           -------------------------------------
                //                           accumulate pleaf-->lca path length in 
                //                           rowcnt(hinbr).
                //                           decrement weight(lca).
                //                           -------------------------------------
                rc[hinbr-1] += level[lownbr] - level[lca];
                --weight[lca];
              }
              //                       ----------------------------------------------
              //                       lownbr now becomes ``previous leaf'' of hinbr.
              //                       ----------------------------------------------
              prvlf[hinbr-1] = lownbr;
              lflag = 1;
            }
            //                   --------------------------------------------------
            //                   lownbr now becomes ``previous neighbor'' of hinbr.
            //                   --------------------------------------------------
            prvnbr[hinbr-1] = lownbr;
          }
        }
        //           ----------------------------------------------------
        //           decrement weight ( parent(lownbr) ).
        //           set ( p(lownbr) ) <-- set ( p(lownbr) ) + set(xsup).
        //           ----------------------------------------------------
        int parent = tree.PostParent(lownbr-1);
        --weight[parent];


        //merge the sets
        if  ( lflag == 1  || nchild[lownbr] >= 2 ) {
          xsup = lownbr;
        }
        set[xsup-1] = parent;
      }

        for(int k = 1; k<=size; ++k){
            int temp = cc[k-1] + weight[k];
            cc[k-1] = temp;
            int parent = tree.PostParent(k-1);
            if  ( parent != 0 ) {
                cc[parent-1] += temp;
            }
        }

  }


void displayMatrix(vector<int> & xadj, vector<int> & adj){
  for(int i = 1; i< xadj.size();++i){
    int fi = xadj[i-1];
    int li = xadj[i]-1;
    
    int prevrow = 0;
    int diagpassed = 0;
    for(int rowi = fi;rowi<=li;++rowi){
      int row = adj[rowi-1];
      if(row>i && !diagpassed){
        diagpassed=1;
        for(int prow = prevrow+1; prow<i; prow++){ 
          cout<<0<<" ";
        }
        cout<<1<<" ";
        prevrow = i;
      }   
   
      for(int prow = prevrow+1; prow<row; prow++){ 
        cout<<0<<" ";
      }
      cout<<1<<" ";
      prevrow = row;
    }
    if(!diagpassed){
      for(int prow = prevrow+1; prow<i; prow++){ cout<<0<<" ";}
      cout<<1<<" ";
      for(int prow = i+1; prow<xadj.size(); prow++){ cout<<0<<" ";}
    }
    else{
      for(int prow = prevrow+1; prow<xadj.size(); prow++){ cout<<0<<" ";}
    }
    cout<<endl;
  }
}

