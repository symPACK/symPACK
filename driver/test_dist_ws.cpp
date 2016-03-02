/// @file run_sparse_fanboth.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2014-01-23
//#define _DEBUG_


#include <mpi.h>
#include <time.h>
#include <omp.h>
#include <upcxx.h>

#include <chrono>
#include <random>
#include <cassert>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
#include <list>

using namespace std;

int iam, np;
ofstream ofs;


class task{
  public:
    int id;
    int weight; //dummy wait time
    size_t dataSize;
    upcxx::global_ptr<char> data;

    task(int pid, int pw){
      id = pid;
      weight = pw;
    }
    task():id(0),weight(0){}
    void work(){
      this_thread::sleep_for(chrono::seconds(weight));
    }
};

//setup upcxx structures




bool busy;
upcxx::shared_array<int64_t> workloads;
upcxx::global_ptr<task> taskArray;
task * pTaskArray = NULL;
list<task*> localTasks;
list<task*> stolenTasks;
int local_task_count = 0;
int64_t * pLocalLoad = NULL;

void reply_steal_task(bool success, upcxx::global_ptr<task> ptr){
  if( success){
    //get the task description
    upcxx::global_ptr<task> lPtr = upcxx::allocate<task>(iam,1);
    upcxx::copy(ptr,lPtr,1);

    //get the task data
    task * newTask = (task*)lPtr;
    upcxx::global_ptr<char> bak = newTask->data;
    newTask->data = upcxx::allocate<char>(iam,newTask->dataSize); 
    upcxx::copy(bak,newTask->data,newTask->dataSize);

    stolenTasks.push_back(newTask);
#ifdef VERBOSE
    ofs<<"P"<<iam<<" stole a task from P"<<ptr.where()<<endl;
#endif
  }
  busy = false;
}

void try_steal_task(int caller_rank){
#ifdef VERBOSE
  ofs<<"P"<<caller_rank<<" is trying to steal some work"<<endl;
#endif
  upcxx::global_ptr<task> ptr = upcxx::global_ptr<task>(NULL);
  bool success = false;
  if(*pLocalLoad>0 && localTasks.size()>0)
  {
    task * curTask = localTasks.front();
#ifdef VERBOSE
    assert(curTask!=NULL);
    ofs<<"WS Removing T"<<curTask->id<<endl;
#endif

    ptr = upcxx::global_ptr<task>(curTask); 
    
//    ofs<<(task*)ptr<<" "<<(task*)&pTaskArray[curTask->id]<<" "<<curTask<<endl;

    localTasks.pop_front();
    *pLocalLoad -= curTask->weight;
#ifdef VERBOSE
    ofs<<"After being robbed of T"<<curTask->id<<", local Load is now "<<*pLocalLoad<<endl;
#endif
    success = true;
  }
  upcxx::async(caller_rank)(reply_steal_task,success,ptr);
}




int main(int argc, char **argv) 
{
  MPI_Init(&argc,&argv);
  upcxx::init(&argc, &argv);

  MPI_Comm worldcomm;
  MPI_Comm_dup(MPI_COMM_WORLD,&worldcomm);
  MPI_Comm_size(worldcomm,&np);
  MPI_Comm_rank(worldcomm,&iam);

  stringstream sstr;
  sstr<<"logTest"<<iam;
  ofs.open(sstr.str()); 
  ofs<<"********* LOGFILE OF P"<<iam<<" *********"<<endl;
  ofs<<"**********************************"<<endl;

  if(iam==0){ cout<<"STARTING SPLITS"<<endl; }
  MPI_Comm workcomm;
  MPI_Comm_split(worldcomm,iam<np,iam,&workcomm);
  ofs<<"MPI SPLITS ARE DONE"<<endl;
  if(iam==0){ cout<<"MPI SPLITS ARE DONE"<<endl; }

  upcxx::team * workteam;
  int new_rank = (iam<np)?iam:iam-np;
  upcxx::team_all.split(iam<np,new_rank, workteam);

  ofs<<"ALL SPLITS ARE DONE"<<endl;

  if(iam==0){ cout<<"ALL SPLITS ARE DONE"<<endl; }

  if(iam<np){
    //  int np, iam;
    MPI_Comm_size(workcomm,&np);
    MPI_Comm_rank(workcomm,&iam);


    //initialize the workload shared array

    workloads.init(np);
    pLocalLoad = (int64_t*)workloads[iam].raw_ptr();
    int64_t & localLoad = *(pLocalLoad);
    //allocate local task lists
    int total_task_count = atoi(argv[1]);
    int dataSize = atoi(argv[2]);
    //int imbalance = 0;
    //if(argc>2){
    //  imbalance = atoi(argv[2]);
    //}   

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //    std::mt19937 generator(seed);
    //    double ratio = generator();

    std::default_random_engine generator;
    generator.seed(seed);
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double ratio = distribution(generator);
    ofs<< (double) total_task_count << endl;
    ofs<< ratio << endl;
    local_task_count = (int)(((double)total_task_count)*ratio); //iam<np-1?total_task_count/np:(total_task_count - iam*total_task_count/np);

    ofs<< local_task_count << endl;

    //int first_task = iam*total_task_count/np;

    //total_task_count += imbalance;
    //if(iam==np-1){
    //  local_task_count += imbalance;
    //}

    if(local_task_count>0){
      taskArray = upcxx::allocate<task>(iam,local_task_count);
      pTaskArray = (task*)taskArray;
      for(int i =0; i<local_task_count;i++){
        pTaskArray[i].dataSize = dataSize * distribution(generator);
        pTaskArray[i].data = upcxx::allocate<char>(iam,pTaskArray[i].dataSize); 
      }
    }



    localLoad = 0;
    for(int i =0; i<local_task_count;i++){
      localTasks.push_back(&pTaskArray[i]); 
      localTasks.back()->id = i;
      localTasks.back()->weight = 1;
      localLoad+= localTasks.back()->weight; 
    }

    ofs<<"Local Load is "<<localLoad<<endl;
    ofs<<"Local task size is "<<localTasks.size()<<endl;

    workteam->barrier(); 

    std::chrono::high_resolution_clock::time_point tstart = std::chrono::high_resolution_clock::now();
    //do the work
    while(localTasks.size()>0)
    {
      task * curTask = localTasks.front();
#ifdef VERBOSE
      assert(curTask!=NULL);
#endif
      localTasks.pop_front();
      curTask->work();
      localLoad -= curTask->weight;     
    }
    workteam->barrier(); 
    chrono::high_resolution_clock::time_point tend = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(tend - tstart);

    if(iam==0){
      cout << "It took " << time_span.count() << " seconds."<<endl;
    }  

localTasks.clear();
#ifdef VERBOSE
    assert(localTasks.empty());
#endif

    localLoad = 0;
    for(int i =0; i<local_task_count;i++){
      localTasks.push_back(&pTaskArray[i]); 
      localTasks.back()->id = i;
      localTasks.back()->weight = 1;
      localLoad+= localTasks.back()->weight; 
    }

    ofs<<"Local Load is "<<localLoad<<endl;
    ofs<<"Local task size is "<<localTasks.size()<<endl;

    workteam->barrier(); 

    tstart = std::chrono::high_resolution_clock::now();

    //do the work

    while(localTasks.size()>0)
    {
      upcxx::progress();
      if(localTasks.size()>0)
      {
        task * curTask = localTasks.front();
    //ofs<<"R "<<(task*)&pTaskArray[curTask->id]<<" "<<curTask<<endl;
        localTasks.pop_front();
#ifdef VERBOSE
        assert(curTask!=NULL);
        ofs<<"Removing T"<<curTask->id<<endl;
#endif
        curTask->work();
        localLoad -= curTask->weight;     
      }
    }

#ifdef VERBOSE
    ofs<<"ALL LOCAL TASKS ARE DONE"<<endl;
    ofs<<"Local task size is now "<<localTasks.size()<<endl;
    ofs<<"Local Load is now "<<localLoad<<endl;
    assert(localLoad==0);
#endif

#if 1
    //find a victim
    int64_t minLoad = -1;
    int victim = -1;
    do{
      minLoad = -1;
      victim = -1;
      for(int i = 0; i<np; i++){
        int64_t load = workloads[i];
        if(i!=iam && load>0 && (minLoad==-1 || load<minLoad) ){
          victim = i;
          minLoad = load;
        }
      }

      if(minLoad>0){
#ifdef VERBOSE
        ofs<<"Trying to steal from P"<<victim<<endl;
#endif
        busy = true;
        upcxx::async(victim)(try_steal_task,iam);
      }

      while(busy){upcxx::progress();}

      //now work with the stolen tasks I already have
      while(stolenTasks.size()>0){
        task * curTask = stolenTasks.front();
#ifdef VERBOSE
        assert(curTask!=NULL);
#endif
        stolenTasks.pop_front();
#ifdef VERBOSE
        ofs<<"Processing stolen task T"<<curTask->id<<endl;
#endif
        curTask->work();

        upcxx::deallocate(curTask->data);
        upcxx::deallocate(upcxx::global_ptr<task>(curTask));
      }
    }while(victim!=-1);
#endif

#ifdef VERBOSE
    ofs<<"ALL TASKS ARE DONE"<<endl;
#endif


    workteam->barrier(); 

    tend = chrono::high_resolution_clock::now();

    time_span = chrono::duration_cast<chrono::duration<double>>(tend - tstart);
    if(iam==0){
      cout << "It took " << time_span.count() << " seconds with WS."<<endl;
    }


    workteam->barrier(); 

    if(local_task_count>0){

      for(int i =0; i<local_task_count;i++){
        upcxx::deallocate(pTaskArray[i].data);
      }
      upcxx::deallocate(taskArray);
    }
    pTaskArray = NULL;

  }

  MPI_Barrier(workcomm);
  MPI_Comm_free(&workcomm);

  MPI_Barrier(worldcomm);
  MPI_Comm_free(&worldcomm);

  upcxx::finalize();
#if 0
  int finalized = 0;
  MPI_Finalized(&finalized);
  if(!finalized){ 
    MPI_Finalize();
  }
#endif
  ofs.close();
  return 0;
}


