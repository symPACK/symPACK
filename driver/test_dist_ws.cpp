/// @file run_sparse_fanboth.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2014-01-23
//#define _DEBUG_


#include <mpi.h>
#include <time.h>
#include <omp.h>
#include <upcxx.h>

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
upcxx::shared_array<int64_t> workloads;
upcxx::global_ptr<task> taskArray;
task * pTaskArray = NULL;
list<task*> localTasks;
list<task*> stolenTasks;
int head = 0;
int local_task_count = 0;
int64_t * pLocalLoad = NULL;

void reply_steal_task(bool success, upcxx::global_ptr<task> ptr){
  if(success){
    //get the task description
    upcxx::global_ptr<task> lPtr = upcxx::allocate<task>(iam,1);
    upcxx::copy(ptr,lPtr,1);
    stolenTasks.push_back((task*)lPtr);
    ofs<<"P"<<iam<<" stole a task from P"<<ptr.where()<<endl;
  }
}

void try_steal_task(int caller_rank){
  ofs<<"P"<<caller_rank<<" is trying to steal some work"<<endl;
  upcxx::global_ptr<task> ptr = upcxx::global_ptr<task>(NULL);
  bool success = false;
  if(*pLocalLoad>0){
    task * curTask = localTasks.front();
    *pLocalLoad -= curTask->weight;
    ofs<<"After being robbed, local Load is now "<<pLocalLoad<<endl;
    localTasks.pop_front();
    ptr = upcxx::global_ptr<task>(curTask); 
    bool success = true;
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
    int total_task_count = 100;
    
    local_task_count = iam<np-1?total_task_count/np:(total_task_count - iam*total_task_count/np);
    int first_task = iam*total_task_count/np;

    taskArray = upcxx::allocate<task>(iam,local_task_count);
    pTaskArray = (task*)taskArray;
    localLoad = 0;
    for(int i =0; i<total_task_count;i++){
      if(i>=first_task && i<first_task+local_task_count){
        localTasks.push_back(&pTaskArray[i]); 
        localTasks.front()->id = i;
        localTasks.front()->weight = 1;
        localLoad+= localTasks.front()->weight; 
      }
    }

    ofs<<"Local Load is "<<localLoad<<endl;

    workteam->barrier(); 
    //do the work
    while(localTasks.size()>0){
      upcxx::progress();

      task * curTask = localTasks.front();
      localTasks.pop_front();
      curTask->work();
      localLoad -= curTask->weight;      
    }

    
      ofs<<"ALL LOCAL TASKS ARE DONE"<<endl;
    ofs<<"Local Load is now "<<localLoad<<endl;
    assert(localLoad==0);

      
   //   
   //   upcxx::async_task<caller_rank>(reply_steal_task,success,ptr);
 
      //find a victim
      int64_t minLoad = -1;
      int victim = -1;
      do{
        minLoad = -1;
        victim = -1;
        for(int i = 0; i<np; i++){
          int64_t load = workloads[i];
          if(minLoad==-1 || load<minLoad){
            victim = i;
            minLoad = load;
          }
        }

        if(minLoad>0){
          ofs<<"Trying to steal from P"<<victim<<endl;
          upcxx::async(victim)(try_steal_task,iam);
        }
        upcxx::progress();

        //now work with the stolen tasks I already have
        while(stolenTasks.size()>0){
          upcxx::progress();

          task * curTask = stolenTasks.front();
          ofs<<"Processing stolen task "<<curTask->id<<endl;
          stolenTasks.pop_front();
          curTask->work();
          upcxx::deallocate(upcxx::global_ptr<task>(curTask));
        }

      }while(victim!=-1);

      ofs<<"ALL TASKS ARE DONE"<<endl;

    

    workteam->barrier(); 
    upcxx::deallocate(taskArray);
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


