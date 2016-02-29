/// @file run_sparse_fanboth.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Mathias Jacquelin
/// @date 2014-01-23
//#define _DEBUG_


#include <mpi.h>
#include <time.h>
#include <omp.h>
#include <upcxx.h>

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
#include <list>

using namespace std;




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
upcxx::global_ptr<task> taskArray = upcxx::global_ptr<task>(NULL);
task * pTaskArray = NULL;
list<task*> localTasks;
int head = 0;
int local_task_count = 0;

void try_steal_task(){

}

void reply_steal_task(){
}



int main(int argc, char **argv) 
{
  MPI_Init(&argc,&argv);
  upcxx::init(&argc, &argv);

  MPI_Comm worldcomm;
  MPI_Comm_dup(MPI_COMM_WORLD,&worldcomm);
  int np, iam;
  MPI_Comm_size(worldcomm,&np);
  MPI_Comm_rank(worldcomm,&iam);

  stringstream sstr;
  sstr<<"logTest"<<iam;
  ofstream ofs(sstr.str()); 
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
    int64_t localLoad = 0;
    //allocate local task lists
    head = 0; 
    local_task_count = 100;
    taskArray = upcxx::allocate<task>(iam,local_task_count);
    pTaskArray = (task*)taskArray;
    for(int i =0; i<local_task_count;i++){ localTasks.push_back(&pTaskArray[i]); localLoad+= localTasks.front()->weight; }
    workloads[iam] = localLoad;

    workteam->barrier(); 
    //do the work
    while(localTasks.size()>0){
      upcxx::progress();

      task * curTask = localTasks.front();
      localTasks.pop_front();

      curTask->work();
      
    }

    //if nothing else to do, try to steal something from remote processes
    


    workteam->barrier(); 
    upcxx::deallocate(taskArray);
    pTaskArray = 0;

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
  return 0;
}


