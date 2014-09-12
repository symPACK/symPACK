#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <time.h>
#include <string>
#include <assert.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <sstream>

#include <omp.h>
#include "ngchol/LogFile.hpp"
#include "ngchol/Environment.hpp"
#include "ngchol/timer.hpp"

#define MAX_NAME_LENGTH 38

//#include <upcxx.h>

namespace LIBCHOLESKY{

  int main_argc = 0;
  char * const * main_argv;
  MPI_Comm comm;
  //double excl_time;
  double complete_time;
  int set_contxt = 0;
  int output_file_counter = 0;

  std::stringstream outstream;

  std::vector<double> arr_excl_time;
  std::vector<double> arr_complete_time;

  class function_timer{
    public:
      char name[MAX_NAME_LENGTH];
      double start_time;
      double start_excl_time;
      double acc_time;
      double acc_excl_time;
      int calls;
      int numcore;

      std::vector<double> arr_start_time;
      std::vector<double> arr_start_excl_time;
      std::vector<double> arr_acc_time;
      std::vector<double> arr_acc_excl_time;
      std::vector<int> arr_calls;


      double total_time;
      double total_excl_time;
      int total_calls;

    public: 
      function_timer(char const * name_, 
          double const start_time_,
          double const start_excl_time_){
        sprintf(name, "%s", name_);
        start_time = start_time_;
        start_excl_time = start_excl_time_;

        numcore = omp_get_num_threads();

        arr_start_time.resize(numcore,0.0); 
        arr_start_excl_time.resize(numcore,0.0); 
        arr_acc_time.resize(numcore,0.0); 
        arr_acc_excl_time.resize(numcore,0.0); 
        arr_calls.resize(numcore,0); 

        if (strlen(name) > MAX_NAME_LENGTH) {
          printf("function name must be fewer than %d characters\n",MAX_NAME_LENGTH);
          assert(0);
        }
        acc_time = 0.0;
        acc_excl_time = 0.0;
        calls = 0;
      }

      void compute_totals(MPI_Comm comm){ 
        if(set_contxt){
          MPI_Allreduce(&acc_time, &total_time, 1, 
              MPI_DOUBLE, MPI_SUM, comm);
          MPI_Allreduce(&acc_excl_time, &total_excl_time, 1, 
              MPI_DOUBLE, MPI_SUM, comm);
          MPI_Allreduce(&calls, &total_calls, 1, 
              MPI_INT, MPI_SUM, comm);
        }
        else{
          total_time=0.0;
          total_excl_time=0.0;
          total_calls=0;
          for(int i =0;i<arr_acc_time.size();i++){
            total_time += arr_acc_time[i];
            total_excl_time += arr_acc_excl_time[i];
            total_calls += arr_calls[i];
          }
        }
      }

      bool operator<(function_timer const & w) const {
        return total_time > w.total_time;
      }

      void print( 
          MPI_Comm const comm, 
          int const      rank,
          int const      np
          ){
        int i;
        if (rank == 0){
          outstream.write(name,(strlen(name))*sizeof(char));
          char * space = (char*)malloc(MAX_NAME_LENGTH-strlen(name)+1);
          for (i=0; i<MAX_NAME_LENGTH-(int)strlen(name); i++){
            space[i] = ' ';
          }
          space[i] = '\0';

          outstream.write(space,(strlen(space))*sizeof(char));
          char outstr[100];
          sprintf(outstr,"%5d    %7.7lg  %3d.%02d  %7.7lg  %3d.%02d\n",
              total_calls/(np*numcore),
              (double)(total_time/(np*numcore)),
              (int)(100.*(total_time)/complete_time),
              ((int)(10000.*(total_time)/complete_time))%100,
              (double)(total_excl_time/(np*numcore)),
              (int)(100.*(total_excl_time)/complete_time),
              ((int)(10000.*(total_excl_time)/complete_time))%100);

          outstream.write(outstr,(strlen(outstr))*sizeof(char));

          free(space);
        } 
      }
  };

  bool comp_name(function_timer const & w1, function_timer const & w2) {
    return strcmp(w1.name, w2.name)>0;
  }

  std::vector<function_timer> function_timers;

  CTF_timer::CTF_timer(const char * name){
#ifdef PROFILE
    int i;
    if (function_timers.size() == 0) {
      if (name[0] == 'M' && name[1] == 'P' && 
          name[2] == 'I' && name[3] == '_'){
        exited = 1;
        original = 0;
        return;
      }
      original = 1;
      index = 0;

      int numcore = omp_get_num_threads();
      arr_excl_time.resize(numcore,0.0);

      double time = get_time();
      function_timers.push_back(function_timer(name, time, 0.0)); 
    } 
    else{
      for (i=0; i<(int)function_timers.size(); i++){
        if (strcmp(function_timers[i].name, name) == 0){
          /*function_timers[i].start_time = MPI_Wtime();
            function_timers[i].start_excl_time = excl_time;*/
          break;
        }
      }
      index = i;
      original = (index==0);
    }
    if (index == (int)function_timers.size()) {
      double time = get_time();
      int core = omp_get_thread_num();

      if(arr_excl_time.size()<core+1){
        arr_excl_time.resize(core+1,0.0);
      }
      function_timers.push_back(function_timer(name, time, arr_excl_time[core])); 
    }
    timer_name = name;
    exited = 0;
#endif
  }

  void CTF_timer::start(){
#ifdef PROFILE
    if (!exited){
      int core = omp_get_thread_num();
      double time = get_time();
      if(function_timers[index].arr_start_time.size()<core+1){
        function_timers[index].arr_start_time.resize(core+1);
        function_timers[index].arr_start_excl_time.resize(core+1);
      }
      function_timers[index].arr_start_time[core] = time;

      if(arr_excl_time.size()<core+1){
        arr_excl_time.resize(core+1,0.0);
      }

      function_timers[index].arr_start_excl_time[core] = arr_excl_time[core];
      function_timers[index].start_time = time;
    }
#endif
  }

  void CTF_timer::stop(){
#ifdef PROFILE
    if (!exited){
      int core = omp_get_thread_num();
      double time = get_time();
      if(function_timers[index].arr_acc_time.size()<core+1){
        function_timers[index].arr_acc_time.resize(core+1);
        function_timers[index].arr_acc_excl_time.resize(core+1);
      }

      if(arr_excl_time.size()<core+1){
        arr_excl_time.resize(core+1,0.0);
      }
      double delta_time = time - function_timers[index].start_time;
      function_timers[index].arr_acc_time[core] += delta_time;
      function_timers[index].arr_acc_excl_time[core] += delta_time - 
        (arr_excl_time[core]- function_timers[index].arr_start_excl_time[core]); 

      arr_excl_time[core] = function_timers[index].arr_start_excl_time[core] + delta_time;

      function_timers[index].arr_calls[core]++;


      exit();
      exited = 1;
    }
#endif
  }

  CTF_timer::~CTF_timer(){
  }





  void CTF_timer::exit(){
#ifdef PROFILE
    if (original && !exited) {
      int core, nc;

      //#ifndef UPCXX
      //#pragma omp parallel
      //    {
      //      nc = omp_get_num_threads();
      //      core = omp_get_thread_num();
      //    }
      //#else
      //    nc = 1;
      //    core=0;
      //#endif

//      int ismpi=0;
//      MPI_Initialized( &ismpi);
      int ismpi=1;


//      if (ismpi && set_contxt){ 
//        int rank, np, i, j, p, len_symbols;
//
//        if(ismpi){
//          MPI_Comm_rank(comm, &rank);
//          MPI_Comm_size(comm, &np);
//        }
//        //      else{
//        //        rank=MYTHREAD;
//        //        np=THREADS;
//        //      }
//
//
//
//        char all_symbols[10000];
//
//        if (rank == 0){
//          char part[300];
//
//          char heading[MAX_NAME_LENGTH+200];
//          for (i=0; i<MAX_NAME_LENGTH; i++){
//            part[i] = ' ';
//          }
//          part[i] = '\0';
//          sprintf(heading,"%s",part);
//          sprintf(part,"       inclusive         exclusive\n");
//          strcat(heading,part);
//
//          outstream.write(heading,(strlen(heading))*sizeof(char));
//
//          for (i=0; i<MAX_NAME_LENGTH; i++){
//            part[i] = ' ';
//          }
//          part[i] = '\0';
//          sprintf(heading,"%s",part);
//          sprintf(part, "calls        sec       %%"); 
//          strcat(heading,part);
//          sprintf(part, "       sec       %%\n"); 
//          strcat(heading,part);
//          outstream.write(heading,(strlen(heading))*sizeof(char));
//
//          len_symbols = 0;
//          for (i=0; i<(int)function_timers.size(); i++){
//            sprintf(all_symbols+len_symbols, "%s", function_timers[i].name);
//            len_symbols += strlen(function_timers[i].name)+1;
//          }
//
//        }
//        if (np > 1){
//          for (p=0; p<np; p++){
//            if (rank == p){
//              MPI_Send(&len_symbols, 1, MPI_INT, (p+1)%np, 1, comm);
//              MPI_Send(all_symbols, len_symbols, MPI_CHAR, (p+1)%np, 2, comm);
//            }
//            if (rank == (p+1)%np){
//              MPI_Status stat;
//              MPI_Recv(&len_symbols, 1, MPI_INT, p, 1, comm, &stat);
//              MPI_Recv(all_symbols, len_symbols, MPI_CHAR, p, 2, comm, &stat);
//              for (i=0; i<(int)function_timers.size(); i++){
//                j=0;
//                while (j<len_symbols && strcmp(all_symbols+j, function_timers[i].name) != 0){
//                  j+=strlen(all_symbols+j)+1;
//                }
//
//                if (j>=len_symbols){
//                  sprintf(all_symbols+len_symbols, "%s", function_timers[i].name);
//                  len_symbols += strlen(function_timers[i].name)+1;
//                }
//              }
//            }
//          }
//          MPI_Bcast(&len_symbols, 1, MPI_INT, 0, comm);
//          MPI_Bcast(all_symbols, len_symbols, MPI_CHAR, 0, comm);
//          j=0;
//          while (j<len_symbols){
//            CTF_timer t(all_symbols+j);
//            j+=strlen(all_symbols+j)+1;
//          }
//        }
//
//        std::sort(function_timers.begin(), function_timers.end(),comp_name);
//        for (i=0; i<(int)function_timers.size(); i++){
//          function_timers[i].compute_totals(comm);
//        }
//        std::sort(function_timers.begin(), function_timers.end());
//        complete_time = function_timers[0].total_time;
//        for (i=0; i<(int)function_timers.size(); i++){
//          function_timers[i].print(comm,rank,np);
//        }
//
//        function_timers.clear();
//      }  
//      else{
        int i, j, p, len_symbols;
        int np, rank;

        np = 1;
        rank = 0;

        char all_symbols[10000];

        char part[300];

        char heading[MAX_NAME_LENGTH+200];
        for (i=0; i<MAX_NAME_LENGTH; i++){
          part[i] = ' ';
        }
        part[i] = '\0';
        sprintf(heading,"%s",part);
        sprintf(part,"       inclusive         exclusive\n");
        strcat(heading,part);
        outstream.write(heading,(strlen(heading))*sizeof(char));

        for (i=0; i<MAX_NAME_LENGTH; i++){
          part[i] = ' ';
        }
        part[i] = '\0';
        sprintf(heading,"%s",part);
        sprintf(part, "calls        sec       %%"); 
        strcat(heading,part);
        sprintf(part, "       sec       %%\n"); 
        strcat(heading,part);
        outstream.write(heading,(strlen(heading))*sizeof(char));

        len_symbols = 0;
        for (i=0; i<(int)function_timers.size(); i++){
          sprintf(all_symbols+len_symbols, "%s", function_timers[i].name);
          len_symbols += strlen(function_timers[i].name)+1;
        }



        std::sort(function_timers.begin(), function_timers.end(),comp_name);
        for (i=0; i<(int)function_timers.size(); i++){
          function_timers[i].compute_totals(comm);
        }
        std::sort(function_timers.begin(), function_timers.end());
        complete_time = function_timers[0].total_time;
        for (i=0; i<(int)function_timers.size(); i++){
          function_timers[i].print(comm,rank,np);
        }

        function_timers.clear();
//      }

//      int iam=0;
//      if(!ismpi){
//        //      iam = MYTHREAD;
//        abort();
//      }
//      else{
//        MPI_Comm_rank(MPI_COMM_WORLD,&iam);
//      }
      char suffix[50];
      sprintf(suffix,"%d",iam);
      profileptr = new LogFile("profile",suffix);

      profileptr->OFS() << outstream.str();


      delete profileptr; 

      //logfileptr->OFS()<<"CALLED"<<std::endl;





    }
#endif
  }

  void CTF_set_main_args(int argc, char * const * argv){
    main_argv = argv;
    main_argc = argc;
  }

  void CTF_set_context(MPI_Comm ctxt){
    set_contxt = 1;
    comm = ctxt;
  }
}
