#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <time.h>
#include <string>
#include <assert.h>
#include <iostream>
#include <vector>
#include <deque>
#include <cstring>
#include <sstream>

#include <omp.h>
#include "sympack.hpp"
#include "sympack/LogFile.hpp"
#include "sympack/Environment.hpp"
#include "sympack/timer.hpp"

#define MAX_NAME_LENGTH 38

#include <chrono>

namespace symPACK{
  using Clock = std::chrono::high_resolution_clock;
  using Duration = Clock::duration;
  using SecondDuration = std::chrono::duration<double, std::ratio<1,1> >;
  using Tick = Clock::time_point;

  int64_t main_argc = 0;
  char * const * main_argv;
  MPI_Comm comm;
#ifdef CUDA_MODE
  cublasHandle_t cublas_handler;
  cusolverDnHandle_t cusolver_handler;
  std::vector<cudaStream_t> streams;
  upcxx::device_allocator<upcxx::cuda_device> gpu_allocator;
  size_t gpu_alloc_size, gpu_block_limit, trsm_limit, potrf_limit, gemm_limit, syrk_limit;
  bool gpu_solve;
  FallbackType fallback_type;
  bool gpu_verbose;
#endif
  std::map<std::string, int> cpu_ops, gpu_ops;
  SecondDuration complete_time;
  int64_t set_contxt = 0;
  int64_t output_file_counter = 0;

  std::stringstream outstream;

  std::vector< Tick > arr_start_excl_tick;
  std::vector< Tick > arr_start_tick;
  std::vector< Tick > arr_excl_tick;
  std::vector< Duration > arr_complete_time;

  class function_timer;
  std::deque<function_timer> function_timers;

  class function_timer{
    public:
      char name[MAX_NAME_LENGTH];
      uint64_t calls;
      int numcore;

      std::vector<Tick> arr_start_tick;
      std::vector<Tick> arr_start_excl_tick;
      std::vector<std::vector<Duration> > arr_acc_time;
      std::vector<std::vector<Duration> > arr_acc_excl_time;
      std::vector<uint64_t> arr_calls;


      Duration total_time;
      Duration total_excl_time;
      uint64_t total_calls;

    public: 
      function_timer(char const * name_, 
          Tick const start_tick_,
          Tick const start_excl_tick_){
        sprintf(name, "%s", name_);

        numcore = omp_get_num_threads();

        arr_start_tick.assign(numcore,start_tick_); 
        arr_start_excl_tick.assign(numcore,start_excl_tick_); 
        arr_acc_time.resize(numcore); 
        arr_acc_excl_time.resize(numcore); 
        arr_calls.resize(numcore,0); 

        if (strlen(name) > MAX_NAME_LENGTH) {
          printf("function name must be fewer than %d characters\n",MAX_NAME_LENGTH);
          //assert(0);
        }
      }

      void compute_totals(MPI_Comm comm){ 
        if(set_contxt){
        }
        else{
          total_time=Duration(0);
          total_excl_time=Duration(0);
          total_calls=0;
          for(int64_t i =0;i<arr_acc_time.size();i++){
            total_time += std::accumulate(arr_acc_time[i].begin(),arr_acc_time[i].end(),Duration(0));
            total_excl_time += std::accumulate(arr_acc_excl_time[i].begin(),arr_acc_excl_time[i].end(),Duration(0));
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
        int64_t i;
        if (rank == 0){
          outstream.write(name,(strlen(name))*sizeof(char));
          char * space = (char*)malloc(MAX_NAME_LENGTH-strlen(name)+1);
          for (i=0; i<MAX_NAME_LENGTH-(int64_t)strlen(name); i++){
            space[i] = ' ';
          }
          space[i] = '\0';

          outstream.write(space,(strlen(space))*sizeof(char));
          char outstr[100];

          SecondDuration ttotal_time(total_time);
          SecondDuration ttotal_excl_time(total_excl_time);

          size_t num_samples = arr_acc_time.front().size();
          auto max_duration = [](std::vector<Duration> & vec)->Duration{ auto it = std::max_element(vec.begin(),vec.end()); return (it==vec.end())?Duration(0):*it; };
          auto it_max_acc_times = std::max_element(arr_acc_time.begin(),arr_acc_time.end(), [&](std::vector<Duration> & veca, std::vector<Duration> & vecb){ return max_duration(veca) > max_duration(vecb);});
          auto max_acc_time = it_max_acc_times!=arr_acc_time.end()?max_duration(*it_max_acc_times):Duration(0);

          sprintf(outstr,"%5ld    %lg  %3lu.%02lu  %lg  %3lu.%02lu %lg %lg %lg\n",
              total_calls/(np*numcore),
              (double)ttotal_time.count()/(np*numcore),
              (uint64_t)((100.*ttotal_time.count())/complete_time.count()),
              (uint64_t)((10000.*ttotal_time.count())/complete_time.count())%100,
              (double)(ttotal_excl_time.count())/(np*numcore),
              (uint64_t)((100.*ttotal_excl_time.count())/complete_time.count()),
              (uint64_t)((10000.*ttotal_excl_time.count())/complete_time.count())%100,
              (double)ttotal_time.count()/(num_samples*np*numcore),
              (double)(SecondDuration(max_acc_time).count()),
              (double)(ttotal_excl_time.count())/(num_samples*np*numcore));





          outstream.write(outstr,(strlen(outstr))*sizeof(char));

          free(space);



        } 
      }
  };


  bool comp_name(function_timer const & w1, function_timer const & w2) {
    return strcmp(w1.name, w2.name)>0;
  }


  symPACK_timer::symPACK_timer(const char * name){
#ifdef SPROFILE
    int64_t i;
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
      arr_start_tick.resize(numcore);
      arr_start_excl_tick.resize(numcore);
      arr_excl_tick.resize(numcore);

      int core = omp_get_thread_num();

      Tick cur_tick = Clock::now();
      arr_start_tick[core] = cur_tick; 
      arr_start_excl_tick[core] = cur_tick; 
      arr_excl_tick[core] = cur_tick; 

      function_timers.push_back(function_timer(name, cur_tick, cur_tick)); 
    } 
    else{
      for (i=0; i<(int64_t)function_timers.size(); i++){
        if (strcmp(function_timers[i].name, name) == 0){
          break;
        }
      }
      index = i;
      original = (index==0);
    }

    if (index == (int64_t)function_timers.size()) {
      Tick cur_tick = Clock::now();
      int core = omp_get_thread_num();

      if(arr_excl_tick.size()<core+1){ arr_excl_tick.resize(core+1); }

      arr_excl_tick[core] = cur_tick; 

      function_timers.push_back(function_timer(name, cur_tick, arr_excl_tick[core])); 
    }
    timer_name = name;
    exited = 0;
#endif
  }

  void symPACK_timer::start(){
#ifdef SPROFILE
    if (!exited){
      int64_t core = omp_get_thread_num();
      Tick cur_tick = Clock::now();
      if(function_timers[index].arr_start_tick.size()<core+1){
        function_timers[index].arr_start_tick.resize(core+1);
        function_timers[index].arr_start_excl_tick.resize(core+1);
      }
      function_timers[index].arr_start_tick[core] = cur_tick;

      if(arr_excl_tick.size()<core+1){
        arr_excl_tick.resize(core+1,cur_tick);
      }

      function_timers[index].arr_start_excl_tick[core] = arr_excl_tick[core];
    }
    else{
      abort();
    }
#endif
  }

  void symPACK_timer::stop(){
#ifdef SPROFILE
    if (!exited){
      int core = omp_get_thread_num();
      Tick cur_tick = Clock::now();

      if(function_timers[index].arr_acc_time.size()<core+1){
        function_timers[index].arr_acc_time.resize(core+1);
        function_timers[index].arr_acc_excl_time.resize(core+1);
      }

      if(arr_excl_tick.size()<core+1){
        arr_excl_tick.resize(core+1,cur_tick);
      }
      

      Duration delta = cur_tick - function_timers[index].arr_start_tick[core];

      function_timers[index].arr_acc_time[core].push_back(delta);
      function_timers[index].arr_acc_excl_time[core].push_back(delta - (arr_excl_tick[core]- function_timers[index].arr_start_excl_tick[core])); 

      arr_excl_tick[core] = function_timers[index].arr_start_excl_tick[core];
      arr_excl_tick[core] = arr_excl_tick[core] + Clock::duration(delta);

      function_timers[index].arr_calls[core]++;


      exit();
      exited = 1;
    }
    else{
      abort();
    }
#endif
  }

  symPACK_timer::~symPACK_timer(){
  }





  void symPACK_timer::exit(){
#ifdef SPROFILE
    if(exited){
      abort();
    }

    if (original && !exited) {
      int core, nc;

      int ismpi=1;

     
      core = omp_get_thread_num();
      int64_t i, j, p, len_symbols;
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
      sprintf(part,"       inclusive         exclusive        avg(incl)    max(incl)     avg(excl)\n");
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
      for (i=0; i<(int64_t)function_timers.size(); i++){
        sprintf(all_symbols+len_symbols, "%s", function_timers[i].name);
        len_symbols += strlen(function_timers[i].name)+1;
      }

      std::sort(function_timers.begin(), function_timers.end(),comp_name);
      for (i=0; i<(int64_t)function_timers.size(); i++){
        function_timers[i].compute_totals(comm);
      }
      std::sort(function_timers.begin(), function_timers.end());


      //TODO

      if(arr_start_tick.size()<core+1){
        arr_start_tick.resize(core+1);
      }
      if(arr_complete_time.size()<core+1){
        arr_complete_time.resize(core+1);
      }
      arr_complete_time[core] = Clock::now() - arr_start_tick[core];
      complete_time = SecondDuration(arr_complete_time[0]);

      for (i=0; i<(int64_t)function_timers.size(); i++){
        function_timers[i].print(comm,rank,np);
      }

      function_timers.clear();

      rank = -1;
      symPACK_Rank(&rank);
      char suffix[50];
      sprintf(suffix,"%d",rank);
      profileptr = new LogFile("sprofile",suffix);

      profileptr->OFS() << outstream.str();


      delete profileptr; 
    }
#endif
  }

  void symPACK_set_main_args(int64_t argc, char * const * argv){
    main_argv = argv;
    main_argc = argc;
  }

  void symPACK_set_context(MPI_Comm ctxt){
    set_contxt = 1;
    comm = ctxt;
  }


}

