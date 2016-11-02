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
#include "sympack/LogFile.hpp"
#include "sympack/Environment.hpp"
#include "sympack/timer.hpp"

#define MAX_NAME_LENGTH 38

#if 0
//#include <upcxx.h>
//upcxx timers
#include <gasnet.h>
#include <gasnet_tools.h>

namespace symPACK{

  int64_t main_argc = 0;
  char * const * main_argv;
  MPI_Comm comm;
  //double excl_time;
  double complete_time;
  int64_t set_contxt = 0;
  int64_t output_file_counter = 0;

  std::stringstream outstream;

  std::vector<gasnett_tick_t> arr_excl_tick;
  std::vector<gasnett_tick_t> arr_complete_tick;

  class function_timer{
    public:
      char name[MAX_NAME_LENGTH];
      gasnett_tick_t start_tick;
      gasnett_tick_t start_excl_tick;
      gasnett_tick_t acc_tick;
      gasnett_tick_t acc_excl_tick;
      int64_t calls;
      int64_t numcore;

      std::vector<gasnett_tick_t> arr_start_tick;
      std::vector<gasnett_tick_t> arr_start_excl_tick;
      std::vector<gasnett_tick_t> arr_acc_tick;
      std::vector<gasnett_tick_t> arr_acc_excl_tick;
      std::vector<int64_t> arr_calls;


      gasnett_tick_t total_tick;
      gasnett_tick_t total_excl_tick;
      int64_t total_calls;

    public: 
      function_timer(char const * name_, 
          gasnett_tick_t const start_tick_,
          gasnett_tick_t const start_excl_tick_){
        sprintf(name, "%s", name_);
        start_tick = start_tick_;
        start_excl_tick = start_excl_tick_;

        numcore = omp_get_num_threads();

        arr_start_tick.resize(numcore,0.0); 
        arr_start_excl_tick.resize(numcore,0.0); 
        arr_acc_tick.resize(numcore,0.0); 
        arr_acc_excl_tick.resize(numcore,0.0); 
        arr_calls.resize(numcore,0); 

        if (strlen(name) > MAX_NAME_LENGTH) {
          printf("function name must be fewer than %d characters\n",MAX_NAME_LENGTH);
          assert(0);
        }
        acc_tick = 0.0;
        acc_excl_tick = 0.0;
        calls = 0;
      }

      void compute_totals(MPI_Comm comm){ 
        if(set_contxt){
          //          MPI_Allreduce(&acc_time, &total_time, 1, 
          //              MPI_DOUBLE, MPI_SUM, comm);
          //          MPI_Allreduce(&acc_excl_time, &total_excl_time, 1, 
          //              MPI_DOUBLE, MPI_SUM, comm);
          //          MPI_Allreduce(&calls, &total_calls, 1, 
          //              MPI_INT, MPI_SUM, comm);
        }
        else{
          total_tick=0.0;
          total_excl_tick=0.0;
          total_calls=0;
          for(int64_t i =0;i<arr_acc_tick.size();i++){
            total_tick += arr_acc_tick[i];
            total_excl_tick += arr_acc_excl_tick[i];
            total_calls += arr_calls[i];
          }
        }
      }

      bool operator<(function_timer const & w) const {
        return total_tick > w.total_tick;
      }

      void print( 
          MPI_Comm const comm, 
          int64_t const      rank,
          int64_t const      np
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
          double ttotal = gasnett_ticks_to_ns(total_tick)/1.0E9;
          double texcl = gasnett_ticks_to_ns(total_excl_tick)/1.0E9;
          logfileptr->OFS()<<ttotal<<" "<< (double)(ttotal)/(np*numcore)<<" ";
          logfileptr->OFS()<<texcl<<" "<< (double)(texcl)/(np*numcore)<<std::endl;

          sprintf(outstr,"%5d    %lg  %3ld.%02ld  %lg  %3ld.%02ld\n",
              total_calls/(np*numcore),
              (double)(ttotal)/(np*numcore),
              (int64_t)(100.*((gasnett_ticks_to_ns(total_tick) / 1.0E9))/complete_time),
              ((int64_t)(10000.*((gasnett_ticks_to_ns(total_tick) / 1.0E9))/complete_time))%100,
              (double)(texcl)/(np*numcore),
              (int64_t)(100.*((gasnett_ticks_to_ns(total_excl_tick) / 1.0E9))/complete_time),
              ((int64_t)(10000.*((gasnett_ticks_to_ns(total_excl_tick) / 1.0E9))/complete_time))%100);

          outstream.write(outstr,(strlen(outstr))*sizeof(char));

          free(space);
        } 
      }
  };

  bool comp_name(function_timer const & w1, function_timer const & w2) {
    return strcmp(w1.name, w2.name)>0;
  }

  std::deque<function_timer> function_timers;

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

      int64_t numcore = omp_get_num_threads();
      arr_excl_tick.resize(numcore,0.0);

      gasnett_tick_t cur_tick;
      cur_tick = gasnett_ticks_now();

      function_timers.push_back(function_timer(name, cur_tick, 0.0)); 
    } 
    else{
      for (i=0; i<(int64_t)function_timers.size(); i++){
        if (strcmp(function_timers[i].name, name) == 0){
          /*function_timers[i].start_time = MPI_Wtime();
            function_timers[i].start_excl_time = excl_time;*/
          break;
        }
      }
      index = i;
      original = (index==0);
    }
    if (index == (int64_t)function_timers.size()) {
      gasnett_tick_t cur_tick;
      cur_tick = gasnett_ticks_now();
      int64_t core = omp_get_thread_num();

      if(arr_excl_tick.size()<core+1){
        arr_excl_tick.resize(core+1,0.0);
      }
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
      gasnett_tick_t cur_tick;
      cur_tick = gasnett_ticks_now();
      if(function_timers[index].arr_start_tick.size()<core+1){
        function_timers[index].arr_start_tick.resize(core+1);
        function_timers[index].arr_start_excl_tick.resize(core+1);
      }
      function_timers[index].arr_start_tick[core] = cur_tick;

      if(arr_excl_tick.size()<core+1){
        arr_excl_tick.resize(core+1,0.0);
      }

      function_timers[index].arr_start_excl_tick[core] = arr_excl_tick[core];
      function_timers[index].start_tick = cur_tick;
    }
    else{
      abort();
    }
#endif
  }

  void symPACK_timer::stop(){
#ifdef SPROFILE
    if (!exited){
      int64_t core = omp_get_thread_num();
      gasnett_tick_t cur_tick;
      cur_tick = gasnett_ticks_now();

      if(function_timers[index].arr_acc_tick.size()<core+1){
        function_timers[index].arr_acc_tick.resize(core+1);
        function_timers[index].arr_acc_excl_tick.resize(core+1);
      }

      if(arr_excl_tick.size()<core+1){
        arr_excl_tick.resize(core+1,0.0);
      }
      gasnett_tick_t delta_tick = cur_tick - function_timers[index].start_tick;
      logfileptr->OFS()<<"cur_tick "<<cur_tick<<" delta "<<delta_tick<<std::endl;
      function_timers[index].arr_acc_tick[core] += delta_tick;
      function_timers[index].arr_acc_excl_tick[core] += delta_tick - 
        (arr_excl_tick[core]- function_timers[index].arr_start_excl_tick[core]); 

      arr_excl_tick[core] = function_timers[index].arr_start_excl_tick[core] + delta_tick;

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
      int64_t core, nc;

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

      //      int64_t ismpi=0;
      //      MPI_Initialized( &ismpi);
      int64_t ismpi=1;


      //      if (ismpi && set_contxt){ 
      //        int64_t rank, np, i, j, p, len_symbols;
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
      //          for (i=0; i<(int64_t)function_timers.size(); i++){
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
      //              for (i=0; i<(int64_t)function_timers.size(); i++){
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
      //            symPACK_timer t(all_symbols+j);
      //            j+=strlen(all_symbols+j)+1;
      //          }
      //        }
      //
      //        std::sort(function_timers.begin(), function_timers.end(),comp_name);
      //        for (i=0; i<(int64_t)function_timers.size(); i++){
      //          function_timers[i].compute_totals(comm);
      //        }
      //        std::sort(function_timers.begin(), function_timers.end());
      //        complete_time = function_timers[0].total_time;
      //        for (i=0; i<(int64_t)function_timers.size(); i++){
      //          function_timers[i].print(comm,rank,np);
      //        }
      //
      //        function_timers.clear();
      //      }  
      //      else{
      int64_t i, j, p, len_symbols;
      int64_t np, rank;

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
      for (i=0; i<(int64_t)function_timers.size(); i++){
        sprintf(all_symbols+len_symbols, "%s", function_timers[i].name);
        len_symbols += strlen(function_timers[i].name)+1;
      }



      std::sort(function_timers.begin(), function_timers.end(),comp_name);
      for (i=0; i<(int64_t)function_timers.size(); i++){
        function_timers[i].compute_totals(comm);
      }
      std::sort(function_timers.begin(), function_timers.end());
      complete_time = gasnett_ticks_to_ns(function_timers[0].total_tick)/1.0E9;
      for (i=0; i<(int64_t)function_timers.size(); i++){
        function_timers[i].print(comm,rank,np);
      }

      function_timers.clear();
      //      }

      //      int64_t iam=0;
      //      if(!ismpi){
      //        //      iam = MYTHREAD;
      //        abort();
      //      }
      //      else{
      //        MPI_Comm_rank(MPI_COMM_WORLD,&iam);
      //      }

      char suffix[50];
      sprintf(suffix,"%d",iam);
      profileptr = new LogFile("sprofile",suffix);

      profileptr->OFS() << outstream.str();


      delete profileptr; 

      //logfileptr->OFS()<<"CALLED"<<std::endl;





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

#else

#include <chrono>

namespace symPACK{
  using Clock = std::chrono::high_resolution_clock;
  //using Duration = std::chrono::duration<double, std::ratio<1,1> >;
  using Duration = Clock::duration;
  using SecondDuration = std::chrono::duration<double, std::ratio<1,1> >;
  using Tick = Clock::time_point;

  int64_t main_argc = 0;
  char * const * main_argv;
  MPI_Comm comm;
  //double excl_time;
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
      //Tick start_tick;
      //Tick start_excl_tick;
      //Tick acc_tick;
      //Tick acc_excl_tick;
      uint64_t calls;
      int numcore;

      std::vector<Tick> arr_start_tick;
      std::vector<Tick> arr_start_excl_tick;
      std::vector<Duration> arr_acc_time;
      std::vector<Duration> arr_acc_excl_time;
      std::vector<uint64_t> arr_calls;


      Duration total_time;
      Duration total_excl_time;
      uint64_t total_calls;

    public: 
      function_timer(char const * name_, 
          Tick const start_tick_,
          Tick const start_excl_tick_){
        sprintf(name, "%s", name_);
        //start_tick = start_tick_;
        //start_excl_tick = start_excl_tick_;

        numcore = omp_get_num_threads();

        arr_start_tick.assign(numcore,start_tick_); 
        arr_start_excl_tick.assign(numcore,start_excl_tick_); 
        arr_acc_time.resize(numcore,Duration(0)); 
        arr_acc_excl_time.resize(numcore,Duration(0)); 
        arr_calls.resize(numcore,0); 

        if (strlen(name) > MAX_NAME_LENGTH) {
          printf("function name must be fewer than %d characters\n",MAX_NAME_LENGTH);
          assert(0);
        }
        //acc_tick = 0.0;
        //acc_excl_tick = 0.0;
        //calls = 0;
      }

      void compute_totals(MPI_Comm comm){ 
        if(set_contxt){
          //          MPI_Allreduce(&acc_time, &total_time, 1, 
          //              MPI_DOUBLE, MPI_SUM, comm);
          //          MPI_Allreduce(&acc_excl_time, &total_excl_time, 1, 
          //              MPI_DOUBLE, MPI_SUM, comm);
          //          MPI_Allreduce(&calls, &total_calls, 1, 
          //              MPI_INT, MPI_SUM, comm);
        }
        else{
          total_time=Duration(0);
          total_excl_time=Duration(0);
          total_calls=0;
          for(int64_t i =0;i<arr_acc_time.size();i++){
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
          int64_t const      rank,
          int64_t const      np
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
 //         double ttotal = gasnett_ticks_to_ns(total_tick)/1.0E9;
 //         double texcl = gasnett_ticks_to_ns(total_excl_tick)/1.0E9;
          
  //        std::chrono::duration<double, std::milli> fp_ms = t2 - t1;


//          logfileptr->OFS()<<ttotal<<" "<< (double)(ttotal)/(np*numcore)<<" ";
//          logfileptr->OFS()<<texcl<<" "<< (double)(texcl)/(np*numcore)<<std::endl;
          SecondDuration ttotal_time(total_time);
          SecondDuration ttotal_excl_time(total_excl_time);
          sprintf(outstr,"%5d    %lg  %3lu.%02lu  %lg  %3lu.%02lu\n",
              total_calls/(np*numcore),
              (double)ttotal_time.count()/(np*numcore),
              (uint64_t)((100.*ttotal_time.count())/complete_time.count()),
              (uint64_t)((10000.*ttotal_time.count())/complete_time.count())%100,
              (double)(ttotal_excl_time.count())/(np*numcore),
              (uint64_t)((100.*ttotal_excl_time.count())/complete_time.count()),
              (uint64_t)((10000.*ttotal_excl_time.count())/complete_time.count())%100);

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

      //if(arr_start_tick.size()<core+1){ arr_start_tick.resize(core+1); }
      //if(arr_start_excl_tick.size()<core+1){ arr_start_excl_tick.resize(core+1); }
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
//      function_timers[index].start_tick = cur_tick;
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

      function_timers[index].arr_acc_time[core] += delta;
      function_timers[index].arr_acc_excl_time[core] += delta - (arr_excl_tick[core]- function_timers[index].arr_start_excl_tick[core]); 

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

      //      int64_t ismpi=0;
      //      MPI_Initialized( &ismpi);
      int ismpi=1;


      //      if (ismpi && set_contxt){ 
      //        int64_t rank, np, i, j, p, len_symbols;
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
      //          for (i=0; i<(int64_t)function_timers.size(); i++){
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
      //              for (i=0; i<(int64_t)function_timers.size(); i++){
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
      //            symPACK_timer t(all_symbols+j);
      //            j+=strlen(all_symbols+j)+1;
      //          }
      //        }
      //
      //        std::sort(function_timers.begin(), function_timers.end(),comp_name);
      //        for (i=0; i<(int64_t)function_timers.size(); i++){
      //          function_timers[i].compute_totals(comm);
      //        }
      //        std::sort(function_timers.begin(), function_timers.end());
      //        complete_time = function_timers[0].total_time;
      //        for (i=0; i<(int64_t)function_timers.size(); i++){
      //          function_timers[i].print(comm,rank,np);
      //        }
      //
      //        function_timers.clear();
      //      }  
      //      else{
      
      core = omp_get_thread_num();
      int64_t i, j, p, len_symbols;
      int64_t np, rank;

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

      char suffix[50];
      sprintf(suffix,"%d",upcxx::myrank());
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
#endif

