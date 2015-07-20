#ifndef __TIMER_H__
#define __TIMER_H__

#if defined (PROFILE) || defined(PMPI)
  #define TAU
#endif



#include <mpi.h>
#include <sys/time.h>
#include "ngchol/Environment.hpp"
#include "ngchol/LogFile.hpp"
#include <iostream>
#include <string>


namespace LIBCHOLESKY{



class CTF_timer{
  public:
    char const * timer_name;
    int index;
    int exited;
    int original;

    std::vector<int> arr_index;
    std::vector<int> arr_exited;
    std::vector<int> arr_original;
  
  public:
    CTF_timer(char const * name);
    ~CTF_timer();
    void stop();
    void start();
    void exit();
    
};

void CTF_set_main_args(int argc, char * const * argv);
void CTF_set_context(MPI_Comm ctxt);

#ifdef TAU
#define TAU_FSTART(ARG)                                           \
  do { CTF_timer t(#ARG);/* logfileptr->OFS()<<"DEBUG START "<<#ARG<<std::endl;*/ t.start(); } while (0)

#define TAU_FSTOP(ARG)                                            \
  do { CTF_timer t(#ARG);/* logfileptr->OFS()<<"DEBUG STOP "<<#ARG<<std::endl;*/ t.stop(); } while (0)

#define TAU_PROFILE_TIMER(ARG1, ARG2, ARG3, ARG4)                 

#define TAU_PROFILE_INIT(argc, argv)                              \
  CTF_set_main_args(argc, argv);

#define TAU_PROFILE_SET_NODE(ARG)

#define TAU_PROFILE_START(ARG)                                    \
  CTF_timer __CTF_timer##ARG(#ARG);

#define TAU_PROFILE_STOP(ARG)                                     \
 __CTF_timer##ARG.stop();

#define TAU_PROFILE_SET_CONTEXT(ARG)                              \
  if (ARG==0) CTF_set_context(MPI_COMM_WORLD);                    \
  else CTF_set_context((MPI_Comm)ARG);
#endif

#ifdef PMPI
#define MPI_Bcast(...)                                            \
  { CTF_timer __t("MPI_Bcast");                                   \
              __t.start();                                        \
    PMPI_Bcast(__VA_ARGS__);                                      \
              __t.stop(); }
#define MPI_Reduce(...)                                           \
  { CTF_timer __t("MPI_Reduce");                                  \
              __t.start();                                        \
    PMPI_Reduce(__VA_ARGS__);                                     \
              __t.stop(); }
#define MPI_Wait(...)                                             \
  { CTF_timer __t("MPI_Wait");                                    \
              __t.start();                                        \
    PMPI_Wait(__VA_ARGS__);                                       \
              __t.stop(); }
#define MPI_Send(...)                                             \
  { CTF_timer __t("MPI_Send");                                    \
              __t.start();                                        \
    PMPI_Send(__VA_ARGS__);                                       \
              __t.stop(); }
#define MPI_Isend(...)                                            \
  { CTF_timer __t("MPI_Isend");                                   \
              __t.start();                                        \
    PMPI_Isend(__VA_ARGS__);                                      \
              __t.stop(); }
#define MPI_Allreduce(...)                                        \
  { CTF_timer __t("MPI_Allreduce");                               \
              __t.start();                                        \
    PMPI_Allreduce(__VA_ARGS__);                                  \
              __t.stop(); }
#define MPI_Allgather(...)                                        \
  { CTF_timer __t("MPI_Allgather");                               \
              __t.start();                                        \
    PMPI_Allgather(__VA_ARGS__);                                  \
              __t.stop(); }
#define MPI_Scatter(...)                                          \
  { CTF_timer __t("MPI_Scatter");                                 \
              __t.start();                                        \
    PMPI_Scatter(__VA_ARGS__);                                    \
              __t.stop(); }
#define MPI_Alltoall(...)                                         \
  { CTF_timer __t("MPI_Alltoall");                                \
              __t.start();                                        \
    PMPI_Alltoall(__VA_ARGS__);                                   \
              __t.stop(); }
#define MPI_Alltoallv(...)                                        \
  { CTF_timer __t("MPI_Alltoallv");                               \
              __t.start();                                        \
    PMPI_Alltoallv(__VA_ARGS__);                                  \
              __t.stop(); }
#define MPI_Gatherv(...)                                          \
  { CTF_timer __t("MPI_Gatherv");                                 \
              __t.start();                                        \
    PMPI_Gatherv(__VA_ARGS__);                                    \
              __t.stop(); }
#define MPI_Scatterv(...)                                         \
  { CTF_timer __t("MPI_Scatterv");                                \
              __t.start();                                        \
   PMPI_Scatterv(__VA_ARGS__);                                    \
              __t.stop(); }
#define MPI_Waitall(...)                                          \
  { CTF_timer __t("MPI_Waitall");                                 \
              __t.start();                                        \
    PMPI_Waitall(__VA_ARGS__);                                    \
              __t.stop(); }
#define MPI_Barrier(...)                                          \
  { CTF_timer __t("MPI_Barrier");                                 \
              __t.start();                                        \
    PMPI_Barrier(__VA_ARGS__);                                    \
              __t.stop(); }
#endif



#ifdef USE_TAU 
#define TIMER_START(a) TAU_START(TOSTRING(a));
#define TIMER_STOP(a) TAU_STOP(TOSTRING(a));
//#define scope_timer(a)
#define scope_timer(b,a) CTF_scope_timer b(#a)
#elif defined (PROFILE)
#define TIMER_START(a) TAU_FSTART(a);
#define TIMER_STOP(a) TAU_FSTOP(a);
#define scope_timer(b,a) CTF_scope_timer b(#a)
#else
#define TIMER_START(a)
#define TIMER_STOP(a)
#define scope_timer(b,a)
#endif

extern int iam;

class CTF_scope_timer{
  std::string name;
  public:
//  scope_timer(){
//    name = "DEFAULT";
//    TIMER_START(name);
//  }

  CTF_scope_timer(const char * pname){
    name = pname;
  //  TIMER_START(name.c_str());

//    if(iam==0){
//      pid_t pid = getpid();
//      std::cout<<"P"<<iam<<" is locked, pid is "<<pid<<std::endl;
//      volatile int lock = 1;
//      while (lock == 1){ }
//      std::cout<<"P"<<iam<<" is unlocked"<<std::endl;
//    }




    CTF_timer t(name.c_str()); 
    t.start();
  }
  ~CTF_scope_timer(){
  //  TIMER_STOP(name.c_str());

//    if(iam==0){
//      pid_t pid = getpid();
//      std::cout<<"P"<<iam<<" is locked, pid is "<<pid<<std::endl;
//      volatile int lock = 1;
//      while (lock == 1){ }
//      std::cout<<"P"<<iam<<" is unlocked"<<std::endl;
//    }

    CTF_timer t(name.c_str());

    t.stop();
  }
};



}

#endif

