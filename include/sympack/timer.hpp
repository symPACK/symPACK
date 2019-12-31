#ifndef __TIMER_H__
#define __TIMER_H__

#include <mpi.h>
#include <sys/time.h>
#include "sympack/Environment.hpp"
#include "sympack/LogFile.hpp"
#include <iostream>
#include <string>


namespace symPACK{



class symPACK_timer{
  public:
    char const * timer_name;
    int64_t index;
    int64_t exited;
    int64_t original;

    vector<int64_t> arr_index;
    vector<int64_t> arr_exited;
    vector<int64_t> arr_original;
  
  public:
    symPACK_timer(char const * name);
    ~symPACK_timer();
    void stop();
    void start();
    void exit();
    
};

void symPACK_set_main_args(int64_t argc, char * const * argv);
void symPACK_set_context(MPI_Comm ctxt);

//extern Int iam;

class symPACK_scope_timer{
  std::string name;
  public:

  symPACK_scope_timer(const char * pname){
    name = pname;
    symPACK_timer t(name.c_str()); 
    t.start();
  }
  ~symPACK_scope_timer(){
    symPACK_timer t(name.c_str());
    t.stop();
  }
};

#ifdef SPROFILE
#define SYMPACK_FSTART(ARG)                                           \
  do { symPACK::symPACK_timer t(#ARG);/* logfileptr->OFS()<<"DEBUG START "<<#ARG<<std::endl;*/ t.start(); } while (0)

#define SYMPACK_FSTOP(ARG)                                            \
  do { symPACK::symPACK_timer t(#ARG);/* logfileptr->OFS()<<"DEBUG STOP "<<#ARG<<std::endl;*/ t.stop(); } while (0)


#define SYMPACK_SPROFILE_INIT(argc, argv)                              \
  symPACK::symPACK_set_main_args(argc, argv);

#define SYMPACK_PROFILE_SET_CONTEXT(ARG)                              \
  if (ARG==0) symPACK::symPACK_set_context(symPACK::world_comm);                    \
  else symPACK::symPACK_set_context((MPI_Comm)ARG);

#define SYMPACK_TIMER_START(a) SYMPACK_FSTART(a);
#define SYMPACK_TIMER_STOP(a) SYMPACK_FSTOP(a);
#define SYMPACK_TIMER_SPECIAL_START(a) SYMPACK_FSTART(a);
#define SYMPACK_TIMER_SPECIAL_STOP(a) SYMPACK_FSTOP(a);
#define scope_timer(b,a) symPACK::symPACK_scope_timer b(#a)
#define scope_timer_special(b,a) symPACK::symPACK_scope_timer b(#a)

#else
#define SYMPACK_TIMER_START(a)
#define SYMPACK_TIMER_STOP(a)
#define SYMPACK_TIMER_SPECIAL_START(a)
#define SYMPACK_TIMER_SPECIAL_STOP(a)
#define scope_timer(b,a)
#define scope_timer_special(b,a)
#endif



}

#endif //__TIMER_H__

