#ifndef _UPCXX_ADDITIONS_IMPL_HPP_
#define _UPCXX_ADDITIONS_IMPL_HPP_



#include "upcxx_additions.hpp"
#include "LogFile.hpp"

extern LogFile * logfileptr; 

namespace upcxx{

    template <typename T> void new_fctn(global_ptr<T> gptr){
      new ((T*)gptr) T();
    }
    template <typename T> void delete_fctn(global_ptr<T> gptr){
      logfileptr->OFS()<<"Deleting"<<std::endl;
      T * toto = gptr;
      logfileptr->OFS()<<"ptr = "<<(long)toto<<endl;
//      toto->~T();
      ((T*)gptr)->~T();
      logfileptr->OFS()<<"Deleted"<<std::endl;
    }
    template <typename T> global_ptr<T>& Create(int where=MYTHREAD){
      //currently it only works with local assignments
      global_ptr<T> gptr = allocate<T>(where,1);
      async(where)(new_fctn<T>,gptr);
      upcxx::wait();
      return gptr;
    }
    template <typename T> void Destroy(global_ptr<T>& gptr){
      //currently it only works with local assignments
      int where = gptr.tid();
      logfileptr->OFS()<<"Calling delete on P"<<where<<std::endl;
        async(where)(delete_fctn<T>,gptr);
      logfileptr->OFS()<<"Called delete"<<std::endl;
      upcxx::wait();
      deallocate<T>(gptr);
    }

}

#endif
