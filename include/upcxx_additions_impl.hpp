#ifndef _UPCXX_ADDITIONS_IMPL_HPP_
#define _UPCXX_ADDITIONS_IMPL_HPP_



#include "upcxx_additions.hpp"


namespace upcxx{

    template <typename T> void new_fctn(global_ptr<T> gptr){
      new ((T*)gptr) T();
    }
    template <typename T> void delete_fctn(global_ptr<T> gptr){
      ((T*)gptr)->~T();
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
      async(gptr.tid())(delete_fctn<T>,gptr);
      upcxx::wait();
      deallocate<T>(gptr);
    }

}

#endif
