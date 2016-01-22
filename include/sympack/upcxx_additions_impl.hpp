#ifndef _UPCXX_ADDITIONS_IMPL_HPP_
#define _UPCXX_ADDITIONS_IMPL_HPP_



#include "sympack/upcxx_additions.hpp"
#include "sympack/LogFile.hpp"

//extern LogFile * logfileptr; 

namespace upcxx{

    template <typename T> void new_fctn(global_ptr<T> gptr){
      new ((T*)gptr) T();
    }
    template <typename T> void delete_fctn(global_ptr<T> gptr){
      ((T*)gptr)->~T();
    }
    template <typename T> global_ptr<T>& Create(int where=-1){
      if(where==-1){where = upcxx::myrank();}
      //currently it only works with local assignments
      global_ptr<T> gptr = allocate<T>(where,1);
      async(where)(new_fctn<T>,gptr);
      upcxx::wait();
      return gptr;
    }
    template <typename T> void Destroy(global_ptr<T>& gptr){
      //currently it only works with local assignments
      int where = gptr.tid();
      async(where)(delete_fctn<T>,gptr);
      upcxx::wait();
      deallocate<T>(gptr);
    }

  template <typename T> inline void ldacopy(int m, int n, global_ptr<T> A, int lda, global_ptr<T> B, int ldb){
    for(int j=0;j<n;j++){
      copy<T>(A + j*lda,B+j*ldb,m);
    }
  }


  template <typename T> inline void ldacopy_async(int m, int n, global_ptr<T> A, int lda, global_ptr<T> B, int ldb, event * e = NULL){
    for(int j=0;j<n;j++){
      async_copy<T>(A + j*lda,B+j*ldb,m,&e[j]);
    }
  }


}

#endif
