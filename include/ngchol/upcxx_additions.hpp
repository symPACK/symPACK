#ifndef _UPCXX_ADDITIONS_HPP_
#define _UPCXX_ADDITIONS_HPP_


#include <upcxx.h>

namespace upcxx{
    template <typename T> global_ptr<T>& Create(int where=MYTHREAD);
    template <typename T> void Destroy(global_ptr<T>& gptr);
}

#include "ngchol/upcxx_additions_impl.hpp"

#endif
