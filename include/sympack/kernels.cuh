#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

namespace symPACK {
namespace cudaKernels {
    
    void update_tgt_wrapper(int ,
        int , int ,
        int * , int ,
        int * , 
        double * , double * );
}
}