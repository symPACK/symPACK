/* Wrapper functions for each cuBLAS call
 * Author: Julian Bellavita, UC Berkeley
 */

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cublas_v2.h"

namespace symPACK {
    /* Test */
    void test_kernel(int);

    /* TRSM wrappers */
    void cublas_strsm(cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t,
                        Int, Int, const float *, const float *, Int, float *, Int);
    
    void cublas_dtrsm(cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t,
                        Int, Int, const double *, const double *, Int, double *, Int);
    
    void cublas_ctrsm(cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t,
                        Int, Int, const cuComplex *, const cuComplex *, Int, cuComplex *, Int);
    
    void cublas_ztrsm(cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t,
                        Int, Int, const cuDoubleComplex *, const cuDoubleComplex *, Int, cuDoubleComplex *, Int);
}