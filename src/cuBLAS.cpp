

#include "sympack/cuBLAS.hpp"


 /* ===== LEVEL 1 BLAS ===== */

/* AXPY */
cublasStatus_t cublas_axpy(cublasHandle_t handle, Int n,
                           const float           *alpha,
                           const float           *x, Int incx,
                           float                 *y, Int incy) {
    return cublasSaxpy(handle, n,
                           alpha,
                           x, incx,
                           y, incy);
}

cublasStatus_t cublas_axpy(cublasHandle_t handle, Int n,
                           const double           *alpha,
                           const double           *x, Int incx,
                           double                 *y, Int incy) {
    return cublasDaxpy(handle, n,
                           alpha
                           x, incx,
                           y, incy);
}

cublasStatus_t cublas_axpy(cublasHandle_t handle, Int n,
                           const cuComplex           *alpha,
                           const cuComplex           *x, Int incx,
                           cuComplex                 *y, Int incy) {
    return cublasCaxpy(handle, n,
                           alpha,
                           x, incx,
                           y, incy);
}

cublasStatus_t cublas_axpy(cublasHandle_t handle, Int n,
                           const cuDoubleComplex           *alpha,
                           const cuDoubleComplex           *x, Int incx,
                           cuDoubleComplex                 *y, Int incy) {
    return cublasZaxpy(handle, n,
                           alpha,
                           x, incx,
                           y, incy);
}


/* COPY */
cublasStatus_t cublas_copy(cublasHandle_t handle, int n,
                           const float           *x, int incx,
                           float                 *y, int incy) {
    return cublasScopy(handle, n, 
                        x, incx,
                        y, incy);
}

cublasStatus_t cublas_copy(cublasHandle_t handle, int n,
                           const double           *x, int incx,
                           double                 *y, int incy) {
    return cublasDcopy(handle, n, 
                        x, incx,
                        y, incy);
}

cublasStatus_t cublas_copy(cublasHandle_t handle, int n,
                           const cuComplex           *x, int incx,
                           cuComplex                 *y, int incy) {
    return cublasCcopy(handle, n, 
                        x, incx,
                        y, incy);
}

cublasStatus_t cublas_copy(cublasHandle_t handle, int n,
                           const cuDoubleComplex           *x, int incx,
                           cuDoubleComplex                 *y, int incy) {
    return cublasZcopy(handle, n, 
                        x, incx,
                        y, incy);
}


/* SCAL */
cublasStatus_t  cublas_scal(cublasHandle_t handle, int n,
                            const float           *alpha,
                            float           *x, int incx) {
    return cublasSscal(handle, n,
                        alpha, 
                        x, incx);
}

cublasStatus_t  cublas_scal(cublasHandle_t handle, int n,
                            const double           *alpha,
                            double           *x, int incx) {
    return cublasDscal(handle, n,
                        alpha, 
                        x, incx);
}

cublasStatus_t  cublas_scal(cublasHandle_t handle, int n,
                            const cuComplex           *alpha,
                            cuComplex           *x, int incx) {
    return cublasCscal(handle, n,
                        alpha, 
                        x, incx);
}

cublasStatus_t  cublas_scal(cublasHandle_t handle, int n,
                            const float           *alpha,
                            cuComplex           *x, int incx) {
    return cublasCsscal(handle, n,
                        alpha, 
                        x, incx);
}

cublasStatus_t  cublas_scal(cublasHandle_t handle, int n,
                            const cuDoubleComplex           *alpha,
                            cuDoubleComplex           *x, int incx) {
    return cublasZscal(handle, n,
                        alpha, 
                        x, incx);
}

cublasStatus_t  cublas_scal(cublasHandle_t handle, int n,
                            const double           *alpha,
                            cuDoubleComplex           *x, int incx) {
    return cublasZdscal(handle, n,
                        alpha, 
                        x, incx);
}


/* ===== LEVEL 2 BLAS ===== */