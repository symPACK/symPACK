

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

/* GEMV */
cublasStatus_t cublas_gemv(cublasHandle_t handle, cublasOperation_t trans,
                           int m, int n,
                           const float           *alpha,
                           const float           *A, int lda,
                           const float           *x, int incx,
                           const float           *beta,
                           float           *y, int incy) {
    return cublasSgemv(handle, trans,
                        m, n,
                        alpha,
                        A, lda,
                        x, incx,
                        beta,
                        y, incy);
}

cublasStatus_t cublas_gemv(cublasHandle_t handle, cublasOperation_t trans,
                           int m, int n,
                           const double           *alpha,
                           const double           *A, int lda,
                           const double           *x, int incx,
                           const double           *beta,
                           double           *y, int incy) {
    return cublasDgemv(handle, trans,
                        m, n,
                        alpha,
                        A, lda,
                        x, incx,
                        beta,
                        y, incy);
}

cublasStatus_t cublas_gemv(cublasHandle_t handle, cublasOperation_t trans,
                           int m, int n,
                           const cuComplex           *alpha,
                           const cuComplex           *A, int lda,
                           const cuComplex           *x, int incx,
                           const cuComplex           *beta,
                           cuComplex           *y, int incy) {
    return cublasCgemv(handle, trans,
                        m, n,
                        alpha,
                        A, lda,
                        x, incx,
                        beta,
                        y, incy);
}

cublasStatus_t cublas_gemv(cublasHandle_t handle, cublasOperation_t trans,
                           int m, int n,
                           const cuDoubleComplex           *alpha,
                           const cuDoubleComplex           *A, int lda,
                           const cuDoubleComplex           *x, int incx,
                           const cuDoubleComplex           *beta,
                           cuDoubleComplex           *y, int incy) {
    return cublasZgemv(handle, trans,
                        m, n,
                        alpha,
                        A, lda,
                        x, incx,
                        beta,
                        y, incy);
}

/* GER */
cublasStatus_t  cublas_ger(cublasHandle_t handle, int m, int n,
                           const float           *alpha,
                           const float           *x, int incx,
                           const float           *y, int incy,
                           float           *A, int lda) {
    return cublasSger(handle,
                      m, n,
                      alpha,
                      x, incx,
                      y, incy,
                      A, lda);
}

cublasStatus_t  cublas_ger(cublasHandle_t handle, int m, int n,
                           const double           *alpha,
                           const double           *x, int incx,
                           const double           *y, int incy,
                           double           *A, int lda) {
    return cublasDger(handle,
                      m, n,
                      alpha,
                      x, incx,
                      y, incy,
                      A, lda);
}

/* GERU */
cublasStatus_t  cublas_geru(cublasHandle_t handle, int m, int n,
                           const cuComplex           *alpha,
                           const cuComplex           *x, int incx,
                           const cuComplex           *y, int incy,
                           cuComplex           *A, int lda) {
    return cublasCgeru(handle,
                      m, n,
                      alpha,
                      x, incx,
                      y, incy,
                      A, lda);
}

cublasStatus_t  cublas_geru(cublasHandle_t handle, int m, int n,
                           const cuDoubleComplex           *alpha,
                           const cuDoubleComplex           *x, int incx,
                           const cuDoubleComplex           *y, int incy,
                           cuDoubleComplex           *A, int lda) {
    return cublasZgeru(handle,
                      m, n,
                      alpha,
                      x, incx,
                      y, incy,
                      A, lda);
}

cublasStatus_t  cublas_gerc(cublasHandle_t handle, int m, int n,
                           const cuDoubleComplex           *alpha,
                           const cuDoubleComplex           *x, int incx,
                           const cuDoubleComplex           *y, int incy,
                           cuDoubleComplex           *A, int lda) {
    return cublasZgerc(handle,
                      m, n,
                      alpha,
                      x, incx,
                      y, incy,
                      A, lda);
}

cublasStatus_t  cublas_gerc(cublasHandle_t handle, int m, int n,
                           const cuComplex           *alpha,
                           const cuComplex           *x, int incx,
                           const cuComplex           *y, int incy,
                           cuComplex           *A, int lda) {
    return cublasCgerc(handle,
                      m, n,
                      alpha,
                      x, incx,
                      y, incy,
                      A, lda);
}