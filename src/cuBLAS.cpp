
#ifdef CUDA_MODE
#include "sympack/cuBLAS.hpp"
#include "sympack/Environment.hpp"

namespace symPACK {
namespace cublas {
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
                           alpha,
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
cublasStatus_t cublas_copy(cublasHandle_t handle, Int n,
                           const float           *x, Int incx,
                           float                 *y, Int incy) {
    return cublasScopy(handle, n, 
                        x, incx,
                        y, incy);
}

cublasStatus_t cublas_copy(cublasHandle_t handle, Int n,
                           const double           *x, Int incx,
                           double                 *y, Int incy) {
    return cublasDcopy(handle, n, 
                        x, incx,
                        y, incy);
}

cublasStatus_t cublas_copy(cublasHandle_t handle, Int n,
                           const cuComplex           *x, Int incx,
                           cuComplex                 *y, Int incy) {
    return cublasCcopy(handle, n, 
                        x, incx,
                        y, incy);
}

cublasStatus_t cublas_copy(cublasHandle_t handle, Int n,
                           const cuDoubleComplex           *x, Int incx,
                           cuDoubleComplex                 *y, Int incy) {
    return cublasZcopy(handle, n, 
                        x, incx,
                        y, incy);
}


/* SCAL */
cublasStatus_t  cublas_scal(cublasHandle_t handle, Int n,
                            const float           *alpha,
                            float           *x, Int incx) {
    return cublasSscal(handle, n,
                        alpha, 
                        x, incx);
}

cublasStatus_t  cublas_scal(cublasHandle_t handle, Int n,
                            const double           *alpha,
                            double           *x, Int incx) {
    return cublasDscal(handle, n,
                        alpha, 
                        x, incx);
}

cublasStatus_t  cublas_scal(cublasHandle_t handle, Int n,
                            const cuComplex           *alpha,
                            cuComplex           *x, Int incx) {
    return cublasCscal(handle, n,
                        alpha, 
                        x, incx);
}

cublasStatus_t  cublas_scal(cublasHandle_t handle, Int n,
                            const float           *alpha,
                            cuComplex           *x, Int incx) {
    return cublasCsscal(handle, n,
                        alpha, 
                        x, incx);
}

cublasStatus_t  cublas_scal(cublasHandle_t handle, Int n,
                            const cuDoubleComplex           *alpha,
                            cuDoubleComplex           *x, Int incx) {
    return cublasZscal(handle, n,
                        alpha, 
                        x, incx);
}

cublasStatus_t  cublas_scal(cublasHandle_t handle, Int n,
                            const double           *alpha,
                            cuDoubleComplex           *x, Int incx) {
    return cublasZdscal(handle, n,
                        alpha, 
                        x, incx);
}


/* ===== LEVEL 2 BLAS ===== */

/* GEMV */
cublasStatus_t cublas_gemv(cublasHandle_t handle, cublasOperation_t trans,
                           Int m, Int n,
                           const float           *alpha,
                           const float           *A, Int lda,
                           const float           *x, Int incx,
                           const float           *beta,
                           float           *y, Int incy) {
    return cublasSgemv(handle, trans,
                        m, n,
                        alpha,
                        A, lda,
                        x, incx,
                        beta,
                        y, incy);
}

cublasStatus_t cublas_gemv(cublasHandle_t handle, cublasOperation_t trans,
                           Int m, Int n,
                           const double           *alpha,
                           const double           *A, Int lda,
                           const double           *x, Int incx,
                           const double           *beta,
                           double           *y, Int incy) {
    return cublasDgemv(handle, trans,
                        m, n,
                        alpha,
                        A, lda,
                        x, incx,
                        beta,
                        y, incy);
}

cublasStatus_t cublas_gemv(cublasHandle_t handle, cublasOperation_t trans,
                           Int m, Int n,
                           const cuComplex           *alpha,
                           const cuComplex           *A, Int lda,
                           const cuComplex           *x, Int incx,
                           const cuComplex           *beta,
                           cuComplex           *y, Int incy) {
    return cublasCgemv(handle, trans,
                        m, n,
                        alpha,
                        A, lda,
                        x, incx,
                        beta,
                        y, incy);
}

cublasStatus_t cublas_gemv(cublasHandle_t handle, cublasOperation_t trans,
                           Int m, Int n,
                           const cuDoubleComplex           *alpha,
                           const cuDoubleComplex           *A, Int lda,
                           const cuDoubleComplex           *x, Int incx,
                           const cuDoubleComplex           *beta,
                           cuDoubleComplex           *y, Int incy) {
    return cublasZgemv(handle, trans,
                        m, n,
                        alpha,
                        A, lda,
                        x, incx,
                        beta,
                        y, incy);
}

/* GER */
cublasStatus_t  cublas_ger(cublasHandle_t handle, Int m, Int n,
                           const float           *alpha,
                           const float           *x, Int incx,
                           const float           *y, Int incy,
                           float           *A, Int lda) {
    return cublasSger(handle,
                      m, n,
                      alpha,
                      x, incx,
                      y, incy,
                      A, lda);
}

cublasStatus_t  cublas_ger(cublasHandle_t handle, Int m, Int n,
                           const double           *alpha,
                           const double           *x, Int incx,
                           const double           *y, Int incy,
                           double           *A, Int lda) {
    return cublasDger(handle,
                      m, n,
                      alpha,
                      x, incx,
                      y, incy,
                      A, lda);
}

/* GERU */
cublasStatus_t  cublas_geru(cublasHandle_t handle, Int m, Int n,
                           const cuComplex           *alpha,
                           const cuComplex           *x, Int incx,
                           const cuComplex           *y, Int incy,
                           cuComplex           *A, Int lda) {
    return cublasCgeru(handle,
                      m, n,
                      alpha,
                      x, incx,
                      y, incy,
                      A, lda);
}

cublasStatus_t  cublas_geru(cublasHandle_t handle, Int m, Int n,
                           const cuDoubleComplex           *alpha,
                           const cuDoubleComplex           *x, Int incx,
                           const cuDoubleComplex           *y, Int incy,
                           cuDoubleComplex           *A, Int lda) {
    return cublasZgeru(handle,
                      m, n,
                      alpha,
                      x, incx,
                      y, incy,
                      A, lda);
}

cublasStatus_t  cublas_gerc(cublasHandle_t handle, Int m, Int n,
                           const cuDoubleComplex           *alpha,
                           const cuDoubleComplex           *x, Int incx,
                           const cuDoubleComplex           *y, Int incy,
                           cuDoubleComplex           *A, Int lda) {
    return cublasZgerc(handle,
                      m, n,
                      alpha,
                      x, incx,
                      y, incy,
                      A, lda);
}

cublasStatus_t  cublas_gerc(cublasHandle_t handle, Int m, Int n,
                           const cuComplex           *alpha,
                           const cuComplex           *x, Int incx,
                           const cuComplex           *y, Int incy,
                           cuComplex           *A, Int lda) {
    return cublasCgerc(handle,
                      m, n,
                      alpha,
                      x, incx,
                      y, incy,
                      A, lda);
}


/* ===== LEVEL 3 BLAS ===== */

/* SYRK */

cublasStatus_t cublas_syrk(cublasHandle_t handle,
                           cublasFillMode_t uplo, cublasOperation_t trans,
                           Int n, Int k,
                           const float           *alpha,
                           const float           *A, Int lda,
                           const float           *beta,
                           float           *C, Int ldc) {
    
    return cublasSsyrk(handle, uplo, trans,
                        n, k,
                        alpha, A, lda,
                        beta, C, ldc);

}

cublasStatus_t cublas_syrk(cublasHandle_t handle,
                           cublasFillMode_t uplo, cublasOperation_t trans,
                           Int n, Int k,
                           const double           *alpha,
                           const double           *A, Int lda,
                           const double           *beta,
                           double           *C, Int ldc) {
    
    return cublasDsyrk(handle, uplo, trans,
                        n, k,
                        alpha, A, lda,
                        beta, C, ldc);

}

cublasStatus_t cublas_syrk(cublasHandle_t handle,
                           cublasFillMode_t uplo, cublasOperation_t trans,
                           Int n, Int k,
                           const cuComplex           *alpha,
                           const cuComplex           *A, Int lda,
                           const cuComplex           *beta,
                           cuComplex           *C, Int ldc) {
    
    return cublasCsyrk(handle, uplo, trans,
                        n, k,
                        alpha, A, lda,
                        beta, C, ldc);

}

cublasStatus_t cublas_syrk(cublasHandle_t handle,
                           cublasFillMode_t uplo, cublasOperation_t trans,
                           Int n, Int k,
                           const cuDoubleComplex           *alpha,
                           const cuDoubleComplex           *A, Int lda,
                           const cuDoubleComplex           *beta,
                           cuDoubleComplex           *C, Int ldc) {
    
    return cublasZsyrk(handle, uplo, trans,
                        n, k,
                        alpha, A, lda,
                        beta, C, ldc);

}

/* GEMM */
cublasStatus_t cublas_gemm(cublasHandle_t handle,
                           cublasOperation_t transa, cublasOperation_t transb,
                           Int m, Int n, Int k,
                           const float           *alpha,
                           const float           *A, Int lda,
                           const float           *B, Int ldb,
                           const float           *beta,
                           float           *C, Int ldc) {
    
    return cublasSgemm(handle, transa, transb,
                        m, n, k,
                        alpha, A, lda,
                        B, ldb, beta,
                        C, ldc);

}

cublasStatus_t cublas_gemm(cublasHandle_t handle,
                           cublasOperation_t transa, cublasOperation_t transb,
                           Int m, Int n, Int k,
                           const double          *alpha,
                           const double           *A, Int lda,
                           const double           *B, Int ldb,
                           const double           *beta,
                           double           *C, Int ldc) {
    
    return cublasDgemm(handle, transa, transb,
                        m, n, k,
                        alpha, A, lda,
                        B, ldb, beta,
                        C, ldc);

}

cublasStatus_t cublas_gemm(cublasHandle_t handle,
                           cublasOperation_t transa, cublasOperation_t transb,
                           Int m, Int n, Int k,
                           const cuComplex          *alpha,
                           const cuComplex           *A, Int lda,
                           const cuComplex           *B, Int ldb,
                           const cuComplex           *beta,
                           cuComplex           *C, Int ldc) {
    
    return cublasCgemm(handle, transa, transb,
                        m, n, k,
                        alpha, A, lda,
                        B, ldb, beta,
                        C, ldc);

}

cublasStatus_t cublas_gemm(cublasHandle_t handle,
                           cublasOperation_t transa, cublasOperation_t transb,
                           Int m, Int n, Int k,
                           const cuDoubleComplex          *alpha,
                           const cuDoubleComplex           *A, Int lda,
                           const cuDoubleComplex           *B, Int ldb,
                           const cuDoubleComplex           *beta,
                           cuDoubleComplex           *C, Int ldc) {
    
    return cublasZgemm(handle, transa, transb,
                        m, n, k,
                        alpha, A, lda,
                        B, ldb, beta,
                        C, ldc);

}

cublasStatus_t cublas_gemm(cublasHandle_t handle,
                           cublasOperation_t transa, cublasOperation_t transb,
                           Int m, Int n, Int k,
                           const __half          *alpha,
                           const __half           *A, Int lda,
                           const __half           *B, Int ldb,
                           const __half           *beta,
                           __half           *C, Int ldc) {
    
    return cublasHgemm(handle, transa, transb,
                        m, n, k,
                        alpha, A, lda,
                        B, ldb, beta,
                        C, ldc);

}

/* TRSM */
cublasStatus_t cublas_trsm(cublasHandle_t handle,
                           cublasSideMode_t side, cublasFillMode_t uplo,
                           cublasOperation_t trans, cublasDiagType_t diag,
                           Int m, Int n,
                           const float           *alpha,
                           const float           *A, Int lda,
                           float           *B, Int ldb) {
    
    return cublasStrsm(handle, side, uplo, trans, diag,
                        m, n,
                        alpha, A, lda,
                        B, ldb);

}

cublasStatus_t cublas_trsm(cublasHandle_t handle,
                           cublasSideMode_t side, cublasFillMode_t uplo,
                           cublasOperation_t trans, cublasDiagType_t diag,
                           Int m, Int n,
                           const double           *alpha,
                           const double           *A, Int lda,
                           double           *B, Int ldb) {
    
    return cublasDtrsm(handle, side, uplo, trans, diag,
                        m, n,
                        alpha, A, lda,
                        B, ldb);

}

cublasStatus_t cublas_trsm(cublasHandle_t handle,
                           cublasSideMode_t side, cublasFillMode_t uplo,
                           cublasOperation_t trans, cublasDiagType_t diag,
                           Int m, Int n,
                           const cuComplex           *alpha,
                           const cuComplex           *A, Int lda,
                           cuComplex           *B, Int ldb) {
    
    return cublasCtrsm(handle, side, uplo, trans, diag,
                        m, n,
                        alpha, A, lda,
                        B, ldb);

}

cublasStatus_t cublas_trsm(cublasHandle_t handle,
                           cublasSideMode_t side, cublasFillMode_t uplo,
                           cublasOperation_t trans, cublasDiagType_t diag,
                           Int m, Int n,
                           const cuDoubleComplex           *alpha,
                           const cuDoubleComplex           *A, Int lda,
                           cuDoubleComplex           *B, Int ldb) {
    
    return cublasZtrsm(handle, side, uplo, trans, diag,
                        m, n,
                        alpha, A, lda,
                        B, ldb);

}
}
}
#endif //CUDA_MODE
