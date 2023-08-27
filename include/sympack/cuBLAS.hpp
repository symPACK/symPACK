/* cuBLAS interface */
#include  "sympack/Environment.hpp"
#ifdef CUDA_MODE
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cublas_v2.h"

                                                                                               \

namespace symPACK {
namespace cublas {
    // NOTE: This only contains a subset of the complete BLAS, just the operations used in symPACK are included.
    // SYRK, TRSM, GEMM, Axpy, Scal, Gemv, Geru, Copy
    // A lot of these operations are only used in the defunct SuperNode methods that have been replaced with
    // the methods in symPACKMatrix2D.hpp, but they're here just in case they become needed at some point.
    typedef  int                    Int;
    /* ===== LEVEL 1 BLAS ===== */

    /* AXPY */
    cublasStatus_t cublas_axpy(cublasHandle_t, Int,
                           const float           *,
                           const float           *, Int,
                           float                 *, Int);

    cublasStatus_t cublas_axpy(cublasHandle_t, Int,
                           const double           *,
                           const double           *, Int,
                           double                 *, Int);

    cublasStatus_t cublas_axpy(cublasHandle_t, Int,
                           const cuComplex           *,
                           const cuComplex           *, Int,
                           cuComplex                 *, Int);

    cublasStatus_t cublas_axpy(cublasHandle_t, Int,
                           const cuDoubleComplex           *,
                           const cuDoubleComplex           *, Int,
                           cuDoubleComplex                 *, Int);

    template <typename T>
    cublasStatus_t cublas_axpy_wrapper2(Int N,
                           const T           DA,
                           const T           * DX, Int incx,
                           T                 * DY, Int incy) {
        CUBLAS_ERROR_CHECK(cublas_axpy(symPACK::cublas_handler, N, &DA, DX, incx, DY, incy));
        return CUBLAS_STATUS_SUCCESS;
    }

    /* COPY */
    cublasStatus_t cublas_copy(cublasHandle_t, Int,
                           const float           *, Int,
                           float                 *, Int);

    cublasStatus_t cublas_copy(cublasHandle_t, Int,
                           const double           *, Int,
                           double                 *, Int);

    cublasStatus_t cublas_copy(cublasHandle_t, Int,
                           const cuComplex           *, Int,
                           cuComplex                 *, Int);

    cublasStatus_t cublas_copy(cublasHandle_t, Int,
                           const cuDoubleComplex           *, Int,
                           cuDoubleComplex                 *, Int);

    template <typename T>
    cublasStatus_t cublas_copy_wrapper(Int n,
		   			 T * dx, Int incx,
					 T * dy, Int incy) {
        throw std::runtime_error("cuBLAS copy is not currently implemented");
        return CUBLAS_STATUS_EXECUTION_FAILED;
    }

    /* SCAL */
    cublasStatus_t  cublas_scal(cublasHandle_t , Int,
                            const float           *,
                            float           *, Int);

    cublasStatus_t  cublas_scal(cublasHandle_t , Int,
                            const double           *,
                            double           *, Int);

    cublasStatus_t  cublas_scal(cublasHandle_t , Int,
                            const cuComplex           *,
                            cuComplex           *, Int);

    cublasStatus_t  cublas_scal(cublasHandle_t , Int,
                            const float           *,
                            cuComplex           *, Int);

    cublasStatus_t  cublas_scal(cublasHandle_t , Int,
                            const cuDoubleComplex           *,
                            cuDoubleComplex           *, Int);
    
    cublasStatus_t  cublas_scal(cublasHandle_t , Int,
                            const double           *,
                            cuDoubleComplex           *, Int);

    template <typename T>
    cublasStatus_t  cublas_scal_wrapper2(Int N,
                            const T           DA,
                            T           * DX, Int incx) {
        throw std::runtime_error("cuBLAS scal is not currently implemented");
        return CUBLAS_STATUS_EXECUTION_FAILED;
    }

    /* ===== LEVEL 2 BLAS ===== */

    /* GEMV */
    cublasStatus_t cublas_gemv(cublasHandle_t , cublasOperation_t ,
                           Int, Int,
                           const float           *,
                           const float           *, Int ,
                           const float           *, Int ,
                           const float           *,
                           float           *, Int );

    cublasStatus_t cublas_gemv(cublasHandle_t , cublasOperation_t ,
                           Int, Int,
                           const double           *,
                           const double           *, Int ,
                           const double           *, Int ,
                           const double           *,
                           double           *, Int );

    cublasStatus_t cublas_gemv(cublasHandle_t , cublasOperation_t ,
                           Int, Int,
                           const cuComplex           *,
                           const cuComplex           *, Int ,
                           const cuComplex           *, Int ,
                           const cuComplex           *,
                           cuComplex           *, Int );

    cublasStatus_t cublas_gemv(cublasHandle_t , cublasOperation_t ,
                           Int, Int,
                           const cuDoubleComplex           *,
                           const cuDoubleComplex           *, Int ,
                           const cuDoubleComplex           *, Int ,
                           const cuDoubleComplex           *,
                           cuDoubleComplex           *, Int );
    
    template <typename T>
    cublasStatus_t cublas_gemv_wrapper(cublasOperation_t op,
                           Int M, Int N,
                           const T           alpha,
                           const T           * A, Int lda,
                           const T           * X, Int incx,
                           const T           beta,
                           T           * Y, Int incy) {
        throw std::runtime_error("cuBLAS gemv is not currently implemented");
        return CUBLAS_STATUS_EXECUTION_FAILED;
    };


    /* GER */
    cublasStatus_t  cublas_ger(cublasHandle_t , Int, Int,
                           const float           *,
                           const float           *, Int,
                           const float           *, Int,
                           float           *, Int);

    cublasStatus_t  cublas_ger(cublasHandle_t , Int, Int,
                           const double           *,
                           const double           *, Int,
                           const double           *, Int,
                           double           *, Int);


    /* GERU */
    cublasStatus_t  cublas_geru(cublasHandle_t , Int, Int,
                           const cuComplex           *,
                           const cuComplex           *, Int,
                           const cuComplex           *, Int,
                           cuComplex           *, Int);

    cublasStatus_t  cublas_geru(cublasHandle_t , Int, Int,
                           const cuDoubleComplex           *,
                           const cuDoubleComplex           *, Int,
                           const cuDoubleComplex           *, Int,
                           cuDoubleComplex           *, Int);
    
    
    /* GERC */
    cublasStatus_t  cublas_gerc(cublasHandle_t , Int, Int,
                           const cuDoubleComplex           *,
                           const cuDoubleComplex           *, Int,
                           const cuDoubleComplex           *, Int,
                           cuDoubleComplex           *, Int);

    cublasStatus_t  cublas_gerc(cublasHandle_t , Int, Int,
                           const cuComplex           *,
                           const cuComplex           *, Int,
                           const cuComplex           *, Int,
                           cuComplex           *, Int);


    template <typename T>
    cublasStatus_t cublas_ger_wrapper(Int M, Int N,
                           const T           alpha,
                           const T           * X, Int incx,
                           const T           * Y, Int incy,
                           T           * A, Int lda) {

        throw std::runtime_error("cuBLAS ger is not currently implemented");
        return CUBLAS_STATUS_EXECUTION_FAILED;
    }

    /* ===== LEVEL 3 BLAS ===== */
    
    /* SYRK */
    cublasStatus_t cublas_syrk(cublasHandle_t,
                           cublasFillMode_t, cublasOperation_t,
                           Int, Int,
                           const float           *,
                           const float           *, Int ,
                           const float           *,
                           float           *, Int);

    cublasStatus_t cublas_syrk(cublasHandle_t,
                           cublasFillMode_t, cublasOperation_t,
                           Int, Int,
                           const double           *,
                           const double           *, Int ,
                           const double           *,
                           double           *, Int);

    cublasStatus_t cublas_syrk(cublasHandle_t,
                           cublasFillMode_t, cublasOperation_t,
                           Int, Int,
                           const cuComplex           *,
                           const cuComplex           *, Int ,
                           const cuComplex           *,
                           cuComplex           *, Int);

    cublasStatus_t cublas_syrk(cublasHandle_t,
                           cublasFillMode_t, cublasOperation_t,
                           Int, Int,
                           const cuDoubleComplex           *,
                           const cuDoubleComplex           *, Int ,
                           const cuDoubleComplex           *,
                           cuDoubleComplex           *, Int);
    
    template <typename T>
    cublasStatus_t cublas_syrk_wrapper2(cublasFillMode_t uplo, cublasOperation_t trans,
                                        Int N, Int K,
                                        const T alpha,
                                        const T * A, Int lda,
                                        const T beta,
                                        T * C, Int ldc
                                        ) {
        

        CUBLAS_ERROR_CHECK(cublas_syrk(symPACK::cublas_handler, 
                                       uplo, trans,
                                       N, K,
                                       &alpha,
                                       A, lda,
                                       &beta,
                                       C, ldc));

        return CUBLAS_STATUS_SUCCESS;
    }


    /* GEMM */
    cublasStatus_t cublas_gemm(cublasHandle_t,
                           cublasOperation_t, cublasOperation_t,
                           Int, Int, Int,
                           const float           *,
                           const float           *, Int,
                           const float           *, Int,
                           const float           *,
                           float           *, Int);

    cublasStatus_t cublas_gemm(cublasHandle_t,
                           cublasOperation_t, cublasOperation_t,
                           Int, Int, Int,
                           const double           *,
                           const double           *, Int,
                           const double           *, Int,
                           const double           *,
                           double           *, Int);
    
    cublasStatus_t cublas_gemm(cublasHandle_t,
                           cublasOperation_t, cublasOperation_t,
                           Int, Int, Int,
                           const cuComplex           *,
                           const cuComplex           *, Int,
                           const cuComplex           *, Int,
                           const cuComplex           *,
                           cuComplex           *, Int);

    cublasStatus_t cublas_gemm(cublasHandle_t,
                           cublasOperation_t, cublasOperation_t,
                           Int, Int, Int,
                           const cuDoubleComplex           *,
                           const cuDoubleComplex           *, Int,
                           const cuDoubleComplex           *, Int,
                           const cuDoubleComplex           *,
                           cuDoubleComplex           *, Int);

    cublasStatus_t cublas_gemm(cublasHandle_t,
                           cublasOperation_t, cublasOperation_t,
                           Int, Int, Int,
                           const __half           *,
                           const __half           *, Int,
                           const __half           *, Int,
                           const __half           *,
                           __half           *, Int);


    template <typename T>
    cublasStatus_t cublas_gemm_wrapper2(
                           cublasOperation_t opA, cublasOperation_t opB,
                           Int M, Int N, Int K,
                           const T           alpha,
                           const T           * A , Int lda,
                           const T          * B, Int ldb,
                           const T           beta,
                           T           * C, Int ldc) {

        CUBLAS_ERROR_CHECK(cublas_gemm(symPACK::cublas_handler, opA, opB,
                    M, N, K,
                    &alpha, A, lda,
                    B, ldb, &beta,
                    C, ldc));
        return CUBLAS_STATUS_SUCCESS;                     
    }

    /* TRSM */
    cublasStatus_t cublas_trsm(cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t,
                        Int, Int, const float *, const float *, Int, float *, Int);
    
    cublasStatus_t cublas_trsm(cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t,
                        Int, Int, const double *, const double *, Int, double *, Int);
    
    cublasStatus_t cublas_trsm(cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t,
                        Int, Int, const cuComplex *, const cuComplex *, Int, cuComplex *, Int);
    
    cublasStatus_t cublas_trsm(cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t,
                        Int, Int, const cuDoubleComplex *, const cuDoubleComplex *, Int, cuDoubleComplex *, Int);
    

    template <typename T>
    cublasStatus_t cublas_trsm_wrapper2(cublasSideMode_t side, cublasFillMode_t fill, cublasOperation_t op, cublasDiagType_t diag,
                        Int M , Int N, 
                        const T alpha, T * A, Int lda, 
                        T * B, Int ldb) {
    

        CUBLAS_ERROR_CHECK(cublas_trsm(symPACK::cublas_handler, side, fill, op, diag, 
                    M, N, 
                    &alpha, A, lda, 
                    B, ldb));      

        return CUBLAS_STATUS_SUCCESS;              

    }

}
}
#endif //CUDA_MODE

