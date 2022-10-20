/* cuBLAS interface
 * Author: Julian Bellavita, UC Berkeley
 */
#include  "sympack/Environment.hpp"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cublas_v2.h"

#define CUBLAS_ERROR_CHECK(s)                                                                       \
{                                                                                                   \
    cublasStatus_t err;                                                                             \
    if ((err = (s)) != CUBLAS_STATUS_SUCCESS)                                                       \
    {                                                                                               \
        std::cout << "cuBLAS Error " << err << " at " << __FILE__ << ":" << __LINE__ << "\n";       \
        exit(1);                                                                                    \
    }                                                                                               \
}                                                                                                   \

#define CUDA_ERROR_CHECK(s)                                                                         \
{                                                                                                   \
    cudaError_t error = s;                                                                          \
    if (error != cudaSuccess) {                                                                     \
        std::cout << "CUDA Error " << error << " at " << __FILE__ << ":" << __LINE__ << "\n";       \
        std::cout << cudaGetErrorString(error) << "\n";                                             \
        exit(1);                                                                                    \
    }                                                                                               \
}                                                                                                   \

namespace symPACK {
namespace cublas {
    // NOTE: This only contains a subset of the complete BLAS, just the operations used in symPACK are included.
    // SYRK, TRSM, GEMM, Axpy, Scal, Gemv, Geru, Copy
    void test(int);
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
                           const T           * Y, Int incy) {
        
        logfileptr->OFS()<<"DOING GEMV\n";
        
        /* Device buffers */
        T *d_A;
        T *d_X;
        T *d_Y;

        cublasStatus_t status;

        /* Set up device buffers */
        CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_A), lda * N * sizeof(A[0])));
        CUBLAS_ERROR_CHECK(cublasSetMatrix(lda, N, sizeof(A[0]), A, lda, d_A, lda));

        long dimx;
        if (op==CUBLAS_OP_T) 
            dimx = (1 + (M - 1) * abs(incx));
        else 
            dimx = (1 + (N - 1) * abs(incx));
    
        CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_X), dimx));
        CUBLAS_ERROR_CHECK(cublasSetVector(dimx, sizeof(X[0]), X, incx, d_X, incx));

        long dimy;
        if (op==CUBLAS_OP_T) 
            dimy = (1 + (N - 1) * abs(incy));
        else 
            dimy = (1 + (M - 1) * abs(incy));
        
        CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_Y), dimy));
        CUBLAS_ERROR_CHECK(cublasSetVector(dimy, sizeof(Y[0]), Y, incy, d_Y, incy));

        /* do gemv */
        CUBLAS_ERROR_CHECK(cublas_gemv(symPACK::handler, op,
                    M, N,
                    &alpha, d_A, lda,
                    d_X, incx,
                    &beta,
                    d_Y, incy));

        /* copy result to host */
        CUBLAS_ERROR_CHECK(cublasGetVector(dimy, sizeof(Y[0]), d_Y, incy, Y, incy));

        /* Cleanup */
        CUDA_ERROR_CHECK(cudaFree(d_A));
        CUDA_ERROR_CHECK(cudaFree(d_Y));
        CUDA_ERROR_CHECK(cudaFree(d_X));

        status = CUBLAS_STATUS_SUCCESS;
        return status;
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
    cublasStatus_t cublas_gemm_wrapper(cublasHandle_t handler,
                           cublasOperation_t opA, cublasOperation_t opB,
                           Int M, Int N, Int K,
                           const T           alpha,
                           const T           * A , Int lda,
                           const T          * B, Int ldb,
                           const T           beta,
                           T           * C, Int ldc) {
        

        T * d_A;
        T * d_B;
        T * d_C;

        logfileptr->OFS()<<"DOING GEMM\n";
        if (opA==CUBLAS_OP_N) {
            //A is lda*K
            CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_A), lda * K * sizeof(A[0])));
            CUBLAS_ERROR_CHECK(cublasSetMatrix(lda, K, sizeof(A[0]), A, lda, d_A, lda));
        } else if (opA==CUBLAS_OP_T) {
            //A is lda*M
            CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_A), lda * M * sizeof(A[0])));
            CUBLAS_ERROR_CHECK(cublasSetMatrix(lda, M, sizeof(A[0]), A, lda, d_A, lda));
        }

        if (opB==CUBLAS_OP_N) {
            //B is ldb*N
            CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_B), ldb * N * sizeof(B[0])));
            CUBLAS_ERROR_CHECK(cublasSetMatrix(ldb, N, sizeof(B[0]), B, ldb, d_B, ldb));
        } else if (opB==CUBLAS_OP_T) {
            //B is lda*K
            CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_B), ldb * K * sizeof(B[0])));
            CUBLAS_ERROR_CHECK(cublasSetMatrix(ldb, K, sizeof(B[0]), B, ldb, d_B, ldb));
        }

        CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_C), ldc * N * sizeof(C[0])));
        CUBLAS_ERROR_CHECK(cublasSetMatrix(ldc, N, sizeof(C[0]), C, ldc, d_C, ldc));

        /* Do GEMM */
        CUBLAS_ERROR_CHECK(cublas_gemm(symPACK::handler, opA, opB,
                    M, N, K,
                    alpha, d_A, lda,
                    d_B, ldb, beta,
                    d_C, ldc));

        /* Copy matrices to host */
        CUBLAS_ERROR_CHECK(cublasGetMatrix(ldc, N, sizeof(C[0]), d_C, ldc, C, ldc));

        /* Cleanup */
        CUDA_ERROR_CHECK(cudaFree(d_A));
        CUDA_ERROR_CHECK(cudaFree(d_B));
        CUDA_ERROR_CHECK(cudaFree(d_C));

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
    cublasStatus_t cublas_trsm_wrapper(cublasHandle_t handler, cublasSideMode_t side, cublasFillMode_t fill, cublasOperation_t op, cublasDiagType_t diag,
                        Int M , Int N, const T alpha, T * A, Int lda, T * B, Int ldb) {
        
        T * d_A;
        T * d_B;
        
        logfileptr->OFS()<<"DOING TRSM\n";
        /* Setup device matrices */
        if (side==CUBLAS_SIDE_LEFT) {
            CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_A), lda * M * sizeof(A[0])));
            CUBLAS_ERROR_CHECK(cublasSetMatrix(lda, M, sizeof(A[0]), A, lda, d_A, lda));
        } else {
            CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_A), lda * N * sizeof(A[0])));
            CUBLAS_ERROR_CHECK(cublasSetMatrix(lda, N, sizeof(A[0]), A, lda, d_A, lda));
        }
        CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_B), ldb * N * sizeof(B[0])));
        CUBLAS_ERROR_CHECK(cublasSetMatrix(ldb, N, sizeof(B[0]), B, ldb, d_B, ldb));
                
        /* Do TRSM */
        CUBLAS_ERROR_CHECK(cublas_trsm(symPACK::handler, side, fill, op, diag, 
                    M, N, 
                    &alpha, d_A, lda, 
                    d_B, ldb));

        /* Copy matrices to host */
        CUBLAS_ERROR_CHECK(cublasGetMatrix(ldb, N, sizeof(B[0]), d_B, ldb, B, ldb));

        /* Cleanup */
        CUDA_ERROR_CHECK(cudaFree(d_A));
        CUDA_ERROR_CHECK(cudaFree(d_B));

        return CUBLAS_STATUS_SUCCESS;
    }

    /* DO_CUBLAS_L3 */

    enum l3_cublas_ops{OP_SYRK,
                       OP_GEMM,
                       OP_TRSM};
    template <class T>
    inline void do_cublas_l3(l3_cublas_ops op, cublasSideMode_t side, cublasFillMode_t fill, 
                                  cublasOperation_t opA, cublasOperation_t opB,
                                  cublasDiagType_t diag, 
                                  Int M, Int N, Int K,
                                  T *alpha, T *A, Int lda, 
                                  T *beta, T *B, Int ldb, 
                                  T *C, Int ldc) {
        T * d_A;
        T * d_B;
        T * d_C;
        
        /* Create handlers for each L3 op */
        switch(op) {
    
            case OP_TRSM:
                logfileptr->OFS()<<"DOING TRSM\n";
                /* Setup device matrices */
                if (side==CUBLAS_SIDE_LEFT) {
                    cudaMalloc(reinterpret_cast<void **>(&d_A), lda * M * sizeof(A[0]));
                    cublasSetMatrix(lda, M, sizeof(A[0]), A, lda, d_A, lda);
                } else {
                    cudaMalloc(reinterpret_cast<void **>(&d_A), lda * N * sizeof(A[0]));
                    cublasSetMatrix(lda, N, sizeof(A[0]), A, lda, d_A, lda);
                }
                cudaMalloc(reinterpret_cast<void **>(&d_B), ldb * N * sizeof(B[0]));
                cublasSetMatrix(ldb, N, sizeof(B[0]), B, ldb, d_B, ldb);
                
                /* Do TRSM */
                cublas_trsm(symPACK::handler, side, fill, opA, diag, 
                            M, N, 
                            alpha, d_A, lda, 
                            d_B, ldb);

                /* Copy matrices to host */
                cublasGetMatrix(ldb, N, sizeof(B[0]), d_B, ldb, B, ldb);

                /* Cleanup */
                cudaFree(d_A);
                cudaFree(d_B);
                break;
            
            case OP_GEMM:
                logfileptr->OFS()<<"DOING GEMM\n";
                if (opA==CUBLAS_OP_N) {
                    //A is lda*K
                    cudaMalloc(reinterpret_cast<void**>(&d_A), lda * K * sizeof(A[0]));
                    cublasSetMatrix(lda, K, sizeof(A[0]), A, lda, d_A, lda);
                } else if (opA==CUBLAS_OP_T) {
                    //A is lda*M
                    cudaMalloc(reinterpret_cast<void**>(&d_A), lda * M * sizeof(A[0]));
                    cublasSetMatrix(lda, M, sizeof(A[0]), A, lda, d_A, lda);
                }

                if (opB==CUBLAS_OP_N) {
                    //B is ldb*N
                    cudaMalloc(reinterpret_cast<void**>(&d_B), ldb * N * sizeof(B[0]));
                    cublasSetMatrix(ldb, N, sizeof(B[0]), B, ldb, d_B, ldb);
                } else if (opB==CUBLAS_OP_T) {
                    //B is lda*K
                    cudaMalloc(reinterpret_cast<void**>(&d_B), ldb * K * sizeof(B[0]));
                    cublasSetMatrix(ldb, K, sizeof(B[0]), B, ldb, d_B, ldb);
                }

                cudaMalloc(reinterpret_cast<void**>(&d_C), ldc * N * sizeof(C[0]));
                cublasSetMatrix(ldc, N, sizeof(C[0]), C, ldc, d_C, ldc);

                /* Do GEMM */
                cublas_gemm(symPACK::handler, opA, opB,
                            M, N, K,
                            alpha, d_A, lda,
                            d_B, ldb, beta,
                            d_C, ldc);

                /* Copy matrices to host */
                cublasGetMatrix(ldc, N, sizeof(C[0]), d_C, ldc, C, ldc);

                /* Cleanup */
                cudaFree(d_A);
                cudaFree(d_B);
                cudaFree(d_C);
                break;
            
            case OP_SYRK:
                logfileptr->OFS()<<"DOING SYRK\n";
                if (opA==CUBLAS_OP_T) {
                    cudaMalloc(reinterpret_cast<void **>(&d_A), lda * N * sizeof(A[0]));
                    cublasSetMatrix(lda, N, sizeof(A[0]), A, lda, d_A, lda);
                } else {
                    cudaMalloc(reinterpret_cast<void **>(&d_A), lda * K * sizeof(A[0]));
                    cublasSetMatrix(lda, K, sizeof(A[0]), A, lda, d_A, lda);
                }

                cudaMalloc(reinterpret_cast<void **>(&d_C), ldc * N * sizeof(C[0]));

                cublas_syrk(symPACK::handler, fill, opA, 
                            N, K, 
                            alpha, d_A, lda,
                            beta, d_C, ldc);

                cublasGetMatrix(ldc, N, sizeof(C[0]), d_C, N, C, ldc);

                cudaFree(d_A);
                cudaFree(d_C);
                break;
        }
    };


}
}