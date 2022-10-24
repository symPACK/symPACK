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

    template <typename T>
    cublasStatus_t cublas_axpy_wrapper(Int N,
                           const T           DA,
                           const T           * DX, Int incx,
                           T                 * DY, Int incy) {
        
        int rank;
        int gpu_id;
        int n_gpus;

        MPI_Comm_rank(symPACK::world_comm, &rank);
        cudaGetDeviceCount(&n_gpus);

        gpu_id = rank % n_gpus;

        logfileptr->OFS()<<"DOING AXPY on GPU " << gpu_id << "\n";
        cudaSetDevice(gpu_id);

        T *d_X;
        T *d_Y;

        long dimx;
        dimx = (1 + (N-1)*abs(incx));
        CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_X), dimx * sizeof(DX[0])));
        CUBLAS_ERROR_CHECK(cublasSetVector(dimx, sizeof(DX[0]), DX, incx, d_X, incx));

        long dimy;
        dimy = (1 + (N-1)*abs(incy));
        CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_Y), dimy * sizeof(DY[0])));
        CUBLAS_ERROR_CHECK(cublasSetVector(dimy, sizeof(DY[0]), DY, incy, d_Y, incy));

        CUBLAS_ERROR_CHECK(cublas_axpy(symPACK::handlers[gpu_id], N, &DA, d_X, incx, d_Y, incy));

        CUBLAS_ERROR_CHECK(cublasGetVector(dimy, sizeof(DY[0]), d_Y, incy, DY, incy));

        CUDA_ERROR_CHECK(cudaFree(d_X));
        CUDA_ERROR_CHECK(cudaFree(d_Y));

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
    cublasStatus_t  cublas_scal_wrapper(Int N,
                            const T           DA,
                            T           * DX, Int incx) {
        int rank;
        int gpu_id;
        int n_gpus;

        MPI_Comm_rank(symPACK::world_comm, &rank);
        cudaGetDeviceCount(&n_gpus);

        gpu_id = rank % n_gpus;

        logfileptr->OFS()<<"DOING SCAL on GPU " << gpu_id << "\n";
        cudaSetDevice(gpu_id);

        T * d_X;

        long dimx;
        dimx = (1 + (N - 1) * abs(incx));

        CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_X), dimx * sizeof(DX[0])));

        CUBLAS_ERROR_CHECK(cublasSetVector(dimx, sizeof(DX[0]), DX, incx, d_X, incx));

        CUBLAS_ERROR_CHECK(cublas_scal(symPACK::handlers[gpu_id], N, &DA, DX, incx));

        CUBLAS_ERROR_CHECK(cublasGetVector(dimx, sizeof(DX[0]), d_X, incx, DX, incx));

        CUDA_ERROR_CHECK(cudaFree(d_X));

        return CUBLAS_STATUS_SUCCESS;

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
        
        int rank;
        int gpu_id;
        int n_gpus;

        MPI_Comm_rank(symPACK::world_comm, &rank);
        cudaGetDeviceCount(&n_gpus);

        gpu_id = rank % n_gpus;

        logfileptr->OFS()<<"DOING GEMV on GPU " << gpu_id << "\n";
        cudaSetDevice(gpu_id);
        
        /* Device buffers */
        T *d_A;
        T *d_X;
        T *d_Y;

        cublasStatus_t status;

        /* Set up device buffers */
        CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_A), lda * N * sizeof(A[0])));
        CUBLAS_ERROR_CHECK(cublasSetMatrix(lda, N, sizeof(A[0]), A, lda, d_A, lda));

        long dimx;
        dimx = (1 + (M - 1) * abs(incx));
        
    
        CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_X), dimx * sizeof(X[0])));
        CUBLAS_ERROR_CHECK(cublasSetVector(dimx, sizeof(X[0]), X, incx, d_X, incx));

        long dimy;
        dimy = (1 + (N - 1) * abs(incy));
        
        CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_Y), dimy * sizeof(Y[0])));
        CUBLAS_ERROR_CHECK(cublasSetVector(dimy, sizeof(Y[0]), Y, incy, d_Y, incy));

        /* do gemv */
        CUBLAS_ERROR_CHECK(cublas_gemv(symPACK::handlers[gpu_id], op,
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


    template <typename T>
    cublasStatus_t cublas_ger_wrapper(Int M, Int N,
                           const T           alpha,
                           const T           * X, Int incx,
                           const T           * Y, Int incy,
                           T           * A, Int lda) {
        int rank;
        int gpu_id;
        int n_gpus;

        MPI_Comm_rank(symPACK::world_comm, &rank);
        cudaGetDeviceCount(&n_gpus);

        gpu_id = rank % n_gpus;

        logfileptr->OFS()<<"DOING GER on GPU " << gpu_id << "\n";
        cudaSetDevice(gpu_id);

        T * d_X;
        T * d_Y;
        T * d_A;

        CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_A), lda * N * sizeof(A[0])));
        CUBLAS_ERROR_CHECK(cublasSetMatrix(lda, N, sizeof(A[0]), A, lda, d_A, lda));
        
        long dimx;
        dimx = (1 + (M - 1) * abs(incx));

        CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_X), dimx * sizeof(X[0])));
        CUBLAS_ERROR_CHECK(cublasSetVector(dimx, sizeof(X[0]), X, incx, d_X, incx));

        long dimy;
        dimy = (1 + (N - 1) * abs(incy));

        
        CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_Y), dimy * sizeof(Y[0])));
        CUBLAS_ERROR_CHECK(cublasSetVector(dimy, sizeof(Y[0]), Y, incy, d_Y, incy));

        if constexpr (std::is_same_v<T, const float> || std::is_same_v<T, const double>)
            CUBLAS_ERROR_CHECK(cublas_ger(symPACK::handlers[gpu_id], M, N, alpha, d_X, incx, d_Y, incy, d_A, lda));
        if constexpr (std::is_same_v<T, const cuComplex> || std::is_same_v<T, const cuDoubleComplex>)
            CUBLAS_ERROR_CHECK(cublas_geru(symPACK::handlers[gpu_id], M, N, alpha, d_X, incx, d_Y, incy, d_A, lda));

        CUBLAS_ERROR_CHECK(cublasGetMatrix(lda, N, sizeof(A[0]), d_A, lda, A, lda));

        CUDA_ERROR_CHECK(cudaFree(d_X));
        CUDA_ERROR_CHECK(cudaFree(d_Y));
        CUDA_ERROR_CHECK(cudaFree(d_A));

        return CUBLAS_STATUS_SUCCESS;
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
    cublasStatus_t cublas_gemm_wrapper(
                           cublasOperation_t opA, cublasOperation_t opB,
                           Int M, Int N, Int K,
                           const T           alpha,
                           const T           * A , Int lda,
                           const T          * B, Int ldb,
                           const T           beta,
                           T           * C, Int ldc) {
        int rank;
        int gpu_id;
        int n_gpus;

        MPI_Comm_rank(symPACK::world_comm, &rank);
        cudaGetDeviceCount(&n_gpus);

        gpu_id = rank % n_gpus;

        T * d_A;
        T * d_B;
        T * d_C;

        logfileptr->OFS()<<"DOING GEMM on GPU " << gpu_id << "\n";
        cudaSetDevice(gpu_id);
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
        CUBLAS_ERROR_CHECK(cublas_gemm(symPACK::handlers[gpu_id], opA, opB,
                    M, N, K,
                    &alpha, d_A, lda,
                    d_B, ldb, &beta,
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
    cublasStatus_t cublas_trsm_wrapper(cublasSideMode_t side, cublasFillMode_t fill, cublasOperation_t op, cublasDiagType_t diag,
                        Int M , Int N, 
                        const T alpha, T * A, Int lda, 
                        T * B, Int ldb) {
        int rank;
        int gpu_id;
        int n_gpus;

        MPI_Comm_rank(symPACK::world_comm, &rank);
        cudaGetDeviceCount(&n_gpus);

        gpu_id = rank % n_gpus;

        logfileptr->OFS()<<"DOING TRSM on GPU " << gpu_id << "\n";
        cudaSetDevice(gpu_id);

        T * d_A;
        T * d_B;
        
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
        CUBLAS_ERROR_CHECK(cublas_trsm(symPACK::handlers[gpu_id], side, fill, op, diag, 
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

    
            
           /* case OP_SYRK:
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
                break;*/
        


}
}