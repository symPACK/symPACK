/* cuBLAS interface
 * Author: Julian Bellavita, UC Berkeley
 */
#include  "sympack/Environment.hpp"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cublas_v2.h"

                                                                                               \

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


        CUBLAS_ERROR_CHECK(cublasGetVector(dimy, sizeof(DY[0]), d_Y, incy, DY, incy));

        CUDA_ERROR_CHECK(cudaFree(d_X));
        CUDA_ERROR_CHECK(cudaFree(d_Y));

        return CUBLAS_STATUS_SUCCESS;

    }

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
	return CUBLAS_STATUS_SUCCESS;
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


        CUBLAS_ERROR_CHECK(cublasGetVector(dimx, sizeof(DX[0]), d_X, incx, DX, incx));

        CUDA_ERROR_CHECK(cudaFree(d_X));

        return CUBLAS_STATUS_SUCCESS;

    }

    template <typename T>
    cublasStatus_t  cublas_scal_wrapper2(Int N,
                            const T           DA,
                            T           * DX, Int incx) {
        
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
        if constexpr (std::is_same_v<T, const cuComplex> || std::is_same_v<T, const cuDoubleComplex>)

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
    
    template <typename T>
    cublasStatus_t cublas_syrk_wrapper2(cublasFillMode_t uplo, cublasOperation_t trans,
                                        Int N, Int K,
                                        const T alpha,
                                        const T * A, Int lda,
                                        const T beta,
                                        T * C, Int ldc
                                        ) {
        
        //logfileptr->OFS()<<"DOING SYRK"<<std::endl;

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

        auto start = std::chrono::system_clock::now();

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
            //CUBLAS_ERROR_CHECK(cublasSetMatrix(lda, K, sizeof(A[0]), A, lda, d_A, lda));
            CUDA_ERROR_CHECK(cudaMemcpyAsync(d_A, A, lda * K * sizeof(A[0]), cudaMemcpyHostToDevice));
        } else if (opA==CUBLAS_OP_T) {
            //A is lda*M
            CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_A), lda * M * sizeof(A[0])));
            //CUBLAS_ERROR_CHECK(cublasSetMatrix(lda, M, sizeof(A[0]), A, lda, d_A, lda));
            CUDA_ERROR_CHECK(cudaMemcpyAsync(d_A, A, lda * M * sizeof(A[0]), cudaMemcpyHostToDevice));
        }

        if (opB==CUBLAS_OP_N) {
            //B is ldb*N
            CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_B), ldb * N * sizeof(B[0])));
            //CUBLAS_ERROR_CHECK(cublasSetMatrix(ldb, N, sizeof(B[0]), B, ldb, d_B, ldb));
            CUDA_ERROR_CHECK(cudaMemcpyAsync(d_B, B, ldb * N * sizeof(B[0]), cudaMemcpyHostToDevice));
        } else if (opB==CUBLAS_OP_T) {
            //B is lda*K
            CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_B), ldb * K * sizeof(B[0])));
            //CUBLAS_ERROR_CHECK(cublasSetMatrix(ldb, K, sizeof(B[0]), B, ldb, d_B, ldb));
            CUDA_ERROR_CHECK(cudaMemcpyAsync(d_B, B, ldb * K * sizeof(B[0]), cudaMemcpyHostToDevice));
        }
        CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_C), ldc * N * sizeof(C[0])));
        //CUBLAS_ERROR_CHECK(cublasSetMatrix(ldc, N, sizeof(C[0]), C, ldc, d_C, ldc));
        CUDA_ERROR_CHECK(cudaMemcpyAsync(d_C, C, ldc * N * sizeof(C[0]), cudaMemcpyHostToDevice));
        //CUDA_ERROR_CHECK(cudaDeviceSynchronize());

        /* Do GEMM */

        /* Copy matrices to host */
        //CUBLAS_ERROR_CHECK(cublasGetMatrix(ldc, N, sizeof(C[0]), d_C, ldc, C, ldc));
        CUDA_ERROR_CHECK(cudaMemcpyAsync(C, d_C, ldc * N * sizeof(C[0]), cudaMemcpyDeviceToHost));

        /* Cleanup */
        CUDA_ERROR_CHECK(cudaFree(d_A));
        CUDA_ERROR_CHECK(cudaFree(d_B));
        CUDA_ERROR_CHECK(cudaFree(d_C));

        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = end - start;
        statfileptr->OFS() << "GEMM time: " << diff.count() << std::endl;

        return CUBLAS_STATUS_SUCCESS;
    }

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
    cublasStatus_t cublas_trsm_wrapper(cublasSideMode_t side, cublasFillMode_t fill, cublasOperation_t op, cublasDiagType_t diag,
                        Int M , Int N, 
                        const T alpha, T * A, Int lda, 
                        T * B, Int ldb) {
        int rank;
        int gpu_id;
        int n_gpus;

        auto start = std::chrono::system_clock::now();

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
            CUDA_ERROR_CHECK(cudaMemcpyAsync(d_A, A, lda * M * sizeof(A[0]), cudaMemcpyHostToDevice));
        } else {
            CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_A), lda * N * sizeof(A[0])));
            CUDA_ERROR_CHECK(cudaMemcpyAsync(d_A, A, lda * N * sizeof(A[0]), cudaMemcpyHostToDevice));
        }
        CUDA_ERROR_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_B), ldb * N * sizeof(B[0])));
        CUDA_ERROR_CHECK(cudaMemcpyAsync(d_B, B, ldb * N * sizeof(B[0]), cudaMemcpyHostToDevice));
                
        /* Do TRSM */

        /* Copy matrices to host */
        CUDA_ERROR_CHECK(cudaMemcpyAsync(B, d_B, ldb * N * sizeof(B[0]), cudaMemcpyDeviceToHost));

        /* Cleanup */
        CUDA_ERROR_CHECK(cudaFree(d_A));
        CUDA_ERROR_CHECK(cudaFree(d_B));
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = end - start;
        statfileptr->OFS() << "TRSM time: " << diff.count() << std::endl;

        return CUBLAS_STATUS_SUCCESS;
    }

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

    
    /* UTILITY FUNCTIONS */

    /* Expand buffer pointed to by DEV_PTR from size n to m. Return pointer to new buffer. */
    template <typename T, typename U>
    upcxx::global_ptr<T, upcxx::memory_kind::cuda_device> cudaReallocManual(upcxx::global_ptr<T, upcxx::memory_kind::cuda_device>& dev_ptr, 
                                                                            size_t n, size_t m) {
        assert(m>n);
        logfileptr->OFS() << "Calling cudaRealloc\n";
        upcxx::global_ptr<T, upcxx::memory_kind::cuda_device> new_ptr = symPACK::gpu_allocator.allocate<T>(m*sizeof(U));
        upcxx::copy(dev_ptr, new_ptr, n*sizeof(U)).wait();
        symPACK::gpu_allocator.deallocate(dev_ptr);
        return new_ptr;
    }

    /* Dumps the contents of a global pointer to the device into the logfile */
    template <typename T>
    void print_dev_buffer(upcxx::global_ptr<T, upcxx::memory_kind::cuda_device> ptr, size_t n) {

        // Copy to host buffer
        T * h_local_ptr = new T[n];
        upcxx::copy(ptr, h_local_ptr, n).wait();

        logfileptr->OFS()<<"CONTENTS OF DEVICE BUFFER OF SIZE "<<n<<std::endl;
        for (int i=0; i<n; i++) {
            logfileptr->OFS()<<h_local_ptr[i];
            if (i%10==0 && i>0) logfileptr->OFS()<<std::endl; //for readability
        }

    }
        


}
}
