/* Wrapper functions for each cuBLAS call
 * Author: Julian Bellavita, UC Berkeley
 */
#include  "sympack/Environment.hpp"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cublas_v2.h"

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


    /* TRSM */
    cublasStatus_t cublas_trsm(cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t,
                        Int, Int, const float *, const float *, Int, float *, Int);
    
    cublasStatus_t cublas_trsm(cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t,
                        Int, Int, const double *, const double *, Int, double *, Int);
    
    cublasStatus_t cublas_trsm(cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t,
                        Int, Int, const cuComplex *, const cuComplex *, Int, cuComplex *, Int);
    
    cublasStatus_t cublas_trsm(cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t,
                        Int, Int, const cuDoubleComplex *, const cuDoubleComplex *, Int, cuDoubleComplex *, Int);
    
    //template<typename T>
    //cublasStatus_t cublas_trsm(cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t,
    //                    Int, Int, const T *, const T *, Int, T *, Int);

    /* DO_CUBLAS_L1 */

    enum l1_cublas_ops{OP_AXPY,
                       OP_COPY,
                       OP_SCAL};
    

    void do_cublas_l1();
    

    /* DO_CUBLAS_L2 */

    enum l2_cublas_ops{OP_GEMV,
                       OP_GER,
                       OP_GERU,
                       OP_GERC};
    

    void do_cublas_l2();


    /* DO_CUBLAS_L3 */

    enum l3_cublas_ops{OP_SYRK,
                       OP_GEMM,
                       OP_TRSM};
    template <class T>
    void do_cublas_l3(l3_cublas_ops op, cublasSideMode_t side, cublasFillMode_t fill, 
                                  cublasOperation_t opA, cublasOperation_t opB,
                                  cublasDiagType_t diag, 
                                  Int M, Int N, Int K,
                                  T *alpha, T *A, Int lda, 
                                  T *beta, T *B, Int ldb, 
                                  T *C, Int ldc) {
        cublasHandle_t handler;
        cublasStatus_t stat;
        stat = cublasCreate(&handler);
        
        switch(op) {
            case OP_TRSM:

                /* Setup device matrices */
                T * d_A;
                T * d_B;
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
                cublas_trsm(handler, side, fill, opA, diag, M, N, alpha, d_A, lda, d_B, ldb);

                /* Copy matrices to host */
                cublasGetMatrix(ldb, N, sizeof(B[0]), d_B, ldb, B, ldb);

                /* Cleanup */
                cudaFree(d_A);
                cudaFree(d_B);
        }
    };


}
}