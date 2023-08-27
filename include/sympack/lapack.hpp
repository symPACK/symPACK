#ifndef _SYMPACK_LAPACK_HEADER_
#define _SYMPACK_LAPACK_HEADER_

#include  "sympack/Environment.hpp"
#include  "sympack/blas.hpp"
#ifdef CUDA_MODE
#include "cusolverDn.h"
#endif

namespace symPACK {

  /// @namespace lapack
  ///
  /// @brief Thin interface to LAPACK.
  namespace lapack {

    typedef  int                    Int; 
    typedef  std::complex<float>    scomplex;
    typedef  std::complex<double>   dcomplex;

    int Ilaenv( int Ispec, const char * name, const char * opts, int n1, int n2, int n3, int n4);

    // *********************************************************************
    // Cholesky factorization
    // *********************************************************************

    void Potrf( char uplo, Int n, const float* A, Int lda );
    void Potrf( char uplo, Int n, const double* A, Int lda );
    void Potrf( char uplo, Int n, const scomplex* A, Int lda );
    void Potrf( char uplo, Int n, const dcomplex* A, Int lda );
#ifdef CUDA_MODE
    cusolverStatus_t cusolver_potrf(cusolverDnHandle_t handle, char uplo, Int n, double* A, Int lda);
#endif

    void Potrs( char uplo, Int n, Int nrhs, const float* A, Int lda, float* B, Int ldb);
    void Potrs( char uplo, Int n, Int nrhs, const double* A, Int lda, double* B, Int ldb);
    void Potrs( char uplo, Int n, Int nrhs, const scomplex* A, Int lda, scomplex* B, Int ldb);
    void Potrs( char uplo, Int n, Int nrhs, const dcomplex* A, Int lda, dcomplex* B, Int ldb);

    void Posv( char uplo, Int n, Int nrhs, const float* A, Int lda, float* B, Int ldb);
    void Posv( char uplo, Int n, Int nrhs, const double* A, Int lda, double* B, Int ldb);
    void Posv( char uplo, Int n, Int nrhs, const scomplex* A, Int lda, scomplex* B, Int ldb);
    void Posv( char uplo, Int n, Int nrhs, const dcomplex* A, Int lda, dcomplex* B, Int ldb);

    void Potf2( char uplo, Int n, const float* A, Int lda );
    void Potf2( char uplo, Int n, const double* A, Int lda );
    void Potf2( char uplo, Int n, const scomplex* A, Int lda );
    void Potf2( char uplo, Int n, const dcomplex* A, Int lda );


    // *********************************************************************
    // LU factorization (with partial pivoting)
    // *********************************************************************

    void Getrf( Int m, Int n, float* A, Int lda, Int* p );
    void Getrf( Int m, Int n, double* A, Int lda, Int* p );
    void Getrf( Int m, Int n, scomplex* A, Int lda, Int* p );
    void Getrf( Int m, Int n, dcomplex* A, Int lda, Int* p );

    // *********************************************************************
    // For reducing well-conditioned Hermitian generalized-definite EVP's
    // to standard form.
    // *********************************************************************

    void Hegst
      ( Int itype, char uplo, 
        Int n, float* A, Int lda, const float* B, Int ldb );
    void Hegst
      ( Int itype, char uplo,
        Int n, double* A, Int lda, const double* B, Int ldb );
    void Hegst
      ( Int itype, char uplo,
        Int n, scomplex* A, Int lda, const scomplex* B, Int ldb );
    void Hegst
      ( Int itype, char uplo,
        Int n, dcomplex* A, Int lda, const dcomplex* B, Int ldb );

    // *********************************************************************
    // For computing the inverse of a triangular matrix
    // *********************************************************************

    void Trtri
      ( char uplo, char diag, Int n, const float* A, Int lda );
    void Trtri
      ( char uplo, char diag, Int n, const double* A, Int lda );
    void Trtri
      ( char uplo, char diag, Int n, const scomplex* A, Int lda );
    void Trtri
      ( char uplo, char diag, Int n, const dcomplex* A, Int lda );


    // *********************************************************************
    // Compute the SVD of a general matrix using a divide and conquer algorithm
    // *********************************************************************

    void DivideAndConquerSVD
      ( Int m, Int n, float* A, Int lda, 
        float* s, float* U, Int ldu, float* VTrans, Int ldvt );
    void DivideAndConquerSVD
      ( Int m, Int n, double* A, Int lda, 
        double* s, double* U, Int ldu, double* VTrans, Int ldvt );
    void DivideAndConquerSVD
      ( Int m, Int n, scomplex* A, Int lda, 
        float* s, scomplex* U, Int ldu, scomplex* VAdj, Int ldva );
    void DivideAndConquerSVD
      ( Int m, Int n, dcomplex* A, Int lda, 
        double* s, dcomplex* U, Int ldu, dcomplex* VAdj, Int ldva );

    //
    // Compute the SVD of a general matrix using the QR algorithm
    //

    void QRSVD
      ( Int m, Int n, float* A, Int lda, 
        float* s, float* U, Int ldu, float* VTrans, Int ldvt );
    void QRSVD
      ( Int m, Int n, double* A, Int lda, 
        double* s, double* U, Int ldu, double* VTrans, Int ldvt );
    void QRSVD
      ( Int m, Int n, scomplex* A, Int lda, 
        float* s, scomplex* U, Int ldu, scomplex* VAdj, Int ldva );
    void QRSVD
      ( Int m, Int n, dcomplex* A, Int lda, 
        double* s, dcomplex* U, Int ldu, dcomplex* VAdj, Int ldva );


    // *********************************************************************
    // Compute the singular values of a general matrix using the QR algorithm
    // *********************************************************************

    void SingularValues( Int m, Int n, float* A, Int lda, float* s );
    void SingularValues( Int m, Int n, double* A, Int lda, double* s );
    void SingularValues( Int m, Int n, scomplex* A, Int lda, float* s );
    void SingularValues( Int m, Int n, dcomplex* A, Int lda, double* s );

    // *********************************************************************
    // Compute the SVD of a bidiagonal matrix using the QR algorithm
    // *********************************************************************

    void BidiagQRAlg
      ( char uplo, Int n, Int numColsVTrans, Int numRowsU,
        float* d, float* e, float* VTrans, Int ldVTrans, float* U, Int ldU );
    void BidiagQRAlg
      ( char uplo, Int n, Int numColsVTrans, Int numRowsU, 
        double* d, double* e, double* VTrans, Int ldVTrans, double* U, Int ldU );
    void BidiagQRAlg
      ( char uplo, Int n, Int numColsVAdj, Int numRowsU,
        float* d, float* e, scomplex* VAdj, Int ldVAdj, scomplex* U, Int ldU );
    void BidiagQRAlg
      ( char uplo, Int n, Int numColsVAdj, Int numRowsU, 
        double* d, double* e, dcomplex* VAdj, Int ldVAdj, dcomplex* U, Int ldU );

    // *********************************************************************
    // Compute the linear least square problem using SVD
    // *********************************************************************
    void SVDLeastSquare( Int m, Int n, Int nrhs, float * A, Int lda,
        float * B, Int ldb, float * S, float rcond,
        Int* rank );
    void SVDLeastSquare( Int m, Int n, Int nrhs, double * A, Int lda,
        double * B, Int ldb, double * S, double rcond,
        Int* rank );
    void SVDLeastSquare( Int m, Int n, Int nrhs, scomplex * A, Int lda,
        scomplex * B, Int ldb, float * S, float rcond,
        Int* rank );
    void SVDLeastSquare( Int m, Int n, Int nrhs, dcomplex * A, Int lda,
        dcomplex * B, Int ldb, double * S, double rcond,
        Int* rank );


    // *********************************************************************
    // Copy
    // *********************************************************************
    void Lacpy( char uplo, Int m, Int n, const float* A, Int lda,
        float* B, Int ldb	);

    void Lacpy( char uplo, Int m, Int n, const double* A, Int lda,
        double* B, Int ldb	);

    void Lacpy( char uplo, Int m, Int n, const scomplex* A, Int lda,
        scomplex* B, Int ldb	);

    void Lacpy( char uplo, Int m, Int n, const dcomplex* A, Int lda,
        dcomplex* B, Int ldb	);

    // *********************************************************************
    // Inverting a factorized matrix: Getri
    // *********************************************************************


    void Getri ( Int n, double* A, Int lda, const Int* ipiv );

    void Getri ( Int n, dcomplex* A, Int lda, const Int* ipiv );

    double Lange(char norm, Int m, Int n, const double* A, Int lda);
    double Lange(char norm, Int m, Int n, const dcomplex* A, Int lda);


    // SCAL
    void Scal( Int N, float DA,    float* DX,    Int INCX);
    void Scal( Int N, double DA,   double* DX,   Int INCX);
    void Scal( Int N, scomplex DA, scomplex* DX, Int INCX);
    void Scal( Int N, dcomplex DA, dcomplex* DX, Int INCX);


    template<typename T>
      void Potf2_LDL( const char * UPLO, Idx N, T * A,  Idx LDA, T* WORK, Int & INFO ){
        //*
        //*  -- LAPACK-like routine --
        //*     Esmond G. Ng, Oak Ridge National Laboratory
        //*     March 17, 1998
        //*
        //*     .. Scalar Arguments ..
        //      CHARACTER          UPLO
        //      INTEGER            INFO, LDA, N, NDEF
        //      DOUBLE PRECISION   TOL
        //*     ..
        //*     .. Array Arguments ..
        //      INTEGER            IDEF(*)
        //      DOUBLE PRECISION   A(LDA,*), WORK(*)
        //*     ..
        //*
        //*  Purpose
        //*  =======
        //*
        //*  DPOTF2_LDL computes the Cholesky factorization of a real symmetric
        //*  positive definite or semi-definite matrix A.
        //*
        //*  The factorization has the form
        //*     A = U' * D * U ,  if UPLO = 'U', or
        //*     A = L  * D * L',  if UPLO = 'L',
        //*  where U is a unit upper triangular matrix, L is unit lower
        //*  triangular, and D is diagonal.
        //*
        //*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
        //*
        //*  Arguments
        //*  =========
        //*
        //*  UPLO    (input) CHARACTER*1
        //*          Specifies whether the upper or lower triangular part of the
        //*          symmetric matrix A is stored.
        //*          = 'U':  Upper triangular
        //*          = 'L':  Lower triangular
        //*
        //*  N       (input) INTEGER
        //*          The order of the matrix A.  N >= 0.
        //*
        //*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
        //*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
        //*          n by n upper triangular part of A contains the upper
        //*          triangular part of the matrix A, and the strictly lower
        //*          triangular part of A is not referenced.  If UPLO = 'L', the
        //*          leading n by n lower triangular part of A contains the lower
        //*          triangular part of the matrix A, and the strictly upper
        //*          triangular part of A is not referenced.
        //*
        //*          On exit, if INFO = 0, the factor U+D or L+D from the Cholesky
        //*          factorization A = U'*D*U  or A = L*D*L'.
        //*
        //*  LDA     (input) INTEGER
        //*          The leading dimension of the array A.  LDA >= max(1,N).
        //*
        //*  NDEF    (output) INTEGER
        //*          The rank deficiency of the matrix A.  NDEF <= N.
        //*
        //*  IDEF    (output) INTEGER array, dimension NDEF
        //*          Indices of columns for which zero pivots are encountered.
        //*
        //*  TOL     (input) DOUBLE PRECISION
        //*          Tolerance for checking if a pivot is zero.
        //*
        //*  WORK    (input) DOUBLE PRECISION array, dimension N.
        //*          Temporary work array.
        //*
        //*  INFO    (output) INTEGER
        //*          = 0: successful exit
        //*          < 0: if INFO = -k, the k-th argument had an illegal value
        //*          > 0: if INFO = k, the leading minor of order k is not
        //*               positive definite or semi-definite, and the
        //*               factorization could not be completed.
        //*
        //*  =====================================================================
        //*
        //*       .. Parameters ..
        //        DOUBLE PRECISION   ONE, ZERO
        //        PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
        //*       ..
        //*       .. Local Scalars ..
        //        LOGICAL            UPPER
        //        INTEGER            I, J
        //        DOUBLE PRECISION   AJJ
        //*       ..
        //*       .. External Functions ..
        //        LOGICAL            LSAME
        //        DOUBLE PRECISION   DDOT
        //        EXTERNAL           LSAME, DDOT
        //*       ..
        //*       .. External Subroutines ..
        //        EXTERNAL           DGEMV, DSCAL, XERBLA
        //*       ..
        //*       .. Intrinsic Functions ..
        //        INTRINSIC          ABS, MAX
        //*       ..
        //*       .. Executable Statements ..
        //*
        //*       Test the input parameters.
        //*
        T AJJ;
        T ONE = T(1.0);
        T MINUS_ONE = T(-1.0);
        T ZERO = T(0);


        INFO = 0;
        bool UPPER = *UPLO == 'U';
        if ( !UPPER  && *UPLO!= 'L' ) {
          INFO = -1;
        }
        else if ( N > 0x7FFFFFFF ) { // means N would be overflow a signed 32-bit int
          INFO = -2;
        }
        else if ( LDA < std::max(Idx(1),N) ) {
          INFO = -4;
        }
        if( INFO != 0 ) {
          //CALL  XERBLA ( 'DPOTF2_LDL', -INFO )
          return;
        }
        //*
        //*       Quick return if possible
        //*
        if ( N == 0 ) return;
        if ( UPPER ) {
          //*
          //*           Compute the Cholesky factorization A = U'*D*U.
          //*
          for (Idx J = 1; J<=N; J++) {
            //*
            //*               Compute U(J,J) and test for non-positive-definiteness.
            //*
            blas::Copy( J-1, &A[(J-1)*LDA], 1, WORK, 1 );
            for (Idx I = 1; I<=J-1; I++) {
              WORK[I-1] = (WORK[I-1]) * A[I-1+(I-1)*LDA];
            }
            AJJ = A[J-1+(J-1)*LDA] - blas::Dotu( J-1, &A[(J-1)*LDA], 1, WORK, 1 );


            if ( std::abs(AJJ) > 0 ) {
              A[J-1+(J-1)*LDA] = AJJ;
              //*
              //*                   Compute elements J+1:N of row J.
              //*
              if( J < N ){
                blas::Gemv( 'T', J-1, N-J, MINUS_ONE, &A[J*LDA], LDA, WORK, 1,ONE, &A[J-1+J*LDA], LDA );
                blas::Scal( N-J, ONE/AJJ, &A[J-1+J*LDA], LDA );
              }
            }
            else {
              A[(J-1)+(J-1)*LDA] = AJJ;
              INFO = J;
              return;
            }
          }
        }
        else{
          //*
          //*           Compute the Cholesky factorization A = L*D*L'.
          //*
          for (Idx J=1; J<=N;J++){
            //*
            //*               Compute L(J,J) and test for non-positive-definiteness.
            //*
            blas::Copy( J-1, &A[J-1], LDA, WORK, 1 );
            for (Idx I = 1; I<=J-1; I++) {
              WORK[I-1] = (WORK[I-1]) * A[I-1+(I-1)*LDA];
            }
            AJJ = A[J-1+(J-1)*LDA] - blas::Dotu( J-1, &A[J-1], LDA, WORK, 1 );

            if ( std::abs(AJJ) > 0 ) {
              A[(J-1)+(J-1)*LDA] = AJJ;
              //*
              //*                   Compute elements J+1:N of column J.
              //*
              if( J < N ){
                blas::Gemv( 'N', N-J, J-1,MINUS_ONE, &A[J], LDA, WORK, 1,ONE, &A[J+(J-1)*LDA], 1 );
                blas::Scal( N-J, ONE/AJJ, &A[J+(J-1)*LDA], 1 );
              }
            }
            else {
              A[(J-1)+(J-1)*LDA] = AJJ;
              INFO = J;
              return;
            }
          }
        }
      }

    template<typename T>
      void Potrf_LDL( const char * UPLO,  Idx N, T * A,  Idx LDA, TempUpdateBuffers<T> & tmpBuffers, Int & INFO) {
        //*
        //*  -- LAPACK-like routine --
        //*     Esmond G. Ng, Oak Ridge National Laboratory
        //*
        //*     .. Scalar Arguments ..
        //      CHARACTER          UPLO
        //      INTEGER            INFO, LDA, N, NDEF
        //      DOUBLE PRECISION   TOL
        //*     ..
        //*     .. Array Arguments ..
        //      INTEGER            IDEF(*)
        //      DOUBLE PRECISION   A(LDA,*), WORK(*)
        //*     ..
        //*
        //*  Purpose
        //*  =======
        //*
        //*  DPOTRF_LDL computes the Cholesky factorization of a real symmetric
        //*  positive definite matrix A.
        //*
        //*  The factorization has the form
        //*     A = U**T * U,  if UPLO = 'U', or
        //*     A = L  * L**T,  if UPLO = 'L',
        //*  where U is an upper triangular matrix and L is lower triangular.
        //*
        //*  This is the block version of the algorithm, calling Level 3 BLAS.
        //*
        //*  Arguments
        //*  =========
        //*
        //*  UPLO    (input) CHARACTER*1
        //*          = 'U':  Upper triangle of A is stored;
        //*          = 'L':  Lower triangle of A is stored.
        //*
        //*  N       (input) INTEGER
        //*          The order of the matrix A.  N >= 0.
        //*
        //*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
        //*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
        //*          N-by-N upper triangular part of A contains the upper
        //*          triangular part of the matrix A, and the strictly lower
        //*          triangular part of A is not referenced.  If UPLO = 'L', the
        //*          leading N-by-N lower triangular part of A contains the lower
        //*          triangular part of the matrix A, and the strictly upper
        //*          triangular part of A is not referenced.
        //*
        //*          On exit, if INFO = 0, the factor U or L from the Cholesky
        //*          factorization A = U**T*U or A = L*L**T.
        //*
        //*  LDA     (input) INTEGER
        //*          The leading dimension of the array A.  LDA >= max(1,N).
        //*
        //*  NDEF    (output) INTEGER
        //*          The rank deficiency of the matrix A.  NDEF <= N.
        //*
        //*  IDEF    (output) INTEGER array, dimension NDEF
        //*          Indices of columns for which zero pivots are encountered.
        //*
        //*  TOL     (input) DOUBLE PRECISION
        //*          Tolerance for checking if a pivot is zero.
        //*
        //*  WORK    (input) DOUBLE PRECISION array, dimension N.
        //*          Temporary work array.
        //*
        //*  INFO    (output) INTEGER
        //*          = 0:  successful exit
        //*          < 0:  if INFO = -i, the i-th argument had an illegal value
        //*          > 0:  if INFO = i, the leading minor of order i is not
        //*                positive definite, and the factorization could not be
        //*                completed.
        //*
        //*  =====================================================================
        //*
        //*     .. Parameters ..
        //      DOUBLE PRECISION  ONE, ZERO
        //      PARAMETER         ( ONE = 1.0D+0, ZERO = 0.0D+0 )
        //*     ..
        //*     .. Local Scalars ..
        //      LOGICAL           UPPER
        //      INTEGER           COL, I, II, J, JB, NB, NDEF0, PTR
        //*     ..
        //*     .. External Functions ..
        //      LOGICAL           LSAME
        //      INTEGER           ILAENV
        //      EXTERNAL          LSAME, ILAENV
        //*     ..
        //*     .. External Subroutines ..
        //      EXTERNAL          DGEMM, DPOTF2_LDL, DSYRK, DTRSM, XERBLA
        //*     ..
        //*     .. Intrinsic Functions ..
        //      INTRINSIC         MAX, MIN
        //*     ..
        //*     .. Executable Statements ..
        //*
        //*     Test the input parameters.
        //*





        T ONE = T(1.0);
        T MINUS_ONE = T(-1.0);
        T ZERO = T(0.0);
        Idx JB,NB,ROW,COL;
        Ptr PTR;


        INFO = 0;
        bool UPPER = *UPLO == 'U';
        if ( !UPPER  && *UPLO!= 'L' ) {
          INFO = -1;
        }
        else if ( N > 0x7FFFFFFF ) { // means N would be overflow a signed 32-bit int
          INFO = -2;
        }
        else if ( LDA < std::max(Idx(1),N) ) {
          INFO = -4;
        }
        if( INFO != 0 ) {
          //CALL XERBLA( 'DPOTRF_LDL', -INFO )
          return;
        }
        //*
        //*       Quick return if possible
        //*
        if ( N == 0 ) return;

        //*
        //*     Determine the block size for this environment.
        //*
        NB = 64;//lapack::Ilaenv( 1, "DPOTRF", UPLO, N, -1, -1, -1 );

        if ( NB <= 1 || NB >= N ) {
          tmpBuffers.Resize(N,1);
          T * WORK = &tmpBuffers.tmpBuf[0];
          //*
          //*        Use unblocked code.
          //*
          Potf2_LDL( UPLO, N, A, LDA, WORK, INFO );
        }
        else {
          tmpBuffers.Resize(NB*N,1);
          T * WORK = &tmpBuffers.tmpBuf[0];
          //*
          //*        Use blocked code.
          //*
          if ( UPPER ) {
            //*
            //*           Compute the Cholesky factorization A = U'*U.
            //*
            for ( Idx J = 1; J <= N; J+=NB ) {
              //*
              //*              Update and factorize the current diagonal block and test
              //*              for non-positive-definiteness.
              //*
              JB = std::min( NB, N-J+1 );
              PTR = 1;
              for ( Idx I = 1; I<=J-1;I++) {
                blas::Copy( JB, &A[I-1+(J-1)*LDA], LDA, &WORK[PTR-1], 1 );
                blas::Scal( JB, A[I-1+(I-1)*LDA], &WORK[PTR-1], 1 );
                PTR += JB;
              }
              //WORK is already stored in a transposed fashion
              blas::Gemm( 'N', 'N', JB, JB, J-1, MINUS_ONE, WORK, JB, &A[(J-1)*LDA], LDA, ONE, &A[J-1+(J-1)*LDA], LDA );
              lapack::Potf2_LDL( "Upper", JB, &A[J-1+(J-1)*LDA], LDA, &WORK[PTR-1], INFO );
              if  ( INFO != 0 ) {
                INFO = INFO + J - 1;
                return;
              }
              if ( J+JB <= N ) {
                //*
                //*                 Compute the current block row.
                //*
                blas::Gemm( 'N', 'N', JB, N-J-JB+1,J-1, MINUS_ONE, WORK, JB, &A[ (J+JB-1)*LDA ],LDA, ONE, &A[ J-1 + (J+JB-1)*LDA ], LDA );
                blas::Trsm( 'L', 'U', 'T', 'U',JB, N-J-JB+1, ONE, &A[J-1+(J-1)*LDA], LDA, &A[J-1+(J+JB-1)*LDA], LDA );
                for ( Idx I = J; I<=J+JB-1;I++) {
                  blas::Scal( N-J-JB+1, ONE/A[I-1+(I-1)*LDA], &A[I-1+(J+JB-1)*LDA], LDA );
                }
              }
            }
          }
          else {
            //*
            //*               Compute the Cholesky factorization A = L*L'.
            //*
            for ( Idx J = 1; J <= N; J+=NB ) {
              //*
              //*                   Update and factorize the current diagonal block
              //*                   and test for non-positive-definiteness.
              //*
              JB = std::min( NB, N-J+1 );
              PTR = 1;
              for ( Idx I = 1; I<=J-1;I++) {
                blas::Copy( JB, &A[J-1+(I-1)*LDA],1,&WORK[PTR-1],1);
                blas::Scal( JB, A[I-1+(I-1)*LDA], &WORK[PTR-1], 1 );
                PTR += JB;
              }
              blas::Gemm( 'N','T', JB, JB, J-1, MINUS_ONE, &A[J-1], LDA, WORK, JB, ONE, &A[J-1+(J-1)*LDA], LDA );

              lapack::Potf2_LDL( "Lower", JB, &A[J-1+(J-1)*LDA], LDA, &WORK[PTR-1], INFO );
              if  ( INFO != 0 ) {
                INFO = INFO + J - 1;
                return;
              }
              if ( J+JB <= N ) {
                //*
                //*                       Compute the current block column.
                //*
                blas::Gemm( 'N','T',N-J-JB+1, JB, J-1,MINUS_ONE, &A[J+JB-1], LDA, WORK, JB, ONE, &A[J+JB-1+(J-1)*LDA], LDA );
                blas::Trsm( 'R','L','T','U', N-J-JB+1, JB, ONE, &A[J-1+(J-1)*LDA], LDA, &A[J+JB-1+(J-1)*LDA], LDA );
                for ( Idx I = J; I<=J+JB-1;I++) {
                  blas::Scal( N-J-JB+1, ONE/A[I-1+(I-1)*LDA], &A[J+JB-1+(I-1)*LDA], 1 );
                }
              }
            }
          }
        }
        return;
        //*
        //*     End of DPOTRF_LDL
        //*
      }


  } // namespace lapack
} // namespace SYMPACK
#endif //_SYMPACK_LAPACK_HEADER_
