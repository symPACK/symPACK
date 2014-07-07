#ifndef _NGCHOL_SP_BLAS_HEADER_
#define _NGCHOL_SP_BLAS_HEADER_

#include "ngchol/Environment.hpp"
#include "ngchol/blas.hpp"
#include "ngchol/DistSparseMatrix.hpp"
#include "ngchol/SparseMatrixStructure.hpp"

namespace LIBCHOLESKY{

//FROM SUPERLU
template <typename MYSCALAR>  int 
sp_dgemv_dist(char *trans, MYSCALAR alpha, DistSparseMatrix<MYSCALAR> & A, MYSCALAR *x, 
    int incx, MYSCALAR beta, MYSCALAR *y, int incy)
{

  /* Local variables */
  SparseMatrixStructure Astruct = A.GetLocalStructure();
  SparseMatrixStructure GAstruct = A.GetGlobalStructure();
  MYSCALAR *Aval = A.nzvalLocal.Data();
  int info;
  double temp;
  int lenx, leny, i, j, irow;
  int iy, jx, jy, kx, ky;
  int notran;

  notran = true;//lsame_(trans, "N");

  /* Test the input parameters */
  info = 0;
  //if ( !notran && !lsame_(trans, "T") && !lsame_(trans, "C")) info = 1;
  /*else*/ if ( A.size < 0 ) info = 3;
  else if (incx == 0) info = 5;
  else if (incy == 0)	info = 8;
  if (info != 0) {
    //    xerbla_("sp_dgemv_dist ", &info);
    return 0;
  }

  /* Quick return if possible. */
  if (A.size == 0 || alpha == 0. && beta == 1.)
    return 0;

  /* Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
     up the start points in  X  and  Y. */

  if (notran) {
    lenx = A.size;
    leny = lenx;
  } else {
    lenx = A.size;
    leny = lenx;
  }

  if (incx > 0) kx = 0;
  else kx =  - (lenx - 1) * incx;
  if (incy > 0) ky = 0;
  else ky =  - (leny - 1) * incy;

  /* Start the operations. In this version the elements of A are   
     accessed sequentially with one pass through A. */
  /* First form  y := beta*y. */
  MYSCALAR * ytmp = y;//new MYSCALAR[leny];

  if (beta != 1.) {
    if (incy == 1) {
      if (beta == 0.)
        for (i = 0; i < leny; ++i) ytmp[i] = 0.;
      else
        for (i = 0; i < leny; ++i) ytmp[i] = beta * y[i];
    } else {
      iy = ky;
      if (beta == 0.)
        for (i = 0; i < leny; ++i) {
          ytmp[iy] = 0.;
          iy += incy;
        }
      else
        for (i = 0; i < leny; ++i) {
          ytmp[iy] = beta * y[iy];
          iy += incy;
        }
    }
  }

  if (alpha == 0.) return 0;

  if ( notran ) {


    /* Form  y := alpha*A*x + y. */
    jx = kx;
    if (incy == 1) {
      int myfirstcol = 1 + iam*(A.size/np);
      int mylastcol = A.size/np + iam*(A.size/np);
      if(iam ==np-1){ mylastcol = A.size;}
#ifdef _DEBUG_
      logfileptr->OFS()<<"Owns col "<<myfirstcol<<" to "<<mylastcol<<std::endl;
#endif
      for (j = 1; j <= A.size; ++j) {
        if (x[jx] != 0.) {
          temp = alpha * x[jx];

          // only the lower triangular part of A is stored
          //this recreates the upper triangular part
          for (int jloc = 1; jloc < Astruct.colptr.m(); ++jloc) {
            int jcol = jloc + iam*(A.size/np);
            if(jcol<j){
              for (i = Astruct.colptr[jloc-1]; i < Astruct.colptr[jloc]; ++i) {
                irow = Astruct.rowind[i-1];
                if(irow == j ){
#ifdef _DEBUG_
                  logfileptr->OFS()<<"1 y["<<jcol<<"] += x["<<jx+1<<"] * A["<<jcol<<","<<irow<<"]"<<std::endl;
#endif
                  ytmp[jcol-1] += temp * Aval[i-1];
                }
                else if(irow>j){
                  break;
                }
              }
            } 
            else{
              break;
            }
          }

          if(j>=myfirstcol && j<=mylastcol){
            int jloc = j-myfirstcol;
            for (i = Astruct.colptr[jloc]; i < Astruct.colptr[jloc+1]; ++i) {
              irow = Astruct.rowind[i-1];
#ifdef _DEBUG_
              logfileptr->OFS()<<"2 y["<<irow<<"] += x["<<jx+1<<"] * A["<<irow<<","<<j<<"] "<<std::endl;
#endif
              ytmp[irow-1] += temp * Aval[i-1];
            }
          }
        }
        jx += incx;
      }

      if(iam==0){
        MYSCALAR * yrecv = new MYSCALAR[leny];
        for(i=1; i<np;++i){
          MPI_Recv(yrecv,leny*sizeof(MYSCALAR),MPI_BYTE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          blas::Axpy(leny,LIBCHOLESKY::ONE<MYSCALAR>(),yrecv,1,y,1);
        }
        delete [] yrecv;
      }
      else{
        MPI_Send(y,leny*sizeof(MYSCALAR),MPI_BYTE,0,0,MPI_COMM_WORLD);
      }

      MPI_Bcast(y,leny*sizeof(MYSCALAR),MPI_BYTE,0,MPI_COMM_WORLD);

      //we need to reduce the column now
      //      MPI_Allreduce(ytmp,y,leny,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 





    } else {
      //      ABORT("Not implemented.");
    }
  } 
  else {
    //    /* Form  y := alpha*A'*x + y. */
    //    jy = ky;
    //    if (incx == 1) {
    //      for (j = 0; j < A.size; ++j) {
    //        temp = 0.;
    //        for (i = Astruct.colptr[j]; i < Astruct.colptr[j+1]; ++i) {
    //          irow = Astruct.rowind[i-1];
    //          temp += Aval[i-1] * x[irow-1];
    //        }
    //        y[jy] += alpha * temp;
    //        jy += incy;
    //      }
    //    } else {
    ////      ABORT("Not implemented.");
    //    }
  }


  //     delete [] ytmp;


  return 0;
}



template <typename MYSCALAR>  int 
sp_dgemm_dist(char *transa, char *transb, int m, int n, int k, 
    MYSCALAR alpha, DistSparseMatrix<MYSCALAR> & A, MYSCALAR *b, int ldb, 
    MYSCALAR beta, MYSCALAR *c, int ldc)
{

  int    incx = 1, incy = 1;
  int    j;

  for (j = 0; j < n; ++j) {
    sp_dgemv_dist(transa, alpha, A, &b[ldb*j], incx, beta, &c[ldc*j], incy);


  }
  return 0;
}


} //end namespace LIBCHOLESKY

#endif // _NGCHOL_SP_BLAS_HEADER_
