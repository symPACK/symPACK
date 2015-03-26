#ifndef IOHB_H
#define IOHB_H
/**
 * @file iohb.h
 * @author K.A. Remington, edited by Mark Hoemmen
 * @version 1.1
 * @date 18 Jul 2006
 *
 * Harwell-Boeing sparse matrix file I/O routines.  
 * Mark Hoemmen (mfh) made some modest bug fixes to iohb.c.
 */
#include <stdio.h>
#include <stdlib.h>
/*
 * MFH 20 Oct 2005, 18 Jul 2006: On many systems, malloc.h doesn't exist (e.g. MacOS X), 
 * and generally it's not needed anyway, as #including stdlib.h suffices.
 */
/*#include <malloc.h>*/ 
#ifdef __cplusplus
extern "C" {
#endif

int readHB_info(const char* filename, int* M, int* N, int* nz, char** Type, 
                                                      int* Nrhs);

int readHB_header(FILE* in_file, char* Title, char* Key, char* Type, 
                    int* Nrow, int* Ncol, int* Nnzero, int* Nrhs,
                    char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt, 
                    int* Ptrcrd, int* Indcrd, int* Valcrd, int* Rhscrd, 
                    char *Rhstype);

int readHB_mat_double(const char* filename, int colptr[], int rowind[], 
                                                                 double val[]);

int readHB_newmat_double(const char* filename, int* M, int* N, int* nonzeros, 
                         int** colptr, int** rowind, double** val);

int readHB_aux_double(const char* filename, const char AuxType, double b[]);

int readHB_newaux_double(const char* filename, const char AuxType, double** b);

int writeHB_mat_double(const char* filename, int M, int N, 
                        int nz, const int colptr[], const int rowind[], 
                        const double val[], int Nrhs, const double rhs[], 
                        const double guess[], const double exact[],
                        const char* Title, const char* Key, const char* Type, 
                        char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                        const char* Rhstype);

int readHB_mat_char(const char* filename, int colptr[], int rowind[], 
                                           char val[], char* Valfmt);

int readHB_newmat_char(const char* filename, int* M, int* N, int* nonzeros, int** colptr, 
                          int** rowind, char** val, char** Valfmt);

int readHB_aux_char(const char* filename, const char AuxType, char b[]);

int readHB_newaux_char(const char* filename, const char AuxType, char** b, char** Rhsfmt);

int writeHB_mat_char(const char* filename, int M, int N, 
                        int nz, const int colptr[], const int rowind[], 
                        const char val[], int Nrhs, const char rhs[], 
                        const char guess[], const char exact[], 
                        const char* Title, const char* Key, const char* Type, 
                        char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                        const char* Rhstype);

int ParseIfmt(char* fmt, int* perline, int* width);

int ParseRfmt(char* fmt, int* perline, int* width, int* prec, int* flag);

void IOHBTerminate(char* message);
#ifdef __cplusplus
}
#endif

#endif /* IOHB_H */
