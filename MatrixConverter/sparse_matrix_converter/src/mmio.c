/* 
Fri Aug 15 16:29:47 EDT 1997
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Matrix Market I/O library for ANSI C
    Roldan Pozo, NIST (pozo@nist.gov)
 
    See http://math.nist.gov/MatrixMarket for details and sample
    calling programs.

    mfh 1 Mar 2006:  Added extern declaration of strdup in 
    mm_typecode_to_str for systems (e.g. Itanium 2, gcc 3.4.3) 
    for which strdup is not declared in string.h.
 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                NOTICE

 Permission to use, copy, modify, and distribute this software and
 its documentation for any purpose and without fee is hereby granted
 provided that the above copyright notice appear in all copies and
 that both the copyright notice and this permission notice appear in
 supporting documentation.

 Neither the Author nor the Institution (National Institute of Standards
 and Technology) make any representations about the suitability of this 
 software for any purpose. This software is provided "as is" without 
 expressed or implied warranty.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
#include <bebop/util/config.h>
#include <bebop/smc/mmio.h>
#include <bebop/util/util.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h> /* malloc */
#include <ctype.h>




int mm_is_valid(char* matcode)
{
    if (!mm_is_matrix(matcode)) return 0;
    if (mm_is_dense(matcode) && mm_is_pattern(matcode)) return 0;
    if (mm_is_real(matcode) && mm_is_hermitian(matcode)) return 0;
    if (mm_is_pattern(matcode) && (mm_is_hermitian(matcode) || 
                mm_is_skew(matcode))) return 0;
    return 1;
}

int mm_read_banner(FILE *f, char** matcode)
{
    char line[MM_MAX_LINE_LENGTH];
    char banner[MM_MAX_TOKEN_LENGTH];
    char mtx[MM_MAX_TOKEN_LENGTH]; 
    char crd[MM_MAX_TOKEN_LENGTH];
    char data_type[MM_MAX_TOKEN_LENGTH];
    char storage_scheme[MM_MAX_TOKEN_LENGTH];
    char *p;

    /* WITH_DEBUG2(fprintf(stderr,"=== mm_read_banner ===\n")); */

    mm_clear_typecode(matcode);  

    if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL) 
        return MM_PREMATURE_EOF;

    /* WITH_DEBUG2(fprintf(stderr,"\tGot line %s\n", line)); */

    if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type, 
        storage_scheme) != 5)
        return MM_PREMATURE_EOF;

    /* WITH_DEBUG2(fprintf(stderr, "\tSplit line:  %s %s %s %s %s\n", banner, mtx, crd, 
       data_type, storage_scheme)); */

    /* 
     * MFH 11 June 2004: BUG: converting strings to upper case, but the strings
     * in the header file are lower case.
     */
    /* for (p=mtx; *p!='\0'; *p=toupper(*p),p++); */ /* convert to upper case */
    /* 
    for (p=crd; *p!='\0'; *p=toupper(*p),p++);  
    for (p=data_type; *p!='\0'; *p=toupper(*p),p++);
    for (p=storage_scheme; *p!='\0'; *p=toupper(*p),p++);
    */

    for (p=mtx; *p!='\0'; *p=tolower(*p),p++); 
    for (p=crd; *p!='\0'; *p=tolower(*p),p++);  
    for (p=data_type; *p!='\0'; *p=tolower(*p),p++);
    for (p=storage_scheme; *p!='\0'; *p=tolower(*p),p++);
 

    /* check for banner */
    if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
        return MM_NO_HEADER;

    /* WITH_DEBUG2(fprintf(stderr,"\tGot header\n")); */

    /* first field should be "mtx" */
    if (strcmp(mtx, MM_MTX_STR) != 0)
        return  MM_NOT_MTX;

    /* WITH_DEBUG2(fprintf(stderr,"\tGot \"matrix\"\n")); */

    mm_set_matrix(matcode);


    /* second field describes whether this is a sparse matrix (in coordinate
       storage) or a dense array */


    if (strcmp(crd, MM_SPARSE_STR) == 0)
      {
        mm_set_sparse(matcode);
        /* WITH_DEBUG2(fprintf(stderr,"\tGot sparse matrix marker\n")); */
      }
    else
    if (strcmp(crd, MM_DENSE_STR) == 0)
      {
        /* WITH_DEBUG2(fprintf(stderr,"\tGot dense matrix marker\n")); */
            mm_set_dense(matcode);
      }
    else
        return MM_UNSUPPORTED_TYPE;
    

    /* third field */

    if (strcmp(data_type, MM_REAL_STR) == 0)
      {
        /* WITH_DEBUG2(fprintf(stderr,"\tGot real matrix marker\n")); */
        mm_set_real(matcode);
      }
    else
    if (strcmp(data_type, MM_COMPLEX_STR) == 0)
      {
        /* WITH_DEBUG2(fprintf(stderr,"\tGot complex matrix marker\n")); */
        mm_set_complex(matcode);
      }
    else
    if (strcmp(data_type, MM_PATTERN_STR) == 0)
      {
        /* WITH_DEBUG2(fprintf(stderr,"\tGot pattern matrix marker\n")); */
        mm_set_pattern(matcode);
      }
    else
    if (strcmp(data_type, MM_INT_STR) == 0)
      {
        /* WITH_DEBUG2(fprintf(stderr,"\tGot integer matrix marker\n")); */
        mm_set_integer(matcode);
      }
    else
        return MM_UNSUPPORTED_TYPE;
    

    /* fourth field */

    if (strcmp(storage_scheme, MM_GENERAL_STR) == 0)
      {
        /* WITH_DEBUG2(fprintf(stderr,"\tGot general storage scheme marker\n")); */
        mm_set_general(matcode);
      }
    else
    if (strcmp(storage_scheme, MM_SYMM_STR) == 0)
      {
        /* WITH_DEBUG2(fprintf(stderr,"\tGot symmetric storage scheme marker\n")); */
        mm_set_symmetric(matcode);
      }
    else
    if (strcmp(storage_scheme, MM_HERM_STR) == 0)
      {
        /* WITH_DEBUG2(fprintf(stderr,"\tGot hermitian storage scheme marker\n")); */
        mm_set_hermitian(matcode);
      }
    else
    if (strcmp(storage_scheme, MM_SKEW_STR) == 0)
      {
        /* WITH_DEBUG2(fprintf(stderr,"\tGot skew storage scheme marker\n")); */
        mm_set_skew(matcode);
      }
    else
        return MM_UNSUPPORTED_TYPE;
        

    return 0;
}

int mm_write_mtx_crd_size(FILE *f, int M, int N, int nz)
{
    if (fprintf(f, "%d %d %d\n", M, N, nz) != 3)
        return MM_COULD_NOT_WRITE_FILE;
    else 
        return 0;
}

int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz )
{
    char line[MM_MAX_LINE_LENGTH];
    int num_items_read;

    /* set return null parameter values, in case we exit with errors */
    *M = *N = *nz = 0;

    /* now continue scanning until you reach the end-of-comments */
    do 
    {
        if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL) 
            return MM_PREMATURE_EOF;
    }while (line[0] == '%');

    /* line[] is either blank or has M,N, nz */
    if (sscanf(line, "%d %d %d", M, N, nz) == 3)
        return 0;
        
    else
    do
    { 
        num_items_read = fscanf(f, "%d %d %d", M, N, nz); 
        if (num_items_read == EOF) return MM_PREMATURE_EOF;
    }
    while (num_items_read != 3);

    return 0;
}


int mm_read_mtx_array_size(FILE *f, int *M, int *N)
{
    char line[MM_MAX_LINE_LENGTH];
    int num_items_read;

    /* set return null parameter values, in case we exit with errors */
    *M = *N = 0;

    /* now continue scanning until you reach the end-of-comments */
    do 
    {
        if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL) 
            return MM_PREMATURE_EOF;
    }while (line[0] == '%');

    /* line[] is either blank or has M,N, nz */
    if (sscanf(line, "%d %d", M, N) == 2)
        return 0;
        
    else /* we have a blank line */
    do
    { 
        num_items_read = fscanf(f, "%d %d", M, N); 
        if (num_items_read == EOF) return MM_PREMATURE_EOF;
    }
    while (num_items_read != 2);

    return 0;
}

int mm_write_mtx_array_size(FILE *f, int M, int N)
{
    if (fprintf(f, "%d %d\n", M, N) != 2)
        return MM_COULD_NOT_WRITE_FILE;
    else 
        return 0;
}



/*-------------------------------------------------------------------------*/

/******************************************************************/
/* use when I[], J[], and val[]J, and val[] are already allocated */
/******************************************************************/

int mm_read_mtx_crd_data(FILE *f, int M, int N, int nz, int I[], int J[],
        double val[], char* matcode)
{
    int i;
    if (mm_is_complex(matcode))
    {
        for (i=0; i<nz; i++)
            if (fscanf(f, "%d %d %lg %lg", &I[i], &J[i], &val[2*i], &val[2*i+1])
                != 4) return MM_PREMATURE_EOF;
    }
    else if (mm_is_real(matcode))
    {
        for (i=0; i<nz; i++)
        {
            if (fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i])
                != 3) return MM_PREMATURE_EOF;

        }
    }

    else if (mm_is_pattern(matcode))
    {
        for (i=0; i<nz; i++)
            if (fscanf(f, "%d %d", &I[i], &J[i])
                != 2) return MM_PREMATURE_EOF;
    }
    else
        return MM_UNSUPPORTED_TYPE;

    return 0;
        
}

int mm_read_mtx_crd_entry(FILE *f, int *I, int *J,
        double *real, double *imag, char* matcode)
{
    if (mm_is_complex(matcode))
    {
            if (fscanf(f, "%d %d %lg %lg", I, J, real, imag)
                != 4) return MM_PREMATURE_EOF;
    }
    else if (mm_is_real(matcode))
    {
            if (fscanf(f, "%d %d %lg\n", I, J, real)
                != 3) return MM_PREMATURE_EOF;

    }

    else if (mm_is_pattern(matcode))
    {
            if (fscanf(f, "%d %d", I, J) != 2) return MM_PREMATURE_EOF;
    }
    else
        return MM_UNSUPPORTED_TYPE;

    return 0;
        
}


/************************************************************************
    mm_read_mtx_crd()  fills M, N, nz, array of values, and return
                        type code, e.g. 'MCRS'

                        if matrix is complex, values[] is of size 2*nz,
                            (nz pairs of real/imaginary values)
************************************************************************/

int mm_read_mtx_crd(char *fname, int *M, int *N, int *nz, int **I, int **J, 
        double **val, char** matcode)
{
    int ret_code;
    FILE *f;

    if (strcmp(fname, "stdin") == 0) f=stdin;
    else
    if ((f = fopen(fname, "r")) == NULL)
        return MM_COULD_NOT_READ_FILE;


    if ((ret_code = mm_read_banner(f, matcode)) != 0)
        return ret_code;

    if (!(mm_is_valid(*matcode) && mm_is_sparse(*matcode) && 
            mm_is_matrix(*matcode)))
        return MM_UNSUPPORTED_TYPE;

    if ((ret_code = mm_read_mtx_crd_size(f, M, N, nz)) != 0)
        return ret_code;


    *I = (int *)  malloc(*nz * sizeof(int));
    *J = (int *)  malloc(*nz * sizeof(int));
    *val = NULL;

    if (mm_is_complex(*matcode))
    {
        *val = (double *) malloc(*nz * 2 * sizeof(double));
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val, 
                *matcode);
        if (ret_code != 0) return ret_code;
    }
    else if (mm_is_real(*matcode))
    {
        *val = (double *) malloc(*nz * sizeof(double));
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val, 
                *matcode);
        if (ret_code != 0) return ret_code;
    }

    else if (mm_is_pattern(*matcode))
    {
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val, 
                *matcode);
        if (ret_code != 0) return ret_code;
    }

    if (f != stdin) fclose(f);
    return 0;
}

int mm_write_banner(FILE *f, char* matcode)
{
    char *str = mm_typecode_to_str(matcode);

    fprintf(f, "%s ", MatrixMarketBanner);
    fprintf(f, "%s\n", str);
    free(str);

    /* BEGIN MFH 25 June 2004 */
    return 0;
    /* END MFH 25 June 2004 */
}

int mm_write_mtx_crd(char fname[], int M, int N, int nz, int I[], int J[],
        double val[], char* matcode)
{
    FILE *f;
    int i;

    if (strcmp(fname, "stdout") == 0) 
        f = stdout;
    else
    if ((f = fopen(fname, "w")) == NULL)
        return MM_COULD_NOT_WRITE_FILE;
    
    /* print banner followed by typecode */
    fprintf(f, "%s ", MatrixMarketBanner);
    fprintf(f, "%s\n", mm_typecode_to_str(matcode));

    /* print matrix sizes and nonzeros */
    fprintf(f, "%d %d %d\n", M, N, nz);

    /* print values */
    if (mm_is_pattern(matcode))
        for (i=0; i<nz; i++)
            fprintf(f, "%d %d\n", I[i], J[i]);
    else
    if (mm_is_real(matcode))
        for (i=0; i<nz; i++)
            fprintf(f, "%d %d %20.16g\n", I[i], J[i], val[i]);
    else
    if (mm_is_complex(matcode))
        for (i=0; i<nz; i++)
            fprintf(f, "%d %d %20.16g %20.16g\n", I[i], J[i], val[2*i], 
                        val[2*i+1]);
    else
    {
        if (f != stdout) fclose(f);
        return MM_UNSUPPORTED_TYPE;
    }

    if (f !=stdout) fclose(f);

    return 0;
}
    

/*
 * mfh 20 Jul 2008: this causes errors in "release" mode, so I
 * commented it out. 
 */
#if 0
/* mfh 1 Mar 2006: Rich Vuduc pointed out to me that with gcc 3.4.3
 * as set up on the Citris (Itanium 2) cluster, strdup isn't in string.h.
 * So for type safety, we need to declare it here. */
extern char* strdup (const char* s);
#endif /* 0 */

char  *mm_typecode_to_str(char* matcode)
{
    char buffer[MM_MAX_LINE_LENGTH];
    char *types[4];
    int error =0;

    /* mfh 26 May 2005: the compiler complains about types not being 
     * initialized.  Let's do that here. */
    types[0] = NULL;
    types[1] = NULL;
    types[2] = NULL;
    types[3] = NULL;

    /* check for MTX type */
    if (mm_is_matrix(matcode)) 
        types[0] = MM_MTX_STR;
    else
        error=1;

    /* check for CRD or ARR matrix */
    if (mm_is_sparse(matcode))
        types[1] = MM_SPARSE_STR;
    else
    if (mm_is_dense(matcode))
        types[1] = MM_DENSE_STR;
    else
        return NULL;

    /* check for element data type */
    if (mm_is_real(matcode))
        types[2] = MM_REAL_STR;
    else
    if (mm_is_complex(matcode))
        types[2] = MM_COMPLEX_STR;
    else
    if (mm_is_pattern(matcode))
        types[2] = MM_PATTERN_STR;
    else
    if (mm_is_integer(matcode))
        types[2] = MM_INT_STR;
    else
        return NULL;


    /* check for symmetry type */
    if (mm_is_general(matcode))
        types[3] = MM_GENERAL_STR;
    else
    if (mm_is_symmetric(matcode))
        types[3] = MM_SYMM_STR;
    else 
    if (mm_is_hermitian(matcode))
        types[3] = MM_HERM_STR;
    else 
    if (mm_is_skew(matcode))
        types[3] = MM_SKEW_STR;
    else
        return NULL;

    sprintf(buffer,"%s %s %s %s", types[0], types[1], types[2], types[3]);
    return strdup(buffer);

}
