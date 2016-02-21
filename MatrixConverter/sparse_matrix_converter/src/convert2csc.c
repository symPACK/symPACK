/**
 * @file main.c
 * @author Mark Hoemmen
 * @since 26 May 2005
 * @date Time-stamp: <2008-07-16 11:23:53 mhoemmen>
 *
 * Driver program for the sparse matrix converter.
 *
 * Copyright (c) 2008, Regents of the University of California 
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the
 * following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright 
 *   notice, this list of conditions and the following disclaimer in 
 *   the documentation and/or other materials provided with the 
 *   distribution.
 *
 * * Neither the name of the University of California, Berkeley, nor
 *   the names of its contributors may be used to endorse or promote
 *   products derived from this software without specific prior
 *   written permission.  
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include <bebop/util/enumerations.h>
#include <bebop/util/config.h>
#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/csc_matrix.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/sparse_matrix_ops.h>


#include <bebop/util/malloc.h>
#include <bebop/util/util.h>
#include <bebop/util/log.h>





#include <bebop/util/get_options.h>
#include <bebop/util/init.h>
#include <bebop/util/log.h>
#include <bebop/util/malloc.h>
#include <bebop/util/timer.h>
#include <bebop/util/util.h>

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#ifndef _SP_base
#define _SP_base 0
#endif

struct 
cmdlineopts_t
{
  char* input_filename;
  char* input_file_format;
  char* output_filename;
  char* output_file_format;
} opts;

  static void
usage (FILE* out, const struct arginfo* arglist, const struct arginfo* ext_args)
{
  fprintf (out, "Usage:\n");
  fprintf (out, "sparse_matrix_converter <in-filename> <in-format> <out-filename> <out-format> [options]\n");
  fprintf (out, "<in-filename>:   name of file containing the sparse matrix to read in\n");
  fprintf (out, "<in-format>:     format of the input file (\"HB\" for Harwell-Boeing, \"ML\" for\n");
  fprintf (out, "                 Matlab or \"MM\" for MatrixMarket)\n");
  fprintf (out, "<out-filename>:  name of file to which to output the sparse matrix\n");
  fprintf (out, "<out-format>:    format of the output file (\"HB\" for Harwell-Boeing, \"ML\" for\n");
  fprintf (out, "                 Matlab or \"MM\" for MatrixMarket)\n");
  fprintf (out, "[options]: -a -- validate the input matrix only, without outputting anything\n");
  fprintf (out, "           -e -- expand symmetric into unsymmetric storage (this option is\n");
  fprintf (out, "                 unnecessary for output to Matlab format, as Matlab format is already\n");
  fprintf (out, "                 expanded)\n");
  fprintf (out, "           -l  -- store lower triangular part only (for CSC output)\n");
  fprintf (out, "           -k  -- convert to skew symmetric storage (for CSC output)\n");
  fprintf (out, "           -s  -- specify that the matrix is symmetric (for CSC output)\n");
  fprintf (out, "           -c  -- specify that the matrix is complex (for CSC output)\n");
  fprintf (out, "           -d  -- turn the matrix into spd\n");
  fprintf (out, "           -v  -- verbose mode\n");
  fprintf (out, "           -h  -- print this usage message and exit\n\n");
}

/**
 * Perform the matrix validation operation specified by the "-a" command-line option.
 */
  static int
do_validate_matrix (int argc, char *argv[], struct arginfo* arglist)
{
  extern int optind;

  enum sparse_matrix_file_format_t informat = 0;
  struct sparse_matrix_t* A = NULL;
  double seconds = 0.0;
  int errcode = 0;

  if (argc - optind != 2)
  {
    fprintf (stderr, "*** Incorrect number of command-line arguments: %d ar"
        "e specified, but there should be %d ***\n", argc - optind, 2);
    dump_usage (stderr, argv[0], arglist, NULL);
    return -1;
  }

  opts.input_filename = argv[optind];
  opts.input_file_format = argv[optind+1];

  if (strcmp (opts.input_file_format, "HB") == 0 || 
      strcmp (opts.input_file_format, "hb") == 0)
    informat = HARWELL_BOEING;
  else if (strcmp (opts.input_file_format, "MM") == 0 ||
      strcmp (opts.input_file_format, "mm") == 0)
    informat = MATRIX_MARKET;
  else if (strcmp (opts.input_file_format, "ML") == 0 ||
      strcmp (opts.input_file_format, "ml") == 0)
    informat = MATLAB;
  else
  {
    fprintf (stderr, "*** Unsupported input file format \"%s\" ***\n", 
        opts.input_file_format);
    dump_usage (stderr, argv[0], arglist, NULL);
    return -1;
  }

  if (got_arg_p (arglist, 'v'))
  {
    printf ("Loading sparse matrix...");
    fflush (stdout); /* Otherwise the message may not be printed until 
                        after loading is complete */
  }
  seconds = get_seconds();
  A = load_sparse_matrix (informat, opts.input_filename);
  seconds = get_seconds() - seconds;
  if (A == NULL)
  {
    fprintf (stderr, "*** Failed to load input matrix file \"%s\" ***\n", 
        opts.input_filename);
    destroy_arginfo_list (arglist);
    return 1;
  }
  if (got_arg_p (arglist, 'v'))
    printf ("done, in %g seconds\n", seconds);


  if (got_arg_p (arglist, 'v'))
  {
    printf ("Validating sparse matrix...");
    fflush (stdout); 
  }
  seconds = get_seconds();
  errcode = valid_sparse_matrix (A);
  seconds = get_seconds() - seconds;
  if (got_arg_p (arglist, 'v'))
    printf ("done, in %g seconds\n", seconds);

  if (valid_sparse_matrix (A))
    printf ("\n### Sparse matrix is valid ###\n\n");
  else 
    printf ("\n### Invalid sparse matrix! ###\n\n");

  destroy_sparse_matrix (A);
  return 0;
}

int save_sparse_matrix_csc(const char * filename, struct sparse_matrix_t* Ain, int toSM, int toLT )
{

  int i;
  FILE *out_file;
  struct csc_matrix_t* A;


  if (Ain->format == CSR)
  {
    A = (struct csc_matrix_t*)csr_to_csc((struct csr_matrix_t*)Ain->repr);
  }
  else if (Ain->format == CSC)
  {
    A = (struct csc_matrix_t*)Ain->repr;
  }
  else if (Ain->format == COO){
    A = (struct csc_matrix_t*)coo_to_csc((struct coo_matrix_t*)Ain->repr);
  }
  else{
    fprintf(stderr,"NOT SUPPORTED\n");
  }

  fprintf(stderr,"m=%d n=%d\n",A->m,A->n);
  fprintf(stderr,"symmetry=%d\n",A->symmetry_type);

  assert( (A->m == A->n) );//&& (A->symmetry_type == SYMMETRIC));

  if ( filename != NULL ) {
    if ( (out_file = fopen( filename, "wb")) == NULL ) {
      fprintf(stderr,"Error: Cannot open file: %s\n",filename);
      return 0;
    }
  } else out_file = stdout;





  int offset = 1-_SP_base;  /* if base 0 storage is declared (via macro definition), */
  /* then storage entries are offset by 1                  */


  int np1 = A->n+1;


  /*if only LT, recompute everything*/
  if(toLT){
  fwrite(&A->n,sizeof(A->n),1,out_file);
    int nnz2=0;
    for (i=0;i<np1;i++)
    {
      int rowid;
      int fi = A->colptr[i];
      int li = A->colptr[i+1]-1;

      int rc = 0;
      for(rowid = fi; rowid<=li;++rowid){
        int row = A->rowidx[rowid];
        if(row>=i) ++rc;
      }
      nnz2+=rc;
    }
    fwrite(&nnz2,sizeof(nnz2),1,out_file);

    /*  Print column pointers:   */
    /*compute new column pointer*/
    fwrite(&np1,sizeof(int),1,out_file);
    int rc;
    rc = offset;
    for (i=0;i<np1;i++)
    {
      int rowid;
      int fi = A->colptr[i];
      int li = A->colptr[i+1]-1;

      fwrite(&rc,sizeof(int),1,out_file);
      
      for(rowid = fi; rowid<=li;++rowid){
        int row = A->rowidx[rowid];
        if(row>=i) ++rc;
      }
    }
    //fwrite(&rc,sizeof(int),1,out_file);

    /*  Print row indices:       */
    fwrite(&nnz2,sizeof(nnz2),1,out_file);
    for (i=0;i<np1;i++)
    {
      int rowid;
      int fi = A->colptr[i];
      int li = A->colptr[i+1]-1;

      for(rowid = fi; rowid<=li;++rowid){
        int row = A->rowidx[rowid];
        if(row>=i){
          int entry = row+offset;
          fwrite(&entry,sizeof(int),1,out_file);
        };
      }
    }

    /*  Print values:            */

    fwrite(&nnz2,sizeof(nnz2),1,out_file);
    if(A->value_type == REAL){
      double * val = (double*)A->values;
      for (i=0;i<np1;i++)
      {
        int rowid;
        int fi = A->colptr[i];
        int li = A->colptr[i+1]-1;

        for (rowid=fi;rowid<=li;rowid++)
        {
          int row = A->rowidx[rowid];
          if(row>=i){
            fwrite(&val[rowid],sizeof(double),1,out_file);
          }
        }
      }
    }
    else if(A->value_type == COMPLEX){
      double _Complex * val = (double _Complex*)A->values;
      for (i=0;i<np1;i++)
      {
        int rowid;
        int fi = A->colptr[i];
        int li = A->colptr[i+1]-1;

        for (rowid=fi;rowid<=li;rowid++)
        {
          int row = A->rowidx[rowid];
          if(row>=i){
            fwrite(&val[rowid],sizeof(double _Complex),1,out_file);
          }
        }
      }
    }


  }
  else{

  fwrite(&A->n,sizeof(A->n),1,out_file);
  fwrite(&A->nnz,sizeof(A->nnz),1,out_file);

    /*  Print column pointers:   */
    fwrite(&np1,sizeof(int),1,out_file);
    for (i=0;i<np1;i++)
    {
      int entry = A->colptr[i]+offset;
      //   fprintf(stderr,"colptr[%d]=%d\n",i,entry);
      fwrite(&entry,sizeof(int),1,out_file);
    }


    fwrite(&A->nnz,sizeof(A->nnz),1,out_file);
    /*  Print row indices:       */
    for (i=0;i<A->nnz;i++)
    {
      int entry = A->rowidx[i]+offset;
      fwrite(&entry,sizeof(int),1,out_file);
    }

    fwrite(&A->nnz,sizeof(A->nnz),1,out_file);

    /*  Print values:            */

    if(A->value_type == REAL){
      double * val = (double*)A->values;
      if(toSM==1){
        for (i=0;i<np1;i++)
        {
          int rowid;
          int fi = A->colptr[i];
          int li = A->colptr[i+1]-1;

          for (rowid=fi;rowid<=li;rowid++)
          {
            int row = A->rowidx[rowid];
            if(i>row){
              double entry = -val[rowid];
              fwrite(&entry,sizeof(double),1,out_file);
            }
            else{
              fwrite(&val[rowid],sizeof(double),1,out_file);
            }
          }
        }
      }
      else{
        for (i=0;i<A->nnz;i++)
        {
          fwrite(&val[i],sizeof(double),1,out_file);
        }
      }
    }
    else if(A->value_type == COMPLEX){
      double _Complex * val = (double _Complex*)A->values;
      if(toSM==1){
        for (i=0;i<np1;i++)
        {
          int rowid;
          int fi = A->colptr[i];
          int li = A->colptr[i+1]-1;

          for (rowid=fi;rowid<=li;rowid++)
          {
            int row = A->rowidx[rowid];
            if(i>row){
              double _Complex entry = -val[rowid];
              fwrite(&entry,sizeof(double _Complex),1,out_file);
            }
            else{
              fwrite(&val[rowid],sizeof(double _Complex),1,out_file);
            }
          }
        }
      }
      else{
        for (i=0;i<A->nnz;i++)
        {
          fwrite(&val[i],sizeof(double _Complex),1,out_file);
        }
      }    }
  }
  fclose(out_file);

  return 0;
}

struct sparse_matrix_t * read_sparse_matrix_csc(const char * filename, int isSymmetric)
{
  struct sparse_matrix_t * Aout;
  int i;
  FILE *in_file;
  struct csc_matrix_t* A;
  int n, nnz;
  double * values;
  int * colptr,*rowidx;


  if ( filename != NULL ) {
    if ( (in_file = fopen( filename, "rb")) == NULL ) {
      fprintf(stderr,"Error: Cannot open file: %s\n",filename);
      return 0;
    }
  } else in_file = stdout;



  fread(&n,sizeof(int),1,in_file);
  fread(&nnz,sizeof(int),1,in_file);


  fprintf(stderr,"n is %d, nnz is %d\n",n,nnz);

  int offset = 1-_SP_base;  /* if base 0 storage is declared (via macro definition), */
  /* then storage entries are offset by 1                  */

  int np1 = 0;
  fread(&np1,sizeof(int),1,in_file);
  colptr = (int *)malloc(np1*sizeof(int));
  /*  Print column pointers:   */
  for (i=0;i<np1;i++)
  {
    int entry;
    fread(&entry,sizeof(int),1,in_file);
    colptr[i] =  entry -offset;
  }


  fread(&np1,sizeof(int),1,in_file);
  rowidx = (int *)malloc(np1*sizeof(int));
  /*  Print row indices:       */
  for (i=0;i<nnz;i++)
  {
    int entry;
    fread(&entry,sizeof(int),1,in_file);
    rowidx[i]=entry-offset;
  }


  /*  Print values:            */
  fread(&np1,sizeof(int),1,in_file);
  values = (void *)malloc(np1*sizeof(double));

  for (i=0;i<nnz;i++)
  {
    fread(&values[i],sizeof(double),1,in_file);
  }



  fclose(in_file);

  fprintf(stderr,"Matrix in memory\n");
  A = create_csc_matrix (n, n, nnz, 
      values, rowidx, colptr,
      (isSymmetric==1?SYMMETRIC:UNSYMMETRIC),
      0,
      REAL,
      LIBRARY_DEALLOCATES,
      &free, NO_COPY);

  Aout = create_sparse_matrix (CSC, A);
  return Aout;
}





struct sparse_matrix_t * read_complex_sparse_matrix_csc(const char * filename, int isSymmetric)
{
  struct sparse_matrix_t * Aout;
  int i;
  FILE *in_file;
  struct csc_matrix_t* A;
  int n, nnz;
  double _Complex * values;
  int * colptr,*rowidx;


  if ( filename != NULL ) {
    if ( (in_file = fopen( filename, "rb")) == NULL ) {
      fprintf(stderr,"Error: Cannot open file: %s\n",filename);
      return 0;
    }
  } else in_file = stdout;



  fread(&n,sizeof(int),1,in_file);
  fread(&nnz,sizeof(int),1,in_file);


  fprintf(stderr,"n is %d, nnz is %d\n",n,nnz);

  int offset = 1-_SP_base;  /* if base 0 storage is declared (via macro definition), */
  /* then storage entries are offset by 1                  */

  int np1 = 0;
  fread(&np1,sizeof(int),1,in_file);
  colptr = (int *)malloc(np1*sizeof(int));
  /*  Print column pointers:   */
  for (i=0;i<np1;i++)
  {
    int entry;
    fread(&entry,sizeof(int),1,in_file);
    colptr[i] =  entry -offset;
  }


  fread(&np1,sizeof(int),1,in_file);
  rowidx = (int *)malloc(np1*sizeof(int));
  /*  Print row indices:       */
  for (i=0;i<nnz;i++)
  {
    int entry;
    fread(&entry,sizeof(int),1,in_file);
    rowidx[i]=entry-offset;
  }


  /*  Print values:            */
  fread(&np1,sizeof(int),1,in_file);
  values = (void *)malloc(np1*sizeof(double _Complex));

  for (i=0;i<nnz;i++)
  {
    fread(&values[i],sizeof(double _Complex),1,in_file);
  }



  fclose(in_file);

  fprintf(stderr,"Matrix in memory\n");
  A = create_csc_matrix (n, n, nnz, 
      values, rowidx, colptr,
      (isSymmetric==1?SYMMETRIC:UNSYMMETRIC),
      0,
      COMPLEX,
      LIBRARY_DEALLOCATES,
      &free, NO_COPY);

  Aout = create_sparse_matrix (CSC, A);
  return Aout;
}











  int
main (int argc, char *argv[])
{
  extern int optind;

  enum sparse_matrix_file_format_t informat = 0;
  enum sparse_matrix_file_format_t outformat = 0;
  struct sparse_matrix_t* A = NULL;
  struct arginfo *arglist = NULL;
  double seconds = 0.0;
  int errcode = 0;

  bebop_default_initialize (argc, argv, &errcode);
  if (errcode != 0)
  {
    fprintf (stderr, "*** Failed to initialize BeBOP Utility Library "
        "(error code %d) ***\n", errcode);
    bebop_exit (EXIT_FAILURE);
  }

  /* Set the get_options usage function to "usage", instead of using the default
   * usage function.  This is necessary because the command-line arguments include
   * things that are not "options" in the strict sense, because they do not follow
   * a "-[a-z]".
   */
  register_usage_function (usage);

  arglist = register_arginfo (arglist, 'v', NULLARG, NULL, "If specified, ac"
      "tivate verbose mode");
  arglist = register_arginfo (arglist, 'e', NULLARG, NULL, "If specified, ex"
      "pand the input matrix from symmetric storage "
      "into unsymmetric storage");
  arglist = register_arginfo (arglist, 'a', NULLARG, NULL, "If specified, va"
      "lidate the input matrix, without outputting a"
      "nything");

  arglist = register_arginfo (arglist, 'l', NULLARG, NULL, "If specified, "
      "convert the ilower triangular symmetric storage. ONLY TO CSC matrices"
  );
  arglist = register_arginfo (arglist, 'k', NULLARG, NULL, "If specified, "
      "convert the matrix to skew symmetric storage. ONLY TO CSC matrices"
      );
  arglist = register_arginfo (arglist, 's', NULLARG, NULL, "If specified, "
      "specify that the CSC matrix is a symmetric matrix. ONLY TO CSC matrices"
      );
  arglist = register_arginfo (arglist, 'c', NULLARG, NULL, "If specified, "
      "specify that the CSC matrix is a complex matrix. ONLY TO CSC matrices"
      );

  arglist = register_arginfo (arglist, 'd', NULLARG, NULL, "If specified, "
      "specify that the CSC matrix has to be made SPD. ONLY TO CSC matrices"
      );
  get_options (argc, argv, arglist, NULL);

  if (got_arg_p (arglist, 'a'))
  {
    int errcode = do_validate_matrix (argc, argv, arglist);
    destroy_arginfo_list (arglist);
    deinit_timer();
    bebop_exit (errcode); /* stops logging */
  }





  if (argc - optind != 4)
  {
    fprintf (stderr, "*** Incorrect number of command-line arguments: %d ar"
        "e specified, but there should be %d ***\n", argc - optind, 4);
    dump_usage (stderr, argv[0], arglist, NULL);
    bebop_exit (EXIT_FAILURE); /* stops logging */
  }

  opts.input_filename = argv[optind];
  opts.input_file_format = argv[optind+1];
  opts.output_filename = argv[optind+2];
  opts.output_file_format = argv[optind+3];

  int isCSC_in = 0;
  if (strcmp (opts.input_file_format, "HB") == 0 || 
      strcmp (opts.input_file_format, "hb") == 0)
    informat = HARWELL_BOEING;
  else if (strcmp (opts.input_file_format, "MM") == 0 ||
      strcmp (opts.input_file_format, "mm") == 0)
    informat = MATRIX_MARKET;
  else if (strcmp (opts.input_file_format, "ML") == 0 ||
      strcmp (opts.input_file_format, "ml") == 0)
    informat = MATLAB;
  else if (strcmp (opts.input_file_format, "CSC") == 0 ||
      strcmp (opts.input_file_format, "csc") == 0)
    isCSC_in = 1;
  else
  {
    fprintf (stderr, "*** Unsupported input file format \"%s\" ***\n", 
        opts.input_file_format);
    dump_usage (stderr, argv[0], arglist, NULL);
    destroy_arginfo_list (arglist);
    bebop_exit (EXIT_FAILURE); /* stops logging */
  }

  int isCSC = 0;

  if (strcmp (opts.output_file_format, "HB") == 0 || 
      strcmp (opts.output_file_format, "hb") == 0)
    outformat = HARWELL_BOEING;
  else if (strcmp (opts.output_file_format, "MM") == 0 ||
      strcmp (opts.output_file_format, "mm") == 0)
    outformat = MATRIX_MARKET;
  else if (strcmp (opts.output_file_format, "ML") == 0 ||
      strcmp (opts.output_file_format, "ml") == 0)
    outformat = MATLAB;
  else if (strcmp (opts.output_file_format, "CSC") == 0 ||
      strcmp (opts.output_file_format, "csc") == 0)
    isCSC = 1;
  else
  {
    fprintf (stderr, "*** Unsupported output file format \"%s\" ***\n", 
        opts.output_file_format);
    dump_usage (stderr, argv[0], arglist, NULL);
    destroy_arginfo_list (arglist);
    bebop_exit (EXIT_FAILURE); /* stops logging */
  }

  if (got_arg_p (arglist, 'v'))
  {
    printf ("Loading sparse matrix...");
    fflush (stdout); /* Otherwise the message may not be printed until 
                        after loading is complete */
  }
  seconds = get_seconds();


  if(isCSC_in){

    int isSymmetric = 0;
    if(got_arg_p (arglist, 's')){
      isSymmetric = 1;
    }
    if (got_arg_p(arglist, 'c')){
      A = read_complex_sparse_matrix_csc(opts.input_filename, isSymmetric);
    }
    else{
      A = read_sparse_matrix_csc(opts.input_filename, isSymmetric);
    }
  }
  else{
    A = load_sparse_matrix (informat, opts.input_filename);
  }

  fprintf(stderr,"FORMAT is %d\n",A->format); 
  if(A->format == CSC){
    struct csc_matrix_t * B = (struct csc_matrix_t*)A->repr;
    fprintf(stderr,"CSC m=%d n=%d nnz=%d\n",B->m,B->n,B->nnz);
    fprintf(stderr,"CSC colptr[0]=%d\n",B->colptr[0]);
  }
  else if(A->format == CSR){
    struct csr_matrix_t * B = (struct csr_matrix_t*)A->repr;
    fprintf(stderr,"CSR m=%d n=%d nnz=%d\n",B->m,B->n,B->nnz);
    fprintf(stderr,"CSR rowptr[0]=%d\n",B->rowptr[0]);
  }

  seconds = get_seconds() - seconds;
  if (A == NULL)
  {
    fprintf (stderr, "*** Failed to load input matrix file \"%s\" ***\n", 
        opts.input_filename);
    destroy_arginfo_list (arglist);
    bebop_exit (1); /* stops logging */
  }
  if (got_arg_p (arglist, 'v'))
    printf ("done, in %g seconds\n", seconds);

  if (got_arg_p (arglist, 'e') || got_arg_p (arglist, 'k'))
  {
    if (got_arg_p (arglist, 'v'))
    {
      printf ("Expanding sparse matrix into unsymmetric storage...");
      fflush (stdout);
    }

    seconds = get_seconds();
    errcode = sparse_matrix_expand_symmetric_storage (A);
    seconds = get_seconds() - seconds;
    if (errcode != 0)
    {
      fprintf (stderr, "*** Failed to expand matrix into symmetric storage ***\n");
      destroy_sparse_matrix (A);
      destroy_arginfo_list (arglist);
      bebop_exit (2);
    }
    if (got_arg_p (arglist, 'v'))
      printf ("done, in %g seconds\n", seconds);
  }


    if (got_arg_p (arglist, 'l'))
    {

      if (got_arg_p (arglist, 'v'))
      {
        printf ("Converting sparse matrix into symmetric storage...");
        fflush (stdout);
      }

      seconds = get_seconds();
      errcode = sparse_matrix_to_symmetric_storage (A);
      seconds = get_seconds() - seconds;
      if (errcode != 0)
      {
        fprintf (stderr, "*** Failed to convert matrix into symmetric storage ***\n");
        destroy_sparse_matrix (A);
        destroy_arginfo_list (arglist);
        bebop_exit (2);
      }
      if (got_arg_p (arglist, 'v'))
        printf ("done, in %g seconds\n", seconds);

    }




    if (got_arg_p (arglist, 'd'))
    {

      if (got_arg_p (arglist, 'v'))
      {
        printf ("Converting sparse matrix into SPD");
        fflush (stdout);
      }

      seconds = get_seconds();
      errcode = sparse_matrix_to_SPD (A);
      seconds = get_seconds() - seconds;
      if (errcode != 0)
      {
        fprintf (stderr, "*** Failed to convert matrix into SPD ***\n");
        destroy_sparse_matrix (A);
        destroy_arginfo_list (arglist);
        bebop_exit (2);
      }
      if (got_arg_p (arglist, 'v'))
        printf ("done, in %g seconds\n", seconds);

    }











  if (got_arg_p (arglist, 'v'))
  {
    printf ("Converting and saving sparse matrix...");
    fflush (stdout);
  }
  seconds = get_seconds();
  if(!isCSC){
    errcode = save_sparse_matrix (opts.output_filename, A, outformat);
  }
  else{
    if(got_arg_p (arglist, 'l')){
      save_sparse_matrix_csc(opts.output_filename, A,0,1);
    }
    else if(got_arg_p (arglist, 'k')){
      save_sparse_matrix_csc(opts.output_filename, A,1,0);
    }
    else{
      save_sparse_matrix_csc(opts.output_filename, A,0,0);
    }
  }
  seconds = get_seconds() - seconds;
  if (errcode != 0)
  {
    fprintf (stderr, "*** Failed to save output matrix to file \"%s\" ***\n",
        opts.output_filename);
  }
  if (got_arg_p (arglist, 'v'))
    printf ("done, in %g seconds\n", seconds);

  destroy_sparse_matrix (A);
  destroy_arginfo_list (arglist);
  deinit_timer();
  bebop_exit (errcode);
  return errcode; /* sentinel to pacify the compiler */
}



