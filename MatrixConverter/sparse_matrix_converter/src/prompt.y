/**
 * @file prompt.y
 * @author Mark Hoemmen
 * @since 05 Jul 2005
 * @date Time-stamp: <2008-07-16 11:24:27 mhoemmen>
 *
 * Simple shell for testing sparse matrix conversion.  It accepts the following 
 * commands:
 *
 * load <filename> <file-format>;
 *   Assumes that the given sparse matrix file <filename> is in the given 
 *   format <file-format>, and loads the file into the current sparse matrix.
 *   The shell only works with one sparse matrix at a time; it is an implicit
 *   argument of all the shell commands.  Supported file formats are:  
 *   `HARWELL_BOEING', `MATRIX_MARKET' and `MATLAB'.
 *
 * save <filename> <file-format>; 
 *   Save the current sparse matrix to the given filename in the given format.
 *   Supported file formats are: `HARWELL_BOEING', `MATRIX_MARKET' and 
 *   `MATLAB'.
 *
 * format?
 *   Prints the internal storage format of the current sparse matrix.
 *
 * convert <internal-format>;
 *   Converts the current sparse matrix from its current internal format to a 
 *   different internal format.  Supported internal formats are:  `CSC', `CSR',
 *   `COO', `BCOO', `BCSR', `JAD'.
 *
 * exit;
 *   Frees the current sparse matrix and exits the shell.
 *
 * Copyright (c) 2006, Regents of the University of California 
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

%{
#include "config.h"
#include "sparse_matrix.h"
#include "sparse_matrix_ops.h"

#include <get_options.h>
#include <smvm_malloc.h>
#include <smvm_util.h>
#include <timer.h>

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

struct sparse_matrix_t* A = NULL;

int
yylex (void);

void
yyerror (char const * s)
{
  fprintf (stderr, "*** ERROR: %s ***\n", s);
}

double seconds = 0.0;
%}

/* Bison declarations */

%union {
  char* str;  /* For reading strings */
}


/* Terminal symbols */
%token <str> STRING 
%token LOAD SAVE FORMAT CONVERT EXIT

%destructor { fprintf (stderr, "Calling destructor on symbol %s\n", $$); free ($$); } STRING 


%%

input:  /* empty */
	| input line          { printf ("$ "); }

line:   '\n'                  
	| stmtlist '\n'   

stmtlist:  /* empty */
	| stmtlist stmt 
	;

/* filename first, then format */
stmt:   LOAD STRING STRING ';' { if (A) destroy_sparse_matrix (A);
                               seconds = get_seconds ();
                               A = load_sparse_matrix (sparse_matrix_file_format_string_to_enum ($3), $2);
			       seconds = get_seconds () - seconds;
			       printf ("Loading the sparse matrix took %g seconds.\n", seconds);
			     }
/* filename first, then format */
	| SAVE STRING STRING ';' { seconds = get_seconds ();
	                         save_sparse_matrix ($2, A, sparse_matrix_file_format_string_to_enum ($3));
			         seconds = get_seconds () - seconds;
				 printf ("Saving the sparse matrix took %g seconds.\n", seconds);
			       }
	| FORMAT               { if (A == NULL)
	                           printf ("Sparse matrix is NULL!\n");
			         else
				   printf ("%s\n", sparse_matrix_format_string (A));
			       }
	| CONVERT STRING ';'   { seconds = get_seconds ();
	                         sparse_matrix_convert (A, sparse_matrix_storage_format_string_to_enum ($2));
				 seconds = get_seconds () - seconds;
				 printf ("Converting the sparse matrix took %g seconds.\n", seconds);
                               }
	| EXIT ';'             { if (A) destroy_sparse_matrix (A);
	                         printf ("Goodbye!\n");
				 exit (EXIT_SUCCESS);
			       }
	;

%%


int
main (int argc, char** argv)
{
  return yyparse ();
}

