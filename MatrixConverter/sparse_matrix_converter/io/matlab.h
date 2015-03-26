#ifndef _matlab_h
#define _matlab_h
/**
 * @file matlab.h
 * @author mfh
 * @since 09 May 2007
 * @date Time-stamp: <2008-07-16 11:00:17 mhoemmen>
 *
 * Basic I/O for Matlab-format sparse matrix files.
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
#include <stdio.h>
#include <enumerations.h> /* value_type_t */

/**
 * Read a Matlab-format sparse matrix from file "infile" and 
 * store its data.
 *
 * @param info  String containing tuning hints for optimizing
 *   reading in the sparse matrix.  This is a comma-delimited
 *   list of case-insensitive VARIABLE=VALUE pairs.  Valid 
 *   variables are nelts (whose value is a nonnegative integer 
 *   -- guess for number of nonzeros in the sparse matrix) and 
 *   save_memory_flag (which is a boolean -- "True" or "False").
 *   If save_memory_flag=True and nelts is not given, we read 
 *   the matrix file twice -- once to determine the number of 
 *   nonzeros, and again to read in the nonzeros' values.
 *   Otherwise, we read in the matrix file once, expanding the 
 *   arrays as necessary via a doubling scheme that may waste
 *   memory (so that we only need to read the file once).  
 *   (Recall that Matlab format does not explicitly store the 
 *   number of nonzeros -- we have to count them up ourselves.)
 *
 * @return Zero if successful, else nonzero error code.
 */
int
read_matlab_sparse_matrix (int* nelts, 
			   value_type_t* valtype,
			   int** __II, int** __JJ, void** __val, 
			   const char* filename,
			   const char* info);

#endif /* _matlab_h */
