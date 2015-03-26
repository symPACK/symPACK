#ifndef _smvm_malloc_h
#define _smvm_malloc_h
/****************************************************************************
 * @file smvm_malloc.h
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-16 10:18:56 mhoemmen>
 *
 * @brief Wrappers for malloc and realloc.
 *
 * Wrappers for malloc and realloc, that check to ensure successful
 * allocation of memory, and abort the process with a useful error
 * message if allocation was not successful.
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
 ****************************************************************************/

#include <bebop/util/config.h>

#ifdef NO_STDDEF_H
  typedef int size_t;
#else
  #include <stddef.h>  /* size_t */
#endif

/**
 * @brief Checked malloc that aborts if allocation failed.
 *
 * Implementation of checked malloc that aborts if allocation failed.  
 * @warn You should use the macro bebop_malloc instead.
 * 
 * @see bebop_malloc()
 *
 * @param size       How many bytes to allocate.
 * @param srcfile    Pass it __FILE__ so if fails, it can output source file 
 *                   name.
 * @param linenum    Pass it __LINE__ so if fails, it can output the line in
 *                   the source file at which it failed.
 * 
 * @return  Non-NULL, valid pointer to allocated memory.  Client is 
 *          responsible for `free'-ing it.
 */
void*
mfh_malloc (size_t size, char* srcfile, int linenum);

/**
 * @brief Checked malloc that aborts if allocation failed.  
 *
 * A convenient wrapper that should be used in place of mfh_malloc.
 * @see mfh_malloc()
 *
 * @param size       How many bytes to allocate.
 * @return  Valid pointer to allocated memory.  Client is responsible for 
 *          `free'-ing it.
 */
#define bebop_malloc( size )   mfh_malloc (size, __FILE__, __LINE__)


/**
 * @brief Checked calloc that aborts if allocation failed.
 *
 * Checked calloc that aborts if allocation failed.  You should use the macro
 * bebop_calloc instead.
 *
 * @see bebop_calloc()
 *
 * @param num_elements  Number of elements in the array.
 * @param element_size  Size of each element of the array.
 * @param srcfile    Pass it __FILE__ so if fails, it can output source file 
 *                   name.
 * @param linenum    Pass it __LINE__ so if fails, it can output the line in
 *                   the source file at which it failed.
 * 
 * @return  Non-NULL, valid pointer to allocated, zeroed memory.  Client is 
 *          responsible for `free'-ing it.
 */
void*
mfh_calloc (size_t num_elements, size_t element_size, char* srcfile, 
	   int linenum);

/**
 * @brief Checked calloc that aborts if allocation failed.  
 *
 * A convenient wrapper that should be used in place of mfh_calloc.
 * @see mfh_calloc()
 *
 * @param num_elements  Number of elements in the array.
 * @param element_size  Size of each element of the array.
 * @return  Valid pointer to allocated, zeroed memory.  Client is responsible 
 *          for `free'-ing it.
 */
#define bebop_calloc( num_elements, element_size )   \
  mfh_calloc (num_elements, element_size, __FILE__, __LINE__)


/**
 * @brief Checked realloc that aborts if allocation failed.
 *
 * Checked realloc that aborts if allocation failed.  You should use the macro
 * bebop_realloc (@see bebop_realloc()) instead.
 *
 * @param ptr        Memory location at which we wish to reallocate.
 * @param size       How much to allocate.
 * @param srcfile    Pass it __FILE__ so if fails, it can output source file name.
 * @param linenum    Pass it __LINE__ for same reason as above.
 */
void*
mfh_realloc (void* ptr, size_t size, char* srcfile, int linenum);

/**
 * A convenient wrapper that should be used in place of mfh_realloc.
 * @see mfh_realloc()
 * 
 * @param size   How many bytes to allocate.
 * @return  Valid pointer to allocated memory.  Client is responsible for 
 *          `free'-ing it.
 */
#define bebop_realloc( ptr, size )   mfh_realloc (ptr, size, __FILE__, __LINE__)


void
mfh_free (void* ptr, char* srcfile, int linenum);


#define bebop_free( ptr )  mfh_free (ptr, __FILE__, __LINE__)



#endif /* #ifndef _smvm_malloc_h */
