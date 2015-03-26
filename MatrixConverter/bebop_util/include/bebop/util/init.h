#ifndef _init_h
#define _init_h
/****************************************************************************
 * @file init.h
 * @author Mark Hoemmen
 * @since 23 Nov 2007
 * @date Time-stamp: <2008-07-16 10:17:52 mhoemmen>
 *
 * Declaration of BeBOP Utility Library initialization functions.
 *
 * @note Moved out of util.h into this file on 23 Nov 2007.
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
#include <stdio.h>   /* FILE* declaration */

/**
 * Does sensible minimal default initializations, in this order:
 *
 * - bebop_set_debug_level_from_environment()
 * - bebop_start_logging()
 *
 * @note If you use this function, don't call the above functions 
 * in your code.
 */
void
bebop_default_initialize (int argc, char** argv, int* const info);

/** 
 * Reads the value of BEBOP_DEBUG_LEVEL from the environment and sets the
 * BEBOP_DEBUG_LEVEL variable accordingly.  If the debug level is > 0,
 * then WITH_DEBUG is operative.  If the debug level is > 1, then
 * WITH_DEBUG2 is operative, etc...
 */
void
bebop_set_debug_level_from_environment ();

/**
 * Returns the debug level.
 */
int
bebop_debug_level ();

/**
 * Sets the debug level.
 */
void 
bebop_set_debug_level (const int level);

/**
 * Special wrapper for exit(int) that shuts down BeBOP Utility Library services.
 */
void
bebop_exit (const int errcode);


#endif /* #ifndef _init_h */
