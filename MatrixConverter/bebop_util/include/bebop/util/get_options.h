#ifndef _get_options_h
#define _get_options_h
/**
 * @file get_options.h
 * @author Mark Hoemmen
 * @date Time-stamp: <2009-05-16 12:44:20 mhoemmen>
 *
 * Functions for parsing and converting command-line arguments.
 *
 * @note If USE_PTHREADS is defined (see bebop_make/options), these
 * functions will be thread-safe.
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
#include <stdio.h>
#include <string.h>

/* Forward declaration */
struct arginfo;

/**
 * Specifies the type of a command-line argument value.
 *
 * @see arginfo
 */
typedef enum { NULLARG, INTARG, DOUBLEARG, STRINGARG } arg_type_t;


/**
 * Returns a pointer to the arginfo object in the given list with the given 
 * label, or NULL if none such object is found in the list.
 */
struct arginfo*
find_arginfo (struct arginfo *list, const char c);


/**
 * Returns nonzero if the given option c was selected at the command line, 
 * else returns zero.
 */
int
got_arg_p (struct arginfo *list, const char c);


/**
 * Returns nonzero if the given option c was selected at the command line
 * and a value was specified there, else returns zero.
 */
int
got_arg_value_p (struct arginfo *list, const char c);


/**
 * Appends a new command-line argument specification to the given list of 
 * command-line argument specifications, using deep copies of the given
 * parameters.
 *
 * @warn It's critical that the initial value of arginfo_list, before you 
 * register any arguments, is NULL!  Otherwise the list will not be 
 * initialized properly.
 *
 * @param arg [OUT]  list of arginfo objects (NULL if list is empty)
 *
 * @param c [IN]     Command-line option name (e.g. 't' for "-t <blah>")
 *
 * @param type [IN]  Type of the command-line argument
 *
 * @param val [IN]   Pointer to the default value of the command-line argument
 *                   e.g. if type==DOUBLEARG or INTARG, val is a pointer to a 
 *                   variable holding the initial value.  If type==STRINGARG, 
 *                   then val is a char* holding the string.  In all of these 
 *                   cases, value will be copied, and val must be a valid 
 *                   pointer.  (You must always set a default value for any
 *                   parameter.)  If type==NULLARG, val will not be referenced.  
 *
 * @param desc [IN]  English description of the command-line argument.  String
 *                   will be (deeply) copied.
 *
 * @return Pointer to the head of the arginfo list
 */
struct arginfo*
register_arginfo (struct arginfo *arginfo_list, 
                  const char c, 
                  const arg_type_t type, 
                  const void *val, 
                  const char *desc);

/**
 * Registers the given function as the user-specified usage function.
 *
 * @param usage_function [IN]  The usage function
 *
 * @note If this function is not called, a default usage function is used.
 */
void 
register_usage_function (void (usage_function) (FILE*, const struct arginfo*, const struct arginfo*));

/**
 * Returns the integer value stored in the given arginfo object.
 * Aborts if arg's value is not of type integer.  
 */
int 
get_int_argval (const struct arginfo* arg);


/**
 * Returns the double value stored in the given arginfo object.
 * Aborts if arg's value is not of type double.  
 */
double
get_double_argval (const struct arginfo* arg);

/**
 * Returns a pointer to the string stored in the given arginfo object.
 * Aborts if arg's value is not a string.
 */
char*
get_string_argval (const struct arginfo* arg);

/**
 * Sets arg's value to val, interpreted as an integer.
 */
void 
set_int_argval (struct arginfo* arg, const int val);

/**
 * Sets arg's value to val, interpreted as a double.
 */
void
set_double_argval (struct arginfo* arg, const double val);

/**
 * Sets arg's value to val, interpreted as a string (char*).
 *
 * @param arg [OUT] receives a deep copy of the given string
 * @param val [IN]  null-terminated character array
 */
void
set_string_argval (struct arginfo* arg, char* val);


/** 
 * Processes the command line options.  Call this function exactly once.
 * After calling it, "optind" (see <unistd.h>) is the index in argv of the 
 * first argv element that is not an option.
 *
 * @param argc [IN]
 * @param argv [IN/OUT]
 * @param core_args [IN/OUT]  Describes and stores the values of the 
 *                            ``core'' (primary) command-line arguments 
 * @param ext_args [IN/OUT]  If non-NULL, describes and stores the values
 *                           of ``extended'' command-line arguments.
 */
void
get_options (const int argc, 
             char* argv[], 
             struct arginfo *core_args, 
             struct arginfo *ext_args);

/**
 * Dumps the standard "usage message" to the given output stream.
 * If a non-default usage function was specified, control is passed to it.
 *
 * @param out [OUT]      Output stream to which to dump the usage message
 * @param pathname [IN]  Name of the executable file of this program
 * @param core_args [IN] List of "core" (primary) command-line arguments
 * @param ext_args [IN]  List of "extended" command-line arguments
 */
void
dump_usage (FILE* out, 
            const char *const pathname, 
            struct arginfo *core_args, 
            struct arginfo *ext_args);

/**
 * Deallocates all the registered command-line argument specifications in the 
 * given list of them.
 *
 * @param list [IN/OUT]  Head of the arginfo list
 */
void
destroy_arginfo_list (struct arginfo *list);



#endif /* _get_options_h */
