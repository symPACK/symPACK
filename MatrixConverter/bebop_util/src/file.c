/**
 * @file file.c
 * @author Mark Hoemmen
 * @since 23 Nov 2007
 * @date Time-stamp: <2008-07-16 10:09:42 mhoemmen>
 *
 * Implementation of filesystem utility functions.  
 *
 * @note Moved out of util.c and into this file on 23 Nov 2007.
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
 *************************************************************************/
#include <bebop/util/config.h>

#include <bebop/util/log.h>
#include <bebop/util/malloc.h>
#include <bebop/util/string.h>
#include <bebop/util/util.h>

/* 
 * Note (mfh 23 Nov 2007) that errno is usually some kind of construct that
 * takes care of its own thread safety.  Darwin chooses the NeXT approach of
 * having errno be a #define to (*__error()), in which __error() is presumably
 * a function that returns a pointer to the actual value of errno.  Presumably
 * also, __error() provides appropriate serialization guarantees.
 */

#include <assert.h>
#include <errno.h>     /* file stats functions use errno */
#include <libgen.h>    /* basename(), dirname() */
#include <stdio.h>     /* FILE* and other file I/O */
#include <string.h>    /* strlen() */

#include <sys/types.h> /* for file stats */
#include <sys/stat.h>  /* for file stats */


void
split_pathname (char** parentdir, char** namestem, char** extn, const char* const path)
{
  char* path_copy = NULL;
  char* s = NULL;
  int n = 0, k = 0;

  *extn = NULL; /* Will be set to at least the empty string */

  path_copy = bebop_strdup (path);
  *parentdir = bebop_strdup (dirname (path_copy));
  bebop_free (path_copy);

  path_copy = bebop_strdup (path);
  s = bebop_strdup (basename (path_copy));
  bebop_free (path_copy);

  n = strlen (s);
  /* 
   * Extract the extension 
   *
   * FIXME (mfh 23 Nov 2007): not clear what to do for filenames with more
   * than one dot in the extension, for example ".so.1" for different versions 
   * of dynamic libraries.  Are they ".1" files, ".so.1" files, or ".so" files?  
   * Here we read backwards from the end of the string, so we would report this
   * example as a ".1" file type.
   */ 
  for (k = n - 1; k >= 0; k--)
    {
      if (s[k] == '.')
	{
	  s[k] = (char) 0;
	  if (k < n - 1) 
	    {
	      /* Save the extension */
	      *extn = bebop_strdup (&s[k + 1]);
	    }
	  break;
	}
    }
  *namestem = s;
  if (*extn == NULL)
    *extn = bebop_strdup ("");
}


int
directory_p (const char* const path)
{
  int saved_errno = 0;
  int status = 0;
  struct stat s;

  errno = 0;
  status = stat (path, &s);
  saved_errno = errno;

  if (status != 0)
    {
      bebop_log (1, "*** directory_p: stat failed: errno = %d ***\n", saved_errno);
      return 0;
    }
  else
    {
      mode_t mode = s.st_mode;
      if (S_ISDIR (mode))
	return 1; /* is a directory */
      else
	return 0;
    }
}


int
regular_file_p (const char* const path)
{
  int saved_errno = 0;
  int status = 0;
  struct stat s;

  errno = 0;
  status = stat (path, &s);
  saved_errno = errno;

  if (status != 0)
    {
      bebop_log (1, "*** regular_file_p: stat failed: errno = %d ***\n", saved_errno);
      return 0;
    }
  else
    {
      mode_t mode = s.st_mode;
      if (S_ISREG (mode))
	return 1; /* is a regular file */
      else
	return 0;
    }
}

unsigned int
max_linelength_in_file (FILE* file)
{
  unsigned int max_linelen = 0;
  unsigned int cur_count = 0;
  unsigned int prev_count = 0;
  unsigned int line_number = 1;
  long file_stream_pos = -1;

  if (ferror (file))
    bebop_fatal_error ("IO", "invalid file stream!\n");

  /* 
   * Get the current file stream position.  After we're done parsing,
   * we'll return the file stream to that position.
   */
  file_stream_pos = ftell (file);
  if (file_stream_pos == -1)
    bebop_fatal_error ("IO", "file stream is valid but ftell returned -1\n");

  while (! feof (file))
    {
      prev_count = 0;
      cur_count = 0;
      /*
       * POTENTIAL BUG (mfh 28 July 2006): what if your system doesn't
       * represent an endline as '\n' ???
       */
      while (fgetc (file) != '\n' && ! feof (file))
	{
	  prev_count = cur_count;
	  cur_count++;
	  if (prev_count > cur_count)
	    {
	      /* Uh oh, overflow! */
	      bebop_fatal_error ("IO", "Number of characters in line "
				 "%d of given file is longer than the"
				 " largest number representable by an" 
				 " unsigned int! ***\n", line_number);
	    }
	}

      max_linelen = (cur_count > max_linelen) ? cur_count : max_linelen;
      line_number++;
    }

  /* Restore the original file position. */
  if (0 != fseek (file, file_stream_pos, SEEK_SET))
    {
      bebop_fatal_error ("IO", "Failed to seek back to original file "
			 "stream position!\n");
    }

  return max_linelen;
}


