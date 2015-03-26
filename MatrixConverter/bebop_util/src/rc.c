/**
 * @file rc.c
 * @author Mark Hoemmen
 * @date Time-stamp: <2008-07-16 10:11:59 mhoemmen>
 * @since 28 July 2006
 *
 * A very simple initialization file (rc file, in Unix jargon) reader.
 * The RC file should have the following format:
 *
 * <line> ::== <comment>|<assignment>
 * <comment> ::== \s*#.*
 * <assignment> ::== <varname>\s*=\s*<value>\s*
 * <varname> ::== [not-whitespace][anything-not-an-equals-sign]*
 * <varname> ::== <number>|<string>
 * <string> ::== ".*"
 * <number> ::== <complex>|<primitive-number>
 * <primitive-number> ::== <floating-point>|<integer>
 * <floating-point> is the usual floating-point number pattern
 * <complex> ::== <primitive-number>\s*\+\s*<primitive-number>{i|I|j|J}
 *                (Note: complex numbers are always stored as pairs 
 *                 of floating-point numbers)
 * <integer> is the usual integer number pattern; any number not a 
 *           floating-point number
 *
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

/**
 * Regular expression for a floating-point number with optional
 * exponent.
 */
#define FP_REGEX_STRING "[-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)?";

/**
 * Regular expression for a single-line string in which the quote
 * character can appear if it is escaped by a backslash. Note the
 * double backslashing.
 */
#define STR_REGEX_STRING "\"[^\"\\\\\\r\\n]*(\\\\.[^\"\\\\\\r\\n]*)*\"";

/**
 * Regular expression for a complex (floating-point) number.
 */
#define COMPLEX_REGEX_STRING  \
  "\\(" FP_REGEX_STRING "\\)" "\\s*[-+]\\s*" "\\(" FP_REGEX_STRING \
  "\\)" "\\s*[iIjJ]"

#define COMMENT_REGEX_STRING "^\s*#"


typedef struct {
  regex_t fp_regex;
  regex_t complex_regex;
  regex_t str_regex;
  regex_t comment_regex;
} regexes_t;



#define BUILD_REGEX( regex, str, what ) \
  do {									\
    errcode = regcomp (&(regex), (str), REG_EXTENDED | REG_NEWLINE);	\
    if (errcode != 0)							\
      {									\
	/* First figure out how long the error message is.  Then */	\
	/* allocate enough space for it and get the error message */	\
	/* in full. */							\
	char* errbuf = NULL;						\
	size_t size = regerror (errcode, &(regex), errbuf, 0);		\
	if (size > 0)							\
	  {								\
	    errbuf = bebop_malloc ((size+1) * sizeof (char));		\
	    size = regerror (errcode, &(regex), errbuf, size);		\
	  }								\
	bebop_log (0, "*** build_regexes: failed to compile regular " \ 
"expression for %s: error message: %s ***\n",		      \
  what, errbuf);					      \
bebop_free (errbuf);					      \
}							      \
} while(0)

    
void
construct_regexes (regexes_t* regexes)
{
  int errcode = 0;

  BUILD_REGEX( (regexes->str_regex), STR_REGEX_STRING, 
	       "quoted strings" );
  BUILD_REGEX( (regexes->fp_regex), FP_REGEX_STRING, 
	       "floating-point numbers" );
  BUILD_REGEX( (regexes->complex_regex), COMPLEX_REGEX_STRING, 
	       "complex floating-point numbers" );
  BUILD_REGEX( (regexes->comment_regex), COMMENT_REGEX_STRING, 
	       "comments" );
}

void
destroy_regexes (regexes_t* regexes)
{
  regfree (&(regexes->fp_regex));
  regfree (&(regexes->str_regex));
  regfree (&(regexes->complex_regex));
  regfree (&(regexes->comment_regex));
}



typedef 
struct 
{ 
  char* varname;
  char* strvalue;
} assignment_t;

assignment_t* 
construct_assignment (char* varname, char* strvalue)
{
  assignment_t* a = bebop_calloc (1, sizeof (assignment_t));
  a->varname = varname;
  a->strvalue = strvalue;
  return a;
}

void
destroy_assignment (assignment_t* a)
{
  if (a != NULL)
    {
      if (a->varname != NULL)
	free (a->varname);
      if (a->strvalue != NULL)
	free (a->strvalue);

      /* Assign NULL to them, even though we are freeing a,
         to make it more likely that we will catch bugs. */
      a->varname = NULL;
      a->strvalue = NULL;
      free (a);
    }
}



typedef enum { STRING_VALUE, COMPLEX_VALUE, FP_VALUE, INTEGER_VALUE, UNKNOWN_VALUE } valtype_t;
  
valtype_t 
what_valtype (assignment_t* a)
{
  /* Assume that whitespace was stripped from beginning and end */

  /* First look for starting and ending quotes */

  /* If quotes not found, check if matches complex number pattern */

  /* If doesn't match complex number pattern, check if matches 
     floating-point pattern */

  /* If doesn't match floating-point pattern, check if matches 
     integer pattern */

  /* If doesn't match integer pattern, report unknown */
  return UNKNOWN_VALUE;
}


/**
 * Returns a list of assignments
 */
list_t 
parse_rc_file (FILE* file)
{
  char* line = NULL;
  /* 
   * First go through and see how long the longest line is in the file.
   * We do this by counting the number of non-newlines in between the 
   * newline characters.  After counting, the function resets the file
   * pointer to its original position (where it was when this function
   * was called).
   */
  const unsigned int maxlinelen = max_linelength_in_file (file);
  list_t L = list_create ();
  regexes_t regexes;

  if (maxlinelen < 2)
    {
      bebop_log (0, "*** parse_rc_file: the longest line in the file "
		"contains less than %d characters, meaning that no line "
		"in the file has the correct format ***\n", 2);
      return L;
    }
 
  /* Allocate a "line buffer" */
  line = bebop_malloc ((maxlinelen + 2) * sizeof (char));

  /* Construct the needed regular expressions */
  construct_regexes (regexes_t* regexes);

  while (! feof (file))
    {
      /* fgets(s,n,f) reads n-1 non-newline characters. */
      fgets (line, maxlinelen+1, file);

      /* Check if line is a whitespace or comment line; if not, parse it */
      if (! whitespace_p (line) && ! comment_p (line))
	{
	  char* varname = NULL;
	  char* strvalue = NULL;

	  if (0 == extract_assignment (&varname, &strvalue, line))
	    L = list_append_item (L, construct_assignment (varname, strvalue));
	  else
	    {
	      bebop_log (0, "*** parse_rc_file: the following line is "
			"not a comment, not whitespace and not an assi"
			"gnment: ***\n%s\n*** We choose to ignore the "
			"line above and continue parsing. ***\n", line);
	    }
	  

	}
    }
  

  
  /* Parse it line-by-line */

  bebop_free (line);
  return L;
}


