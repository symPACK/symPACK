# Author: Mark Hoemmen
# Created: 04 Feb 2008
# Last modified: 05 Feb 2008
#
# Filters each line of stdin and does substitutions on it, in order to
# generate a sparse matrix code.
#
# template-parser (BeBOP C "templates" parser and code generator)
# Copyright (C) 2008 Mark Hoemmen <mhoemmen@cs.berkeley.edu>
#
# This file is part of the template-parser library, which parses C
# "templates" and generates code for different data types.
# template-parser is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# template-parser is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the source code of this library.  If not, see
# <http://www.gnu.org/licenses/>.
#
######################################################################

import string
import sys

def do_subst(dict, body_template):
    """Do keyword substitution on body_template, using dict.

    Arguments:
    dict -- Python dictionary mapping from variable names (without the
      dollar sign or braces) to their string values.
    body_template -- String template, possibly containing variable
      names on which to do substitution.  These follow the conventions
      of string.Template.
      
    """
    return string.Template(body_template).substitute(dict)


def make_dict(mattype, indextype, indextype_string, datatype, datatype_string):
    """Return a Python dictionary containing the keywords we use."""
    return {'mt' : mattype,
            'it' : indextype,
            'itstring' : indextype_string,
            'dt' : datatype,
            'dtstring' : datatype_string}


dict = make_dict('CSR', 'size_t', 'sizet', 'double', 'r64')
for line in sys.stdin:
	print do_subst(dict, line.rstrip())
