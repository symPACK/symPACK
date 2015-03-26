# Author: Mark Hoemmen
# Time-stamp: <2008-07-04 14:15:57 mhoemmen>
#
# Makefile for the BeBOP Utility Library source code directory.
#
# Copyright (c) 2008, Regents of the University of California 
# All rights reserved.
# Redistribution and use in source and binary forms, with or
# without modification, are permitted provided that the
# following conditions are met:
# 
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright 
#   notice, this list of conditions and the following disclaimer in 
#   the documentation and/or other materials provided with the 
#   distribution.
#
# * Neither the name of the University of California, Berkeley, nor
#   the names of its contributors may be used to endorse or promote
#   products derived from this software without specific prior
#   written permission.  
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE.
##########################################################################

include ../../bebop_make/options
include ../../bebop_make/Makefile.include
include ../../bebop_make/ecl
include ../../bebop_make/setup

##
## Object files for the BeBOP utility library.  Does not include
## ecl_interface.o or other ECL-related files -- see the "ifeq" below.
##
LIBUTILOBJ = complex.o avlstree.o avltree_intpair.o \
	     convert.o extstring.o file.o \
             get_matches.o get_options.o init.o list.o \
             log.o malloc.o merge_sort.o mt19937ar.o \
             random_number.o sort_joint_arrays.o split.o \
             string.o timer.o util.o

ifeq "$(USE_ECL)" "1"
  LIBUTILOBJ += ecl_interface.o
endif

## Filename of static library
LIBUTIL_STATIC = ../libbebop_util$(STATIC_LIB_EXTENSION)

## Filename of dynamic library
LIBUTIL_DYNAMIC = ../libbebop_util$(DYNAMIC_LIB_EXTENSION)

## Location of bebop_util header files
HEADERDIR=../include/bebop/util


#########################################################################
## GENERAL BUILD RULES
#########################################################################

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@


#########################################################################
## OUTPUT RULES
#########################################################################

.PHONY: all dynamic static 

all: $(LIBUTIL_STATIC) $(LIBUTIL_DYNAMIC) 

dynamic: $(LIBUTIL_DYNAMIC)

static: $(LIBUTIL_STATIC) 

$(LIBUTIL_DYNAMIC): $(LIBUTILOBJ)
	$(LINKER) $(DYLIB_CREATE_FLAG) $(CFLAGS) $(LIBS) $(LDFLAGS) $(LIBUTILOBJ) -o $(LIBUTIL_DYNAMIC)

$(LIBUTIL_STATIC): $(LIBUTILOBJ)
	$(AR) $(ARFLAGS) $@ $(LIBUTILOBJ)


#########################################################################
## DEPENDENCIES
#########################################################################

complex.o: complex.c $(HEADERDIR)/complex.h $(HEADERDIR)/config.h

avlstree.o: avlstree.c $(HEADERDIR)/avlstree.h \
	$(HEADERDIR)/malloc.h $(HEADERDIR)/string.h \
	$(HEADERDIR)/util.h $(HEADERDIR)/config.h

avltree_intpair.o: avltree_intpair.c \
	$(HEADERDIR)/avltree_intpair.h $(HEADERDIR)/malloc.h \
	$(HEADERDIR)/util.h $(HEADERDIR)/config.h

convert.o: $(HEADERDIR)/convert.h convert.c $(HEADERDIR)/config.h

extstring.o: extstring.c $(HEADERDIR)/extstring.h \
	$(HEADERDIR)/init.h \
	$(HEADERDIR)/log.h \
	$(HEADERDIR)/malloc.h \
	$(HEADERDIR)/util.h

file.o: file.c $(HEADERDIR)/file.h \
	$(HEADERDIR)/log.h \
	$(HEADERDIR)/malloc.h \
	$(HEADERDIR)/string.h \
	$(HEADERDIR)/util.h

get_matches.o: get_matches.c $(HEADERDIR)/get_matches.h \
	$(HEADERDIR)/list.h

get_options.o: get_options.c $(HEADERDIR)/get_options.h \
	$(HEADERDIR)/malloc.h $(HEADERDIR)/util.h \
	$(HEADERDIR)/log.h $(HEADERDIR)/config.h

init.o: init.c $(HEADERDIR)/init.h \
	$(HEADERDIR)/log.h \
        $(HEADERDIR)/ecl_interface.h \
	$(HEADERDIR)/config.h 

list.o: list.c $(HEADERDIR)/list.h

log.o: log.c $(HEADERDIR)/log.h \
	$(HEADERDIR)/util.h

merge_sort.o: merge_sort.c $(HEADERDIR)/merge_sort.h \
	$(HEADERDIR)/malloc.h \
	$(HEADERDIR)/util.h \
	$(HEADERDIR)/config.h

malloc.o: malloc.c $(HEADERDIR)/malloc.h \
	$(HEADERDIR)/util.h $(HEADERDIR)/config.h

mt19937ar.o: mt19937ar.c $(HEADERDIR)/config.h

random_number.o: random_number.c $(HEADERDIR)/random_number.h \
	$(HEADERDIR)/malloc.h $(HEADERDIR)/util.h \
	$(HEADERDIR)/log.h $(HEADERDIR)/config.h

sort_joint_arrays.o: sort_joint_arrays.c $(HEADERDIR)/sort_joint_arrays.h \
	$(HEADERDIR)/log.h \
	$(HEADERDIR)/malloc.h \
	$(HEADERDIR)/random_number.h \
	$(HEADERDIR)/util.h \
	$(HEADERDIR)/config.h

split.o: split.c $(HEADERDIR)/split.h \
	$(HEADERDIR)/malloc.h \
	$(HEADERDIR)/log.h \
	$(HEADERDIR)/string.h \
	$(HEADERDIR)/util.h \
	$(HEADERDIR)/config.h

string.o: string.c $(HEADERDIR)/string.h \
	$(HEADERDIR)/malloc.h \
	$(HEADERDIR)/config.h

timer.o: timer.c $(HEADERDIR)/timer.h $(HEADERDIR)/config.h

util.o: util.c $(HEADERDIR)/util.h \
	$(HEADERDIR)/random_number.h \
	$(HEADERDIR)/config.h


##########################################################################
## UTILITY STEPS
##########################################################################

ctags: 
	ctags *.[ch]

etags:
	etags *.[ch]

clean:
	$(RM) $(LIBUTILOBJ)
	$(RM) $(LIBUTIL_STATIC)
	$(RM) $(LIBUTIL_DYNAMIC)
