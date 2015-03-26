# Author: Mark Hoemmen
# Time-stamp: <2008-07-04 14:16:09 mhoemmen>
#
# Makefile for the BeBOP Utility Library test suite.
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


## Object files for unit tests.  Each one has a main() routine.
TESTOBJ = test_extstring.o \
	test_list.o \
	test_random_number.o \
	test_split.o \
	test_split_pathname.o \
	test_max_linelength_in_file.o \
	test_complex.o

TESTEXE = test_extstring.x \
	test_list.x \
	test_random_number.x \
	test_split.x \
	test_split_pathname.x \
	test_max_linelength_in_file.x \
	test_complex.x

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

.PHONY: all 

all: $(TESTEXE)

test_extstring.x: test_extstring.o $(LIBUTIL_STATIC)
	$(LINKER) $(CFLAGS) $(LIBS) $(LDFLAGS) -o $@ $< $(LIBUTIL_STATIC)	 $(LIBS)

test_list.x: test_list.o $(LIBUTIL_STATIC) 
	$(LINKER) $(CFLAGS) $(LIBS) $(LDFLAGS) -o $@ $< $(LIBUTIL_STATIC)	 $(LIBS)

test_random_number.x: test_random_number.o $(LIBUTIL_STATIC)
	$(LINKER) $(CFLAGS) $(LIBS) $(LDFLAGS) -o $@ $< $(LIBUTIL_STATIC)	 $(LIBS)

test_split.x: test_split.o $(LIBUTIL_STATIC)
	$(LINKER) $(CFLAGS) $(LIBS) $(LDFLAGS) -o $@ $< $(LIBUTIL_STATIC)	 $(LIBS)

test_split_pathname.x: test_split_pathname.o $(LIBUTIL_STATIC)
	$(LINKER) $(CFLAGS) $(LIBS) $(LDFLAGS) -o $@ $< $(LIBUTIL_STATIC)	 $(LIBS)

test_max_linelength_in_file.x: test_max_linelength_in_file.o $(LIBUTIL_STATIC)
	$(LINKER) $(CFLAGS) $(LIBS) $(LDFLAGS) -o $@ $< $(LIBUTIL_STATIC)	 $(LIBS)

test_complex.x: test_complex.o $(LIBUTIL_STATIC)
	$(LINKER) $(CFLAGS) $(LIBS) $(LDFLAGS) -o $@ $< $(LIBUTIL_STATIC)	 $(LIBS)

#####################################################################

test_extstring.o: test_extstring.c $(HEADERDIR)/extstring.h \
	$(HEADERDIR)/log.h \
	$(HEADERDIR)/malloc.h \
	$(HEADERDIR)/random_number.h \
	$(HEADERDIR)/util.h 

test_random_number.o: test_random_number.c \
	$(HEADERDIR)/random_number.h $(HEADERDIR)/malloc.h \
	$(HEADERDIR)/util.h

test_split.o: test_split.c $(HEADERDIR)/malloc.h \
	$(HEADERDIR)/util.h $(HEADERDIR)/log.h $(HEADERDIR)/split.h \
	$(HEADERDIR)/extstring.h

test_split_pathname.o: test_split_pathname.c \
	$(HEADERDIR)/file.h \
	$(HEADERDIR)/log.h  \
	$(HEADERDIR)/malloc.h \
	$(HEADERDIR)/util.h

test_list.o: test_list.c $(HEADERDIR)/list.h \
	$(HEADERDIR)/log.h \
	$(HEADERDIR)/util.h

test_max_linelength_in_file.o: test_max_linelength_in_file.c \
	$(HEADERDIR)/file.h \
	$(HEADERDIR)/log.h \
	$(HEADERDIR)/random_number.h \
	$(HEADERDIR)/util.h \
	$(HEADERDIR)/config.h 

test_complex.o: test_complex.c $(HEADERDIR)/complex.h \
	$(HEADERDIR)/util.h $(HEADERDIR)/log.h $(HEADERDIR)/config.h

##########################################################################
## UTILITY STEPS
##########################################################################

ctags: 
	ctags *.[ch]

etags:
	etags *.[ch]

clean:
	$(RM) $(TESTEXE) $(TESTOBJ)
