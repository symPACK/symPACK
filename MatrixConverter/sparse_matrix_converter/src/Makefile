# Author: Mark Hoemmen
# Since: 26 May 2005
# Time-stamp: <2008-07-04 14:21:04 mhoemmen>
#
# For a sparse matrix utility that converts between various file formats 
# for sparse matrices, and also between various internal storage formats
# for sparse matrices.  Harwell-Boeing, MatrixMarket and Matlab file
# formats are supported.
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

## Set up library and include paths
CPPFLAGS += -I../../bebop_util/include 
LIBS += -L../../bebop_util -lbebop_util 

#
# Object files for the sparse matrix converter library.
#
LIBSMCOBJ = bcoo_matrix.o bcsr_matrix.o coo_matrix.o coord_elem.o \
	    csc_matrix.o csr_matrix.o iohb.o jad_matrix.o mmio.o \
	    read_mm.o sparse_matrix.o sparse_matrix_ops.o sparse_vector.o \
            csr_matmatmult.o interface.o csr_spmv.o csr_triple_product.o \
            csr_transpose.o csr_expand_to_dense.o csr_weighted_jacobi.o \
            csr_trisolve.o
EXEOBJ = convert2csc.o
 #main.o  

## Filename of static library
LIBSMC_STATIC = ../libsparse_matrix_converter$(STATIC_LIB_EXTENSION)

## Filename of dynamic library
LIBSMC_DYNAMIC = ../libsparse_matrix_converter$(DYNAMIC_LIB_EXTENSION)

## Filename of SMC executable
EXEOUT = ../sparse_matrix_converter

## Location of SMC header files
HEADERDIR = ../include/bebop/smc

#########################################################################
## GENERAL BUILD RULES 
#########################################################################

.c.o: 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

#########################################################################
## OUTPUT RULES
#########################################################################

.PHONY: all dynamic static 

all: $(EXEOUT) dynamic static

dynamic: $(LIBSMC_DYNAMIC)

static: $(LIBSMC_STATIC)

$(EXEOUT): $(LIBSMC_STATIC) $(EXEOBJ)
	$(LINKER) -static $(CPPFLAGS) $(CFLAGS) -o $(EXEOUT) $(EXEOBJ) $(LIBSMC_STATIC) $(LIBS) -lpthread

$(LIBSMC_STATIC): $(LIBSMCOBJ) 
	$(AR) $(ARFLAGS) $@ $(LIBSMCOBJ)

$(LIBSMC_DYNAMIC): $(LIBSMCOBJ)
	$(LINKER) $(DYLIB_CREATE_FLAG) $(CFLAGS) $(LIBSMCOBJ) $(LIBS) $(LDFLAGS) -o $(LIBSMC_DYNAMIC) $(LIBS) -lpthread

smc_prompt: $(LIBSMC_STATIC) prompt.tab.o prompt.lexer.o
	$(LINKER) $(CPPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) -o $@ $< $(LIBSMC_STATIC) -lfl $(LIBS) -lpthread

#########################################################################

prompt.tab.o: prompt.tab.c $(HEADERDIR)/sparse_matrix.h $(HEADERDIR)/sparse_matrix_ops.h 

prompt.lexer.o: prompt.lexer.c

prompt.lexer.c: prompt.lex
	flex prompt.lex
	mv lex.yy.c prompt.lexer.c

prompt.tab.c: prompt.y
	bison -d prompt.y

#interface: interface_wrap.c 
#	$(CC) $(CPPFLAGS) $(CFLAGS) -c interface.c interface_wrap.c -I/usr/apps/python/default/include/python2.3
#	$(CC) -shared $(LIBS) interface.o interface_wrap.o $(LIBSMCOUT) -o interface.so

# FIXME: needs more dependencies (on includes)!
#interface_wrap.c: interface.i $(LIBSMCOUT)
#	swig -python interface.i


#########################################################################
## BUILD RULES 
#########################################################################

convert2csc.o: convert2csc.c $(HEADERDIR)/sparse_matrix.h $(HEADERDIR)/sparse_matrix_ops.h 
main.o: main.c $(HEADERDIR)/sparse_matrix.h $(HEADERDIR)/sparse_matrix_ops.h 

bcoo_matrix.o: bcoo_matrix.c $(HEADERDIR)/bcoo_matrix.h $(HEADERDIR)/bcsr_matrix.h 

bcsr_matrix.o: bcsr_matrix.c $(HEADERDIR)/bcsr_matrix.h 

coo_matrix.o: coo_matrix.c $(HEADERDIR)/coo_matrix.h \
	$(HEADERDIR)/csc_matrix.h $(HEADERDIR)/read_mm.h  

coord_elem.o: coord_elem.c $(HEADERDIR)/coord_elem.h 

csc_matrix.o: csc_matrix.c $(HEADERDIR)/csc_matrix.h \
	$(HEADERDIR)/coo_matrix.h $(HEADERDIR)/iohb.h \
	$(HEADERDIR)/sparse_vector.h 

csr_expand_to_dense.o: csr_expand_to_dense.c $(HEADERDIR)/csr_expand_to_dense.h

csr_trisolve.o: csr_trisolve.c $(HEADERDIR)/csr_trisolve.h

csr_weighted_jacobi.o: csr_weighted_jacobi.c $(HEADERDIR)/csr_weighted_jacobi.h

csr_matrix.o: csr_matrix.c $(HEADERDIR)/csr_matrix.h \
	$(HEADERDIR)/coo_matrix.h $(HEADERDIR)/csc_matrix.h \
	$(HEADERDIR)/iohb.h $(HEADERDIR)/read_mm.h \
	$(HEADERDIR)/sparse_vector.h $(HEADERDIR)/csr_matmatmult.h \
	$(HEADERDIR)/csr_triple_product.h $(HEADERDIR)/csr_expand_to_dense.h \
	$(HEADERDIR)/csr_spmv.h $(HEADERDIR)/csr_trisolve.h

csr_matmatmult.o: csr_matmatmult.c $(HEADERDIR)/csr_matmatmult.h 

csr_spmv.o: csr_spmv.c $(HEADERDIR)/csr_spmv.h

csr_transpose.o: csr_transpose.c $(HEADERDIR)/csr_transpose.h

csr_triple_product.o: csr_triple_product.c $(HEADERDIR)/csr_triple_product.h \
	$(HEADERDIR)/csr_transpose.h

interface.o: interface.c $(HEADERDIR)/interface.h $(HEADERDIR)/sparse_matrix.h \
	$(HEADERDIR)/sparse_matrix_ops.h $(HEADERDIR)/csr_matrix.h

iohb.o: iohb.c $(HEADERDIR)/iohb.h

jad_matrix.o: jad_matrix.c $(HEADERDIR)/jad_matrix.h $(HEADERDIR)/csr_matrix.h 

mmio.o: mmio.c $(HEADERDIR)/mmio.h 

read_mm.o: read_mm.c $(HEADERDIR)/read_mm.h $(HEADERDIR)/coo_matrix.h \
	$(HEADERDIR)/coord_elem.h $(HEADERDIR)/csc_matrix.h \
	$(HEADERDIR)/csr_matrix.h $(HEADERDIR)/mmio.h 

sparse_matrix.o: sparse_matrix.c $(HEADERDIR)/sparse_matrix.h \
	$(HEADERDIR)/bcoo_matrix.h $(HEADERDIR)/bcsr_matrix.h \
	$(HEADERDIR)/coo_matrix.h $(HEADERDIR)/csc_matrix.h \
	$(HEADERDIR)/csr_matrix.h $(HEADERDIR)/jad_matrix.h 

sparse_matrix_ops.o: sparse_matrix_ops.c $(HEADERDIR)/sparse_matrix_ops.h \
	$(HEADERDIR)/bcoo_matrix.h $(HEADERDIR)/bcsr_matrix.h \
	$(HEADERDIR)/coo_matrix.h $(HEADERDIR)/csc_matrix.h \
	$(HEADERDIR)/csr_matrix.h $(HEADERDIR)/jad_matrix.h \
	$(HEADERDIR)/read_mm.h $(HEADERDIR)/sparse_matrix.h

sparse_vector.o: sparse_vector.c $(HEADERDIR)/sparse_vector.h 

##########################################################################
## UTILITY STEPS
##########################################################################

ctags: 
	ctags *.[ch]

etags:
	etags *.[ch]

clean:
	$(RM) $(LIBSMCOBJ) 
	$(RM) $(LIBSMC_STATIC) $(LIBSMC_DYNAMIC)
	$(RM) $(EXEOBJ)
	$(RM) $(EXEOUT)
	$(RM) prompt.tab.c prompt.tab.h prompt.tab.o
	$(RM) prompt.lexer.c prompt.lexer.o
	$(RM) smc_prompt
#	$(RM) interface_wrap.*
#	$(RM) interface.so 
