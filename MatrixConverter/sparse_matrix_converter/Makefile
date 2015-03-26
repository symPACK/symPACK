# Author: Mark Hoemmen
# Time-stamp: <2008-07-04 14:18:54 mhoemmen>
#
# Makefile for a sparse matrix utility that converts between various
# file formats for sparse matrices, and also between various internal
# storage formats for sparse matrices.  Harwell-Boeing, MatrixMarket
# and Matlab file formats are supported.
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

##
## Subdirectories in which to run Make recursively.  Yes, yes, I'm aware
## of the paper "Recursive MAKE Considered Harmful."
##
SUBDIRS = src

#########################################################################
## OUTPUT RULES
#########################################################################

.PHONY: all 

all: 
	@for i in $(SUBDIRS); do $(MAKE) -C $$i; done

ctags: 
	@for i in $(SUBDIRS); do $(MAKE) -C $$i ctags; done

etags:
	@for i in $(SUBDIRS); do $(MAKE) -C $$i etags; done

clean:
	@for i in $(SUBDIRS); do $(MAKE) -C $$i clean; done

