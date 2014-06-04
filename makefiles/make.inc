#/usr/bin/bash
PROFILE   =   1
PROF_MPI   =  1

USE_TAU   =   0
USE_AUTO_TAU   =   0


MODE=release
MODE=debug

SUFFIX        =$(MODE)_${PLATFORM}

ifdef USE_TAU
  ifeq ($(USE_TAU),1)
    SUFFIX=$(MODE)_tau_${PLATFORM}
  endif
endif

  ifeq ($(PROFILE),1)
    SUFFIX=$(MODE)_tau_${PLATFORM}
  endif

  ifeq ($(PROF_MPI),1)
    SUFFIX=$(MODE)_tau_${PLATFORM}
  endif

CHOLESKY_DIR     = ${HOME}/Work/postdoc-lbl/cholesky

UPCXX_INSTALL_DIR     = ${PREFIX}/${MODE}/upcxx_install
#UPCXX_INSTALL_DIR     = ${PREFIX}/release/upcxx_install
## Add the following line in the UPC++ application Makefile
include $(UPCXX_INSTALL_DIR)/include/upcxx.mak



# includes

#UPCXX_INCLUDE = -I${UPCXX_DIR}/include -DGASNET_SEQ -DUSE_GASNET_FAST_SEGMENT -DONLY_MSPACES $(GASNET_CPPFLAGS) $(GASNET_CFLAGS)
#UPCXX_INCLUDE = -I${UPCXX_DIR}/include -DGASNET_SEQ $(GASNET_CPPFLAGS) $(GASNET_CFLAGS)
SPARSELIB_INCLUDE    = -I${HOME}/Work/postdoc-lbl/pexsi/Converter/bebop_util/include -I${HOME}/Work/postdoc-lbl/pexsi/Converter/sparse_matrix_converter/include
CHOLESKY_INCLUDE    = -I${CHOLESKY_DIR}/include 
INCLUDES         = ${CHOLESKY_INCLUDE} ${SPARSELIB_INCLUDE} ${UPCXX_INCLUDE}

# Libraries
GFORTRAN_LIB = -lgfortran
LAPACK_LIB       = -L${HOME}/software/scalapack-2.0.2/build/lib -lscalapack  
BLAS_LIB         = -L${HOME}/software/lapack-3.4.2/lib -llapack -llapacke -ltmglib -lblas 
#BLAS_LIB         = -lopenblas
CHOLESKY_LIB        = ${CHOLESKY_DIR}/src/libcholesky_${SUFFIX}.a
SPARSELIB_LIB    =  ${HOME}/Work/postdoc-lbl/pexsi/Converter/sparse_matrix_converter/libsparse_matrix_converter.a ${HOME}/Work/postdoc-lbl/pexsi/Converter/bebop_util/libbebop_util.a
#UPCXX_LIB        = ${UPCXX_DIR}/lib/libupcxx.a $(GASNET_LDFLAGS) $(GASNET_LIBS)

#UPCXX_FLAGS =  -DGASNET_SEQ -DUSE_GASNET_FAST_SEGMENT -DONLY_MSPACES

LIBS             = ${SPARSELIB_LIB} ${CHOLESKY_LIB} ${UPCXX_LDLIBS} ${LAPACK_LIB} ${BLAS_LIB} ${GFORTRAN_LIB} -lgomp ${IPM} 

CC           = mpicc 
CXX          = mpicxx
FC           = mpif90
LOADER       = mpicxx


AR           = ar 
ARFLAGS      = rvscu
# For System V based machine without ranlib, like Cray and SGI,
# use touch instead.
#RANLIB      = touch
RANLIB       = ranlib

RM           = rm
RMFLAGS      = -f

# Different compiling and linking options.

ifeq ($(MODE), debug)
	COMMONDEFS   = -DRELEASE -DAdd_ 
  CFLAGS       = ${INCLUDES} -O0 -w -g 
  FFLAGS       = ${INCLUDES} -O0 -w -g 
  CXXFLAGS     = $(UPCXX_CXXFLAGS) ${INCLUDES} -O0 -w -g -fpermissive
  LDFLAGS      = $(UPCXX_LDFLAGS)
	CCDEFS       = ${COMMONDEFS}
	CPPDEFS      = ${COMMONDEFS}
  LOADOPTS     = ${LIBS}
  FLOADOPTS    = ${LIBS} -lstdc++

endif

ifeq ($(MODE), release)
	COMMONDEFS   = -DRELEASE -g -DAdd_ 
  CFLAGS       = ${INCLUDES} -O3 -w -g 
  FFLAGS       = ${INCLUDES} -O3 -w -g 
  CXXFLAGS     = $(UPCXX_CXXFLAGS) ${INCLUDES} -O3 -w -g -fpermissive
  LDFLAGS      = $(UPCXX_LDFLAGS)
	CCDEFS       = ${COMMONDEFS}
	CPPDEFS      = ${COMMONDEFS}
  LOADOPTS     = ${LIBS}
  FLOADOPTS    = ${LIBS} -lstdc++
endif
#  COMMONDEFS += -DUPCXX

  #COMMONDEFS += -DMAP=col2D
  #COMMONDEFS += -DMAP=modwrap2Dns
  COMMONDEFS += -DMAP=modwrap2D
  #COMMONDEFS += -DMAP=antichevron2D


  ## THIS ONE HAS TO BE DEBUGGED
  #COMMONDEFS += -DMAP=antidiag2D
  COMMONDEFS += -DUSE_ABORT
  COMMONDEFS += -DASYNC_AGGREGATES
  COMMONDEFS += -DASYNC_UPDATES
  #COMMONDEFS += -DNO_ASYNC_COPY
  #COMMONDEFS += -DASYNC_COPY_FENCE
  #COMMONDEFS += -DDRAW_GRAPH
  COMMONDEFS += -DCALLBACK



#  COMMONDEFS += -DUPDATE_LIST 
#  COMMONDEFS += -D_DEBUG_ 
  COMMONDEFS += -D_CHECK_RESULT_ 
  COMMONDEFS += -D_CHECK_RESULT_SEQ_ 
#  COMMONDEFS += -DPROBE_FIRST 
#  COMMONDEFS += -DNO_INTRA_PROFILE 
#WIP
#  COMMONDEFS += -DRELAXED_SNODE 
  COMMONDEFS += -DSINGLE_BLAS





  ifeq ($(USE_TAU),1)
    COMMONDEFS     += -DUSE_TAU 
  endif

  ifeq ($(PROFILE),1)
    COMMONDEFS += -DPROFILE
  endif

  ifeq ($(PROF_MPI),1)
    COMMONDEFS  += -DPMPI
  endif



	COMMONDEFS   += -fopenmp -std=gnu++11 
  ifeq ($(USE_TAU),1)
    CXXFLAGS     += -DUSE_TAU 
  endif

  ifeq ($(PROFILE),1)
    CXXFLAGS     += -DPROFILE
  endif

  ifeq ($(PROF_MPI),1)
    CXXFLAGS     += -DPMPI
  endif



# Generate auto-dependencies 
%.d: %.c
	@set -e; rm -f $@; \
	$(CC) -M $(CCDEFS) $(CFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;\
	rm -f $@.$$$$

%.d: %.cpp
	@set -e; rm -f $@; \
	$(CXX) -M $(CPPDEFS) $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;\
	rm -f $@.$$$$
