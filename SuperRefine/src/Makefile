CC           = mpiicc 
CXX          = mpiicpc -gxx-name=g++-4.9
FC           = mpiifort
LOADER       = mpiicpc -gxx-name=g++-4.9
#CC           = mpicc 
#CXX          = mpic++
#FC           = mpif90
#LOADER       = mpic++


SWPREFIX = ${HOME}/software

#NGCHOL directory
#NGCHOL_DIR     = ${HOME}/Work/postdoc-lbl/ngchol
NGCHOL_DIR     = ${HOME}/Work/postdoc-lbl/cholesky
METIS_DIR     = $(SWPREFIX)/intel/release/metis-5.1.0/build/Linux-x86_64
PARMETIS_DIR  = $(SWPREFIX)/intel/release/parmetis-4.0.3/build/Linux-x86_64
SCOTCH_DIR    = $(SWPREFIX)/intel/release/scotch_6.0.0_esmumps
#SCOTCH_DIR    = $(PREFIX)/release/scotch_6.0.4
#SCOTCH_DIR    = $(PREFIX)/release/scotch_5.1.12_esmumps

#SPARSELIB_DIR = ${NGCHOL_DIR}/MatrixConverter
#SPARSELIB_INCLUDE = -I${SPARSELIB_DIR}/bebop_util/include -I${SPARSELIB_DIR}/sparse_matrix_converter/include
METIS_INCLUDE        = -I${METIS_DIR}/../../include
SCOTCH_INCLUDE        = -I${SCOTCH_DIR}/include
INCLUDES          = -I${NGCHOL_DIR}/SuperRefine/include ${METIS_INCLUDE} ${SCOTCH_INCLUDE} 
#${SPARSELIB_INCLUDE}

#GFORTRAN_LIB     = -lgfortran
#SPARSELIB_LIB    = ${SPARSELIB_DIR}/sparse_matrix_converter/libsparse_matrix_converter.a  ${SPARSELIB_DIR}/bebop_util/libbebop_util.a

SCOTCH_LIB       = -L${SCOTCH_DIR}/lib -lscotch -lscotcherr
METIS_LIB        = -L${METIS_DIR}/libmetis -lmetis
PARMETIS_LIB     = -L${PARMETIS_DIR}/libparmetis -lparmetis 

COMPILE_FLAG   = -O3 -w -g
LIBS  += ${SCOTCH_LIB} ${GFORTRAN_LIB}
LIBS  += ${METIS_LIB} ${GFORTRAN_LIB}





COMPILE_DEF += -fopenmp -std=c++11 -DAdd_ -DPROFILE 
COMPILE_DEF += -DUSE_METIS
COMPILE_DEF += -DUSE_SCOTCH
COMPILE_DEF += -DFUNDAMENTAL
#COMPILE_DEF += -DRELAXED_SNODE

CFLAGS       = ${COMPILE_FLAG} ${INCLUDES}
FFLAGS       = ${COMPILE_FLAG} ${INCLUDES}
CXXFLAGS     = ${COMPILE_FLAG} ${INCLUDES} 
CCDEFS       = ${COMPILE_DEF} 
CPPDEFS      = ${COMPILE_DEF} 
LOADOPTS     = ${PROFILE_FLAG} ${LIBS} -fopenmp


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

SRCS_C   = iohb.c 
SRCS_CPP = ETree.cpp SparseMatrixStructure.cpp Ordering.cpp main.cpp timer.cpp LogFile.cpp symbol_reordering.cpp  
SRCS_F   = genmmd.f  mmdelm.f  mmdint.f  mmdnum.f  mmdupd.f  ordmmd.f amdbar.f

OBJS     = ${UPCXX_CPP:.cpp=.o} ${SRCS_CPP:.cpp=.o} ${SRCS_C:.c=.o} ${SRCS_F:.f=.o}
DEPS     = ${UPCXX_CPP:.cpp=.d} ${SRCS_CPP:.cpp=.d} ${SRCS_C:.c=.d} ${SRCS_F:.f=.d}

# Compilation replacement rules
%.o: %.c
	${CC} -c ${CFLAGS} ${CCDEFS} $< 
%.o: %.cpp
	${CXX} -c $< ${CXXFLAGS} ${CPPDEFS} 
%.o: %.f
	${FC} -c ${FFLAGS} $<
%.o: %.F
	${FC} -c ${FFLAGS} $<


all: main

-include ${DEPS}

main: ${OBJS}
	($(LOADER) -o $@ ${OBJS} $(LOADOPTS)  )

clean:
	${RM} -f ${OBJS} ${DEPS} *.d.*



