CC           = mpicc 
CXX          = mpicxx
FC           = mpif90
LOADER       = mpicxx


PREFIX = ${HOME}/software

#NGCHOL directory
NGCHOL_DIR     = ${HOME}/Work/postdoc-lbl/cholesky
#SPARSELIB_DIR = ${NGCHOL_DIR}/MatrixConverter
#SPARSELIB_INCLUDE = -I${SPARSELIB_DIR}/bebop_util/include -I${SPARSELIB_DIR}/sparse_matrix_converter/include
INCLUDES          = -I${NGCHOL_DIR}/SuperRefine/include 
#${SPARSELIB_INCLUDE}

GFORTRAN_LIB     = -lgfortran
#SPARSELIB_LIB    = ${SPARSELIB_DIR}/sparse_matrix_converter/libsparse_matrix_converter.a  ${SPARSELIB_DIR}/bebop_util/libbebop_util.a


COMPILE_FLAG   = -O0 -w -g
LIBS  = ${GFORTRAN_LIB}





COMPILE_DEF += -fopenmp -std=c++11 -DAdd_ 

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
SRCS_CPP = ETree.cpp SparseMatrixStructure.cpp Ordering.cpp main.cpp 
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

cleanall:
	${RM} -f ${OBJS} ${DEPS} *.d.*



