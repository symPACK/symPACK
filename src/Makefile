include ../make.inc

SRCS_C   = 
SRCS_CPP = ETree.cpp utility.cpp timer.cpp blas.cpp lapack.cpp global.cpp LogFile.cpp SparseMatrixStructure.cpp IntervalTree.cpp Ordering.cpp SuperNode.cpp CommPull.cpp 
#UPCXX_CPP = FBMatrix.cpp FBMatrix_mpi.cpp FBMatrix_upcxx.cpp
UPCXX_CPP = 
SRCS_F   = genmmd.f  mmdelm.f  mmdint.f  mmdnum.f  mmdupd.f  ordmmd.f amdbar.f boxnd.f gridnd.f

OBJS     = ${UPCXX_CPP:.cpp=.o} ${SRCS_CPP:.cpp=.o} ${SRCS_C:.c=.o} ${SRCS_F:.f=.o}
DEPS     = ${UPCXX_CPP:.cpp=.d} ${SRCS_CPP:.cpp=.d} ${SRCS_C:.c=.d} ${SRCS_F:.f=.d}

# Compilation replacement rules

#FBMatrix.o: FBMatrix.cpp
#	${CXX} -c $< ${UPCXX_CXXFLAGS} ${CXXFLAGS} ${CPPDEFS} 
#FBMatrix_upcxx.o: FBMatrix_upcxx.cpp
#	${CXX} -c $< ${UPCXX_CXXFLAGS} ${CXXFLAGS} ${CPPDEFS} 
#FBMatrix_mpi.o: FBMatrix_mpi.cpp
#	${CXX} -c $< ${UPCXX_CXXFLAGS} ${CXXFLAGS} ${CPPDEFS} 


%.o: %.c
	${CC} -c ${CFLAGS} ${CCDEFS} $< 
%.o: %.cpp
	${CXX} -c $<  ${UPCXX_CXXFLAGS} ${CXXFLAGS} ${CPPDEFS} 
%.o: %.f
	${FC} -c ${FFLAGS} $<
%.o: %.F
	${FC} -c ${FFLAGS} $<


all: libsympack

-include ${DEPS}

libsympack: ${OBJS}
	$(AR) $(ARFLAGS) $(NGCHOL_LIB) ${OBJS}
	$(RANLIB) $(NGCHOL_LIB)

cleanall:
	${RM} -f ${OBJS} ${DEPS} $(NGCHOL_LIB) *.d.*


