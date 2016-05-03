include ../make.inc

SRCS_C   = 
SRCS_CPP = ETree.cpp utility.cpp timer.cpp blas.cpp lapack.cpp global.cpp LogFile.cpp \
           SparseMatrixStructure.cpp IntervalTree.cpp Ordering.cpp SuperNode.cpp CommPull.cpp\
           DistSparseMatrixGraph.cpp
#MemTrack.cpp 
#UPCXX_CPP = FBMatrix.cpp FBMatrix_mpi.cpp FBMatrix_upcxx.cpp
UPCXX_CPP = 
SRCS_F   =  boxnd.f gridnd.f
SRCS_FF   = genmmd.F mmdelm.F  mmdint.F  mmdnum.F  mmdupd.F  ordmmd.F amdbar.F 

OBJS     = ${UPCXX_CPP:.cpp=.o} ${SRCS_CPP:.cpp=.o} ${SRCS_C:.c=.o} ${SRCS_F:.f=.o} ${SRCS_FF:.F=.o}
DEPS     = ${UPCXX_CPP:.cpp=.d} ${SRCS_CPP:.cpp=.d} ${SRCS_C:.c=.d} ${SRCS_F:.f=.d} ${SRCS_FF:.F=.d}

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
	${FC} -c ${FFLAGS} ${FDEFS} $<
%.o: %.F
	${FC} -c ${FFLAGS} ${FDEFS} $<


all: libsympack

-include ${DEPS}

libsympack: ${OBJS}
	$(AR) $(ARFLAGS) $(SYMPACK_LIB) ${OBJS}
	$(RANLIB) $(SYMPACK_LIB)

cleanall:
	${RM} -f ${OBJS} ${DEPS} $(SYMPACK_LIB) *.d.*


