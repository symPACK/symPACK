include ../make.inc

SRCS_C   = 
SRCS_CPP = ETree.cpp utility.cpp timer.cpp blas.cpp lapack.cpp global.cpp LogFile.cpp FBMatrix.cpp 
SRCS_F   =

OBJS     = ${SRCS_CPP:.cpp=.o} ${SRCS_C:.c=.o} ${SRCS_F:.f=.o}
DEPS     = ${SRCS_CPP:.cpp=.d} ${SRCS_C:.c=.d} ${SRCS_F:.f=.d}

# Compilation replacement rules

%.o: %.c
	${CC} -c ${CFLAGS} ${CCDEFS} $< 
%.o: %.cpp
	${CXX} -c $< ${CXXFLAGS} ${CPPDEFS} 
%.o: %.f
	${FC} -c ${FFLAGS} $<
%.o: %.F
	${FC} -c ${FFLAGS} $<


all: libcholesky.a 

-include ${DEPS}

libcholesky.a: ${OBJS}
	$(AR) $(ARFLAGS) $(CHOLESKY_LIB) ${OBJS}
	$(RANLIB) $(CHOLESKY_LIB)

cleanall:
	${RM} -f ${OBJS} ${DEPS} $(CHOLESKY_LIB) *.d.*


