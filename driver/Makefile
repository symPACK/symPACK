# NOTE: This Makefile does NOT support auto-dependency for the .h files.
# If the header files are changed, do "make clean" first.



include ../make.inc


SRCS  = run_cholesky.cpp run_dense_fanboth.cpp 

OBJS = ${SRCS:.cpp=.o}
DEPS = ${SRCS:.cpp=.d}
EXES = ${SRCS:.cpp=}	

# Compilation replacement rules

%.o: %.c
	${CC} -c ${CFLAGS} ${CCDEFS} $< 
%.o: %.cpp
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 
%.o: %.f
	${FC} -c ${FFLAGS} $<
%.o: %.F
	${FC} -c ${FFLAGS} $<

#all: run_cholesky run_dense_fanboth
all: run_dense_fanboth_upcxx run_dense_fanboth


run_cholesky: run_cholesky.o ${CHOLESKY_LIB}
	($(LOADER) -o $@ run_cholesky.o  $(LOADOPTS) $(IPM) )


run_cholesky.o: run_cholesky.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

run_dense_fanboth: run_dense_fanboth.o ${CHOLESKY_LIB}
	($(LOADER) -o $@ run_dense_fanboth.o  $(LOADOPTS) $(IPM) )

run_dense_fanboth.o: run_dense_fanboth.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

run_dense_fanboth_upcxx: run_dense_fanboth_upcxx.o ${CHOLESKY_LIB}
	($(LOADER)  -DGASNET_SEQ -DUSE_GASNET_FAST_SEGMENT -DONLY_MSPACES  -o $@ run_dense_fanboth_upcxx.o  $(LOADOPTS)   $(IPM) )

run_dense_fanboth_upcxx.o: run_dense_fanboth_upcxx.cpp 
	${CXX} -g  -c ${CXXFLAGS} ${CPPDEFS} $< 


-include ${DEPS}

${CHOLESKY_LIB}:
	(cd ${PEXSI_DIR}/src; make all)

cleanall:
	rm -f ${EXES} ${OBJS} ${DEPS} *.d.* *.o


include $(GASNET_MAKEFILE)
