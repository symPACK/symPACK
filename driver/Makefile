# NOTE: This Makefile does NOT support auto-dependency for the .h files.
# If the header files are changed, do "make clean" first.



include ../make.inc


SRCS  = 
#run_cholesky.cpp run_dense_fanboth.cpp 

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
all: run_dense_fanboth_upcxx run_dense_fanboth_mpi run_dense_scalapack run_sparse_fanboth_upcxx unit_test 

#run_sparse_fanboth_upcxx 

#run_dense_fanboth


run_cholesky: run_cholesky.o ${NGCHOL_LIB}
	($(LOADER) -o $@ run_cholesky.o  $(LOADOPTS)  )


run_cholesky.o: run_cholesky.cpp 
	${CXX} -g -c ${CXXFLAGS} ${CPPDEFS} $< 

run_dense_fanboth: run_dense_fanboth.o ${NGCHOL_LIB}
	($(LOADER) -o $@ run_dense_fanboth.o  $(LOADOPTS) )

run_dense_fanboth.o: run_dense_fanboth.cpp 
	${CXX} -g -c ${CXXFLAGS} ${CPPDEFS} $< 

run_dense_fanboth_upcxx: run_dense_fanboth_upcxx.o ${NGCHOL_LIB}
	($(LOADER) $(LDFLAGS)  -o $@ run_dense_fanboth_upcxx.o  $(LOADOPTS)   )

run_dense_fanboth_upcxx.o: run_dense_fanboth_upcxx.cpp 
	${CXX} -g  -c -DUPCXX ${CXXFLAGS} ${CPPDEFS} $< 

run_dense_fanboth_mpi: run_dense_fanboth_mpi.o ${NGCHOL_LIB}
	($(LOADER) $(LDFLAGS)  -o $@ run_dense_fanboth_mpi.o  $(LOADOPTS)  ) 

run_dense_fanboth_mpi.o: run_dense_fanboth_mpi.cpp 
	${CXX} -g  -o $@ -c ${CXXFLAGS} ${CPPDEFS} $< 


run_sparse_fanboth_upcxx: run_sparse_fanboth_upcxx.o ${NGCHOL_LIB}
	($(LOADER) $(LDFLAGS)  -o $@ run_sparse_fanboth_upcxx.o  $(LOADOPTS)   )

run_sparse_fanboth_upcxx.o: run_sparse_fanboth_upcxx.cpp 
	${CXX} -g  -c ${CXXFLAGS} ${CPPDEFS} $< 

unit_test: unit_test.o ${NGCHOL_LIB}
	($(LOADER) $(LDFLAGS)  -o $@ unit_test.o  $(LOADOPTS)    )

unit_test.o: unit_test.cpp 
	${CXX} -g  -c ${CXXFLAGS} ${CPPDEFS} $< 



run_dense_scalapack: run_dense_scalapack.o ${NGCHOL_LIB}
	($(LOADER) -o $@ run_dense_scalapack.o  $(LOADOPTS)  )

run_dense_scalapck.o: run_dense_scalapack.cpp 
	${CXX} -g  -c ${CXXFLAGS} ${CPPDEFS} $< 


run_mumps: run_mumps.o ${NGCHOL_LIB}
	($(LOADER) -o $@ run_mumps.o  $(DMUMPS_LIB) $(ZMUMPS_LIB) $(MUMPS_ORDERINGS) ${LOADOPTS} )

run_mumps.o: run_mumps.cpp
	${CXX} -c $(MUMPS_INCLUDE)  ${CXXFLAGS} ${CPPDEFS} $<



-include ${DEPS}

${NGCHOL_LIB}:
	(cd ${PEXSI_DIR}/src; make all)

cleanall:
	rm -f ${EXES} ${OBJS} ${DEPS} *.d.* *.o


