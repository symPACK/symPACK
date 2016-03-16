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
all: run_dense_fanboth_upcxx run_dense_fanboth_mpi run_dense_scalapack run_sparse_push unit_test 

#run_sparse_fanboth_upcxx 

#run_dense_fanboth


run_cholesky: run_cholesky.o ${SYMPACK_LIB}
	($(LOADER) -o $@ run_cholesky.o  $(LOADOPTS)  )


run_cholesky.o: run_cholesky.cpp 
	${CXX} -g -c ${CXXFLAGS} ${CPPDEFS} $< 

run_dense_fanboth: run_dense_fanboth.o ${SYMPACK_LIB}
	($(LOADER) -o $@ run_dense_fanboth.o  $(LOADOPTS) )

run_dense_fanboth.o: run_dense_fanboth.cpp 
	${CXX} -g -c ${CXXFLAGS} ${CPPDEFS} $< 

run_dense_fanboth_upcxx: run_dense_fanboth_upcxx.o ${SYMPACK_LIB}
	($(LOADER) $(LDFLAGS)  -o $@ run_dense_fanboth_upcxx.o  $(LOADOPTS)   )

run_dense_fanboth_upcxx.o: run_dense_fanboth_upcxx.cpp 
	${CXX} -g  -c -DUPCXX ${CXXFLAGS} ${CPPDEFS} $< 

run_dense_fanboth_mpi: run_dense_fanboth_mpi.o ${SYMPACK_LIB}
	($(LOADER) $(LDFLAGS)  -o $@ run_dense_fanboth_mpi.o  $(LOADOPTS)  ) 

run_dense_fanboth_mpi.o: run_dense_fanboth_mpi.cpp 
	${CXX} -g  -o $@ -c ${CXXFLAGS} ${CPPDEFS} $< 


run_sparse_push: run_sparse_push.o ${SYMPACK_LIB}
	($(LOADER) $(LDFLAGS)  $(UPCXX_LDFLAGS) -o $@${SUFFIX} run_sparse_push.o  $(LOADOPTS)  $(UPCXX_LDLIBS) )

run_sparse_push.o: run_sparse_push.cpp 
	${CXX} -g  -c ${UPCXX_CXXFLAGS} ${CXXFLAGS} ${CPPDEFS} $< 

run_sparse_pull: run_sparse_pull.o ${SYMPACK_LIB}
	($(LOADER) $(LDFLAGS)  $(UPCXX_LDFLAGS) -o $@${SUFFIX} run_sparse_pull.o  $(LOADOPTS)  $(UPCXX_LDLIBS) )

run_sparse_pull.o: run_sparse_pull.cpp 
	${CXX} -c ${UPCXX_CXXFLAGS} ${CXXFLAGS} ${CPPDEFS} $< 


test_ordering: test_ordering.o ${SYMPACK_LIB}
	($(LOADER) $(LDFLAGS)  $(UPCXX_LDFLAGS) -o $@${SUFFIX} test_ordering.o  $(LOADOPTS)  $(UPCXX_LDLIBS) )

test_ordering.o: test_ordering.cpp 
	${CXX} -c ${UPCXX_CXXFLAGS} ${CXXFLAGS} -g  ${CPPDEFS} $< 





run_sparse_pull_static: run_sparse_pull_static.o ${SYMPACK_LIB}
	($(LOADER) $(LDFLAGS)  $(UPCXX_LDFLAGS) -o $@${SUFFIX} run_sparse_pull_static.o  $(LOADOPTS)  $(UPCXX_LDLIBS) )

run_sparse_pull_static.o: run_sparse_pull_static.cpp 
	${CXX} -g  -c ${UPCXX_CXXFLAGS} ${CXXFLAGS} ${CPPDEFS} $< 


test_dist_ws: test_dist_ws.o ${SYMPACK_LIB}
	($(LOADER) $(LDFLAGS)  $(UPCXX_LDFLAGS) -o $@${SUFFIX} test_dist_ws.o  $(LOADOPTS)  $(UPCXX_LDLIBS) )

test_dist_ws.o: test_dist_ws.cpp 
	${CXX} -g  -c ${UPCXX_CXXFLAGS} ${CXXFLAGS} ${CPPDEFS} $< 



test_split: test_split.o ${SYMPACK_LIB}
	($(LOADER) $(UPCXX_LDFLAGS)  $(LDFLAGS) -o $@ test_split.o  $(LOADOPTS)  $(UPCXX_LDLIBS) )

test_split.o: test_split.cpp 
	${CXX} -g  -c ${UPCXX_CXXFLAGS} ${CXXFLAGS} ${CPPDEFS} $< 

run_pull_test: run_pull_test.o ${SYMPACK_LIB}
	($(LOADER) $(UPCXX_LDFLAGS)  $(LDFLAGS) -o $@ run_pull_test.o  $(LOADOPTS)  $(UPCXX_LDLIBS) )

run_pull_test.o: run_pull_test.cpp 
	${CXX} -g  -c ${UPCXX_CXXFLAGS} ${CXXFLAGS} ${CPPDEFS} $< 

run_isendrecv_test: run_isendrecv_test.o ${SYMPACK_LIB}
	($(LOADER) $(UPCXX_LDFLAGS)  $(LDFLAGS) -o $@ run_isendrecv_test.o  $(LOADOPTS)  $(UPCXX_LDLIBS) )

run_isendrecv_test.o: run_isendrecv_test.cpp 
	${CXX} -g  -c ${UPCXX_CXXFLAGS} ${CXXFLAGS} ${CPPDEFS} $< 





unit_test: unit_test.o ${SYMPACK_LIB}
	($(LOADER) $(LDFLAGS)  -o $@ unit_test.o  $(LOADOPTS)    )

unit_test.o: unit_test.cpp 
	${CXX} -g  -c ${CXXFLAGS} ${CPPDEFS} $< 



run_dense_scalapack: run_dense_scalapack.o ${SYMPACK_LIB}
	($(LOADER) -o $@ run_dense_scalapack.o  $(LOADOPTS)  )

run_dense_scalapack.o: run_dense_scalapack.cpp 
	${CXX} -g  -c ${CXXFLAGS} ${CPPDEFS} $< 






run_mumps: run_mumps.o ${SYMPACK_LIB}
	($(LOADER) -o $@${SUFFIX} run_mumps.o  $(DMUMPS_LIB) $(ZMUMPS_LIB) $(MUMPS_ORDERINGS) ${MUMPSLOADOPTS} )

run_mumps.o: run_mumps.cpp
	${CXX} -c $(MUMPS_INCLUDE)  ${CXXFLAGS} ${CPPDEFS} -O3 $<


run_superlu: run_superlu3.o ${SYMPACK_LIB}
	($(LOADER) $(LDFLAGS) -o $@ run_superlu3.o  $(SUPERLU_LIB) $(PEXSI_LIB) $(PARMETIS_LIB) ${LOADOPTS} )

run_superlu3.o: run_superlu3.cpp
	${CXX} -c $(SUPERLU_INCLUDE)  ${CXXFLAGS} ${CPPDEFS} $<

run_pastix: run_pastix.o ${SYMPACK_LIB}
	($(LOADER) $(LDFLAGS) -o $@${SUFFIX} run_pastix.o  $(PASTIX_LIB) $(PEXSI_LIB) $(PARMETIS_LIB) ${LOADOPTS} )

run_pastix.o: run_pastix.cpp
	${CXX} -c $(PASTIX_INCLUDE)  ${CXXFLAGS} ${CPPDEFS} $< 


-include ${DEPS}

${SYMPACK_LIB}:
	(cd ${PEXSI_DIR}/src; make all)

cleanall:
	rm -f ${OBJS} ${DEPS} *.d.* *.o


