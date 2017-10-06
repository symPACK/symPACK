set( ENV{UPCXX_DIR} "/usr/common/ftg/upcxx/default/PrgEnv-intel--17.0.2.174")
set( ENV{GASNET_LIBRARIES} "-L/opt/cray/pe/pmi/default/lib64 -L/opt/cray/ugni/default/lib64 -L/opt/cray/udreg/default/lib64 -lpmi -lugni -ludreg -lxpmem -lhugetlbfs -lm")
set( ENV{GASNET_DIR} "$ENV{UPCXX_DIR}/gasnet" )
set( ENV{GASNET_CONDUIT} "$ENV{GASNET_DIR}/include/aries-conduit/aries-seq.mak")
set( ENV{GASNET_DEFINES} "-DGASNET_CONDUIT_ARIES -D_GNU_SOURCE=1 -DGASNET_SEQ -DGASNET_ALLOW_OPTIMIZED_DEBUG")

set(ENABLE_ARIES ON CACHE BOOL "...")
#need to set METIS_DIR accordingly
#set(ENABLE_METIS ON CACHE BOOL "...")
#need to set PARMETIS_DIR accordingly
#set(ENABLE_PARMETIS ON CACHE BOOL "...")

#set(ENABLE_MKL ON CACHE BOOL "...")

set(CMAKE_CXX_COMPILER CC)
set(CMAKE_C_COMPILER cc)
set(CMAKE_Fortran_COMPILER ftn)
set(MPI_C_COMPILER cc)
set(MPI_CXX_COMPILER CC)
set(MPI_CXX_COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS} -mkl=sequential")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=sequential")

add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-mkl=sequential>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>")

set(SCOTCH_DIR $ENV{SWPREFIX}/release/scotch_install)
set(ENABLE_SCOTCH ON CACHE BOOL "...")
