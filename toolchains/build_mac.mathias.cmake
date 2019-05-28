#set(BLAS_DIR  /usr/local/Cellar/lapack/3.7.0_1)
#set(LAPACK_DIR  /usr/local/Cellar/lapack/3.7.0_1)
set(BLAS_DIR /usr/local/Cellar/openblas/0.3.5)
set(LAPACK_DIR /usr/local/Cellar/openblas/0.3.5)
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${BLAS_DIR})
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${LAPACK_DIR})

#set(SCOTCH_DIR /home/mjacquel/software/gnu/release/scotch_6.0.4)
#set(ENABLE_SCOTCH ON CACHE BOOL "...")

#set(METIS_DIR /home/mjacquel/software/gnu/release/parmetis_install)
set(ENABLE_METIS ON CACHE BOOL "...")

#set(PARMETIS_DIR /home/mjacquel/software/gnu/release/parmetis_install)
set(ENABLE_PARMETIS ON CACHE BOOL "...")

set(CMAKE_CXX_COMPILER g++-8)
set(CMAKE_C_COMPILER gcc-8)
set(CMAKE_Fortran_COMPILER gfortran-8)
set(MPI_C_COMPILER mpicc)
set(MPI_CXX_COMPILER mpic++)

add_compile_options("-Wl,-no_pie")
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -Wl,-no_pie" CACHE STRING "" FORCE )

add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>")
STRING( TOLOWER "${CMAKE_BUILD_TYPE}" config_type )
if(config_type STREQUAL "debug")
#add_compile_options("$<$<COMPILE_LANGUAGE:C>:-g>")
#add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-g>")
#add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-g>")
#set(MPI_CXX_COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS} -g" CACHE STRING "" FORCE)
endif()




