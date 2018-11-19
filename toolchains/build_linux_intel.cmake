set(BLAS_DIR /opt/intel/mkl)
set(LAPACK_DIR /opt/intel/mkl)
#Required for mkl
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${BLAS_DIR}/lib/intel64)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${LAPACK_DIR}/lib/intel64)

set(CMAKE_PREFIX_PATH ${CMAKE_LIBRARY_PATH} ${BLAS_DIR})
set(CMAKE_PREFIX_PATH ${CMAKE_LIBRARY_PATH} ${LAPACK_DIR})

set(SCOTCH_DIR /home/mjacquel/software/gnu/release/scotch_6.0.4_install)
set(ENABLE_SCOTCH ON CACHE BOOL "...")

set(METIS_DIR /home/mjacquel/software/intel/release/parmetis_install)
set(ENABLE_METIS ON CACHE BOOL "...")

set(PARMETIS_DIR /home/mjacquel/software/intel/release/parmetis_install)
set(ENABLE_PARMETIS ON CACHE BOOL "...")

set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_C_COMPILER icc)
set(CMAKE_Fortran_COMPILER ifort )
set(MPI_C_COMPILER mpiicc )
set(MPI_CXX_COMPILER impiicpc)

add_compile_options("-mkl=sequential")
#add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-gxx-name=g++-4.9>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>")


option(ENABLE_VTUNE "Enable VTUNE" OFF)

STRING( TOLOWER "${CMAKE_BUILD_TYPE}" config_type )
if(ENABLE_VTUNE)
if(NOT config_type STREQUAL "debug")
add_compile_options("$<$<COMPILE_LANGUAGE:C>:-g>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-g>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-g>")
add_compile_options("$<$<COMPILE_LANGUAGE:C>:-debug\ inline-debug-info>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-debug\ inline-debug-info>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-debug\ inline-debug-info>")
#SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -g -dynamic" CACHE STRING "" FORCE )
#else()
  #SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -dynamic" CACHE STRING "" FORCE )
endif()

#add_compile_options("$<$<COMPILE_LANGUAGE:C>:-dynamic>")
#add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-dynamic>")
#add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-dynamic>")

#add_compile_options("$<$<COMPILE_LANGUAGE:C>:-opt-report=5>")
#add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-opt-report=5>")
#add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-opt-report=5>")
#
#add_compile_options("$<$<COMPILE_LANGUAGE:C>:-opt-report-phase=vec>")
#add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-opt-report-phase=vec>")
#add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-opt-report-phase=vec>")
#
#add_compile_options("$<$<COMPILE_LANGUAGE:C>:-opt-report-file=stderr>")
#add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-opt-report-file=stderr>")
#add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-opt-report-file=stderr>")
endif()


