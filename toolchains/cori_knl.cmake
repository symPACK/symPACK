
option(ENABLE_VTUNE "Enable VTUNE" OFF)

set(ENABLE_ARIES ON CACHE BOOL "...")
set(ENABLE_KNL ON CACHE BOOL "...")


set(SCOTCH_DIR $ENV{SWPREFIX}/release/scotch_install)
set(METIS_DIR $ENV{SWPREFIX}/release/metis_install)
set(PARMETIS_DIR $ENV{SWPREFIX}/release/parmetis_install)


#need to set METIS_DIR accordingly
set(ENABLE_METIS ON CACHE BOOL "...")
#need to set PARMETIS_DIR accordingly
set(ENABLE_PARMETIS ON CACHE BOOL "...")
#need to set SCOTCH_DIR accordingly
set(ENABLE_SCOTCH ON CACHE BOOL "...")

set(CMAKE_CXX_COMPILER CC)
#set(MPI_CXX_COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS} -mkl=sequential -axMIC-AVX512,AVX2 -g -dynamic" CACHE STRING "" FORCE)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=sequential -axMIC-AVX512,AVX2 -g -dynamic" CACHE STRING "" FORCE)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mkl=sequential -axMIC-AVX512,AVX2 -g -dynamic" CACHE STRING "" FORCE)
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl=sequential -axMIC-AVX512,AVX2 -g -dynamic" CACHE STRING "" FORCE)
set(CMAKE_C_COMPILER cc)
set(CMAKE_Fortran_COMPILER ftn)
set(MPI_C_COMPILER cc)
set(MPI_CXX_COMPILER CC)


if(ENABLE_VTUNE)
add_compile_options("$<$<COMPILE_LANGUAGE:C>:-g>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-g>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-g>")


add_compile_options("$<$<COMPILE_LANGUAGE:C>:-inline-debug-info>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-inline-debug-info>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-inline-debug-info>")




add_compile_options("$<$<COMPILE_LANGUAGE:C>:-dynamic>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-dynamic>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-dynamic>")

SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -g -dynamic" CACHE STRING "" FORCE )


add_compile_options("$<$<COMPILE_LANGUAGE:C>:-opt-report=5>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-opt-report=5>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-opt-report=5>")

add_compile_options("$<$<COMPILE_LANGUAGE:C>:-opt-report-phase=vec>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-opt-report-phase=vec>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-opt-report-phase=vec>")

add_compile_options("$<$<COMPILE_LANGUAGE:C>:-opt-report-file=stderr>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-opt-report-file=stderr>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-opt-report-file=stderr>")
endif()

add_compile_options("$<$<COMPILE_LANGUAGE:C>:-axMIC-AVX512,AVX2>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-axMIC-AVX512,AVX2>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-axMIC-AVX512,AVX2>")

add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-mkl>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>")

