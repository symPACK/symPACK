set(ENABLE_ARIES ON CACHE BOOL "...")
#need to set METIS_DIR accordingly
set(ENABLE_METIS ON CACHE BOOL "...")
#need to set PARMETIS_DIR accordingly
set(ENABLE_PARMETIS ON CACHE BOOL "...")

set(CMAKE_CXX_COMPILER CC)
set(CMAKE_CXX_COMPILER CC)
set(CMAKE_CXX_COMPILER CC)
set(MPI_CXX_COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS} -gxx-name=g++-4.7 -mkl=sequential" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gxx-name=g++-4.7 -mkl=sequential" CACHE STRING "" FORCE)
set(CMAKE_C_COMPILER cc)
set(CMAKE_Fortran_COMPILER ftn)
set(MPI_C_COMPILER cc)
set(MPI_CXX_COMPILER CC)

add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-gxx-name=g++-4.7>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-mkl=sequential>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>")

