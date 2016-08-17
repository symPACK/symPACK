
# this one is important
#SET(CMAKE_SYSTEM_NAME Linux)
#this one not so much
#SET(CMAKE_SYSTEM_VERSION 1)

set(CMAKE_CXX_COMPILER icpc)
set(MPI_CXX_COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS} -gxx-name=g++-4.9" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gxx-name=g++-4.9" CACHE STRING "" FORCE)
set(CMAKE_C_COMPILER icc)
set(CMAKE_Fortran_COMPILER ifort)
set(MPI_C_COMPILER mpiicc)
set(MPI_CXX_COMPILER mpiicpc)
#add_definitions("-gxx-name=g++-4.9")
#add_definitions("-std=c++11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" CACHE STRING "" FORCE)
#add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-gxx-name=g++-4.9>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>")


