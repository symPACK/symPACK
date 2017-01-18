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
#set(MPI_CXX_COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS} -mkl=sequential -axMIC-AVX-512,AVX" CACHE STRING "" FORCE)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=sequential -axMIC-AVX-512,AVX" CACHE STRING "" FORCE)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mkl=sequential -axMIC-AVX-512,AVX" CACHE STRING "" FORCE)
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl=sequential -axMIC-AVX-512,AVX" CACHE STRING "" FORCE)
set(CMAKE_C_COMPILER cc)
set(CMAKE_Fortran_COMPILER ftn)
set(MPI_C_COMPILER cc)
set(MPI_CXX_COMPILER CC)

add_compile_options("$<$<COMPILE_LANGUAGE:C>:-axMIC-AVX-512,AVX>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-axMIC-AVX-512,AVX>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-axMIC-AVX-512,AVX>")

add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-mkl=sequential>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>")

