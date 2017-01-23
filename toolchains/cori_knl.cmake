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


add_compile_options("$<$<COMPILE_LANGUAGE:C>:-g>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-g>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-g>")


add_compile_options("$<$<COMPILE_LANGUAGE:C>:-debug inline-debug-info>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-debug inline-debug-info>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-debug inline-debug-info>")




add_compile_options("$<$<COMPILE_LANGUAGE:C>:-dynamic>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-dynamic>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-dynamic>")

add_compile_options("$<$<COMPILE_LANGUAGE:C>:-axMIC-AVX512,AVX2>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-axMIC-AVX512,AVX2>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-axMIC-AVX512,AVX2>")

add_compile_options("$<$<COMPILE_LANGUAGE:C>:-opt-report=5>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-opt-report=5>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-opt-report=5>")

add_compile_options("$<$<COMPILE_LANGUAGE:C>:-opt-report-phase=vec>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-opt-report-phase=vec>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-opt-report-phase=vec>")

#add_compile_options("$<$<COMPILE_LANGUAGE:C>:-opt-report-file=$ENV{HOME}/vecreport>")
#add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-opt-report-file=$ENV{HOME}/vecreport>")
#add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-opt-report-file=$ENV{HOME}/vecreport>")
add_compile_options("$<$<COMPILE_LANGUAGE:C>:-opt-report-file=stderr>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-opt-report-file=stderr>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-opt-report-file=stderr>")

add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-mkl=sequential>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>")

