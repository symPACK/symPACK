cmake_minimum_required(VERSION 3.1...3.12)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()


project(symPACK VERSION 2.0
                DESCRIPTION "symPACK is a direct linear solver for sparse symmetric matrices"
                LANGUAGES C CXX Fortran)

              #include(ExternalProject)



############################## Options ######################################

#all options
option(GASNET_PAR "Use Gasnet Parallel conduits" OFF)
option(ENABLE_ARIES "Enable ARIES" OFF)
option(ENABLE_KNL "Cross compiling" OFF)
option(ENABLE_KNL_ONLY "Cross compiling" OFF)
option(ENABLE_PROFILE "Enable profiling" OFF)
option(ENABLE_OPENMP "Enable OpenMP" OFF)
option(TRACK_MEMORY "Enable (partial) Memory allocation tracking" OFF)
option(ENABLE_PTSCOTCH "Enable PTSCOTCH" OFF)
option(ENABLE_SCOTCH "Enable SCOTCH" OFF)
option(ENABLE_PARMETIS "Enable PARMETIS" OFF)
option(ENABLE_METIS "Enable METIS" OFF)
option(ENABLE_MKL "Enable MKL" OFF)
option(AMD_IDX_64 "Use 64 bit integers for AMD" OFF)
option(MMD_IDX_64 "Use 64 bit integers for MMD" OFF)
option(RCM_IDX_64 "Use 64 bit integers for RCM" OFF)



############################## Targets ######################################

add_executable(run_sympack EXCLUDE_FROM_ALL driver/run_sympack.cpp)
add_executable(run_sympack2D EXCLUDE_FROM_ALL driver/run_sympack2D.cpp)
add_executable(run_sympack_c EXCLUDE_FROM_ALL driver/run_c_interface.cpp)
add_custom_target(examples DEPENDS run_sympack run_sympack_c)

target_include_directories( run_sympack2D PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/driver> )


########################### Dependencies ####################################




#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")
list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")

#include( PEXSICompileFlags ) # Compile Flags
include( symPACKDepends ) # Basic Dependencies


# libsympack
add_subdirectory( src )


# examples
#add_subdirectory( driver )



string( TOLOWER "${CMAKE_BUILD_TYPE}" config_type )
IF( config_type STREQUAL "debug")
  add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-g>")
  add_compile_options("$<$<COMPILE_LANGUAGE:C>:-g>")
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-g>")
  add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-O0>")
  add_compile_options("$<$<COMPILE_LANGUAGE:C>:-O0>")
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-O0>")
  message("Debug mode")
ENDIF()  




target_link_libraries(run_sympack PUBLIC   symPACK::sympack)
target_link_libraries(run_sympack2D PUBLIC symPACK::sympack)
target_link_libraries(run_sympack_c PUBLIC symPACK::sympack)

find_package( Threads )
target_link_libraries( run_sympack PUBLIC ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries( run_sympack2D PUBLIC ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries( run_sympack_c PUBLIC ${CMAKE_THREAD_LIBS_INIT} )


###########################################################################

