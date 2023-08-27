cmake_minimum_required(VERSION 3.11...3.12)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()


if (ENABLE_CUDA)
    project(symPACK VERSION 2.0
                DESCRIPTION "symPACK is a direct linear solver for sparse symmetric matrices"
                LANGUAGES C CXX Fortran CUDA)
else()
    project(symPACK VERSION 2.0
                DESCRIPTION "symPACK is a direct linear solver for sparse symmetric matrices"
                LANGUAGES C CXX Fortran)
endif()

              #include(ExternalProject)



############################## Options ######################################

#all options
option(ENABLE_PTSCOTCH "Enable PTSCOTCH" OFF)
option(ENABLE_SCOTCH "Enable SCOTCH" OFF)
option(ENABLE_PARMETIS "Enable PARMETIS" OFF)
option(ENABLE_METIS "Enable METIS" OFF)
option(ENABLE_MKL "Enable MKL" OFF)
option(AMD_IDX_64 "Use 64 bit integers for AMD" OFF)
option(MMD_IDX_64 "Use 64 bit integers for MMD" OFF)
option(RCM_IDX_64 "Use 64 bit integers for RCM" OFF)
option(ENABLE_THREADS "Enable Threading" OFF)



############################## Targets ######################################

add_executable(run_sympack EXCLUDE_FROM_ALL driver/run_sympack.cpp)
#add_executable(run_sympack2D EXCLUDE_FROM_ALL driver/run_sympack2D.cpp)
add_custom_target(examples DEPENDS run_sympack run_sympack2D )

#target_include_directories( run_sympack2D PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/driver> )
target_include_directories( run_sympack PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/driver> )


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
#target_link_libraries(run_sympack2D PUBLIC symPACK::sympack)
IF(ENABLE_CUDA)
	target_link_libraries(run_sympack PUBLIC symPACK::sympack -lcublas -lcusolver)
	#target_link_libraries(run_sympack2D  PUBLIC symPACK::sympack -lcublas -lcusolver)
ENDIF()
find_package( Threads )
target_link_libraries( run_sympack PUBLIC ${CMAKE_THREAD_LIBS_INIT} )
#target_link_libraries( run_sympack2D PUBLIC ${CMAKE_THREAD_LIBS_INIT} )


###########################################################################

