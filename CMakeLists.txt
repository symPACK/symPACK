cmake_minimum_required(VERSION 3.1...3.12)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()


project(symPACK VERSION 2.0
                DESCRIPTION "symPACK is a direct linear solver for sparse symmetric matrices"
                LANGUAGES C CXX Fortran)

              #include(ExternalProject)


############################## Source files ######################################

set(SYMPACKSRC
src/ETree.cpp
src/utility.cpp
src/timer.cpp
src/blas.cpp
src/lapack.cpp
src/global.cpp
src/LogFile.cpp
src/IntervalTree.cpp
src/Ordering.cpp
src/CommPull.cpp
src/symPACK_c_interface.cpp
src/DistSparseMatrixGraph.cpp
src/amdbar.F
src/genmmd.F
src/mmdelm.F
src/mmdint.F
src/mmdnum.F
src/mmdupd.F
src/ordmmd.F
src/betree.f
src/del_heap.f
src/ins_heap.f
src/mcssup.f
src/move_down.f
src/ordsup.f
src/boxnd.f
src/gridnd.f
src/degree.F
src/fnroot.F
src/genrcm.F
src/rcm.F
src/rootls.F
src/TSPB/compress_skel2.f
src/TSPB/fnskel2.f
src/TSPB/near_ins2.f
src/TSPB/ordsup_ind_tsp_paths2.f
src/TSPB/ordsup_tsp_paths2.f
src/TSPB/post_suptree2.f
src/TSPB/sup_etpost2.f
src/TSPB/fnskel.f
src/TSPB/near_ins.f
src/TSPB/ordsup_ind_tsp_paths.f
src/TSPB/ordsup_tsp_paths.f
src/TSPB/post_suptree.f
src/TSPB/sup_etpost.f
src/invinv.f
src/mod_heap.f
src/move_up.f
)



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
add_library(sympack ${SYMPACKSRC})
target_include_directories(sympack PUBLIC include)

add_executable(run_sympack EXCLUDE_FROM_ALL driver/run_sympack.cpp)
add_executable(run_sympack2D EXCLUDE_FROM_ALL driver/run_sympack2D.cpp)
add_executable(run_sympack_c EXCLUDE_FROM_ALL driver/run_c_interface.cpp)
add_custom_target(examples DEPENDS run_sympack run_sympack_c)
########################### Dependencies ####################################




set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

  
find_package(MPI REQUIRED)
if(MPI_CXX_FOUND)
  message(STATUS "Found CXX MPI compiler: ${MPI_CXX_COMPILER}")
  message( STATUS "MPI_CXX_LINK_FLAGS: " ${MPI_CXX_LINK_FLAGS} )
  message( STATUS "MPI_CXX_LIBRARIES: " ${MPI_CXX_LIBRARIES} )
  target_link_libraries(sympack PUBLIC MPI::MPI_CXX)
endif()


include(FortranCInterface)
FortranCInterface_HEADER("${PROJECT_BINARY_DIR}/sympack_config.hpp" MACRO_NAMESPACE "FC_")
FortranCInterface_VERIFY(CXX)


#######   ORDERING LIBRARIES ######

if(ENABLE_PARMETIS)
  find_package(ParMETIS REQUIRED)
  target_link_libraries(sympack PUBLIC ParMETIS::parmetis)
  target_compile_definitions( sympack PUBLIC "-DUSE_PARMETIS")
  set(ENABLE_METIS ON)
endif()

if(ENABLE_METIS)
  if (NOT METIS_FOUND)
    find_package(METIS REQUIRED)
    target_link_libraries(sympack PUBLIC METIS::metis)
  endif()
  target_compile_definitions( sympack PUBLIC "-DUSE_METIS")
endif()

if(ENABLE_PTSCOTCH)
  find_package(PTSCOTCH REQUIRED)
  target_link_libraries(sympack PUBLIC PTSCOTCH::ptscotch)
  target_compile_definitions( sympack PUBLIC "-DUSE_PTSCOTCH")
  set(ENABLE_SCOTCH ON)
endif()

if(ENABLE_SCOTCH)
  if (NOT SCOTCH_FOUND)
    find_package(SCOTCH REQUIRED)
    target_link_libraries(sympack PUBLIC SCOTCH::scotch)
  endif()
  target_compile_definitions( sympack PUBLIC "-DUSE_SCOTCH")
endif()


if(ENABLE_MKL)
  find_package(MKL)
  if(MKL_FOUND)
    message(STATUS "Found MKL library: " ${MKL_LIBRARIES})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MKL_LINKER_FLAGS}")
    set(LIBS ${LIBS} ${MKL_LIBRARIES})
  else()
    message(WARNING "Warning: could not find a MKL library")
    add_compile_options("-mkl=sequential")
    SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -mkl=sequential" )
  endif()
else()
  set(MKL_FOUND false)
  get_directory_property( ALL_OPTIONS DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_OPTIONS )
  foreach( d ${ALL_OPTIONS} )
    string(FIND ${d} "-mkl" pos)
    if(NOT pos EQUAL "-1")
      set(MKL_FOUND true)
    endif()
  endforeach()
  
  if(MKL_FOUND)
    SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -mkl=sequential" )
    message(STATUS "Found MKL in flags")
  else()
    find_package(BLAS )
    if(BLAS_FOUND)
      message(STATUS "Found BLAS library: " ${BLAS_LIBRARIES})
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")
      set(LIBS ${LIBS} ${BLAS_LIBRARIES})
    else()
      message(WARNING "Warning: could not find a BLAS library")
    endif()
    
    find_package(LAPACK )
    if(LAPACK_FOUND)
      message(STATUS "Found LAPACK library: " ${LAPACK_LIBRARIES})
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
      set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
    else()
      message(WARNING "Warning: could not find a LAPACK library")
    endif()
  endif()
endif()


#######   SYMPACK OPTIONS ######

if(TRACK_MEMORY)
  target_compile_definitions( sympack PUBLIC "-D_TRACK_MEMORY_")
endif()

if(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  if(OPENMP_FOUND)
    message(STATUS "OpenMP found")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

if(ENABLE_PROFILE)
  if(NOT ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OPENMP_FOUND)
      message(STATUS "OpenMP found")
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
  endif()

  target_compile_definitions( sympack PUBLIC "-DSPROFILE")
  target_compile_definitions( sympack PUBLIC "-DNO_INTRA_PROFILE")
endif()


if(RCM_IDX_64)
  target_compile_definitions( sympack PUBLIC  "-DRCMIDXSIZE=64")
else()
  target_compile_definitions( sympack PUBLIC  "-DRCMIDXSIZE=32")
endif()

if(MMD_IDX_64)
  target_compile_definitions( sympack PUBLIC  "-DMMDIDXSIZE=64")
else()
  target_compile_definitions( sympack PUBLIC  "-DMMDIDXSIZE=32")
endif()

if(AMD_IDX_64)
  target_compile_definitions( sympack PUBLIC  "-DAMDIDXSIZE=64")
else()
  target_compile_definitions( sympack PUBLIC  "-DAMDIDXSIZE=32")
endif()






find_package(UPCXX REQUIRED)
target_link_libraries(sympack PUBLIC UPCXX::upcxx)
  


target_link_libraries(sympack PUBLIC ${LIBS})
target_link_libraries(sympack PUBLIC ${CMAKE_DL_LIBS})

############################ Compile Options ################################

target_compile_features(sympack PUBLIC cxx_std_11)
set_target_properties(sympack PROPERTIES POSITION_INDEPENDENT_CODE ON)

#find_library(MATH_LIBRARY m)
#if(MATH_LIBRARY)
#   target_link_libraries(sympack PUBLIC ${MATH_LIBRARY})
#endif()

target_include_directories(sympack PRIVATE ${PROJECT_BINARY_DIR})

target_compile_definitions(sympack PUBLIC "ITREE")
target_compile_definitions(sympack PUBLIC "_DSW_ITREE_")
#add_definitions("-DNOTRY")


###############################################
##  generate a config file
###############################################
get_target_property( ALL_DEFINITIONS sympack COMPILE_DEFINITIONS )
file(WRITE "sympack_definitions.hpp.cmake" "#ifndef _HEADER_SYMPACK_DEFINITIONS_\n#define _HEADER_SYMPACK_DEFINITIONS_\n" )
foreach( d ${ALL_DEFINITIONS} )
  string(REPLACE "=" " " d ${d})
  string(REPLACE " " ";" DEF_LIST ${d})
  list(GET DEF_LIST 0 def_name)

  file(APPEND "sympack_definitions.hpp.cmake" "\n#ifndef ${def_name}\n" )
  file(APPEND "sympack_definitions.hpp.cmake" "#define ${d}\n" )
  file(APPEND "sympack_definitions.hpp.cmake" "#endif // ${def_name}\n" )
endforeach()
file(APPEND "sympack_definitions.hpp.cmake" "\n#endif //_HEADER_SYMPACK_DEFINITIONS_\n" )

CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/sympack_definitions.hpp.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/sympack_definitions.hpp
)

####################################################
####### generate a makefile include file
####################################################
#####file(WRITE "sympack.mak.cmake" "SYMPACK_DIR = ${CMAKE_INSTALL_PREFIX}\n\n")
#####file(APPEND "sympack.mak.cmake" "SYMPACK_INCLUDE = -I\$(SYMPACK_DIR)/include\n")
#####
#####if(ENABLE_SCOTCH)
#####  file(APPEND "sympack.mak.cmake" "SYMPACK_INCLUDE += -I${SCOTCH_INCLUDE_DIRS}\n")
#####endif()
#####if(ENABLE_METIS)
#####  file(APPEND "sympack.mak.cmake" "SYMPACK_INCLUDE += -I${METIS_INCLUDE_DIRS}\n")
#####endif()
#####if(ENABLE_PARMETIS)
#####  file(APPEND "sympack.mak.cmake" "SYMPACK_INCLUDE += -I${PARMETIS_INCLUDE_DIRS}\n")
#####endif()
#####
#####file(APPEND "sympack.mak.cmake" "SYMPACK_LIB = -L\$(SYMPACK_DIR)/lib -lsympack\n")
#####
#####foreach( d ${LIBS} )
#####  file(APPEND "sympack.mak.cmake" "SYMPACK_LIB += ${d}\n")
#####endforeach()
#####
#####IF( NOT CMAKE_DL_LIBS STREQUAL "")
#####file(APPEND "sympack.mak.cmake" "SYMPACK_LIB += -l${CMAKE_DL_LIBS}\n")
#####endif()
#####
#####unset(TMPUPCXXDIR)
#####
#####CONFIGURE_FILE(
#####    ${CMAKE_CURRENT_SOURCE_DIR}/sympack.mak.cmake
#####    ${CMAKE_CURRENT_BINARY_DIR}/sympack.mak
#####)
#####
#####



      #set(CMAKE_CXX_FLAGS "${UPCXX_DEFINES} ${CMAKE_CXX_FLAGS}")
      #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${UPCXX_LDFLAGS}")




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




###########################################################################


 
if(UNIX AND NOT APPLE)
  target_link_libraries(sympack PUBLIC rt)
endif()


install(TARGETS sympack LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
install(FILES 
  ${PROJECT_BINARY_DIR}/sympack_config.hpp 
  ${PROJECT_BINARY_DIR}/sympack_definitions.hpp
  ${PROJECT_BINARY_DIR}/sympack.mak
  include/sympack.hpp
  DESTINATION include)

install(FILES 
  include/sympack/SuperNode.hpp
  include/sympack/SuperNodeInd.hpp
  include/sympack/CommPull.hpp
  include/sympack/CommTypes.hpp
  include/sympack/DistSparseMatrix.hpp
  include/sympack/symPACKMatrix.hpp
  include/sympack/symPACKMatrixBase.hpp
  include/sympack/symPACKMatrix2D.hpp
  include/sympack/DistSparseMatrixGraph.hpp
  include/sympack/ETree.hpp
  include/sympack/Environment.hpp
  include/sympack/IntervalTree.hpp
  include/sympack/LoadBalancer.hpp 
  include/sympack/Task.hpp
  include/sympack/LogFile.hpp
  include/sympack/Types.hpp
  include/sympack/Mapping.hpp 
  include/sympack/blas.hpp
  include/sympack/datatypes.hpp
  include/sympack/lapack.hpp
  include/sympack/mpi_interf.hpp
  include/sympack/timer.hpp
  include/sympack/Ordering.hpp
  include/sympack/Scheduler.hpp
  include/sympack/utility.hpp 
  include/sympack/supernodalTaskGraph.hpp
  DESTINATION include/sympack)

install(FILES 
  include/sympack/impl/symPACKMatrixBase_impl.hpp
  include/sympack/impl/symPACKMatrix_impl.hpp
  include/sympack/impl/symPACKMatrix2D_impl.hpp
  include/sympack/impl/SuperNode_impl.hpp
  include/sympack/impl/SuperNodeInd_impl.hpp
  include/sympack/impl/DistSparseMatrix_impl.hpp
  include/sympack/impl/symPACKMatrix_impl_FB_pull.hpp
  include/sympack/impl/symPACKMatrix_impl_solve.hpp
  include/sympack/impl/Environment_impl.hpp 
  include/sympack/impl/IntervalTree_impl.hpp 
  include/sympack/impl/supernodalTaskGraph_impl.hpp
  include/sympack/impl/Scheduler_impl.hpp
  DESTINATION include/sympack/impl)


#target_link_libraries(run_sympack sympack libupcxx libgasnet-conduit ${LIBS}  ${GASNET_LIBRARIES} ${CMAKE_DL_LIBS})



#add_dependencies(run_sympack sympack)# gasnet upcxx) 
#add_dependencies(run_sympack2D sympack)# gasnet upcxx) 
#add_dependencies(run_sympack_c sympack)# gasnet upcxx) 

target_link_libraries(run_sympack PUBLIC sympack)
target_link_libraries(run_sympack2D PUBLIC sympack)
target_link_libraries(run_sympack_c PUBLIC sympack)

find_package( Threads )
target_link_libraries( run_sympack PUBLIC ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries( run_sympack2D PUBLIC ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries( run_sympack_c PUBLIC ${CMAKE_THREAD_LIBS_INIT} )

#install(TARGETS run_sympack RUNTIME DESTINATION bin)

###########################################################################

set(SYMPACK_VERSION_MAJOR 1)
set(SYMPACK_VERSION_MINOR 0)
set(SYMPACK_VERSION_PATCH 0)
set(SYMPACK_VERSION_FULL "${SYMPACK_VERSION_MAJOR}.${SYMPACK_VERSION_MINOR}")
set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION_MAJOR ${SYMPACK_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${SYMPACK_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${SYMPACK_VERSION_PATCH})
#SET(CPACK_PACKAGE_FILE_NAME "sympack-"${SYMPACK_VERSION_FULL})
include(CPack)


