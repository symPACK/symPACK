cmake_minimum_required(VERSION 2.8)
project(SYMPACK C CXX Fortran)
include(ExternalProject)

#all options
option(GASNET_PAR "Use Gasnet Parallel conduits" OFF)
option(ENABLE_ARIES "Enable ARIES" OFF)
option(ENABLE_KNL "Cross compiling" OFF)
option(ENABLE_KNL_ONLY "Cross compiling" OFF)
option(ENABLE_PROFILE "Enable profiling" OFF)
option(ENABLE_OPENMP "Enable OpenMP" OFF)
option(TRACK_MEMORY "Enable (partial) Memory allocation tracking" OFF)
option(ENABLE_SCOTCH "Enable SCOTCH" OFF)
option(ENABLE_PARMETIS "Enable PARMETIS" OFF)
option(ENABLE_METIS "Enable METIS" OFF)
option(ENABLE_MKL "Enable MKL" OFF)
option(AMD_IDX_64 "Use 64 bit integers for AMD" OFF)
option(MMD_IDX_64 "Use 64 bit integers for MMD" OFF)
option(RCM_IDX_64 "Use 64 bit integers for RCM" OFF)

option(NEW_UPCXX "Use new UPCXX" OFF)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)

string( TOLOWER "${CMAKE_BUILD_TYPE}" config_type )
IF( "${config_type}" STREQUAL "debug")
  message("Debug mode")
ENDIF()  
  
find_package(MPI REQUIRED)
if(MPI_CXX_FOUND)
  message(STATUS "Found CXX MPI compiler: ${MPI_CXX_COMPILER}")
  include_directories(${MPI_CXX_INCLUDE_PATH})
  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS})
  set(LIBS ${LIBS} ${MPI_CXX_LIBRARIES})
endif()

if(MPI_C_FOUND)
  message(STATUS "Found C MPI compiler: ${MPI_C_COMPILER}")
  include_directories(${MPI_C_INCLUDE_PATH})
  set(CMAKE_CXX_FLAGS ${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS})
  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS})
  set(LIBS ${LIBS} ${MPI_C_LIBRARIES})
endif()



include(FortranCInterface)
FortranCInterface_HEADER("${PROJECT_BINARY_DIR}/sympack_config.hpp" MACRO_NAMESPACE "FC_")
FortranCInterface_VERIFY(CXX)

if(RCM_IDX_64)
  add_definitions( "-DRCMIDXSIZE=64")
else()
  add_definitions( "-DRCMIDXSIZE=32")
endif()

if(MMD_IDX_64)
  add_definitions( "-DMMDIDXSIZE=64")
else()
  add_definitions( "-DMMDIDXSIZE=32")
endif()

if(AMD_IDX_64)
  add_definitions( "-DAMDIDXSIZE=64")
else()
  add_definitions( "-DAMDIDXSIZE=32")
endif()

if(ENABLE_MKL)
  find_package(MKL)
  if(MKL_FOUND)
    message(STATUS "Found MKL library: " ${MKL_LIBRARIES})
    set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${MKL_LINKER_FLAGS})
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
      set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS})
      set(LIBS ${LIBS} ${BLAS_LIBRARIES})
    else()
      message(WARNING "Warning: could not find a BLAS library")
    endif()
    
    find_package(LAPACK )
    if(LAPACK_FOUND)
      message(STATUS "Found LAPACK library: " ${LAPACK_LIBRARIES})
      set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
      set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
    else()
      message(WARNING "Warning: could not find a LAPACK library")
    endif()
  endif()
endif()


#######   ORDERING LIBRARIES ######

if(ENABLE_METIS)
  find_package(METIS REQUIRED)
  set(LIBS ${LIBS} ${METIS_LIBRARIES})
  include_directories(${METIS_INCLUDE_DIRS})
  add_definitions("-DUSE_METIS")
endif()

if(ENABLE_PARMETIS)
  find_package(ParMETIS REQUIRED)
  set(LIBS ${LIBS} ${PARMETIS_LIBRARIES})
  include_directories(${PARMETIS_INCLUDE_DIRS})
  add_definitions("-DUSE_METIS")
  add_definitions("-DUSE_PARMETIS")
endif()

if(ENABLE_SCOTCH)
  find_package(Scotch REQUIRED)
  set(LIBS ${LIBS} ${SCOTCH_LIBRARIES})
  include_directories(${SCOTCH_INCLUDE_DIRS})
  add_definitions("-DUSE_SCOTCH")
  add_definitions("-DUSE_PTSCOTCH")
endif()


#######   SYMPACK OPTIONS ######

if(TRACK_MEMORY)
  add_definitions("-D_TRACK_MEMORY_")
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


  add_definitions("-DSPROFILE")
  add_definitions("-DNO_INTRA_PROFILE")
endif()


if(NEW_UPCXX)
  #set(upcxx-meta "/home/mjacquel/Work/lbl/pagoda/upcxx_dev/upcxx_install/bin/upcxx-meta") 
  set(upcxx-meta "upcxx-meta") 
  execute_process( COMMAND ${upcxx-meta} PPFLAGS OUTPUT_VARIABLE UPCXX_PPFLAGS)
  message( STATUS "UPCXX_PPFLAGS: " ${UPCXX_PPFLAGS} )
  execute_process( COMMAND ${upcxx-meta} LDFLAGS OUTPUT_VARIABLE UPCXX_LDFLAGS)
  message( STATUS "UPCXX_LDFLAGS: " ${UPCXX_LDFLAGS} )
  execute_process( COMMAND ${upcxx-meta} LIBFLAGS OUTPUT_VARIABLE UPCXX_LIBFLAGS)
  message( STATUS "UPCXX_LIBFLAGS: " ${UPCXX_LIBFLAGS} )
  string(REPLACE "\n" " " UPCXX_LIBFLAGS ${UPCXX_LIBFLAGS})
  string(REPLACE "\n" " " UPCXX_PPFLAGS ${UPCXX_PPFLAGS})
  string(REPLACE "\n" " " UPCXX_LDFLAGS ${UPCXX_LDFLAGS})

  string(STRIP ${UPCXX_LIBFLAGS} UPCXX_LIBFLAGS)
  string(STRIP ${UPCXX_PPFLAGS}  UPCXX_PPFLAGS)
  string(STRIP ${UPCXX_LDFLAGS}  UPCXX_LDFLAGS)
  set(UPCXX_DEFINES "${UPCXX_PPFLAGS} ${UPCXX_LDFLAGS} -DNEW_UPCXX")
  set(UPCXX_LIBRARIES "${UPCXX_LIBFLAGS}")

  #string(REPLACE " " ";" UPCXX_LIBRARIES ${UPCXX_LIBRARIES})
  #message( STATUS "UPCXX_LIBRARIES: " ${UPCXX_LIBRARIES} )
  #foreach( d ${UPCXX_LIBRARIES} )
  #  string(STRIP ${d} d)
  #  message( STATUS "Found UPCXX Define: " ${d} )
  #  add_definitions("${d}")
  #endforeach()

  #endif()


  
else()
#######   GASNET AND UPCXX OPTIONS ######

  include(cmake/gasnet.cmake)
  include(cmake/upcxx.cmake)


include_directories(${UPCXX_INCLUDE_DIR})
link_directories(${UPCXX_LIBRARY_PATH})

if(GASNET_DEFINES)
  string(REPLACE " " ";" GASNET_DEFINES ${GASNET_DEFINES})
  message( STATUS "GASNET_DEFINES: " ${GASNET_DEFINES} )
  foreach( d ${GASNET_DEFINES} )
    string(STRIP ${d} d)
    message( STATUS "Found GASNET Define: " ${d} )
    add_definitions("${d}")
  endforeach()
endif()

#this work with both old and new upcxx
if(UPCXX_DEFINES)
  string(REPLACE " " ";" UPCXX_DEFINES ${UPCXX_DEFINES})
  message( STATUS "UPCXX_DEFINES: " ${UPCXX_DEFINES} )
  foreach( d ${UPCXX_DEFINES} )
    string(STRIP ${d} d)
    message( STATUS "Found UPCXX Define: " ${d} )
    add_definitions("${d}")
  endforeach()
endif()

endif()

set(SYMPACKSRC src/ETree.cpp src/utility.cpp src/timer.cpp src/blas.cpp
  src/lapack.cpp src/global.cpp src/LogFile.cpp 
  src/IntervalTree.cpp src/Ordering.cpp src/CommPull.cpp
  src/symPACK_c_interface.cpp
  src/DistSparseMatrixGraph.cpp src/amdbar.F src/genmmd.F src/mmdelm.F
  src/mmdint.F src/mmdnum.F src/mmdupd.F src/ordmmd.F
  src/betree.f src/del_heap.f src/ins_heap.f src/mcssup.f
  src/move_down.f src/ordsup.f src/boxnd.f src/gridnd.f
  src/degree.F src/fnroot.F src/genrcm.F src/rcm.F src/rootls.F
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
  src/invinv.f src/mod_heap.f src/move_up.f)

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_BINARY_DIR})

add_definitions("-DITREE")
add_definitions("-D_DSW_ITREE_")
#add_definitions("-DNOTRY")


###############################################
##  generate a config file
###############################################
get_directory_property( ALL_DEFINITIONS DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )
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

###############################################
## generate a makefile include file
###############################################
file(WRITE "sympack.mak.cmake" "SYMPACK_DIR = ${CMAKE_INSTALL_PREFIX}\n\n")
file(APPEND "sympack.mak.cmake" "SYMPACK_INCLUDE = -I\$(SYMPACK_DIR)/include\n")

if(ENABLE_SCOTCH)
  file(APPEND "sympack.mak.cmake" "SYMPACK_INCLUDE += -I${SCOTCH_INCLUDE_DIRS}\n")
endif()
if(ENABLE_METIS)
  file(APPEND "sympack.mak.cmake" "SYMPACK_INCLUDE += -I${METIS_INCLUDE_DIRS}\n")
endif()
if(ENABLE_PARMETIS)
  file(APPEND "sympack.mak.cmake" "SYMPACK_INCLUDE += -I${PARMETIS_INCLUDE_DIRS}\n")
endif()



if(NEW_UPCXX)
else()
list(GET UPCXX_INCLUDE_DIR 0 TMPUPCXXDIR)

set(UPCXX_PREFIX $ENV{UPCXX_DIR})
if (UPCXX_PREFIX)
  file(APPEND "sympack.mak.cmake" "\ninclude ${TMPUPCXXDIR}/upcxx.mak\n")
  file(APPEND "sympack.mak.cmake" "VPATH=${TMPUPCXXDIR}\n\n")
else()
  file(APPEND "sympack.mak.cmake" "\ninclude \$(SYMPACK_DIR)/upcxx_install/include/upcxx.mak\n")
  file(APPEND "sympack.mak.cmake" "VPATH=\$(SYMPACK_DIR)/upcxx_install\n\n")
endif()

file(APPEND "sympack.mak.cmake" "SYMPACK_INCLUDE += \$(UPCXX_CXXFLAGS)\n\n")
endif()

file(APPEND "sympack.mak.cmake" "SYMPACK_LIB = -L\$(SYMPACK_DIR)/lib -lsympack\n")

foreach( d ${LIBS} )
  file(APPEND "sympack.mak.cmake" "SYMPACK_LIB += ${d}\n")
endforeach()

if(NEW_UPCXX)
else()
file(APPEND "sympack.mak.cmake" "SYMPACK_LIB += \$(UPCXX_LDLIBS)\n")
endif()

file(APPEND "sympack.mak.cmake" "SYMPACK_LIB += -l${CMAKE_DL_LIBS}\n")

unset(TMPUPCXXDIR)

CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/sympack.mak.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/sympack.mak
)


add_library(sympack ${SYMPACKSRC})



if(NEW_UPCXX) 
      set(CMAKE_CXX_FLAGS "${UPCXX_DEFINES} ${CMAKE_CXX_FLAGS}")
      #      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${UPCXX_LIBFLAGS}")
else()
add_dependencies(sympack gasnet upcxx)
set_target_properties(gasnet PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(upcxx PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

set_target_properties(sympack PROPERTIES POSITION_INDEPENDENT_CODE ON)
 
if(UNIX AND NOT APPLE)
    target_link_libraries(sympack rt)
endif()


if(NEW_UPCXX)
target_link_libraries(sympack ${LIBS} ${CMAKE_DL_LIBS})
else()
target_link_libraries(sympack libupcxx libgasnet-conduit ${LIBS} ${UPCXX_LIBRARIES} ${GASNET_LIBRARIES} ${CMAKE_DL_LIBS})
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
  include/sympack/impl/symPACKMatrix_impl.hpp
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

if(NEW_UPCXX)
else()
install (SCRIPT "${CMAKE_SOURCE_DIR}/cmake/gasnetPostInstall.cmake")
install (SCRIPT "${CMAKE_SOURCE_DIR}/cmake/upcxxPostInstall.cmake")
endif()

add_executable(run_sympack EXCLUDE_FROM_ALL driver/run_sympack.cpp)
add_executable(run_sympack2D EXCLUDE_FROM_ALL driver/run_sympack2D.cpp)
add_executable(run_sympack_c EXCLUDE_FROM_ALL driver/run_c_interface.cpp)
add_custom_target(examples DEPENDS run_sympack run_sympack_c)

add_dependencies(run_sympack sympack)# gasnet upcxx) 
add_dependencies(run_sympack2D sympack)# gasnet upcxx) 
add_dependencies(run_sympack_c sympack)# gasnet upcxx) 
#target_link_libraries(run_sympack sympack libupcxx libgasnet-conduit ${LIBS}  ${GASNET_LIBRARIES} ${CMAKE_DL_LIBS})




if(NEW_UPCXX)
  target_link_libraries(run_sympack sympack ${UPCXX_LIBFLAGS})
  target_link_libraries(run_sympack2D sympack ${UPCXX_LIBFLAGS})
target_link_libraries(run_sympack_c sympack)
else()
target_link_libraries(run_sympack sympack libupcxx libgasnet-conduit)
target_link_libraries(run_sympack2D sympack libupcxx libgasnet-conduit)
target_link_libraries(run_sympack_c sympack libupcxx libgasnet-conduit)
endif()


find_package( Threads )
target_link_libraries( run_sympack ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries( run_sympack2D ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries( run_sympack_c ${CMAKE_THREAD_LIBS_INIT} )

#install(TARGETS run_sympack RUNTIME DESTINATION bin)

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


