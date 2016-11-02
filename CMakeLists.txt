cmake_minimum_required(VERSION 3.3)
project(SYMPACK CXX C Fortran)
include(ExternalProject)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)

#SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0" CACHE STRING "" FORCE)
#SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0" CACHE STRING "" FORCE)
#SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0" CACHE STRING "" FORCE)
add_compile_options("$<$<CONFIG:Debug>:-O0>")
add_compile_options("$<$<CONFIG:Debug>:-g>")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
#include_directories(/path/to/upcxx)
#set(LIBS ${LIBS} )



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

option(RCM_IDX_64 "Use 64 bit integers for RCM" OFF)
if(RCM_IDX_64)
add_definitions( "-DRCMIDXSIZE=64")
else()
add_definitions( "-DRCMIDXSIZE=32")
endif()

option(MMD_IDX_64 "Use 64 bit integers for MMD" OFF)
if(MMD_IDX_64)
add_definitions( "-DMMDIDXSIZE=64")
else()
add_definitions( "-DMMDIDXSIZE=32")
endif()

option(AMD_IDX_64 "Use 64 bit integers for AMD" OFF)
if(AMD_IDX_64)
add_definitions( "-DAMDIDXSIZE=64")
else()
add_definitions( "-DAMDIDXSIZE=32")
endif()

option(ENABLE_MKL "Enable MKL" OFF)
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
    if(pos EQUAL "0")
      set(MKL_FOUND true)
    endif()
  endforeach()
  
  
  if(MKL_FOUND)
    SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -mkl=sequential" )
    message(STATUS "Found MKL in flags")
  else()
    find_package(BLAS)
    if(BLAS_FOUND)
      message(STATUS "Found BLAS library: " ${BLAS_LIBRARIES})
      set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS})
      set(LIBS ${LIBS} ${BLAS_LIBRARIES})
    else()
      message(WARNING "Warning: could not find a BLAS library")
    endif()
    
    find_package(LAPACK)
    if(LAPACK_FOUND)
      message(STATUS "Found LAPACK library: " ${LAPACK_LIBRARIES})
      set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
      set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
    else()
      message(WARNING "Warning: could not find a LAPACK library")
    endif()
  endif()
endif()


#find_package(MKL)
#if(MKL_FOUND)
#  message(STATUS "Found MKL library: " ${MKL_LIBRARIES})
#  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${MKL_LINKER_FLAGS})
#  set(LIBS ${LIBS} ${MKL_LIBRARIES})
#else()
#  message(WARNING "Warning: could not find a MKL library")
#endif()

#find_package(LIBSCI)
#if(LIBSCI_FOUND)
#  message(STATUS "Found LIBSCI library: " ${LIBSCI_LIBRARIES})
#  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${LIBSCI_LINKER_FLAGS})
#  set(LIBS ${LIBS} ${LIBSCI_LIBRARIES})
#else()
#  message(WARNING "Warning: could not find a LIBSCI library")
#endif()





#find_package(GASNET)
#if(GASNET_FOUND)
#  message(STATUS "Found GASNET library: " ${GASNET_LIBRARIES})
#  include_directories(${GASNET_INCLUDE_DIR})
#  set(LIBS ${LIBS} ${GASNET_LIBRARIES})
#endif()


#option(ENABLE_GASNET "Enable GASNET and dependencies" ON)
#cmake_dependent_option(ENABLE_GASNET "Enable gasnet" ON
#                       "ENABLE_GASNET" OFF)


##MACRO( CONFIGURE_UPCXX )
##  find_package(OpenMP REQUIRED)
#find_package(UPCXX REQUIRED)
#if(UPCXX_FOUND)
#  message(STATUS "Found UPC++ library: " ${UPCXX_LIBRARIES})
#  add_definitions(${UPCXX_DEFINES})
#  include_directories(${UPCXX_INCLUDE_DIR})
#  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
#  set(LIBS ${LIBS} ${UPCXX_LIBRARIES})
##LINK_DIRECTORIES(${UPCXX_LIBRARY_PATH})
#endif()
##ENDMACRO( CONFIGURE_UPCXX )



#######   ORDERING LIBRARIES ######

option(ENABLE_METIS "Enable METIS" OFF)
if(ENABLE_METIS)
  find_package(METIS REQUIRED)
  set(LIBS ${LIBS} ${METIS_LIBRARIES})
  include_directories(${METIS_INCLUDE_DIRS})
  add_definitions("-DUSE_METIS")
endif()

option(ENABLE_PARMETIS "Enable PARMETIS" OFF)
if(ENABLE_PARMETIS)
  find_package(ParMETIS REQUIRED)
  set(LIBS ${LIBS} ${PARMETIS_LIBRARIES})
  include_directories(${PARMETIS_INCLUDE_DIRS})
  add_definitions("-DUSE_METIS")
  add_definitions("-DUSE_PARMETIS")
endif()

option(ENABLE_SCOTCH "Enable SCOTCH" OFF)
if(ENABLE_SCOTCH)
  find_package(Scotch REQUIRED)
  set(LIBS ${LIBS} ${SCOTCH_LIBRARIES})
  include_directories(${SCOTCH_INCLUDE_DIRS})
  add_definitions("-DUSE_SCOTCH")
  add_definitions("-DUSE_PTSCOTCH")
endif()


#######   SYMPACK OPTIONS ######

option(TRACK_MEMORY "Enable (partial) Memory allocation tracking" OFF)
if(TRACK_MEMORY)
  add_definitions("-D_TRACK_MEMORY_")
endif()


option(ENABLE_PROFILE "Enable profiling" OFF)
if(ENABLE_PROFILE)
  find_package(OpenMP REQUIRED)
  if(OPENMP_FOUND)
    message(STATUS "OpenMP found")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:${OpenMP_CXX_FLAGS}>")
#add_compile_options("$<$<COMPILE_LANGUAGE:C>:${OpenMP_C_FLAGS}>")
#add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:${OpenMP_Fortran_FLAGS}>")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
  #  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
  #  set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
  endif()


  add_definitions("-DSPROFILE")
  add_definitions("-DNO_INTRA_PROFILE")
endif()


#######   GASNET AND UPCXX OPTIONS ######

option(ENABLE_ARIES "Enable ARIES" OFF)

include(cmake/gasnet.cmake)
include(cmake/upcxx.cmake)

include_directories(${UPCXX_INCLUDE_DIR})
link_directories(${UPCXX_LIBRARY_PATH})

#list( APPEND CMAKE_CXX_FLAGS -std=c++11)

#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${UPCXX_DEFINES})
#add_definitions("-std=c++11")
#add_definitions(${UPCXX_DEFINES})

if(GASNET_DEFINES)
  string(REPLACE " " ";" GASNET_DEFINES ${GASNET_DEFINES})
  message( STATUS "GASNET_DEFINES: " ${GASNET_DEFINES} )
  foreach( d ${GASNET_DEFINES} )
    string(STRIP ${d} d)
    message( STATUS "Found GASNET Define: " ${d} )
    add_definitions("${d}")
  endforeach()
endif()

if(UPCXX_DEFINES)
  string(REPLACE " " ";" UPCXX_DEFINES ${UPCXX_DEFINES})
  message( STATUS "UPCXX_DEFINES: " ${UPCXX_DEFINES} )
  foreach( d ${UPCXX_DEFINES} )
    string(STRIP ${d} d)
    message( STATUS "Found UPCXX Define: " ${d} )
    add_definitions("${d}")
  endforeach()
endif()


set(SYMPACKSRC src/ETree.cpp src/utility.cpp src/timer.cpp src/blas.cpp
  src/lapack.cpp src/global.cpp src/LogFile.cpp 
  src/IntervalTree.cpp src/Ordering.cpp src/CommPull.cpp
  src/DistSparseMatrixGraph.cpp src/amdbar.F src/genmmd.F src/mmdelm.F
  src/mmdint.F src/mmdnum.F src/mmdupd.F src/ordmmd.F
  src/betree.f src/del_heap.f src/ins_heap.f src/mcssup.f
  src/move_down.f src/ordsup.f src/boxnd.f src/gridnd.f
  src/degree.F src/fnroot.F src/genrcm.F src/rcm.F src/rootls.F
  src/invinv.f src/mod_heap.f src/move_up.f)

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_BINARY_DIR})

add_definitions("-DITREE")
add_definitions("-D_DSW_ITREE_")
add_definitions("-DNOTRY")
get_directory_property( ALL_DEFINITIONS DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )


###############################################
##  generate a config file
###############################################
file(WRITE "sympack_definitions.hpp.cmake" "#ifndef _HEADER_SYMPACK_DEFINITIONS_\n#define _HEADER_SYMPACK_DEFINITIONS_\n" )
#file(APPEND "sympack_definitions.hpp.cmake" "#cmakedefine _CONFIG_HEADER_H_\n" )
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

# generate a makefile include file

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

list(GET UPCXX_INCLUDE_DIR 0 TMPUPCXXDIR)
file(APPEND "sympack.mak.cmake" "\ninclude ${TMPUPCXXDIR}/upcxx.mak\n")
file(APPEND "sympack.mak.cmake" "VPATH=${TMPUPCXXDIR}\n\n")
file(APPEND "sympack.mak.cmake" "SYMPACK_INCLUDE += \$(UPCXX_CXXFLAGS)\n\n")
file(APPEND "sympack.mak.cmake" "SYMPACK_LIB = -L\$(SYMPACK_DIR)/lib -lsympack\n")

foreach( d ${LIBS} )
  file(APPEND "sympack.mak.cmake" "SYMPACK_LIB += ${d}\n")
endforeach()

#foreach( d ${GASNET_LIBRARIES} )
#  file(APPEND "sympack.mak.cmake" "SYMPACK_LIB += ${d}\n")
#endforeach()

file(APPEND "sympack.mak.cmake" "SYMPACK_LIB += \$(UPCXX_LDLIBS)\n")
file(APPEND "sympack.mak.cmake" "SYMPACK_LIB += -l${CMAKE_DL_LIBS}\n")

unset(MAKLIBS)
unset(TMPUPCXXDIR)

#if(ENABLE_SCOTCH)
#  file(APPEND "sympack.mak" "SYMPACK_LIB += ${SCOTCH_LIBRARIES}\n")
#endif()
#if(ENABLE_METIS)
#  file(APPEND "sympack.mak" "SYMPACK_LIB += ${METIS_LIBRARIES}\n")
#endif()
#if(ENABLE_PARMETIS)
#  file(APPEND "sympack.mak" "SYMPACK_LIB += ${PARMETIS_LIBRARIES}\n")
#endif()

CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/sympack.mak.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/sympack.mak
)


message(STATUS "ALL DEFINES: " ${ALL_DEFINITIONS})
 
add_library(sympack ${SYMPACKSRC})
add_dependencies(sympack gasnet upcxx)
 
if(UNIX AND NOT APPLE)
    target_link_libraries(sympack rt)
endif()
target_link_libraries(sympack libupcxx libgasnet-conduit ${LIBS}  ${GASNET_LIBRARIES} ${CMAKE_DL_LIBS})

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
  include/sympack/SuperNode_impl.hpp
  include/sympack/SuperNodeInd_impl.hpp
  include/sympack/CommTypes.hpp
  include/sympack/DistSparseMatrix.hpp
  include/sympack/symPACKMatrix.hpp
  include/sympack/DistSparseMatrixGraph.hpp
  include/sympack/symPACKMatrix_impl.hpp
  include/sympack/DistSparseMatrix_impl.hpp
  include/sympack/symPACKMatrix_impl_FB_pull.hpp
  include/sympack/ETree.hpp
  #include/sympack/SupernodalMatrixBase.hpp
  include/sympack/Environment.hpp
  include/sympack/Environment_impl.hpp 
  include/sympack/IntervalTree.hpp
  include/sympack/IntervalTree_impl.hpp 
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
  DESTINATION include/sympack)



add_executable(run_sympack EXCLUDE_FROM_ALL driver/run_sparse_pull.cpp)
add_custom_target(examples DEPENDS run_sympack)

add_dependencies(run_sympack sympack gasnet upcxx) 
#target_link_libraries(run_sympack sympack libupcxx libgasnet-conduit ${LIBS}  ${GASNET_LIBRARIES} ${CMAKE_DL_LIBS})
target_link_libraries(run_sympack sympack libupcxx libgasnet-conduit)

#install(TARGETS run_sympack RUNTIME DESTINATION bin)

set(SYMPACK_VERSION_MAJOR 0)
set(SYMPACK_VERSION_MINOR 1)
set(SYMPACK_VERSION_PATCH 0)
set(SYMPACK_VERSION_FULL "${SYMPACK_VERSION_MAJOR}.${SYMPACK_VERSION_MINOR}")
set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION_MAJOR ${SYMPACK_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${SYMPACK_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${SYMPACK_VERSION_PATCH})
#SET(CPACK_PACKAGE_FILE_NAME "sympack-"${SYMPACK_VERSION_FULL})
include(CPack)


