
# Compile

############################## Source files ######################################

set(SYMPACK_SRC_CXX
ETree.cpp
utility.cpp
timer.cpp
blas.cpp
lapack.cpp
global.cpp
LogFile.cpp
IntervalTree.cpp
Ordering.cpp
CommPull.cpp
DistSparseMatrixGraph.cpp
cuBLAS.cpp
)

set(SYMPACK_SRC_Fortran
amdbar.F
genmmd.F
mmdelm.F
mmdint.F
mmdnum.F
mmdupd.F
ordmmd.F
betree.f
del_heap.f
ins_heap.f
mcssup.f
move_down.f
ordsup.f
boxnd.f
gridnd.f
degree.F
fnroot.F
genrcm.F
rcm.F
rootls.F
TSPB/compress_skel2.f
TSPB/fnskel2.f
TSPB/near_ins2.f
TSPB/ordsup_ind_tsp_paths2.f
TSPB/ordsup_tsp_paths2.f
TSPB/post_suptree2.f
TSPB/sup_etpost2.f
TSPB/fnskel.f
TSPB/near_ins.f
TSPB/ordsup_ind_tsp_paths.f
TSPB/ordsup_tsp_paths.f
TSPB/post_suptree.f
TSPB/sup_etpost.f
invinv.f
mod_heap.f
move_up.f
)


add_library(sympack ${SYMPACK_SRC_CXX} ${SYMPACK_SRC_Fortran})
add_library( symPACK::sympack ALIAS sympack )


include(FortranCInterface)
FortranCInterface_HEADER("${PROJECT_BINARY_DIR}/sympack_config.hpp" MACRO_NAMESPACE "FC_")
FortranCInterface_VERIFY(CXX)

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

# Handle OpenMP
if( ENABLE_OPENMP )
  find_package( OpenMP REQUIRED )
  target_link_libraries( symPACK::parallel_cxx     INTERFACE OpenMP::OpenMP_CXX     )   
  target_link_libraries( symPACK::parallel_c       INTERFACE OpenMP::OpenMP_C       )
  target_link_libraries( symPACK::parallel_fortran INTERFACE OpenMP::OpenMP_Fortran )
endif()

if(ENABLE_PROFILE)
  if(NOT ENABLE_OPENMP)
    find_package( OpenMP REQUIRED )
    target_link_libraries( symPACK::parallel_cxx     INTERFACE OpenMP::OpenMP_CXX     )   
    target_link_libraries( symPACK::parallel_c       INTERFACE OpenMP::OpenMP_C       )
    target_link_libraries( symPACK::parallel_fortran INTERFACE OpenMP::OpenMP_Fortran )
  endif()

  target_compile_definitions( sympack PUBLIC "SPROFILE")
  target_compile_definitions( sympack PUBLIC "NO_INTRA_PROFILE")
endif()


if(RCM_IDX_64)
  target_compile_definitions( sympack PUBLIC  "RCMIDXSIZE=64")
else()
  target_compile_definitions( sympack PUBLIC  "RCMIDXSIZE=32")
endif()

if(MMD_IDX_64)
  target_compile_definitions( sympack PUBLIC  "MMDIDXSIZE=64")
else()
  target_compile_definitions( sympack PUBLIC  "MMDIDXSIZE=32")
endif()

if(AMD_IDX_64)
  target_compile_definitions( sympack PUBLIC  "AMDIDXSIZE=64")
else()
  target_compile_definitions( sympack PUBLIC  "AMDIDXSIZE=32")
endif()

if(ENABLE_THREADS)
  get_target_property(UPCXX_DEFINITIONS UPCXX::upcxx INTERFACE_COMPILE_DEFINITIONS)
  if("GASNET_PAR" IN_LIST UPCXX_DEFINITIONS)
    target_compile_definitions( sympack PUBLIC  "SP_THREADS")
  else()
    message(FATAL_ERROR "UPCXX_THREADMODE=par is required when using ENABLE_THREADS=ON\n(GASNET_PAR must be in the compiler definitions)")
  endif()
endif()



target_compile_definitions(sympack PUBLIC "ITREE")
target_compile_definitions(sympack PUBLIC "_DSW_ITREE_")
  

target_link_libraries(sympack PUBLIC 
  symPACK::ordering
  symPACK::parallel_cxx
  symPACK::math
  ${LIBS} 
  ${CMAKE_DL_LIBS}
)

if(UNIX AND NOT APPLE)
  target_link_libraries(sympack PUBLIC rt)
endif()

############################ Compile Options ################################

target_compile_features(sympack PUBLIC cxx_std_11)
set_target_properties(sympack PROPERTIES POSITION_INDEPENDENT_CODE ON)



target_include_directories( sympack  
  PUBLIC
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  )


###############################################
##  generate a config file
###############################################
get_target_property( ALL_DEFINITIONS sympack COMPILE_DEFINITIONS )
file(WRITE "${CMAKE_BINARY_DIR}/sympack_definitions.hpp.cmake" "#ifndef _HEADER_SYMPACK_DEFINITIONS_\n#define _HEADER_SYMPACK_DEFINITIONS_\n" )
foreach( d ${ALL_DEFINITIONS} )
  string(REPLACE "=" " " d ${d})
  string(REPLACE " " ";" DEF_LIST ${d})
  list(GET DEF_LIST 0 def_name)

  file(APPEND "${CMAKE_BINARY_DIR}/sympack_definitions.hpp.cmake" "\n#ifndef ${def_name}\n" )
  file(APPEND "${CMAKE_BINARY_DIR}/sympack_definitions.hpp.cmake" "#define ${d}\n" )
  file(APPEND "${CMAKE_BINARY_DIR}/sympack_definitions.hpp.cmake" "#endif // ${def_name}\n" )
endforeach()
file(APPEND "${CMAKE_BINARY_DIR}/sympack_definitions.hpp.cmake" "\n#endif //_HEADER_SYMPACK_DEFINITIONS_\n" )

CONFIGURE_FILE(
    ${CMAKE_BINARY_DIR}/sympack_definitions.hpp.cmake
    ${CMAKE_BINARY_DIR}/sympack_definitions.hpp
)

###############################################
## generate a makefile include file
###############################################
message(STATUS "install path " ${CMAKE_INSTALL_PREFIX})

function( get_target_prop_internal targetname propname top)
  #message(STATUS "INPUT ${${top}_libraries}")
  if ( TARGET ${targetname} )
    get_target_property( tgt_LIBS ${targetname} ${propname} )
    foreach( d ${tgt_LIBS} )
      get_target_prop_internal(${d} ${propname} ${targetname})  
    endforeach()
  elseif( targetname STREQUAL "dl" )
    set( ${targetname}_libraries "-ldl")
  elseif( targetname STREQUAL "rt" )
    set( ${targetname}_libraries "-lrt")
  else()
    set( ${targetname}_libraries "${targetname}")
  endif()
  #message(STATUS "OUTPUT from within ${${targetname}_libraries}")

  if (${top}_libraries) 
    set( ${top}_libraries "${${top}_libraries}" "${${targetname}_libraries}"  PARENT_SCOPE)
  else()
    set( ${top}_libraries "${${targetname}_libraries}"  PARENT_SCOPE)
  endif()
endfunction(get_target_prop_internal)

function( get_target_prop varname targetname propname)
   get_target_prop_internal( ${targetname} ${propname} ${targetname})
   #message(STATUS "OUTPUT ${${targetname}_libraries}")
   list(REVERSE ${targetname}_libraries)
   list(REMOVE_DUPLICATES ${targetname}_libraries)
   list(REVERSE ${targetname}_libraries)
   set( ${varname} "${${targetname}_libraries}"  PARENT_SCOPE) 
endfunction(get_target_prop)

get_target_prop(sympack_LIBRARIES sympack INTERFACE_LINK_LIBRARIES)

string(REPLACE ";" " " sympack_LIBRARIES "${sympack_LIBRARIES}")
message(STATUS "SYMPACK LIBRARIES ${sympack_LIBRARIES}")
#doesnt work
#get_target_prop(sympack_INCDIRS sympack INTERFACE_INCLUDE_DIRECTORIES)
#message(STATUS "SYMPACK INCDIRS ${sympack_INCDIRS}")


# Deal with @...@ configurable items:
configure_file(
  ${CMAKE_SOURCE_DIR}/cmake/sympack.mak.in
  ${CMAKE_BINARY_DIR}/sympack.mak_at_configured
  @ONLY
  )

# Deal with items configured with generator expressions:
file(GENERATE
  OUTPUT
  ${CMAKE_BINARY_DIR}/sympack.mak_at_generated
  INPUT
  ${CMAKE_BINARY_DIR}/sympack.mak_at_configured
  )

# Build driver
add_executable(run_sympack2D ../driver/run_sympack2D.cpp)
target_include_directories( run_sympack2D PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/driver> )
target_link_libraries(run_sympack2D PUBLIC symPACK::sympack)
if (ENABLE_CUDA)
    target_link_libraries(run_sympack2D PUBLIC -lcublas -lcusolver)
endif()
set_target_properties(run_sympack2D PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)

#### INSTALL rules
###
###include( GNUInstallDirs )
###set( INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/symPACK )
###
#### Targets
###install(TARGETS sympack
###  EXPORT sympack-targets
###  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
###  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
###  )
###
###set_target_properties( sympack PROPERTIES EXPORT_NAME symPACK )
###
#### Install Headers
###install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )
###
#### Export target to script
###install(EXPORT sympack-targets
###  FILE
###  symPACKTargets.cmake
###  NAMESPACE
###  symPACK::
###  DESTINATION
###  ${INSTALL_CONFIGDIR}
###  )
###
###
####Create a ConfigVersion.cmake file
###include(CMakePackageConfigHelpers)
###write_basic_package_version_file(
###  ${CMAKE_CURRENT_BINARY_DIR}/symPACKConfigVersion.cmake
###  VERSION ${PROJECT_VERSION}
###  COMPATIBILITY AnyNewerVersion
###  )
###
#### Setup PEXSIConfig.cmake
###configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/symPACKConfig.cmake.in
###  ${CMAKE_CURRENT_BINARY_DIR}/symPACKConfig.cmake
###  INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
###  )
###
####Install the config, configversion and custom find modules
###install(DIRECTORY 
###  ${PROJECT_SOURCE_DIR}/cmake/ 
###  DESTINATION ${INSTALL_CONFIGDIR} 
###  FILES_MATCHING PATTERN "*.cmake"
###  )
###
###install(FILES
###  ${CMAKE_CURRENT_BINARY_DIR}/symPACKConfig.cmake
###  ${CMAKE_CURRENT_BINARY_DIR}/symPACKConfigVersion.cmake
###  DESTINATION ${INSTALL_CONFIGDIR}
###  )

###########################################################################

include( GNUInstallDirs )
set( INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/symPACK )

# Targets
install(TARGETS sympack
  EXPORT sympack-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

set_target_properties( sympack PROPERTIES EXPORT_NAME symPACK )

# Install Headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )

install(FILES 
  ${PROJECT_BINARY_DIR}/sympack_config.hpp 
  ${PROJECT_BINARY_DIR}/sympack_definitions.hpp
  #${PROJECT_BINARY_DIR}/sympack.mak
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})


# Export target to script
install(EXPORT sympack-targets
  FILE
  symPACKTargets.cmake
  NAMESPACE
  symPACK::
  DESTINATION
  ${INSTALL_CONFIGDIR}
  )

#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/symPACKConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
  )

# Setup symPACKConfig.cmake
configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/symPACKConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/symPACKConfig.cmake
  INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
  )

#Install the config, configversion and custom find modules
install(DIRECTORY 
  ${PROJECT_SOURCE_DIR}/cmake/ 
  DESTINATION ${INSTALL_CONFIGDIR} 
  FILES_MATCHING PATTERN "*.cmake"
  )

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/symPACKConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/symPACKConfigVersion.cmake
  DESTINATION ${INSTALL_CONFIGDIR}
  )

# Deal with items configured with generator expressions:
#install(CODE "MESSAGE(\"Sample install message.\")")
#install( CODE "FILE(GENERATE OUTPUT \${CMAKE_INSTALL_INCLUDEDIR}/sympack.mak INPUT \${CMAKE_BINARY_DIR}/sympack.mak_at_configured )" )
configure_file(
  ${CMAKE_SOURCE_DIR}/cmake/install_mak.cmake.in
  ${CMAKE_BINARY_DIR}/install_mak.cmake
  @ONLY
  )


install(SCRIPT ${CMAKE_BINARY_DIR}/install_mak.cmake )



