cmake_minimum_required( VERSION 3.10 ) # Require CMake 3.10+

get_filename_component(SYMPACK_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(CMakeFindDependencyMacro)

list(APPEND CMAKE_MODULE_PATH ${SYMPACK_CMAKE_DIR})
list(APPEND CMAKE_MODULE_PATH ${SYMPACK_CMAKE_DIR}/Modules)

include( symPACKDepends ) # Basic Dependencies

list(REMOVE_AT CMAKE_MODULE_PATH -1) 

if(NOT TARGET symPACK::sympack)
  include("${SYMPACK_CMAKE_DIR}/symPACKTargets.cmake")
endif()

set(SYMPACK_LIBRARIES symPACK::sympack)

