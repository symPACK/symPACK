cmake_minimum_required( VERSION 3.11 ) # Require CMake 3.11+
# Set up some auxillary vars if hints have been set



if( DEFINED ENV{UPCXX_INSTALL} )
  find_program( UPCXX_EXECUTABLE upcxx-meta HINTS $ENV{UPCXX_INSTALL}/bin NO_DEFAULT_PATH )
else()
  find_program( UPCXX_EXECUTABLE upcxx-meta )
endif()


if( UPCXX_EXECUTABLE )
  execute_process( COMMAND ${UPCXX_EXECUTABLE} CPPFLAGS OUTPUT_VARIABLE UPCXX_PPFLAGS)
  #execute_process( COMMAND ${UPCXX_EXECUTABLE} LDFLAGS OUTPUT_VARIABLE UPCXX_LDFLAGS)
  execute_process( COMMAND ${UPCXX_EXECUTABLE} LIBS OUTPUT_VARIABLE UPCXX_LIBFLAGS)
  #execute_process( COMMAND ${UPCXX_EXECUTABLE} CXXFLAGS OUTPUT_VARIABLE UPCXX_CXX_FLAGS)

  #message( STATUS "UPCXX_PPFLAGS: " ${UPCXX_PPFLAGS} )
  #message( STATUS "UPCXX_LDFLAGS: " ${UPCXX_LDFLAGS} )
  #message( STATUS "UPCXX_LIBFLAGS: " ${UPCXX_LIBFLAGS} )
  #message( STATUS "UPCXX_CXX_FLAGS: " ${UPCXX_CXX_FLAGS} )

  string(REPLACE "\n" " " UPCXX_LIBFLAGS ${UPCXX_LIBFLAGS})
  string(REPLACE "\n" " " UPCXX_PPFLAGS ${UPCXX_PPFLAGS})
  #string(REPLACE "\n" " " UPCXX_LDFLAGS ${UPCXX_LDFLAGS})
  #string(REPLACE "\n" " " UPCXX_CXX_FLAGS ${UPCXX_CXX_FLAGS})

  string(STRIP ${UPCXX_LIBFLAGS} UPCXX_LIBFLAGS)
  string(STRIP ${UPCXX_PPFLAGS}  UPCXX_PPFLAGS)
  #string(STRIP ${UPCXX_LDFLAGS}  UPCXX_LDFLAGS)
  #string(STRIP ${UPCXX_CXX_FLAGS}  UPCXX_CXX_FLAGS)

  #list( APPEND UPCXX_LIBRARIES ${UPCXX_LDFLAGS})
  list( APPEND UPCXX_LIBRARIES ${UPCXX_LIBFLAGS})

  #now separate include dirs from flags
  if(UPCXX_PPFLAGS)
    string(REPLACE " " ";" UPCXX_PPFLAGS ${UPCXX_PPFLAGS})
    foreach( option ${UPCXX_PPFLAGS} )
      string(STRIP ${option} option)
      string(REGEX MATCH "^-I" UPCXX_INCLUDE ${option})
      if( UPCXX_INCLUDE )
        string( REGEX REPLACE "^-I" "" option ${option} )
        list( APPEND UPCXX_INCLUDE_DIR ${option})
      else()
        string(REGEX MATCH "^-D" UPCXX_DEFINE ${option})
        if( UPCXX_DEFINE )
          string( REGEX REPLACE "^-D" "" option ${option} )
          list( APPEND UPCXX_DEFINITIONS ${option})
          #else()
          #list( APPEND UPCXX_CXX_FLAGS ${option})
        endif()
      endif()
    endforeach()
  endif()

  unset( UPCXX_LIBFLAGS )
  unset( UPCXX_PPFLAGS )
  #unset( UPCXX_LDFLAGS )
  unset( UPCXX_INCLUDE )
  unset( UPCXX_DEFINE )


endif()




foreach( dir ${UPCXX_UPCXX_INCLUDE_DIR} )

  if( EXISTS ${dir}/upcxx/upcxx.h )
  set( version_pattern 
    "^#define[\t ]+UPCXX_VERSION[\t ]+([0-9]+)$"
  )
  file( STRINGS ${dir}/upcxx/upcxx.h upcxx_version
        REGEX ${version_pattern} )
  
  foreach( match ${upcxx_version} )
    #string(REGEX REPLACE ${version_pattern} 
    #  "${UPCXX_VERSION_STRING}\\1" UPCXX_VERSION_STRING ${match}
    #)
    set(UPCXX_VERSION_STRING ${CMAKE_MATCH_2})
  endforeach()
  
  unset( upcxx_version )
  unset( version_pattern )

  endif()

endforeach()

if(UPCXX_VERSION_STRING)
  message( STATUS "UPCXX VERSION: " ${UPCXX_VERSION_STRING} )
endif()

# Determine if we've found UPCXX
mark_as_advanced( UPCXX_FOUND UPCXX_EXECUTABLE UPCXX_INCLUDE_DIR UPCXX_LIBRARIES UPCXX_DEFINITIONS ) #UPCXX_CXX_FLAGS

#message( STATUS "UPCXX_LIBRARIES: " ${UPCXX_LIBRARIES} )
#message( STATUS "UPCXX_INCLUDE_DIR: " ${UPCXX_INCLUDE_DIR} )
#message( STATUS "UPCXX_DEFINITIONS: " ${UPCXX_DEFINITIONS} )
#message( STATUS "UPCXX_CXX_FLAGS: " ${UPCXX_CXX_FLAGS} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( UPCXX
  REQUIRED_VARS UPCXX_EXECUTABLE UPCXX_LIBRARIES UPCXX_INCLUDE_DIR UPCXX_DEFINITIONS #UPCXX_CXX_FLAGS
  VERSION_VAR UPCXX_VERSION_STRING
  HANDLE_COMPONENTS
)


# Export target
if( UPCXX_FOUND AND NOT TARGET UPCXX::upcxx )

  add_library( UPCXX::upcxx INTERFACE IMPORTED )
  set_target_properties( UPCXX::upcxx PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${UPCXX_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES      "${UPCXX_LIBRARIES}" 
    INTERFACE_COMPILE_DEFINITIONS "${UPCXX_DEFINITIONS}" 
    INTERFACE_COMPILE_FEATURES    cxx_std_11
    #INTERFACE_COMPILE_OPTIONS     "${UPCXX_CXX_FLAGS}" 
    #INTERFACE_LINK_OPTIONS      "${UPCXX_LDFLAGS}" 
  )
endif()
