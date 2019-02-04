add_library( symPACK::ordering     INTERFACE IMPORTED )
add_library( symPACK::parallel_cxx     INTERFACE IMPORTED )
add_library( symPACK::parallel_c       INTERFACE IMPORTED )
add_library( symPACK::parallel_fortran INTERFACE IMPORTED )


# Handle MPI and UPCXX 
find_package( MPI )
find_package( UPCXX REQUIRED )
#######   ORDERING LIBRARIES ######

if(ENABLE_METIS)
  if (NOT METIS_FOUND)
    find_package(METIS REQUIRED)
    target_link_libraries(symPACK::ordering INTERFACE METIS::metis)
  endif()
  target_compile_definitions( symPACK::ordering INTERFACE "-DUSE_METIS")
endif()

if(ENABLE_PARMETIS)
  find_package(ParMETIS REQUIRED)
  target_link_libraries(symPACK::ordering INTERFACE ParMETIS::parmetis)
  target_compile_definitions( symPACK::ordering INTERFACE "-DUSE_PARMETIS")
  set(ENABLE_METIS ON)
endif()


if(ENABLE_SCOTCH)
  if (NOT SCOTCH_FOUND)
    find_package(SCOTCH REQUIRED)
    target_link_libraries(symPACK::ordering INTERFACE SCOTCH::scotch)
  endif()
  target_compile_definitions( symPACK::ordering INTERFACE "-DUSE_SCOTCH")
endif()


if(ENABLE_PTSCOTCH)
  find_package(PTSCOTCH REQUIRED)
  target_link_libraries(symPACK::ordering INTERFACE PTSCOTCH::ptscotch)
  target_compile_definitions( symPACK::ordering INTERFACE "-DUSE_PTSCOTCH")
  set(ENABLE_SCOTCH ON)
endif()



target_link_libraries( symPACK::parallel_cxx     INTERFACE MPI::MPI_CXX  UPCXX::upcxx )   
target_link_libraries( symPACK::parallel_c       INTERFACE MPI::MPI_C       )
target_link_libraries( symPACK::parallel_fortran INTERFACE MPI::MPI_Fortran )

# Handle OpenMP
if( symPACK_ENABLE_OPENMP )

  find_package( OpenMP REQUIRED )

  target_link_libraries( symPACK::parallel_cxx     INTERFACE OpenMP::OpenMP_CXX     )   
  target_link_libraries( symPACK::parallel_c       INTERFACE OpenMP::OpenMP_C       )
  target_link_libraries( symPACK::parallel_fortran INTERFACE OpenMP::OpenMP_Fortran )

endif()

