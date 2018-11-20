set(GASNET_PREFIX $ENV{GASNET_DIR})
if (GASNET_PREFIX)
    message("GASNET has been provided, not compiling local version.")
    set(GASNET_INCLUDE_DIR ${GASNET_PREFIX}/include)
    set(GASNET_LIBRARY_PATH ${GASNET_PREFIX}/lib)

    add_library(libgasnet-conduit STATIC IMPORTED)

    #make sure GASNET_CONDUIT is in the environment
    set(GASNET_CONDUIT $ENV{GASNET_CONDUIT})
    get_filename_component(GASNET_CONDUIT_INCLUDE_DIR ${GASNET_CONDUIT} DIRECTORY)

    set(GASNET_DEFINES "")
    set(GASNET_DEFINES $ENV{GASNET_DEFINES})
    #set(GASNET_DEFINES ${GASNET_DEFINES} "-DGASNET_CONDUIT_ARIES -DGASNET_SEQ -DUSE_GASNET_FAST_SEGMENT -DONLY_MSPACES -DGASNET_ALLOW_OPTIMIZED_DEBUG")


    get_filename_component(GASNET_CONDUIT_NAME ${GASNET_CONDUIT} NAME_WE)

    set_property(TARGET libgasnet-conduit PROPERTY IMPORTED_LOCATION ${GASNET_PREFIX}/lib/libgasnet-${GASNET_CONDUIT_NAME}.a)

    set(GASNET_LIBRARIES $ENV{GASNET_LIBRARIES} ${MPI_CXX_LIBRARIES})

    add_library(gasnet STATIC IMPORTED)
else()

  set(GASNET_NAME gasnet)
  set(GASNET_URL https://gasnet.lbl.gov/download/)
  set(GASNET_GZ  GASNet-1.30.0.tar.gz)
  #set(GASNET_MD5 "4c6255469c9d2922ad8e2c8ab757eff1")
  set(GASNET_CFLAGS "${CMAKE_C_FLAGS} -DGASNETI_PSHM_BARRIER_HIER=0")
 
  if(ENABLE_ARIES)
    ExternalProject_Add(${GASNET_NAME}
    URL ${GASNET_URL}/${GASNET_GZ}
    URL_MD5 ${GASNET_MD5}
    #URL ${PROJECT_SOURCE_DIR}/tarballs/${GASNET_GZ}
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --enable-pshm-hugetlbfs --enable-pshm-xpmem --enable-pshm --disable-pshm-posix --enable-aries --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} MPI_CC=${MPI_C_COMPILER} CFLAGS=${GASNET_CFLAGS} CXXFLAGS=${CMAKE_CXX_FLAGS}

        )
  
  
    set(GASNET_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install/include)
    set(GASNET_LIBRARY_PATH ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install/lib)
    ExternalProject_Get_Property(${GASNET_NAME} install_dir)
    add_library(libgasnet-conduit STATIC IMPORTED)
    set(GASNET_DEFINES "")
    
    
    set(GASNET_CONDUIT_INCLUDE_DIR ${GASNET_INCLUDE_DIR}/aries-conduit)
    set(GASNET_CONDUIT ${GASNET_CONDUIT_INCLUDE_DIR}/aries-seq.mak)
    set(GASNET_DEFINES ${GASNET_DEFINES} "-DGASNET_CONDUIT_ARIES -DGASNET_SEQ -DUSE_GASNET_FAST_SEGMENT -DONLY_MSPACES -DGASNET_ALLOW_OPTIMIZED_DEBUG")
    set_property(TARGET libgasnet-conduit PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libgasnet-aries-seq.a)
    
    if(APPLE)
      SET(GASNET_LIBRARIES -lammpi ${MPI_CXX_LIBRARIES} -lpthread )
    else()
      SET(GASNET_LIBRARIES -lammpi ${MPI_CXX_LIBRARIES} -lpthread -lrt )
    endif()
    
    add_dependencies(libgasnet-conduit ${GASNET_NAME})
  
  else()
    ExternalProject_Add(${GASNET_NAME}
        URL ${GASNET_URL}/${GASNET_GZ}
        URL_MD5 ${GASNET_MD5}
  #      URL ${PROJECT_SOURCE_DIR}/tarballs/${GASNET_GZ}
        INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --disable-aligned-segments --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} MPI_CC=${MPI_C_COMPILER} CFLAGS=${GASNET_CFLAGS} CXXFLAGS=${CMAKE_CXX_FLAGS}
        )
  
    set(GASNET_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install/include)
    set(GASNET_LIBRARY_PATH ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install/lib)
    
    ExternalProject_Get_Property(${GASNET_NAME} install_dir)
    add_library(libgasnet-conduit STATIC IMPORTED)
    set(GASNET_DEFINES "")
  
    set(GASNET_CONDUIT_INCLUDE_DIR ${GASNET_INCLUDE_DIR}/mpi-conduit)
    set(GASNET_CONDUIT ${GASNET_CONDUIT_INCLUDE_DIR}/mpi-seq.mak)
    set_property(TARGET libgasnet-conduit PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libgasnet-mpi-seq.a)
    set(GASNET_DEFINES ${GASNET_DEFINES} "-DGASNET_SEQ -DGASNETT_USE_GCC_ATTRIBUTE_MAYALIAS=1")
  
    if(APPLE)
      SET(GASNET_LIBRARIES -lammpi ${MPI_CXX_LIBRARIES} -lpthread )
    else()
      SET(GASNET_LIBRARIES -lammpi ${MPI_CXX_LIBRARIES} -lpthread -lrt )
    endif()
  
    add_dependencies(libgasnet-conduit ${GASNET_NAME})
  endif()

endif()


#set(GASNET_CONDUIT_INCLUDE_DIR ${GASNET_INCLUDE_DIR}/ibv-conduit)
#set(GASNET_CONDUIT ${GASNET_CONDUIT_INCLUDE_DIR}/ibv-seq.mak)
#set_property(TARGET libgasnet-conduit PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libgasnet-ibv-seq.a)
#add_dependencies(libgasnet-conduit ${GASNET_NAME})
#set(GASNET_DEFINES "-DGASNET_SEQ -DGASNETT_USE_GCC_ATTRIBUTE_MAYALIAS=1 -D_REENTRANT -D_GNU_SOURCE")
#SET(GASNET_LIBRARIES -lammpi -libverbs ${MPI_CXX_LIBRARIES} -lpthread -lrt )
