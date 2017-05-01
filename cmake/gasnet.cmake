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


    get_filename_component(GASNET_CONDUIT_NAME ${GASNET_CONDUIT} NAME_WE)

    set_property(TARGET libgasnet-conduit PROPERTY IMPORTED_LOCATION ${GASNET_PREFIX}/lib/libgasnet-${GASNET_CONDUIT_NAME}.a)

    set(GASNET_LIBRARIES $ENV{GASNET_LIBRARIES} ${MPI_CXX_LIBRARIES})

    add_library(gasnet STATIC IMPORTED)
else()

  set(GASNET_NAME gasnet)
  set(GASNET_URL https://gasnet.lbl.gov/)
  set(GASNET_GZ  GASNet-1.28.2.tar.gz)
  set(GASNET_MD5 "6ca0463dc2430570e40646c4d1e97b36")
  set(GASNET_CFLAGS "${CMAKE_C_FLAGS} -DGASNETI_PSHM_BARRIER_HIER=0")
 
  if(ENABLE_ARIES)
    if(ENABLE_KNL)
    if(ENABLE_KNL_ONLY)
ExternalProject_Add(${GASNET_NAME}
    URL ${GASNET_URL}/${GASNET_GZ}
    URL_MD5 ${GASNET_MD5}
    #URL ${PROJECT_SOURCE_DIR}/tarballs/${GASNET_GZ}
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install
    #PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gasnet_prefix
    #INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/gasnet_install
    #INSTALL_COMMAND ""
    CONFIGURE_COMMAND <SOURCE_DIR>/cross-configure-intel-knl-gasnet-knl-only --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} MPI_CC=${MPI_C_COMPILER} CFLAGS=${GASNET_CFLAGS} CXXFLAGS=${CMAKE_CXX_FLAGS}
        )

  ExternalProject_Add_Step(${GASNET_NAME} cross-compile
      DEPENDEES patch update patch download
      DEPENDERS configure
      COMMAND ln -s ${PROJECT_SOURCE_DIR}/contrib/config/cross-configure-intel-knl-gasnet-knl-only <SOURCE_DIR> 
      WORKING_DIRECTORY <SOURCE_DIR>
#ALWAYS 1
      )
    else()
    ExternalProject_Add(${GASNET_NAME}
    URL ${GASNET_URL}/${GASNET_GZ}
    URL_MD5 ${GASNET_MD5}
    #URL ${PROJECT_SOURCE_DIR}/tarballs/${GASNET_GZ}
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install
    #PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gasnet_prefix
    #INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/gasnet_install
    #INSTALL_COMMAND ""
    CONFIGURE_COMMAND <SOURCE_DIR>/cross-configure-intel-knl-gasnet --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} MPI_CC=${MPI_C_COMPILER} CFLAGS=${GASNET_CFLAGS} CXXFLAGS=${CMAKE_CXX_FLAGS}
        )

  ExternalProject_Add_Step(${GASNET_NAME} cross-compile
      DEPENDEES patch update patch download
      DEPENDERS configure
      COMMAND ln -s ${PROJECT_SOURCE_DIR}/contrib/config/cross-configure-intel-knl-gasnet <SOURCE_DIR> 
      WORKING_DIRECTORY <SOURCE_DIR>
#ALWAYS 1
      )
    endif()

    set(GASNET_DEFINES "")

    if(GASNET_PAR)
      set(GASNET_DEFINES ${GASNET_DEFINES} "-DGASNET_CONDUIT_ARIES -DGASNET_PAR -DGASNET_ALLOW_OPTIMIZED_DEBUG -D_REENTRANT")
      set(GASNET_LIBRARIES ${MPI_CXX_LIBRARIES} -L/opt/cray/pe/pmi/5.0.10-1.0000.11050.0.0.ari/lib64 -lpmi -L/opt/cray/ugni/6.0.12-2.1/lib64 -lugni -L/opt/cray/udreg/2.3.2-4.6/lib64 -ludreg -L/opt/cray/xpmem/0.1-4.5/lib64 -lxpmem -lhugetlbfs -lm )
    else()
      set(GASNET_DEFINES ${GASNET_DEFINES} "-DGASNET_CONDUIT_ARIES -DGASNET_SEQ -DGASNET_ALLOW_OPTIMIZED_DEBUG")
      set(GASNET_LIBRARIES ${MPI_CXX_LIBRARIES} -L/opt/cray/pe/pmi/5.0.10-1.0000.11050.0.0.ari/lib64 -lpmi -L/opt/cray/ugni/6.0.12-2.1/lib64 -lugni -L/opt/cray/udreg/2.3.2-4.6/lib64 -ludreg  -L/opt/cray/xpmem/0.1-4.5/lib64 -lxpmem -lhugetlbfs -lm )
    endif()





    else()
     ExternalProject_Add(${GASNET_NAME}
    URL ${GASNET_URL}/${GASNET_GZ}
    URL_MD5 ${GASNET_MD5}
    #URL ${PROJECT_SOURCE_DIR}/tarballs/${GASNET_GZ}
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install
    #PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gasnet_prefix
    #INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/gasnet_install
    #INSTALL_COMMAND ""
    CONFIGURE_COMMAND <SOURCE_DIR>/configure  --enable-aries --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} MPI_CC=${MPI_C_COMPILER} CFLAGS=${GASNET_CFLAGS} CXXFLAGS=${CMAKE_CXX_FLAGS}
    #CONFIGURE_COMMAND <SOURCE_DIR>/configure  --enable-pshm-hugetlbfs --enable-pshm-xpmem --enable-pshm --disable-pshm-posix --enable-aries --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} MPI_CC=${MPI_C_COMPILER} CFLAGS=${GASNET_CFLAGS} CXXFLAGS=${CMAKE_CXX_FLAGS}
        )

    if(GASNET_PAR)
      set(GASNET_DEFINES "")
      set(GASNET_DEFINES ${GASNET_DEFINES} "-DGASNET_CONDUIT_ARIES -DGASNET_PAR -DUSE_GASNET_FAST_SEGMENT -DONLY_MSPACES -DGASNET_ALLOW_OPTIMIZED_DEBUG -D_REENTRANT")
    else()
      set(GASNET_DEFINES "")
      set(GASNET_DEFINES ${GASNET_DEFINES} "-DGASNET_CONDUIT_ARIES -DGASNET_SEQ -DUSE_GASNET_FAST_SEGMENT -DONLY_MSPACES -DGASNET_ALLOW_OPTIMIZED_DEBUG")
    endif()
    
    if(APPLE)
      SET(GASNET_LIBRARIES -lammpi ${MPI_CXX_LIBRARIES} -lhugetlbfs -lpthread )
    else()
      SET(GASNET_LIBRARIES -lammpi ${MPI_CXX_LIBRARIES} -lhugetlbfs -lpthread -lrt )
    endif()

    endif()
  
  
    set(GASNET_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install/include)
    set(GASNET_LIBRARY_PATH ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install/lib)
    ExternalProject_Get_Property(${GASNET_NAME} install_dir)
    add_library(libgasnet-conduit STATIC IMPORTED)
    
    
    set(GASNET_CONDUIT_INCLUDE_DIR ${GASNET_INCLUDE_DIR}/aries-conduit)
    if(GASNET_PAR)
    set(GASNET_CONDUIT ${GASNET_CONDUIT_INCLUDE_DIR}/aries-par.mak)
    set_property(TARGET libgasnet-conduit PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libgasnet-aries-par.a)
    else()
    set(GASNET_CONDUIT ${GASNET_CONDUIT_INCLUDE_DIR}/aries-seq.mak)
    set_property(TARGET libgasnet-conduit PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libgasnet-aries-seq.a)
    endif()

    
    add_dependencies(libgasnet-conduit ${GASNET_NAME})
  
  else()
    ExternalProject_Add(${GASNET_NAME}
        URL ${GASNET_URL}/${GASNET_GZ}
        URL_MD5 ${GASNET_MD5}
  #      URL ${PROJECT_SOURCE_DIR}/tarballs/${GASNET_GZ}
        INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install
    #PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gasnet_prefix
    #    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/gasnet_install
    #    INSTALL_COMMAND ""
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --disable-aligned-segments --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} MPI_CC=${MPI_C_COMPILER} CFLAGS=${GASNET_CFLAGS} CXXFLAGS=${CMAKE_CXX_FLAGS}
        )
  
    set(GASNET_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install/include)
    set(GASNET_LIBRARY_PATH ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install/lib)
    
    ExternalProject_Get_Property(${GASNET_NAME} install_dir)
    add_library(libgasnet-conduit STATIC IMPORTED)
    set(GASNET_DEFINES "")
  
    set(GASNET_CONDUIT_INCLUDE_DIR ${GASNET_INCLUDE_DIR}/mpi-conduit)
    if(GASNET_PAR)
      set(GASNET_CONDUIT ${GASNET_CONDUIT_INCLUDE_DIR}/mpi-par.mak)
      set_property(TARGET libgasnet-conduit PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libgasnet-mpi-par.a)
      set(GASNET_DEFINES ${GASNET_DEFINES} "-DGASNET_PAR -D_REENTRANT -DGASNETT_USE_GCC_ATTRIBUTE_MAYALIAS=1")
    else()
      set(GASNET_CONDUIT ${GASNET_CONDUIT_INCLUDE_DIR}/mpi-seq.mak)
      set_property(TARGET libgasnet-conduit PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libgasnet-mpi-seq.a)
      set(GASNET_DEFINES ${GASNET_DEFINES} "-DGASNET_SEQ -DGASNETT_USE_GCC_ATTRIBUTE_MAYALIAS=1")
    endif()

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
