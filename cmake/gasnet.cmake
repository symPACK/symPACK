set(GASNET_NAME gasnet)
set(GASNET_URL https://gasnet.lbl.gov/)
set(GASNET_GZ  GASNet-1.26.3.tar.gz)
set(GASNET_MD5 "93ab985881ca19de319bb73583c385d2")

ExternalProject_Add(${GASNET_NAME}
 URL ${GASNET_URL}/${GASNET_GZ}
 URL_MD5 ${GASNET_MD5}
 INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install
 CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> CFLAGS=-DGASNETI_PSHM_BARRIER_HIER=0
)


set(GASNET_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install/include)
#TODO handle others conduit

set(GASNET_LIBRARY_PATH ${CMAKE_CURRENT_BINARY_DIR}/external/gasnet_install/lib)


#if (APPLE)
#  SET(GASNET_LIBRARIES -gasnet-mpi-seq -ammpi ${MPI_CXX_LIBRARIES} CACHE INTERNAL "gasnet libs")
#else ()
#  SET(GASNET_LIBRARIES -gasnet-mpi-seq -ammpi ${MPI_CXX_LIBRARIES} -pthread -dl -hugetlbfs CACHE INTERNAL "gasnet libs")
#endif ()

ExternalProject_Get_Property(${GASNET_NAME} install_dir)
add_library(libgasnet-conduit STATIC IMPORTED)

set(GASNET_CONDUIT_INCLUDE_DIR ${GASNET_INCLUDE_DIR}/mpi-conduit)
set(GASNET_CONDUIT ${GASNET_CONDUIT_INCLUDE_DIR}/mpi-seq.mak)
set_property(TARGET libgasnet-conduit PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libgasnet-mpi-seq.a)
SET(GASNET_LIBRARIES -lammpi ${MPI_CXX_LIBRARIES} -lpthread -lrt )
add_dependencies(libgasnet-conduit ${GASNET_NAME})
set(GASNET_DEFINES "-DGASNET_SEQ -DGASNETT_USE_GCC_ATTRIBUTE_MAYALIAS=1")

#set(GASNET_CONDUIT_INCLUDE_DIR ${GASNET_INCLUDE_DIR}/ibv-conduit)
#set(GASNET_CONDUIT ${GASNET_CONDUIT_INCLUDE_DIR}/ibv-seq.mak)
#set_property(TARGET libgasnet-conduit PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libgasnet-ibv-seq.a)
#add_dependencies(libgasnet-conduit ${GASNET_NAME})
#set(GASNET_DEFINES "-DGASNET_SEQ -DGASNETT_USE_GCC_ATTRIBUTE_MAYALIAS=1 -D_REENTRANT -D_GNU_SOURCE")
#SET(GASNET_LIBRARIES -lammpi -libverbs ${MPI_CXX_LIBRARIES} -lpthread -lrt )
