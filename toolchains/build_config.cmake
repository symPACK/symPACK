#BLAS and Lapack if not in standard path
set(BLAS_DIR /path/to/blas CACHE BOOL "...")
set(LAPACK_DIR /path/to/lapack CACHE BOOL "...")
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${BLAS_DIR})
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${LAPACK_DIR})

#If BLAS or LAPACK cannot be found, especially when using latest cmake versions,
#please uncomment the following lines and comment the lines above

#set(BLAS_LIBRARIES /path/to/libblas.a CACHE BOOL "...")
#set(LAPACK_LIBRARIES /path/to/liblapack.a CACHE BOOL "...")


#optional Ordering libraries
#set(SCOTCH_DIR /path/to/scotch)
#set(ENABLE_SCOTCH ON CACHE BOOL "...")

#set(METIS_DIR /path/to/metis)
#set(ENABLE_METIS ON CACHE BOOL "...")

#set(PARMETIS_DIR /path/to/parmetis)
#set(ENABLE_PARMETIS ON CACHE BOOL "...")
