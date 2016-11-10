#BLAS and Lapack if not in standard path
set(BLAS_DIR /path/to/blas)
set(LAPACK_DIR /path/to/lapack)
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${BLAS_DIR})
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${LAPACK_DIR})

#optional Ordering libraries
#set(SCOTCH_DIR /path/to/scotch)
#set(ENABLE_SCOTCH ON CACHE BOOL "...")

#set(METIS_DIR /path/to/metis)
#set(ENABLE_METIS ON CACHE BOOL "...")

#set(PARMETIS_DIR /path/to/parmetis)
#set(ENABLE_PARMETIS ON CACHE BOOL "...")
