#UPCXX if not in standard path
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${UPCXX_INSTALL})

#BLAS and Lapack if not in standard path
set(BLAS_DIR /path/to/blas)
set(LAPACK_DIR /path/to/lapack)
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${BLAS_DIR})
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${LAPACK_DIR})

#optional Ordering libraries
#set(scotch_PREFIX /path/to/scotch)
#set(ENABLE_SCOTCH ON CACHE BOOL "...")

#set(ptscotch_PREFIX /path/to/ptscotch)
#set(ENABLE_PTSCOTCH ON CACHE BOOL "...")

#setmetis_PREFIX /path/to/metis)
#set(ENABLE_METIS ON CACHE BOOL "...")

#set(parmetis_PREFIX /path/to/parmetis)
#set(ENABLE_PARMETIS ON CACHE BOOL "...")

