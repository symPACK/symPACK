set(BLAS_DIR /opt/intel/mkl)
set(LAPACK_DIR /opt/intel/mkl)
#Required for mkl
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${BLAS_DIR}/lib/intel64)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${LAPACK_DIR}/lib/intel64)

#set(BLAS_DIR /home/mjacquel/software/gnu/release/OpenBLAS/install)
#set(LAPACK_DIR /home/mjacquel/software/gnu/release/OpenBLAS/install)
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${BLAS_DIR})
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${LAPACK_DIR})

set(ptscotch_PREFIX /home/mjacquel/software/gnu/release/scotch_6.0.4_install)
set(ENABLE_PTSCOTCH ON CACHE BOOL "...")

set(scotch_PREFIX /home/mjacquel/software/gnu/release/scotch_6.0.4_install)
#set(SCOTCH_DIR /home/mjacquel/software/gnu/release/scotch_6.0.4)
set(ENABLE_SCOTCH ON CACHE BOOL "...")


set(metis_PREFIX /home/mjacquel/software/gnu/release/parmetis_install)
#set(METIS_DIR /home/mjacquel/software/gnu/release/parmetis_install)
set(ENABLE_METIS ON CACHE BOOL "...")

set(parmetis_PREFIX /home/mjacquel/software/gnu/release/parmetis_install)
#set(PARMETIS_DIR /home/mjacquel/software/gnu/release/parmetis_install)
set(ENABLE_PARMETIS ON CACHE BOOL "...")

