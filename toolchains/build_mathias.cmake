set(BLAS_DIR /home/mjacquel/software/gnu/release/OpenBLAS/install)
set(LAPACK_DIR /home/mjacquel/software/gnu/release/OpenBLAS/install)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${BLAS_DIR}/lib)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${LAPACK_DIR}/lib)

set(SCOTCH_DIR /home/mjacquel/software/gnu/release/scotch_6.0.4_install)
set(METIS_DIR /home/mjacquel/software/intel/release/parmetis_install)
set(PARMETIS_DIR /home/mjacquel/software/intel/release/parmetis_install)

set(ENABLE_SCOTCH ON CACHE BOOL "...")
set(ENABLE_METIS ON CACHE BOOL "...")
set(ENABLE_PARMETIS ON CACHE BOOL "...")
