set(CMAKE_CXX_COMPILER CC)
set(CMAKE_C_COMPILER cc)
set(CMAKE_Fortran_COMPILER ftn)
set(MPI_C_COMPILER cc)
set(MPI_CXX_COMPILER CC)

set(ptscotch_PREFIX /home/mjacquel/software/gnu/release/scotch_6.0.4_install_openmpi)
set(ptscotch_PREFIX /home/mjacquel/software/gnu/release/scotch_6.0.4_install)
set(ENABLE_PTSCOTCH ON CACHE BOOL "...")

set(scotch_PREFIX /home/mjacquel/software/gnu/release/scotch_6.0.4_install_openmpi)
set(scotch_PREFIX /home/mjacquel/software/gnu/release/scotch_6.0.4_install)
set(SCOTCH_DIR /home/mjacquel/software/gnu/release/scotch_6.0.4)
set(ENABLE_SCOTCH ON CACHE BOOL "...")


set(metis_PREFIX /home/mjacquel/software/gnu/release/parmetis_install)
set(METIS_DIR /home/mjacquel/software/gnu/release/parmetis_install)
set(ENABLE_METIS ON CACHE BOOL "...")

set(parmetis_PREFIX /home/mjacquel/software/gnu/release/parmetis_install)
set(PARMETIS_DIR /home/mjacquel/software/gnu/release/parmetis_install)
set(ENABLE_PARMETIS ON CACHE BOOL "...")

