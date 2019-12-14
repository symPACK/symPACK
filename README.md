# Configuring and building symPACK
--------------------------


**symPACK** is a sparse symmetric matrix direct linear solver. It can be built using the following procedure.

First, download **symPACK** as follows:


```
git clone git@github.com:symPACK/symPACK.git  /path/to/sympack

```

## External dependencies
---------------------------

### UPC++

**SymPACK** requires the **UPC++ v1.0** library to be installed. The minimum supported release version is 2019.9.0. 
It can be downloaded at [upcxx.lbl.gov](upcxx.lbl.gov).
UPC++ contains a CMake config file which **symPACK** is using to link to the library. The install path
needs to be provided to CMake as follows:
```
-DCMAKE_PREFIX_PATH=/path/to/upcxx-2019.9.0
``` 

### Ordering libraries
Several environment variables can be optionally set before configuring the build:

- `metis_PREFIX` = Installation directory for **MeTiS**

- `parmetis_PREFIX` = Installation directory for **ParMETIS**

- `scotch_PREFIX` = Installation directory for **SCOTCH**

- `ptscotch_PREFIX` = Installation directory for **PT-SCOTCH**

Then, create a build directory, enter that directory and type:

```
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install/sympack -DCMAKE_PREFIX_PATH=/path/to/upcxx/install
 ...OPTIONS... /path/to/sympack

```

The `...OPTIONS...` can be one of the following:

- `-DENABLE_METIS=ON|OFF` to make **MeTiS** ordering available in **symPACK** (`metis_PREFIX` must be set in the environment)

- `-DENABLE_PARMETIS=ON|OFF` to make **ParMETIS** ordering available in **symPACK** (`parmetis_PREFIX` must be set in the environment, `metis_PREFIX` is required as well)

- `-DENABLE_SCOTCH=ON|OFF` to make **SCOTCH** ordering available in **symPACK** (`scotch_PREFIX` must be set in the environment)

- `-DENABLE_PTSCOTCH=ON|OFF` to make **PT-SCOTCH** ordering available in **symPACK** (`ptscotch_PREFIX` must be set in the environment, `scotch_PREFIX` is required as well)

### Notes for specific platforms

Some platforms have preconfigured toolchain files which can be used by adding the following option to the `cmake` command:
```
-DCMAKE_TOOLCHAIN_FILE=/path/to/sympack/toolchains/cori.cmake     
(To build on NERSC Cori for instance)

```

A sample toolchain file can be found in `/path/to/sympack/toolchains/build_config.cmake` and customized for the target platform.

## Building symPACK
---------------------------

The `cmake` command will configure the build process, which can now start by typing:
```
make
make install
```

Additionally, standalone drivers for **symPACK** can be built by typing `make examples`

### Build options

**SymPACK** has several optional features. We have already mentioned the options related to sparse matrix ordering. The following list is the complete set of options available when building **symPACK**:

- `-DENABLE_METIS` to enable **MeTiS** ordering (`ON|OFF`)

- `-DENABLE_PARMETIS` to enable **ParMETIS** ordering (`ON|OFF`)

- `-DENABLE_SCOTCH` to enable **SCOTCH** ordering (`ON|OFF`)

- `-DENABLE_PTSCOTCH` to enable **PT-SCOTCH** ordering (`ON|OFF`)

- `-DAMD_IDX_64` to use 64 bit integers for AMD ordering (`ON|OFF`)

- `-DMMD_IDX_64` to use 64 bit integers for MMD ordering (`ON|OFF`)

- `-DRCM_IDX_64` to use 64 bit integers for RCM ordering (`ON|OFF`)

- `-DENABLE_MKL` to enable Intel MKL through the `-mkl` flag (`ON|OFF`)

- `-DENABLE_THREADS` to enable multithreading (`ON|OFF`). `UPCXX_THREADMODE=par` is required during cmake configuration. **SymPACK** implements its own multithreading and as such should be linked with **sequential** BLAS/LAPACK libraries.

# Running symPACK
---------------------------

**SymPACK** uses both **UPC++** and **MPI**, therefore, some preliminary environment variables may need to be defined on most platforms.
More details can be found in the **UPC++** [documentation](https://bitbucket.org/berkeleylab/upcxx/wiki/docs/mpi-hybrid):
```
export UPCXX_NETWORK=udp #good option on most workstations, but NOT on high performance networks 
export GASNET_SPAWNFN='C'
export GASNET_CSPAWN_CMD='mpirun -np %N %C'
```

To run the standalone **symPACK** driver, for instance on a sample problem from the [Suite Sparse matrix collection](https://sparse.tamu.edu),
the procedure below can be followed:
```
#first download the input matrix (we are using nasa2146 here)
wget https://sparse.tamu.edu/RB/Nasa/nasa2146.tar.gz
#extract the matrix
tar xzf nasa2146.tar.gz

#run sympack
upcxx-run -n 4 -shared-heap 3GB -- ./run_sympack -in nasa2146/nasa2146.rb -ordering MMD -nrhs 1
```
For larger problems, it is strongly advised to use more efficient orderings (such as **METIS**).
Note that to run **symPACK**, `upcxx-run` is used in the example above, but on some platforms, such as NERSC Cori,
other launchers may be used to both spawn **MPI** and **UPC++**, such as `srun` if the system is using SLURM.
Moreover, for larger problems, the `-shared-heap XX` argument to `upcxx-run` may be needed (Please refer to **UPC++** documentation).
