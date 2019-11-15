# Configuring and building symPACK
--------------------------


**symPACK** is a sparse symmetric matrix direct linear solver. It can be built using the following procedure.

First, download **symPACK** as follows:


```
#!

git clone git@github.com:symPACK/symPACK.git  /path/to/sympack

```

## External dependencies
---------------------------

### UPC++

**SymPACK** requires the UPC++ library to be installed. It can be downloaded at <upcxx.lbl.gov>.
UPC++ contains a CMake config file which **symPACK** is using to link to the library. The install path
needs to be provided to CMake as follows:
```
#!

-DCMAKE_PREFIX_PATH=/home/mjacquel/Work/lbl/pagoda/upcxx_dev/upcxx_install_develop
``` 

### Ordering libraries
Several environment variables can be optionally set before configuring the build:

- `metis_PREFIX` = Installation directory for **MeTiS**

- `parmetis_PREFIX` = Installation directory for **ParMETIS**

- `scotch_PREFIX` = Installation directory for **SCOTCH**

- `ptscotch_PREFIX` = Installation directory for **PT-SCOTCH**

Then, create a build directory, enter that directory and type:

```
#!

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install/sympack -DCMAKE_PREFIX_PATH=/path/to/upcxx/install
 ...OPTIONS... /path/to/sympack

```

The `...OPTIONS...` can be one of the following:

* `-DENABLE_METIS=ON|OFF` to make **MeTiS** ordering available in **symPACK** (`metis_PREFIX` must be set in the environment)

* `-DENABLE_PARMETIS=ON|OFF` to make **ParMETIS** ordering available in **symPACK** (`parmetis_PREFIX` must be set in the environment, `metis_PREFIX` is required as well)

* `-DENABLE_SCOTCH=ON|OFF` to make **SCOTCH** ordering available in **symPACK** (`scotch_PREFIX` must be set in the environment)

* `-DENABLE_PTSCOTCH=ON|OFF` to make **PT-SCOTCH** ordering available in **symPACK** (`ptscotch_PREFIX` must be set in the environment, `scotch_PREFIX` is required as well)

### Notes for specific platforms

Some platforms have preconfigured toolchain files which can be used by adding the following option to the `cmake` command:
```
#!

-DCMAKE_TOOLCHAIN_FILE=/path/to/sympack/toolchains/cori.cmake     
(To build on NERSC Cori for instance)

```

A sample toolchain file can be found in `/path/to/sympack/toolchains/build_config.cmake` and customized for the target platform.

## Building symPACK
---------------------------

The `cmake` command will configure the build process, which can now start by typing:
```
#!

make
make install
```

Additionally, standalone drivers for **symPACK** can be built by typing `make examples`

# Running symPACK
---------------------------

**SymPACK** uses both UPC++ and MPI, therefore, some preliminary environment variables may need to be defined on most platforms.
More details can be found in the UPC++ [documentation](https://bitbucket.org/berkeleylab/upcxx/wiki/docs/mpi-hybrid):
```
#!

export UPCXX_NETWORK=udp #good option on most workstations, but NOT on high performance networks 
export GASNET_SPAWNFN='C'
export GASNET_CSPAWN_CMD='mpirun -np %N %C'
```

To run the standalone **symPACK** driver, for instance on a sample problem from the [Suite Sparse matrix collection](https://sparse.tamu.edu),
the procedure below can be followed:
```
#!

#first download the input matrix (we are using audikw_1 here)
wget https://sparse.tamu.edu/RB/GHS_psdef/audikw_1.tar.gz
#extract the matrix
tar xzf audikw_1.tar.gz

#run sympack
upcxx-run -n 4 -shared-heap 3GB -- ./run_sympack -in audikw_1/audikw_1.rb -ordering METIS -nrhs 1
```

Note that to run **symPACK**, `upcxx-run` is used in the example above, but on some platforms, such as NERSC Cori,
other launchers may be used to both spawn MPI and UPC++, such as `srun` if the system is using SLURM.
