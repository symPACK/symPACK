# Building symPACK
--------------------------


**symPACK** is a sparse symmetric matrix direct linear solver. It can be built using the following procedure.

First, download **symPACK** as follows:


```
git clone https://github.com/symPACK/symPACK.git  /path/to/sympack
```

Several environment variables can be optionally set before configuring the build:

- `METIS_DIR` = Installation directory for **MeTiS**

- `PARMETIS_DIR` = Installation directory for **ParMETIS**

- `SCOTCH_DIR` = Installation directory for **SCOTCH** and **PT-SCOTCH**

Then, create a build directory, enter that directory and type:

```
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install/sympack
 ...OPTIONS... /path/to/sympack
```

The `...OPTIONS...` can be one of the following:

* `-DENABLE_METIS=ON|OFF`   to make **MeTiS** ordering available in **symPACK** (`METIS_DIR` must be set in the environment)

* `-DENABLE_PARMETIS=ON|OFF`   to make **ParMETIS** ordering available in **symPACK** (`PARMETIS_DIR` must be set in the environment, `METIS_DIR` is required as well)

* `-DENABLE_SCOTCH=ON|OFF`   to make **SCOTCH** / **PT-SCOTCH** orderings available in **symPACK** (`SCOTCH_DIR` must be set in the environment)



Some platforms have preconfigured toolchain files which can be used by adding the following option to the `cmake` command:
```
-DCMAKE_TOOLCHAIN_FILE=/path/to/sympack/toolchains/edison.cmake     
(To build on NERSC Edison for instance)

```

A sample toolchain file can be found in `/path/to/sympack/toolchains/build_config.cmake` and customized for the target platform.


The `cmake` command will configure the build process, which can now start by typing:
```
make
make install
```

Additionally, a standalone driver for **symPACK** can be built by typing `make examples`


# Troubleshooting

On some systems such as OS X, the build may fail and report that the **libtoolize** command is not found. 
A workaround is to create a symbolic link to the **glibtoolize**, assuming that it is available (it can be installed using Homebrew on OS X for instance).
```
ln -s /usr/local/bin/glibtoolize /usr/local/bin/libtoolize
```
