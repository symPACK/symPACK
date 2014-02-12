#! /bin/bash
export GASNET_BACKTRACE=1;
mpirun -np $1 ./run_dense_scalapack -in $2 -inf HARWELL_BOEING -b $3
