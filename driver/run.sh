#! /bin/bash
mpirun -np $1 ./run_dense_fanboth_upcxx -in nasa2146/nasa2146.rb -inf HARWELL_BOEING -b 1 -pref 0
