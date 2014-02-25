#! /bin/bash
export GASNET_BACKTRACE=1;

if [[ -z "$5" ]]
then
s=12
else
s=$5
fi



#mpirun -np $1 ./run_dense_fanboth_mpi -in $2 -inf HARWELL_BOEING -b $3
aprun -n $1 ./run_dense_fanboth_mpi -in $2 -inf HARWELL_BOEING -b $3
