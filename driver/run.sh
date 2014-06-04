#! /bin/bash

if [[ -z "$6" ]]
then
s=12
else
s=$6
fi

export GASNET_BACKTRACE=1;
#mpirun -np $1 ./run_dense_fanboth_upcxx -in $2 -inf HARWELL_BOEING -b $3 -pref $4
#aprun -n $1 -S $s ./run_dense_fanboth_upcxx -in $2 -inf HARWELL_BOEING -b $3 -pref $4
echo "aprun -n $1 -S $s ./run_dense_fanboth_upcxx -in $2 -inf HARWELL_BOEING -b $3 -pref $4 -loop $5"
echo "mpirun -np $1 ./run_dense_fanboth_upcxx -in $2 -inf HARWELL_BOEING -b $3 -pref $4 -loop $5"
