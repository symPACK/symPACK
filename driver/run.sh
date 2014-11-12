#! /bin/bash

if [[ -z "$5" ]]
then
s=12
else
s=$5
fi

export GASNET_BACKTRACE=1;
#mpirun -np $1 ./run_dense_fanboth_upcxx -in $2 -inf HARWELL_BOEING -b $3 -pref $4
#aprun -n $1 -S $s ./run_dense_fanboth_upcxx -in $2 -inf HARWELL_BOEING -b $3 -pref $4
echo "aprun -n $1 -S $s ./run_dense_fanboth_upcxx -in $2 -inf HARWELL_BOEING -b $3 -pref $4"
echo "mpirun -np $1 ./run_dense_fanboth_upcxx -in $2 -inf HARWELL_BOEING -b $3 -pref $4"



for p in 512; do echo "-------------$p-----------"; for i in `seq 1 5`; do rm core; rm logTest*; rm progress*; aprun -n $p ./run_sparse_fanboth_upcxx
 -in ../matrices/boneS10/boneS10.csc -inf CSC -map AntiDiag2D -b 100 -ir 0 -is 0 -fb 1 | grep "Factorization time" ; done; done >>Results.txt
