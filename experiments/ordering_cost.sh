#!/bin/bash

ofile=tmp2
tmpfile=tmp3
tmpfile2=tmp4

p=2
b="-1"
map="COL2D"
stencil=27
#b=20
#b=40

orderings=( "SCOTCH" "NDBOX" "METIS" "MMD" "AMD" )

rm -rf $tmpfile
echo -n "k, " >> $tmpfile;
for ordering in "${orderings[@]}"
do 
    echo -n "$ordering, " >> $tmpfile;
done

cat $tmpfile | sed -e '{
N
N
N
N
N
N
s/\n/ /g
}' >> $ofile

echo " " >> $ofile;

#for k in 3 9 16 24 32 48 64 96
for k in 128 192 256
do
rm -rf $tmpfile2
echo -n "$k, " >> $tmpfile2;


problem="../utilities/${k}x${k}x${k}_${stencil}.rb"
pushd ../utilities
./3d_adjacency $k $k $k $stencil > tmp.adj; python adj_to_hb.py tmp.adj > ${k}x${k}x${k}_${stencil}.rb
popd

for ordering in "${orderings[@]}"
do 
    mpirun -n $p ./run_sparse_pull_nocomp -in $problem -inf HARWELL_BOEING -map $map -ir 0 -is 0 -fb 1 \
    -ordering $ordering -lb SUBCUBE-FO -scheduler DL -relax -1 \
    | grep  "NNZ" \
    | awk '{print $5","}' | tr '\n' ' ' >> $tmpfile2
done

cat $tmpfile2 | sed -e '{
N
N
N
N
N
s/\n/ /g
}' >> $ofile
echo " " >> $ofile;


done

#rm -rf $tmpfile
#rm -rf $tmpfile2

