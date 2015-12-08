# !/bin/sh

ofile=tmp
problem="../utilities/16x16x16_27.rb"
#problem="../utilities/24x24x24_27.rb"
b=20
#b=20
#b=40
LB=NOLB
ordering=SCOTCH
#ordering=METIS
#ordering=NDBOX
#ordering=MMD
#ordering=AMD

export GASNET_MAX_SEGSIZE=512M

#for p in 2 4 8 16 24 32 64
#for p in 2 4 8 16 24 64 81 100 169
for p in 81
do 
  #for map in "ROW2D" "COL2D" "MODWRAP2D"
  for map in "COL2D" i
  do
    echo -n "$p, $b, $map, " >> $ofile;
    mpirun -n $p ./run_sparse_pull -in $problem -inf HARWELL_BOEING -map $map -b $b -ir 0 -is 0 -fb 1 \
    -ordering $ordering -lb $LB -scheduler DL -relax -1 \
    | grep -e "Total\|NNZ" \
    | awk '{print $5","}' | sed -e '{
N
N
N
s/\n/ /g
}' >> $ofile
  done
done
