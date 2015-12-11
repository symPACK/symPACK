# !/bin/sh

ofile=tmp
#problem="../utilities/16x16x16_7.rb"
#problem="../utilities/16x16x16_27.rb"
#problem="../utilities/24x24x24_27.rb"
problem="../utilities/24x24x24_7.rb"
#problem="../utilities/3x3x3_27.rb"
b=20
#b=20
#b=40
LB=NOLB

export GASNET_MAX_SEGSIZE=2G
orderings=( "SCOTCH" "NDBOX" "METIS" "MMD" "AMD" )
mappings=( "MODWRAP2D" "COL2D" "ROW2D" )
mappings=( "COL2D" )

for ordering in "${orderings[@]}"
do 
#for p in 2 4 8 16 24 32 64
#for p in 2 4 8 16 24 64 100 169
#for p in 2 4 8 16 24 64 100 169
#for p in 64
#for p in 2 4 8 16 24 64 100 169
for p in 4
do 
  #for map in "COL2D" 
  for map in "${mappings[@]}"
  do
    echo -n "$p, $b, $ordering, $map, " >> $ofile;
    mpirun -n $p ./run_sparse_pull -in $problem -inf HARWELL_BOEING -map $map -b $b -ir 0 -is 0 -fb 1 \
    -ordering $ordering -lb $LB -scheduler DL -relax -1 \
    | grep -e "Total\|NNZ" \
    | awk '{print $5","}' | tr '\n' ' ' | sed -e '{
N
N
N
s/\n/ /g
}' >> $ofile
    echo " " >> $ofile
  done
done
done
