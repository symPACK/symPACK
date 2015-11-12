#!/bin/sh


NGCHOLDIR=$HOME/ngchol_scratch
#EXECS=($NGCHOLDIR/driver/run_sparse_pull $NGCHOLDIR/driver/run_sparse_fanboth_upcxx $NGCHOLDIR/driver/run_mumps)
EXECS=($NGCHOLDIR/driver/run_sparse_pull $NGCHOLDIR/driver/run_sparse_fanboth_upcxx)
#EXECS=($NGCHOLDIR/driver/run_sparse_pull $NGCHOLDIR/driver/run_mumps)
#EXECS=($NGCHOLDIR/driver/run_sparse_pull)

#EXEC=$NGCHOLDIR/driver/run_sparse_fanboth_upcxx_nolb


#MATRIX=$NGCHOLDIR/matrices/boneS10/boneS10.csc
#MATRIX=$NGCHOLDIR/matrices/boneS10/boneS10.csc
#MATRIX=$HOME/pexsi_project/UFL-Matrices/pwtk/pwtk.rb
MATRIX=$HOME/pexsi_project/UFL-Matrices/audikw_1/audikw_1.rb
#MATRIX=$HOME/pexsi_project/UFL-Matrices/ngchol/boneS10/boneS10.rb
#MATRIX=$HOME/pexsi_project/DNA_715_64cell/H.symm.csc
#MATRIX=$HOME/pexsi_project/LU_C_BN_C_1by1/H.symm.csc
#MATRIX=$HOME/pexsi_project/LU_DNA/H.symm.csc
FORMAT="HARWELL_BOEING"

MAPPING=ROW2D
#TREE
ORDER=SCOTCH
LB=SUBCUBE-FO
NUMTRIES=1
FB=1
ISEND=100
IRECV=0
#BLOCKS=(58)
BLOCKS=(100)
#BLOCKS=(100)
#RELAX="-relax -1"
SCHEDULER=DL

#BLOCK=30
#PROCS=(24)
#PROCS=(4 16 24)
#PROCS=(32 48)
#PROCS=(512 768 1024)
#PROCS=(25 169 121 361)
#PROCS=(484 729 1024)
#PROCS=(1 4 16 24 32 48 64)
#PROCS=(96 128 192 256 384)
#PROCS=(1 4 16 24 32 48 64 96 128 192 256 384)
#PROCS=(24)
#ALLOCPROCS=(24 24 24 24 72 96 144 192 264 384 576)
PROCS=(1 4 16 24 64 96 128 192 256 384 576)

#ALLOCPROCS=(768 1032 2136 4104)
PROCS=(768 1024 2116 4096)
#PROCS=(128)

#PROCS=(1 4 16 24 25 48 64 81 121 169 256 361)
#PROCS=(192 256)
#PROCS=(1024 768 512 384 256 192 128 96 64)
#PROCS=(1225 1444 1600 1764 1936 2025)

QUEUE=regular
#WALLTIME=01:20:00
#QUEUE=debug
WALLTIME=00:30:00

for BLOCK in "${BLOCKS[@]}"
do
for p in "${PROCS[@]}"
do
JOBNAME="NGCHOL_`basename $MATRIX '.csc'`_${p}_FB_${FB}_${MAPPING}_${BLOCK}_$SCHEDULER_${LB}"
ALLOCPROC=$(echo "$(printf "%.0f\n" $(echo "scale=2;$p/24 +0.5" | bc))*24" | bc)
rm tmp.pbs
echo "#PBS -N $JOBNAME" >> tmp.pbs
echo "#PBS -q $QUEUE" >> tmp.pbs
echo "#PBS -l mppwidth=$ALLOCPROC" >> tmp.pbs
echo "#PBS -l walltime=$WALLTIME" >> tmp.pbs
echo "#PBS -j eo" >> tmp.pbs
echo "#PBS -V" >> tmp.pbs
echo "#PBS -S /bin/sh" >> tmp.pbs
echo "#PBS -A mp127" >> tmp.pbs
echo "#PBS -m ae" >> tmp.pbs
if [ $p -gt 24 ]
then
echo "export GASNET_MAX_SEGSIZE=1.2G" >> tmp.pbs
#echo "export GASNET_MAX_SEGSIZE=1.8G" >> tmp.pbs
else
echo "export GASNET_MAX_SEGSIZE=3.5G" >> tmp.pbs
fi
echo "echo \"-------------$p-----------\"" >> tmp.pbs

echo "for EXEC in ${EXECS[@]};" >> tmp.pbs
echo "do " >> tmp.pbs

echo "echo \"-------------\$EXEC-----------\"" >> tmp.pbs
echo "echo \"aprun -n $p \$EXEC -in $MATRIX  -inf $FORMAT -map $MAPPING -b $BLOCK -ir $IRECV -is $ISEND -fb $FB -ordering $ORDER -lb $LB -scheduler $SCHEDULER $RELAX\"" >> tmp.pbs
echo "for i in \`seq 1 $NUMTRIES\`;" >> tmp.pbs
echo "do " >> tmp.pbs
echo "workdir=$NGCHOLDIR/experiments/$JOBNAME/\`basename \$EXEC\`/\$i" >> tmp.pbs
echo "echo \$workdir" >> tmp.pbs
echo "mkdir --parents \$workdir" >> tmp.pbs
echo "cd \$workdir" >> tmp.pbs
echo "rm -f core;" >> tmp.pbs
echo "rm -f logTest*;" >> tmp.pbs
echo "rm -f progress*;" >> tmp.pbs
echo "rm -f profile*;" >> tmp.pbs
#echo "aprun -n $p \$EXEC -in $MATRIX -inf HARWELL_BOEING -map $MAPPING -b $BLOCK -ir $IRECV -is $ISEND -fb $FB -ordering $ORDER -lb $LB -scheduler $SCHEDULER" >> tmp.pbs
#echo "aprun -n $p \$EXEC -in $MATRIX -inf CSC -map $MAPPING -b $BLOCK -ir $IRECV -is $ISEND -fb $FB -ordering $ORDER -lb $LB -scheduler $SCHEDULER $RELAX" >> tmp.pbs
echo "aprun -n $p \$EXEC -in $MATRIX -inf $FORMAT -map $MAPPING -b $BLOCK -ir $IRECV -is $ISEND -fb $FB -ordering $ORDER -lb $LB -scheduler $SCHEDULER $RELAX -R 50 -Fact 1 -Ord $ORDER" >> tmp.pbs
echo "done" >> tmp.pbs
echo "done" >> tmp.pbs

qsub tmp.pbs



done
done
