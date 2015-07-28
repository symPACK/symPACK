#!/bin/sh


NGCHOLDIR=$HOME/ngchol_scratch
#EXEC=$NGCHOLDIR/driver/run_sparse_fanboth_upcxx
EXECS=($NGCHOLDIR/driver/run_sparse_pull $NGCHOLDIR/driver/run_sparse_fanboth_upcxx)

#EXEC=$NGCHOLDIR/driver/run_sparse_fanboth_upcxx_nolb


#MATRIX=$NGCHOLDIR/matrices/boneS10/boneS10.csc
#MATRIX=$HOME/pexsi_project/UFL-Matrices/pwtk/pwtk.rb
#MATRIX=$HOME/pexsi_project/DNA_715_64cell/H.symm.csc
#MATRIX=$HOME/pexsi_project/LU_C_BN_C_1by1/H.symm.csc
MATRIX=$HOME/pexsi_project/LU_DNA/H.symm.csc


MAPPING=Row2D
ORDER=AMD
LB=SUBCUBE
NUMTRIES=1
FB=1
ISEND=100
IRECV=-1
#BLOCKS=(58)
BLOCKS=(50)
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
PROCS=(64)
#PROCS=(1 4 16 24 25 48 64 81 121 169 256 361)
#PROCS=(192 256)
#PROCS=(1024 768 512 384 256 192 128 96 64)
#PROCS=(1225 1444 1600 1764 1936 2025)

#QUEUE=regular
#WALLTIME=01:20:00
QUEUE=debug
WALLTIME=00:30:00

for BLOCK in "${BLOCKS[@]}"
do
for p in "${PROCS[@]}"
do
JOBNAME="NGCHOL_`basename $MATRIX '.csc'`_${p}_FB_${FB}_${MAPPING}_${BLOCK}_$SCHEDULER"

rm tmp.pbs
echo "#PBS -N $JOBNAME" >> tmp.pbs
echo "#PBS -q $QUEUE" >> tmp.pbs
echo "#PBS -l mppwidth=$p" >> tmp.pbs
echo "#PBS -l walltime=$WALLTIME" >> tmp.pbs
echo "#PBS -j eo" >> tmp.pbs
echo "#PBS -V" >> tmp.pbs
echo "#PBS -S /bin/sh" >> tmp.pbs
echo "#PBS -A mp127" >> tmp.pbs
echo "#PBS -m ae" >> tmp.pbs
if [ $p -gt 8 ]
then
#cho "export GASNET_MAX_SEGSIZE=1.5G" >> tmp.pbs
echo "export GASNET_MAX_SEGSIZE=1.8G" >> tmp.pbs
else
echo "export GASNET_MAX_SEGSIZE=3.5G" >> tmp.pbs
fi
echo "echo \"-------------$p-----------\"" >> tmp.pbs

echo "for EXEC in ${EXECS[@]};" >> tmp.pbs
echo "do " >> tmp.pbs

echo "echo \"-------------\$EXEC-----------\"" >> tmp.pbs
echo "echo \"aprun -n $p \$EXEC -in $MATRIX  -inf CSC -map $MAPPING -b $BLOCK -ir $IRECV -is $ISEND -fb $FB -ordering $ORDER -lb $LB -scheduler $SCHEDULER\"" >> tmp.pbs
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
echo "aprun -n $p \$EXEC -in $MATRIX -inf CSC -map $MAPPING -b $BLOCK -ir $IRECV -is $ISEND -fb $FB -ordering $ORDER -lb $LB -scheduler $SCHEDULER" >> tmp.pbs
echo "done" >> tmp.pbs
echo "done" >> tmp.pbs

qsub tmp.pbs



done
done
