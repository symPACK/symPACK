#!/bin/sh


NGCHOLDIR=$HOME/ngchol_scratch
EXEC=$NGCHOLDIR/driver/run_sparse_fanboth_upcxx
#EXEC=$NGCHOLDIR/driver/run_sparse_fanboth_upcxx_nolb


MATRIX=$NGCHOLDIR/matrices/boneS10/boneS10.csc


MAPPING=Col2D
ORDER=AMD
LB=SUBCUBE
NUMTRIES=5
FB=1
ISEND=400
IRECV=0
BLOCK=50
#PROCS=(24)
#PROCS=(4 16 24)
#PROCS=(32 48)
PROCS=(1 4 16 24 32 48 64 96 128 192 256 384 512)
#PROCS=(512 768 1024)
#PROCS=(1024 768 512 384 256 192 128 96 64)
#PROCS=(1225 1444 1600 1764 1936 2025)

#QUEUE=regular
#WALLTIME=01:20:00
QUEUE=debug
WALLTIME=00:30:00

for p in "${PROCS[@]}"
do
JOBNAME="NGCHOL_`basename $MATRIX '.csc'`_${p}_FB_${FB}_${MAPPING}"

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
echo "cd $NGCHOLDIR" >> tmp.pbs
echo "echo \"-------------$p-----------\"" >> tmp.pbs
echo "echo \"aprun -n $p $EXEC -in $MATRIX  -inf CSC -map $MAPPING -b $BLOCK -ir $IRECV -is $ISEND -fb $FB -ordering $ORDER -lb $LB\"" >> tmp.pbs
echo "for i in \`seq 1 $NUMTRIES\`;" >> tmp.pbs
echo "do " >> tmp.pbs
echo "rm -f core;" >> tmp.pbs
echo "rm -f logTest*;" >> tmp.pbs
echo "rm -f progress*;" >> tmp.pbs
echo "rm -f profile*;" >> tmp.pbs
echo "aprun -n $p $EXEC -in $MATRIX -inf CSC -map $MAPPING -b $BLOCK -ir $IRECV -is $ISEND -fb $FB -ordering $ORDER -lb $LB" >> tmp.pbs
echo "done" >> tmp.pbs

qsub tmp.pbs



done
