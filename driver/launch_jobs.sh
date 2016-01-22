#!/bin/sh


NGCHOLDIR=$HOME/sympack_devel
#EXECS=($NGCHOLDIR/driver/run_sparse_pull $NGCHOLDIR/driver/run_sparse_fanboth_upcxx $NGCHOLDIR/driver/run_mumps)
#EXECS=($NGCHOLDIR/driver/run_sparse_pull $NGCHOLDIR/driver/run_sparse_fanboth_upcxx)
#EXECS=($NGCHOLDIR/driver/run_sparse_pull $NGCHOLDIR/driver/run_mumps)
EXECS=($NGCHOLDIR/driver/run_sparse_pull_time $NGCHOLDIR/driver/run_mumps)

#EXEC=$NGCHOLDIR/driver/run_sparse_fanboth_upcxx_nolb



#MATRIX=$NGCHOLDIR/matrices/BenElechi1/BenElechi1.rb
#MATRIX=$NGCHOLDIR/matrices/Flan_1565/Flan_1565.rb
#MATRIX=$NGCHOLDIR/matrices/G3_circuit/G3_circuit.rb
#MATRIX=$NGCHOLDIR/matrices/StocF-1465/StocF-1465.rb
#MATRIX=$NGCHOLDIR/matrices/af_shell7/af_shell7.rb
MATRIX=$NGCHOLDIR/matrices/audikw_1/audikw_1.rb
#MATRIX=$NGCHOLDIR/matrices/bone010/bone010.rb
#MATRIX=$NGCHOLDIR/matrices/boneS10/boneS10.rb
#MATRIX=$NGCHOLDIR/matrices/dielFilterV3real/dielFilterV3real.rb
#MATRIX=$NGCHOLDIR/matrices/ex5/ex5.rb
#MATRIX=$NGCHOLDIR/matrices/fv2/fv2.rb
#MATRIX=$NGCHOLDIR/matrices/inline_1/inline_1.rb
#MATRIX=$NGCHOLDIR/matrices/nasa2146/nasa2146.rb
#MATRIX=$NGCHOLDIR/matrices/offshore/offshore.rb
#MATRIX=$NGCHOLDIR/matrices/thermal2/thermal2.rb










#MATRIX=$NGCHOLDIR/matrices/pwtk/pwtk.rb
#MATRIX=$HOME/pexsi_project/DNA_715_64cell/H.symm.csc
#MATRIX=$HOME/pexsi_project/LU_C_BN_C_1by1/H.symm.csc
#MATRIX=$HOME/pexsi_project/LU_DNA/H.symm.csc
FORMAT="HARWELL_BOEING"

#MAPPING=ROW2D
#MAPPING=COL2D
MAPPINGS=(MODWRAP2D COL2D ROW2D)
#TREE
ORDER=SCOTCH
LB=SUBCUBE-FO
NUMTRIES=2
FB=1
ISEND=10
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
#PROCS=(1 4 16 24 64 96 128 192) 
PROCS=(1 4 16 24 64 96 192 256 384 576)

#ALLOCPROCS=(768 1032 2136 4104)
#PROCS=(768 1024 2116 4096)
#PROCS=(24)

#PROCS=(1 4 16 24 25 48 64 81 121 169 192 256 361)
#PROCS=(192 256)
#PROCS=(1024 768 512 384 256 192 128 96 64)
#PROCS=(1225 1444 1600 1764 1936 2025)

QUEUE=regular
#WALLTIME=01:20:00
QUEUE=debug
WALLTIME=00:30:00

NODE_MEMORY=40
NODE_CORES=24

for BLOCK in "${BLOCKS[@]}"
do
for p in "${PROCS[@]}"
do
#JOBNAME="NGCHOL_`basename $MATRIX '.csc'`_${p}_FB_${FB}_${MAPPING}_${BLOCK}_$SCHEDULER_${LB}"
JOBNAME="NGCHOL_`basename $MATRIX '.rb'`_${p}_FB_${FB}_${BLOCK}_$SCHEDULER_${LB}"
ALLOCNODES=$(echo "$(printf "%.0f\n" $(echo "scale=2;$p/$NODE_CORES +0.5" | bc))" | bc)
rm tmp.slurm

echo "#!/bin/bash -l" >> tmp.slurm
echo "#SBATCH -p $QUEUE " >> tmp.slurm
echo "#SBATCH -N $ALLOCNODES " >> tmp.slurm
echo "#SBATCH -t $WALLTIME" >> tmp.slurm
echo "#SBATCH -J $JOBNAME" >> tmp.slurm
echo "#SBATCH --output=$JOBNAME-%J" >> tmp.slurm
#echo "#SBATCH --output=$JOBNAME-%J" >> tmp.slurm

#echo "#PBS -N $JOBNAME" >> tmp.slurm

if [ $p -gt $NODE_CORES ]
then
MAXSEGSIZE=$(echo "scale=2;($NODE_MEMORY/$NODE_CORES)*0.8" | bc)
elif [ $p -eq 1 ]
then
MAXSEGSIZE=40
else
MAXSEGSIZE=$(echo "scale=2;($NODE_MEMORY/$p)*0.8" | bc)
fi

echo $MAXSEGSIZE
echo "export GASNET_MAX_SEGSIZE=${MAXSEGSIZE}G" >> tmp.slurm
#if [ $p -gt 24 ]
#then
#echo "export GASNET_MAX_SEGSIZE=1.2G" >> tmp.slurm
##echo "export GASNET_MAX_SEGSIZE=1.8G" >> tmp.slurm
#else
#echo "export GASNET_MAX_SEGSIZE=3.5G" >> tmp.slurm
#fi
echo "echo \"-------------$p-----------\"" >> tmp.slurm

echo "EXECS=(${EXECS[@]})" >> tmp.slurm
echo "MAPPINGS=(${MAPPINGS[@]})" >> tmp.slurm

echo "for EXEC in \${EXECS[@]};" >> tmp.slurm
echo "do " >> tmp.slurm

echo "echo \"-------------\$EXEC-----------\"" >> tmp.slurm


echo "for MAPPING_IDX in \${!MAPPINGS[@]};" >> tmp.slurm
echo "do " >> tmp.slurm

echo "if [ \"\$MAPPING_IDX\" != \"1\" ] && [ \"\$EXEC\" == \"/global/homes/m/mjac/sympack_devel/driver/run_mumps\" ]" >> tmp.slurm
echo "then" >> tmp.slurm
echo "continue" >> tmp.slurm
echo "fi" >> tmp.slurm

echo "MAPPING=\${MAPPINGS[\$MAPPING_IDX]}" >> tmp.slurm

echo "echo \"-------------\$MAPPING-----------\"" >> tmp.slurm

echo "for i in \`seq 1 $NUMTRIES\`;" >> tmp.slurm
echo "do " >> tmp.slurm
echo "workdir=$NGCHOLDIR/experiments/logs/$JOBNAME/\`basename \$EXEC\`/\$i" >> tmp.slurm
echo "echo \$workdir" >> tmp.slurm
echo "mkdir --parents \$workdir" >> tmp.slurm
echo "cd \$workdir" >> tmp.slurm
echo "rm -f core;" >> tmp.slurm
echo "rm -f logTest*;" >> tmp.slurm
echo "rm -f progress*;" >> tmp.slurm
echo "rm -f profile*;" >> tmp.slurm
#echo "aprun -n $p \$EXEC -in $MATRIX -inf HARWELL_BOEING -map $MAPPING -b $BLOCK -ir $IRECV -is $ISEND -fb $FB -ordering $ORDER -lb $LB -scheduler $SCHEDULER" >> tmp.slurm
#echo "aprun -n $p \$EXEC -in $MATRIX -inf CSC -map $MAPPING -b $BLOCK -ir $IRECV -is $ISEND -fb $FB -ordering $ORDER -lb $LB -scheduler $SCHEDULER $RELAX" >> tmp.slurm

echo "echo \"srun -n $p \$EXEC -in $MATRIX  -inf $FORMAT -map \$MAPPING -b $BLOCK -ir $IRECV -is $ISEND -fb $FB -ordering $ORDER -lb $LB -scheduler $SCHEDULER $RELAX\"" >> tmp.slurm
echo "srun -n $p \$EXEC -in $MATRIX -inf $FORMAT -map \$MAPPING -b $BLOCK -ir $IRECV -is $ISEND -fb $FB -ordering $ORDER -lb $LB -scheduler $SCHEDULER $RELAX -R 50 -Fact 1 -Ord $ORDER" >> tmp.slurm
echo "done" >> tmp.slurm

echo "echo \"------------------------------\"" >> tmp.slurm
echo "done" >> tmp.slurm
echo "echo \"------------------------------\"" >> tmp.slurm
echo "done" >> tmp.slurm

sbatch tmp.slurm



done
done
