#!/bin/sh


NGCHOLDIR=$HOME/sympack_devel
#EXECS=($NGCHOLDIR/driver/run_sparse_pull $NGCHOLDIR/driver/run_sparse_fanboth_upcxx $NGCHOLDIR/driver/run_mumps)
#EXECS=($NGCHOLDIR/driver/run_sparse_pull $NGCHOLDIR/driver/run_sparse_fanboth_upcxx)
#EXECS=($NGCHOLDIR/driver/run_sparse_pull $NGCHOLDIR/driver/run_mumps)

#EXECS=($NGCHOLDIR/driver/run_sparse_pull_debug $NGCHOLDIR/driver/run_sparse_push_debug $NGCHOLDIR/driver/run_mumps )
#EXECS=( $NGCHOLDIR/driver/run_pastix ) 
#EXECS=( $NGCHOLDIR/driver/run_sparse_pull_debug ) # $NGCHOLDIR/driver/run_sparse_push_debug    )
#EXECS=( $NGCHOLDIR/driver/run_superlu_cori )
#EXECS=( $NGCHOLDIR/driver/run_mumps_cori )
#EXECS=( $NGCHOLDIR/driver/run_elemental_cori )
#EXECS=( $NGCHOLDIR/driver/run_pastix_cori )
#EXECS=( $NGCHOLDIR/driver/run_sparse_pull_cori )
#EXECS=( $NGCHOLDIR/driver/run_sparse_push_cori )
#EXECS=( $NGCHOLDIR/driver/run_sparse_push_cori_mkl_icomm )
#EXECS=( $NGCHOLDIR/driver/run_sparse_pull_cori_mkl )
#EXECS=( $NGCHOLDIR/driver/run_sparse_pull_cori_mkl_icomm )
declare -a EXECS=( "\"$NGCHOLDIR/driver/run_sparse_pull_cori_mkl_icomm -fb dyn\"" $NGCHOLDIR/driver/run_mumps_cori $NGCHOLDIR/driver/run_superlu_cori "\"$NGCHOLDIR/driver/run_sparse_pull_cori_mkl_icomm -fb static\"" "\"$NGCHOLDIR/driver/run_sparse_push_cori_mkl_icomm -fb 1\"" $NGCHOLDIR/driver/run_pastix_cori  $NGCHOLDIR/driver/run_elemental_cori )
declare -a EXECS=( "\"$NGCHOLDIR/driver/run_sparse_pull_cori -fb dyn\"" $NGCHOLDIR/driver/run_elemental_cori $NGCHOLDIR/driver/run_pastix_cori $NGCHOLDIR/driver/run_mumps_cori  "\"$NGCHOLDIR/driver/run_sparse_push_cori -fb 1\""  $NGCHOLDIR/driver/run_superlu_cori )
#declare -a EXECS=( "\"$NGCHOLDIR/driver/run_sparse_pull_cori_mkl_icomm -fb dyn\"" $NGCHOLDIR/driver/run_mumps_cori $NGCHOLDIR/driver/run_pastix_cori  $NGCHOLDIR/driver/run_elemental_cori $NGCHOLDIR/driver/run_superlu_cori)

#declare -a EXECS=( "\"$NGCHOLDIR/driver/run_sparse_pull_cori_mkl_icomm -fb dyn\"" $NGCHOLDIR/driver/run_mumps_cori $NGCHOLDIR/driver/run_superlu_cori "\"$NGCHOLDIR/driver/run_sparse_push_cori_mkl_icomm -fb 1\"" $NGCHOLDIR/driver/run_pastix_cori  $NGCHOLDIR/driver/run_elemental_cori )
#declare -a EXECS=( $NGCHOLDIR/driver/run_superlu_cori "\"$NGCHOLDIR/driver/run_sparse_push_cori_mkl_icomm -fb 1\"" $NGCHOLDIR/driver/run_pastix_cori  $NGCHOLDIR/driver/run_elemental_cori )
#declare -a EXECS=( "\"$NGCHOLDIR/driver/run_sparse_pull_cori_mkl_icomm -fb dyn\""  )


#MATRIX=$NGCHOLDIR/matrices/BenElechi1/BenElechi1.rb
#MATRIX=$NGCHOLDIR/matrices/Flan_1565/Flan_1565.rb
MATRIX=$NGCHOLDIR/matrices/G3_circuit/G3_circuit.rb
#MATRIX=$NGCHOLDIR/matrices/StocF-1465/StocF-1465.rb
#MATRIX=$NGCHOLDIR/matrices/af_shell7/af_shell7.rb
#MATRIX=$NGCHOLDIR/matrices/audikw_1/audikw_1.rb
#MATRIX=$NGCHOLDIR/matrices/bone010/bone010.rb
#MATRIX=$NGCHOLDIR/matrices/boneS10/boneS10.rb
#MATRIX=$NGCHOLDIR/matrices/Serena/Serena.rb
#MATRIX=$NGCHOLDIR/matrices/nd24k/nd24k.rb

#MATRIX=$NGCHOLDIR/matrices/dielFilterV3real/dielFilterV3real.rb
#MATRIX=$NGCHOLDIR/matrices/ex5/ex5.rb
#MATRIX=$NGCHOLDIR/matrices/fv2/fv2.rb
#MATRIX=$NGCHOLDIR/matrices/inline_1/inline_1.rb
#MATRIX=$NGCHOLDIR/matrices/nasa2146/nasa2146.rb
#MATRIX=$NGCHOLDIR/matrices/offshore/offshore.rb
#MATRIX=$NGCHOLDIR/matrices/thermal2/thermal2.rb

#MATRIX=/global/homes/m/mjac/pexsi_project/DNA_715_64cell/DNA_715_64cell.rb




MATRIXCSC=$(echo $MATRIX | sed 's/\.rb/.extend.csc/')
if [ ! -f $MATRIXCSC ]; then
  sparse_matrix_converter $MATRIX HB $MATRIXCSC CSC -e
fi


#MATRIX=$NGCHOLDIR/matrices/pwtk/pwtk.rb
#MATRIX=$HOME/pexsi_project/DNA_715_64cell/H.symm.csc
#MATRIX=$HOME/pexsi_project/LU_C_BN_C_1by1/H.symm.csc
#MATRIX=$HOME/pexsi_project/LU_DNA/H.symm.csc
FORMAT="HARWELL_BOEING"

#MAPPING=ROW2D
#MAPPING=COL2D
MAPPINGS=(ROW2D MODWRAP2D COL2D)
#MAPPINGS=(ROW2D COL2D)
#MAPPINGS=( MODWRAP2D )
#TREE
ORDER=SCOTCH
LB=SUBCUBE-FO
NUMTRIES=1
FB=1
#SUPERLU_METIS
#DYN_MKL
#DYN_LIBSCI
#SUPERLU_METIS
#ELEMENTAL
#dynamic
#MUMPS
#PASTIX
#SUPERLU
#1
#static
ISEND=50
IRECV=0
#BLOCKS=(58)
BLOCKS=(150)
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
#PROCS=(4 16 24 64 96)
#PROCS=(4 16 24 64 96)
#ALLOCPROCS=(24 24 24 24 72 96 144 192 264 384 576)
#PROCS=(1 4 16 24 64 96 128 192) 

PROCS=(1 4 16 24 32 64 96 192 256 384 576 784 1024 2048)
PROCS=(24 32 64 96 192 256 384 576 784 1024 2048)
#PROCS=(1 4 16)

#PROCS=(64 96 192 256 384 576 784 1024)
#PROCS=(1024 1600)
#PROCS=(192 256 384 576)
#PROCS=(384 576 784 1024)

#PROCS=(768 1024 2116 4096)
#PROCS=(64 96 192 256)
#PROCS=(784 1024)
#PROCS=(1 32 96)
#PROCS=(4 16 24 64 192 256 384 576) # 784 1024)

#PROCS=(1 4 16 24 25 48 64 81 121 169 192 256 361)
#PROCS=(192 256)
#PROCS=(1024 768 512 384 256 192 128 96 64)
#PROCS=(1225 1444 1600 1764 1936 2025)

QUEUE=regular
WALLTIME=01:00:00
#QUEUE=debug
#WALLTIME=00:30:00

NODE_MEMORY=40
NODE_CORES=24
NODE_MEMORY=65
NODE_CORES=32

for p in "${PROCS[@]}"
do
#ALLOCNODE=$(echo "$(printf "%.0f\n" $(echo "scale=2;$p/$NODE_CORES +0.5" | bc))" | bc)
#ALLOCNODE=$(echo "($p/$NODE_CORES +0.5)/1" | bc)
ALLOCNODE=$(printf "%.0f\n" $(echo "scale=2;$p/$NODE_CORES +0.49" | bc))
#echo "$p $ALLOCNODE"
ALLOCNODES=( ${ALLOCNODES[@]} $ALLOCNODE )
done

ALLOCNODES=( $(echo "${ALLOCNODES[@]}" | tr ' ' '\n' | sort -u -n | tr '\n' ' ') )
echo ${ALLOCNODES[@]}



for BLOCK in "${BLOCKS[@]}"
do

    for ALLOCNODE in "${ALLOCNODES[@]}"
    do
    
        MINPROC=$(echo "scale=0; ($ALLOCNODE-1)*$NODE_CORES+1" | bc)
        MAXPROC=$(echo "scale=0; $ALLOCNODE*$NODE_CORES" | bc)
        unset CURPROCS 
        for p in "${PROCS[@]}"
        do
#          echo "$p $MINPROC $MAXPROC"
          if [ $p -ge $MINPROC ] && [ $p -le $MAXPROC ]
          then
            CURPROCS=( ${CURPROCS[@]} $p )
          fi
        done
        
        echo ${CURPROCS[@]}
       
#        continue 
        
        
        JOBNAME="NGCHOL_`basename $MATRIX '.rb'`_${ALLOCNODE}_${BLOCK}_$SCHEDULER_${LB}"
        
        rm tmp.slurm
        
        echo "#!/bin/bash -l" >> tmp.slurm
        echo "#SBATCH -p $QUEUE " >> tmp.slurm
        echo "#SBATCH -N $ALLOCNODE " >> tmp.slurm
        echo "#SBATCH -t $WALLTIME" >> tmp.slurm
        echo "#SBATCH -J $JOBNAME" >> tmp.slurm
        echo "#SBATCH --output=$JOBNAME-%J" >> tmp.slurm
        
        echo "declare -a EXECS=(${EXECS[@]})" >> tmp.slurm
        echo "MAPPINGS=(${MAPPINGS[@]})" >> tmp.slurm
        echo "PROCS=(${CURPROCS[@]})" >> tmp.slurm
        echo "for p in \${PROCS[@]};" >> tmp.slurm
        echo "do " >> tmp.slurm
        echo "  SQRTPROC=\$(echo \"scale=0; sqrt(\$p)\" | bc)" >> tmp.slurm
        echo "  if [ \$p -gt $NODE_CORES ]" >> tmp.slurm
        echo "  then" >> tmp.slurm
        echo "    MAXSEGSIZE=\$(echo \"scale=2;($NODE_MEMORY/$NODE_CORES)*0.8\" | bc)" >> tmp.slurm
        echo "  elif [ \$p -eq 1 ]" >> tmp.slurm
        echo "  then" >> tmp.slurm
        echo "    MAXSEGSIZE=$NODE_MEMORY" >> tmp.slurm
        echo "  else" >> tmp.slurm
        echo "    MAXSEGSIZE=\$(echo \"scale=2;($NODE_MEMORY/\$p)*0.8\" | bc)" >> tmp.slurm
        echo "  fi" >> tmp.slurm
        echo "  export GASNET_MAX_SEGSIZE=\${MAXSEGSIZE}G" >> tmp.slurm
        echo "  echo \"-------------\$p-----------\"" >> tmp.slurm
        echo "  echo \"GASNET_MAX_SEGSIZE=\${GASNET_MAX_SEGSIZE}\"" >> tmp.slurm
        echo "  for EXEC in \"\${EXECS[@]}\";" >> tmp.slurm
        echo "  do " >> tmp.slurm
        echo "      read -r -a tmparr <<< \"\$EXEC\" " >> tmp.slurm
        echo "      BASEEXEC=\${tmparr[0]}">> tmp.slurm
        echo "      echo \"-------------\"\$EXEC\"-----------\"" >> tmp.slurm
        echo "      for MAPPING_IDX in \${!MAPPINGS[@]};" >> tmp.slurm
        echo "      do " >> tmp.slurm
        echo "          if [ \"\$MAPPING_IDX\" != \"0\" ] && ( [[ \$(basename \"\$BASEEXEC\") == \"run_mumps\"* ]] || [[ \$(basename \"\$BASEEXEC\") == \"run_pastix\"* ]]  || [[ \$(basename \"\$BASEEXEC\") == \"run_superlu\"* ]]  || [[ \$(basename \"\$BASEEXEC\") == \"run_elemental\"* ]])" >> tmp.slurm
        echo "          then" >> tmp.slurm
        echo "              continue" >> tmp.slurm
        echo "          fi" >> tmp.slurm
        echo "          MAPPING=\${MAPPINGS[\$MAPPING_IDX]}" >> tmp.slurm
        echo "          echo \"-------------\$MAPPING-----------\"" >> tmp.slurm
        echo "          for i in \`seq 1 $NUMTRIES\`;" >> tmp.slurm
        echo "          do " >> tmp.slurm
        echo "              workdir=$NGCHOLDIR/experiments/logs/$JOBNAME-\${SLURM_JOB_ID}/\`basename \"\$BASEEXEC\"\`/\$i" >> tmp.slurm
        echo "              echo \$workdir" >> tmp.slurm
        echo "              mkdir --parents \$workdir" >> tmp.slurm
        echo "              cd \$workdir" >> tmp.slurm
        echo "              rm -f core;" >> tmp.slurm
        echo "              rm -f logTest*;" >> tmp.slurm
        echo "              rm -f progress*;" >> tmp.slurm
        echo "              rm -f profile*;" >> tmp.slurm
        echo "              if [[ \$(basename \"\$BASEEXEC\") == \"run_pastix\"* ]]" >> tmp.slurm
        echo "              then" >> tmp.slurm
        echo "                  echo \"srun -n \$p \$EXEC -t 1 -hb $MATRIX -ord scotch\"" >> tmp.slurm
        echo "                  srun -n \$p \$EXEC -t 1 -hb $MATRIX -ord scotch |grep -v \"vs\" | grep -v \"Factorization time\"" >> tmp.slurm
        echo "              elif [[ \$(basename \"\$BASEEXEC\") == \"run_elemental\"* ]]" >> tmp.slurm
        echo "              then" >> tmp.slurm
        echo "                  echo \"srun -n \$p \$EXEC -in $MATRIX\" " >> tmp.slurm
        echo "                  srun -n \$p \$EXEC -in $MATRIX " >> tmp.slurm
        echo "              elif [[ \$(basename \"\$BASEEXEC\") == \"run_superlu\"* ]]" >> tmp.slurm
        echo "              then" >> tmp.slurm
        #echo "                echo \"srun -n $p \$EXEC -in $MATRIXCSC -inf CSC -colperm PARMETIS -r $SQRTPROC -c $SQRTPROC\"" >> tmp.slurm
        #echo "                srun -n $p \$EXEC -in $MATRIXCSC -inf CSC -colperm PARMETIS -r $SQRTPROC -c $SQRTPROC " >> tmp.slurm
        echo "                  echo \"srun -n \$p \$EXEC -in $MATRIXCSC -inf CSC -colperm METIS_AT_PLUS_A -r \$SQRTPROC -c \$SQRTPROC\"" >> tmp.slurm
        echo "                  srun -n \$p \$EXEC -in $MATRIXCSC -inf CSC -colperm METIS_AT_PLUS_A -r \$SQRTPROC -c \$SQRTPROC  | grep -v \"dQuery_Space\"" >> tmp.slurm
        echo "              else" >> tmp.slurm
        echo "                  echo \"srun -n \$p \$EXEC -in $MATRIX  -inf $FORMAT -map \$MAPPING -b $BLOCK -ir $IRECV -is $ISEND -ordering $ORDER -lb $LB -scheduler $SCHEDULER $RELAX\"" >> tmp.slurm
        echo "                  srun -n \$p \$EXEC -in $MATRIX -inf $FORMAT -map \$MAPPING -b $BLOCK -ir $IRECV -is $ISEND -ordering $ORDER -lb $LB -scheduler $SCHEDULER $RELAX -R 50 -Fact 1 -Ord $ORDER" >> tmp.slurm
        echo "              fi" >> tmp.slurm
        echo "          done" >> tmp.slurm
        echo "          echo \"------------------------------\"" >> tmp.slurm
        echo "      done" >> tmp.slurm
        echo "      echo \"------------------------------\"" >> tmp.slurm
        echo "  done" >> tmp.slurm
        echo "done " >> tmp.slurm
        
        sbatch tmp.slurm
        
    done
done
exit


