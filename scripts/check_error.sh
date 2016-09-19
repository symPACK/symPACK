output=$1
numexp=$2
problem=$3
p=4
b=12
rm $output
rm supernodes_*
rm matrix_output_*
for i in `seq 1 $numexp`
do 
  echo "--------------$i-----------">>$output; 
  mpirun -n $p ./run_sympack -in $problem -inf HARWELL_BOEING -lb SUBCUBE-FI -ordering METIS -fb 1 -b $b -ir 0 -is 0 -nrhs 1 -map ROW2D -refine NO -fact LDL 2>/dev/null >> $output

  rm matrix_output_$i

  for cp in `seq 0 $(($p-1))`
  do
    cat logTest$cp | grep sparse | head -n 1 >> matrix_output_$i
  done

  for cp in `seq 0 $(($p-1))`
  do
    cat logTest$cp | grep sparse | tail -n 1 >> matrix_output_$i
  done

  cat logTest* | grep "Supernode partition" > supernodes_$i
  cat logTest0 | grep "RHS" > solution_$i
  cat logTest0 | grep "X" >> solution_$i
done
diff -q --from-file supernodes_* ; echo $?
cat $output | grep Norm | sort -k 7 -g | tail
cat $output | grep Norm | sort -k 7 -g -r | tail
