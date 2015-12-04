# !/bin/bash

file=$1

pushd `dirname $0` > /dev/null
SCRIPTPATH=`pwd`
popd > /dev/null
ofile=`basename $file .rb`.mtx
$SCRIPTPATH/../MatrixConverter/sparse_matrix_converter/sparse_matrix_converter $file HB $ofile MM
python $SCRIPTPATH/spy.py $ofile & 
#rm tmp.mtx
