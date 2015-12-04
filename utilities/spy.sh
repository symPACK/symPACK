# !/bin/bash

file=$1

pushd `dirname $0` > /dev/null
SCRIPTPATH=`pwd`
popd > /dev/null

$SCRIPTPATH/../MatrixConverter/sparse_matrix_converter/sparse_matrix_converter $file HB tmp.mtx MM
python $SCRIPTPATH/spy.py tmp.mtx
rm tmp.mtx
