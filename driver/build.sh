#! /bin/bash
cd ../src; make cleanall; make -j; cd ../driver; rm $1 $1.o; make $1 
