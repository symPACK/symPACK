#! /bin/bash
cd ../src; make cleanall; make -j; cd ../driver; make cleanall; make $1 
