#!/bin/bash

if [ -e $1 ]
  then
    rm -r $1
fi

mkdir $1
cp qm IN qm.f fit.f90 $1
cp mod.f90 lattice-file-180 $1
#cp *.py $1
cp derivs.f $1
cp darter $1
cd $1
qsub darter
