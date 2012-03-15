#!/bin/bash

export wt="21:00:00";

# for dul-cores machine, set 32 for 64 processor
export np=32

for d in $( ls -d -1  *mA ); do
 cd $d/ 
 ln -s ../run-cycl-1 .
 FIRST=`qsub -l size=$np,walltime=$wt run-cycl-1`
 FIRST=${FIRST%.*}
 sleep 5
 echo "Run in " `pwd` " with ID= "$FIRST
 echo `pwd` " with ID= "$FIRST > .status 
 cd ..
done
exit 0


