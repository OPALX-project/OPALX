#!/bin/bash

export wt="00:10:00";

# for dul-cores machine, set 32 for 64 processor
for n in 32  
do
 export np=$n
 rm -rf Ring-$np\P
 mkdir Ring-$np\P
 cd Ring-$np\P
 ln -s ../NewRing-multiBunch.in .
 ln -s ../run-cycl-1 .
 FIRST=`qsub -l size=$np,walltime=$wt run-cycl-1`
 FIRST=${FIRST%.*}
 sleep 5
 echo "Run in " `pwd` " with ID= "$FIRST
 echo `pwd` " with ID= "$FIRST > .status 
 cd ..
done
exit 0



