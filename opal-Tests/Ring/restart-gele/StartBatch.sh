#!/bin/bash

export wt="01:00:00";

for n in 16 
do
 export np=$n
 rm -rf Ring-$np\P
 mkdir Ring-$np\P
 cd Ring-$np\P
 ln -s ../testcycl-5.in .
 ln -s ../run-cycl-1 .
 FIRST=`qsub -l size=$np,walltime=$wt run-cycl-1`
 FIRST=${FIRST%.*}
 sleep 5
 echo "Run in " `pwd` " with ID= "$FIRST
 echo `pwd` " with ID= "$FIRST > .status 
 cd ..
done
exit 0
# SECOND=`qsub -Wafterany:i$FIRST -l size=$np,walltime=$wt run-cycl-2`


