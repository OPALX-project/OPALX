#!/bin/bash
# do the analysis
#
# usage

sleepTime=600
export wt="01:00:00";
export np=16
export inpn=testcycl-5-P16

curd=/scratch/yang/Grid-16-gele-Ring/
while [ 1 ] 
  do
  qstat -u yang
  sleep $sleepTime
  
  export flag=`qstat | grep yang | awk '{print $5}'` 
  echo qflag=$flag
  # check if running job is waiting not in the queue 
  if [ "$flag" != 'Q' ]; then
      for d in $( ls -d -1 Ring-* ); do 
	  cd $curd/$d
	  echo check in `pwd`
	  tail --lines=1 *.out > A
	  sleep 30
	  tail --lines=1 *.out > B
	  diff A B > q
	  if [ "$?" -ne "-0" ]; then
	      echo the former job is running, please wait ...
	      rm -rf A B q 
	  else
	      grep "Wall" *.out > A
	      if [ "$?" -ne "0" ]; then
		  tail --lines=1 *.out | grep "testcycl-5.in" >A
		  if [ "$?" -ne "0" ]; then
		      MYPID=`cat .status | awk '{print $4}'`
		      qstat $MYPID > t
		      tail -1 t >> tt
		      grep Q tt > t
		      if [ "$?" -ne "0" ]; then
			  echo Dead process in `pwd` kill PID $MYPID
			  qdel $MYPID
			  sleep 5
			  rm .status
		      fi
		      rstep=`grep "Dump step" $inpn.out | tail -1 | awk '{print $5}'`
		      echo restart from dump step $rstep
		      rm -rf run-cycl-2
		      cp ../run-cycl-2 .
		      perl -pi -e "s|STEP|$rstep|g" run-cycl-2
		      cat run-cycl-2
		      FIRST=`qsub -l size=$np,walltime=$wt run-cycl-2`
		      FIRST=${FIRST%.*}
		      sleep 5
		      echo "Run in " `pwd` " with ID= "$FIRST
		      echo `pwd` " with ID= "$FIRST > .status 
		  else 
		      echo "formor job is Reading h5 file,please sleep ... "
		      sleep $sleepTime
		  fi
	      else
		  echo Simulation D O N E
		  exit
	      fi
	      rm -rf A B q t tt
	  fi 
	  cd ..
	  cd $curd
      done
  else
      echo process is waiting at Queue    
  fi
done
