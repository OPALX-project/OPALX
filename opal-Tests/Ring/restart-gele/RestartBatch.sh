#!/bin/bash
# do the analysis
#
# usage

sleepTime=120
export wt="01:59:00";
export np=32
export outn=NewRing-multiBunch-P64.out
export inpn=NewRing-multiBunch
export workdirn=Ring
export myjobn=run-cycl-2
curd=/scratch/yang/Ring_3Bunch
while [ 1 ] 
  do
  qstat -u yang
  sleep $sleepTime
  
  export flag=`qstat | grep yang | awk '{print $5}'` 
  echo qflag=$flag
  # check if running job is waiting not in the queue 
  if [ "$flag" != 'Q' ]; then
      for d in $( ls -d -1  $workdirn-* ); do 
	  cd $curd/$d
	  echo check in `pwd`
	  tail --lines=1 $outn > A
	  sleep 60
	  tail --lines=1 $outn > B
	  diff A B > q
	  if [ "$?" -ne "-0" ]; then
	      echo the former job is running, please wait ...
	      rm -rf A B q 
	  else
	      grep "Wall" $outn > A
	      if [ "$?" -ne "0" ]; then
		  tail --lines=5 $outn | grep "phase space dumping" >A
		  if [ "$?" -eq "0" ]; then
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
		      rstep=`grep "Dump step" $outn | tail -1 | awk '{print $5}'`
		      echo restart from dump step $rstep
		      rm -rf $myjobn
		      cp ../$myjobn .
   
		      h5fn=`ls -rt1 $inpn*.h5|tail -n 1`
		      perl -pi -e "s|STEP|$rstep -restartfn $h5fn|g" $myjobn
		      cat $myjobn
		      FIRST=`qsub -l size=$np,walltime=$wt $myjobn`
		      FIRST=${FIRST%.*}
		      sleep 5
		      echo "Run in " `pwd` " with ID= "$FIRST
		      echo `pwd` " with ID= "$FIRST > .status 
		  else 
		      echo "job is stuck before tracking, please check ... "
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
