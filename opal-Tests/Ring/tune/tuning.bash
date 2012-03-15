#!/bin/bash
rm -rf data.dat
rm -f plotdata.dat 
rm -f tuningresult

exec 6<FIXPO_SEO
N="257"
j="1"
while read -u 6 E1 r  pr
# read in Energy, initil R and intial Pr of each SEO from FIXPO output
  do
  rm -rf data.dat
  echo -n j = 
  echo "$j"
  echo -n energy= > data.dat 
  echo -n "$E1"  >> data.dat
  echo  ";" >> data.dat

  echo -n r=  >> data.dat
  echo -n "$r" >> data.dat
  echo ";" >> data.dat

  echo -n pr= >> data.dat
  echo -n "$pr" >> data.dat
  echo  ";" >> data.dat
# execute OPAL to calculate tuning frquency and store
  opal tune.in --commlib mpi --info 0 | grep "Max" >>tuningresult
   j=$[$j+1]
done
exec 6<&-
rm -rf data.dat

# post porcess
# read in the result of FIXPO and OPAL together into file plotdata.dat  
exec 8<tuningresult
exec 9<FIXPOUdata 
rm -f plotdata.dat 
i="0"

echo "# Energy   nu_r(FIXPO)   nu_z(FIXPO)   nu_r(OPAL-CYCL)   nu_z(OPAL-CYCL)" >>plotdata.dat 
while [ $i -lt $N ]
  do 
  read -u 8 a b ur1 d 
  read -u 8 aa bb uz1 dd 
  read -u 9 E ur2 uz2
# format : eneryg, ur_fixpo, uz_fixpo, ur_opal-cycl, ur_opal-cycl
  echo "$E   $ur2   $uz2   $ur1   $uz1" >>plotdata.dat 

  i=$[$i+1]
done

exec 8<&-
exec 9<&-
rm -f tuningresult
