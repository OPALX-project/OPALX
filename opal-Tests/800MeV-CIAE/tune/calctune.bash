#!/bin/bash
rm -f data.dat
rm -f plotdata.dat 
rm -f tuneresult

exec 6<800MeV9S-seo.dat
N="71"
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
#  opal 800MeV9S-tune.in --commlib mpi --info 0 | grep "Max" >>tuneresult
  opal 800MeV9S-tune.in --commlib mpi --info 0 > out
#  opal-20091208 800MeV9S-tune.in --commlib mpi --info 0 > out
  grep "Max" out  >>tuneresult
  j=$[$j+1]
done
exec 6<&-
#rm -rf data.dat

# post porcess
# read in the result of xxxx and OPAL together into file plotdata.dat  
exec 8<tuneresult
#exec 9<refdata
rm -f plotdata.dat 
i="0"

echo "# nu_r(OPAL-CYCL)   nu_z(OPAL-CYCL)" >>plotdata.dat 
while [ $i -lt $N ]
  do 
  read -u 8 a b ur1 d 
  read -u 8 aa bb uz1 dd 
# read -u 9 E ur2 uz2
# format : ur_opal-cycl, ur_opal-cycl
  echo "$ur1   $uz1" >>plotdata.dat 
  i=$[$i+1]
done

exec 8<&-
exec 9<&-
rm -f tunegresult
