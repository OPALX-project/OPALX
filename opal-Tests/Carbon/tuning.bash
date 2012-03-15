#!/bin/bash
rm -rf data.dat
rm -f tunedata 
rm -f tuningresult

exec 6<SEO.data
N="200"
j="1"
while read -u 6 E1 r  pr
# read in Energy, initil R and intial Pr of each SEO from FIXPO output
  do
  rm -rf data.dat
  echo -n j = 
  echo "$j"

  echo E= "$E1"
  echo -n energy= > data.dat 
  echo -n "$E1"  >> data.dat
  echo  ";" >> data.dat

  echo r= "$r"
  echo -n r=  >> data.dat
  echo -n "$r" >> data.dat
  echo ";" >> data.dat
  
  echo pr= "$pr"
  echo -n pr= >> data.dat
  echo -n "$pr" >> data.dat
  echo  ";" >> data.dat
# execute OPAL to calculate tuning frquency and store
  opal tuning.in --commlib mpi --info 0 | grep "Max" >>tuningresult
# opal tuning.in --commlib mpi --info 0 | tee tuning.out
   j=$[$j+1]
done
exec 6<&-
rm -rf data.dat

# post porcess
exec 8<tuningresult
exec 6<SEO.data

rm -f tunedata
i="0"
echo "#energy  R[mm]  Pr  nur   nuz  " >>tunedata 
while [ $i -lt $N ]
  do 
  read -u 8 a b ur1 d 
  read -u 8 aa bb uz1 dd 
  read -u 6 E R Pr 
  echo "$E   $R   $Pr  $ur1   $uz1" >>tunedata 
  i=$[$i+1]
done

exec 8<&-
exec 6<&-
rm -f tuningresult
