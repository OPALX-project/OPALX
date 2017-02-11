  #!/bin/bash
  rm -rf tempfile
  rm -f plotdata
  rm -f tuningresult

  exec 6<FIXPO_SEO
  N="260"
  j="1"
  while read -u 6 E1 r  pr
  # read in Energy, initil R and intial Pr of each
  # SEO from FIXPO output
  do
  rm -rf tempfile
  echo -n j =
  echo "$j"
  echo -n energy= > tempfile
  echo -n "$E1"  >> tempfile
  echo  ";" >> tempfile

  echo -n r=  >> tempfile
  echo -n "$r" >> tempfile
  echo ";" >> tempfile

  echo -n pr= >> tempfile
  echo -n "$pr" >> tempfile
  echo  ";" >> tempfile
  # execute OPAL to calculate tuning frquency and store
  opal testcycl.in --commlib mpi --info 0 | \
     grep "Max" >>tuningresult
  j=$[$j+1]
  done
  exec 6<&-
  rm -rf tempfile

  # post porcess
  exec 8<tuningresult
    rm -f plotdata
  i="0"

  while [ $i -lt $N ]
  do
  read -u 8 a b ur1 d
  read -u 8 aa bb uz1 dd
  echo "$ur1   $uz1" >>plotdata
  i=$[$i+1]
  done

  exec 8<&-
  rm -f tuningresult