scp -r gele.cscs.ch:/scratch/adelmann/RegressionTests/DriftTest/OPAL/figs/*.* Drift
cd Drift
for f  in `ls -1 *.eps`
do
 echo convert to pdf $f
 epstopdf $f
 rm $f
done
cd ..
