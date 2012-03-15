
echo on
start /w automesh soll2.am
start /w poisson  soll2.t35
start /w WSFplot  soll2.t35  3

echo off
echo The next command runs postprocessor SF7 with input 
echo from file SOLL2.IN7. The code creates file SOLL2.T7, 
echo which contains the input data for Parmela, and it 
echo creates file BFLD.QKP, which contains representative 
echo data from file SOLL2.T7. After SF7 ends, Quikplot 
echo displays the fields in file BFLD.QKP.
pause

start /w sf7 soll2.in7
start /w quikplot bfld.qkp  3

echo -------------------------------------------
Echo File SOLL2.T7 contains field-map data for Parmela.   
Echo Batch file ends after next entry. 
pause

