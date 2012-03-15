Rem Delete old 1300c1xx.t7 files 
del /q   2piovr??.t7  >nul

start /w autofish   2piovr3.af
start /w sf7        2piovr3.in7
start /w wsfplot    2piovr3.t35  3
start /w quikplot   rffld.qkp    3
echo off 

Echo Accept the default settings for prompts in EFLDTR except 
Echo for the cell length, which is 3.333 cm.  The 11 Fourier 
Echo coefficients produced by this code would appear on the 
Echo first TRWAVE line in the a Parmela input file.
pause

start /w efldtr 2piovr01.t7

Echo Results are in file FORCOEFF.TXT.
Echo Batch file ends after next entry.
pause

