Rem  Uncomment this line to skip directly to the Parmela calculation. 
Rem  goto RunParmela


echo off
Echo This file generates the sine and cosine solutions needed for 
Echo the and TRWCfield line in Parmela for the SLAC accelerator.
Echo Fields in the cosine solution are normalized to 1 MV/m over
Echo the last half cell of the problem.  The sine solution uses 
Echo the normalization that produces the same stored energy as the 
Echo cosine solution. Program SF7 is run with "NoNormalization" 
Echo on the command line.  Note: Poisson Superfish version 7.11 
Echo and earlier do not recognize this keyword.   

Echo Delete old SF7 input files, if any, before starting. 
del /q   slac_???.t7   >nul
del /q   slac1??.t7   >nul

Echo First generate the cosine solution consisting Superfish fields 
Echo for the 2pi/3 mode of the traveling-wave structure with Neumann
Echo boundary conditions at left and right edges.
pause

start /w autofish   slac_c.af
start /w WSFplot    slac_c.t35    3
start /w sf7        slac_c.in7
start /w quikplot   rffld.qkp       3

Echo Next generate the sine solution consisting Superfish fields 
Echo for the 2pi/3 mode with Dirichlet boundary conditions at left 
Echo and right edges.
pause

start /w autofish   slac_s.af
start /w WSFplot    slac_s.t35    3
start /w sf7        slac_s.in7
start /w quikplot   rffld.qkp       3



Echo Next generate the Superfish fields for a half cell of the structure
Echo with a beam tube attached.
pause

start /w autofish   slac1.af
start /w WSFplot    slac1.t35    3
start /w sf7        slac1.in7
start /w quikplot   rffld.qkp    3



:RunParmela

start /w Parmela travelwv.acc


Echo Use the program menu to cycle through the displays in PARGRAF.
pause

start /w pargraf travelwv.pgf

Echo Next run Tablplot on the .TBL files to view the fields.

pause

start /w  tablplot RFFLD000.TBL 3
start /w  tablplot RFFLD030.TBL 3
start /w  tablplot RFFLD060.TBL 3
start /w  tablplot RFFLD090.TBL 3
start /w  tablplot RFFLD120.TBL 3
start /w  tablplot RFFLD150.TBL 3
start /w  tablplot RFFLD180.TBL 3
start /w  tablplot RFFLD210.TBL 3
start /w  tablplot RFFLD240.TBL 3
start /w  tablplot RFFLD270.TBL 3
start /w  tablplot RFFLD300.TBL 3
start /w  tablplot RFFLD330.TBL 3

Echo Batch file ends after next entry.
pause

