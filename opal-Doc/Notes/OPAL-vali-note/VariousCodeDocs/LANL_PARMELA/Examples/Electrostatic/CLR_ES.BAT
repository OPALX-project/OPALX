Rem  This procedure deletes all files created by Run_ES.BAT in 
Rem  each directory.  Run from a Windows tool such as My Computer.
Rem  Running from DOS windows will not work because of long filenames.                  

Rem  Copyright 1998, by the University of California. 
Rem  Unauthorized commercial use is prohibited. 

cd "%SFDIR%Examples\Electrostatic"
del /s /q  *.inf
del /s /q  *.t35
del /s /q  out*.*
del /s /q  *.log
del /s /q  *.sfo
del /s /q  *.hgl
del /s /q  tape*.*
del /s /q  *.prf
del /s /q  *.tbl
del /s /q  *.qkp
del /s /q  *.t7
del /s /q  wsfprf.txt
del /s /q  ?plotprf.txt

