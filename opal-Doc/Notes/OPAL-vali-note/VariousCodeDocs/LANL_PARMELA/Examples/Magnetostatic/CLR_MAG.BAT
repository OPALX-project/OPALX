Rem  This procedure deletes all OUTxxx and TAPExx files created by the
Rem  Run_MAG.BAT run in each directory.  Run from a Windows tool such
Rem  as My Computer.  Running from DOS windows will not work because of
Rem  long filenames.

Rem  Copyright 1998, by the University of California.
Rem  Unauthorized commercial use is prohibited.


cd  "%SFDIR%Examples\Magnetostatic"
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

