Rem  This procedure deletes files created by Run_TC.BAT in the
Rem  CavityTuning directories.  Run from a Windows tool such as
Rem  My Computer.  Running from DOS windows will not work because
Rem  of long filenames.

Rem  Copyright 1998, by the University of California.
Rem  Unauthorized commercial use is prohibited.

cd "%SFDIR%Examples\CavityTuning"
del /s /q  *.inf
del /s /q  *.t35
del /s /q  out*.*
del /s /q  *.log
del /s /q  *.am
del /s /q  *.seg
del /s /q  *.sfo
del /s /q  *.pmi
del /s /q  *.hgl
del /s /q  tape*.*
del /s /q  *.prf
del /s /q  *.tbl
del /s /q  *.qkp
del /s /q  *.t7
del /s /q  *.rdp
del /s /q  wsfprf.txt
del /s /q  ?plotprf.txt
del /s /q  fishscan.*

del /s /q  ccdtl\?????0?.cdt

del /s /q  ccl\????0?.ccl
del /s /q  ccl\?????0?.ccl
del /s /q  ccl\??????0?.ccl

del /s /q  dtl\?????0?.dtl
del /s /q  dtl\??????0?.dtl
del /s /q  dtl\??????~?.dtl
del /s /q  dtl\*.in7

del /s /q EllipticalCavity\???0?.ell
del /s /q EllipticalCavity\????0?.ell
del /s /q EllipticalCavity\*.af
del /s /q EllipticalCavity\filein.ctl

del /s /q MulticellDTL\??????~?.mdt

del /s /q  rfq\??????0?.rfq
del /s /q  rfq\??????~?.rfq

del /s /q  SideCouplingCell\??0?.scc
del /s /q  SideCouplingCell\?????0?.scc
del /s /q  SideCouplingCell\*.0??

