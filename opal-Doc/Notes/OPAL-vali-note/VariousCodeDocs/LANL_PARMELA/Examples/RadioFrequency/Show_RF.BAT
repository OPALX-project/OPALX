Rem  File Show_RF.BAT displays the solutions calculated by the 
Rem  batch file Run_RF.BAT. 

Rem  Copyright 1998, by the University of California. 
Rem  Unauthorized commercial use is prohibited. 

cd "%SFDIR%Examples\RadioFrequency\CFISH"
start /W  " "  "%SFDIR%WSFplot"   dirod.t35     3
start /W  " "  "%SFDIR%WSFplot"   coaxwg.t35    3
start /W  " "  "%SFDIR%Tablplot"  coaxwg1.tbl  3
start /W  " "  "%SFDIR%WSFplot"   splitter.t35  3
start /W  " "  "%SFDIR%WSFplot"   qztube.t35    3
start /W  " "  "%SFDIR%Tablplot"  qztube1.tbl  3


cd "%SFDIR%Examples\RadioFrequency\FerriteCavity"
start /W  " "  "%SFDIR%WSFplot" ftest.t35  3 


cd "%SFDIR%Examples\RadioFrequency\FrequencyScan"
start /W  " "  "%SFDIR%WSFplot"   mode825.t35   3
start /W  " "  "%SFDIR%WSFplot"   mode763.t35   3
start /W  " "  "%SFDIR%Tablplot"  FishScan.tbl  3
start /W  " "  "%SFDIR%Tablplot"  Segs.tbl      3


cd "%SFDIR%Examples\RadioFrequency\PillboxCavities"
start /W  " "  "%SFDIR%WSFplot"  pill1.t35    3
start /W  " "  "%SFDIR%WSFplot"  pill2.t35    3
start /W  " "  "%SFDIR%WSFplot"  modpill.t35  3


cd "%SFDIR%Examples\RadioFrequency\QuarterWaveResonator"
start /W  " "  "%SFDIR%WSFplot"  qwave.t35  3


cd "%SFDIR%Examples\RadioFrequency\RFQCavity"
start /W  " "  "%SFDIR%WSFplot"  4thrfq.t35 3
start /W  " "  "%SFDIR%WSFplot"  8thrfq.t35 3


cd "%SFDIR%Examples\RadioFrequency\SphericalCavity"
start /W  " "  "%SFDIR%WSFplot"  Sphere.t35 3


cd "%SFDIR%Examples\RadioFrequency\Waveguides"
start /W  " "  "%SFDIR%WSFplot"  wr2300.t35 3
start /W  " "  "%SFDIR%WSFplot"  wghex.t35  3


cd "%SFDIR%Examples\RadioFrequency"