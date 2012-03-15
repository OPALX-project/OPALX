Rem  File Show_ES.BAT displays the solutions calculated 
Rem  by the batch file Run_ES.BAT. 

Rem  Copyright 1998, by the University of California. 
Rem  Unauthorized commercial use is prohibited. 

cd "%SFDIR%Examples\Electrostatic\ChargeDensity
START /W  " "  "%SFDIR%WSFplot"   c1q.T35   3
START /W  " "  "%SFDIR%Tablplot"  c1q1.tbl 3
START /W  " "  "%SFDIR%WSFplot"   c1qc.T35   3
START /W  " "  "%SFDIR%Tablplot"  c1qc1.tbl 3
START /W  " "  "%SFDIR%WSFplot"   c1qmm.T35   3
START /W  " "  "%SFDIR%Tablplot"  c1qmm1.tbl 3
START /W  " "  "%SFDIR%WSFplot"   c4q.T35   3
START /W  " "  "%SFDIR%Tablplot"  c4q01.tbl 3
START /W  " "  "%SFDIR%WSFplot"   c4qtri.T35   3
START /W  " "  "%SFDIR%Tablplot"  c4qtri1.tbl 3
START /W  " "  "%SFDIR%WSFplot"   c4qpnt.T35   3
START /W  " "  "%SFDIR%Tablplot"  c4qpnt1.tbl 3


cd "%SFDIR%Examples\Electrostatic\CoaxCapacitor
START /W  " "  "%SFDIR%WSFplot"  caprz.t35    3
START /W  " "  "%SFDIR%WSFplot"  capxy.t35    3
START /W  " "  "%SFDIR%Tablplot" outpoi2.tbl  3


cd "%SFDIR%Examples\Electrostatic\CoaxCylinder
START /W  " "  "%SFDIR%WSFplot" coaxtest.t35 3


cd "%SFDIR%Examples\Electrostatic\Dielectrics
START /W  " "  "%SFDIR%WSFplot" 1epsilon.t35 3
START /W  " "  "%SFDIR%WSFplot" 2epsilon.t35 3


cd "%SFDIR%Examples\Electrostatic\IonSource
START /W  " "  "%SFDIR%WSFplot" spsource.t35 3


cd "%SFDIR%Examples\Electrostatic\Plates
START /W  " "  "%SFDIR%WSFplot"  Plate.t35    3
START /W  " "  "%SFDIR%Tablplot" plate1.tbl  3


cd "%SFDIR%Examples\Electrostatic\QuadLens
START /W  " "  "%SFDIR%WSFplot" esquad.t35   3
START /W  " "  "%SFDIR%WSFplot" esquad8.t35  3


cd "%SFDIR%Examples\Electrostatic\StoredEnergy
START /W  " "  "%SFDIR%WSFplot" ese1.t35     3
START /W  " "  "%SFDIR%WSFplot" ese2.t35     3


