Rem  File Run_ES.BAT runs Poisson Superfish in each of the
Rem  Electrostatic subdirectories.    Run from a Windows tool
Rem  such as My Computer.  Running from DOS windows will not
Rem  work because of long filenames.

Rem  Copyright 1998, by the University of California.
Rem  Unauthorized commercial use is prohibited.

Rem  Remove "Rem" from the next line to prevent codes from opening display windows.
Rem set screenoutput=no

cd "%SFDIR%Examples\Electrostatic\ChargeDensity"
del /q  c1q*.tbl
del /q  c4q*.tbl
start /W  " "  "%SFDIR%automesh"  c1q.am
START /W  " "  "%SFDIR%poisson"   c1q.t35
copy outaut.txt outaut1.txt
copy outpoi.txt outpoi1.txt
START /W  " "  "%SFDIR%SF7"       c1q.t35
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"   c1q.T35   3
if (%1)==(p) START /W  " "  "%SFDIR%Tablplot"  c1q01.tbl 3

start /W  " "  "%SFDIR%automesh"  c1qc.am
START /W  " "  "%SFDIR%poisson"   c1qc.t35
copy outaut.txt outaut2.txt
copy outpoi.txt outpoi2.txt
START /W  " "  "%SFDIR%SF7"       c1qc.t35
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"   c1qc.T35   3
if (%1)==(p) START /W  " "  "%SFDIR%Tablplot"  c1qc01.tbl 3

start /W  " "  "%SFDIR%automesh"  c1qmm.am
START /W  " "  "%SFDIR%poisson"   c1qmm.t35
copy outaut.txt outaut3.txt
copy outpoi.txt outpoi3.txt
START /W  " "  "%SFDIR%SF7"       c1qmm.t35
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"   c1qmm.T35   3
if (%1)==(p) START /W  " "  "%SFDIR%Tablplot"  c1qmm1.tbl 3


start /W  " "  "%SFDIR%automesh"  c4q.am
START /W  " "  "%SFDIR%poisson"   c4q.t35
copy outaut.txt outaut4.txt
copy outpoi.txt outpoi4.txt
START /W  " "  "%SFDIR%SF7"       c4q.t35
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"   c4q.T35   3
if (%1)==(p) START /W  " "  "%SFDIR%Tablplot"  c4q1.tbl 3


start /W  " "  "%SFDIR%automesh"  c4qtri.am
START /W  " "  "%SFDIR%poisson"   c4qtri.t35
copy outaut.txt outaut5.txt
copy outpoi.txt outpoi5.txt
START /W  " "  "%SFDIR%SF7"       c4qtri.t35
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"   c4qtri.T35   3
if (%1)==(p) START /W  " "  "%SFDIR%Tablplot"  c4qtri1.tbl 3


start /W  " "  "%SFDIR%automesh"  c4qpnt.am
START /W  " "  "%SFDIR%poisson"   c4qpnt.t35
copy outaut.txt outaut6.txt
copy outpoi.txt outpoi6.txt
START /W  " "  "%SFDIR%SF7"       c4qpnt.t35
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"   c4qpnt.T35   3
if (%1)==(p) START /W  " "  "%SFDIR%Tablplot"  c4qpnt1.tbl 3


cd "%SFDIR%Examples\Electrostatic\CoaxCapacitor"
START /W  " "  "%SFDIR%automesh"  caprz
START /W  " "  "%SFDIR%poisson"   caprz
copy outaut.txt outaut1.txt
copy outpoi.txt outpoi1.txt
copy outpoi.tbl outpoi1.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" caprz.t35  3

START /W  " "  "%SFDIR%automesh"  capxy
START /W  " "  "%SFDIR%poisson"   capxy
copy outaut.txt outaut2.txt
copy outpoi.txt outpoi2.txt
copy outpoi.tbl outpoi2.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" capxy.t35  3


cd "%SFDIR%Examples\Electrostatic\CoaxCylinder"
START /W  " "  "%SFDIR%automesh"  coaxtest
START /W  " "  "%SFDIR%poisson"   coaxtest
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" coaxtest.t35  3


cd "%SFDIR%Examples\Electrostatic\Dielectrics"
START /W  " "  "%SFDIR%automesh"  1epsilon
START /W  " "  "%SFDIR%poisson"   1epsilon
copy outaut.txt outaut1.txt
copy outpoi.txt outpoi1.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" 1epsilon.t35  3

START /W  " "  "%SFDIR%automesh"  2epsilon
START /W  " "  "%SFDIR%poisson"   2epsilon
copy outaut.txt outaut2.txt
copy outpoi.txt outpoi2.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"  2epsilon.t35  3


cd "%SFDIR%Examples\Electrostatic\IonSource"
START /W  " "  "%SFDIR%automesh"   spsource
START /W  " "  "%SFDIR%poisson"    spsource
START /W  " "  "%SFDIR%SFO"        spsource
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"  spsource.t35  3


cd "%SFDIR%Examples\Electrostatic\Plates"
START /W  " "  "%SFDIR%automesh"   plate
START /W  " "  "%SFDIR%poisson"    plate
START /W  " "  "%SFDIR%SF7"        plate
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"  plate.t35  3
if (%1)==(p) START /W  " "  "%SFDIR%Tablplot" plate1.tbl  3


cd "%SFDIR%Examples\Electrostatic\QuadLens"
START /W  " "  "%SFDIR%automesh"  esquad
START /W  " "  "%SFDIR%poisson"   esquad
START /W  " "  "%SFDIR%SFO"       esquad.t35
copy outaut.txt outaut1.txt
copy outpoi.txt outpoi1.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" esquad.t35  3

START /W  " "  "%SFDIR%automesh"  esquad8
START /W  " "  "%SFDIR%poisson"   esquad8
START /W  " "  "%SFDIR%SFO"       esquad8
copy outaut.txt outaut2.txt
copy outpoi.txt outpoi2.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot esquad8.t35  3


cd "%SFDIR%Examples\Electrostatic\StoredEnergy"
START /W  " "  "%SFDIR%automesh"  ese1
START /W  " "  "%SFDIR%pandira"   ese1
copy outaut.txt outaut1.txt
copy outpan.txt outpan1.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"  ese1.t35  3

START /W  " "  "%SFDIR%automesh"  ese2
START /W  " "  "%SFDIR%pandira"   ese2
copy outaut.txt outaut2.txt
copy outpan.txt outpan2.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"  ese2.t35  3


cd "%SFDIR%Examples\Electrostatic"

