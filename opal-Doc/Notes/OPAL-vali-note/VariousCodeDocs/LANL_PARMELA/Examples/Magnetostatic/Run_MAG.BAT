Rem  File Run_MAG.BAT runs Poisson Superfish in each of the
Rem  Magnetostatic subdirectories.  Run from a Windows tool
Rem  such as My Computer.  Running from DOS windows will not
Rem  work because of long filenames.

Rem  Copyright 1998, by the University of California.
Rem  Unauthorized commercial use is prohibited.

Rem  Remove "Rem" from the next line to prevent codes from opening display windows.
Rem set screenoutput=no


cd "%SFDIR%Examples\Magnetostatic\Dipoles"
START /W  " "  "%SFDIR%automesh"   dipole1
START /W  " "  "%SFDIR%pandira"    dipole1
copy outaut.txt outaut1.txt
copy outpan.txt outpan1.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" dipole1.t35 3

START /W  " "  "%SFDIR%automesh"   dipole2
START /W  " "  "%SFDIR%poisson"    dipole2
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" dipole2.t35 3


cd "%SFDIR%Examples\Magnetostatic\ForceOnWire"
START /W  " "  "%SFDIR%automesh"      1wire
rem START /W  " "  "%SFDIR%poisson"   1wire
START /W  " "  "%SFDIR%pandira"       1wire
START /W  " "  "%SFDIR%force"         1wire
copy outaut.txt outaut1.txt
rem  Use OUTPOI.TXT if you run Poisson  above
rem  copy outpoi.txt outpoi1.txt
copy outpan.txt outpan1.txt
copy outfor.txt outfor1.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" 1wire.t35 3

START /W  " "  "%SFDIR%automesh"      1w90
rem START /W  " "  "%SFDIR%poisson"   1w90
START /W  " "  "%SFDIR%pandira"       1w90
START /W  " "  "%SFDIR%force"         1w90
copy outaut.txt outaut2.txt
rem  Use OUTPOI.TXT if you run Poisson  above
rem  copy outpoi.txt outpoi2.txt
copy outpan.txt outpan2.txt
copy outfor.txt outfor2.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" 1w90.t35 3

START /W  " "  "%SFDIR%automesh"      1w45
rem START /W  " "  "%SFDIR%poisson"   1w45
START /W  " "  "%SFDIR%pandira"       1w45
START /W  " "  "%SFDIR%force"         1w45
copy outaut.txt outaut3.txt
rem  Use OUTPOI.TXT if you run Poisson above
rem  copy outpoi.txt outpoi3.txt
copy outpan.txt outpan3.txt
copy outfor.txt outfor3.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" 1w45.t35 3

START /W  " "  "%SFDIR%automesh"     2wire
START /W  " "  "%SFDIR%pandira"      2wire
START /W  " "  "%SFDIR%force"        2wire
copy outfor.txt outfor4.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" 2wire.t35 3

START /W  " "  "%SFDIR%automesh"   loop1
START /W  " "  "%SFDIR%poisson"    loop1
START /W  " "  "%SFDIR%force"      loop1
copy outfor.txt outfor5.txt
copy outpoi.txt outpoi5.txt
copy outfor.txt outfor5.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" loop1.t35 3

START /W  " "  "%SFDIR%automesh"   loop2
START /W  " "  "%SFDIR%poisson"    loop2
START /W  " "  "%SFDIR%force"      loop2
copy outfor.txt outfor6.txt
copy outpoi.txt outpoi6.txt
copy outfor.txt outfor6.txt
del outfor.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" loop2.t35 3


cd "%SFDIR%Examples\Magnetostatic\ForceTest"
START /W  " "  "%SFDIR%automesh"     iron1
rem START /W  " "  "%SFDIR%poisson"  iron1
START /W  " "  "%SFDIR%pandira"      iron1
START /W  " "  "%SFDIR%force"        iron1
copy outaut.txt outaut1.txt
rem  Use OUTPOI.TXT if you run Poisson above
rem  copy outpoi.txt outpoi1.txt
copy outpan.txt outpan1.txt
copy outfor.txt outfor1.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" iron1.t35 3

START /W  " "  "%SFDIR%automesh"     iron2
rem START /W  " "  "%SFDIR%poisson"  iron2
START /W  " "  "%SFDIR%pandira"      iron2
START /W  " "  "%SFDIR%force"        iron2
copy outaut.txt outaut2.txt
rem  Use OUTPOI.TXT if you run Poisson above
rem  copy outpoi.txt outpoi2.txt
copy outpan.txt outpan2.txt
copy outfor.txt outfor2.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" iron2.t35 3



cd "%SFDIR%Examples\Magnetostatic\H-Magnet"
START /W  " "  "%SFDIR%automesh"   htest1
START /W  " "  "%SFDIR%poisson"    htest1
copy outaut.txt outaut1.txt
copy outpoi.txt outpoi1.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" htest1.t35 3

START /W  " "  "%SFDIR%automesh"   htest2
START /W  " "  "%SFDIR%poisson"    htest2
copy outaut.txt outaut2.txt
copy outpoi.txt outpoi2.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" htest2.t35 3

START /W  " "  "%SFDIR%automesh"   hfull
START /W  " "  "%SFDIR%pandira"    hfull
START /W  " "  "%SFDIR%force"      hfull
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" hfull.t35 3



cd "%SFDIR%Examples\Magnetostatic\PMCylinders"
START /W  " "  "%SFDIR%automesh"   pmcylin
START /W  " "  "%SFDIR%pandira"    pmcylin
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" pmcylin.t35 3


cd "%SFDIR%Examples\Magnetostatic\PMDipoles"
START /W  " "  "%SFDIR%automesh"   pmdipole
START /W  " "  "%SFDIR%pandira"    pmdipole
START /W  " "  "%SFDIR%SF7"        pmdipole
copy outaut.txt outaut1.txt
copy outpan.txt outpan1.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"  pmdipole.t35 3
if (%1)==(p) START /W  " "  "%SFDIR%Tablplot" pmdipo2.tbl 3

START /W  " "  "%SFDIR%automesh"   pmdvar
START /W  " "  "%SFDIR%pandira"    pmdvar
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" pmvar.t35 3


cd "%SFDIR%Examples\Magnetostatic\PMQuads"
START /W  " "  "%SFDIR%automesh"   pmquad1
START /W  " "  "%SFDIR%pandira"    pmquad1
copy outaut.txt outaut1.txt
copy outpan.txt outpan1.txt
copy outpan.tbl outpan1.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" pmquad1.t35 3

START /W  " "  "%SFDIR%automesh"   pmquad2
START /W  " "  "%SFDIR%pandira"    pmquad2
copy outaut.txt outaut2.txt
copy outpan.txt outpan2.txt
copy outpan.tbl outpan2.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" pmquad2.t35 3

START /W  " "  "%SFDIR%automesh"   pmquad3
START /W  " "  "%SFDIR%pandira"    pmquad3
copy outaut.txt outaut3.txt
copy outpan.txt outpan3.txt
copy outpan.tbl outpan3.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" pmquad3.t35 3

del outaut.txt
del outpan.txt
del outpan.tbl



cd "%SFDIR%Examples\Magnetostatic\POAsample"
START /W  " "  "%SFDIR%automesh"   poa1
copy outaut.txt outaut1.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" poa1.t35 3

START /W  " "  "%SFDIR%automesh"   poa2
START /W  " "  "%SFDIR%poisson"    poa2
copy outaut.txt outaut2.txt
copy outpoi.txt outpoi2.txt
copy outpoi.tbl outpoi2.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" poa2.t35 3

del outaut.txt
del outpoi.txt
del outpoi.tbl



cd "%SFDIR%Examples\Magnetostatic\Quadrupoles"
START /W  " "  "%SFDIR%automesh"   qtest1
START /W  " "  "%SFDIR%poisson"    qtest1
copy outaut.txt outaut1.txt
copy outpoi.txt outpoi1.txt
copy outpoi.tbl outpoi1.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" qtest1.t35 3

START /W  " "  "%SFDIR%automesh"   qtest2
START /W  " "  "%SFDIR%poisson"    qtest2
copy outaut.txt outaut2.txt
copy outpoi.txt outpoi2.txt
copy outpoi.tbl outpoi2.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" qtest2.t35 3

del outaut.txt
del outpoi.txt
del outpoi.tbl


cd "%SFDIR%Examples\Magnetostatic\Septum"
START /W  " "  "%SFDIR%automesh"   septum
START /W  " "  "%SFDIR%poisson"    septum
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" septum.t35 3


cd "%SFDIR%Examples\Magnetostatic\Shielding"
START /W  " "  "%SFDIR%automesh"   vecp
START /W  " "  "%SFDIR%poisson"    vecp
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" vecp.t35 3


cd "%SFDIR%Examples\Magnetostatic\Solenoids"
START /W  " "  "%SFDIR%automesh"   soll1
START /W  " "  "%SFDIR%pandira"    soll1
copy outaut.txt outaut1.txt
copy outpan.txt outpan1.txt
copy outpan.tbl outpan1.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" soll1.t35 3

START /W  " "  "%SFDIR%automesh"   soll2
START /W  " "  "%SFDIR%poisson"    soll2
START /W  " "  "%SFDIR%force"      soll2
START /W  " "  "%SFDIR%SF7"        soll2
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" soll2.t35 3



cd "%SFDIR%Examples\Magnetostatic\StoredEnergy"
START /W  " "  "%SFDIR%automesh"   mse1
START /W  " "  "%SFDIR%pandira"    mse1
copy outaut.txt outaut1.txt
copy outpan.txt outpan1.txt
copy outpan.tbl outpan1.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" mse1.t35 3

START /W  " "  "%SFDIR%automesh"   mse2
START /W  " "  "%SFDIR%pandira"    mse2
copy outaut.txt outaut2.txt
copy outpan.txt outpan2.txt
copy outpan.tbl outpan2.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" mse2.t35 3

START /W  " "  "%SFDIR%automesh"   mse3
START /W  " "  "%SFDIR%pandira"    mse3
copy outaut.txt outaut3.txt
copy outpan.txt outpan3.txt
copy outpan.tbl outpan3.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" mse3.t35 3

del outaut.txt
del outpan.txt
del outpan.tbl


cd "%SFDIR%Examples\Magnetostatic"
