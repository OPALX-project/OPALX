Rem  File Show_MAG.BAT displays the solutions calculated by the
Rem  batch file Run_MAG.BAT.

Rem  Copyright 1998, by the University of California.
Rem  Unauthorized commercial use is prohibited.

cd "%SFDIR%Examples\Magnetostatic\Dipoles
START /W  " "  "%SFDIR%WSFplot"  dipole1.t35 3
START /W  " "  "%SFDIR%Tablplot" outpan.tbl  3
START /W  " "  "%SFDIR%WSFplot"  dipole2.t35 3


cd "%SFDIR%Examples\Magnetostatic\ForceOnWire
START /W  " "  "%SFDIR%WSFplot"  1wire.t35 3
START /W  " "  "%SFDIR%WSFplot"  1w90.t35  3
START /W  " "  "%SFDIR%WSFplot"  1w45.t35  3
START /W  " "  "%SFDIR%WSFplot"  2wire.t35 3
START /W  " "  "%SFDIR%WSFplot"  loop1.t35 3
START /W  " "  "%SFDIR%WSFplot"  loop2.t35 3


cd "%SFDIR%Examples\Magnetostatic\ForceTest
START /W  " "  "%SFDIR%WSFplot" iron1.t35  3
START /W  " "  "%SFDIR%WSFplot" iron2.t35  3


cd "%SFDIR%Examples\Magnetostatic\H-Magnet
START /W  " "  "%SFDIR%WSFplot"  htest1.t35 3
START /W  " "  "%SFDIR%Tablplot" outpoi.tbl 3
START /W  " "  "%SFDIR%WSFplot"  htest2.t35 3
START /W  " "  "%SFDIR%WSFplot"  hfull.t35  3
START /W  " "  "%SFDIR%Tablplot" outpan.tbl 3


cd "%SFDIR%Examples\Magnetostatic\PMCylinders
START /W  " "  "%SFDIR%WSFplot" PMCYLIN.T35 3
START /W  " "  "%SFDIR%Tablplot" outpan.tbl 3


cd "%SFDIR%Examples\Magnetostatic\PMDipoles
START /W  " "  "%SFDIR%WSFplot"  pmdipole.t35 3
START /W  " "  "%SFDIR%Tablplot" pmdipo2.tbl 3
START /W  " "  "%SFDIR%WSFplot"  pmdvar.t35   3
START /W  " "  "%SFDIR%Tablplot" outpan.tbl   3


cd "%SFDIR%Examples\Magnetostatic\PMQuads
START /W  " "  "%SFDIR%WSFplot"  pmquad1.t35 3
START /W  " "  "%SFDIR%Tablplot" outpan1.tbl 3
START /W  " "  "%SFDIR%WSFplot"  pmquad2.t35 3
START /W  " "  "%SFDIR%Tablplot" outpan2.tbl 3
START /W  " "  "%SFDIR%WSFplot"  pmquad3.t35 3
START /W  " "  "%SFDIR%Tablplot" outpan3.tbl 3



cd "%SFDIR%Examples\Magnetostatic\POAsample
START /W  " "  "%SFDIR%WSFplot"  poa2.t35  3
START /W  " "  "%SFDIR%Tablplot" outpoi2.tbl 3


cd "%SFDIR%Examples\Magnetostatic\Quadrupoles
START /W  " "  "%SFDIR%WSFplot"  qtest1.t35  3
START /W  " "  "%SFDIR%Tablplot" outpoi1.tbl 3
START /W  " "  "%SFDIR%WSFplot"  qtest2.t35  3
START /W  " "  "%SFDIR%Tablplot" outpoi2.tbl 3


cd "%SFDIR%Examples\Magnetostatic\Septum
START /W  " "  "%SFDIR%WSFplot" septum.t35  3
START /W  " "  "%SFDIR%Tablplot" outpoi.tbl 3


cd "%SFDIR%Examples\Magnetostatic\Shielding
START /W  " "  "%SFDIR%WSFplot" vecp.t35    3
START /W  " "  "%SFDIR%Tablplot" outpoi.tbl 3


cd "%SFDIR%Examples\Magnetostatic\Solenoids
START /W  " "  "%SFDIR%WSFplot"  soll1.t35   3
START /W  " "  "%SFDIR%Tablplot" outpan1.tbl 3
START /W  " "  "%SFDIR%WSFplot"  soll2.t35   3
START /W  " "  "%SFDIR%Tablplot" soll21.tbl 3
START /W  " "  "%SFDIR%Quikplot" bfld.qkp    3


cd "%SFDIR%Examples\Magnetostatic\StoredEnergy
START /W  " "  "%SFDIR%WSFplot"  mse1.t35    3
START /W  " "  "%SFDIR%Tablplot" outpan1.tbl 3
START /W  " "  "%SFDIR%WSFplot"  mse2.t35    3
START /W  " "  "%SFDIR%Tablplot" outpan2.tbl 3
START /W  " "  "%SFDIR%WSFplot"  mse3.t35    3
START /W  " "  "%SFDIR%Tablplot" outpan3.tbl 3

