Rem  File SHW_TC.BAT displays the solutions calculated by the 
Rem  batch file Run_TC.BAT. 

Rem  Copyright 1998, by the University of California. 
Rem  Unauthorized commercial use is prohibited. 


cd "%SFDIR%Examples\CavityTuning\CCDTL"
START /W  " "  "%SFDIR%WSFplot" 1dtfc1.t35  3
START /W  " "  "%SFDIR%WSFplot" 2dthc1.t35  3


cd "%SFDIR%Examples\CavityTuning\CCL"
START /W  " "  "%SFDIR%WSFplot" 805ccl1.t35  3
START /W  " "  "%SFDIR%WSFplot" 805er1.t35   3
START /W  " "  "%SFDIR%WSFplot" 805wr1.t35   3
START /W  " "  "%SFDIR%WSFplot" 3ghz1.t35    3
START /W  " "  "%SFDIR%WSFplot" 425ew1.t35   3


cd "%SFDIR%Examples\CavityTuning\DTL"
START /W  " "  "%SFDIR%WSFplot"  rgdtl1.t35  3
START /W  " "  "%SFDIR%WSFplot"  test1.t35   3
START /W  " "  "%SFDIR%Tablplot" rgdtl.tbl   3


cd "%SFDIR%Examples\CavityTuning\EllipticalCavity"
START /W  " "  "%SFDIR%WSFplot"  82b1.t35    3
START /W  " "  "%SFDIR%WSFplot"  82be1.t35   3
START /W  " "  "%SFDIR%WSFplot"  8C82BE1.t35 3
START /W  " "  "%SFDIR%Tablplot" 8CFIELD.TBL 3


cd "%SFDIR%Examples\CavityTuning\MulticellDTL"
START /W  " "  "%SFDIR%WSFplot" test_1.t35  3


cd "%SFDIR%Examples\CavityTuning\RFQ"
START /W  " "  "%SFDIR%WSFplot" rfq_a1.t35 3
START /W  " "  "%SFDIR%WSFplot" rfq_b1.t35 3
START /W  " "  "%SFDIR%WSFplot" rfq_c1.t35 3
START /W  " "  "%SFDIR%WSFplot" rfq_d1.t35 3
START /W  " "  "%SFDIR%WSFplot" rfq_d2.t35 3
START /W  " "  "%SFDIR%WSFplot" rfq_d3.t35 3
START /W  " "  "%SFDIR%WSFplot" rfq_e1.t35 3
START /W  " "  "%SFDIR%WSFplot" rfq_f1.t35 3

START /W  " "  "%SFDIR%WSFplot" la1.t35 3
START /W  " "  "%SFDIR%WSFplot" la2.t35 3
START /W  " "  "%SFDIR%WSFplot" la3.t35 3
START /W  " "  "%SFDIR%WSFplot" la4.t35 3


cd "%SFDIR%Examples\CavityTuning\SideCouplingCell"
START /W  " "  "%SFDIR%WSFplot" 805cc1.t35 3
START /W  " "  "%SFDIR%WSFplot" lf1.t35    3

cd "%SFDIR%Examples\CavityTuning"

