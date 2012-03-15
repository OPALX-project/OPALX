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
