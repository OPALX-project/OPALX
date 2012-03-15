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
