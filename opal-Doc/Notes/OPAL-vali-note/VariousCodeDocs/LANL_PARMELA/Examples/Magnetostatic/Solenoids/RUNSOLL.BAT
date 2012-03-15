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

