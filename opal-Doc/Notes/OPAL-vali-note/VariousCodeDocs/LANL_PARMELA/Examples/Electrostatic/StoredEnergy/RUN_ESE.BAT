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
