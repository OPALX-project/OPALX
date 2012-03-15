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

