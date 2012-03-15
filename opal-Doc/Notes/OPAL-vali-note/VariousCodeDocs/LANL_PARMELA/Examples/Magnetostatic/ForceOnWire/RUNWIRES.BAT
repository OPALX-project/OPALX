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

