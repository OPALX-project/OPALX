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

