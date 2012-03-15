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

