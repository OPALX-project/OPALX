START /W  " "  "%SFDIR%automesh"  caprz
START /W  " "  "%SFDIR%poisson"   caprz
copy outaut.txt outaut1.txt
copy outpoi.txt outpoi1.txt
copy outpoi.tbl outpoi1.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" caprz.t35  3

START /W  " "  "%SFDIR%automesh"  capxy
START /W  " "  "%SFDIR%poisson"   capxy
copy outaut.txt outaut2.txt
copy outpoi.txt outpoi2.txt
copy outpoi.tbl outpoi2.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" capxy.t35  3
