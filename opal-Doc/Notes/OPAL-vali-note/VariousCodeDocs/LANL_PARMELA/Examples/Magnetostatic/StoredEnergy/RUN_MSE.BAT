START /W  " "  "%SFDIR%automesh"   mse1
START /W  " "  "%SFDIR%pandira"    mse1
copy outaut.txt outaut1.txt
copy outpan.txt outpan1.txt
copy outpan.tbl outpan1.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" mse1.t35 3

START /W  " "  "%SFDIR%automesh"   mse2
START /W  " "  "%SFDIR%pandira"    mse2
copy outaut.txt outaut2.txt
copy outpan.txt outpan2.txt
copy outpan.tbl outpan2.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" mse2.t35 3

START /W  " "  "%SFDIR%automesh"   mse3
START /W  " "  "%SFDIR%pandira"    mse3
copy outaut.txt outaut3.txt
copy outpan.txt outpan3.txt
copy outpan.tbl outpan3.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" mse3.t35 3

del outaut.txt
del outpan.txt
del outpan.tbl

