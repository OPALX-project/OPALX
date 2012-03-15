START /W  " "  "%SFDIR%automesh"   pmquad1
START /W  " "  "%SFDIR%pandira"    pmquad1
copy outaut.txt outaut1.txt
copy outpan.txt outpan1.txt
copy outpan.tbl outpan1.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" pmquad1.t35 3

START /W  " "  "%SFDIR%automesh"   pmquad2
START /W  " "  "%SFDIR%pandira"    pmquad2
copy outaut.txt outaut2.txt
copy outpan.txt outpan2.txt
copy outpan.tbl outpan2.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" pmquad2.t35 3

START /W  " "  "%SFDIR%automesh"   pmquad3
START /W  " "  "%SFDIR%pandira"    pmquad3
copy outaut.txt outaut3.txt
copy outpan.txt outpan3.txt
copy outpan.tbl outpan3.tbl
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" pmquad3.t35 3

del outaut.txt
del outpan.txt
del outpan.tbl

