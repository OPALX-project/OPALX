START /W  " "  "%SFDIR%automesh"  esquad
START /W  " "  "%SFDIR%poisson"   esquad
START /W  " "  "%SFDIR%SFO"       esquad.t35
copy outaut.txt outaut1.txt
copy outpoi.txt outpoi1.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" esquad.t35  3

START /W  " "  "%SFDIR%automesh"  esquad8
START /W  " "  "%SFDIR%poisson"   esquad8
START /W  " "  "%SFDIR%SFO"       esquad8
copy outaut.txt outaut2.txt
copy outpoi.txt outpoi2.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot esquad8.t35  3
