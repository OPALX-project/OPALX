START /W  " "  "%SFDIR%automesh"   dipole1
START /W  " "  "%SFDIR%pandira"    dipole1
copy outaut.txt outaut1.txt
copy outpan.txt outpan1.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" dipole1.t35 3

START /W  " "  "%SFDIR%automesh"   dipole2
START /W  " "  "%SFDIR%poisson"    dipole2
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" dipole2.t35 3
