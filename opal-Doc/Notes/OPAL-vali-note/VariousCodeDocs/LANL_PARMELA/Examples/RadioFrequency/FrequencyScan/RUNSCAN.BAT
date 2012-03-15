Rem Delete any old plot file, then do frequency scan.
if exist fishscan.tbl del fishscan.tbl
start /W  " "  "%SFDIR%automesh"  scan
start /W  " "  "%SFDIR%fish"      scan
copy outfis.txt outfis1.txt
if (%1)==(p) start /W  " "  "%SFDIR%Tablplot" fishscan.tbl 3

Rem Run Autofish to find the higher-frequecy mode near 825 MHz.
start /W  " "  "%SFDIR%autofish"  mode825
copy outfis.txt outfis2.txt
if (%1)==(p) start /W  " "  "%SFDIR%WSFplot" mode825.t35 3

Rem Run Autofish to find the lower-frequecy mode near 763 MHz.
START /W  " " "%SFDIR%autofish"  mode763
copy outfis.txt outfis3.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" mode763.t35 3

Rem Plot fields on segments 2 through 7 for the 763-MHz mode.
start /W  " "  "%SFDIR%segfield"   mode763
if (%1)==(p) start /W  " "  "%SFDIR%Tablplot" Segs.tbl 3

