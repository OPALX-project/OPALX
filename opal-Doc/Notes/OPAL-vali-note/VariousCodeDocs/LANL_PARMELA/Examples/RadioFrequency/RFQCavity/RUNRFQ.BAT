start /W  " "  "%SFDIR%autofish"  4thrfq
copy outaut.txt outaut1.txt
copy outfis.txt outfis1.txt
if (%1)==(p) start /W  " "  "%SFDIR%WSFplot" 4thrfq.t35 3

start /W  " "  "%SFDIR%autofish"  8thrfq
if (%1)==(p) start /W  " "  "%SFDIR%WSFplot" 8thrfq.t35 3
