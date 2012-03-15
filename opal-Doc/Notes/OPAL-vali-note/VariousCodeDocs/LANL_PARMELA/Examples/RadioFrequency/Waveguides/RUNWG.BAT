start /W  " "  "%SFDIR%automesh"  wr2300
start /W  " "  "%SFDIR%cfish"     wr2300
copy outaut.txt outaut1.txt
copy outfis.txt outfis1.txt
if (%1)==(p) start /W  " "  "%SFDIR%WSFplot" wr2300.t35 3

start /W  " "  "%SFDIR%automesh"  wghex
start /W  " "  "%SFDIR%cfish"     wghex
copy outaut.txt outaut2.txt
copy outfis.txt outfis2.txt
if (%1)==(p) start /W  " "  "%SFDIR%WSFplot" wghex.t35 3

