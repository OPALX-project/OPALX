start /W  " "  "%SFDIR%autofish"  pill1
copy outaut.txt outaut1.txt
copy outfis.txt outfis1.txt
copy pill1.sfo  pill1a.sfo 

Rem Here we use different names for the SFO"  input files
Rem It is usually simpler to make sure the SEG file has the
Rem same name as the T35 solution file because then SFO"  can
Rem be started with either one of the files.  It will find
Rem the other file automatically.

start /W  " "  "%SFDIR%SFO"  pill1.t35 pill1b.seg
copy pill1.sfo  pill1b.sfo 
start /W  " "  "%SFDIR%SFO"  pill1.t35 pill1c.seg
copy pill1.sfo  pill1c.sfo 
start /W  " "  "%SFDIR%SFO"  pill1.t35 pill1d.seg
copy pill1.sfo  pill1d.sfo 
start /W  " "  "%SFDIR%SFO"  pill1.t35 pill1e.seg
copy pill1.sfo  pill1e.sfo 
if (%1)==(p) start /W  " "  "%SFDIR%WSFplot" pill1.t35 3

start /W  " "  "%SFDIR%autofish"  pill2
copy outaut.txt outaut2.txt
copy outfis.txt outfis2.txt
if (%1)==(p) start /W  " "  "%SFDIR%WSFplot" pill2.t35 3

start /W  " "  "%SFDIR%autofish"  modpill
if (%1)==(p) start /W  " "  "%SFDIR%WSFplot" modpill.t35 3

