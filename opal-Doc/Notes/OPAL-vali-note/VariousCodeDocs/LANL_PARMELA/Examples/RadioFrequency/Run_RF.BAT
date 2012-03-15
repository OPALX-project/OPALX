Rem  File Run_RF.BAT runs Poisson Superfish in each of the
Rem  RadioFrequency subdirectories.  Run from a Windows tool
Rem  such as My Computer.  Running from DOS windows will not
Rem  work because of long filenames.

Rem  Copyright 1998, by the University of California.
Rem  Unauthorized commercial use is prohibited.

Rem  Remove "Rem" from the next line to prevent codes from opening display windows.
Rem set screenoutput=no

cd "%SFDIR%Examples\RadioFrequency\Cfish"
START /W  " "  "%SFDIR%automesh"  dirod
START /W  " "  "%SFDIR%cfish"     dirod
copy outaut.txt outaut1.txt
copy outfis.txt outfis1.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" dirod.t35 3

START /W  " "  "%SFDIR%automesh"  coaxwg
START /W  " "  "%SFDIR%cfish"     coaxwg
START /W  " "  "%SFDIR%SFO"       coaxwg
START /W  " "  "%SFDIR%SF7"       coaxwg
copy outaut.txt outaut2.txt
copy outfis.txt outfis2.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" coaxwg.t35 3

START /W  " "  "%SFDIR%automesh"  splitter
START /W  " "  "%SFDIR%cfish"     splitter
copy outaut.txt outaut3.txt
copy outfis.txt outfis3.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" splitter.t35 3

START /W  " "  "%SFDIR%automesh"  qztube
START /W  " "  "%SFDIR%cfish"     qztube
START /W  " "  "%SFDIR%SFO"       qztube
START /W  " "  "%SFDIR%SF7"       qztube
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" qztube.t35 3


cd "%SFDIR%Examples\RadioFrequency\FerriteCavity"
START /W  " "  "%SFDIR%automesh"  ftest
START /W  " "  "%SFDIR%fish"      ftest
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" ftest.t35 3


cd "%SFDIR%Examples\RadioFrequency\FrequencyScan"
Rem Delete any old plot file, then do frequency scan.
if exist fishscan.tbl del fishscan.tbl
START /W  " "  "%SFDIR%automesh"  scan
START /W  " "  "%SFDIR%fish"      scan
copy outfis.txt outfis1.txt
if (%1)==(p) START /W  " "  "%SFDIR%Tablplot" fishscan.tbl 3

Rem Run Autofish to find the higher-frequecy mode near 825 MHz.
START /W  " "  "%SFDIR%autofish"  mode825
copy outfis.txt outfis2.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" mode825.t35 3

Rem Run Autofish to find the lower-frequecy mode near 763 MHz.
START /W  " " "%SFDIR%autofish"  mode763
copy outfis.txt outfis3.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" mode763.t35 3

Rem Plot fields on segments 2 through 7 for the 763-MHz mode.
START /W  " "  "%SFDIR%segfield"   mode763
if (%1)==(p) START /W  " "  "%SFDIR%Tablplot" Segs.tbl 3


cd "%SFDIR%Examples\RadioFrequency\PillboxCavities"
START /W  " "  "%SFDIR%autofish"  pill1
copy outaut.txt outaut1.txt
copy outfis.txt outfis1.txt
copy pill1.sfo  pill1a.sfo 

Rem Here we use different names for the SFO  input files
Rem It is usually simpler to make sure the SEG file has the
Rem same name as the T35 solution file because then SFO"  can
Rem be STARTed with either one of the files.  It will find
Rem the other file automatically.

START /W  " "  "%SFDIR%SFO"  pill1.t35 pill1b.seg
copy pill1.sfo  pill1b.sfo 
START /W  " "  "%SFDIR%SFO"  pill1.t35 pill1c.seg
copy pill1.sfo  pill1c.sfo 
START /W  " "  "%SFDIR%SFO"  pill1.t35 pill1d.seg
copy pill1.sfo  pill1d.sfo 
START /W  " "  "%SFDIR%SFO"  pill1.t35 pill1e.seg
copy pill1.sfo  pill1e.sfo 
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" pill1.t35 3

START /W  " "  "%SFDIR%autofish"  pill2
copy outaut.txt outaut2.txt
copy outfis.txt outfis2.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" pill2.t35 3

START /W  " "  "%SFDIR%autofish"  modpill
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" modpill.t35 3


cd "%SFDIR%Examples\RadioFrequency\QuarterWaveResonator"
START /W  " "  "%SFDIR%autofish"  qwave
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" qwave.t35 3


cd "%SFDIR%Examples\RadioFrequency\RFQCavity"
START /W  " "  "%SFDIR%autofish"  4thrfq
copy outaut.txt outaut1.txt
copy outfis.txt outfis1.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" 4thrfq.t35 3

START /W  " "  "%SFDIR%autofish"  8thrfq
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" 8thrfq.t35 3


cd "%SFDIR%Examples\RadioFrequency\SphericalCavity"
START /W  " "  "%SFDIR%autofish"  sphere
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" sphere.t35 3


cd "%SFDIR%Examples\RadioFrequency\Waveguides"
START /W  " "  "%SFDIR%automesh"  wr2300
START /W  " "  "%SFDIR%cfish"     wr2300
copy outaut.txt outaut1.txt
copy outfis.txt outfis1.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" wr2300.t35 3

START /W  " "  "%SFDIR%automesh"  wghex
START /W  " "  "%SFDIR%cfish"     wghex
copy outaut.txt outaut2.txt
copy outfis.txt outfis2.txt
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" wghex.t35 3


cd "%SFDIR%Examples\RadioFrequency"
