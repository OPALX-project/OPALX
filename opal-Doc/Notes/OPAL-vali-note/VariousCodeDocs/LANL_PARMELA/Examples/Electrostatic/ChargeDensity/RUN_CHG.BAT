del /q  c1q*.tbl
del /q  c4q*.tbl
start /W  " "  "%SFDIR%automesh"  c1q.am
START /W  " "  "%SFDIR%poisson"   c1q.t35
copy outaut.txt outaut1.txt
copy outpoi.txt outpoi1.txt
START /W  " "  "%SFDIR%SF7"       c1q.t35
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"   c1q.T35   3
if (%1)==(p) START /W  " "  "%SFDIR%Tablplot"  c1q1.tbl 3

start /W  " "  "%SFDIR%automesh"  c1qc.am
START /W  " "  "%SFDIR%poisson"   c1qc.t35
copy outaut.txt outaut2.txt
copy outpoi.txt outpoi2.txt
START /W  " "  "%SFDIR%SF7"       c1qc.t35
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"   c1qc.T35   3
if (%1)==(p) START /W  " "  "%SFDIR%Tablplot"  c1qc1.tbl 3

start /W  " "  "%SFDIR%automesh"  c1qmm.am
START /W  " "  "%SFDIR%poisson"   c1qmm.t35
copy outaut.txt outaut3.txt
copy outpoi.txt outpoi3.txt
START /W  " "  "%SFDIR%SF7"       c1qmm.t35
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"   c1qmm.T35   3
if (%1)==(p) START /W  " "  "%SFDIR%Tablplot"  c1qmm1.tbl 3


start /W  " "  "%SFDIR%automesh"  c4q.am
START /W  " "  "%SFDIR%poisson"   c4q.t35
copy outaut.txt outaut4.txt
copy outpoi.txt outpoi4.txt
START /W  " "  "%SFDIR%SF7"       c4q.t35
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"   c4q.T35   3
if (%1)==(p) START /W  " "  "%SFDIR%Tablplot"  c4q1.tbl 3


start /W  " "  "%SFDIR%automesh"  c4qtri.am
START /W  " "  "%SFDIR%poisson"   c4qtri.t35
copy outaut.txt outaut5.txt
copy outpoi.txt outpoi5.txt
START /W  " "  "%SFDIR%SF7"       c4qtri.t35
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"   c4qtri.T35   3
if (%1)==(p) START /W  " "  "%SFDIR%Tablplot"  c4qtri1.tbl 3


start /W  " "  "%SFDIR%automesh"  c4qpnt.am
START /W  " "  "%SFDIR%poisson"   c4qpnt.t35
copy outaut.txt outaut6.txt
copy outpoi.txt outpoi6.txt
START /W  " "  "%SFDIR%SF7"       c4qpnt.t35
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"   c4qpnt.T35   3
if (%1)==(p) START /W  " "  "%SFDIR%Tablplot"  c4qpnt1.tbl 3
