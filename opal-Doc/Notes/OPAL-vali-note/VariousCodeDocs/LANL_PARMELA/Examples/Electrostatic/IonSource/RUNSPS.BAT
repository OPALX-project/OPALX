START /W  " "  "%SFDIR%automesh"   spsource
START /W  " "  "%SFDIR%poisson"    spsource
START /W  " "  "%SFDIR%SFO"        spsource
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot"  spsource.t35  3
