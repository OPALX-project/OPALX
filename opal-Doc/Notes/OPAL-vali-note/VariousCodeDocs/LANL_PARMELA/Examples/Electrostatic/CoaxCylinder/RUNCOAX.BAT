START /W  " "  "%SFDIR%automesh"  coaxtest
START /W  " "  "%SFDIR%poisson"   coaxtest
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" coaxtest.t35  3
