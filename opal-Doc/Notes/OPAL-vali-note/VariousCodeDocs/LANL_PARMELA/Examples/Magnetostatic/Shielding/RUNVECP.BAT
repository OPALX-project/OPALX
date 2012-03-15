START /W  " "  "%SFDIR%automesh"   vecp
START /W  " "  "%SFDIR%poisson"    vecp
if (%1)==(p) START /W  " "  "%SFDIR%WSFplot" vecp.t35 3
