#ifndef OPAL_TpsWerrf_HH
#define OPAL_TpsWerrf_HH 1

// ------------------------------------------------------------------------
// $RCSfile: TpsWerrf.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Function: double TpsWerrf(double)
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:48 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/FTps.h"


/// Complex error function.
//  The argument is given by two DA variables.
//  The algorithms is based on:
//  [center]
//  Walter Gautschi, Efficient Computation of the Complex Error Function,
//  [br]
//  SIAM J. Numer. Anal., Vol 7, No. 1, March 1970, pp. 187-198.
//  [/center]
//  The argument is [tt]xx + i*yy[/tt], the result is [tt]wx + i*wy[/tt].

extern void TpsWerrf(const FTps<double, 6> &xx, const FTps<double, 6> &yy,
                     FTps<double, 6> &wx, FTps<double, 6> &wy);

#endif
