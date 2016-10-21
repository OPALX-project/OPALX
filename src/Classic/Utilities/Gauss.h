// ------------------------------------------------------------------------
// $RCSfile: Gauss.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Gaussian density function.
//
// ------------------------------------------------------------------------
// Class category: Utilities
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:38 $
// $Author: fci $
//
// ------------------------------------------------------------------------


/// Gaussian density function.
//  Arguments:
//  [DL]
//  [DT][b]sigx[/b]
//  [DD]the standard deviation in x-direction.
//  [DT][b]sigy[/b]
//  [DD]the standard deviation in y-direction.
//  [DT][b]x[/b]
//  [DD]the bias in x-direction.
//  [DT][b]y[/b]
//  [DD]the bias in y-direction.
//  [/DL]
//  Return value: Gaussian density[br][br]
//  exp(-(x/sigx)^2 - (y/sigy)^2) / (2 pi * sigx * sigy)

extern double Gauss(double sigx, double sigy, double dx, double dy);
