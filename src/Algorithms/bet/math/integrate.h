/*
   Project: Beam Envelope Tracker (BET)
   Author:  Rene Bakker et al.
   Created: 07-03-2006

   integration routines using Rombergs method

   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
*/

#ifndef _BET_INTEGRATE_H
#define _BET_INTEGRATE_H

/* qromb()
   Returns the integral of the function func from a to b
   with error tolerance eps.
   Integration is performed by Romberg's method. The order
   is set in in the program source code (default 5).
*/
double qromb(
    double( *)(double), // function
    double, double,     // boundaries
    double = 1.0e-4);   // error (eps)

#endif

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
