/* root.C
   root finding routines

   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    ------------
   09-03-2006    Copied from "numerical recipies"                Rene Bakker

   Last Revision:
   $Id: root.C 29 2007-04-14 17:03:18Z l_bakker $
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Algorithms/bet/math/root.h"
#include "Algorithms/bet/BetError.h"

/*
   ======================================================================= */
#define MAXIT 100


/* findRoot()
   Originally rtsafe()
   Using a combination of Newton-Raphson and bisection, find the root of a
   function bracketed between x1 and x2. The root, returned as the function
   value rtsafe, will be refined until its accuracy is known within ±xacc.
   funcd is a user-supplied routine that returns both the function value
   and the first derivative of the function. */
double findRoot
(
    void (*funcd)(double, double *, double *), // function to find root
    double x1, double x2,                      // boundaries
    double xacc) {                             // accuracy
    int j;
    double df, dx, dxold, f, fh, fl;
    double temp, xh, xl, rts;

    (*funcd)(x1, &fl, &df);
    (*funcd)(x2, &fh, &df);
    if((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
        //writeBetError(errModeAll,errFatal,"Root must be bracketed in findRoot()");
        writeBetError("Root must be bracketed in frindRoot");
    }
    if(fl == 0.0) return x1;
    if(fh == 0.0) return x2;
    if(fl < 0.0) {
        xl = x1;
        xh = x2;
    } else {
        xh = x1;
        xl = x2;
    }
    rts = 0.5 * (x1 + x2);
    dxold = fabs(x2 - x1);
    dx = dxold;
    (*funcd)(rts, &f, &df);
    for(j = 1; j <= MAXIT; j++) {
        if((((rts - xh)*df - f) * ((rts - xl)*df - f) >= 0.0)
           || (fabs(2.0 * f) > fabs(dxold * df))) {
            dxold = dx;
            dx = 0.5 * (xh - xl);
            rts = xl + dx;
            if(xl == rts) return rts;
        } else {
            dxold = dx;
            dx = f / df;
            temp = rts;
            rts -= dx;
            if(temp == rts) return rts;
        }
        if(fabs(dx) < xacc) return rts;
        (*funcd)(rts, &f, &df);
        if(f < 0.0)
            xl = rts;
        else
            xh = rts;
    }
    //writeBetError(errModeAll,errFatal,"Maximum number of iterations (%d) exceeded in findRoot()",MAXIT);
    writeBetError("Maximum number of iterations exeeded in findRoot");

    return 0.0;
} /* findroot() */

#undef MAXIT

