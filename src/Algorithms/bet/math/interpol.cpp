/*
   Project: Beam Envelope Tracker (BET)
   Author:  Rene Bakker et al.
   Created: 09-03-2006

   Spline interpolation routines

   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
*/

#include "Utilities/OpalException.h"
#include "Algorithms/bet/math/interpol.h"

#include <cmath>
#include <vector>

/* Internal functions
   ========================================================================= */

/*
   Given arrays x[0..n-1] and y[0..n-1] containing a tabulated function,
   i.e., yi = f(xi), with x1 < x2 < .. . < xN, this routine returns an
   array y2[0..n-1] that contains the second derivatives of the
   interpolating function at the tabulated points x[i].
   Derrivate at bounderies set to zero !
*/
void spline(double x[], double y[], int n, double y2[]) {

    std::vector<double> u(n - 1);
    y2[0] = u[0] = 0.0;
    // :FIXME: why commented out? Impact of this change?
    //  y2[0] = -0.5;
    // u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0]));

    for(auto i = 1; i <= n - 2; i++) {
        auto sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
        auto p = sig * y2[i-1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        u[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]) - (y[i] - y[i-1]) / (x[i] - x[i-1]);
        u[i] = (6.0 * u[i] / (x[i+1] - x[i-1]) - sig * u[i-1]) / p;
    }

    double qn = 0.0;
    double un = 0.0;
    // :FIXME: why commented out?  Impact of this change?
    // qn=0.5;
    // un=(3.0/(x[n-1]-x[n-2]))*((y[n-1]-y[n-2])/(x[n-1]-x[n-2]));

    y2[n-1] = (un - qn * u[n-2]) / (qn * y2[n-2] + 1.0);
    for(auto k = n - 2; k >= 0; k--)
        y2[k] = y2[k] * y2[k+1] + u[k];
}

/*
   Given the arrays xa[0..n-1] and ya[0..n-1], which tabulate a function
   (with the xa[i] in order), and given the array y2a[0..n-1], which
   is the output from spline above, and given a value of x, this routine
   returns a cubic-spline interpolated value y, which is limited between
   the y-values of the adjacent points.
*/
double lsplint(double xa[], double ya[], double y2a[], int n, double x) {

    auto xb = x;
    if(xb < xa[0]) {
        //xb = xa[0];
        throw OpalException(
            "BetMath_interpol::lsplint()",
            "called with x < xMin.");
    }
    if(xb > xa[n-1]) {
        //xb = xa[n-1];
        throw OpalException(
            "BetMath_interpol::lsplint()",
            "called with x > xMax.");
    }

    auto klo = 0;
    auto khi = n - 1;
    while(khi - klo > 1) {
        auto k = (khi + klo) >> 1;
        if(xa[k] > xb)
            khi = k;
        else
            klo = k;
    }
    auto h = xa[khi] - xa[klo];
    if(h == 0.0)
        throw OpalException(
            "BetMath_interpol::lsplint()",
            "Bad xa input.");
    auto a = (xa[khi] - xb) / h;
    auto b = (xb - xa[klo]) / h;

    auto y = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
    if(ya[khi] > ya[klo]) {
        if(y > ya[khi])
            y = ya[khi];
        else if(y < ya[klo])
            y = ya[klo];
    } else {
        if(y > ya[klo])
            y = ya[klo];
        else if(y < ya[khi])
            y = ya[khi];
    }
    return y;
}

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
