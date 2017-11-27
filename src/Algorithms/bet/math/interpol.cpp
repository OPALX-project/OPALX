/* interpol.C

   Spline interpolation routines
   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   09-03-06      Created                                         Rene Bakker

   Last Revision:
   $Id: interpol.C 29 2007-04-14 17:03:18Z l_bakker $
*/

#include "Algorithms/bet/BetError.h"
#include "Algorithms/bet/math/interpol.h"

#include <cstdlib>
#include <cstddef>
#include <cmath>
#include <vector>

using namespace std;





/* Internal functions
   ========================================================================= */

#if 0
/* vector
   allocate an array of doubles [0..n-1] and check memory */
static double *vector(int n) {
    double *b;

    b = (double *) malloc(sizeof(double) * n);
    if(!b) {
        //writeBetError(errModeAll,errFatal,
        //    "Insufficient memory malloc %d bytes (interpol.C)",sizeof(double)*n);
        writeBetError("Insufficient memory (interpol.C)");
    }
    return b;
} /* vector */
#endif

/* Internal functions
   ========================================================================= */

/* spline()
   Given arrays x[0..n-1] and y[0..n-1] containing a tabulated function,
   i.e., yi = f(xi), with x1 < x2 < .. . < xN, this routine returns an
   array y2[0..n-1] that contains the second derivatives of the
   interpolating function at the tabulated points x[i].
   Derrivate at bounderies set to zero !
*/
void spline(double x[], double y[], int n, double y2[]) {
    int i, k;
    double p, qn, sig, un;

    std::vector<double> u(n - 1);
    y2[0] = u[0] = 0.0;
    //  y2[0] = -0.5;
    // u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0]));

    for(i = 1; i <= n - 2; i++) {
        sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
        p = sig * y2[i-1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        u[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]) - (y[i] - y[i-1]) / (x[i] - x[i-1]);
        u[i] = (6.0 * u[i] / (x[i+1] - x[i-1]) - sig * u[i-1]) / p;
    }

    qn = un = 0.0;
    // qn=0.5;
    // un=(3.0/(x[n-1]-x[n-2]))*((y[n-1]-y[n-2])/(x[n-1]-x[n-2]));

    y2[n-1] = (un - qn * u[n-2]) / (qn * y2[n-2] + 1.0);
    for(k = n - 2; k >= 0; k--)
        y2[k] = y2[k] * y2[k+1] + u[k];
} /* spline() */


/* splint()
   Given the arrays xa[0..n-1] and ya[0..n-1], which tabulate a function
   (with the xa[i] in order), and given the array y2a[0..n-1], which
   is the output from spline above, and given a value of x, this routine
   returns a cubic-spline interpolated value y.
*/
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y) {
    int klo, khi, k;
    double xb, h, b, a;

    xb = x;
    if(xb < xa[0]) {
        xb = xa[0];
        //writeBetError(errModeAll,errGeneral,"slint called with x < xMin: %lf < %lf",x,xa[0]);
        writeBetError("slint called with x < xMin");
    }
    if(xb > xa[n-1]) {
        xb = xa[n-1];
        //writeBetError(errModeAll,errGeneral,"slint called with x > xMax: %lf > %lf",x,xa[n-1]);
        writeBetError("slint called with x > xMax");
    }

    klo = 0;
    khi = n - 1;
    while(khi - klo > 1) {
        k = (khi + klo) >> 1;
        if(xa[k] > xb) khi = k;
        else klo = k;
    }
    h = xa[khi] - xa[klo];
    if(h == 0.0)
        //writeBetError(errModeAll,errFatal,
        //      "Bad xa input to routine splint: x=%lf [%lf,%lf]",
        //         x,xa[0],xa[n-1]);
        writeBetError("Bad xa input to routine splint");
    a = (xa[khi] - xb) / h;
    b = (xb - xa[klo]) / h;
    *y = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
} /* splint() */


/* splint()
   Given the arrays xa[0..n-1] and ya[0..n-1], which tabulate a function
   (with the xa[i] in order), and given the array y2a[0..n-1], which
   is the output from spline above, and given a value of x, this routine
   returns a cubic-spline interpolated value y, which is limited between
   the y-values of the adjacent points.
*/
void lsplint(double xa[], double ya[], double y2a[], int n, double x, double *y) {
    int klo, khi, k;
    double xb, h, b, a, yy;

    xb = x;
    if(xb < xa[0]) {
        xb = xa[0];
        //writeBetError(errModeAll,errGeneral,"lsplint called with x < xMin: %lf < %lf",x,xa[0]);
        writeBetError("lslplint x < xMin");
    }
    if(xb > xa[n-1]) {
        xb = xa[n-1];
        //writeBetError(errModeAll,errGeneral,"lsplint called with x > xMax: %lf > %lf",x,xa[n-1]);
        writeBetError("lsplint x > xMax");
    }

    klo = 0;
    khi = n - 1;
    while(khi - klo > 1) {
        k = (khi + klo) >> 1;
        if(xa[k] > xb) khi = k;
        else klo = k;
    }
    h = xa[khi] - xa[klo];
    if(h == 0.0)
        //writeBetError(errModeAll,errFatal,
        //         "Bad xa input to routine lsplint: x=%lf [%lf,%lf]",
        //         x,xa[0],xa[n-1]);
        writeBetError("Bad xa input to routine lsplint");
    a = (xa[khi] - xb) / h;
    b = (xb - xa[klo]) / h;

    yy = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
    if(ya[khi] > ya[klo]) {
        if(yy > ya[khi]) yy = ya[khi];
        else if(yy < ya[klo]) yy = ya[klo];
    } else {
        if(yy > ya[klo]) yy = ya[klo];
        else if(yy < ya[khi]) yy = ya[khi];
    }
    *y = yy;
} /* splint() */


