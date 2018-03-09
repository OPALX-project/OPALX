/* rk.C

   Runge-Kutta with adaptive stepsize control
   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   07-03-06      Created                                         Rene Bakker

   Last Revision:
   $Id: rk.C 29 2007-04-14 17:03:18Z l_bakker $
*/

#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "Algorithms/bet/BetError.h"
#include "Algorithms/bet/math/rk.h"

/* Internal functions
   ========================================================================= */
static double maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

/* vector
   allocate an array of doubles [0..n-1] and check memory */
static double *vector(int n) {
    double *b;

    b = (double *) malloc(sizeof(double) * n);
    if(!b) {
        //writeBetError(1,errModeAll,errFatal,"Insufficient memory malloc %d bytes (rk.C)",sizeof(double)*n);
        writeBetError("Insufficient memory malloc in rk.C");
    }
    return b;
} /* vector */

static double **matrix(int n1, int n2)
/* allocate a double matrix with subscript range [0..n1-1][0..n2-1] */
{
    int    i;
    double **m;

    /* allocate pointers to rows */
    m = (double **) malloc(sizeof(double *) * (n1 + 1));
    if(!m) {
        //  writeBetError(1,errModeAll,errFatal,"Allocation failure1 matrix(%d,%d) in rk.C\n",n1,n2);
        writeBetError("Allocation failure in matrix rk.C");
    }

    for(i = 0; i < n1; i++) {
        m[i] = (double *) malloc(sizeof(double) * n2);
        if(!m[i]) {
            // writeBetError(1,errModeAll,errFatal,"Allocation failure2 matrix(%d,%d) in rk.C at i=%d\n",n1,n2,i);
            writeBetError("allocation failure in matrix rk.C");
        }
    }
    m[n1] = NULL;
    return m;
} /* matrix */

static void free_matrix(double **m)
/* free a double matrix allocated by matrix() */
{
    int i = 0;

    while(m[i] != NULL) {
        free(m[i++]);
    }
    free(m);
} /* free_matrix */


/* Global variables + init for odeint
   ================================== */
static int kmax = 0; // do not store intermediate results for inspection
static int kount;
static double
*xp  = NULL,       // must only be set if kmax > 0
 **yp = NULL, // must only be set if kmax > 0
   dxsav;

static void rkSetBuffer(int nvar) {
    if(nvar < 1) {
        //writeBetError(1,errModeAll,errFatal,"odeint nvar (%d) < 1 (impossible)",nvar);
        writeBetError("odeint nvar < 1");
    }
    if(kmax > 0) {
        if(xp) free(xp);
        xp = vector(kmax);

        if(yp) free_matrix(yp);
        yp = matrix(nvar, kmax);
    }
} /* rkSetBuffer */

/* rkck
   Given values for n variables y[0..n-1] and their derivatives dydx[0..n-1]
   known at x, use the fifth-order Cash-Karp Runge-Kutta method
   to advance the solution over an interval h and return the incremented
   variables as yout[0..n-1].

   Also return an estimate of the local truncation error in yout
   using the embedded fourth-order method.

   The user supplies the routine derivs(x,y,dydx),
   which returns derivatives dydx at x.
*/
void rkck
(
    double y[],
    double dydx[],
    int    n,
    double x,
    double h,
    double yout[],
    double yerr[],
    void (*derivs)(double, double [], double [])) {
    int i;
    static double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875, b21 = 0.2,
                  b31 = 3.0 / 40.0, b32 = 9.0 / 40.0, b41 = 0.3, b42 = -0.9, b43 = 1.2,
                  b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0, b54 = 35.0 / 27.0,
                  b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0, b63 = 575.0 / 13824.0,
                  b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0, c1 = 37.0 / 378.0,
                  c3 = 250.0 / 621.0, c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0,
                  dc5 = -277.0 / 14336.0;
    double dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0,
           dc4 = c4 - 13525.0 / 55296.0, dc6 = c6 - 0.25;
    double *ak2, *ak3, *ak4, *ak5, *ak6, *ytemp;

    ak2  = vector(n);
    ak3  = vector(n);
    ak4  = vector(n);
    ak5  = vector(n);
    ak6  = vector(n);
    ytemp = vector(n);

    for(i = 0; i < n; i++)
        ytemp[i] = y[i] + b21 * h * dydx[i];
    (*derivs)(x + a2 * h, ytemp, ak2);
    for(i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i]);
    (*derivs)(x + a3 * h, ytemp, ak3);
    for(i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
    (*derivs)(x + a4 * h, ytemp, ak4);
    for(i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);
    (*derivs)(x + a5 * h, ytemp, ak5);
    for(i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);
    (*derivs)(x + a6 * h, ytemp, ak6);
    for(i = 0; i < n; i++)
        yout[i] = y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
    for(i = 0; i < n; i++)
        yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);

    free(ytemp);
    free(ak6);
    free(ak5);
    free(ak4);
    free(ak3);
    free(ak2);
} /* rkck */


#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
// The value ERRCON equals (5/SAFETY) raised to the power (1/PGROW), see use below.

/* rkqs
   Fifth-order Runge-Kutta step with monitoring of local truncation error
   to ensure accuracy and adjust stepsize.
   Input are:
   - the dependent variable vector y[0..n-1]
   - its derivative dydx[0..n-1]
   at the starting value of the independent variable x.

   Also input are:
   - the stepsize to be attempted htry
   - the required accuracy eps
   - the vector yscal[0..n-1] against which the error is scaled.

   RETURNS 1 on error and 0 otherwise

   On output, y and x are replaced by their new values,
   hdid is the stepsize that was actually accomplished,
   and hnext is the estimated next stepsize.

   derivs is the user-supplied routine that computes the
   right-hand side derivatives.
*/

static int rkqs
(
    double y[], double dydx[], int n, double *x, double htry, double eps,
    double yscal[], double *hdid, double *hnext,
    void (*derivs)(double, double [], double [])) {
    int i, nLoop;
    double errmax, h, xnew, *yerr, *ytemp;

    yerr = vector(n);
    ytemp = vector(n);
    h = htry;
    nLoop = 0;
    for(;;) {
        ++nLoop;
        rkck(y, dydx, n, *x, h, ytemp, yerr, derivs);
        errmax = 0.0;
        for(i = 0; i < n; i++) errmax = FMAX(errmax, fabs(yerr[i] / yscal[i]));
        errmax /= eps;
        if(errmax > 1.0) {
            h = SAFETY * h * pow(errmax, PSHRNK);
            if(h < 0.1 * h) h *= 0.1;
            xnew = (*x) + h;
            if(xnew == *x) {
                //writeBetError(1,errModeAll,errWarning,
                //     "stepsize underflow in rkqs (rk.C) h = %le, n = %d",
                //     h,nLoop);
                writeBetError("stepsize underflow rkgs");
                free(ytemp);
                free(yerr);
                return 1;
            }
            continue;
        } else {
            if(errmax > ERRCON) *hnext = SAFETY * h * pow(errmax, PGROW);
            else *hnext = 5.0 * h;
            *x += (*hdid = h);
            for(i = 0; i < n; i++) y[i] = ytemp[i];
            break;
        }
    }
    free(ytemp);
    free(yerr);
    return 0;
} /* rkqs */

#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef NRANSI

#undef FMAX

/* External functions
   ========================================================================= */

/* rk4()
   Given values for the variables y[0..n-1] use the fourth-order Runge-Kutta
   method to advance the solution over an interval h and return the
   incremented variables in y[0..n-1]. The user supplies the routine
   derivs(x,y,dydx), which returns derivatives dydx at x.t
*/
void rk4(double y[], int n, double x, double h,
         void (*derivs)(double, double [], double [])) {
    int i;
    double xh, hh, h6, *dydx, *dym, *dyt, *yt;

    dydx = vector(n);
    dym = vector(n);
    dyt = vector(n);
    yt = vector(n);

    (*derivs)(x, y, dydx);
    hh = h * 0.5;
    h6 = h / 6.0;
    xh = x + hh;
    for(i = 0; i < n; i++) yt[i] = y[i] + hh * dydx[i];
    (*derivs)(xh, yt, dyt);
    for(i = 0; i < n; i++) yt[i] = y[i] + hh * dyt[i];
    (*derivs)(xh, yt, dym);
    for(i = 0; i < n; i++) {
        yt[i] = y[i] + h * dym[i];
        dym[i] += dyt[i];
    }
    (*derivs)(x + h, yt, dyt);
    for(i = 0; i < n; i++)
        y[i] += h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
    free(yt);
    free(dyt);
    free(dym);
    free(dydx);
} /* rk4() */


#define MAXSTP 10000
#define TINY 1.0e-30

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


/* odeint
   Runge-Kutta driver with adaptive stepsize control.

   Integrate starting values ystart[0..nvar-1]
   from x1 to x2 with accuracy eps,

   storing intermediate results in global variables if kmax > 0
   (see above).
   h1 should be set as a guessed first stepsize,
   hmin as the minimum allowed stepsize (can be zero).

   On output:
   - nok and nbad are the number of good and bad (but retried and fixed)
     steps taken
   - ystart is replaced by values at the end of the integration interval.
   - RETURNS 0 on success and 1 on error

   derivs is the user-supplied routine for calculating the right-hand side
   derivative.

   rkqs is the name of the stepper routine (see above).
*/
int odeint(
    double ystart[],  // initial condition
    int    nvar,      // number of equations
    double x1,
    double x2,
    double eps,
    double h1,
    double hmin,
    int    *nok,
    int    *nbad,
    void (*derivs)(double, double [], double [])) {
    int nstp, i;
    double xsav = 0.0, x, hnext, hdid, h;
    double *yscal, *y, *dydx;

    if(kmax > 0) {  // set buffer
        rkSetBuffer(nvar);
    }

    yscal = vector(nvar);
    y     = vector(nvar);
    dydx  = vector(nvar);
    x = x1;
    h = SIGN(h1, x2 - x1);
    *nok = (*nbad) = kount = 0;
    for(i = 0; i < nvar; i++) y[i] = ystart[i];
    if(kmax > 0) xsav = x - dxsav * 2.0;
    for(nstp = 0; nstp < MAXSTP; nstp++) {
        (*derivs)(x, y, dydx);
        for(i = 0; i < nvar; i++) {
            yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + TINY;
        }
        if(kmax > 0 && kount < kmax - 1 && fabs(x - xsav) > fabs(dxsav)) {
            xp[++kount] = x;
            for(i = 0; i < nvar; i++) yp[i][kount] = y[i];
            xsav = x;
        }
        if((x + h - x2) * (x + h - x1) > 0.0) h = x2 - x;
        if(rkqs(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, derivs) == 1) {
            // error message already written
            free(dydx);
            free(y);
            free(yscal);
            return 1;
        }
        if(hdid == h) ++(*nok);
        else ++(*nbad);
        if((x - x2) * (x2 - x1) >= 0.0) {
            for(i = 0; i < nvar; i++) ystart[i] = y[i];
            if(kmax) {
                xp[++kount] = x;
                for(i = 0; i < nvar; i++) yp[i][kount] = y[i];
            }
            free(dydx);
            free(y);
            free(yscal);
            return 0;
        }
        if(fabs(hnext) <= hmin) {
            // writeBetError(1,errModeAll,errWarning,"Step size too small in odeint (%le < %le)",
            //     hnext,hmin);
            writeBetError("Step size too small in odeint");
            return 1;
        }
        h = hnext;
    }
    //writeBetError(1,errModeAll,errWarning,"Too many steps in routine odeint (%d)",MAXSTP);
    writeBetError("Too many steps in routine odeint");
    return 1;
} /* odeint */


void rkActivateBuffer(int max) {
    if(max < 1) {
        if(xp) free(xp);
        xp = NULL;
        if(yp) free_matrix(yp);
        yp = NULL;

        kmax = 0;
    } else {
        kmax = max;
    }
} /* setActiveBuffer */

void rkPrintBuffer(FILE *f) {
    int i, j;

    fprintf(f, "odeint output buffer, kmax = %d/%d\n", kount, kmax);

    fprintf(f, "%-4s %-10s", "k", "xp");
    j = 0;
    while(yp[j]) fprintf(f, "    [%2d]   ", j++);
    fprintf(f, "\n");

    for(i = 0; i < kount; i++) {
        fprintf(f, "%4d %10.4f", i, xp[i]);
        j = 0;
        while(yp[j]) fprintf(f, " %10.4f", yp[j++][i]);
        fprintf(f, "\n");
    }
} /* rkPrintBuffer */

#undef MAXSTP
#undef TINY
#undef NRANSI

#undef SIGN


