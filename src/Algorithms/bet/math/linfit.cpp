/* linfit.C
   linear fitting routines

   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    ------------
   30/04/1992    Copied from "numerical recipies"                Rene Bakker

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Algorithms/bet/BetError.h"
#include "Algorithms/bet/math/linfit.h"

#include <math.h>
/* Internal functions
   ========================================================================= */

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

/* gcf()
 */
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

/* gammln()
   Returns the ln-value of the gamma function ln[G(xx)] for xx > 0.
*/
static double gammln(double xx) {
    double x, y, tmp, ser;
    static double cof[6] = {76.18009172947146, -86.50532032941677,
                            24.01409824083091, -1.231739572450155,
                            0.1208650973866179e-2, -0.5395239384953e-5
                           };
    int j;

    y = x = xx;
    tmp  = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser  = 1.000000000190015;
    for(j = 0; j <= 5; j++) ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / x);
} /* gammln */


/* gcf()

   Returns the incomplete gamma function Q(a, x) evaluated by its
   continued fraction representation as gammcf. Also returns Gamma(a) as
   gln.
*/

static void gcf(double *gammcf, double a, double x, double *gln) {
    int i;
    double an, b, c, d, del, h;

    *gln = gammln(a);
    b = x + 1.0 - a;
    c = 1.0 / FPMIN;
    d = 1.0 / b;
    h = d;
    for(i = 1; i <= ITMAX; i++) {
        an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if(fabs(d) < FPMIN) d = FPMIN;
        c = b + an / c;
        if(fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if(fabs(del - 1.0) < EPS) break;
    }
    if(i > ITMAX) {
        // writeBetError(errModeAll,errFatal,
        //       "linfit-->gfc(): a (%le) too large, ITMAX (%d) too small",
        //       a,ITMAX);
        writeBetError();
    }
    *gammcf = exp(-x + a * log(x) - (*gln)) * h;
} /* gcf */


/* gser() Returns the incomplete gamma function P(a,x) evaluated by
   its series representation as gamser.  Also returns ln Gamma(a) as gln.
*/
static void gser(double *gamser, double a, double x, double *gln) {
    int n;
    double sum, del, ap;

    *gln = gammln(a);
    if(x <= 0.0) {
        if(x < 0.0) {
            //writeBetError(errModeAll,errFatal,
            //         "linfit-->gser(): x (%lf) less than 0",
            //         x);
            writeBetError();
        }
        *gamser = 0.0;
        return;
    } else {
        ap = a;
        del = sum = 1.0 / a;
        for(n = 1; n <= ITMAX; n++) {
            ++ap;
            del *= x / ap;
            sum += del;
            if(fabs(del) < fabs(sum)*EPS) {
                *gamser = sum * exp(-x + a * log(x) - (*gln));
                return;
            }
        }
        //writeBetError(errModeAll,errFatal,
        //       "linfit-->gser(): a (%le) too large, ITMAX (%d) too small",
        //       a,ITMAX);
        writeBetError();
        return;
    }
} /* gser */

#undef ITMAX
#undef EPS
#undef FPMIN


/* gammq()
   Returns the incomplete gamma function Q(a, x) = 1 - P(a, x). */

static double gammq(double a, double x) {
    double gamser, gammcf, gln;

    if(x < 0.0 || a <= 0.0) {
        //writeBetError(errModeAll,errFatal,
        //       "lintfit-->gammq(): Invalid arguments in routine gammq");
        writeBetError();
    }
    if(x < (a + 1.0)) {
        gser(&gamser, a, x, &gln);
        return 1.0 - gamser;
    } else {
        gcf(&gammcf, a, x, &gln);
        return gammcf;
    }
} /* gammq */

/* external functions
   ========================================================================= */

/* linfit() Given a set of data points x[0..ndata-1],y[0..ndata-1] with
   individual standard deviations sig[0..ndata-1], fit them to a
   straight line y = a + bx by minimizing chi2. Returned are a,b and
   their respective probable uncertainties siga and sigb, the
   chi-square chi2, and the goodness-of-fit probability q (that the
   fit would have chi2 this large or larger). If mwt=0 on input, then
   the standard deviations are assumed to be unavailable: q is
   returned as 1.0 and the normalization of chi2 is to unit standard
   deviation on all points.
*/

void linfit(double x[], double y[], int ndata,
            double sig[], int mwt, double *a, double *b, double *siga,
            double *sigb, double *chi2, double *q) {

    int
    i;
    double
    wt, t, sxoss, sx = 0.0, sy = 0.0, st2 = 0.0, ss, sigdat;

    *b = 0.0;
    if(mwt) {
        ss = 0.0;
        for(i = 0; i < ndata; i++) {
            wt = 1.0 / SQR(sig[i]);
            ss += wt;
            sx += x[i] * wt;
            sy += y[i] * wt;
        }
    } else {
        for(i = 0; i < ndata; i++) {
            sx += x[i];
            sy += y[i];
        }
        ss = ndata;
    }
    sxoss = sx / ss;
    if(mwt) {
        for(i = 0; i < ndata; i++) {
            t = (x[i] - sxoss) / sig[i];
            st2 += t * t;
            *b += t * y[i] / sig[i];
        }
    } else {
        for(i = 0; i < ndata; i++) {
            t = x[i] - sxoss;
            st2 += t * t;
            *b += t * y[i];
        }
    }
    *b /= st2;
    *a = (sy - sx * (*b)) / ss;
    *siga = sqrt((1.0 + sx * sx / (ss * st2)) / ss);
    *sigb = sqrt(1.0 / st2);
    *chi2 = 0.0;
    if(mwt == 0) {
        for(i = 0; i < ndata; i++)
            *chi2 += SQR(y[i] - (*a) - (*b) * x[i]);
        *q = 1.0;
        sigdat = sqrt((*chi2) / (ndata - 2));
        *siga *= sigdat;
        *sigb *= sigdat;
    } else {
        for(i = 0; i < ndata; i++)
            *chi2 += SQR((y[i] - (*a) - (*b) * x[i]) / sig[i]);
        *q = gammq(0.5 * (ndata - 2), 0.5 * (*chi2));
    }
} /* linfit */

/* linfit()
   Given a set of data points x[0..ndata-1],y[0..ndata-1] with, fit
   them to a straight line y = a + bx by minimizing chi2. Returned are
   a,b and their respective probable uncertainties siga and sigb, and
   the chi-square chi2.
*/
void linfit(double x[], double y[], int ndata,
            double *a, double *b, double *siga,
            double *sigb, double *chi2) {

    int
    i;
    double
    t, sxoss, sx, sy, st2, sigdat;

    *b = 0.0;
    sx = 0.0;
    sy = 0.0;
    for(i = 0; i < ndata; i++) {
        sx += x[i];
        sy += y[i];
    }

    sxoss = sx / ndata;

    st2 = 0.0;
    for(i = 0; i < ndata; i++) {
        t = x[i] - sxoss;
        st2 += t * t;
        *b += t * y[i];
    }
    *b /= st2;
    *a = (sy - sx * (*b)) / ndata;
    *siga = sqrt((1.0 + sx * sx / (ndata * st2)) / ndata);
    *sigb = sqrt(1.0 / st2);
    *chi2 = 0.0;

    for(i = 0; i < ndata; i++)
        *chi2 += SQR(y[i] - (*a) - (*b) * x[i]);

    sigdat = sqrt((*chi2) / (ndata - 2));
    *siga *= sigdat;
    *sigb *= sigdat;
} /* linfit */


