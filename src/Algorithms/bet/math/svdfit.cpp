/* svdfit.C
   SVD fitting routines

   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    ------------
   10/03/2006    Copied from "numerical recipies"                Rene Bakker

   Last Revision:
   $Id: svdfit.C 29 2007-04-14 17:03:18Z l_bakker $
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Algorithms/bet/BetError.h"
#include "Algorithms/bet/math/svdfit.h"

/* Internal functions
   ========================================================================= */

/* vector
   allocate an array of doubles [0..n-1] and check memory */
static double *vector(int n) {
    double *b;

    b = (double *) malloc(sizeof(double) * n);
    if(!b) {
        //writeBetError(errModeAll,errFatal,"Insufficient memory malloc %d bytes (svdfit.C)",sizeof(double)*n);
        writeBetError("Insufficient memory in malloc svdfit.C");
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
        //writeBetError(errModeAll,errFatal,"Allocation failure1 matrix(%d,%d) in svdfit.C\n",n1,n2);
        writeBetError("allocation failure 1 in svdfit.C");
    }

    for(i = 0; i < n1; i++) {
        m[i] = (double *) malloc(sizeof(double) * n2);
        if(!m[i]) {
            //writeBetError(errModeAll,errFatal,"Allocation failure2 matrix(%d,%d) in svdfit.C at i=%d\n",n1,n2,i);
            writeBetError("allocation failure 2 in svdfit.C");
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

/* svdksb()
   Solves A·X = B for a vector X, where A is specified by the arrays
   u[1..m][1..n], w[1..n], v[1..n][1..n] as returned by svdcmp.
   m and n are the dimensions of a, and will be equal for square matrices.
   b[1..m] is the input right-hand side. x[1..n] is the output solution vector.
   No input quantities are destroyed, so the routine may be called sequentially
   with different b's.
*/
static void svbksb
(
    double **u, double w[], double **v, int m, int n, double b[], double x[]) {
    int jj, j, i;
    double s, *tmp;

    tmp = vector(n + 1);
    for(j = 1; j <= n; j++) {
        s = 0.0;
        if(w[j]) {
            for(i = 1; i <= m; i++) s += u[i][j] * b[i];
            s /= w[j];
        }
        tmp[j] = s;
    }
    for(j = 1; j <= n; j++) {
        s = 0.0;
        for(jj = 1; jj <= n; jj++) s += v[j][jj] * tmp[jj];
        x[j] = s;
    }
    free(tmp);
} /* svdksb() */


/* svdcmp
   Given a matrix a[1..m][1..n], this routine computes its singular value
   decomposition, A =U·W·V^T. The matrix U replaces a on output. The diagonal
   matrix of singular values W is output as a vector w[1..n]. The matrix V
   (not the transpose V^T ) is output as v[1..n][1..n].
*/

static double maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double pythag(double a, double b) {
    double absa, absb;
    absa = fabs(a);
    absb = fabs(b);
    if(absa > absb) return absa * sqrt(1.0 + SQR(absb / absa));
    else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}

static void svdcmp(double **a, int m, int n, double w[], double **v) {
    int flag, i, its, j, jj, k, l, nm;
    double anorm, c, f, g, h, s, scale, x, y, z, *rv1;

    rv1 = vector(n + 1);
    g = scale = anorm = 0.0;
    for(i = 1; i <= n; i++) {
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if(i <= m) {
            for(k = i; k <= m; k++) scale += fabs(a[k][i]);
            if(scale) {
                for(k = i; k <= m; k++) {
                    a[k][i] /= scale;
                    s += a[k][i] * a[k][i];
                }
                f = a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = f - g;
                for(j = l; j <= n; j++) {
                    for(s = 0.0, k = i; k <= m; k++) s += a[k][i] * a[k][j];
                    f = s / h;
                    for(k = i; k <= m; k++) a[k][j] += f * a[k][i];
                }
                for(k = i; k <= m; k++) a[k][i] *= scale;
            }
        }
        w[i] = scale * g;
        g = s = scale = 0.0;
        if(i <= m && i != n) {
            for(k = l; k <= n; k++) scale += fabs(a[i][k]);
            if(scale) {
                for(k = l; k <= n; k++) {
                    a[i][k] /= scale;
                    s += a[i][k] * a[i][k];
                }
                f = a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = f - g;
                for(k = l; k <= n; k++) rv1[k] = a[i][k] / h;
                for(j = l; j <= m; j++) {
                    for(s = 0.0, k = l; k <= n; k++) s += a[j][k] * a[i][k];
                    for(k = l; k <= n; k++) a[j][k] += s * rv1[k];
                }
                for(k = l; k <= n; k++) a[i][k] *= scale;
            }
        }
        anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }
    for(i = n; i >= 1; i--) {
        if(i < n) {
            if(g) {
                for(j = l; j <= n; j++)
                    v[j][i] = (a[i][j] / a[i][l]) / g;
                for(j = l; j <= n; j++) {
                    for(s = 0.0, k = l; k <= n; k++) s += a[i][k] * v[k][j];
                    for(k = l; k <= n; k++) v[k][j] += s * v[k][i];
                }
            }
            for(j = l; j <= n; j++) v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
    for(i = IMIN(m, n); i >= 1; i--) {
        l = i + 1;
        g = w[i];
        for(j = l; j <= n; j++) a[i][j] = 0.0;
        if(g) {
            g = 1.0 / g;
            for(j = l; j <= n; j++) {
                for(s = 0.0, k = l; k <= m; k++) s += a[k][i] * a[k][j];
                f = (s / a[i][i]) * g;
                for(k = i; k <= m; k++) a[k][j] += f * a[k][i];
            }
            for(j = i; j <= m; j++) a[j][i] *= g;
        } else for(j = i; j <= m; j++) a[j][i] = 0.0;
        ++a[i][i];
    }
    for(k = n; k >= 1; k--) {
        for(its = 1; its <= 30; its++) {
            flag = 1;
            for(l = k; l >= 1; l--) {
                nm = l - 1;
                if((double)(fabs(rv1[l]) + anorm) == anorm) {
                    flag = 0;
                    break;
                }
                if((double)(fabs(w[nm]) + anorm) == anorm) break;
            }
            if(flag) {
                c = 0.0;
                s = 1.0;
                for(i = l; i <= k; i++) {
                    f = s * rv1[i];
                    rv1[i] = c * rv1[i];
                    if((double)(fabs(f) + anorm) == anorm) break;
                    g = w[i];
                    h = pythag(f, g);
                    w[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for(j = 1; j <= m; j++) {
                        y = a[j][nm];
                        z = a[j][i];
                        a[j][nm] = y * c + z * s;
                        a[j][i] = z * c - y * s;
                    }
                }
            }
            z = w[k];
            if(l == k) {
                if(z < 0.0) {
                    w[k] = -z;
                    for(j = 1; j <= n; j++) v[j][k] = -v[j][k];
                }
                break;
            }
            if(its == 30) {
                //writeBetError(errModeAll,errFatal,"svdfit(): no convergence in 30 svdcmp iterations");
                writeBetError("svdfit() no convergence in 30 svdcmp iterations");
            }
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
            c = s = 1.0;
            for(j = l; j <= nm; j++) {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for(jj = 1; jj <= n; jj++) {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = x * c + z * s;
                    v[jj][i] = z * c - x * s;
                }
                z = pythag(f, h);
                w[j] = z;
                if(z) {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for(jj = 1; jj <= m; jj++) {
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = y * c + z * s;
                    a[jj][i] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }
    free(rv1);
}

#undef FMAX
#undef IMIN
#undef SIGN

/* svdfit
   Given a set of data points x[1..ndata],y[1..ndata] with individual standard
   deviations sig[1..ndata], use X^2 minimization to determine the coefficients
   a[1..ma] of the fitting function y = Sum a[i]*afunci(x). Here we solve the
   fitting equations using singular value decomposition of the ndata by ma matrix,
   as in paragraph 2.6. Arrays u[1..ndata][1..ma], v[1..ma][1..ma], and w[1..ma]
   provide workspace on input; on output they define the singular value decomposition,
   and can be used to obtain the covariance matrix. The program returns values for the
   ma fit parameters a, and X^2, chisq. The user supplies a routine funcs(x,afunc,ma)
   that returns the ma basis functions evaluated at x = x in the array afunc[1..ma].
*/
#define TOL 1.0e-12

static void svdfit1
(
    double x[], double y[], double sig[], int ndata, double a[], int ma,
    double **u, double **v, double w[], double *chisq,
    void (*funcs)(double, double [], int)) {
    int j, i;
    double wmax, tmp, thresh, sum, *b, *afunc;

    b = vector(ndata + 1);
    afunc = vector(ma + 1);
    for(i = 1; i <= ndata; i++) {
        (*funcs)(x[i], afunc, ma);
        tmp = 1.0 / sig[i];
        for(j = 1; j <= ma; j++) u[i][j] = afunc[j] * tmp;
        b[i] = y[i] * tmp;
    }
    svdcmp(u, ndata, ma, w, v);
    wmax = 0.0;
    for(j = 1; j <= ma; j++)
        if(w[j] > wmax) wmax = w[j];
    thresh = TOL * wmax;
    for(j = 1; j <= ma; j++)
        if(w[j] < thresh) w[j] = 0.0;
    svbksb(u, w, v, ndata, ma, b, a);
    *chisq = 0.0;
    for(i = 1; i <= ndata; i++) {
        (*funcs)(x[i], afunc, ma);
        for(sum = 0.0, j = 1; j <= ma; j++) sum += a[j] * afunc[j];
        *chisq += (tmp = (y[i] - sum) / sig[i], tmp * tmp);
    }
    free(afunc);
    free(b);
}


/* svdvar()
   To evaluate the covariance matrix cvm[1..ma][1..ma] of the fit for ma
   parameters obtained by svdfit, call this routine with matrices
   v[1..ma][1..ma], w[1..ma] as returned from svdfit.
*/
void svdvar(double **v, int ma, double w[], double **cvm) {
    int k, j, i;
    double sum, *wti;

    wti = vector(ma + 1);
    for(i = 1; i <= ma; i++) {
        wti[i] = 0.0;
        if(w[i]) wti[i] = 1.0 / (w[i] * w[i]);
    }
    for(i = 1; i <= ma; i++) {
        for(j = 1; j <= i; j++) {
            for(sum = 0.0, k = 1; k <= ma; k++) sum += v[i][k] * v[j][k] * wti[k];
            cvm[j][i] = cvm[i][j] = sum;
        }
    }
    free(wti);
}

/* External functions
   ========================================================================= */

/* svdfit
   Given a set of data points x[0..ndata-1],y[0..ndata-1], use X^2 minimization
   to determine the coefficients a[0..ma-1] of the fitting function y = Sum a[i]*afunci(x).
   Here we solve the fitting equations using singular value decomposition of the
   ndata by ma matrix, as in paragraph 2.6.

   The program returns values for the  ma fit parameters a, the estimated error,
   and X^2, chisq. The user supplies a routine funcs(x,afunc,ma) that returns
   the ma basis functions evaluated at x = x in the array afunc[0..ma-1].
*/

static void (*theFuncs)(double, double *, int);

static void locFuncs(double x, double *p, int n) {
    theFuncs(x, p + 1, n);
}

void svdfit
(
    double *x, double *y, double *sig, int ndata, double *a, int ma,
    double *err, double *chisq,
    void (*funcs)(double, double *, int)) {
    int i;

    double
    *w   = vector(ndata + 1),
     **u  = matrix(ndata + 1, ma + 1),
       **v  = matrix(ndata + 1, ma + 1),
         **cvm = matrix(ma + 1, ma + 1);

    /* Array coversion [0..n-1] to [1..n]
       just lazy: instead of changing the routines above
       I've just changed the array index on calling */

    theFuncs = funcs;
    svdfit1(x - 1, y - 1, sig - 1, ndata, a - 1, ma, u, v, w, chisq, locFuncs);

    if(err) {
        svdvar(v, ma, w, cvm);

        for(i = 0; i < ma; i++) {
            err[i] = cvm[i+1][i+1];
        }
    }

    free_matrix(cvm);
    free_matrix(u);
    free_matrix(v);
    free(w);
} /* svdfit */

/* svdfit with sig[i] = 0.0, with i [0..ndata-1] */
void svdfit
(
    double *x, double *y, int ndata, double *a, int ma,
    double *err, double *chisq,
    void (*funcs)(double, double *, int)) {
    int i;
    double
    *sig = (double *) malloc(sizeof(double) * ndata);

    for(i = 0; i < ndata; i++) sig[i] = 1.0;

    svdfit(x, y, sig, ndata, a, ma, err, chisq, funcs);
    free(sig);
} /* svdfit */


/* svdfit with sig[i] = 0.0 and i [0..ndata-1]
   fit polynominal of order ma */

static void fpoly(double x, double p[], int np) {
    int j;

    p[0] = 1.0;
    for(j = 1; j < np; j++) p[j] = p[j-1] * x;
}


void svdfitP
(
    double *x, double *y, int ndata, double *a, int ma,
    double *err, double *chisq) {
    svdfit(x, y, ndata, a, ma, err, chisq, fpoly);
}

#undef TOL
