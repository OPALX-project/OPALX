/* integration.C
   integration routines using Rombergs method

   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    ------------
   10/03/2006    Copied from "numerical recipies"                Rene Bakker

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Algorithms/bet/BetError.h"
#include "Algorithms/bet/math/integrate.h"

/* internal functions
   ====================================================================== */

/* vector
   allocate an array of doubles [0..n-1] and check memory */
static double *vector(int n) {
    double *b;

    b = (double *) malloc(sizeof(double) * n);
    if(!b) {
        writeBetError("Insufficient memory malloc bytes (rk.C)"); //,sizeof(double)*n);
    }
    return b;
} /* vector */

#define FUNC(x) ((*func)(x))

/* tranzd()
   This routine computes the nth stage of refinement of an extended
   trapezoidal rule. func is input as a pointer to the function to
   be integrated between limits a and b, also input. When called with
   n=1, the routine returns the crudest estimate of integral a->b f(x)dx.
   Subsequent calls with n=2,3,... (in that sequential order) will improve
   the accuracy by adding 2n-2 additional interior points.
*/
static double trapzd(double(*func)(double), double a, double b, int n) {
    double x, tnm, sum, del;
    static double s;
    int it, j;

    if(n == 1) {
        return (s = 0.5 * (b - a) * (FUNC(a) + FUNC(b)));
    } else {
        for(it = 1, j = 1; j < n - 1; j++) it <<= 1;
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for(sum = 0.0, j = 1; j <= it; j++, x += del) sum += FUNC(x);
        s = 0.5 * (s + (b - a) * sum / tnm);
        return s;
    }
} /* trapzd() */

#undef FUNC

/* polint
   Given arrays xa[1..n] and ya[1..n], and given a value x,
   this routine returns a value y, and an error estimate dy.
   If P(x) is the polynomial of degree N - 1 such that
   P(xa[i]) = ya[i], i = 1,..n, then the returned value y = P(x).
 */
static void polint
(
    double xa[], double ya[], int n,
    double x, double *y, double *dy) {
    int i, m, ns = 1;
    double den, dif, dift, ho, hp, w;
    double *c, *d;

    dif = fabs(x - xa[1]);
    c = vector(n + 1);
    d = vector(n + 1);
    for(i = 1; i <= n; i++) {
        if((dift = fabs(x - xa[i])) < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    *y = ya[ns--];
    for(m = 1; m < n; m++) {
        for(i = 1; i <= n - m; i++) {
            ho = xa[i] - x;
            hp = xa[i+m] - x;
            w = c[i+1] - d[i];
            if((den = ho - hp) == 0.0) {
                int j;

                fprintf(stderr, "Polint: n=%d, i=%d\n", n, i);
                for(j = 0; j < n; j++) fprintf(stderr, "%5d %20.12e %20.12e\n", j, xa[j], ya[j]);
                writeBetError("Error in routine polint (singular point)");
            }
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        *y += (*dy = (2 * ns < (n - m) ? c[ns+1] : d[ns--]));
    }
    free(d);
    free(c);
} /* polint() */

/* external functions
   ======================================================================= */

#define JMAX 30
#define JMAXP (JMAX+1)
#define K 5

/* qromb()
   Returns the integral of the function func from a to b.
   Integration is performed by Romberg's method of order 2K,
   where, e.g., K=2 is Simpson's rule.
*/
double qromb
(
    double(*func)(double),
    double a, double b,
    double eps) {
    double aeps, ss, dss;
    double s[JMAXP+1], h[JMAXP+1];
    int j;

    aeps = fabs(eps);
    h[1] = 1.0;
    for(j = 1; j <= JMAX; j++) {
        s[j] = trapzd(func, a, b, j);
        if(j >= K) {
            polint(&h[j-K], &s[j-K], K, 0.0, &ss, &dss);
            if(fabs(dss) < aeps * fabs(ss)) return ss;
        }
        s[j+1] = s[j];
        h[j+1] = 0.25 * h[j];
    }
    //writeBetError(errModeAll,errFatal,"Too many steps in routine qromb (%d) [%le/%le]",
    //         JMAX,fabs(dss),fabs(ss));
    writeBetError("Too many steps in routine gromb");
    return 0.0;
}
#undef JMAX
#undef JMAXP
#undef K
