/*
   Project: Beam Envelope Tracker (BET)
   Author:  Rene Bakker et al.
   Created: 07-03-2006

   integration routines using Rombergs method

   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
*/

#include <cmath>

#include "Algorithms/bet/math/integrate.h"
#include "Utilities/OpalException.h"

/* internal functions
   ====================================================================== */

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

    static double s;

    if(n == 1) {
        s = 0.5 * (b - a) * (FUNC(a) + FUNC(b));
    } else {
        auto it = 1;
        for(auto j = 1; j < n - 1; j++)
            it <<= 1;
        auto tnm = double(it);
        auto del = (b - a) / tnm;
        auto x = a + 0.5 * del;
        auto sum = 0.0;
        for(auto j = 1; j <= it; j++, x += del)
            sum += FUNC(x);
        s = 0.5 * (s + (b - a) * sum / tnm);
    }
    return s;
}

/* polint
   Given arrays xa[1..n] and ya[1..n], and given a value x,
   this routine returns a value y, and an error estimate dy.
   If P(x) is the polynomial of degree N - 1 such that
   P(xa[i]) = ya[i], i = 1,..n, then the returned value y = P(x).
 */
static void polint (double xa[], double ya[], int n, double x, double &y, double &dy) {

    auto ns = 1;
    double c[n + 1];
    double d[n + 1];
    auto dif = std::fabs(x - xa[1]);
    for(auto i = 1; i <= n; i++) {
        auto dift = std::fabs(x - xa[i]); 
        if(dift < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    y = ya[ns--];
    for(auto m = 1; m < n; m++) {
        for(auto i = 1; i <= n - m; i++) {
            auto ho = xa[i] - x;
            auto hp = xa[i+m] - x;
            auto w = c[i+1] - d[i];
            auto den = ho - hp; 
            if(den == 0.0) {
                fprintf(stderr, "Polint: n=%d, i=%d\n", n, i);
                for(auto j = 0; j < n; j++)
                    fprintf(stderr, "%5d %20.12e %20.12e\n", j, xa[j], ya[j]);
                throw OpalException(
                    "BetMath_integrate::polint()",
                    "singular point error.");
            }
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        dy = (2 * ns) < (n - m) ? c[ns+1] : d[ns--];
        y += dy;
    }
}

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
double qromb (double(*func)(double), double a, double b, double eps) {
    
    double s[JMAXP+1], h[JMAXP+1];
    auto aeps = std::fabs(eps);
    h[1] = 1.0;
    for(int j = 1; j <= JMAX; j++) {
        s[j] = trapzd(func, a, b, j);
        if(j >= K) {
            double ss, dss;
            polint(&h[j-K], &s[j-K], K, 0.0, ss, dss);
            if(std::fabs(dss) < aeps * std::fabs(ss))
                return ss;
        }
        s[j+1] = s[j];
        h[j+1] = 0.25 * h[j];
    }
    throw OpalException (
        "BetMath_integrate::qromb()",
        "Too many steps.");
    return 0.0;
}

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
