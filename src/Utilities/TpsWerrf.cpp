#include "Utilities/TpsWerrf.h"
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FTpsMath.h"
#include <cmath>

// Complex error function
// The algorithms is based on:
//   Walter Gautschi,
//   Efficient Computation of the Complex Error Function,
//   SIAM J. Numer. Anal., Vol 7, No. 1, March 1970, pp. 187-198.
// ------------------------------------------------------------------------


void TpsWerrf(const FTps<double, 6> &xx, const FTps<double, 6> &yy,
              FTps<double, 6> &wx, FTps<double, 6> &wy) {
    static const double cc   = 1.12837916709551;  // 2 / sqrt(pi)
    static const double xlim = 5.33;
    static const double ylim = 4.29;
    FTps<double, 6> x = (xx[0] < 0.0) ? -xx : xx;
    FTps<double, 6> y = (yy[0] < 0.0) ? -yy : yy;

    if(y[0] < ylim  &&  x[0] < xlim) {
        // Inside limit rectangle: equation (3.8): q = s(z);
        const FTps<double, 6> q =
            (1.0 - y / ylim) * sqrt(1.0 - (x * x) / (xlim * xlim));

        // Equations (3.11).
        const int nc = 6 + int(23.0 * q[0]);
        const int nu = 9 + int(21.0 * q[0]);
        const FTps<double, 6> h = 1.0 / (3.2 * q);
        const FTps<double, 6> zzx = 1.6 * q + y;    // h - i*z
        const FTps<double, 6> zzy = x;              // h - i*z
        FTps<double, 6> xl = pow(h, 1 - nc);

        // Equations (3.12) for n = nu, nu - 1, ..., N.

        FTps<double, 6> rx[33];
        FTps<double, 6> ry[33];
        rx[nu+1] = 0.0;                     // r_{nu}
        ry[nu+1] = 0.0;                     // r_{nu}
        for(int n = nu; n >= 0; n--) {             // for n = nu, nu-1, ..., 0
            FTps<double, 6> tx = zzx + double(n + 1) * rx[n+1];
            FTps<double, 6> ty = zzy - double(n + 1) * ry[n+1];
            FTps<double, 6> tn = tx * tx + ty * ty;
            rx[n] = 0.5 * tx / tn;            // r_{n-1}
            ry[n] = 0.5 * ty / tn;            // r_{n-1}
        }

        // Equations (3.12) for n = N, N - 1, ..., 0.
        FTps<double, 6> sx;                 // s_{N}
        FTps<double, 6> sy;                 // s_{N}
        for(int n = nc; n >= 0; n--) {      // for n = N, N-1, ..., 0
            FTps<double, 6> saux = sx + xl;
            sx = rx[n] * saux - ry[n] * sy;   // s_{n-1}
            sy = rx[n] * sy + ry[n] * saux;   // s_{n-1}
            xl *= h;
        }

        // Equation (3.13).
        wx = cc * sx;
        wy = cc * sy;
    } else {
        // Outside limit rectangle: equations (3.12) for h = 0.
        const FTps<double, 6> zzx = y;      // - i*z
        const FTps<double, 6> zzy = x;      // - i*z
        FTps<double, 6> rx;                 // r_{nu}
        FTps<double, 6> ry;                 // r_{nu}
        for(int n = 9; n >= 1; n--) {       // for n = nu, nu-1, ..., 0
            FTps<double, 6> tx = zzx + double(n) * rx;
            FTps<double, 6> ty = zzy - double(n) * ry;
            FTps<double, 6> tn = tx * tx + ty * ty;
            rx = 0.5 * tx / tn;               // r_{n-1}
            ry = 0.5 * ty / tn;               // r_{n-1}
        }

        // Equation (3.13).
        wx = cc * rx;
        wy = cc * ry;
    }

    // Equations (3.1).
    if(yy[0] < 0.0) {
        wx =  2.0 * exp(y * y - x * x) * cos(2.0 * x * y) - wx;
        wy = -2.0 * exp(y * y - x * x) * sin(2.0 * x * y) - wy;
        if(xx[0] > 0.0) wy = -wy;
    } else {
        if(xx[0] < 0.0) wy = -wy;
    }
}
