#include "TrimCoils/TrimCoilFit.h"

#include <cmath>
#include <iostream>

TrimCoilFit::TrimCoilFit(double bmax,
                         double rmin,
                         double rmax,
                         const std::vector<double>& coefnum,
                         const std::vector<double>& coefdenom):
    TrimCoil(),
    coefnum_m(coefnum),
    coefdenom_m(coefdenom)
{
    // convert to m
    const double mm2m = 0.001;
    rmin_m = rmin * mm2m;
    rmax_m = rmax * mm2m;
    // convert to kG
    bmax_m   = bmax * 10.0;

    // normal polynom if no denominator coefficients (denominator = 1)
    if (coefdenom_m.empty())
      coefdenom_m.push_back(1.0);
}

void TrimCoilFit::doApplyField(const double r, const double z, double *br, double *bz)
{
    if (std::abs(bmax_m) < 1e-20) return;
    // check range, go up to 1 meter outside trim coil
    if (r < rmin_m - 1 || r > rmax_m + 1) return;

    double num  = 0.0; // numerator
    double powr = 1.0; // power of r
    for (std::size_t i = 0; i < coefnum_m.size(); ++i) {
        num  += coefnum_m[i] * powr;
        powr *= r;
    }

    double denom = 0.0; // denominator
    powr         = 1.0; // power of r
    for (std::size_t i = 0; i < coefdenom_m.size(); ++i) {
        denom += coefdenom_m[i] * powr;
        powr  *= r;
    }

    double dr  = num / denom;

    // FIXME
    // integral of dr
    double btr = 0.0;

    //std::cout << "r " << r << " dr " <<  dr << std::endl;

    *bz += bmax_m * btr;
    *br += bmax_m * dr * z;
}