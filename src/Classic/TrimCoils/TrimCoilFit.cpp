#include "TrimCoils/TrimCoilFit.h"

#include <cmath>

TrimCoilFit::TrimCoilFit(double bmax,
                         double rmin,
                         double rmax,
                         const std::vector<double>& coefnum,
                         const std::vector<double>& coefdenom):
    TrimCoil(),
    bmax_m(bmax),
    coefnum_m(coefnum),
    coefdenom_m(coefdenom)
{
    // convert to m
    const double mm2m = 0.001;
    rmin_m = rmin * mm2m;
    rmax_m = rmax * mm2m;
}

void TrimCoilFit::doApplyField(const double r, const double z, double *br, double *bz)
{
    if (std::abs(bmax_m) < 1e-20) return;
    // check range, go up to 1 meter outside trim coil
    if (r < rmin_m - 1 || r > rmax_m + 1) return;

    double num     = 0.0; // numerator
    double dnum_dr = 0.0; // derivative of numerator
    // add constant
    num += coefnum_m[0];
    for (std::size_t i = 1; i < coefnum_m.size(); ++i) {
        double dummy = coefnum_m[i] * std::pow(r, i-1);
        num     += dummy * r;
        dnum_dr += dummy * i;
    }

    double denom     = 0.0; // denominator
    double ddenom_dr = 0.0; // derivate of denominator
    // add constant
    denom += coefdenom_m[0];
    for (std::size_t i = 1; i < coefdenom_m.size(); ++i) {
        double dummy = coefdenom_m[i] * std::pow(r, i-1);
        num     += dummy * r;
        dnum_dr += dummy * i;
    }

    double dr  = num / denom;
    // derivative of dr with quotient rule
    double ddr = (dnum_dr * denom - ddenom_dr * num) / (denom*denom);

    *bz -= ddr;
    *br -= dr * z;
}