#include "TrimCoils/TrimCoilFit.h"

#include <cmath>

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
    // check range
    if (r < rmin_m || r > rmax_m) return;

    double num     = 0.0; // numerator
    double dnum_dr = 0.0; // derivative of numerator
    double powr    = 1.0; // power of r

    // add constant
    num += coefnum_m[0];
    for (std::size_t i = 1; i < coefnum_m.size(); ++i) {
        dnum_dr += coefnum_m[i] * powr * i;
        powr    *= r;
        num     += coefnum_m[i] * powr;
    }

    double denom     = 0.0; // denominator
    double ddenom_dr = 0.0; // derivative of denominator
    powr             = 1.0; // power of r

    // add constant
    denom += coefdenom_m[0];
    for (std::size_t i = 0; i < coefdenom_m.size(); ++i) {
        ddenom_dr += coefdenom_m[i] * powr * i;
        powr      *= r;
        denom     += coefdenom_m[i] * powr;
    }

    double btr = num / denom;

    // derivative of dr with quotient rule
    double dr = (dnum_dr * denom - ddenom_dr * num) / (denom*denom);

    //std::cout << "r " << r << " dr " <<  dr << std::endl;

    *bz -= bmax_m * btr;
    *br -= bmax_m * dr * z;
}